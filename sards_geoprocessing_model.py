import os
import sys
import ntpath
import arcpy
import json
import csv
import jinja2
import multiprocessing
import logging

# multi processing for geometric operations
cpu_cnt = multiprocessing.cpu_count()

arcpy.env.overwriteOutput = 1
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(3005)

script_path = sys.path[0]
arcpy.AddMessage("Script Path: " + script_path)
acsv = ".\\index.csv"
tcsv = ".\\threats.csv"
indicatorCsv = os.path.join(script_path,'index.csv')
#input geometry as Type = Featureset
#if not on 10.2 server then use testing data
if arcpy.GetParameterAsText(0):
    tAOI = arcpy.GetParameterAsText(0)
    arcpy.AddMessage("AOI is: " + tAOI)
    areaOfInterest = arcpy.FeatureSet(tAOI)
else:
    arcpy.AddWarning('WARNING -- NO INPUT POLYGON DEFAULTING')
    tAOI = os.path.join(script_path,'kb_sards.gdb','aoi_test')
    areaOfInterest = arcpy.FeatureSet(tAOI)
arcpy.AddMessage("Starting SARDS Reporting ---")
#workspace related
src_gdb = os.path.join(script_path,'kb_sards.gdb')
outfile =  os.path.join(arcpy.env.scratchFolder,"thisOutput.html")
#wildlife info
sar_range_draft_wm= os.path.join(src_gdb,'KB_DRAFT_SAR_WM')
biot_occr_non_sens_area_svw_wm='biot_occr_non_sens_area_svw_wm'
wcp_critical_habitat_sp_kb_wm='wcp_critical_habitat_sp_kb_wm'
threat_table = os.path.join(script_path,'threats.csv')

#Threat table looks like this
#'OBJECTID', 'ENG_NAME', 'scientific_name', 'climate', 'cooridor', 'nat_syst_mod', 'agriculture', 'bio_resource_use',
# 'pollution', 'developement', 'human_intrusion', 'invasive_genes', 'climate_weather', 'energy_mining', 'geo_event'

class Aoi:
    #establish an output file
    wrkspace = 'in_memory'
        #class for the input AOI polygon.  Deals with reprojecting from geographic for area calc
    class Poly:
        def __init__(self,featureset):
            arcpy.AddMessage("Building AOI Polygon")
            self.source = featureset
            self.area = self._Build_area()
        def GetHectares(self):
            if self.area:
                return self.area / 10000
            else:
                return None
        def _Build_area(self):
            arcpy.AddMessage("Calculating Area")
            #define area with bc albers projection
            a = self.source
            sr_out = arcpy.SpatialReference(3005)
            sr_in = arcpy.SpatialReference(int(json.loads(a.JSON)['spatialReference']['wkid']))
            trans = arcpy.ListTransformations(sr_in, sr_out)
            area_metres = 0
            for r in arcpy.da.SearchCursor (a, ["SHAPE@"]):
                geom = r[0]
                if sr_in.factoryCode <> sr_out.factoryCode:
                # if the geographic coordinate systems are different, there will be transformations, and we will use the first one
                    if len(trans):
                        area_metres += geom.projectAs(sr_out, trans[0]).getArea('GEODESIC','METERS')
                    # if the geographic coordinate systems are the same, we don't need a transformation, so run the function without
                    else:
                        self.area += geom.getArea('GEODESIC','METERS')
                else:
                    self.area += geom.getArea('GEODESIC','METERS')
            return area_metres

    class IndicatorInputs:
        def __init__(self,csvfile):
            arcpy.AddMessage("Inputs from: " + csvfile)
            self.inputs = []#list of indicator objects
            #build dictionary of inputs from fields indicator,source,threat

            with open(csvfile,'rb') as csvFile:
                ireader = csv.DictReader(csvFile)
                for r in ireader:
                    arcpy.AddMessage("Adding: " + str(r))
                    i = Indicator(r['indicator'],r['source'],r['threat'])
                    self.inputs.append(i)
        def Remove(self,Indicator):
            #remove first indicators with matching source
            for i in self.inputs:
                if i == Indicator:
                    self.inputs.remove(i)
        def ListInputNames(self):#list names of indicators
            l = []
            for i in self.inputs:
                l.append(i.name)
            return list(set(l))
        def ListInputSources(self):#list sources of indicators
            l = []
            for i in self.inputs:
                l.append(i.source)
            return list(set(l))
        def ListInputThreats(self):#list threats of indicators
            l = []
            for i in self.inputs:
                l.append(i.threat)
            return list(set(l))
        def ListInputs(self,searchParam=None):#find a input with matching name,source,threat parameter
            l = []
            if searchParam:
                for i in self.inputs:
                    if searchParam in [i.name,i.source,i.threat]:
                        l.append(i)
            else:
                l = self.inputs
            return l

    def __init__(self,aoi_featureset,workspace,rangefc, indicatorInputCSV,occurance_fc, critical_habitat_fc):
        arcpy.AddMessage("Constructing AOI...")
        arcpy.env.workspace = workspace
        self.input = self.Poly(aoi_featureset) #aoi_poly class
        self.hectares = self.input.GetHectares()
        self.indicatorInputs = self.IndicatorInputs(indicatorInputCSV)
        self.__filterIndicators()
        self.species = []
        self.summary = {} # dictionary of indicator: area/length
        self.summary_tables_info = {} #indicator: summary
        self.clipedIndicatorList = self.__Clip_fcList(self.indicatorInputs.ListInputSources()) #list of indicator fc inside AOI

        self.occurance_fc = occurance_fc
        self.critical_habitat_fc = critical_habitat_fc
        self.range = self.__Clip_fcList([rangefc])[0]
        #self.__Species_Maker()
        self.__Species_Maker()

        #self.__GenerateResults()
    def __Species_Maker(self):
        #unique species from clipped range data
        values = [row[0] for row in arcpy.da.SearchCursor(self.range, "ENG_NAME")]
        uniqueValues = set(values)
        #for each species make an object
        for s in uniqueValues:
            print 'Create Species: ' + s
            sfix = filter(str.isalnum, str(s))
            if '\'' in s:
                sLayer = arcpy.MakeFeatureLayer_management(self.range,'lyr','"ENG_NAME" =' + '\''+ s.replace('\'','\'\'') + '\'' )
            else:
                sLayer = arcpy.MakeFeatureLayer_management(self.range,'lyr','"ENG_NAME" =' + '\''+ s + '\'' )
            sfc = arcpy.CopyFeatures_management(sLayer,"in_memory/range_"+sfix)

            nSpecies = Species(s,sfc.getOutput(0),self.indicatorInputs,threat_table)
            self.species.append(nSpecies)
    def __mp_Species_Maker(self):
        #multiprocessing enabled
        pool = multiprocessing.Pool(cpu_cnt-2)
        print "Using multiprocessing with " + str(cpu_cnt-2) + " CPU CORES"
        def merge_params(a, b, c, d):
            return '{} & {}'.format(a, b, c, d)

        def merge_params_unpack(args):
            return merge_params(*args)
        #unique species from clipped range data
        values = [row[0] for row in arcpy.da.SearchCursor(self.range, "ENG_NAME")]
        uniqueValues = set(values)
        #for each species make an object
        il = []
        for s in uniqueValues:
            print 'Create Species: ' + s
            sfix = filter(str.isalnum, str(s))
            if '\'' in s:
                sLayer = arcpy.MakeFeatureLayer_management(self.range,'lyr','"ENG_NAME" =' + '\''+ s.replace('\'','\'\'') + '\'' )
            else:
                sLayer = arcpy.MakeFeatureLayer_management(self.range,'lyr','"ENG_NAME" =' + '\''+ s + '\'' )
            sfc = arcpy.CopyFeatures_management(sLayer,"in_memory/range_"+sfix)
            il.append([s,sfc.getOutput(0),self.indicatorInputs,threat_table])
        extResults = pool.map(self.__mp_process,il)
        pool.close()
        pool.join()
        for mp_s in extResults:
            self.species.append(mp_s)
    def __mp_process(self,inputList):
        mp_s = Species(inputList[0],inputList[1],inputList[2],inputList[3])
        return mp_s
    def __filterIndicators(self): #filters the input indicators to just those that overlap AOI
        srcs = self.indicatorInputs.inputs
        arcpy.AddMessage(str(srcs))
        for i in srcs:
            l= arcpy.MakeFeatureLayer_management(i.source,'lyr')
            arcpy.SelectLayerByLocation_management(l,'intersect',self.input.source)
            cnt = int(arcpy.GetCount_management(l).getOutput(0))
            #rl = self.__Clip_fcList([i.source])
            if cnt == 0:
                #no intersection with AOI
                self.indicatorInputs.Remove(i)
    def __getArea__(self,fcLayer):
        fArea = 0.0 #area in hectares
        d = arcpy.Describe(fcLayer)
        if d.shapeType == 'Polygon':
            sc = arcpy.da.SearchCursor(fcLayer,['SHAPE@AREA'])
            for r in sc:
                fArea = fArea + r[0]/10000
        return fArea
    def __getLenght(self,fcLayer):
        fLength = 0.0 #Lenght KM
        d = arcpy.Describe(fcLayer)
        if d.shapeType == 'Polyline':
            sc = arcpy.da.SearchCursor(fcLayer,['SHAPE@LENGTH'])
            for r in sc:
                fLength = fLength + r[0]/1000
        return fLength
    def __Clip_fcList(self, fclist):
        clpList = []

        for fc in fclist:
            output = 'rslt_'+ ntpath.basename(fc)

            arcpy.Clip_analysis(fc,self.input.source, self.wrkspace + '/'+ output)
            cr = arcpy.GetCount_management(self.wrkspace + '/'+ output)
            if int(cr.getOutput(0))>0:
                clpList.append(self.wrkspace + '/'+ output)
        return clpList
class Indicator:

    def __init__(self, name, source, threat):
        arcpy.AddMessage("Constructing indicator: " + name)
        self.name = name
        self.source = None
        self.threat = threat
        self.value = None
        self.units = None
        self.UpdateSource(source)
    def UpdateSource(self, newSource):
        #update source, trigger area calc
        if arcpy.Exists(newSource):
            self.source = newSource
            d = arcpy.Describe(newSource)
            if d.shapeType == 'Polygon':
                self.value = self.__getArea__(newSource)
                self.units = 'ha'
            elif d.shapeType == 'Polyline':
                self.value = self.__getLength__(newSource)
                self.units = 'km'
        else:
            print "Source " + newSource + " does not exist"
    def __getArea__(self,fcLayer):
        fArea = 0.0 #area in hectares
        d = arcpy.Describe(fcLayer)
        if d.shapeType == 'Polygon':
            sc = arcpy.da.SearchCursor(fcLayer,['SHAPE@AREA'])
            for r in sc:
                fArea = fArea + r[0]/10000
        return fArea
    def __getLength__(self,fcLayer):
        fLength = 0.0 #Lenght KM
        d = arcpy.Describe(fcLayer)
        if d.shapeType == 'Polyline':
            sc = arcpy.da.SearchCursor(fcLayer,['SHAPE@LENGTH'])
            for r in sc:
                fLength = fLength + r[0]/1000
        return fLength
#class deals with individual species / species specific reporting
#species name, area of range(ha), threat_lu_table(threat severity by species), indicator threat lookup dictionary
class Species:
    def __init__(self, species_eng_name, range_featureclass, AOIindicators, threat_csv):
        arcpy.AddMessage("Creating species: " + species_eng_name)
        #name of species
        self.name = species_eng_name
        self.source = range_featureclass
        #total range area (ha)
        self.area_hectares = self.__getArea__(self.source)
        self.indicators = [] # list of indicators
        self.threats = {} #{threat:severity}
        self.threat_lookup_dictionary = {} #all possible threat:severity values for species
        self.threatFactory(threat_csv)
        self.indicatorFactory(AOIindicators.inputs)
    def __getArea__(self,fcLayer):
        fArea = 0.0 #area in hectares
        d = arcpy.Describe(fcLayer)
        if d.shapeType == 'Polygon':
            sc = arcpy.da.SearchCursor(fcLayer,['SHAPE@AREA'])
            for r in sc:
                fArea = fArea + r[0]/10000
        return fArea
    #Add an indicator intersecting species range eng_name,area value, area units
    def indicatorFactory(self, inidcatorList):
        sfix = filter(str.isalnum, str(self.name))
        #collect indicators for each species
        for i in inidcatorList:
            basename = ntpath.basename(i.source)
            output =  'ind_' + basename + '_'+ sfix
            rslt = arcpy.Clip_analysis(i.source,self.source, os.path.join("in_memory",output))
            output = rslt.getOutput(0)
            cr = arcpy.GetCount_management(rslt)
            if int(cr.getOutput(0))>0:
                speciesIndicator = Indicator(i.name,output,i.threat)
                self.indicators.append(speciesIndicator)
                if i.threat in self.threat_lookup_dictionary:
                    self.threats[i.threat] = self.threat_lookup_dictionary[i.threat]


    def threatFactory(self,threat_csv):
        #builds dictionary of threats:severity for species
        with open(threat_csv,'rb') as csvfile:
            ireader = csv.DictReader(csvfile)
            threatList = ['climate','cooridor','nat_syst_mod','agriculture','bio_resource_use','pollution',\
            'developement','human_intrusion','invasive_genes','climate_weather','energy_mining','geo_event']
            threatList = ['BC vulnerability to climate change','Transportation & service corridors',\
            'Natural system modifications','Agriculture & aquaculture','Biological resource use','Pollution',\
            'Residential & commercial development','Human intrusions and disturbance',\
            'Invasive & other problematic species & genes','Climate change & severe weather',\
            'Energy production & mining','Geologic events']
            for r in ireader:

                if str(r['ENG_NAME'])==self.name:#only use the row for matching species
                    for k in r:#look at each column
                        if k in threatList:#field is a threat category
                            if r[k]:
                                self.threat_lookup_dictionary[k]=r[k]

def main():
    
    a = Aoi(areaOfInterest,script_path,sar_range_draft_wm,indicatorCsv,biot_occr_non_sens_area_svw_wm,wcp_critical_habitat_sp_kb_wm)
    arcpy.AddMessage("Starting to write report...")
    #build summary
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(
        searchpath=script_path)
        )
    template = env.get_template('result.txt')
    # build up some summary info
    ol = {} # species range info species:areaha
    iDict = {}# indicator: max value
    threatCounts = {} # {threat:{severity:count}} all species
    for s in a.species:
        ol[s.name] = s.area_hectares
        for i in s.indicators:
            if i.name in iDict:
                if iDict[i.name][0] < round(i.value,2):
                    iDict[i.name] = [round(i.value,2), i.units]
            else:
                iDict[i.name] = [round(i.value,2), i.units]
        for t in s.threats:#{threat:severity}
            #each threat

            if t in threatCounts:#threat has severities
                sevDict = threatCounts[t] #{severity:counts}
                sv = s.threats[t]
                if sv in sevDict: #severity has count
                    sevDict[s.threats[t]] = sevDict[s.threats[t]] + 1 #add count to severity for this threat
                else:#severity has no count
                    sevDict[sv] = 1
            else: # threat hasn't been counted yet
                ns = s.threats[t]
                sevDict = {ns:1}
            threatCounts[t] = sevDict

    threatSummary =[]
    for t in threatCounts:
        outstr = '<p><u>'+t+'</u></p>\n<p>Threat Level   |   Count </p>'
        for s, c in threatCounts[t].iteritems():
            outstr = outstr +  '\n<p>' + s +'\t|\t' + str(c) + '</p>'
        threatSummary.append(outstr)
    forHtml = '\n'.join(threatSummary)


    iList = []
    sList =[]
    for k in iDict:
        iList.append(k + ' ' + str(iDict[k][0])+ ' ' + iDict[k][1])
    for s in ol:
        sList.append(s + ' ' + str(ol[s]) + ' ha')
    ahtml = template.render(aoi_hectares = round(a.hectares,1),threatString=forHtml, myIndicators = iList, mylist = ol.keys() )
    with open(outfile, 'w') as f:
        f.write(ahtml)
    #the last hurah!
    arcpy.SetParameterAsText(1, outfile)

if __name__== "__main__":
    main()


