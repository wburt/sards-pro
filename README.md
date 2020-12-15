# Automated Status Tool
![img](https://img.shields.io/badge/Lifecycle-Dormant-ff7f2a)
## Usage

sards_geoprocessing_module AOI
AOI : ESRI FEATURESET or ESRIJSON

### Description

This script is reports on threats to species at risk that intersect species range within
a user defined area of interest.

### Dependencies
index.csv  -- csv table linking gis data source to an indicator and threat
threats.csv -- csv table table identifying threat rating for each species at risk
result.txt  -- jinga2 html template for results reporting

## Requirements
arcpy

## License
    Copyright 2015-2016 Province of British Columbia

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at 

       http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
