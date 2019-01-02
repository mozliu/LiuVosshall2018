
# LiuVosshall2019

This repository contains all data presented in Liu & Vosshall, 2019 ("General visual and contingent thermal cues interact to elicit attraction in female Aedes aegypti mosquitoes"), and the scripts used to process the data. Each folder corresponds to a group of figures in the paper sharing an assay and processing steps. Within each folder, raw data are stored as LiuVosshall_Raw_&ast;.zip. Raw data were processed using &ast;.m, &ast;.py, and &ast;.ipynb scripts, and outputted processed data are stored as LiuVosshall_Data_&ast;.csv. Processed data can be graphed to replicate the plots within the paper using Fig&ast;.ipynb notebooks.

All &ast;.py and &ast;.ipynb files depend on Python3 modules numpy, matplotlib, and csv. All &ast;.m files depend on the MATLAB Image Processing toolbox. Specific dependencies are noted per script below, and versions of all software are listed at the bottom of this file.

## Master files in root:

* [LiuVosshall_Data.xlsx](LiuVosshall_Data.xlsx): All processed data, collated from LiuVosshall_Data_&ast;.csv files
* [LiuVosshall_Raw_Dir.xlsx](LiuVosshall_Raw_Dir.xlsx): Directory of all raw data, connecting raw file names to conditions and details for each experiment

## Folders and their contents

### [Fig1+S1_magnotether](Fig1+S1_magnotether)

Mosquito and fly data obtained using the magnotether, graphed in Figures 1 and S1.

* [traxio.py](/Fig1+S1_magnotether/traxio.py)
  * Reads in magnotether data
  * [Source](https://github.com/motmot/flytrax/blob/master/motmot/flytrax/traxio.py) (Straw & Dickinson, 2009)
* [magno.py](/Fig1+S1_magnotether/magno.py)
  * Python3 functions for processing fly data
  * Depends on: [neo.io.WinEdrIO](https://github.com/NeuralEnsemble/python-neo/blob/0.3.0/neo/io/winedrio.py)
* [mozmagno.py](/Fig1+S1_magnotether/mozmagno.py)
  * Python3 functions for processing mosquito data
  * Depends on: [traxio.py](traxio.py), pandas(https://pandas.pydata.org/)
* [Fig1KL+S1DEFG.ipynb](/Fig1+S1_magnotether/Fig1KL+S1DEFG.ipynb)
  * Python3 notebook to interactively plot mosquito magnotether data for dark and light shapes under air and CO2
  * Input: [LiuVosshall_Raw_1KL+S1DEFG.zip](/Fig1+S1_magnotether/LiuVosshall_Raw_1KL+S1DEFG.zip): Python3 pickle files
  * Output:
    * [LiuVosshall_Data_1K.csv](/Fig1+S1_magnotether/LiuVosshall_Data_1K.csv)
    * [LiuVosshall_Data_1L.csv](/Fig1+S1_magnotether/LiuVosshall_Data_1L.csv)
    * [LiuVosshall_Data_S1D.csv](/Fig1+S1_magnotether/LiuVosshall_Data_S1D.csv)
    * [LiuVosshall_Data_S1E.csv](/Fig1+S1_magnotether/LiuVosshall_Data_S1E.csv)
    * [LiuVosshall_Data_S1E.csv](/Fig1+S1_magnotether/LiuVosshall_Data_S1E.csv)
    * [LiuVosshall_Data_S1F.csv](/Fig1+S1_magnotether/LiuVosshall_Data_S1F.csv)
    * [LiuVosshall_Data_S1G.csv](/Fig1+S1_magnotether/LiuVosshall_Data_S1G.csv)
  * Depends on: [mozmagno.py](/Fig1+S1_magnotether/mozmagno.py)
* [Fig1MNOP.ipynb](Fig1MNOP.ipynb)
  * Python3 notebook to interactively plot mosquito magnotether data for mosquitoes fed blood and saline
  * Input: [LiuVosshall_Raw_1MNOP.zip](LiuVosshall_Raw_1MNOP.zip): Python3 pickle files
  * Output:
    * [LiuVosshall_Data_1M.csv](LiuVosshall_Data_1M.csv)
    * [LiuVosshall_Data_1N.csv](LiuVosshall_Data_1N.csv)
    * [LiuVosshall_Data_1O.csv](LiuVosshall_Data_1O.csv)
    * [LiuVosshall_Data_1P.csv](LiuVosshall_Data_1P.csv)
  * Depends on: [mozmagno.py](mozmagno.py)
* [FigS1ABC.ipynb](FigS1ABC.ipynb)
  * Python3 notebook to interactively plot fly magnotether data
  * Input: [LiuVosshall_Raw_S1ABC.zip](LiuVosshall_Raw_S1ABC.zip): WinEDR files
  * Output:
    * [LiuVosshall_Data_S1A.csv](LiuVosshall_Data_S1A.csv)
    * [LiuVosshall_Data_S1B.csv](LiuVosshall_Data_S1B.csv)
    * [LiuVosshall_Data_S1C.csv](LiuVosshall_Data_S1C.csv)
  * Depends on: [magno.py](magno.py)

### [Fig2+3+S2+S3_heatseeking](Fig2+3+S2+S3_heatseeking)

Mosquito data obtained using the heat-seeking assay, graphed in Figures 2, 3, S2, and S3.

#### Producing timeseries and summary data of occupancy on Peltier

* count_pelt_model.m
  * MATLAB function to count mosquitoes in a single frame from the heat-seeking assay using a bagged tree classifier
  * Depends on: hsfit.mat (model including predictor, training data, and test data)
* heatseeking_assay_count.m
  * MATLAB function to count mosquitoes over time of a heat-seeking assay
  * Input: Directory of .tiff images, as contained in Fig2+3+S2+S3_heatseeking/LiuVosshall_Raw_&ast;.zip files
  * Output: &ast;.mat data of counts within directory of .tiff images
  * Depends on: count_pelt_model.m 
* manual_heatseeking_count.m
  * MATLAB function to manually count mosquitoes on the dark dot
  * Input: &ast;.mat data produced by heatseeking_assay_count.m
  * Output: &ast;.mat data incorporating manually generated counts
* Fig2.ipynb
  * Python3 notebook to interactively plot dot vs. blank heat-seeking data
  * Input: LiuVosshall_Raw_2.zip: &ast;topy.mat data within image directories
  * Output:
    * LiuVosshall_Data_2DE.csv
    * LiuVosshall_Data_2F.csv
    * LiuVosshall_Data_2G.csv
    * LiuVosshall_Data_2I.csv
* Fig3.ipynb
  * Python3 notebook to interactively plot *Gr3* mutant heat-seeking data
  * Input: LiuVosshall_Raw_3.zip: &ast;topy.mat data within image directories
  * Output:
    * LiuVosshall_Data_3AB.csv
    * LiuVosshall_Data_3C.csv
    * LiuVosshall_Data_3D.csv
    * LiuVosshall_Data_3F.csv
* 20171021_heatseeking.ipynb
  * Python3 notebook to interactively plot shuffled dot vs. blank heat-seeking data
  * Input: LiuVosshall_Raw_S2.zip: &ast;topy.mat data within image directories
  * Output:
    * LiuVosshall_Data_S2AB.csv
    * LiuVosshall_Data_S2C.csv
    * LiuVosshall_Data_S2D.csv
    * LiuVosshall_Data_S2F.csv

#### Producing heatmaps of spatial occupancy on Peltier

* gather_pos.m
  * MATLAB function to aggregate position data within an experiment
* gather_all_pos.m
  * MATLAB function to aggregate position data across experiments
  * Depends on: gather_pos.m
* heatseeking_heatmap.m
  * MATLAB function to plot heatmap of aggregated position data
* orl_heatmap.m
  * Input: LiuVosshall_Raw_2.zip: &ast;.mat data within directories
  * Output: LiuVosshall_Data_2H.csv
  * Depends on: gather_all_pos.m, heatseeking_heatmap.m
* gr3_heatmap.m
  * Input: LiuVosshall_Raw_3.zip: &ast;.mat data within directories
  * Output: LiuVosshall_Data_3E.csv
  * Depends on: gather_all_pos.m, heatseeking_heatmap.m
* shuffle_heatmap.m
  * Input: LiuVosshall_Raw_S2.zip: &ast;.mat data within directories
  * Output: LiuVosshall_Data_S2E.csv
  * Depends on: gather_all_pos.m, heatseeking_heatmap.m

### Fig3_dwelling

Computing dwell times of mosquitoes in the heat-seeking assay, graphed in Figure 3.

* manual_dwell.m
  * MATLAB function to count landings and takeoffs from heat-seeking assay images
  * Input: LiuVosshall_Raw_3.zip: &ast;.tiff files within image directories
  * Output: LiuVosshall_Data_3I.csv
* manual_dwell_count.m
  * MATLAB function to compute dwell time and landings from landing and takeoff data
* manual_dwell_wrapper.m
  * MATLAB function to compute dwell times and landings across entire experiment
  * Input: LiuVosshall_Data_3I.csv
  * Output:
    * LiuVosshall_Data_3F.csv
    * LiuVosshall_Data_3G.csv
* dwell_heatmap.m
  * MATLAB function to plot heatmap of landings and takeoffs
  * Input: LiuVosshall_Data_3I.csv
  * Depends on: heatseeking_heatmap.m

### Software versions

Software | Version | Source
--- | --- | ---
neo | 0.3.0 | https://pypi.org/project/neo/
pandas | 0.20.1 | https://pandas.pydata.org/


