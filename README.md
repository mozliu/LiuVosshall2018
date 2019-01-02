
# LiuVosshall2019

This repository contains all data presented in Liu & Vosshall, 2019 (visual and thermal mosquito behavior), and the scripts used to process the data. Each folder corresponds to a group of figures in the paper sharing an assay and processing steps. Within each folder, raw data are stored as LiuVosshall_Raw_&ast;.zip. Raw data were processed using &ast;.m, &ast;.py, and &ast;.ipynb scripts, and outputted processed data are stored as LiuVosshall_Data_&ast;.csv. Processed data can be graphed to replicate the plots within the paper using Fig&ast;.ipynb notebooks.

All &ast;.py and &ast;.ipynb files depend on Python3 modules numpy, matplotlib, and csv. All &ast;.m files depend on the MATLAB Image Processing toolbox. Specific dependencies are noted per script below, and versions of all software are listed at the bottom of this file.

## Master files in root:

LiuVosshall_Data.xlsx		All processed data, collated from LiuVosshall_Data_*.csv files
LiuVosshall_Raw_Dir.xlsx	Directory of all raw data, connecting raw file names to conditions and details for each experiment

## Folders and their contents

### Fig1+S1_magnotether

Mosquito and fly data obtained using the magnotether, graphed in Figures 1 and S1.

#### Raw data

* LiuVosshall_Raw_1KL+S1DEFG.zip
  * Mosquito magnotether data for dark and light shapes under air and CO2 conditions
  * Format: Python pickle files
* LiuVosshall_Raw_1MNOP.zip
  * Mosquito magnotether data for mosquitoes fed blood and saline
  * Format: Python pickle files
* LiuVosshall_Raw_S1ABC.zip
  * Fly magnotether data for dark shapes under air conditions
  * Format: WinEDR files

#### Processing scripts

* traxio.py
  * Reads in magnotether data
  * Source: https://github.com/motmot/flytrax/blob/master/motmot/flytrax/traxio.py (Straw & Dickinson, 2009)
* magno.py
  * Python3 functions for processing fly data
  * Depends on: neo.io.WinEdrIO
* mozmagno.py
  * Python3 functions for processing mosquito data
  * Depends on: traxio, pandas
* 20160811_fourshapes.ipynb
  * Python3 notebook for mosquito magnotether data for dark and light shapes under air and CO2
  * Input: LiuVosshall_Raw_1KL+S1DEFG.zip
  * Output:
    * LiuVosshall_Data_1K.csv
    * LiuVosshall_Data_1L.csv
    * LiuVosshall_Data_S1D.csv
    * LiuVosshall_Data_S1E.csv
    * LiuVosshall_Data_S1E.csv
    * LiuVosshall_Data_S1F.csv
    * LiuVosshall_Data_S1G.csv
  * Depends on: mozmagno
* 20160811_fourshapes.ipynb
  * Python3 notebook for mosquito magnotether data for mosquitoes fed blood and saline
  * Input: LiuVosshall_Raw_1MNOP.zip
  * Output:
    * LiuVosshall_Data_1M.csv
    * LiuVosshall_Data_1N.csv
    * LiuVosshall_Data_1O.csv
    * LiuVosshall_Data_1P.csv
  * Depends on: mozmagno
* magno_20150624.ipynb
  * Python3 notebook for fly magnotether data
  * Input: LiuVosshall_Raw_S1ABC.zip
  * Output:
    * LiuVosshall_Data_S1A.csv
    * LiuVosshall_Data_S1B.csv
    * LiuVosshall_Data_S1C.csv

### Software versions

neo 0.3.0	https://github.com/NeuralEnsemble/python-neo/blob/0.3.0/neo/io/winedrio.py


