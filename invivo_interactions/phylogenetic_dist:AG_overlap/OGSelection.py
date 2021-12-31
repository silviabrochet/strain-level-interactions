#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 09:50:47 2021

@author: garancesarton-loheac
"""


# extract orthogroups from a clade
#%%
import numpy as np
import pandas as pd
from os import path
import sys



#%%
# Set the variables


pathToOG=sys.argv[1] # orthogroups.tsv file
PathToSAVE=sys.argv[2] # where to save result

#selectedClade='stinglessbee'
#pathToTsv='~/Documents/Teaching/SAGE_2021/SAGE_II/GenomesDataset_tsv.txt'
#pathToOG='~/Documents/Teaching/SAGE_2021/SAGE_II/OrthoFinder/200412_Orthofinder_Result/Orthogroups/Orthogroups.tsv'#
#PathToSAVE='~/Documents/Teaching/SAGE_2021/SAGE_II/test.csv'

#%%
## prepare OG table

OGdf = pd.read_csv(pathToOG, sep='\t', header=0) # read the table


NumberOfSamples=OGdf.shape[1]-1 #First column is OG name = number of sample = dim-1
NumberofEmptyOG=np.array(OGdf.apply(lambda x : x.isnull().sum(), axis=1)) # count number of empty orthogroups
selectEmpty=np.where(NumberofEmptyOG >= 1)[0].tolist() #select rows with selected number of genomes

# Subsample dataframe, keep all columns, select rows with choosen number of sp in OG
selectedOGs=OGdf.iloc[selectEmpty]
selectedOGs.columns = OGdf.columns
# save file
selectedOGs.to_csv(PathToSAVE, sep='\t', index=False, header=True)


#%% Print useful data

print('Orthogroup Table dimension \n * Number of Genomes: %d ; \n * Number of Orthogroups:  %d' %(OGdf.shape[1]-1,OGdf.shape[0]))
print('Selected Orthogroup table \n * Number of Genomes: %d ; \n * Number of Orthogroups:  %d' %(selectedOGs.shape[1]-1, selectedOGs.shape[0]))
