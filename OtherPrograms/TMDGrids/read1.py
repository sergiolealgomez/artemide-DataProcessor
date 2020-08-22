#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 22 14:51:04 2020

Reading the grids for TMDs

@author: vla18041
"""


testFile='/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/OtherPrograms/TMDGrids/Grids/InitalExample/grid_pdf_0001.yaml'

import yaml

with open(testFile, 'r') as stream:
    try:
        res=yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)
    
TMDs=res["TMDs"]
#%%

print(len(res["Qg"]))
print(len(res["xg"]))
print(len(res["qToQg"]))