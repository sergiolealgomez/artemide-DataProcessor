#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 20:15:46 2019

Program that parse various DY data files to ADP-frendly formal

@author: vla18041
"""

import sys
sys.path.append("/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/")

import DataProcessor.Point
import DataProcessor.DataSet

path_to_save="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/OtherPrograms/LHC_impact_Study/"

M_Z=91.### mass of Z-boson

#%%
Qbins=[[30.,45.],[45.,70.],[70.,110.],[110.,150.],[150.,200.],[200.,270.],[270.,350.],[350.,440.],[440.,540.]]
ybins=[[0.,0.4],[0.4,0.8],[0.8,1.2],[1.2,1.6],[1.6,2.0],[2.0,2.4]]
ptbins=[[0.01,1.],[1.,2.],[2.,3.],[3.,4.],[4.,5.],[5.,6.],[6.,7.],[7.,8.],[8.,10.],[10.,12.],[12.,14.],[14.,16.],[16.,20.],[20.,24.],
        [24.,28.],[28.,36.],[36.,42.],[42.,50.],[50.,60.],[60.,75.],[75.,90.]]

#%%
DataCurrent=DataProcessor.DataSet.DataSet('LHC13TeV-Empty',"DY")
DataCurrent.comment="Empty set for LHC at 13 TeV"
DataCurrent.reference="AV"

DataCurrent.isNormalized=False
proc_current=[1,1,5]
s_current=13000.**2
Q_current=[66.,116.]
y_current=[0.,0.4]
incCut=True
cutParam=[20.,20.,-2.4,2.4]

i=0
for QQ in Qbins:
    for yy in ybins:
        for pp in ptbins:
            # makeup a point
            p=DataProcessor.Point.CreateDYPoint(str(i))
            p["process"]=proc_current
            p["s"]=s_current
            p["qT"]=pp
            p["thFactor"]=2/(p["qT"][1]-p["qT"][0])#devide by bin size## since we have to symmetrize y
            p["Q"]=QQ
            p["y"]=yy
            p["includeCuts"]=incCut
            p["cutParams"]=cutParam
            p["xSec"]=100.
            p["uncorrErr"].append(1.)
            DataCurrent.AddPoint(p)
            
            i+=1

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
