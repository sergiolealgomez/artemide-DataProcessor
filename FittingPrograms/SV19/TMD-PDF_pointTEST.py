#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 10:03:22 2020

@author: vla18041
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: vla18041
"""

#######################################
# Paths
#######################################
PathToHarpy="/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/harpy"
PathToDataProcessor="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/"
PathToDataLibrary=PathToDataProcessor+"DataLib/unpolDY/"
PathToLog=PathToDataProcessor+"FittingPrograms/SV19/LOGS/"
PathToConstantsFile=PathToDataProcessor+"FittingPrograms/SV19/Constants-files/DY_nnlo/const-NNPDF31_NNLO"

import sys
#sys.path.append(PathToHarpy)
sys.path.append(PathToDataProcessor)

#### Starting and final replica (included)
StartReplica=1
FinalReplica=3

#%%
#######################################
# importing libraries
#######################################
import numpy
import time

import DataProcessor.harpyInterface
import DataProcessor.DataMultiSet

#%%
#######################################
# Output paths
#######################################
import socket
PCname=socket.gethostname()

replicaFile=PathToLog+"SV19_NNPDF31_pointsTEST.txt"
logFile=PathToLog+PCname+".log"

#%%
#######################################
#Initialize artemide
#######################################
import harpy

harpy.initialize(PathToConstantsFile)
initializationArray=[2.0340, 0.0299, 0.2512, 7.7572,334.6108, 2.4543,-4.8203, 0.1000,  0.0000]
harpy.setNPparameters(initializationArray)

#%%
#######################################
# read the list of files and return the list of DataSets
#######################################
def loadThisData(listOfNames):    
    import DataProcessor.DataSet
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(PathToDataLibrary+name+".csv")
        dataCollection.append(loadedData)   

    return dataCollection

#%%
#######################################
# Data cut function
#######################################
def cutFunc(p):    
    par=1.0
    if p["type"]=="DY":
        if(p["xSec"]>0):
            err=numpy.sqrt(sum([i**2 for i in p["uncorrErr"]]))/p["xSec"]
        else:
            err=100.
        delta=p["<qT>"]/p["<Q>"]
        
        if(p["id"][0] == "E"):
            delta=p["<qT>"]/p["Q"][1] 
        
        
        if(p["id"][0:4] == "E605"):
            if(p["Q"][0]==10.5):#UPSILON resonance-bin
                return False , p
        elif(p["id"][0:4] == "E772"):
            if(p["Q"][0]<10):#these bins seems broken
                return False , p
        elif(p["id"][0:4] == "E615"):
            if(9<p["<Q>"]<11.2):#UPSILON resonance-bin
                return False , p
        elif(p["id"][0:4] == "E228"):
            if(9<p["<Q>"]<11):#UPSILON resonance-bin
                return False , p
        else:
            if(9<p["<Q>"]<11):#UPSILON resonance-bin
                return False , p
    
    if p["type"]=="SIDIS":        
        if p["<z>"]>0.8:
            return False , p
        ## bins with low z drop
        if p["<z>"]<0.2:
            return False , p
        
        par=1.0
        if p["xSec"]<0.00000001:
            err=1
            delta=1
        else:
            ##############3 I MULTIPLY THE ERROR BY 100 (so it does not affect the cuts)
            err=10000#*numpy.sqrt(p.uncorrErrorsSquare)/p.xSec    
            gamma2=(2.0*p["M_target"]*p["<x>"]/p["<Q>"])**2
            rho2=(p["M_product"]/p["<z>"]/(p["<Q>"]))**2
            qT=p["<pT>"]/p["<z>"]*numpy.sqrt((1+gamma2)/(1-gamma2*rho2))
            delta=qT/(p["<Q>"])
            
            ### compute the largest possible qT (approximate)
            gamma2WORST=(2.0*p["M_target"]*p["x"][1]/p["<Q>"])**2
            # it is definitely not a TMD point
            if gamma2WORST*rho2>1:
                return False , p
            qTWORST=p["pT"][1]/p["z"][0]*numpy.sqrt((1+gamma2WORST)/(1-gamma2WORST*rho2))
    
            ## drop if qT>Q/2
            if qTWORST>p["<Q>"]/2:
                return False , p
    
        ### drop Q<2
        if p["<Q>"]<2 :
            return False , p
    
#    return delta<0.5 and p.qT_avarage<80
    return (delta<0.1 or (delta<0.25 and par/err*delta**2<1)) , p

#%%
#######################################
# Loading the data set
#######################################
theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisData(['CDF1', 'CDF2', 'D01', 'D02', 'D02m', 
                      'A7-00y10', 'A7-10y20','A7-20y24', 
                      'A8-00y04', 'A8-04y08','A8-08y12',
                      'A8-12y16', 'A8-16y20','A8-20y24',
                      'A8-16y20','A8-116Q150',
                      'CMS7', 'CMS8', 
                      'LHCb7', 'LHCb8', 'LHCb13', 
                      'PHE200', 'E228-200', 'E228-300', 'E228-400', 
                      'E772', 'E605']))

setDY=theData.CutData(cutFunc) 

print('Loaded ', setDY.numberOfSets, 'data sets with ', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
print('Loaded experiments are', [i.name for i in setDY.sets])

#%%
for s in setDY.sets:
    s.isNormalized=True
    s.normalizationMethod='bestChi2'

#%%
f = open("/misc/data2/braun/vla18041/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/LOGS/NNPDF31/allPDF_centralEXP_full.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

#%%
import ast

replicas=[]
for l in data_from_f:
    replicas.append(ast.literal_eval(l))


#%%
f = open("/misc/data2/braun/vla18041/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/LOGS/NNPDF31/allPDF_centralEXP_full_POINTStest.dat",'a')
f.write(str([p["id"] for p in setDY.points]))
f.write("\n")
f.close()
#%%

for j in range(655,867):    
    r=replicas[j]
    print('>>> DOING >>>> ',j,'  = ',r[3])    
    harpy.setPDFreplica(r[3])
    harpy.setNPparameters(r[4])
    YY=DataProcessor.harpyInterface.ComputeXSec(setDY)
    RR=[]
    for i in range(len(YY)):
        p=setDY.points[i]
        RR.append((p["xSec"]-YY[i])**2/numpy.dot(p["uncorrErr"],p["uncorrErr"]))
    
    f = open("/misc/data2/braun/vla18041/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/LOGS/NNPDF31/allPDF_centralEXP_full_POINTStest.dat",'a')
    f.write(str(RR))
    f.write("\n")
    f.close()
    