#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 16:40:59 2019

@author: vla18041
"""

#######################################
# importing libraries
#######################################

import sys
import time
import numpy
sys.path.append("/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor")
import DataProcessor.harpyInterface
import DataProcessor.DataMultiSet

MAINPATH="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/FittingPrograms/Vpion19/"
#%%
#######################################
#Initialize artemide with a replica -2
#######################################
import harpy
path_to_constants=MAINPATH+"Constants-files/"
#harpy.initialize(path_to_constants+"DY_nnlo/const-Vpion19_nnlo")
harpy.initialize(path_to_constants+"const-Vpion19_nnlo_all=0")

#### All=0 Case
harpy.setNPparameters_TMDR([2., 0.0398333])
harpy.setNPparameters_uTMDPDF([0.184739, 6.22437, 588.193, 2.44327, -2.51106, 0.,  0.17, 0.48, 2.15])


#%%
### read the list of files and return the list of DataSets
def loadThisData(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/DataLib/piDY/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   

    return dataCollection



#%%
##################Cut function
def cutFunc(p):
    
    delta=p["<qT>"]/p["<Q>"]
    
    if(p["id"][0] == "E"):
        delta=p["<qT>"]/p["Q"][1] 
    
    return (delta<0.1 or delta<0.3) , p

#%%
### Loading the data set
theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisData([
                                'E615(pi)-dxF-0.0','E615(pi)-dxF-0.1','E615(pi)-dxF-0.2','E615(pi)-dxF-0.3',
                                'E615(pi)-dxF-0.4','E615(pi)-dxF-0.5','E615(pi)-dxF-0.6','E615(pi)-dxF-0.7'
                                ]))

setDY=theData.CutData(cutFunc) 

print('Loaded ', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
print('Loaded experiments are', [i.name for i in setDY.sets])

#%%
DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)
    
#%%
#######################################
# Minimisation
#######################################
totalN=setDY.numberOfPoints

def chi_2(x):
    startT=time.time()
    
    nn=[0.184739, 6.22437, 588.193, 2.44327, -2.51106, 0.,  x[0],x[1],x[2]]    
    harpy.setNPparameters_uTMDPDF(nn)
    print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
    
    ccDY2,cc3=DataProcessor.harpyInterface.ComputeChi2(setDY)
    
    cc=ccDY2/totalN
    endT=time.time()
    print(':->',cc,'       t=',endT-startT)
    return ccDY2

#%%

from iminuit import Minuit

initialValues=(0.01, 0.40, 3.44)

initialErrors=(0.1,       0.3,     1. )
searchLimits=((0,None), (0,None), (0,None))
parametersToMinimize=(False, False,  False)

m = Minuit.from_array_func(chi_2, initialValues,
      error=initialErrors, limit=searchLimits, fix=parametersToMinimize, errordef=1)

#m.get_param_states()

m.tol=0.0001*totalN*10000 ### the last 0.0001 is to compensate MINUIT def
m.strategy=1

#%%


# m.tol=0.0001*totalN*10000 ### the last 0.0001 is to compensate MINUIT def
# m.strategy=1

# m.migrad()

# print(m.params)


#%%
def MinForReplica():
    
    
    def repchi_2(x):        
        startT=time.time()
        nn=[0.184739, 6.22437, 588.193, 2.44327, -2.51106, 0.,  x[0],x[1],x[2]]    
        harpy.setNPparameters_uTMDPDF(nn)
        print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
        
        ccDY2,cc3=DataProcessor.harpyInterface.ComputeChi2(repDataDY)
        
        cc=ccDY2/totalNnew
        endT=time.time()
        print(':->',cc,'       t=',endT-startT)
        return ccDY2
    
    repDataDY=setDY.GenerateReplica()
    totalNnew=repDataDY.numberOfPoints
    
    localM = Minuit.from_array_func(repchi_2, initialValues,
      error=initialErrors, limit=searchLimits, fix=parametersToMinimize, errordef=1)
    
    localM.tol=0.0001*totalNnew*10000 ### the last 0.0001 is to compensate MINUIT def
    m.strategy=1

    localM.migrad()
    
    return [localM.fval,localM.values.values()]

#%%
#
# Generate pseudo data and minimise   100 times
#
numOfReplicas=50
REPPATH=MAINPATH+"LOGS/"+"Vpion19_REPLICAS_all=0.txt"
for i in range(numOfReplicas):
    print('---------------------------------------------------------------')
    print('------------REPLICA ',i,'/',numOfReplicas,'--------------------')
    print('---------------------------------------------------------------')
    repRes=MinForReplica()
    print(repRes)
    f=open(REPPATH,"a+")
    print('SAVING >>  ',f.name)
    f.write(str(repRes)+"\n")
    f.close()