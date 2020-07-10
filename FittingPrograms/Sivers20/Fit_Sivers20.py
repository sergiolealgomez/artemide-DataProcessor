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
#sys.path.append("/home/m/Github/artemide-DataProcessor")
sys.path.append("/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor")
import DataProcessor.harpyInterface
import DataProcessor.DataMultiSet

#MAINPATH="/home/m/Github/artemide-DataProcessor/"
MAINPATH="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/"

#%%
#######################################
#Initialize artemide
#######################################
import harpy
path_to_constants=MAINPATH+"FittingPrograms/Sivers20/Constants-files/"
#harpy.initialize(path_to_constants+"const-Sivers20_lo")
harpy.initialize(path_to_constants+"const-Sivers20_nnlo")
harpy.setNPparameters_TMDR([1.93, 0.0434])
harpy.setNPparameters_uTMDPDF([0.253434, 9.04351, 346.999, 2.47992, -5.69988, 0.1, 0.17, 0.48, 2.15])
harpy.setNPparameters_uTMDFF([0.264,0.479,0.459,0.539])
harpy.setNPparameters_SiversTMDPDF([5.2, 0., -0.6, 15.9, 0.5, -0.2, 21.6, -0.5, -0.1, 0.4, -1.1]) 
#%%
### read the list of files and return the list of DataSets
def loadThisData(listOfNames):    
    import DataProcessor.DataSet.
    
    #path_to_data="/home/m/Github/artemide-DataProcessor/DataLib/Sivers/"
    path_to_data="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/DataLib/Sivers/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   

    return dataCollection

#%%
#################### LOG save function
LOGPATH=MAINPATH+"FittingPrograms/LOGS/"+"Sivers20["+time.ctime()+"].log"
def SaveToLog(logTitle,text):
    with open(LOGPATH, 'a') as file:
        file.write(time.ctime())
        file.write(' --> '+logTitle+'\n')
        file.write(text)
        file.write('\n \n \n')

#%%
##################Cut function
def cutFunc(p):
    import copy
    if p["type"]=="DY":
        delta=p["<qT>"]/p["<Q>"]        
        
        if(9<p["<Q>"]<11):#UPSILON resonance-bin
            return False , p
    
    if p["type"]=="SIDIS":  
        
        delta=p["<pT>"]/p["<z>"]/p["<Q>"]        
    
    deltaTEST=0.3
    
    if delta<deltaTEST:
        pNew=copy.deepcopy(p)    
        pNew["process"]=pNew["weightProcess"]
        if p["type"]=="SIDIS":
            normX, normX1=DataProcessor.harpyInterface.ComputeXSec(pNew,method="central")        
        elif p["type"]=="DY":
            normX, normX1=DataProcessor.harpyInterface.ComputeXSec(pNew)        
        else:
            print("Are you crazy?")
        p["thFactor"]=p["thFactor"]/normX        
    
#    return delta<0.5 and p.qT_avarage<80
    return delta<deltaTEST, p

#%%
### Loading the data set
# theData=DataProcessor.DataMultiSet.DataMultiSet("SIDISset",loadThisData([
#                     'compass.sivers.pi+.dpt', 'compass.sivers.pi-.dpt',
#                     'compass.sivers.k+.dpt', 'compass.sivers.k-.dpt',
#                     'compass.sivers.pi+.dx', 'compass.sivers.pi-.dx',
#                     'compass.sivers.k+.dx', 'compass.sivers.k-.dx',
#                     'compass.sivers.pi+.dz', 'compass.sivers.pi-.dz',
#                     'compass.sivers.k+.dz', 'compass.sivers.k-.dz',                    
#                     'hermes.sivers.pi+.Qint.dpt','hermes.sivers.k+.Qint.dpt',
#                     'hermes.sivers.pi-.Qint.dpt','hermes.sivers.k-.Qint.dpt',
#                     'hermes.sivers.pi+.Qint.dx','hermes.sivers.k+.Qint.dx',
#                     'hermes.sivers.pi-.Qint.dx','hermes.sivers.k-.Qint.dx',
#                     'hermes.sivers.pi+.Qint.dz','hermes.sivers.k+.Qint.dz',
#                     'hermes.sivers.pi-.Qint.dz','hermes.sivers.k-.Qint.dz',
#                     'hermes.sivers.pi+.Q<2.dpt','hermes.sivers.k+.Q<2.dpt',
#                     'hermes.sivers.pi+.Q>2.dpt','hermes.sivers.k+.Q>2.dpt',                    
#                     'hermes.sivers.pi+.Q<2.dz','hermes.sivers.k+.Q<2.dz',
#                     'hermes.sivers.pi+.Q>2.dz','hermes.sivers.k+.Q>2.dz',                    
#                     'hermes.sivers.pi+.3d','hermes.sivers.pi-.3d']))

# theData=DataProcessor.DataMultiSet.DataMultiSet("SIDISset",loadThisData([
# 'compass16.sivers.h+.1<z<2.dpt',
# 'compass16.sivers.h-.1<z<2.dpt',
# 'compass16.sivers.h+.1<z<2.dx',
# 'compass16.sivers.h-.1<z<2.dx',
# 'compass16.sivers.h+.1<z<2.dz',
# 'compass16.sivers.h-.1<z<2.dz',
# 'compass16.sivers.h+.z>1.dpt' ,
# 'compass16.sivers.h-.z>1.dpt' ,
# 'compass16.sivers.h+.z>1.dx'  ,
# 'compass16.sivers.h-.z>1.dx'  ,
# 'compass16.sivers.h+.z>1.dz'  ,
# 'compass16.sivers.h-.z>1.dz'  ,
# 'compass16.sivers.h+.z>2.dpt' ,
# 'compass16.sivers.h-.z>2.dpt' ,
# 'compass16.sivers.h+.z>2.dx'  ,
# 'compass16.sivers.h-.z>2.dx'  ,
# 'compass16.sivers.h+.z>2.dz'  ,
# 'compass16.sivers.h-.z>2.dz']))

# theData=DataProcessor.DataMultiSet.DataMultiSet("SIDISset",loadThisData([
#     "jlab.sivers.pi+","jlab.sivers.pi-","jlab.sivers.k+","jlab.sivers.k-"]))
    

theData=DataProcessor.DataMultiSet.DataMultiSet("SIDISset",loadThisData([
                    'compass.sivers.pi+.dpt', 'compass.sivers.pi-.dpt',
                    'compass.sivers.k+.dpt', 'compass.sivers.k-.dpt',
                    'compass16.sivers.h+.1<z<2.dpt','compass16.sivers.h-.1<z<2.dpt',
                    'compass16.sivers.h+.z>2.dpt' ,'compass16.sivers.h-.z>2.dpt',
                    'hermes.sivers.pi+.3d','hermes.sivers.pi-.3d',
                    'jlab.sivers.pi+','jlab.sivers.pi-','jlab.sivers.k+']))

setSIDIS=theData.CutData(cutFunc) 

theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisData([
                    'star.sivers.W+.dqT','star.sivers.W-.dqT',
                    #'star.sivers.W+.dy','star.sivers.W-.dy',
                    'star.sivers.Z',
                    'compass.sivers.piDY.dqT'
                    #'compass.sivers.piDY.dQ','compass.sivers.piDY.dxF'
                    ]))

setDY=theData.CutData(cutFunc) 

print('Loaded (SIDIS)', setSIDIS.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setSIDIS.sets]), 'points.')
print('Loaded SIDIS experiments are', [i.name for i in setSIDIS.sets])

print('Loaded (DY)', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
print('Loaded DY experiments are', [i.name for i in setDY.sets])

#%%
harpy.setNPparameters_TMDR([1.92819, 0.0390534])
harpy.setNPparameters_uTMDPDF([0.198279, 9.29836, 431.647, 2.11829, -4.44162, 0., 0.17, 0.48, 2.15])
harpy.setNPparameters_uTMDFF([0.259499, 0.476235, 0.477143, 0.482977])
#harpy.setNPparameters_SiversTMDPDF([0.23, 0., 0.5, 7, -0.1, -0.2, 6, -0.1, -0.03, 8, -0.2])
harpy.setNPparameters_SiversTMDPDF([4.8, 0.0, 0.4, 18, -1.2, 1, 5.1, 0.5, -0.15, 1.09, -0.96])

#%%
DataProcessor.harpyInterface.PrintChi2Table(setSIDIS,method="central",printSysShift=False)

DataProcessor.harpyInterface.PrintChi2Table(setDY)

    
#%%
#######################################
# Minimisation
#######################################
totalN=setSIDIS.numberOfPoints+setDY.numberOfPoints

def chi_2(x):
    startT=time.time()
    harpy.setNPparameters_SiversTMDPDF(x)
    print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
    
    ccSIDIS2,cc3=DataProcessor.harpyInterface.ComputeChi2(setSIDIS,method="central")
    ccDY2,cc3=DataProcessor.harpyInterface.ComputeChi2(setDY)
    
    cc=(ccSIDIS2+ccDY2)/totalN
    endT=time.time()
    print(':->',cc,'       t=',endT-startT)
    return ccSIDIS2+5*ccDY2

#%%

from iminuit import Minuit


initialValues=(4.8, 0.0, 0.4, 18, -1.2, 1, 5.1, 0.5, -0.15, 1.09, -0.96)

initialErrors=(0.1, 0.001, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,0.1, 0.1, 0.1)
searchLimits=((0,30),(0,200), (-25,25),(0.01,20),(-50,50), (-25,25),(0.01,20),(-50,50), (-25,25),(0.01,20),(-50,50))
parametersToMinimize=(False,False, False, False, False, False, False, False,False, False, False)

m = Minuit.from_array_func(chi_2, initialValues,
      error=initialErrors, limit=searchLimits, fix=parametersToMinimize, errordef=1)

#m.get_param_states()

m.tol=0.0001*totalN*10000 ### the last 0.0001 is to compensate MINUIT def
m.strategy=1

SaveToLog("MINIMIZATION STARTED",str(m.params))
#%%


m.tol=0.0001*totalN*10000 ### the last 0.0001 is to compensate MINUIT def
m.strategy=1

m.migrad()

print(m.params)

SaveToLog("MINIMIZATION FINISHED",str(m.params))
SaveToLog("CORRELATION MATRIX",str(m.matrix(correlation=True)))
#%%

m.hesse()

print(m.params)

SaveToLog("HESSE FINISHED",str(m.params))
SaveToLog("CORRELATION MATRIX",str(m.matrix(correlation=True)))
#%%

m.minos()

print(m.params)

SaveToLog("MINOS FINISHED",str(m.params))
SaveToLog("CORRELATION MATRIX",str(m.matrix(correlation=True)))

#%%

print("{")
for j in range(len(setSIDIS.sets)):
    s=setSIDIS.sets[j]
    YY=DataProcessor.harpyInterface.ComputeXSec(s,method="central")
    print('{"'+s.name+'",{')
    for i in range(s.numberOfPoints):
        print("{"+"{:2.4f},{:2.4f},{:2.4f},{:2.4f},{:2.4f},{:2.4f},{:12.6f},{:12.6f},{:12.6f}".format(
            s.points[i]["x"][0],s.points[i]["x"][1],
            s.points[i]["z"][0],s.points[i]["z"][1],
            s.points[i]["pT"][0],s.points[i]["pT"][1],
            s.points[i]["xSec"],numpy.sqrt(numpy.sum(numpy.array(s.points[i]["uncorrErr"])**2)),
            YY[i]
            ),end="")
        if i==s.numberOfPoints-1:
            if j==len(setSIDIS.sets)-1:
                print("}}}}")
            else:
                print("}}},")
        else:
            print("},")
#%%
print("{")
for j in range(len(setDY.sets)):
    s=setDY.sets[j]
    YY=DataProcessor.harpyInterface.ComputeXSec(s)
    print('{"'+s.name+'",{')
    for i in range(s.numberOfPoints):
        print("{"+"{:2.4f},{:2.4f},{:2.4f},{:2.4f},{:12.6f},{:12.6f},{:12.6f}".format(
            s.points[i]["y"][0],s.points[i]["y"][1],
            s.points[i]["qT"][0],s.points[i]["qT"][1],
            s.points[i]["xSec"],numpy.sqrt(numpy.sum(numpy.array(s.points[i]["uncorrErr"])**2)),
            YY[i]
            ),end="")
        if i==s.numberOfPoints-1:
            if j==len(setSIDIS.sets)-1:
                print("}}}}")
            else:
                print("}}},")
        else:
            print("},")
