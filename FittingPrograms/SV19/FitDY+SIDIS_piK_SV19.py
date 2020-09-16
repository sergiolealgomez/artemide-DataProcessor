#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 16:40:59 2019

SIDIS only fit with separation of pi and kaon TMDFF

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

MAINPATH="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/"
logFile=MAINPATH+"FittingPrograms/LOGS/"+"SV19[all=0_piK].log"

#%%
#######################################
# LOG save function
#######################################
savedTime=time.time()
def SaveToLog(text):
    global savedTime,logFile
    newTime=time.time()
    
    import socket
    PCname=socket.gethostname()
    
    passedTime=newTime-savedTime
    hours=int(passedTime/3600)
    minutes=int((passedTime-hours*3600)/60)
    seconds=int(passedTime-hours*3600-minutes*60)
    
    with open(logFile, 'a') as file:
        file.write(PCname+ ' : ' + time.ctime()+' :  [+'+str(hours)+':'+str(minutes)+':'+str(seconds)+' ]\n')
        file.write(' --> '+text+'\n')
        file.write('\n')
    savedTime=time.time()
#%%
#######################################
#Initialize artemide
#######################################
import harpy
path_to_constants=MAINPATH+"FittingPrograms/SV19/Constants-files/"
harpy.initialize(path_to_constants+"DY+SIDIS_nnlo_piK_all=0/const-DY+SIDIS_NNPDF31+DSS_nnlo_piK_all=0")
harpy.setNPparameters_TMDR([1.93, 0.0434])
harpy.setNPparameters_uTMDPDF([0.253434, 9.04351, 346.999, 2.47992, -5.69988, 0.1, 0.])
harpy.setNPparameters_uTMDFF([0.264,0.479,0.459,0.539,0.264,0.479,0.459,0.539]) 

#%%
### read the list of files and return the list of DataSets
def loadThisDataDY(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/DataLib/unpolDY/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   

    return dataCollection

def loadThisDataSIDIS(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/DataLib/unpolSIDIS/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   

    return dataCollection

#%%
##################Cut function
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
### Loading the DY data set
theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY([
                          'CDF1', 'CDF2', 'D01', 'D02', 'D02m', 
                          'A7-00y10', 'A7-10y20','A7-20y24', 
                          'A8-00y04', 'A8-04y08', 'A8-08y12', 'A8-12y16', 'A8-16y20', 'A8-20y24', 
                          'A8-46Q66', 'A8-116Q150', 
                          'CMS7', 'CMS8', 
                          'LHCb7', 'LHCb8', 'LHCb13', 
                          'PHE200', 'E228-200', 'E228-300', 'E228-400', 
                          'E772',
                          'E605']))

setDY=theData.CutData(cutFunc) 

### Loading the SIDIS data set
theData=DataProcessor.DataMultiSet.DataMultiSet("SIDISset",loadThisDataSIDIS([
                      'hermes.p.vmsub.zxpt.pi+','hermes.p.vmsub.zxpt.pi-',
                      'hermes.d.vmsub.zxpt.pi+','hermes.d.vmsub.zxpt.pi-',
                      'hermes.p.vmsub.zxpt.k+','hermes.p.vmsub.zxpt.k-',
                      'hermes.d.vmsub.zxpt.k+','hermes.d.vmsub.zxpt.k-',
                      'compass.d.h+','compass.d.h-']))

setSIDIS=theData.CutData(cutFunc) 

print('Loaded ', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setSIDIS.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setSIDIS.sets]), 'points.')
print('Loaded experiments are', [i.name for i in setSIDIS.sets])

#%%
for s in setDY.sets:
    if s.isNormalized:
        #s.normalizationMethod="bestChi2"
        s.normalizationMethod="integral"

#%%

#harpy.setNPparameters([1.92819, 0.0390534, 0.198279, 9.29836, 431.647, 2.11829, -4.44162, 0., 0., 0.259499, 0.476235, 0.477143, 0.482977]) ##NNPDF+DSS (paper)
#harpy.setNPparameters([1.92516, 0.0426578, 0.223809, 9.23868, 375.888, 2.14611, -4.97177, 0., 0., 0.233382, 0.478562, 0.47218, 0.511187]) ##NNPDF+DSS n3lo (paper)

##NNPDF+DSS (paper) pi/K
harpy.setNPparameters([2., 0.0397287, 0.180468, 4.33024, 423.481, 2.22637, -0.0862019, 0., 0.,
                       0.291174, 0.461295, 0.444895, 0.581916, 
                       0.428639, 0.457899, 0.295924, 1.22788]) 

# harpy.setNPparameters([2., 0.0396753, 0.185239, 6.22706, 580.946, 2.44166, -2.53161, 0., 0., 
#                        0.279443, 0.460015, 0.435955, 0.551302, 
#                        0.279443, 0.460015, 0.435955, 0.551302 ]) 

DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)
#DataProcessor.harpyInterface.PrintChi2Table(setSIDIS,printDecomposedChi2=True)
    
#%%
#######################################
# Minimisation
#######################################
totalN=setDY.numberOfPoints+setSIDIS.numberOfPoints

def chi_2(x):
    startT=time.time()
    harpy.setNPparameters(x)
    print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
    
    ccDY2,cc3=DataProcessor.harpyInterface.ComputeChi2(setDY)
    ccSIDIS2,cc3=DataProcessor.harpyInterface.ComputeChi2(setSIDIS)
    
    cc=(ccDY2+ccSIDIS2)/totalN
    endT=time.time()
    print(':->',cc,'       t=',endT-startT)
    return ccDY2+ccSIDIS2

#%%

from iminuit import Minuit

initialValues=( 2.000, 0.044, 0.172, 3.964, 528.650, 2.294, -0.011, 0.000, 0.000,
               0.303, 0.467, 0.432, 0.609, 
               0.425, 0.475, 0.355, 0.916)#pi/K-separation

initialErrors=(0.1,  0.05,  0.01,  1.0,   100,  0.1,   0.5,    1., 1., 0.1,   0.1,  0.1,    0.1,0.1,   0.1,  0.1,    0.1)

searchLimits=((1.,5.),   (0.,4.),  (0.,10.), (0.,10.),(0., 1000.),(0.,10.),None,None,None,
              (0.,5.),(0.,5.),(0.,5.),None,(0.,5.),(0.,5.),(0.,5.),None)
# True= FIX
parametersToMinimize=(True, False, False,False,False,False,False,True,True, False, False, False, False, False, False, False, False)


m = Minuit.from_array_func(chi_2, initialValues,
      error=initialErrors, limit=searchLimits, fix=parametersToMinimize, errordef=1)

#m.get_param_states()

m.tol=0.0001*totalN*10000 ### the last 0.0001 is to compensate MINUIT def
m.strategy=1

# SaveToLog("MINIMIZATION STARTED",str(m.params))
#%%
####################################
# Search for minimum
###################################

# m.migrad()

# print(m.params)

# DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)
# DataProcessor.harpyInterface.PrintChi2Table(setSIDIS,printDecomposedChi2=True)

# sys.exit()

# SaveToLog("MINIMIZATION FINISHED",str(m.params))
# SaveToLog("CORRELATION MATRIX",str(m.matrix(correlation=True)))

#%%
####################################
# Hesse matrix
###################################

# m.hesse()

# print(m.params)

# SaveToLog("HESSE FINISHED",str(m.params))
# SaveToLog("CORRELATION MATRIX",str(m.matrix(correlation=True)))

# sys.exit()
#%%
####################################
# Minos error-band estimation (VERY long)
###################################

# m.minos()

# print(m.params)

# SaveToLog("MINOS FINISHED",str(m.params))
# SaveToLog("CORRELATION MATRIX",str(m.matrix(correlation=True)))

#%%
def MinForReplica():
    
    
    def repchi_2(x):        
        startT=time.time()
        harpy.setNPparameters(x)
        print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
        
        ccDY2,cc3=DataProcessor.harpyInterface.ComputeChi2(repDataDY)
        ccSIDIS2,cc3=DataProcessor.harpyInterface.ComputeChi2(repDataSIDIS)
        
        cc=(ccDY2+ccSIDIS2)/totalNnew
        endT=time.time()
        print(':->',cc,'       t=',endT-startT)
        return ccDY2+ccSIDIS2
    
    repDataDY=setDY.GenerateReplica()
    repDataSIDIS=setSIDIS.GenerateReplica()
    totalNnew=repDataDY.numberOfPoints+repDataSIDIS.numberOfPoints
    
    localM = Minuit.from_array_func(repchi_2, initialValues,
      error=initialErrors, limit=searchLimits, fix=parametersToMinimize, errordef=1)
    
    localM.tol=0.0001*totalNnew*10000 ### the last 0.0001 is to compensate MINUIT def
    localM.set_strategy(1)

    localM.migrad()
    
    return [localM.fval,localM.values.values()]

#%%
#
# Generate pseudo data and minimise   100 times
#
numOfReplicas=25
REPPATH=MAINPATH+"FittingPrograms/LOGS/"+"SV19_REPLICAS_all=0_piK.txt"
for i in range(numOfReplicas):
    print('---------------------------------------------------------------')
    print('------------REPLICA ',i,'/',numOfReplicas,'--------------------')
    print('---------------------------------------------------------------')
    SaveToLog("Start computation of replica "+str(i) +"/"+ str(numOfReplicas))
    repRes=MinForReplica()
    SaveToLog("Minimization for replica "+str(i) +"/"+ str(numOfReplicas)+" finished.")
    print(repRes)
    f=open(REPPATH,"a+")
    print('SAVING >>  ',f.name)
    f.write(str(repRes)+"\n")
    f.close()
    
#%%
SaveToLog("Computation finished correctly ("+int(numOfReplicas)+" replicas computed) +\n "
          +"----------------------------------------------------------------------------------")