#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: vla18041
"""

#%%
#######################################################################
# Global parameter of a run
#######################################################################

#PDFinUse="HERA20"
#PDFinUse="NNPDF31"
PDFinUse="CT18"

## if this trigger is ON the LHC data will be fit only by shape
useNormalizedLHCdata=False
## Include ATLAS 7TeV?
useA7data=False
## Split the low-energy experiment <Upsilon and >Upsilon
splitUpsilon=False

#### Starting and final replica (included)
StartReplica=1
FinalReplica=3

## automatic name generation for the run
runName="model2.2_"+PDFinUse+"_PDFrep_"
if (not useA7data): runName+="noA7_"
if (splitUpsilon): runName+="spltUPS_"
if (useNormalizedLHCdata): runName+="norm_"
if (runName[-1]=="_"): runName=runName[0:-1]

print(" RUN: "+runName)

#%%
#######################################
# Paths
#######################################
PathToHarpy="/home/vla18041/LinkData2/arTeMiDe_Repository/artemide-ForPDF/harpy"
PathToDataProcessor="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/"
PathToDataLibrary=PathToDataProcessor+"DataLib/unpolDY/"
PathToLog=PathToDataProcessor+"FittingPrograms/PDF-TMD/LOGS/"
PathToConstantsFile=PathToDataProcessor+"/FittingPrograms/PDF-TMD/Constants-files/const-"+PDFinUse+"_NNLO_7p"

import sys
sys.path.remove('/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/harpy')
sys.path.append(PathToHarpy)
sys.path.append(PathToDataProcessor)

#%%
#######################################
# Output paths
#######################################
import socket
PCname=socket.gethostname()

replicaFile=PathToLog+runName+".txt"
logFile=PathToLog+PCname+".log"

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
# LOG save function
#######################################
savedTime=time.time()
def SaveToLog(text):
    global savedTime,logFile
    newTime=time.time()
    
    passedTime=newTime-savedTime
    hours=int(passedTime/3600)
    minutes=int((passedTime-hours*3600)/60)
    seconds=int(passedTime-hours*3600-minutes*60)
    
    with open(logFile, 'a') as file:
        file.write(time.ctime()+' :  [+'+str(hours)+':'+str(minutes)+':'+str(seconds)+' ]\n')
        file.write(' run '+runName+'\n')
        file.write(' --> '+text+'\n')
        file.write('\n')
    savedTime=time.time()

#%%
#######################################
#Initialize artemide
#######################################
import harpy

SaveToLog('Initialization with : \n'+PathToConstantsFile)

harpy.initialize(PathToConstantsFile)
initializationArray=[2.0340, 0.0299, 0.2512, 7.7572, 0.2512, 7.7572, 0.2512, 7.7572, 10000.]
harpy.setNPparameters(initializationArray)

#%%
#######################################
# read the list of files and return the list of DataSets
#######################################
def loadThisData(listOfNames):    
    import DataProcessor.DataSet
    
    dataCollection=[]
    for name in listOfNames:
        if( name==''): continue
        loadedData=DataProcessor.DataSet.LoadCSV(PathToDataLibrary+name+".csv")
        dataCollection.append(loadedData)   

    return dataCollection

#%%
#######################################
# Data cut function
#######################################
def cutFunc(p):    
    par=1.0

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
        
    if(p["id"][-2:]=="<u" and p["<Q>"]>10.5):
        return False,p
    if(p["id"][-2:]==">u" and p["<Q>"]<10.5):
        return False,p
    
#    return delta<0.5 and p.qT_avarage<80
    return (delta<0.1 or (delta<0.25 and par/err*delta**2<1)) , p

#%%
#######################################
# Loading the data set
#######################################

setHE=loadThisData(['CDF1', 'CDF2', 'D01', 'D02', 'D02m', 
                      'A7-00y10' if useA7data else '',
                      'A7-10y20' if useA7data else '',
                      'A7-20y24' if useA7data else '', 
                      'A8-00y04-norm' if useNormalizedLHCdata else 'A8-00y04',
                      'A8-04y08-norm' if useNormalizedLHCdata else 'A8-04y08',
                      'A8-08y12-norm' if useNormalizedLHCdata else 'A8-08y12',
                      'A8-12y16-norm' if useNormalizedLHCdata else 'A8-12y16',
                      'A8-16y20-norm' if useNormalizedLHCdata else 'A8-16y20',
                      'A8-20y24-norm' if useNormalizedLHCdata else 'A8-20y24',
                      'A8-46Q66-norm' if useNormalizedLHCdata else 'A8-46Q66',
                      'A8-116Q150-norm' if useNormalizedLHCdata else 'A8-116Q150',
                      'CMS7', 'CMS8', 
                      'LHCb7', 'LHCb8', 'LHCb13'])

#### If I separate data above and below UPSILON, I create two copies of LE data with different names
#### the data to be split only  'E228-300', 'E228-400' and E605
if(splitUpsilon):
    setLE1=loadThisData(['PHE200', 'E228-200','E772'])
    setLE2=loadThisData(['E228-300', 'E228-400','E605'])
    setLE3=loadThisData(['E228-300', 'E228-400','E605'])
    for s in setLE2:
        s.name+="-blwUPS"
        for p in s.points:
            p["id"]+="<u"
    for s in setLE3:
        s.name+="-abvUPS"
        for p in s.points:
            p["id"]+=">u"
    setLE=[setLE1[0],setLE1[1],setLE2[0],setLE3[0],setLE2[1],setLE3[1],setLE1[2],setLE2[2],setLE3[2]]
else:
    setLE=loadThisData(['PHE200', 'E228-200', 'E228-300', 'E228-400','E772','E605'])

theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",setHE+setLE)

setDY=theData.CutData(cutFunc) 

print('Loaded ', setDY.numberOfSets, 'data sets with ', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
print('Loaded experiments are', [i.name for i in setDY.sets])

SaveToLog('Loaded '+ str(setDY.numberOfSets) + ' data sets with '+str(sum([i.numberOfPoints for i in setDY.sets])) + ' points. \n'
+'Loaded experiments are '+str([i.name for i in setDY.sets]))

#%%
if useNormalizedLHCdata:
    for s in setDY.sets:
        if s.name[0:4]=='LHCb':
            s.isNormalized=True
            
        if s.isNormalized:
            s.normalizationMethod='bestChi2'
#%%

if(PDFinUse=="HERA20"):
    #initialValues=(2.000,  0.033, 0.230, 5.609, 0.252, 8.021, 570.223, 0.000, 0.000) #model 1.0 HERA
    initialValues=(2.000,  0.033, 0.230, 5.609, 0.252, 8.021, 0.252, 8.021, 570.223) #model 2.0 HERA
if(PDFinUse=="CT18"):
    #initialValues=(2.000,  0.039, 0.161, 7.904, 0.212, 5.301, 700.642, 0. , 0.)      #model 1.0 CT18
    initialValues=(2.000,  0.039, 0.161, 7.904, 0.212, 5.301, 0.212, 5.301, 700.642) #model 2.0 CT18
if(PDFinUse=="NNPDF31"):
    #initialValues=(2.000,  0.029, 0.345, 2.587, 0.152, 7.561, 232.544, 0. , 0.)      #model 1.0 NNPDF31
    initialValues=(2.000,  0.029, 0.345, 2.587, 0.152, 7.561, 0.152, 7.561, 232.544) #model 2.0 NNPDF31

harpy.setNPparameters(list(initialValues))

DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)
#%%
#######################################
# Main chi2 formula
#######################################
totalN=setDY.numberOfPoints

def chi_2(x):
    startT=time.time()
    harpy.setNPparameters(x)
    print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
    
    ccDY2,cc3=DataProcessor.harpyInterface.ComputeChi2(setDY)
    
    cc=ccDY2/totalN
    endT=time.time()
    print(':->',cc,'       t=',endT-startT)
    return ccDY2

#%%
#######################################
# Setting-up Minuit
#######################################
from iminuit import Minuit

##### model 1.0
#initialErrors=(0.1,       0.002,      0.05,    0.2,  0.05,   0.2,    10.0,   0.1,   0.1)
#searchLimits=((1.4,4.5), (0.0001,5.0),(0.0,2.0),(0.,25.0),(0.0,2.0),(0.,25.0),(0.,1000),None, None)
#parametersToMinimize=(True,     False,    False,    False,    False,     False,  False, True, True)

##### model 2.0
initialErrors=(0.1,       0.002,      0.05,    0.2,  0.05,    0.2,  0.05,   0.2,    10.0)
searchLimits=((1.4,4.5), (0.0001,5.0),(0.0,10.0),(0.,50.0),(0.0,10.0),(0.,50.0),(0.0,10.0),(0.,50.0),(0.,2500))
parametersToMinimize=(True,     False,    False,    False,    False,     False,  False, False, False)

# ##### model 3.0
# initialErrors=(0.1,       0.002,      0.05,    0.2,  0.05,   0.2,    10.0,   0.1,   0.1)
# searchLimits=((1.4,4.5), (0.0001,5.0),(0.0,2.0),(0.,25.0),(0.0,2.0),(0.,25.0),(0.,2500),None, (-0.5, 10.))
# parametersToMinimize=(True,     False,    False,    False,    False,     False,  False, False, False)


#%%

m = Minuit.from_array_func(chi_2, initialValues,
      error=initialErrors, limit=searchLimits, fix=parametersToMinimize, errordef=1)

#m.get_param_states()

m.tol=0.0001*totalN*10000 ### the last 0.0001 is to compensate MINUIT def
m.strategy=1
m.migrad()

## print parameters
print(m.params)

## print chi^2 table
harpy.setNPparameters(m.values.values())
DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)

sys.exit()
#%%
#######################################
# Generate replica of data and compute chi2
#######################################
def MinForReplica():
    global totalN,setDY,initialValues,initialErrors,searchLimits,parametersToMinimize
        
    def repchi_2(x):        
        global totalN
        startT=time.time()
        harpy.setNPparameters(x)
        print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
        
        ccDY2,cc2=DataProcessor.harpyInterface.ComputeChi2(repDataDY)
        
        cc=(ccDY2)/totalN
        endT=time.time()
        print(':->',cc,'       t=',endT-startT)
        return ccDY2
    
    repDataDY=setDY
    
    localM = Minuit.from_array_func(repchi_2, initialValues,
      error=initialErrors, limit=searchLimits, fix=parametersToMinimize, errordef=1)
    
    localM.tol=0.0001*totalN*10000 ### the last 0.0001 is to compensate MINUIT def
    localM.strategy=1

    localM.migrad()
    
    ### [chi^2, NP-parameters]
    return [localM.fval,localM.values.values()]

#%%
#######################################
# This is the main cicle. 
# It generates replica of data take random PDF and minimize it
# Save to log.
#######################################
for i in range(StartReplica,FinalReplica+1):
    print('---------------------------------------------------------------')
    print('------------REPLICA ',i,' from [',StartReplica,' , ',FinalReplica,']------------------')
    print('---------------------------------------------------------------')
    
    ## reset PDF
    harpy.setNPparameters(initializationArray)    
    harpy.setPDFreplica(i)
    SaveToLog("Start computation of replica "+str(i) +"in ["+ str(StartReplica)+','+str(FinalReplica)+"]")
    
    ## got to pseudo-data and minimization
    repRes=MinForReplica()
    print(repRes)
    SaveToLog("Minimization for replica "+str(i) +"in ["+ str(StartReplica)+','+str(FinalReplica)+"] finished.")    
    
    ## compute the chi2 for true data
    mainDY, mainDY2 =DataProcessor.harpyInterface.ComputeChi2(setDY)    
    SaveToLog("Central chi^2 for "+str(i) +"in ["+ str(StartReplica)+','+str(FinalReplica)+" computed. \n Saving to log >> "+replicaFile)
    
    ## save to file
    f=open(replicaFile,"a+")
    print('SAVING >>  ',f.name)
    ### [total chi^2(cenral), total chi^2 (pseudo data), list of chi^2 for experiments(central), number of PDF, list of NP-parameters]
    f.write(str([mainDY,repRes[0],mainDY2,i,repRes[1]])+"\n")
    f.close()
    
#%%
#######################################
# Finalizing log
#######################################
print('Computation finished')
SaveToLog("Computation finished correctly (["+ str(StartReplica)+','+str(FinalReplica)+"] range of replicas computed) +\n "
          +"----------------------------------------------------------------------------------")