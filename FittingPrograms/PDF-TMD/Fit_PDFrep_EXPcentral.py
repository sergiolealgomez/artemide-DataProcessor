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
splitUpsilon=True
## Use the reduced set of the data
useReducedSet=True


#### Starting and final replica (included)
StartReplica=1
FinalReplica=3

## automatic name generation for the run
runName="model2.2_"+PDFinUse+"_PDFrep_"
if (not useA7data): runName+="noA7_"
if (splitUpsilon): runName+="spltUPS_"
if (useNormalizedLHCdata): runName+="norm_"
if (useReducedSet): runName+="reducedSet_"
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
dropPoints=["A7-10y20.7", "A7-10y20.8", "A7-10y20.9", "A7-20y24.4", "A7-20y24.5", \
"A7-20y24.6", "A7-20y24.7", "A7-20y24.8", "A7-20y24.9", "A8-08y12.8", \
"A8-116Q150.1", "A8-116Q150.2", "A8-116Q150.3", "A8-116Q150.4", \
"A8-116Q150.5", "A8-116Q150.6", "A8-116Q150.7", "A8-116Q150.8", \
"A8-116Q150.9", "A8-12y16.8", "A8-16y20.6", "A8-16y20.7", \
"A8-16y20.8", "A8-20y24.5", "A8-20y24.6", "A8-20y24.7", "A8-20y24.8", \
"CDF1.10", "CDF1.11", "CDF1.12", "CDF1.13", "CDF1.14", "CDF1.15", \
"CDF1.16", "CDF1.17", "CDF1.18", "CDF1.19", "CDF1.20", "CDF1.21", \
"CDF1.22", "CDF1.23", "CDF1.24", "CDF1.25", "CDF1.26", "CDF1.27", \
"CDF1.28", "CDF1.29", "CDF1.3", "CDF1.30", "CDF1.31", "CDF1.32", \
"CDF1.4", "CDF1.5", "CDF1.6", "CDF1.7", "CDF1.8", "CDF1.9", \
"CDF2.15", "CDF2.16", "CDF2.17", "CDF2.18", "CDF2.19", "CDF2.20", \
"CDF2.21", "CDF2.22", "CDF2.23", "CDF2.24", "CDF2.25", "CDF2.26", \
"CDF2.27", "CDF2.28", "CDF2.29", "CDF2.30", "CDF2.31", "CDF2.32", \
"CDF2.33", "CDF2.34", "CDF2.35", "CDF2.36", "CDF2.37", "CDF2.38", \
"CDF2.39", "CDF2.40", "CDF2.41", "CDF2.42", "CDF2.43", "CDF2.44", \
"CMS7.1", "CMS7.2", "CMS7.3", "CMS7.4", "CMS7.5", "CMS7.6", "CMS7.7", \
"CMS8.0", "CMS8.1", "CMS8.2", "CMS8.3", "CMS8.4", "CMS8.5", "CMS8.6", \
"CMS8.7", "D01.10", "D01.11", "D01.12", "D01.13", "D01.14", "D01.15", \
"D01.2", "D01.3", "D01.4", "D01.5", "D01.6", "D01.7", "D01.8", \
"D01.9", "D02.1", "D02.2", "D02.3", "D02.4", "D02.5", "D02.6", \
"D02.7", "D02.8", "D02m.3", "D02m.4", "E228-200.8Q9.10", \
"E228-200.8Q9.5", "E228-200.8Q9.6", "E228-200.8Q9.9", \
"E228-300.11Q12.11", "E228-300.11Q12.3", "E228-300.11Q12.4", \
"E228-300.11Q12.5", "E228-300.11Q12.6", "E228-300.11Q12.7", \
"E228-300.11Q12.8", "E228-300.8Q9.10", "E228-400.13Q14.11", \
"E228-400.13Q14.12", "E228-400.13Q14.6", "E228-400.13Q14.7", \
"E228-400.13Q14.9", "E605.7Q8.7", "E605.7Q8.8", "E605.7Q8.9", \
"E772.12Q13.10", "E772.12Q13.6", "E772.12Q13.8", "E772.12Q13.9", \
"E772.13Q14.4", "E772.13Q14.5", "E772.14Q15.0", "E772.14Q15.3", \
"E772.14Q15.4", "E772.14Q15.5", "E772.14Q15.6", "LHCb13.1", \
"LHCb13.10", "LHCb13.2", "LHCb13.3", "LHCb13.4", "LHCb13.5", \
"LHCb13.6", "LHCb13.7", "LHCb13.8", "LHCb13.9", "LHCb7.10", \
"LHCb7.2", "LHCb7.3", "LHCb7.4", "LHCb7.7", "LHCb7.8", "LHCb7.9", \
"LHCb8.10", "LHCb8.8", "LHCb8.9", "PHE200.2",\
##### extra points by hands
"A8-116Q150.0","CDF1.0","CDF1.1","CDF1.2","D01.0", "D01.1", "D02.0","CMS7.0","LHCb13.0"]


def cutFunc(p):    
    
    #### check against the presence in the reduced set.
    if(useReducedSet and p["id"] in dropPoints):
        return False,p
    
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
    #initialValues=(2.000,  0.033, 0.230, 5.609, 0.252, 8.021, 0.252, 8.021, 570.223) #model 2.0 HERA
    initialValues=(2.000, 0.034, 0.145, 10.393, 0.333, 0.000, 0.366, 10.261, 677.434) #model 2.2 HERA
if(PDFinUse=="NNPDF31"):
    #initialValues=(2.000,  0.029, 0.345, 2.587, 0.152, 7.561, 232.544, 0. , 0.)      #model 1.0 NNPDF31
    #initialValues=(2.000,  0.029, 0.345, 2.587, 0.152, 7.561, 0.152, 7.561, 232.544) #model 2.0 NNPDF31
    initialValues=(2.000, 0.030, 0.188, 5.542, 0.200, 4.375, 0.486, 0.009, 232.793) #model 2.2 NNPDF31
if(PDFinUse=="CT18"):
    #initialValues=(2.000,  0.039, 0.161, 7.904, 0.212, 5.301, 700.642, 0. , 0.)      #model 1.0 CT18
    #initialValues=(2.000,  0.039, 0.161, 7.904, 0.212, 5.301, 0.212, 5.301, 700.642) #model 2.0 CT18
    initialValues=(2.000, 0.042, 0.094, 12.534, 0.293, 0.004, 0.003, 16.568, 819.267) #model 2.2 CT18
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

m.get_param_states()

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