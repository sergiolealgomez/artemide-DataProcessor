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
import DataProcessor.ArtemideReplicaSet

#MAINPATH="/home/m/Github/artemide-DataProcessor/"
MAINPATH="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/"

useOrder="nnlo"
#useOrder="n3lo"

#%%
#######################################
#Initialize artemide
#######################################
import harpy
path_to_constants=MAINPATH+"FittingPrograms/Sivers20/Constants-files/"
#harpy.initialize(path_to_constants+"const-Sivers20_lo")

# harpy.initialize(path_to_constants+"const-Sivers20_nnlo_piK")

if(useOrder=="nnlo"):
    harpy.initialize(path_to_constants+"const-Sivers20_nnlo")
    
    #### All=0 Case
    harpy.setNPparameters_TMDR([2., 0.0398333])
    harpy.setNPparameters_uTMDPDF([0.185239, 6.22706, 580.946, 2.44166, -2.53161, 0.,  0.0014, 0.442, 4.14])
    harpy.setNPparameters_uTMDFF([0.279443, 0.460015, 0.435955, 0.551302])
    
    # # #### All=0 case piK
    # harpy.setNPparameters_TMDR([2., 0.0394095])
    # harpy.setNPparameters_uTMDPDF([0.180718, 4.38119, 426.208, 2.22347, -0.0646396, 0., 0.17, 0.48, 2.15])
    # harpy.setNPparameters_uTMDFF([0.293548, 0.462093, 0.442867, 0.590596, 0.427915, 0.462578, 0.304421,1.18113])
elif(useOrder=="n3lo"):
    harpy.initialize(path_to_constants+"const-Sivers20_n3lo")
    #### All=0 Case n3lo
    harpy.setNPparameters_TMDR([2., 0.0442327])
    harpy.setNPparameters_uTMDPDF([0.17975, 3.9081, 453.883, 2.07199, 1.60774, 0,  0.0014, 0.442, 4.14])
    harpy.setNPparameters_uTMDFF([0.26994, 0.456091, 0.423312, 0.615092])
    

harpy.setNPparameters_SiversTMDPDF([5.2, 0.,0.,0.,0., -0.6, 15.9, 0.5, -0.2, 21.6, -0.5, -0.1, 0.4, -1.1]) 
#%%
### read the list of files and return the list of DataSets
def loadThisData(listOfNames):    
    import DataProcessor.DataSet
    
    #path_to_data="/home/m/Github/artemide-DataProcessor/DataLib/Sivers/"
    path_to_data="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/DataLib/Sivers/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   

    return dataCollection

#%%
##################Cut function
def cutFunc(p):
    import copy
    
    if p["type"]=="DY":
        deltaTEST=0.3
        delta=p["<qT>"]/p["<Q>"]        
        
        
        if(9<p["<Q>"]<11):#UPSILON resonance-bin
            return False , p
    
    if p["type"]=="SIDIS":   
        deltaTEST=0.3
        delta=p["<pT>"]/p["<z>"]/p["<Q>"]        
    
    
    if delta<deltaTEST:
        pNew=copy.deepcopy(p)    
        pNew["process"]=pNew["weightProcess"]
        if p["type"]=="SIDIS":
            normX=DataProcessor.harpyInterface.ComputeXSec(pNew,method="central")        
        elif p["type"]=="DY":
            normX=DataProcessor.harpyInterface.ComputeXSec(pNew)        
        else:
            print("Are you crazy?")
        p["thFactor"]=p["thFactor"]/normX        
    
    #### This is because star measures AN
    if p["id"][0:4]=="star":
        p["thFactor"]=-p["thFactor"]        
        
#    return delta<0.5 and p.qT_avarage<80
    return delta<deltaTEST, p

#%%
### Loading the data set
# theData=DataProcessor.DataMultiSet.DataMultiSet("SIDISset",loadThisData([
#                     'compass.sivers.pi+.dpt', 'compass.sivers.pi-.dpt',
#                     'compass.sivers.k+.dpt', 'compass.sivers.k-.dpt',                  
#                     'hermes.sivers.pi+.Qint.dpt','hermes.sivers.k+.Qint.dpt',
#                     'hermes.sivers.pi-.Qint.dpt','hermes.sivers.k-.Qint.dpt',
#                     'hermes.sivers.pi+.Q<2.dpt','hermes.sivers.k+.Q<2.dpt',
#                     'hermes.sivers.pi+.Q>2.dpt','hermes.sivers.k+.Q>2.dpt',                  
#                     'hermes.sivers.pi+.3d','hermes.sivers.pi-.3d',
#                     'hermes.sivers.k+.3d','hermes.sivers.k-.3d',
#                     'compass16.sivers.h+.1<z<2.dpt',
#                     'compass16.sivers.h-.1<z<2.dpt',
#                     'compass16.sivers.h+.z>1.dpt' ,
#                     'compass16.sivers.h-.z>1.dpt' ,
#                     'compass16.sivers.h+.z>2.dpt' ,
#                     'compass16.sivers.h-.z>2.dpt' ,
#                     "jlab.sivers.pi+","jlab.sivers.pi-",
#                     "jlab.sivers.k+"]))
    

theData=DataProcessor.DataMultiSet.DataMultiSet("SIDISset",loadThisData([
                    'compass.sivers.pi+.dpt', 'compass.sivers.pi-.dpt',
                    'compass.sivers.k+.dpt', 'compass.sivers.k-.dpt',
                    'compass16.sivers.h+.1<z<2.dpt','compass16.sivers.h-.1<z<2.dpt',
                    'compass16.sivers.h+.z>2.dpt' ,'compass16.sivers.h-.z>2.dpt',
                    'hermes.sivers.pi+.3d','hermes.sivers.pi-.3d',
                    'hermes.sivers.k+.3d','hermes.sivers.k-.3d',
                    'jlab.sivers.pi+','jlab.sivers.pi-','jlab.sivers.k+'
                    ]))

setSIDIS=theData.CutData(cutFunc) 

theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisData([
                    'star.sivers.W+.dqT','star.sivers.W-.dqT',
                    #'star.sivers.W+.dy','star.sivers.W-.dy',
                    'star.sivers.Z',
                    'compass.sivers.piDY.dqT'
                    #,'compass.sivers.piDY.dQ','compass.sivers.piDY.dxF'
                    ]))

setDY=theData.CutData(cutFunc) 

print('Loaded (SIDIS)', setSIDIS.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setSIDIS.sets]), 'points.')
print('Loaded SIDIS experiments are', [i.name for i in setSIDIS.sets])

print('Loaded (DY)', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
print('Loaded DY experiments are', [i.name for i in setDY.sets])

print('Total number of points:',setSIDIS.numberOfPoints+setDY.numberOfPoints)

#%%
rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile("/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/Sivers20/REPS/"+
                                                  "Sivers20_model9case1.rep")
rSet.SetReplica()

DataProcessor.harpyInterface.PrintChi2Table(setSIDIS,method="central",printSysShift=False)

DataProcessor.harpyInterface.PrintChi2Table(setDY)

#%%
###########################
### Computation of chi^2 for each replica
##########################
cccSIDIS=[]
cccDY=[]
cccZ=[]
cccTotal=[]
for i in range(1,rSet.numberOfReplicas+1):
    rSet.SetReplica(i)
    
    ccSIDIS2,cc3=DataProcessor.harpyInterface.ComputeChi2(setSIDIS,method="central")
    ccDY2,cc3=DataProcessor.harpyInterface.ComputeChi2(setDY)
    cccSIDIS.append(ccSIDIS2)
    cccDY.append(ccDY2)
    cccZ.append(cc3[0]+cc3[1]+cc3[3])
    cccTotal.append(ccSIDIS2+ccDY2)

#%%
### SIDIS nnlo
#chiSIDIS=0.865
#chiDY=1.248
#chiZ=2.873

### SIDIS+DY nnlo
#chiSIDIS=0.880
#chiDY=0.782
#chiZ=1.624

### SIDIS n3lo
#chiSIDIS=0.848
#chiDY=1.236
#chiZ=2.841

### SIDIS+DY n3lo
# chiSIDIS=0.883
# chiDY=0.792
# chiZ=1.559

# rr=ComputeParameters(cccSIDIS)
# print("SIDIS: ",chiSIDIS,rr[2]/63-chiSIDIS,rr[3]/63-chiSIDIS)

# rr=ComputeParameters(cccDY)
# print("DY: ",chiDY,rr[2]/12-chiDY,rr[3]/12-chiDY)

# chii=(chiDY*12-chiZ)/11
# rr=ComputeParameters(cccZ)
# print("DY/Z: ",chii,rr[2]/11-chii,rr[3]/11-chii)

# chii=(chiDY*12+chiSIDIS*63)/75
# rr=ComputeParameters(cccTotal)
# print("Total: ",chii,rr[2]/75-chii,rr[3]/75-chii)


#%%
# rSet.SetReplica(0)
# for x in [0.01,0.02, 0.04, 0.06,0.08,0.1,0.2,0.4,0.6,0.8,0.9]:
#     for kT in [0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0]:
#         tmd=harpy.get_uTMDPDF_kT(x, kT, 1,mu=100.)
#         #tmd=harpy.get_SiversTMDPDF_kT(x, kT, 1, mu=2.)
#         print("{"+"{:12.6f},{:12.6f},{:12.6f}".format(x,kT,tmd[1+5]) +"},")
        
#%%
# for x in [0.01,0.02, 0.04, 0.06,0.08,0.1,0.2,0.4,0.6,0.8,0.9]:
#     for kT in [0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0]:
#         tmd=getSiversRowP(x,kT,1,mu=100.)
#         pp=ComputeParameters(tmd)
#         pp2=numpy.max([numpy.abs(pp[2]),numpy.abs(pp[3])])
#         print("{"+"{:12.6f},{:12.6f},{:12.6f}".format(x,kT,pp2) +"},")

#%%
#######################################
# Minimisation
#######################################
includeDY=True
if includeDY:
    totalN=setSIDIS.numberOfPoints+setDY.numberOfPoints
else :
    totalN=setSIDIS.numberOfPoints

def chi_2(x):
    startT=time.time()
    harpy.setNPparameters_SiversTMDPDF(x)
    print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
    
    ccSIDIS2,cc3=DataProcessor.harpyInterface.ComputeChi2(setSIDIS,method="central")
    
    if includeDY:
        ccDY2,cc3=DataProcessor.harpyInterface.ComputeChi2(setDY)
    else:
        ccDY2,cc3=0,0
    
    cc=(ccSIDIS2+ccDY2)/totalN
    endT=time.time()
    print(':->',cc,'       t=',endT-startT)
    return ccSIDIS2+ccDY2
#%%
from iminuit import Minuit

if(useOrder=="nnlo"):
    #initialValues=(0.02617, 1.65313, 29.5239, 0.0, 0.0, -0.04842, -0.39636, -2.61426, 0.18038, -0.45053, 14.0649, 0.45734, 2.5715, -0.25895)
    #initialValues=(0.131196, 3.29739, 64.8559, 0., 0., -0.0350288, -0.360633, -3.48127, 0.278947, -0.543139, 11.7035, 0.749227, 2.6025, -0.443391)
    initialValues=(0.15856, 3.45784, 69.1464, 0.0, 0.0, -0.0352, -0.35094, -3.49901, 0.27941, -0.58113, 7.50691, 0.70375, 2.59306, -0.42195)
elif(useOrder=="n3lo"):
    initialValues=(0.04915, 3.285, 51.92, 0.0, 0.0, -0.03589, -0.3689, -3.5055, 0.2822, -0.59204, 9.572, 0.8442, 2.744, -0.4627)

initialErrors=(0.1, 10., 10. , 0.1,0.1,
               0.3, 0.5, 1., 
               1., 1., 1., 
               1., 1., 0.1)
# searchLimits=((0.0001,2.),(0.01,80), (0.0,400), None,None,
#               (-2.,2.),(-0.99,2.5), (-10.,10.),
#               (-20.,20.),(-0.99,6.), (-15.,15.),
#               (-5.,5.),(-0.99,8), None)
searchLimits=((0,None),(0,None), (0,None), None, None,
              (-5.,5.),(-0.99,30), (-30.,30.),
              (-5.,5.),(-0.99,60.),(-30.,30.),
              (-5.,15.),(-0.99,30.), (-5,5))

parametersToMinimize=(False,False,False,True,True, 
                      False, False, False, 
                      False, False, False,
                      False, False, False)


m = Minuit.from_array_func(chi_2, initialValues,
      error=initialErrors, limit=searchLimits, fix=parametersToMinimize, errordef=1)

#m.get_param_states()

m.tol=0.0001*totalN*10000 ### the last 0.0001 is to compensate MINUIT def
m.strategy=1

#%%


# m.tol=0.0001*totalN*10000 ### the last 0.0001 is to compensate MINUIT def
# m.strategy=1

m.migrad(ncall=100)

print(m.params)
# sys.exit()
# SaveToLog("MINIMIZATION FINISHED",str(m.params))
# SaveToLog("CORRELATION MATRIX",str(m.matrix(correlation=True)))
#%%

# m.hesse()

# print(m.params)

# print(m.matrix(correlation=True))

# SaveToLog("HESSE FINISHED",str(m.params))
# SaveToLog("CORRELATION MATRIX",str(m.matrix(correlation=True)))
#%%

# m.minos()

# print(m.params)

# SaveToLog("MINOS FINISHED",str(m.params))
# SaveToLog("CORRELATION MATRIX",str(m.matrix(correlation=True)))

 
#%%
# # # #### JOINED PLOT without bins over replicas
# print("{")
# #rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile("/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/Sivers20/REPS/Sivers20_model9case1(noDY).rep")

# for j in range(len(setSIDIS.sets)):
#     s=setSIDIS.sets[j]
#     YYlist=[]
#     for r in range(rSet.numberOfReplicas):
#         rSet.SetReplica(r)
#         YY0=DataProcessor.harpyInterface.ComputeXSec(s,method="central")
#         YYlist.append(YY0)
        
#     YY=numpy.mean(YYlist,axis=0)
#     YYstd=numpy.std(YYlist,axis=0)
                
#     print('{"'+s.name+'",{')
#     for i in range(s.numberOfPoints):
#         print("{"+"{:2.4f},{:2.4f},{:2.4f},{:12.6f},{:12.6f},{:12.6f},{:12.6f}".format(
#             s.points[i]["<x>"],
#             s.points[i]["<z>"],
#             s.points[i]["<pT>"],
#             s.points[i]["xSec"],numpy.sqrt(numpy.sum(numpy.array(s.points[i]["uncorrErr"])**2)),
#             YY[i],YYstd[i]
#             ),end="")
#         if i==s.numberOfPoints-1:
#             print("}}},")                
#         else:
#             print("},")
# for j in range(len(setDY.sets)):
#     s=setDY.sets[j]
#     YYlist=[]
#     for r in range(rSet.numberOfReplicas):
#         rSet.SetReplica(r)
#         YY0=DataProcessor.harpyInterface.ComputeXSec(s)
#         YYlist.append(YY0)
        
#     YY=numpy.mean(YYlist,axis=0)
#     YYstd=numpy.std(YYlist,axis=0)
        
#     print('{"'+s.name+'",{')
#     for i in range(s.numberOfPoints):
#         print("{"+"{:4d},{:4d},{:2.4f},{:12.6f},{:12.6f},{:12.6f},{:12.6f}".format(
#             -1,
#             -1,
#             s.points[i]["<qT>"],
#             s.points[i]["xSec"],numpy.sqrt(numpy.sum(numpy.array(s.points[i]["uncorrErr"])**2)),
#             YY[i],YYstd[i]
#             ),end="")
#         if i==s.numberOfPoints-1:
#             if j==len(setDY.sets)-1:
#                 print("}}}}")
#             else:
#                 print("}}},")                
#         else:
#             print("},")

#%%

def MinForReplica():
    
    
    def repchi_2(x):        
        startT=time.time()
        harpy.setNPparameters_SiversTMDPDF(x)
        print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
        
        
        if(includeDY):
            ccDY2,cc3=DataProcessor.harpyInterface.ComputeChi2(repDataDY)
        else:
            ccDY2,cc3=0,0
            
        ccSIDIS2,cc3=DataProcessor.harpyInterface.ComputeChi2(repDataSIDIS,method="central")
        cc=(ccDY2+ccSIDIS2)/totalNnew
        
        endT=time.time()
        print(':->',cc,'       t=',endT-startT)
        return ccSIDIS2+ccDY2
    
    repDataSIDIS=setSIDIS.GenerateReplica()
    if(includeDY):
        repDataDY=setDY.GenerateReplica()
        totalNnew=repDataSIDIS.numberOfPoints+repDataDY.numberOfPoints
    else:
        totalNnew=repDataSIDIS.numberOfPoints    
    
    localM = Minuit.from_array_func(repchi_2, initialValues,
      error=initialErrors, limit=searchLimits, fix=parametersToMinimize, errordef=1)
    
    localM.tol=0.0001*totalNnew*10000 ### the last 0.0001 is to compensate MINUIT def
    localM.strategy=1

    localM.migrad()
    
    chi2Central=chi_2(localM.values.values())
    
    return [localM.fval,chi2Central,localM.values.values()]

#%%
#
# Generate pseudo data and minimise   100 times
#
numOfReplicas=450
REPPATH=MAINPATH+"FittingPrograms/Sivers20/LOGS/"+"model9case1(noDY-n3lo)EXTRA-replicas.txt"
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