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

#harpy.setNPparameters_TMDR([1.93, 0.0434])
#harpy.setNPparameters_uTMDPDF([0.253434, 9.04351, 346.999, 2.47992, -5.69988, 0.1, 0.17, 0.48, 2.15])
#harpy.setNPparameters_uTMDFF([0.264,0.479,0.459,0.539])
#harpy.setNPparameters_SiversTMDPDF([5.2, 0., -0.6, 15.9, 0.5, -0.2, 21.6, -0.5, -0.1, 0.4, -1.1]) 
# #### M=0 Case
# harpy.setNPparameters_TMDR([2., 0.0402402])
# harpy.setNPparameters_uTMDPDF([0.187436, 7.53637, 513.438, 2.24779, -2.52609, 0.,  0.17, 0.48, 2.15])
# harpy.setNPparameters_uTMDFF([0.192616, 0.474368, 0.504594, 0.419496])
# harpy.setNPparameters_SiversTMDPDF([5.2, 0.,0.,0.,0., -0.6, 15.9, 0.5, -0.2, 21.6, -0.5, -0.1, 0.4, -1.1]) 
#### All=0 Case
harpy.setNPparameters_TMDR([2., 0.0398333])
harpy.setNPparameters_uTMDPDF([0.184739, 6.22437, 588.193, 2.44327, -2.51106, 0.,  0.17, 0.48, 2.15])
harpy.setNPparameters_uTMDFF([0.277974, 0.459238, 0.43427, 0.55001])
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
#################### LOG save function
LOGPATH=MAINPATH+"FittingPrograms/Sivers20/LOGS/"+"Sivers20["+time.ctime()+"].log"
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
        deltaTEST=0.3
        
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
                    'hermes.sivers.k+.3d','hermes.sivers.k-.3d'
                    ,'jlab.sivers.pi+','jlab.sivers.pi-','jlab.sivers.k+'
                    ]))

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

print('Total number of points:',setSIDIS.numberOfPoints+setDY.numberOfPoints)

#%%
#all=0 
harpy.setNPparameters_TMDR([2., 0.0398333])
harpy.setNPparameters_uTMDPDF([0.184739, 6.22437, 588.193, 2.44327, -2.51106, 0.,  0.17, 0.48, 2.15])
harpy.setNPparameters_uTMDFF([0.277974, 0.459238, 0.43427, 0.55001])
#harpy.setNPparameters_SiversTMDPDF([0.23, 0., 0.5, 7, -0.1, -0.2, 6, -0.1, -0.03, 8, -0.2])
#harpy.setNPparameters_SiversTMDPDF([4.8, 0.0,0.0,0.0,0.0, 0.4, 18, -1.2, 1, 5.1, 0.5, -0.15, 1.09, -0.96])
## cosh model
#harpy.setNPparameters_SiversTMDPDF([0.942, -2.234, 0.000, 0.000, 0.000, 1.022, 9.954, 0.130, -0.611, 5.921, 0.082, 0.000,0.010, -1.086])
## hermes s=0
#harpy.setNPparameters_SiversTMDPDF([3.26491, -11.3924, 0., 0., 0., 2.93554, 21.1682, -1.42077, -1.33249,11.7114, -0.0910189, 0., 0.01, -0.25])
## hermes 
#harpy.setNPparameters_SiversTMDPDF([1.44543, -4.45842, 0., 0., 0., 1.84186, 11.4789, -0.219184, -1.51701, 10.3272, -0.799934, 0.0341125, 11.8556, -0.836211])
## hermes+compass s=0
#harpy.setNPparameters_SiversTMDPDF([0.493437, -0.584008, 0., 0., 0., 2.20891, 6.37475, 0.342108, -4.76376, 9.88546, 0.889865, 0., 0.01, -0.25])
## hermes+compass 
#harpy.setNPparameters_SiversTMDPDF([0.222126, -0.364747, 0., 0., 0., 2.53923, 5.55035, 0.761213, -7.17387, 8.79212, 0.887841, 0.308428, 6.10136, 0.750209])
## hermes+compass alpha=1
#harpy.setNPparameters_SiversTMDPDF([-0.119032, 2.53877, 0., 0., 0., 0.206096, 1., 0.458724, -5.75387, 1., 2.30951, -0.00791191, 1., -0.17452])
## hermes+compass alpha=2
#harpy.setNPparameters_SiversTMDPDF([-0.144051, 2.59002, 0., 0., 0., 0.604819, 2., 0.647903, -4.464, 2., 1.68727, 0.0655617, 2., -0.042426])
## hermes+compass; cosh2; alpha=1
#harpy.setNPparameters_SiversTMDPDF([0., -0.269222, 17.6731, 0., 0., 0.895419, 1., 0.773234, -3.95582, 1., 1.19785, 0.0639664, 1., 0.090659])
## hermes+compass model 2
harpy.setNPparameters_SiversTMDPDF([-0.032751, 1.95094, 0., 0., 0., 1.92258, -1.80177, 1.07152, -5.39127, -0.21834, 1.53353, 0.143266, -2.89196, 0.314361])

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
    ccDY2,cc3=0,0#DataProcessor.harpyInterface.ComputeChi2(setDY)
    
    cc=(ccSIDIS2+ccDY2)/totalN
    endT=time.time()
    print(':->',cc,'       t=',endT-startT)
    return ccSIDIS2+ccDY2

#%%

from iminuit import Minuit


#initialValues=(7.861, 4.929, 0.000, 0.000, 0.000, 0.268, 0.367, -2.386, 0.974, 0.1517, -1.5301, -0.4305, 0.010, -1.0857)
initialValues=(-0.119032, 2.53877, 0.0 , 0., 0., 
               2.53923, 0., 0.761213, 
               -7.17387,0., 0.887841, 
               0.308428,0., 0.750209)

initialErrors=(0.1, 0.1, 0.1,0.1,0.1,0.1, 0.1, 0.1, 0.1, 0.1, 0.1,0.1, 0.1, 0.1)
searchLimits=(None,None, None, None,None,
             (-25,25),None,(-25,25), 
             (-25,25),None,(-25,25), 
             (-25,25),None,(-25,25))
parametersToMinimize=(False,False,True,True,True, 
                      False, False, False, 
                      False, False, False,
                      False, False, False)


m = Minuit.from_array_func(chi_2, initialValues,
      error=initialErrors, limit=searchLimits, fix=parametersToMinimize, errordef=1)

#m.get_param_states()

m.tol=0.0001*totalN*10000 ### the last 0.0001 is to compensate MINUIT def
m.strategy=1

# SaveToLog("MINIMIZATION STARTED",str(m.params))
#%%


# m.tol=0.0001*totalN*10000 ### the last 0.0001 is to compensate MINUIT def
# m.strategy=1

# m.migrad()

# print(m.params)
# sys.exit()
# SaveToLog("MINIMIZATION FINISHED",str(m.params))
# SaveToLog("CORRELATION MATRIX",str(m.matrix(correlation=True)))
#%%

# m.hesse()

# print(m.params)

# SaveToLog("HESSE FINISHED",str(m.params))
# SaveToLog("CORRELATION MATRIX",str(m.matrix(correlation=True)))
#%%

# m.minos()

# print(m.params)

# SaveToLog("MINOS FINISHED",str(m.params))
# SaveToLog("CORRELATION MATRIX",str(m.matrix(correlation=True)))

#%%
##### SIDIS PLOT with bins
# print("{")
# for j in range(len(setSIDIS.sets)):
#     s=setSIDIS.sets[j]
#     YY=DataProcessor.harpyInterface.ComputeXSec(s,method="central")
#     print('{"'+s.name+'",{')
#     for i in range(s.numberOfPoints):
#         print("{"+"{:2.4f},{:2.4f},{:2.4f},{:2.4f},{:2.4f},{:2.4f},{:12.6f},{:12.6f},{:12.6f}".format(
#             s.points[i]["x"][0],s.points[i]["x"][1],
#             s.points[i]["z"][0],s.points[i]["z"][1],
#             s.points[i]["pT"][0],s.points[i]["pT"][1],
#             s.points[i]["xSec"],numpy.sqrt(numpy.sum(numpy.array(s.points[i]["uncorrErr"])**2)),
#             YY[i]
#             ),end="")
#         if i==s.numberOfPoints-1:
#             if j==len(setSIDIS.sets)-1:
#                 print("}}}}")
#             else:
#                 print("}}},")
#         else:
#             print("},")
#%%
##### DY PLOT with bins
# print("{")
# for j in range(len(setDY.sets)):
#     s=setDY.sets[j]
#     YY=DataProcessor.harpyInterface.ComputeXSec(s)
#     print('{"'+s.name+'",{')
#     for i in range(s.numberOfPoints):
#         print("{"+"{:2.4f},{:2.4f},{:12.6f},{:12.6f},{:12.6f}".format(
#             s.points[i]["qT"][0],s.points[i]["qT"][1],
#             s.points[i]["xSec"],numpy.sqrt(numpy.sum(numpy.array(s.points[i]["uncorrErr"])**2)),
#             YY[i]
#             ),end="")
#         if i==s.numberOfPoints-1:
#             if j==len(setSIDIS.sets)-1:
#                 print("}}}}")
#             else:
#                 print("}}},")
#         else:
#             print("},")

#%%
##### JOINED PLOT without bins
# print("{")
# for j in range(len(setSIDIS.sets)):
#     s=setSIDIS.sets[j]
#     YY=DataProcessor.harpyInterface.ComputeXSec(s,method="central")
#     print('{"'+s.name+'",{')
#     for i in range(s.numberOfPoints):
#         print("{"+"{:2.4f},{:2.4f},{:2.4f},{:12.6f},{:12.6f},{:12.6f}".format(
#             s.points[i]["<x>"],
#             s.points[i]["<z>"],
#             s.points[i]["<pT>"],
#             s.points[i]["xSec"],numpy.sqrt(numpy.sum(numpy.array(s.points[i]["uncorrErr"])**2)),
#             YY[i]
#             ),end="")
#         if i==s.numberOfPoints-1:
#             print("}}},")                
#         else:
#             print("},")
# for j in range(len(setDY.sets)):
#     s=setDY.sets[j]
#     YY=DataProcessor.harpyInterface.ComputeXSec(s)
#     print('{"'+s.name+'",{')
#     for i in range(s.numberOfPoints):
#         print("{"+"{:4d},{:4d},{:2.4f},{:12.6f},{:12.6f},{:12.6f}".format(
#             -1,
#             -1,
#             s.points[i]["<qT>"],
#             s.points[i]["xSec"],numpy.sqrt(numpy.sum(numpy.array(s.points[i]["uncorrErr"])**2)),
#             YY[i]
#             ),end="")
#         if i==s.numberOfPoints-1:
#             if j==len(setDY.sets)-1:
#                 print("}}}}")
#             else:
#                 print("}}},")                
#         else:
#             print("},")
            
#%%
##### JOINED PLOT without bins over replicas
print("{")
import DataProcessor.ArtemideReplicaSet
rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile("/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/Sivers20/REPS/Sivers20_model2.rep")

for j in range(len(setSIDIS.sets)):
    s=setSIDIS.sets[j]
    YYlist=[]
    for r in range(rSet.numberOfReplicas):
        rSet.SetReplica(r)
        YY0=DataProcessor.harpyInterface.ComputeXSec(s,method="central")
        YYlist.append(YY0)
        
    YY=numpy.mean(YYlist,axis=0)
    YYstd=numpy.std(YYlist,axis=0)
                
    print('{"'+s.name+'",{')
    for i in range(s.numberOfPoints):
        print("{"+"{:2.4f},{:2.4f},{:2.4f},{:12.6f},{:12.6f},{:12.6f},{:12.6f}".format(
            s.points[i]["<x>"],
            s.points[i]["<z>"],
            s.points[i]["<pT>"],
            s.points[i]["xSec"],numpy.sqrt(numpy.sum(numpy.array(s.points[i]["uncorrErr"])**2)),
            YY[i],YYstd[i]
            ),end="")
        if i==s.numberOfPoints-1:
            print("}}},")                
        else:
            print("},")
for j in range(len(setDY.sets)):
    s=setDY.sets[j]
    YYlist=[]
    for r in range(rSet.numberOfReplicas):
        rSet.SetReplica(r)
        YY0=DataProcessor.harpyInterface.ComputeXSec(s)
        YYlist.append(YY0)
        
    YY=numpy.mean(YYlist,axis=0)
    YYstd=numpy.std(YYlist,axis=0)
        
    print('{"'+s.name+'",{')
    for i in range(s.numberOfPoints):
        print("{"+"{:4d},{:4d},{:2.4f},{:12.6f},{:12.6f},{:12.6f},{:12.6f}".format(
            -1,
            -1,
            s.points[i]["<qT>"],
            s.points[i]["xSec"],numpy.sqrt(numpy.sum(numpy.array(s.points[i]["uncorrErr"])**2)),
            YY[i],YYstd[i]
            ),end="")
        if i==s.numberOfPoints-1:
            if j==len(setDY.sets)-1:
                print("}}}}")
            else:
                print("}}},")                
        else:
            print("},")

#%%
def MinForReplica():
    
    
    def repchi_2(x):        
        startT=time.time()
        harpy.setNPparameters_SiversTMDPDF(x)
        print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
        
        ccDY2,cc3=0,0#DataProcessor.harpyInterface.ComputeChi2(repDataDY)
        ccSIDIS2,cc3=DataProcessor.harpyInterface.ComputeChi2(repDataSIDIS,method="central")
        
        cc=(ccDY2+ccSIDIS2)/totalNnew
        endT=time.time()
        print(':->',cc,'       t=',endT-startT)
        return ccSIDIS2+5*ccDY2
    
    #repDataDY=setDY.GenerateReplica()
    repDataSIDIS=setSIDIS.GenerateReplica()
    totalNnew=repDataSIDIS.numberOfPoints#+repDataDY.numberOfPoints
    
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
numOfReplicas=100
REPPATH=MAINPATH+"FittingPrograms/Sivers20/LOGS/"+"COSH-MODEL2-hermes+compass-replicas+++.txt"
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