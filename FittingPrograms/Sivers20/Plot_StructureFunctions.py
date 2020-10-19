#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 09:04:35 2020

@author: vla18041
"""

#######################################
# importing libraries
#######################################

import sys
import numpy
#sys.path.append("/home/m/Github/artemide-DataProcessor")
sys.path.append("/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor")
import DataProcessor.harpyInterface
import DataProcessor.DataMultiSet
import DataProcessor.DataSet
import DataProcessor.ArtemideReplicaSet

MAINPATH="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/OtherPrograms/EIC_impact_Study/"
#%%
rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile("/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/Replicas/DY+SIDIS/SV19_nnlo_all=0.rep")

#%%
#######################################
#Initialize artemide
#######################################
import harpy
harpy.initialize(MAINPATH+"/Constants-files/"+"const-NNPDF31+DSS_nnlo_reduce")
rSet.SetReplica(0)

#%%
def cutFuncFUU(p):
    
    from math import log10, floor
    def round_2dg(x):
        return round(x, 1-int(floor(log10(abs(x)))))
    #p["<x>"]=numpy.round(numpy.mean(p["x"]),5)
    #p["<z>"]=numpy.round(numpy.mean(p["z"]),3)
    #p["<pT>"]=numpy.round(numpy.mean(p["pT"]),2)
    
    p["<x>"]=round_2dg(numpy.mean(p["x"]))
    p["<z>"]=round_2dg(numpy.mean(p["z"]))
    p["<pT>"]=round_2dg(numpy.mean(p["pT"]))
    ### for it I recompute normalization (process changes as p[1,1,a]->p[1,2,a])
    p["process"][1]=2
    
    ### compute the theory prediciton
    th=DataProcessor.harpyInterface.ComputeXSec(p,method="central")      
    ### rescale sys.error by new theory
    p["uncorrErr"]=[u*th/numpy.abs(p["xSec"]) for u in p["uncorrErr"]]
    ### replace xSec by theory prediciton 
    p["xSec"]=th
    
    p["thFactor"]=1.
    return True, p


#%%
################################################
### Load all unpolarized data
###############################################
import os
##load
path_to_load="/home/vla18041/LinkData2/WorkingFiles/TMD/YR_Studies/Repository/EIC_YR_TMD/Unpolarized/PseudoData/Data4_cut"
##### for all directories
loadedSets=[]
listOfDir=os.listdir(path_to_load)
for dirName in listOfDir:
    ## all directories strat with `ep'
    if(dirName[0:2]!='ep'):
        continue
    listOfFiles=os.listdir(path_to_load+'/'+dirName)
    for fileName in listOfFiles:
        ### Ralf's data files have names like
        ### 'ep_noradcor.10x100_pim_ACC_opt3.txt'
        if(fileName[0:12]!='ep_noradcor.'):
            print("File <",fileName,"> skipped.")
            continue
        if(fileName.split("_")[2]!="pip"):
            #print("File <",fileName,"> skipped.")
            continue
        if(fileName.split("_")[3]!="ACC"):
            #print("File <",fileName,"> skipped.")
            continue
        if(fileName.split("_")[4]!="opt5.csv"):
            #print("File <",fileName,"> skipped.")
            continue
        
        print("Working file: "+dirName+'/'+fileName)
        loadedSets.append(DataProcessor.DataSet.LoadCSV(path_to_load+'/'+dirName+'/'+fileName))

set0=DataProcessor.DataMultiSet.DataMultiSet("SIDISset",loadedSets)
setFUU=set0.CutData(cutFuncFUU,computeCovarianceMatrix=False) 

#%%%
################################################
### Load empty unpolarized data
###############################################
import os
##load
path_to_load="/home/vla18041/LinkData2/WorkingFiles/TMD/YR_Studies/Repository/EIC_YR_TMD/Unpolarized/PseudoData/DataEmpty_cut"
##### for all directories
loadedSets=[]
listOfFiles=os.listdir(path_to_load)
for fileName in listOfFiles:        
    print("Working file: "+'/'+fileName)
    loadedSets.append(DataProcessor.DataSet.LoadCSV(path_to_load+'/'+fileName))

set0=DataProcessor.DataMultiSet.DataMultiSet("SIDISset",loadedSets)
setFUUempty=set0.CutData(cutFuncFUU,computeCovarianceMatrix=False) 

#%%
#######################
# Joining empty and HB set, to fill the points which are not covered by acceptance
#######################
import copy
newPoints=[]

def comparepoints(p1,p2):
    # return (numpy.abs(p1["x"][0]-p2["x"][0])/p1["<x>"]<0.001 and numpy.abs(p1["x"][1]-p2["x"][1])/p1["<x>"]<0.001
    #    and numpy.abs(p1["pT"][0]-p2["pT"][0])/p1["<pT>"]<0.001 and numpy.abs(p1["pT"][1]-p2["pT"][1])/p1["<pT>"]<0.001
    #     and numpy.abs(p1["z"][0]-p2["z"][0])/p1["<z>"]<0.001 and numpy.abs(p1["z"][1]-p2["z"][1])/p1["<z>"]<0.001
    #     and numpy.abs(p1["Q"][0]-p2["Q"][0])/p1["<Q>"]<0.001 and numpy.abs(p1["Q"][1]-p2["Q"][1])/p1["<Q>"]<0.001)
    return (p1["<x>"]==p2["<x>"] and  p1["<z>"]==p2["<z>"] and p1["<pT>"]==p2["<pT>"]
         and numpy.abs(p1["Q"][0]-p2["Q"][0])/p1["<Q>"]<0.001 and numpy.abs(p1["Q"][1]-p2["Q"][1])/p1["<Q>"]<0.001)
def isPointInSet(s,p):
    for pp in s.points:
        if(comparepoints(p, pp)):
            return True
    return False

for p in setFUUempty.points:
    if(not isPointInSet(setFUU,p)):
        newPoints.append(p)

newSet=DataProcessor.DataSet.DataSet("addendum","SIDIS")
for p in newPoints:
    newSet.AddPoint(p)

newSet.FinalizeSet(computeCovarianceMatrix=False)    
newSets=copy.deepcopy(setFUU.sets)
newSets.append(newSet)

completeFUU=DataProcessor.DataMultiSet.DataMultiSet("SIDISset",newSets)
#%%
xztcases = [] 
[xztcases.append(x) for x in [[p["<x>"],p["<z>"],p["<pT>"]] for p in setFUU.points] if x not in xztcases]

#%%
xcases = [] 
[xcases.append(x) for x in [p[0] for p in xztcases] if x not in xcases]

zcases = [] 
[zcases.append(x) for x in [p[1] for p in xztcases] if x not in zcases]

ptcases = [] 
[ptcases.append(x) for x in [p[2] for p in xztcases] if x not in ptcases]


#%%
def ExtractXZTcase(xzt,ss):
    points=[]
    for p in ss.points:
        if([p["<x>"],p["<z>"],p["<pT>"]]==xzt):
            points.append(p)
    return points

def ExtractXZTcase2(xzt,ss):
    points=[]
    for p in ss:
        if([p[0],p[1],p[2]]==xzt):
            points.append(p)
    return points

#%%
from matplotlib import pyplot
#fig, ax = pyplot.subplots()
fig, ax = pyplot.subplots(figsize=(5, 5*2))
for c in xztcases:
    if(c[2]==0.025 and c[1]==0.35):
        pp=ExtractXZTcase(c,setFUU)
        ppE=ExtractXZTcase(c,completeFUU)
        if(len(pp)>1):
            curve=numpy.transpose([[p["<Q>"]**2 for p in ppE],[p["xSec"]/p["<x>"] for p in ppE]])
            curve=numpy.transpose(curve[curve[:,0].argsort()])
            #print(curve)
            #pyplot.plot([p["<Q>"]**2 for p in pp],[p["xSec"] for p in pp],label=c[0])
            pyplot.plot(curve[0],curve[1],label="x="+str(c[0]))
            
            pyplot.errorbar([p["<Q>"]**2 for p in pp], [p["xSec"]/p["<x>"] for p in pp], 
                            yerr=[numpy.sqrt(p["uncorrErr"][0]**2+p["uncorrErr"][1]**2)/p["<x>"] for p in pp], fmt='o')

pyplot.yscale("log")
pyplot.xscale("log")
pyplot.xlim(2,1800)
pyplot.legend(loc="upper right")
pyplot.title("FUU z=0.35 pT=0.02" )
pyplot.xlabel("Q^2" )
pyplot.ylabel("FUU/x" )
pyplot.show()

#%%
def cutToDraw(p):
    if(p["<z>"]==0.35 and p["<pT>"]==0.025):
        return True, p
    return False, p

subsetToDrawTH=completeFUU.CutData(cutToDraw,computeCovarianceMatrix=False) 
subsetToDrawEX=setFUU.CutData(cutToDraw,computeCovarianceMatrix=False) 

#%%
from matplotlib import pyplot
#fig, ax = pyplot.subplots()
fig, ax = pyplot.subplots(figsize=(5, 5*2))
for c in xztcases:    
    pp=ExtractXZTcase(c,subsetToDrawEX)
    ppE=ExtractXZTcase(c,subsetToDrawTH)
    if(len(pp)>1):
        curve=numpy.transpose([[p["<Q>"]**2 for p in ppE],[p["xSec"]/p["<x>"] for p in ppE]])
        curve=numpy.transpose(curve[curve[:,0].argsort()])
        #print(curve)
        #pyplot.plot([p["<Q>"]**2 for p in pp],[p["xSec"] for p in pp],label=c[0])
        pyplot.plot(curve[0],curve[1],label="x="+str(c[0]))
        
        pyplot.errorbar([p["<Q>"]**2 for p in pp], [p["xSec"]/p["<x>"] for p in pp], 
                        yerr=[numpy.sqrt(p["uncorrErr"][0]**2+p["uncorrErr"][1]**2)/p["<x>"] for p in pp], fmt='o')

pyplot.yscale("log")
pyplot.xscale("log")
pyplot.xlim(2,1800)
pyplot.legend(loc="upper right")
pyplot.title("FUU z=0.35 pT=0.02" )
pyplot.xlabel("Q^2" )
pyplot.ylabel("FUU/x" )
pyplot.show()

#%%
central0=DataProcessor.harpyInterface.ComputeXSec(subsetToDrawTH,method="central")
ss=[]
for i in range(rSet.numberOfReplicas):
    rSet.SetReplica(i+1)
    ss.append(DataProcessor.harpyInterface.ComputeXSec(subsetToDrawTH,method="central"))
rSet.SetReplica(0)
central=numpy.mean(ss,axis=0)
std=numpy.std(ss,axis=0)
THpoints=numpy.transpose([[p["<x>"] for p in subsetToDrawTH.points],[p["<z>"] for p in subsetToDrawTH.points],[p["<pT>"] for p in subsetToDrawTH.points],
          [p["<Q>"] for p in subsetToDrawTH.points],central0,std])

#%%
from matplotlib import pyplot
#fig, ax = pyplot.subplots()
fig, ax = pyplot.subplots(figsize=(5, 5*2))
for c in xztcases:    
    pp=ExtractXZTcase(c,subsetToDrawEX)
    ppE=ExtractXZTcase2(c,THpoints)
    if(len(pp)>1):
        curve=numpy.transpose([[p[3]**2 for p in ppE],[p[4]/p[0] for p in ppE]])
        curvePLUS=numpy.transpose([[p[3]**2 for p in ppE],[(p[4]/p[0]+p[5]/p[0]) for p in ppE]])
        curveMINUS=numpy.transpose([[p[3]**2 for p in ppE],[(p[4]/p[0]-p[5]/p[0]) for p in ppE]])
        
        curve=numpy.transpose(curve[curve[:,0].argsort()])
        curvePLUS=numpy.transpose(curvePLUS[curvePLUS[:,0].argsort()])
        curveMINUS=numpy.transpose(curveMINUS[curveMINUS[:,0].argsort()])
        
        
                
        #print(curve)
        #pyplot.plot([p["<Q>"]**2 for p in pp],[p["xSec"] for p in pp],label=c[0])
        pyplot.fill_between(curve[0],curveMINUS[1],curvePLUS[1],color="orange")
        pyplot.plot(curve[0],curve[1],label="x="+str(c[0]))
        
        
        
        pyplot.errorbar([p["<Q>"]**2 for p in pp], [p["xSec"]/p["<x>"] for p in pp], 
                        yerr=[numpy.sqrt(p["uncorrErr"][0]**2+p["uncorrErr"][1]**2)/p["<x>"] for p in pp], fmt='o')

pyplot.yscale("log")
pyplot.xscale("log")
pyplot.xlim(2,1800)
pyplot.ylim(0.1,2000)
pyplot.legend(loc="upper right")
pyplot.title("FUU z=0.35 pT=0.02" )
pyplot.xlabel("Q^2" )
pyplot.ylabel("FUU/x" )
pyplot.show()


#%%
for p in THpoints:
    #print("{"+"{:12.8f} ,{:12.8f} ,{:12.8f},{:12.8f}".format(p[0],p[3],p[4]/p[0],p[5]/p[0])+"},")
    print("{"+"{:12.8f} ,{:12.8f} ,{:12.8f},{:12.8f}".format(p[0],p[3],p[4],p[5])+"},")
    
#%%
for p in subsetToDrawEX.points:
    print("{"+"{:12.8f} ,{:12.8f} ,{:12.8f} ,{:12.8f}".format(
        p["<x>"],p["<Q>"],p["xSec"]/p["<x>"],numpy.sqrt(p["uncorrErr"][0]**2+p["uncorrErr"][1]**2)/p["<x>"])+"},")
        #p["<x>"],p["<Q>"],p["xSec"],numpy.sqrt(p["uncorrErr"][0]**2+p["uncorrErr"][1]**2))+"},")