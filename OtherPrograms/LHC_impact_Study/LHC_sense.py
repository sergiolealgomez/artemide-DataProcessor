#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 20:20:15 2020

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
import DataProcessor.ArtemideReplicaSet

MAINPATH="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/"

#%%
#######################################
#Initialize artemide
#######################################
import harpy
path_to_constants=MAINPATH+"FittingPrograms/SV19/Constants-files/"+"DY+SIDIS_nnlo/const-DY+SIDIS_NNPDF31+DSS_nnlo"

harpy.initialize(path_to_constants)

harpy.setNPparameters_TMDR([1.93, 0.0434])
harpy.setNPparameters_uTMDPDF([0.253434, 9.04351, 346.999, 2.47992, -5.69988, 0.1, 0.])
harpy.setNPparameters_uTMDFF([0.264,0.479,0.459,0.539]) 

#%%
##################Cut function
def cutFunc(p):
    if p["type"]=="DY":
        delta=p["<qT>"]/p["<Q>"]    
    return delta<0.25, p

#%%
### Loading the DY data set
theData=DataProcessor.DataSet.LoadCSV("/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/OtherPrograms/LHC_impact_Study/LHC13TeV-Empty.csv")

setDY=theData.CutData(cutFunc) 
#%%
rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
    "/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/Replicas/DY+SIDIS/DY+SIDIS_NNPDF31+DSS_nnlo.rep")
rSet.SetReplica()
#%%
######################
### This function compute correlaiton coefficient for the vector of replicas for the given set
### It test 4 functions (given by F1,F2,F3,F4)
### Returns [corr1,corr2,corr3,corr4, stdX, deltaO(rescaled), avarageXsec, Ralf's xSec]
### where corr is the correlation number ,stdX is (theory) uncertanty of the xSec, deltaO is experimental unvertanty of xSec
#####################
def CorrelationParameters(inputSet):
    
    vectorO=[]
    vectorF1=[]
    vectorF2=[]
    vectorF3=[]
    vectorF4=[]    
    vectorF5=[]
    vectorF6=[]      
    
    n1=0
    n2=0
    ## compute vectors
    for n in range(1,rSet.numberOfReplicas+1):
        rSet.SetReplica(n)
        NPparams=rSet.GetReplica(n)
        ### vector of cross-section
        XX1=DataProcessor.harpyInterface.ComputeXSec(inputSet)    
        ### vector of test function
        F1=[]
        F2=[]
        F3=[]
        F4=[]
        F5=[]
        F6=[]
        for p in inputSet.points:               
            
            F1.append(NPparams[1])
            F2.append(NPparams[2])
            F3.append(NPparams[3])            
            F4.append(NPparams[4])
            F5.append(NPparams[5])
            F6.append(NPparams[6])
            
        vectorO.append(XX1)
        vectorF1.append(F1)
        vectorF2.append(F2)
        vectorF3.append(F3)
        vectorF4.append(F4)
        vectorF5.append(F5)
        vectorF6.append(F6)
        
        n1+=1
        n2+=1
        if(n1>(rSet.numberOfReplicas/10.)):
            n1=0
            print(n2/(rSet.numberOfReplicas+1)*100.,' %')
            
    
        
    vectorO=numpy.array(vectorO)
    vectorF1=numpy.array(vectorF1)
    vectorF2=numpy.array(vectorF2)
    vectorF3=numpy.array(vectorF3)
    vectorF4=numpy.array(vectorF4)
    vectorF5=numpy.array(vectorF5)
    vectorF6=numpy.array(vectorF6)
    
    avOF1=numpy.mean(numpy.multiply(vectorO,vectorF1),axis=0)
    avOF2=numpy.mean(numpy.multiply(vectorO,vectorF2),axis=0)
    avOF3=numpy.mean(numpy.multiply(vectorO,vectorF3),axis=0)
    avOF4=numpy.mean(numpy.multiply(vectorO,vectorF4),axis=0)
    avOF5=numpy.mean(numpy.multiply(vectorO,vectorF5),axis=0)
    avOF6=numpy.mean(numpy.multiply(vectorO,vectorF6),axis=0)
    
    avO=numpy.mean(vectorO,axis=0)
    avF1=numpy.mean(vectorF1,axis=0)
    avF2=numpy.mean(vectorF2,axis=0)
    avF3=numpy.mean(vectorF3,axis=0)
    avF4=numpy.mean(vectorF4,axis=0)
    avF5=numpy.mean(vectorF5,axis=0)
    avF6=numpy.mean(vectorF6,axis=0)
    
    
    stdO=numpy.std(vectorO,axis=0)    
    stdF1=numpy.std(vectorF1,axis=0)
    stdF2=numpy.std(vectorF2,axis=0)
    stdF3=numpy.std(vectorF3,axis=0)
    stdF4=numpy.std(vectorF4,axis=0)
    stdF5=numpy.std(vectorF5,axis=0)
    stdF6=numpy.std(vectorF6,axis=0)
    
    rho1=numpy.divide(avOF1-numpy.multiply(avO,avF1),numpy.multiply(stdO,stdF1))
    rho2=numpy.divide(avOF2-numpy.multiply(avO,avF2),numpy.multiply(stdO,stdF2))
    rho3=numpy.divide(avOF3-numpy.multiply(avO,avF3),numpy.multiply(stdO,stdF3))
    rho4=numpy.divide(avOF4-numpy.multiply(avO,avF4),numpy.multiply(stdO,stdF4))
    rho5=numpy.divide(avOF5-numpy.multiply(avO,avF5),numpy.multiply(stdO,stdF5))
    rho6=numpy.divide(avOF6-numpy.multiply(avO,avF6),numpy.multiply(stdO,stdF6))
    
    ### relative experimental error
    deltaOth=numpy.divide(stdO,avO)
    
    return numpy.transpose([rho1,rho2,rho3,rho4,rho5,rho6,deltaOth,avO])

#%%
### Compute the correlation matrix and save to file
def computeAndSaveSense(setIN,path):
    rho=CorrelationParameters(setIN)
    f=open(path,"w")        
    for i in range(setIN.numberOfPoints):    
        p=setIN.points[i]
        line=("{ "
              +"{:f} , {:f} , {:f} , ".format(numpy.mean(p["Q"]),numpy.mean(p["y"]),numpy.mean(p["qT"]))
              +"{:f} , {:f} , {:f} , {:f} , {:f} , {:f} , {:f} , {:f}".format(
                  rho[i][0],rho[i][1],rho[i][2],rho[i][3],rho[i][4],rho[i][5],rho[i][6],rho[i][7])
              +"}\n")
        f.write(line)
    
    f.close()
#%%
computeAndSaveSense(setDY,"/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/OtherPrograms/LHC_impact_Study/LOGS/LHC_sense.dat")