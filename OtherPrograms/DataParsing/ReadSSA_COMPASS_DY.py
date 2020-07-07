#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 20:15:46 2019
This program collect all the data on the SIDIS and save in  "SIDISdata_uncut.pkl"
@author: vla18041
"""
import sys
sys.path.append("/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/")
import DataProcessor.Point
import DataProcessor.DataSet
import numpy

########### Alexey desktop
path_to_data="/home/vla18041/LinkData2/arTeMiDe_Repository/data"
path_to_COMPASS="/COMPASS/1704.00488/"
path_to_save="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/DataLib/Sivers/"

########## Marcin laptop
#path_to_data="/home/m/Dropbox/Sivers/Data"
#path_to_COMPASS="/COMPASS08/"

M_proton=0.938
m_pion=0.139

### Scale uncertanty
targetPolarizationUncertanty=0.05
dilutionUncertanty=0.08

#%%
###############################################################################
f = open(path_to_data+path_to_COMPASS+"A_T_sin_phiS_M.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
######    d Q            ######
data_current=data_from_f[4:7]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";","\t")
    data_current[i]=data_current[i].split("\t")    
    data_current[i]=[float(j) for j in data_current[i]]
    
print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('compass.sivers.piDY.dQ',"DY")
DataCurrent.comment="COMPASS A_T pionDY (differential in Q)"
DataCurrent.reference="1704.00488"

proc_current=[2,2,10101]
proc_denominator=[2,2,101]
s_current=M_proton**2+m_pion**2+2.*M_proton*190.
includeCuts=False
cutParameters=[25.,25.,-1.,1.] #pT1,pT2,yMin,yMax
DataCurrent.normErr.append(targetPolarizationUncertanty)
DataCurrent.normErr.append(dilutionUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))    
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<qT>"]=data_current[i][5]    
    p1["qT"]=[0.0000001,5.]    
    p1["<Q>"]=data_current[i][6]    
    p1["Q"]=[data_current[i][0],data_current[i][1]] ### I have checked that +- 40 gives ~10^{-4} error in comparison to the total volume
    p1["<y>"]=data_current[i][4]
    p1["y"]=[-1.,1.]
    p1["xSec"]=data_current[i][7]
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][9])
    p1["weightProcess"]=proc_denominator
    #
    DataCurrent.AddPoint(p1)    

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
f = open(path_to_data+path_to_COMPASS+"A_T_sin_phiS_x_F.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
######    d Q            ######
data_current=data_from_f[4:7]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";","\t")
    data_current[i]=data_current[i].split("\t")    
    data_current[i]=[float(j) for j in data_current[i]]
    
print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('compass.sivers.piDY.dxF',"DY")
DataCurrent.comment="COMPASS A_T pionDY (differential in xF)"
DataCurrent.reference="1704.00488"

proc_current=[2,2,10101]
proc_denominator=[2,2,101]
s_current=M_proton**2+m_pion**2+2.*M_proton*190.
includeCuts=False
cutParameters=[25.,25.,-1.,1.] #pT1,pT2,yMin,yMax
DataCurrent.normErr.append(targetPolarizationUncertanty)
DataCurrent.normErr.append(dilutionUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))    
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<qT>"]=data_current[i][5]    
    p1["qT"]=[0.0000001,5.]    
    p1["<Q>"]=data_current[i][6]    
    p1["Q"]=[4.3,8.5]
    p1["<y>"]=data_current[i][4]
    p1["y"]=[data_current[i][0],data_current[i][1]]
    p1["xSec"]=data_current[i][7]
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][9])
    p1["weightProcess"]=proc_denominator
    #
    DataCurrent.AddPoint(p1)    

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
f = open(path_to_data+path_to_COMPASS+"A_T_sin_phiS_q_T.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
######    d Q            ######
data_current=data_from_f[4:7]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";","\t")
    data_current[i]=data_current[i].split("\t")    
    data_current[i]=[float(j) for j in data_current[i]]
    
print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('compass.sivers.piDY.dqT',"DY")
DataCurrent.comment="COMPASS A_T pionDY (differential in qT)"
DataCurrent.reference="1704.00488"

proc_current=[2,2,10101]
proc_denominator=[2,2,101]
s_current=M_proton**2+m_pion**2+2.*M_proton*190.
includeCuts=False
cutParameters=[25.,25.,-1.,1.] #pT1,pT2,yMin,yMax
DataCurrent.normErr.append(targetPolarizationUncertanty)
DataCurrent.normErr.append(dilutionUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))    
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<qT>"]=data_current[i][5]    
    p1["qT"]=[data_current[i][0],data_current[i][1]]
    p1["<Q>"]=data_current[i][6]    
    p1["Q"]=[4.3,8.5]
    p1["<y>"]=data_current[i][4]
    p1["y"]=[-1.,1.]
    p1["xSec"]=data_current[i][7]
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][9])
    p1["weightProcess"]=proc_denominator
    #
    DataCurrent.AddPoint(p1)    

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")