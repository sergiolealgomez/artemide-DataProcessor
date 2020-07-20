#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 12:06:14 2020

Class ArtemideReplicaSet incorporates the information 
about a particular replica set, encode the to parametrize TMDs in artemide,
it reads the replica files and allow to operate with it.

The object ArtemideReplicaSet has following fields:
    
numberOfReplicas:  the number of replicas


@author: vla18041
"""

def ReadRepFile(path):
    import os.path
    if not os.path.exists(path):
        raise FileNotFoundError('replica-file '+path+' NOT FOUND')
    
    
    with open(path,"r") as file:
        listFromF=file.readlines()
    
    ### I read line by line and delete lines fron the list
    ### it corresponds to the logical way I saved it in the fortran
    ### the line start with *?? indcate that the next line is the variable
    
    ### search for name entry
    while not listFromF[0].startswith("*A   "):
        listFromF.pop(0)
    listFromF.pop(0)
        
    name=(listFromF.pop(0)).replace("\n","")
    rSet=ArtemideReplicaSet()
    rSet.name=name
    
    ### search for indexing parameters
    while not listFromF[0].startswith("*B   "):
        listFromF.pop(0)
    ## total size
    while not listFromF[0].startswith("*0   "):
        listFromF.pop(0)
    listFromF.pop(0)
    
    rSet._totalLength=int(listFromF.pop(0))
    
    ## TMDR size
    while not listFromF[0].startswith("*3   "):
        listFromF.pop(0)
    listFromF.pop(0)
    
    line=(listFromF.pop(0)).split(",")
    
    rSet._TMDRstart=int(line[0])-2
    rSet._TMDRend=int(line[1])-1
    
    ## TMDPDF size
    while not listFromF[0].startswith("*4   "):
        listFromF.pop(0)
    listFromF.pop(0)
    
    line=(listFromF.pop(0)).split(",")
    rSet._uTMDPDFstart=int(line[0])-2
    rSet._uTMDPDFend=int(line[1])-1
    
    ## TMDFF size
    while not listFromF[0].startswith("*5   "):
        listFromF.pop(0)
    listFromF.pop(0)
    
    line=(listFromF.pop(0)).split(",")
    rSet._uTMDFFstart=int(line[0])-2
    rSet._uTMDFFend=int(line[1])-1
    
    ## TMDFF size
    while not listFromF[0].startswith("*10  "):
        listFromF.pop(0)
    listFromF.pop(0)
    
    line=(listFromF.pop(0)).split(",")
    rSet._lpTMDPDFstart=int(line[0])-2
    rSet._lpTMDPDFend=int(line[1])-1
    
    ### search for number of replicas
    while not listFromF[0].startswith("*C   "):
        listFromF.pop(0)
    listFromF.pop(0)
    
    rSet.numberOfReplicas=int(listFromF.pop(0))
    
    ### search for technical replicas
    while not listFromF[0].startswith("*D   "):
        listFromF.pop(0)
    listFromF.pop(0)
    
    line=listFromF.pop(0).split(",")
    rSet.initialReplica=[float(x) for x in line[1:]]
    
    line=listFromF.pop(0).split(",")
    rSet.meanReplica=[float(x) for x in line[1:]]
    
    ### search for the start of the replica lst
    while not listFromF[0].startswith("*R   "):
        listFromF.pop(0)
    listFromF.pop(0)
    
    rSet.replicaList=[]
    for line in listFromF[:rSet.numberOfReplicas]:
        lineX=line.split(",")
        rSet.replicaList.append([float(x) for x in lineX[1:]])
    
    return rSet
        

class ArtemideReplicaSet:
    def ArtemideReplicaSet(self):
        
        self.name="NONAME"
        self.numberOfReplicas=0
        
        self._totalLength=0
        
        self._TMDRstart=0
        self._uTMDPDFstart=0
        self._uTMDFFstart=0
        self._lpTMDPDFstart=0
        
        self._TMDRend=0
        self._uTMDPDFend=0
        self._uTMDFFend=0
        self._lpTMDPDFend=0
        
        ### replica suggested for the initialization
        self.initialReplica=[]
        ### mean replica
        self.meanReplica=[]
        ### all other replicas
        self.replicaList=[]
        
    def __repr__ (self):
        return "<ArtemideReplicaSet: %s with %s replicas>" % (self.name, self.numberOfReplicas)
    
    
    def SetReplica(self,num=0):
        """
        Set the replica according to the list

        Parameters
        ----------
        num : -1 = initla replica
              0  = mean replica
              1.... = replica from list
            The default is 0.

        Returns
        -------
        None.

        """        
        import harpy
        
        if(num==0):
            r=self.meanReplica
        elif(num==-1):
            r=self.initialReplica
        else:
            r=self.replicaList[num-1]
        
        if(self._TMDRend>self._TMDRstart+1):
            harpy.setNPparameters_TMDR(r[self._TMDRstart:self._TMDRend])
        if(self._uTMDPDFend>self._uTMDPDFstart+1):
            harpy.setNPparameters_uTMDPDF(r[self._uTMDPDFstart:self._uTMDPDFend])
        if(self._uTMDFFend>self._uTMDFFstart+1):
            harpy.setNPparameters_uTMDFF(r[self._uTMDFFstart:self._uTMDFFend])
        if(self._lpTMDPDFend>self._lpTMDPDFstart+1):
            harpy.setNPparameters_lpTMDPDF(r[self._lpTMDPDFstart:self._lpTMDPDFend])
            
    def GetReplica(self,num,part="full"):
        if(num==0):
            r=self.meanReplica
        elif(num==-1):
            r=self.initialReplica
        else:
            r=self.replicaList[num-1]
        
        if(part=="full"):
            return r
        elif(part=="TMDR"):
            return r[self._TMDRstart:self._TMDRend]
        elif(part=="uTMDPDF"):
            return r[self._uTMDPDFstart:self._uTMDPDFend]
        elif(part=="uTMDFF"):
            return r[self._uTMDFFstart:self._uTMDFFend]
        elif(part=="lpTMDPDF"):
            return r[self._lpTMDPDFstart:self._lpTMDPDFend]
        else:
            raise ValueError("part should correspond to a TMD")