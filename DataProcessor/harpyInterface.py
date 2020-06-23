#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May  8 14:52:27 2020

@author: vla18041
"""

import harpy
from . import DataSet
from . import DataMultiSet


def ComputeXSec(data,method="default"):
    """ Computes the cross-section values for the given data (matched to data)
    data can be: DataSet or DataMultiSet or a point-like dictionary
    methods are:
        default = usual one
        binless = evaluated with avarage values of bins only
    """
    if isinstance(data,DataSet.DataSet) or isinstance(data,DataMultiSet.DataMultiSet):
        YY=_ComputeXSec_Data(data,method)
    elif  isinstance(data,dict):
        YY=_ComputeXSec_Point(data,method)
    else:
        raise Exception('ComputeXSec works only with DataSet or DataMultiSet or a point-like dictionary')
        
    return YY

def _ComputeXSec_Data(data,method="default"):
    """ Computes the cross-section values for the given data (matched to data)
    data can be: DataSet or DataMultiSet
    methods are:
        default = usual one
        binless = evaluated with avarage values of bins only
    """
    
    if method=="default":        
        if data.processType == "DY":
            XX=harpy.DY.xSecList([d["process"] for d in data.points],
                                    [d["s"] for d in data.points],
                                    [d["qT"] for d in data.points],
                                    [d["Q"] for d in data.points],
                                    [d["y"] for d in data.points],
                                    [d["includeCuts"] for d in data.points],
                                    [d["cutParams"] for d in data.points])
            
        elif data.processType=="SIDIS":
            XX=harpy.SIDIS.xSecList([d["process"] for d in data.points],
                                    [d["s"] for d in data.points],
                                    [d["pT"] for d in data.points],
                                    [d["z"] for d in data.points],
                                    [d["x"] for d in data.points],
                                    [d["Q"] for d in data.points],
                                    [d["includeCuts"] for d in data.points],
                                    [d["cutParams"] for d in data.points],
                                    [[d["M_target"],d["M_product"]] for d in data.points])
    elif method=="binless":
        if data.processType == "DY":
            XX1=harpy.DY.xSecListBINLESS([d["process"] for d in data.points],
                                    [d["s"] for d in data.points],
                                    [d["<qT>"] for d in data.points],
                                    [d["<Q>"] for d in data.points],
                                    [d["<y>"] for d in data.points],
                                    [d["includeCuts"] for d in data.points],
                                    [d["cutParams"] for d in data.points])
            XX=[]
            ## I weight the binless computation by the width of the bin
            for i in range(len(XX1)):
                XX.append(XX1[i]*
                          (data.points[i]["Q"][1]**2-data.points[i]["Q"][0]**2)*
                          (data.points[i]["qT"][1]**2-data.points[i]["qT"][0]**2)*
                          (data.points[i]["y"][1]-data.points[i]["y"][0]))
            
        elif data.processType=="SIDIS":
            XX1=harpy.SIDIS.xSecListBINLESS([d["process"] for d in data.points],
                                    [d["s"] for d in data.points],
                                    [d["<pT>"] for d in data.points],
                                    [d["<z>"] for d in data.points],
                                    [d["<x>"] for d in data.points],
                                    [d["<Q>"] for d in data.points],
                                    [[d["M_target"],d["M_product"]] for d in data.points])
            XX=[]
            for i in range(len(XX1)):
                XX.append(XX1[i]*
                          (data.points[i]["Q"][1]**2-data.points[i]["Q"][0]**2)*
                          (data.points[i]["pT"][1]**2-data.points[i]["pT"][0]**2)*
                          (data.points[i]["x"][1]-data.points[i]["x"][0])*
                          (data.points[i]["z"][1]-data.points[i]["z"][0]))
     
    YY=data.MatchWithData(XX)
    return YY

def _ComputeXSec_Point(p,method="default",processType="default"):
    """ Computes the cross-section values for the given point
    methods are:
        default = usual one
        binless = evaluated with avarage values of bins only
    processType are:
        default = the process number is taken for process
        weight = the process number is take from weightProcess
    """
    
    if method=="default":        
        if p["type"] == "DY":
            XX1=harpy.DY.xSecList([p["process"]],[p["s"]],[p["qT"]],[p["Q"]],
                                    [p["y"]],[p["includeCuts"]],[p["cutParams"]])
            
        elif p["type"]=="SIDIS":
            XX1=harpy.SIDIS.xSecList([p["process"]],[p["s"]],[p["pT"]],[p["z"]],
                                    [p["x"]],[p["Q"]],[p["includeCuts"]],
                                    [p["cutParams"]],[[p["M_target"],p["M_product"]]])
        
        XX=XX1[0]*p["thFactor"]
        
    elif method=="binless":
        if p["type"] == "DY":
            XX1=harpy.DY.xSecListBINLESS([p["process"]],[p["s"]],[p["<qT>"]],[p["<Q>"]],
                                    [p["<y>"]],[p["includeCuts"]],[p["cutParams"]])
            ### Matching with data (note that the normalization is ignored)
            XX=XX1[0]*(p["Q"][1]**2-p["Q"][0]**2)*(p["qT"][1]**2-p["qT"][0]**2)* \
                          (p["y"][1]-p["y"][0])*p["thFactor"]
            
        elif p["type"]=="SIDIS":
            XX1=harpy.SIDIS.xSecListBINLESS([p["process"]],[p["s"]],[p["<pT>"]],[p["<z>"]],
                                    [p["<x>"]],[p["<Q>"]],[[p["M_target"],p["M_product"]]])            

            ### Matching with data (note that the normalization is ignored)
            XX=XX1[0]*(p["Q"][1]**2-p["Q"][0]**2)*(p["pT"][1]**2-p["pT"][0]**2)* \
                (p["x"][1]-p["x"][0])*(p["z"][1]-p["z"][0])*p["thFactor"]
    else:
        raise Exception('The dictionary is not a point.')
       
    return XX,[XX]

def ComputeChi2(data,method="default"):
    """ Computes the chi^2 values for the given data    
    data can be: DataSet or DataMultiSet
    """
    YY=ComputeXSec(data,method)
    
    ZZ=data.chi2(YY)
    
    if isinstance(data,DataSet.DataSet):
        return ZZ, [ZZ]
    else:
        return ZZ

def PrintChi2Table(data,method="default",processType="default"):
    """ Compute and print the table of chi^2 for experiments
    """
    import time
    
    startT=time.time()
    
    chi2T, chi2Part = ComputeChi2(data,method,processType)
    
    endT=time.time()
    
    if isinstance(data,DataSet.DataSet):
        maxLength=len(data.name)
    elif isinstance(data,DataMultiSet.DataMultiSet):
        maxLength=max([len(d.name) for d in data.sets])
    else:
        maxLength=10
    
    print("{:{width}}".format("name",width=maxLength)+' | '+' chi^2/N  '+' | ')
    print("{:-<{width}}".format("",width=maxLength)+'-|-'+'----------'+'-|-')
    #Only one set in MultiSet of just Set
    if len(chi2Part)==1:
        print("{:{width}} | {:10.3f} |".format(data.name,chi2T/data.numberOfPoints,width=maxLength))
    else:
        for i in range(len(chi2Part)):
            print("{:{width}} | {:10.3f} |".format(
                    data.sets[i].name,chi2Part[i]/data.sets[i].numberOfPoints,width=maxLength))
        print("{:-<{width}}".format("",width=maxLength)+'-|-'+'----------'+'-|-')
        print("{:{width}} | {:10.3f} |".format('Total',chi2T/data.numberOfPoints,width=maxLength))
        
    print("Computation time = ",endT-startT,' sec.')
        