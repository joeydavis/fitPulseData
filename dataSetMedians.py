# -*- coding: utf-8 -*-
"""
Created on Mon May 13 22:19:48 2013

@author: jhdavis
"""
import qMS

def readCSV(filename, pulse=False):
    data = qMS.readIsoCSV(filename, pulse)
    return data

def calcValue(data, num, den, offset=0.0):
    ns = [data.loc[x] for x in num]
    valNum = reduce(lambda x,y:x+y, ns)
    ds = [data.loc[x] for x in den]
    valDen = reduce(lambda x,y:x+y, ds)
    return float(valNum)/float(valDen) + offset
    
def transformValue(inputData, transformData, toTransform, normalizeTo):
    inputNorm = float(inputData.loc[toTransform])/float(inputData.loc[normalizeTo])
    transformNorm = float(transformData.loc[toTransform])/float(transformData.loc[normalizeTo])
    inputNormOffset = inputNorm - transformNorm
    inputTransformed = inputNormOffset*float(inputData.loc[normalizeTo])
    return inputTransformed
    
if __name__ == '__main__':
    listOfFiles = ['S'+str(i) for i in range(17,39)]
    datapath = '/home/jhdavis/data/2013_03_11-NCMGradients/SingleSpike/'
    
    datafile = datapath+listOfFiles[0]+'/'+listOfFiles[0]+'_iso.csv'
    dataFrame = readCSV(datafile, pulse=True)
    calcNum = ['AMP_U', 'AMP_L']
    calcDen = ['AMP_U', 'AMP_L', 'AMP_S']