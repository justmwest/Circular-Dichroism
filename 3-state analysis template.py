#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 12 07:35:58 2019

@author: Justin
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob

# Pull all csv filenames in the folder containing the script and print them
extension = 'csv'
Fnames = [i for i in glob.glob('*.{}'.format(extension))]
Fnames=sorted(Fnames)
print(("Found ")+str(len(Fnames))+(" files:"))     
for filename in Fnames:
    print(filename)

#import mdeg values
def import_mdeg(file_number,high_wavelength = 260,low_wavelength = 195):
    '''
    input is a file number
    output is a list of mdeg
    '''
    file = Fnames[file_number]
    nRows=(high_wavelength-low_wavelength)+1
    global wavelengths
    wavelengths = np.array(pd.read_csv(file,usecols=[0], skiprows=18,nrows=nRows))
    mdeg = np.array(pd.read_csv(file,usecols=[1], skiprows=18,nrows=nRows))
    return mdeg

#plot function will take mdeg as input
def plot(mdeg,plot_title = "",mdegErr=None):
    '''
    Plots data from file number only
    input is an integer within the len(files)
    output is a plt.plot set up, user just needs to call plt.show()
    ''' 
    plt.title(plot_title)
    plt.xlabel('λ (nm)')
    plt.ylabel('mdeg')
    plt.errorbar(wavelengths,mdeg,mdegErr)
    #plt.show() #uncomment this to show automatically when function is called

#import all the data
hiλ = 260
loλ = 198
mdegList = []
for i in range(len(Fnames)):
    mdegList.append(import_mdeg(i,hiλ,loλ))
        
#Show all the raw data together with plot
for i in range(len(mdegList)):
    plot(mdegList[i],Fnames[i])
plt.legend(Fnames)
plt.title("Raw")
plt.show()
 
#Define and average the data
lipidBlank=mdegList[0]
lipidBlankerr=None
bufferBlank=mdegList[1]
bufferBlankerr=None
highpH=np.mean(mdegList[5:7],axis=0)
highpHerr=np.std(mdegList[5:7],axis=0)
lowpH=np.mean(mdegList[3:5],axis=0)
lowpHerr=np.std(mdegList[3:5],axis=0)
peptideOnly=mdegList[2]
peptideOnlyerr=None
plot_legends = ["State I","State II","State III"]

#Plot raw averaged three-state mdeg
data_to_plot=[peptideOnly,highpH,lowpH]
errors_to_plot=[peptideOnlyerr,highpHerr,lowpHerr]
#for i in range(len(data_to_plot)):
#    plot(data_to_plot[i],None,errors_to_plot[i])
#plt.legend(plot_legends)
#plt.title("Averaged three-state")
#plt.show()

#Function to subtract the data, ignore None
def sub(sig,blk):
    if sig is None: return None
    elif blk is None: return sig
    else: return sig-blk

#Function to propagate the errors, ignore None
def suberr(e1,e2):
    if e1 is None: 
        if e2 is None: return None
        else: return e2
    elif e2 is None: return e1
    else: return ((e1**2)+(e2**2))**0.5

#Define the blanks    
blanks_to_plot=[bufferBlank,lipidBlank,lipidBlank]
blank_errors=[bufferBlankerr,lipidBlankerr,lipidBlankerr]

#Subtract the data and propagate errors
subtractedData=[]
subtractedErr=[]
for i in range(len(data_to_plot)):
    subtractedData.append(sub(data_to_plot[i],blanks_to_plot[i]))
    subtractedErr.append(suberr(errors_to_plot[i],blank_errors[i]))

#Plot subtracted data
data_to_plot=subtractedData
errors_to_plot=subtractedErr
for i in range(len(data_to_plot)):
    if data_to_plot[i] is not None:
        plot(data_to_plot[i],None,errors_to_plot[i])
plt.legend(plot_legends)
plt.title("Subtracted and Averaged")
plt.show()

#Plot blanks
#data_to_plot=[lipidBlank,bufferBlank]
#for i in range(len(data_to_plot)):
#    plot(data_to_plot[i])
#plt.legend(["Lipid Blank","bufferBlank"])
#plt.title("Blanks")
#plt.show()

#Basline Subtraction
def bsub(sig):
    if sig is None: return None
    else: return sig-sig[0]

#Baseline subtract the data
bsmdeg=[]
for curve in subtractedData:
    bsmdeg.append(bsub(curve))

#plot baseline subtracted data
#data_to_plot=bsmdeg
#for i in range(len(data_to_plot)): 
#    plot(data_to_plot[i],None,errors_to_plot[i])
#plt.legend(plot_legends)
#plt.title("Baseline Subtracted Data")
#plt.show()

#MRE conversion function
def MRE(mdeg,path=1,nAA=41,M=5):
    if mdeg is None: return None
    else: return mdeg*abs(1/(path*(nAA-1)*M)/1000)

#Convert to MRE
concList=[7,7,7]
concList=[x/1000000 for x in concList]
path=2
nAA=41

MRE_list=[]
for i in range(len(bsmdeg)):
    MRE_list.append(MRE(bsmdeg[i],path,nAA,concList[i]))

MRE_err_list=[]
for i in range(len(errors_to_plot)):
    if errors_to_plot[i] is None:
        MRE_err_list.append(None)
    else:
        MRE_err_list.append(MRE(errors_to_plot[i],path,nAA,concList[i]))

#plot MRE data        
data_to_plot=MRE_list
errors_to_plot=MRE_err_list
for i in range(len(data_to_plot)):
    plot(data_to_plot[i],None,errors_to_plot[i])
plt.legend(plot_legends)
plt.title("Mean Residue Ellipticity (output)")
plt.ylabel('MRE')
plt.show()

#Export function. Don't run it, I've already done this. 
def export(data_list,file_name="untitled.xlsx"):
    '''
    Exports data to excel file in state 1 2 3 format.
    Works even if a value is none!
    '''
    list_for_export=[]
    for each in data_list:
        if each is None:
            list_for_export.append(np.ndarray([]).flatten())
        else:
            list_for_export.append(each.flatten())
    df = pd.DataFrame(list_for_export, index=None)
    df = df.T
    df.to_excel(file_name,index=None,header=None)

export(MRE_list,"MRE_list.xlsx")
export(MRE_err_list,"MRE_err_list.xlsx")