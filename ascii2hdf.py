

import h5py 
import numpy as np
import datetime
import time
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit  
import warnings
warnings.filterwarnings("ignore")
import glob 
import sys 

#from DAQscripts.general.HDFwriterGen2DOMPMTinspection import HDFwriterGen2DOMPMTinspection
from HDFwriterGen2DOMPMTinspection import HDFwriterGen2DOMPMTinspection_Ex

if(__name__ == "__main__"):


    PMTID    = 9883
    MCUID    = 1234
    pwmfreq  = 12345
    temperature = 25
    dacvalue = 600
    
    wubaseID = -1
    userdict = {"Lasersetting":-1, "position":np.array([0,0,0]), "B-field":np.array([0,0,0])}


    # 2.0V/4096bits/1670ohm/60MHz
    conversion_ch1 = 2.0*1.85/4096/1670.0/60e6*1e12 # pC/1 ADC sample calibrated by micro base
    conversion_ch2 = conversion_ch1/0.0115 # measured by Chiba
    yourname = "Chiba"
    description = "Data collected with no injection."


    filename    = sys.argv[1]
    outfilename = sys.argv[2]
    PMTID = sys.argv[3]
    print(f"Reading {filename}")

    PMTIDstr = "BB{0}".format(PMTID)


    WAIT_FOR_HIT=0
    NSAMPLES_NEXT=1
    TIMESTAMP_NEXT=2
    TDCWORD_NEXT=3
    CH1_NEXT=4
    CH1_IN_PROGRESS=5
    CH2_IN_PROGRESS=6
    UID_NEXT = 7
    FREQ_NEXT = 8
    TEMP_NEXT = 9
    readstate=WAIT_FOR_HIT
    Nsample_max = 100
    Nwfm = 0

    hitData = []
    vol=[]
    Timestamp = []
    Tcword    = []
    Nsample = []

    f = open(filename, "r")

    for line in f:
        linestrip=line.rstrip()
        tokens=linestrip.split()
        if readstate==WAIT_FOR_HIT:
            if tokens[0]=='Ubaud':
                Ubaud = tokens[1]
            if tokens[0]=='get_Uid':
                readstate=UID_NEXT
            if tokens[0]=='get_freq':
                readstate=FREQ_NEXT
                #freq = tokens[1]
            if tokens[0]=='get_avg_temp':
                readstate=TEMP_NEXT
                #temperature = tokens[1]
            if tokens[0]=='voltage':
                vdata = tokens[1]
                vol.append(vdata)
            if tokens[0]=='dac':
                dacvalue = tokens[2]
            if tokens[0]=='fpgatrig':
                runtype = tokens[1]
            if len(tokens)!=1:continue  #to the next data 
            if tokens[0]!="V1" and tokens[0]!="V1X": continue
            
            readstate=NSAMPLES_NEXT   #len(tokens)==1 and V1X
            continue 
        elif readstate==UID_NEXT:
            if len(tokens)!=3:
                readstate=WAIT_FOR_HIT #unexpected input
            else:
                MCUID    = tokens[0]+tokens[1]+tokens[2]
                readstate=WAIT_FOR_HIT 
            continue
        elif readstate==FREQ_NEXT:
            if len(tokens)!=1:
                readstate=WAIT_FOR_HIT #unexpected input
            else:
                pwmfreq    = tokens[0]
                readstate=WAIT_FOR_HIT 
            continue
        elif readstate==TEMP_NEXT:
            if len(tokens)!=1:
                readstate=WAIT_FOR_HIT #unexpected input
            else:
                temperature    = tokens[0]
                readstate=WAIT_FOR_HIT 
            continue
        elif readstate==NSAMPLES_NEXT:
            if len(tokens)!=1:
                readstate=WAIT_FOR_HIT #unexpected input
            else:
                nsamples=int(tokens[0],16)  #use
                readstate=TIMESTAMP_NEXT
            if nsamples>Nsample_max:
                readstate=WAIT_FOR_HIT # unexpected input
            continue
        elif readstate==TIMESTAMP_NEXT:
            if len(tokens)!=1:
                readstate=WAIT_FOR_HIT #unexpected input
            else:
                timestamp =int(tokens[0],16)   #use
                readstate=TDCWORD_NEXT
            continue

        elif readstate==TDCWORD_NEXT:

            if len(tokens)!=1:
                readstate=WAIT_FOR_HIT #unexpected input
            else:
                tcword = int(tokens[0],16)   #use
                readstate=CH1_NEXT
            continue

        elif readstate==CH1_NEXT:
            ch1=[]
            ch2=[]
            discraw=[]
            discsync=[]
            for hexvalue in tokens:   #the first value of ch1
                raw=int(hexvalue,16)&0xfff  #one bin value
                if len(ch1)>=nsamples:
                    readstate=WAIT_FOR_HIT #more values than expected
                    continue
                ch1.append(raw)
                if len(ch1)==nsamples:
                    readstate=CH2_IN_PROGRESS
                else:
                    readstate=CH1_IN_PROGRESS
                continue

        elif readstate==CH1_IN_PROGRESS:
            for hexvalue in tokens:
                raw=int(hexvalue,16)&0xfff
                if len(ch1)>=nsamples:
                    readstate=WAIT_FOR_HIT #more values than expected
                    continue
                ch1.append(raw)
            if len(ch1)==nsamples:  #complated ch1
                while len(ch1)!=Nsample_max:
                    ch1.append(0)
                if len(ch1)==Nsample_max:
                    Ch1 = np.asarray(ch1).astype(np.uint16)
                    Nwfm  += 1
                if Nwfm == 1:
                    adcs_ch1 = Ch1
                elif Nwfm > 1:
                    adcs_ch1 = np.vstack((adcs_ch1,Ch1))
                    
                readstate=CH2_IN_PROGRESS
            else:
                readstate=CH1_IN_PROGRESS
            continue
        elif readstate==CH2_IN_PROGRESS:
            for hexvalue in tokens:
                raw=int(hexvalue,16)
                if len(ch2)>=nsamples:
                    readstate=WAIT_FOR_HIT #more values than expected
                    continue
                ch2.append(raw&0xfff)
                discsync.append((raw>>13)&1)
                discraw.append((raw>>12)&1)
                
            if len(ch2)==nsamples:
                while len(ch2)!=Nsample_max:
                    ch2.append(0)

                if len(ch2)==Nsample_max:
                    Ch2 = np.asarray(ch2).astype(np.uint16)
                    #Nwfm  += 1
                if Nwfm == 1:
                    adcs_ch2 = Ch2
                elif Nwfm > 1:
                    adcs_ch2 = np.vstack((adcs_ch2,Ch2))

                Nsample.append(nsamples)
                Timestamp.append(timestamp)
                Tcword.append(tcword)
                readstate=WAIT_FOR_HIT # that's next after pPMTID

    
    #runtype = 0
    voltage10 = vol[0]
    

    # a class to save HDF file
    hdf = HDFwriterGen2DOMPMTinspection_Ex(Nsample_max)

    # set the metadata
    hdf.fill_metadata(PMTID, PMTIDstr, wubaseID, MCUID, yourname, runtype, voltage10, \
         pwmfreq, dacvalue, temperature, conversion_ch1, conversion_ch2, description, userdict)

    # number of samples. in this dummy creation, this is fixed, but real data would dynamically change.
    Nsample = np.array(Nsample).astype(np.uint16)
    
    # FPGA time
    FPGA_time = np.array(Timestamp).astype(np.uint64)
    FPGA_tcword = np.array(Tcword).astype(np.uint64)

    ## ------------- end of dummy data creation ------------------##

    # save the waveform and their basic values
    hdf.fill(Nsample, FPGA_time, FPGA_tcword, adcs_ch1, adcs_ch2)
    
    # fit Makino function to precisely evealuate charge/time/peak.
    hdf.fit_v0()

    hdf.write("{0}".format(outfilename))
