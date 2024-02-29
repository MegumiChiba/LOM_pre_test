
import h5py 
import numpy as np
import datetime
import time
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit  
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import glob 
import sys 
from matplotlib.backends.backend_pdf import PdfPages

file = sys.argv[1]
Nfile = int(sys.argv[2])	

timelist = np.array([])
ratelist = np.array([])
countlist = np.array([])

#pdfname = f'{file}_monitor.pdf'
pdfname = f'/home/icecube603/Chiba/wuBaseV1.4/new_pdf/new.pdf'
pdf = PdfPages(f"{pdfname}")

for n in range(Nfile):
    filename = f'{file}{n}.hd5'
    try:
        f = h5py.File(filename, mode="r")
        NEXTDATA = 0
    except:
        print(filename, 'is not found')
        NEXTDATA = 1

    if NEXTDATA == 0:
        for key in f["data"].keys():
            print("{0}\t{1}".format(key, f["data"][key].dtype))
        #for key in f["metadata"].keys():
        #	print("{0}\t{1}".format(key, f["metadata"][key].dtype))


        nsample          = f["data"]["nsample"][()]
        FPGAtime         = f["data"]["FPGAtime"][()]
        FPGAtcword       = f["data"]["FPGAtcword"][()]
        charge_ch1       = f["data"]["charge_ch1"][()]        
        peak_ch1         = f["data"]["peak_ch1"][()]          
        time_ch1         = f["data"]["time_ch1"][()]          
        charge_fit_ch1   = f["data"]["charge_fit_ch1"][()]    
        peak_fit_ch1     = f["data"]["peak_fit_ch1"][()]      
        time_fit_ch1     = f["data"]["time_fit_ch1"][()]      
        pedestal_ch1     = f["data"]["pedestal_ch1"][()]     
        ADC_ch1          = f["data"]["ADC_ch1"][()]               
        charge_ch2       = f["data"]["charge_ch2"][()]        
        peak_ch2         = f["data"]["peak_ch2"][()]          
        time_ch2         = f["data"]["time_ch2"][()]          
        charge_fit_ch2   = f["data"]["charge_fit_ch2"][()]    
        peak_fit_ch2     = f["data"]["peak_fit_ch2"][()]      
        time_fit_ch2     = f["data"]["time_fit_ch2"][()]      
        pedestal_ch2     = f["data"]["pedestal_ch2"][()]     
        ADC_ch2          = f["data"]["ADC_ch2"][()]               
        Nwfm             = f["metadata"]["Nwfm"][()]
        Nsample_max      = f["metadata"]["Nsample_max"][()]          
        #voltage          = f["metadata"]["voltage"][()]       
        DACvalue         = f["metadata"]["DACvalue"][()]       
        conversion_ch1   = f["metadata"]["conversion_ch1"][()]
        conversion_ch2   = f["metadata"]["conversion_ch2"][()]
        date             = f["metadata"]["date"][()]          
        PMTID            = f["metadata"]["PMTID"][()]   
        PMTIDstr         = f["metadata"]["PMTIDstr"][()]
        wubaseID         = f["metadata"]["wubaseID"][()]      
        runtype          = f["metadata"]["runtype"][()]       
        #OSCfreq          = f["metadata"]["OSCfreq"][()]       
        creatorname      = f["metadata"]["creatorname"][()]#.decode("utf-8")       
        description      = f["metadata"]["description"][()]#.decode("utf-8")   
        fitversion       = f["metadata"]["fitversion"][()]   
        temperature      = f["metadata"]["temperature"][()]

        #description = 9877
        print(Nwfm,'data')
        try:
            print(int(description))
            description = int(description)
            if description<50000:
                description = description + 240000
        except:
            descri=0

        def Gauss1(x,Aspe,Vspe,Mspe):
            # Mspe=mean
            gaussian = Aspe *np.exp(-(x-Mspe)**2/(2*Vspe**2))
            return gaussian

        def Darkrate(dt,tau,A):
            entries = A*np.exp(-dt/tau)
            return entries
        
        bins_set = 100
        range_set = (0,200)

        hist, bins = np.histogram(peak_ch1, bins = bins_set,range = range_set)
        x_list = []
        for i in range(len(hist)):
            x_list.append((bins[i] + bins[i+1]) / 2)

        x_data = np.array(x_list)
        y_data = hist
        limit = int(bins_set*0.15)
        limit2 = int(bins_set*0.85)
        mean = np.mean(peak_ch1[limit:])
        mean = 90

        try:
            popt, pcov = curve_fit(Gauss1,x_data[limit:limit2],y_data[limit:limit2],p0=[1000,20,mean])
            ADCspe = popt[2]
            print(popt)
            print('SPE  ',popt[2],'counts')
            print('0.2PE',popt[2]*0.2,'counts')
            print('0.25PE',popt[2]*0.25,'counts')
        except RuntimeError:
            popt=[1,1,1]
            pcov= [[1, 1, 1],[1, 1, 1],[1, 1, 1]]
            print( 'Optimal parameters not found')
        '''
        plt.hist(peak_ch1,bins = bins_set,log = True,range = range_set,histtype="step")
        plt.plot(x_data[limit:limit2],Gauss1(x_data[limit:limit2], popt[0],popt[1],popt[2]),'-')   #result line
        plt.vlines(popt[2]*0.2,0.8,100000)
        plt.vlines(popt[2]*0.25,0.8,100000)
        plt.ylim(0.8,np.max(hist)*2)
        plt.title("Peak hight distribution")
        plt.show()
        '''
        nwf = Nwfm-1
        #nwf = 50000
        threshold = 0.25
        deltat = np.array([])
        peak_dark = np.array([])
        #peak_dark=np.append(peak_dark,1)
        timestamp = FPGAtime[0]   # *6e7[s]

        print(nwf,'calcuration')

        for i in range(nwf-11):
            if ADC_ch1[i][37]>=int(pedestal_ch1[i])-2:  # multi hit
                sigindex = np.where(ADC_ch1[i][:]>=ADCspe*threshold+pedestal_ch1[i]) 
                sigindex = np.asarray(sigindex)
                #plt.plot(ADC_ch1[i][:50])
                #plt.scatter(sigindex,ADC_ch1[i][sigindex])
                
                for j in range(len(sigindex[0,:])):
                    if int(ADC_ch1[i][sigindex[0,j]])-int(ADC_ch1[i][sigindex[0,j]-1])>0:
                        if int(ADC_ch1[i][sigindex[0,j]])-int(ADC_ch1[i][sigindex[0,j]+1])>0:
                            deltat = np.append(deltat,(((FPGAtime[i]+sigindex[0,j])-timestamp)/6*1e-7))
                            
                            #print('deltat',((FPGAtime[i]+sigindex[0,j])-timestamp)/6*1e-7)
                            #print('index',sigindex[0,j],' t[s]',sigindex[0,j]/6e7,' Tfpga[s]',FPGAtime[i]/6*1e-7+sigindex[0,j]/6e7)
                            timestamp = FPGAtime[i]+sigindex[0,j]
            
            elif peak_ch1[i] >= ADCspe*threshold:        # single fit
                deltat = np.append(deltat,(((FPGAtime[i]+time_fit_ch1[i])-timestamp)/6*1e-7))
                
                #print((((FPGAtime[i]+time_fit_ch1[i])-timestamp)/6*1e-7))
                timestamp = FPGAtime[i]+time_fit_ch1[i] # *6e7[s]
        bins_set = 100
        range_set = (0,0.01)

        hist, bins = np.histogram(deltat, bins = bins_set,range = range_set)
        x_list = []
        for i in range(len(hist)):
            x_list.append((bins[i] + bins[i+1]) / 2)

        x_data = np.array(x_list)
        y_data = hist

        start = 30
        popt, pcov = curve_fit(Darkrate,x_data[start:], y_data[start:], p0=[0.001,1000])
        print(popt)
        print('threshold {} PE'.format(threshold))
        print('darkrate {:.5g} Hz'.format(1/popt[0]))
        '''
        plt.hist(deltat,bins = bins_set,range = range_set,histtype="step")
        plt.plot(x_data[start:],Darkrate(x_data[start:], popt[0],popt[1]),label='tau={:.5g}'.format(popt[0]))
        
        plt.title(f' {filename}', fontsize=14)
        plt.ylabel('entries')
        plt.xlabel('deltat [s]')
        plt.yscale('log')
        plt.legend()
        plt.ylim(2,np.max(hist))
        plt.show()

        plt.hist(np.log10(deltat),bins=bins_set,log = False,range = (-8,2),histtype="step")
        plt.ylabel('entries')
        plt.xlabel('log10(deltat) [s]')
        plt.show()
        '''
        nwf = Nwfm-10
        #nwf = 50000
        threshold = 0.25
        deltat = np.array([])
        deltat = np.array([])
        timestamp = FPGAtime[0]   # *6e7[s]

        print(nwf,'calcuration')

        for i in range(nwf):
            if ADC_ch1[i][35]>=int(pedestal_ch1[i])-2:
                sigindex = np.where(ADC_ch1[i][:]>=ADCspe*threshold+pedestal_ch1[i]) 
                sigindex = np.asarray(sigindex)
                #plt.plot(ADC_ch1[i][:50])
                #plt.scatter(sigindex,ADC_ch1[i][sigindex])
                
                for j in range(len(sigindex[0,:])):
                    #print(i)
                    if int(ADC_ch1[i][sigindex[0,j]])-int(ADC_ch1[i][sigindex[0,j]-1])>0:
                        if int(ADC_ch1[i][sigindex[0,j]])-int(ADC_ch1[i][sigindex[0,j]+1])>0:
                            deltat = np.append(deltat,(((FPGAtime[i]+sigindex[0,j])-timestamp)/6*1e-7))
                            #print(deltat[-1])
                            peak_dark = np.append(peak_dark,ADC_ch1[i][sigindex[0,j]]-pedestal_ch1[i])
                            #print('deltat',((FPGAtime[i]+sigindex[0,j])-timestamp)/6*1e-7)
                            #print('index',sigindex[0,j],' t[s]',sigindex[0,j]/6e7,' Tfpga[s]',FPGAtime[i]/6*1e-7+sigindex[0,j]/6e7)
                            timestamp = FPGAtime[i]+sigindex[0,j]
                            
            
            elif peak_ch1[i] >= ADCspe*threshold:     
                deltat = np.append(deltat,(((FPGAtime[i]+time_fit_ch1[i])-timestamp)/6*1e-7))
                peak_dark = np.append(peak_dark,peak_ch1[i])
                #print((((FPGAtime[i]+time_fit_ch1[i])-timestamp)/6*1e-7))
                timestamp = FPGAtime[i]+time_fit_ch1[i] # *6e7[s]
        print(len(deltat[np.where(deltat<1)]),sum(deltat[np.where(deltat<1)]))
        print((sum(deltat[np.where(deltat<1)])/len(deltat[np.where(deltat<1)]))**(-1),'Hz')
        darkcount = (sum(deltat[np.where(deltat<1)])/len(deltat[np.where(deltat<1)]))**(-1)
        
        timelist = np.append(timelist,int(description))
        ratelist = np.append(ratelist,1/popt[0])
        countlist = np.append(countlist,darkcount)

meanrate = np.mean(countlist)
print('mean of countrate: {:.5g} Hz'.format(np.mean(countlist[np.where(countlist>=meanrate-10)])))
print('sigma of countrate: {:.5g} Hz'.format(np.std(countlist[np.where(countlist>=meanrate-10)])))
fig, ax = plt.subplots()
#plt.annotate('mean of countrate: {:.5g} Hz'.format(np.mean(countlist[np.where(countlist>=meanrate-10)])),xy=(0.5,2e3), fontsize=14)
plt.scatter(timelist,countlist)
plt.ylim(np.mean(countlist[np.where(countlist>=meanrate-10)])-20,np.mean(countlist[np.where(countlist>=meanrate-10)])+20)
plt.title(f'BB{PMTID}  '+'mean of countrate: {:.5g} Hz'.format(np.mean(countlist[np.where(countlist>=meanrate-10)])),fontsize = 15)
plt.xlabel('Time',fontsize = 15)
plt.ylabel('Dark Rate[Hz]', fontsize = 15)
#ax.set_aspect()
fig.tight_layout()
pdf.savefig(fig)
plt.show()

pdf.close()
'''
import matplotlib.cm as cm
import matplotlib.colors

fig = plt.figure()
ax = fig.add_subplot(111)

x = np.log10(deltat)
y = peak_dark[:len(x)]*2.35*conversion_ch1

H = ax.hist2d(x,y, bins=[np.linspace(-8,0,100),np.linspace(0,25,100)], norm=matplotlib.colors.LogNorm(), cmap=cm.jet)
ax.set_title(f'BB{PMTID}')
ax.set_xlabel('log10(deltat)')
ax.set_ylabel('Charge')
fig.colorbar(H[3],ax=ax)
plt.grid()
plt.show()

Z, X, Y = np.histogram2d(x,y, bins=[np.linspace(-8,0,101),np.linspace(0,25,101)])
bins_set = 100
hist, bins = np.histogram(np.log10(deltat),bins=bins_set,range = (-8,0))
h1d = hist[:,np.newaxis]
newZ = Z/h1d

plt.pcolormesh(X,Y,np.transpose(newZ), norm=matplotlib.colors.LogNorm(),cmap="jet")
plt.title(f'BB{PMTID}  2Dhist/1Dhist')
plt.xlabel('log10(deltat)')
plt.ylabel('Charge')
plt.colorbar()
plt.grid()
plt.show()
'''
