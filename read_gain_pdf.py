
import h5py 
import numpy as np
import datetime
import time
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit  
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages     
import sys 

file    = sys.argv[1]
pdfname = sys.argv[2]
PMTID = sys.argv[3]
vlist = [80,82,84,86,88,90,92,94,96,98]
vlist = [76,78,80,82,84,86,88,90,92,94]
#vlist = [74,76,78,80,82,84,86,88,90,92]
#vlist = [90,92,94,96,98,100,102,104,106,108]
#vlist = [85]
blank = []
bins_set = 80
pdfname = f'{pdfname}_gainplots.pdf'
#pdf = PdfPages(f"data/lom_test/BB{PMTID}/{pdfname}")
pdf = PdfPages(f"{pdfname}")
#pdf = PdfPages(f"data/BB{PMTID}/{pdfname}")

gain_list = np.array([])

def Gauss1(x,Aspe,Vspe,Mspe):
	# Mspe=mean
	gaussian = Aspe *np.exp(-(x-Mspe)**2/(2*Vspe**2))
	return gaussian
def GainFit(Vct,Again):
    k = 0.75
    n = 10
    gainfit = Again *Vct**(k*n)
    return gainfit

for i in range(len(vlist)):
    filename    = file + str(vlist[i]) + '.hd5'
    try:
        f = h5py.File(filename, mode="r")
    except:
        blank = np.append(blank,i)
        continue

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
   # voltage          = f["metadata"]["voltage"][()]       
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

    print(vlist[i],'V')
    print(Nwfm,'data')

    range_set = (-0, 1.5 + (vlist[i]-vlist[0] )/8)
    lim1 = int(bins_set*0.15)
    lim2 = int(bins_set*0.9)
    
    d = charge_fit_ch1*conversion_ch1
    hist, bins = np.histogram(d, bins = bins_set,range = range_set)
    x_list = []
    for j in range(len(hist)):
        x_list.append((bins[j] + bins[j+1]) / 2)

    x_data = np.array(x_list)
    y_data = hist
    mean = np.mean(d[lim1:lim2])

    try:
        popt, pcov = curve_fit(Gauss1,x_data[lim1:lim2],y_data[lim1:lim2],p0=[100,0.3,mean])
        print(popt)
    except RuntimeError:
        popt=[1,1,1]
        pcov= [[1, 1, 1],[1, 1, 1],[1, 1, 1]]
        print( 'Optimal parameters not found')
    
    gain = popt[2]*10**(-12)/(1.6*10**(-19))
    gain_list = np.append(gain_list,gain)
    print('SPE mean',popt[2])
    print('Gain {:.2g}'.format(gain))

    fig, ax = plt.subplots()
    ax.hist(charge_ch1*conversion_ch1, bins=bins_set, log=True,range = range_set, histtype="step",label="charge")
    ax.hist(charge_fit_ch1*conversion_ch1, bins=bins_set,range = range_set, histtype="step",label="charge_fit")
    ax.plot(x_data[lim1:lim2],Gauss1(x_data[lim1:lim2], popt[0],popt[1],popt[2]),'-',label = 'gaussian')   #result line
    ax.set_ylim(0.8,1000)
    ax.set_title(f'BB{PMTID} Charge distribution {vlist[i]}V')
    plt.legend()
    fig.tight_layout()
    pdf.savefig(fig)
    plt.show()

#pdf.close()
try:
    del vlist[int(blank)]
except:
    nblank=0
x_data = np.asarray(vlist)
y_data = gain_list#*0.82
err = np.sqrt(y_data)#np.where(y_data==0, 1000, np.sqrt(y_data))
x_plot = np.linspace(min(x_data),max(x_data),100)

try:
    popt, pcov = curve_fit(GainFit,x_data,y_data, sigma=err)#,p0=)
    print(popt, pcov)
except RuntimeError:
    popt=[1,1,1]
    pcov= [[1, 1, 1],[1, 1, 1],[1, 1, 1]]
    print( 'Optimal parameters not found')

gaincurve = GainFit(x_plot, popt[0])
deltamin = np.min(abs(gaincurve-5e6))
index5e6 = np.where(abs(gaincurve-5e6)==deltamin)
vc_result = x_plot[index5e6]
print('Control voltage of 5e6 gain is {:.4g} V'.format(vc_result[0]))
print('HV ~ {:.4g} V'.format(vc_result[0]*12))

#plt.title(f'Charge distribution {vlist[i]}V')
fig, ax = plt.subplots()
ax.annotate('Control voltage of 5e6 gain is {:.4g} V'.format(vc_result[0]),xy=(76,5e6), fontsize=14)
ax.scatter(vlist,y_data)
ax.plot(x_plot,gaincurve,'-',label = 'gainfit')   #result line
ax.hlines(5e6,np.min(vlist),np.max(vlist),linestyle = "dotted")
ax.set_title(f'BB9869', fontsize = 15)
ax.set_xlabel('Vcontrl [V]', fontsize = 15)
ax.set_ylabel('Gain', fontsize = 15)
plt.grid()
plt.show()
fig.tight_layout()
#pdf.savefig(f"data/lom_test/BB{PMTID}/{fig}")
pdf.savefig(fig)
pdf.close()
'''
    binpeak = np.linspace(int(np.min(peak_ch1)), int(np.max(peak_ch1)), 100)
    plt.title("Peak distribution")
    plt.hist(peak_ch1, bins=binpeak, histtype="step",label = 'peak')
    plt.hist(peak_fit_ch1, bins=binpeak, histtype="step",label = 'peak_fit')
    plt.legend()
    plt.show()
    '''
