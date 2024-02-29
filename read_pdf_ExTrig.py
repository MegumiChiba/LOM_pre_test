
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

filename = sys.argv[1]	
f = h5py.File(filename, mode="r")

for key in f["data"].keys():
	print("{0}\t{1}".format(key, f["data"][key].dtype))
#for key in f["metadata"].keys():
#	print("{0}\t{1}".format(key, f["metadata"][key].dtype))
pdfname = f'{filename}_results.pdf'
pdf = PdfPages(f"{pdfname}")

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



print(Nwfm,'data')
#plt.hist(time_fit_ch1,bins=100)
#plt.show()
#plt.hist(time_fit_ch2,bins=100)
#plt.show()

def Gauss1(x,Aspe,Vspe,Mspe):
	# Mspe=mean
	gaussian = Aspe *np.exp(-(x-Mspe)**2/(2*Vspe**2))
	return gaussian
bins_set = 80

spetime = np.array([])
spetime2 = np.array([])
for i in range(len(charge_fit_ch1)):
    if charge_fit_ch1[i]*conversion_ch1>=0.4 and charge_fit_ch1[i]*conversion_ch1<=1.5:
        #print('spe',i)
        if time_fit_ch1[i]<=26 and time_fit_ch1[i]>=15:
            spetime = np.append(spetime,time_fit_ch1[i])
            spetime2 = np.append(spetime2, time_fit_ch2[i])
#print(spetime[:10],spetime2[:10],(spetime[:10]+spetime2[:10])*1/60*1000)

d = (spetime+spetime2)*1/60*1000
range_set = (np.mean(d)-20,np.mean(d)+20)
hist, bins = np.histogram(d, bins = bins_set,range =range_set)
x_list = []
for i in range(len(hist)):
    x_list.append((bins[i] + bins[i+1]) / 2)

x_data = np.array(x_list)
y_data = hist
limit = int(bins_set*0.1)
mean = np.mean(d)

try:
    popt, pcov = curve_fit(Gauss1,x_data[limit:],y_data[limit:],p0=[100,3,mean])
    print(popt)
except RuntimeError:
    popt=[1,1,1]
    pcov= [[1, 1, 1],[1, 1, 1],[1, 1, 1]]
    print( 'Optimal parameters not found')
print('Time resolution',popt[1])#*1/60*1000)
print(len(spetime),'wf was counted')
#print('mean',mean,'max',max)
fig, ax = plt.subplots()
ax.hist(d,bins = bins_set,range =range_set)#,log = True)
ax.plot(x_data,Gauss1(x_data, popt[0],popt[1],popt[2]),'-')   #result line
#plt.ylim(0.8,np.max(hist)*2)
ax.set_title(f'Time distribution BB{PMTID}', fontsize=16)
ax.set_xlabel('Tch1 - Tch2 [ns]', fontsize=16)
ax.set_ylabel('entries', fontsize=16)
plt.annotate('Time resolution {:.4g} ns'.format(popt[1]),xy=(popt[2]+2,popt[0]-20), fontsize=14)
plt.tick_params(labelsize=14)
#plt.xlim(popt[2]-40,popt[2]+40)
fig.tight_layout()
pdf.savefig(fig)
plt.show()


binchg = np.linspace(-0.5#int(np.min(charge_ch1*conversion_ch1))
					, int(np.max(charge_ch1*conversion_ch1)), 100)
bins_chargeplot = 100
range_chargeplot = (0,3)
fig, ax = plt.subplots()
#plt.title("Voltage = {0:.1f} V, RunType={1}, PMT={2}, Nwfm={3}".format(voltage, runtype, PMTID, Nwfm))
ax.hist(charge_ch1*conversion_ch1, bins=bins_chargeplot, log=True,range=range_chargeplot, histtype="step",label="charge")
ax.hist(charge_fit_ch1*conversion_ch1, bins=bins_chargeplot,range=range_chargeplot, histtype="step",label="charge_fit")
ax.set_title("Charge distribution")
plt.legend()
fig.tight_layout()
pdf.savefig(fig)
plt.show()

fig, ax = plt.subplots()
for i in range(100):
    ax.plot(ADC_ch1[i][:30])
ax.set_title('CH1 100 waves')
ax.set_xlabel('bins')
ax.set_ylabel('ADC')
fig.tight_layout()
pdf.savefig(fig)
plt.show()

fig, ax = plt.subplots()
for i in range(100):
    ax.plot(ADC_ch2[i][:30])
ax.set_title('CH2 100 waves')
ax.set_xlabel('bins')
ax.set_ylabel('ADC')
#plt.grid()
fig.tight_layout()
pdf.savefig(fig)
plt.show()

def FitCharge(q,q0,q1,N0,Sped,Sspe,lam,eta,tau):
	# Sspe:Sigma_spe
    fit_ped = N0/(np.sqrt(2*np.pi)*Sped) * np.exp(-(q-q0)**2/(2*Sped**2))
    fit_exp = N0*eta/tau * np.exp(-q/tau) * 1*(q>Sped)
    fit_spe = N0*lam*np.exp(-lam)/(np.sqrt(2*np.pi)*Sspe) * np.exp(-(q-q1-q0)**2/(2*Sspe**2))
    fit_2pe = N0*lam**2*np.exp(-lam)/(2*np.sqrt(2*2*np.pi)*Sspe) * np.exp(-(q-2*q1-q0)**2/(2*2*Sspe**2))
    fit = fit_ped + fit_exp + fit_spe + fit_2pe
    return fit

def Fitipe(q,q0,q1,N0,Sped,Sspe,lam,i,tau):
	# Sspe:Sigma_spe
    fit_2pe = N0*lam**2*np.exp(-lam)/(i*np.sqrt(i*2*np.pi)*Sspe) * np.exp(-(q-i*q1-q0)**2/(i*2*Sspe**2))
    return fit_2pe

def FitPed(q,q0,q1,N0,Sped,Sspe,lam,eta,tau):
	# Sspe:Sigma_spe
    fit_ped = N0/(np.sqrt(2*np.pi)*Sped) * np.exp(-(q-q0)**2/(2*Sped**2))
    return fit_ped

def FitSpe(q,q0,q1,N0,Sped,Sspe,lam,eta,tau):
	# Sspe:Sigma_spe
    fit_spe = N0*lam*np.exp(-lam)/(np.sqrt(2*np.pi)*Sspe) * np.exp(-(q-q1-q0)**2/(2*Sspe**2))
    return fit_spe

def FitExp(q,q0,q1,N0,Sped,Sspe,lam,eta,tau):
	# Sspe:Sigma_spe
    fit_exp = N0*eta/tau * np.exp(-q/tau)* 1*(q>0)#ped)
    return fit_exp

def Fit2pe(q,q0,q1,N0,Sped,Sspe,lam,eta,tau):
	# Sspe:Sigma_spe
    fit_2pe = N0*lam**2*np.exp(-lam)/(2*np.sqrt(2*2*np.pi)*Sspe) * np.exp(-(q-2*q1-q0)**2/(2*2*Sspe**2))
    return fit_2pe

bins_set = 60 #Changed 20240227
range_set = (-0.08,2.2)
hist, bins = np.histogram(charge_fit_ch1*conversion_ch1, bins = bins_set, range=range_set)
x_list = []

for i in range(len(hist)):
    x_list.append((bins[i] + bins[i+1]) / 2)

x_data = np.array(x_list)
y_data = hist
y_data[np.where(y_data==0)] = 1

limit_min = int(len(y_data)*0)
limit_max = int(len(y_data)*1)
sigma = np.sqrt(y_data[limit_min:limit_max])

p0_q1 = x_data[np.where(y_data==np.max(y_data[np.where(x_data>0.2)]))]
#print(p0_q1,p0_Sspe)
p0_N0 = np.max(y_data)/10
p0_Sped = x_data[np.where(y_data == np.min(y_data[np.where(x_data<p0_q1[0])]))] /5
p0_Sspe = x_data[np.where(x_data ==p0_q1[0])] /2.5
print(p0_Sped,p0_Sspe)
p0=[0, p0_q1[0],p0_N0,abs(p0_Sped[0]),p0_Sspe[0], 0.5,0.03,0.3 ]
#  (q0, q1,       N0,Sped,Sspe,lam(spe counts),eta(hight),tau(width))

try:
	popt, pcov = curve_fit(FitCharge,x_data[limit_min:limit_max],y_data[limit_min:limit_max],p0=p0,sigma=sigma)
	sigma_param = np.sqrt([pcov[0,0],pcov[1,1],pcov[2,2]])
except RuntimeError:
    popt=[1,1,1,1,1,1,1,1]
    pcov= [[1, 1, 1,1, 1, 1,1],[1, 1, 1,1, 1, 1,1],[1, 1, 1,1, 1, 1,1],[1, 1, 1,1, 1, 1,1]
    ,[1, 1, 1,1, 1, 1,1],[1, 1, 1,1, 1, 1,1],[1, 1, 1,1, 1, 1,1]]
    print('Optimal parameters not found')

x_fitplot = np.linspace(x_data[limit_min],x_data[limit_max-1],10000)
y_fitplot = FitCharge(x_fitplot, popt[0],popt[1],popt[2], popt[3],popt[4],popt[5], popt[6],popt[7])
valley = np.min(y_fitplot[np.min(np.where(x_fitplot>0)):np.min(np.where(x_fitplot>popt[1]))])
peak = np.max(y_fitplot[np.min(np.where(x_fitplot>popt[1]-popt[4])):np.min(np.where(x_fitplot>popt[1]+popt[4]))])
PtoV = peak/valley    #Vに対する比

print('popt',popt)
print('resolution',popt[4]/popt[1])

fig, ax = plt.subplots()
plt.hist(charge_fit_ch1*conversion_ch1,bins = bins_set,log = True,histtype="step",range = range_set,label='data')

plt.plot(x_fitplot,
            FitCharge(x_fitplot, popt[0],popt[1],popt[2], popt[3],popt[4],popt[5], popt[6],popt[7]),
                    '-',label='fit',c='b')   #result line
plt.plot(x_fitplot,
            FitPed(x_fitplot, 
                popt[0],popt[1],popt[2], popt[3],popt[4],popt[5], popt[6],popt[7]),
                    '--',c='gray')   #result line
plt.plot(x_data[limit_min:limit_max],
            FitSpe(x_data[limit_min:limit_max], 
                popt[0],popt[1],popt[2], popt[3],popt[4],popt[5], popt[6],popt[7]),
                '--',c='gray')   #result line
plt.plot(x_data[limit_min:limit_max],
            FitExp(x_data[limit_min:limit_max], 
                popt[0],popt[1],popt[2], popt[3],popt[4],popt[5], popt[6],popt[7]),
                '--',c='gray')   #result line
plt.plot(x_data[limit_min:limit_max],
            Fit2pe(x_data[limit_min:limit_max], 
                popt[0],popt[1],popt[2], popt[3],popt[4],popt[5], popt[6],popt[7]),
                '--',c='gray')   #result line
plt.annotate('Sigma_ped : '+f"{popt[3]:.3f}"+'pC',xy=(0.5,2e3), fontsize=14)
plt.annotate('Mean_spe  : '+f"{popt[1]:.3f}"+'pC',xy=(0.5,1e3), fontsize=14)
plt.annotate('Resolution : '+f"{popt[4]/popt[1]*100:.1f}"+'%',xy=(0.5,5e2), fontsize=14)
plt.annotate('P/V ratio   : '+f"{PtoV:.2f}"+'',xy=(0.5,2.5e2), fontsize=14)

plt.ylim(0.8,20000)
plt.ylim(0.8,np.max(hist)*2)
plt.ylabel('entries', fontsize=16)
plt.xlabel('charge [pC]', fontsize=16)
plt.title(f'BB{PMTID}', fontsize=16)
#plt.title(f'{basename}_mod {IVconv} Ohm', fontsize=20)
plt.tick_params(labelsize=14)
plt.legend(ncol=2,loc='upper center')
plt.grid()
ax.set_aspect(0.4)
fig.tight_layout()
pdf.savefig(fig)
plt.show()

pdf.close()
