from distutils.command.build_scripts import first_line_re
import numpy as np
import matplotlib.pyplot as plt
import sys 
from scipy import integrate
from scipy.optimize import curve_fit

WAIT_FOR_HIT=0
NSAMPLES_NEXT=1
TIMESTAMP_NEXT=2
TDCWORD_NEXT=3
CH1_NEXT=4
CH1_IN_PROGRESS=5
CH2_IN_PROGRESS=6
readstate=WAIT_FOR_HIT

hitData = []	

filename = sys.argv[1]

f = open(filename, "r")

readingCh1=False   # will be True when line is expected to have ch1 data
readingCh2=True   # similar for ch2
#waveformCount = 0from scipy import integrate
ch1charge = np.array([])
ch2charge = np.array([])
deltach1_array = np.array([])
deltach2_array = np.array([])
timelist = np.array([])
max_Ch1 = np.array([])
max_Ch2 = np.array([])
onedatatime = 15*10**(-9)   #[s]  15ns
count = 0
counts = 0
skipped = 0
conv = 2.65 #micto to wu

##### setting of histgram
bins_set = 80
range_set = (0,100)#pC
range_set2 = (0,100)#ADC

fig, axs = plt.subplots(1,3, figsize=(18,3.7), constrained_layout=True)

for line in f:
	linestrip=line.rstrip()    #文字列.rstrip(除去する文字)
	tokens=linestrip.split()    #空白文字で分割
	if readstate==WAIT_FOR_HIT:
		if len(tokens)!=1:continue
		if tokens[0]!="V1" and tokens[0]!="V1X": continue
		readstate=NSAMPLES_NEXT
		continue
	elif readstate==NSAMPLES_NEXT:
		if len(tokens)!=1:
			readstate=WAIT_FOR_HIT #unexpected input
		else:
			nsamples=int(tokens[0],16)
			readstate=TIMESTAMP_NEXT
		if nsamples>256:
			readstate=WAIT_FOR_HIT # unexpected input
		continue
	elif readstate==TIMESTAMP_NEXT:
		if len(tokens)!=1:
			readstate=WAIT_FOR_HIT #unexpected input
		else:
			time=int(tokens[0],16)
			if time < 2923270018009:
				timelist=np.append(timelist,time)
			readstate=TDCWORD_NEXT
		continue

	elif readstate==TDCWORD_NEXT:

		if len(tokens)!=1:

			readstate=WAIT_FOR_HIT #unexpected input
			readstate=WAIT_FOR_HIT #unexpected input
		else:
			readstate=CH1_NEXT
		continue

	elif readstate==CH1_NEXT:
		ch1=[]
		ch2=[]
		discraw=[]
		discsync=[]
		for hexvalue in tokens:
			raw=int(hexvalue,16)&0xfff
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
			#count += 1
		#ch1long=np.append(ch1long,ch1)
		if len(ch1)==nsamples:
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
			readstate=WAIT_FOR_HIT # that's next after processing data
			
			ymax=max(ch1)
			if ymax==0: continue   # empty waveform = not good
			
			if(len(ch2)>30 ):  ################ threshold ##############
				Ch1 = np.asarray(ch1)
				Ch2 = np.asarray(ch2)
				count += 1
				#axs[0].plot(Ch1)
				if np.min(Ch2[5:10])>3770:
					baseline1 = np.mean(Ch1[5:10])
					baseline2 = np.mean(Ch2[5:10])
					Ch1 = Ch1 - baseline1
					Ch2 = Ch2 - baseline2
					# axs[0].plot(Ch1[:30])
				elif np.min(Ch2[10:15])>3770:
					baseline1 = np.mean(Ch1[10:15])
					baseline2 = np.mean(Ch2[10:15])
					Ch1 = Ch1 - baseline1
					Ch2 = Ch2 - baseline2
					# baseline = np.mean(Ch1[10:15])
					# Ch1 = Ch1 - baseline
					# axs[0].plot(Ch1[:30])
				else:
					baseline1 = np.mean(Ch1[5:15])
					baseline2 = np.mean(Ch2[5:15])
					Ch1 = Ch1 - baseline1
					Ch2 = Ch2 - baseline2
					# baseline = np.mean(Ch1[5:15])
					# Ch1 = Ch1 - baseline
					# axs[0].plot(Ch1[:30])
				#print(baseline1)
				#print(baseline2)
				if np.max(Ch2)> 100:
					skipped += 1
				elif np.min(Ch2)<-4000:
					skipped += 1
				else:
					if counts<100:
						axs[0].plot(Ch2)
					counts += 1
					# if np.min(Ch1)<-5:
					# 	skipped+=1
					max_ch2 = np.min(Ch2)
					max_ch1 = np.max(Ch1)
					max_Ch2 = np.append(max_Ch2,max_ch2)
					max_Ch1 = np.append(max_Ch1,max_ch1)
					max2 = ch2.index(min(ch2))
					max1 = ch1.index(min(ch1))
					Ch2_masked1=Ch2[max2-4 : max2+5]
					Ch1_masked1=Ch1[max1-4 : max1+5]
					if len(Ch1_masked1) == 9:
						if len(Ch2_masked1) == 9:
							Ch1_masked = Ch1_masked1*(-1)
							Ch2_masked = Ch2_masked1*(-1)
							# axs[0].plot(Ch1_masked)
						
							# onedatatime = 15*10**(-9)   #[s]  15ns
							voltage1 = Ch1_masked * 0.44 * 10**(-3) #[V]  0.44[mV]
							voltage2 = Ch2_masked * 0.44 * 10**(-3) #[V]  0.44[mV]
							current1 = voltage1 / (3600*1) * conv#[A]
							current2 = voltage2 / (3600*0.015) * conv#[A]
							# max_Ch2i = np.min(current2)
							# max_Ch2 = np.append(max_Ch2,max_Ch2i)
							integrated_ch1 = np.sum(current1)
							integrated_ch2 = np.sum(current2)
							charge1 = integrated_ch1 * onedatatime *10**12  #[pC]
							charge2 = integrated_ch2 * onedatatime *10**12  #[pC]
							ch1charge = np.append(ch1charge, charge1)
							ch2charge = np.append(ch2charge, charge2)
							
							#counts += 1
					else:
						skipped +=1
				
		else:
			readstate=CH2_IN_PROGRESS
		continue


max_Ch2i = -max_Ch2 #* -0.44 * 10**(-3) / (3600*0.015) * conv *10**3 #[uA]
################################## fit #######################################

hist, bins = np.histogram(ch2charge, bins = bins_set, range=range_set)
hist2, bins2 = np.histogram(max_Ch2i, bins = bins_set, range=range_set2)
x_list = []
x_list2 = []

for i in range(len(hist)):
    x_list.append((bins[i] + bins[i+1]) / 2)

for i in range(len(hist2)):
    x_list2.append((bins2[i] + bins2[i+1]) / 2)

x_data = np.array(x_list)
y_data = hist
x_data2 = np.array(x_list2)
y_data2 = hist2

def Gauss1(x,Aspe,Vspe,Mspe):
	# Mspe=mean
	gaussian = Aspe *np.exp(-(x-Mspe)**2/(2*Vspe**2))
	return gaussian

def Gauss2(x,Aped,Vped,Mped,Aspe,Vspe,Mspe):
    doublegaussian = Aped *np.exp(-(x-Mped)**2/(2*Vped**2)) + Aspe *np.exp(-(x-Mspe)**2/(2*Vspe**2))
    return doublegaussian
limit = 0
err = np.where(y_data[limit:]==0, 1000, np.sqrt(y_data[limit:]))
err2 = np.where(y_data2[limit:]==0, 1000, np.sqrt(y_data2[limit:]))

ini1 = np.mean(ch2charge)  #initial value
ini2 = np.mean(max_Ch2i)

try:
	popt, pcov = curve_fit(Gauss1,x_data[limit:],y_data[limit:], sigma=err,p0=[0,1,ini1])
except:
	popt=[1,1,1]
	pcov= [[1, 1, 1],[1, 1, 1],[1, 1, 1]]
	print( 'Optimal parameters not found')
try:
	popt2, pcov2 = curve_fit(Gauss1,x_data2[limit:],y_data2[limit:], sigma=err2,p0=[0,1,ini2])
except RuntimeError:
	popt2=[1,1,1]
	pcov2= [[1, 1, 1],[1, 1, 1],[1, 1, 1]]
	print( 'Optimal parameters not found')


# axs[1].plot(x_data[limit],y_data[limit],'o',color='k')
#########################################################################

# cut=ch2charge[np.where(ch2charge<range_fit[1])]
# fit_charge = cut[np.where(cut>range_fit[0])]

axs[1].hist(ch2charge,bins = bins_set,log = True,range = range_set)
axs[2].hist(max_Ch2i,bins = bins_set,log = True,range = range_set2)
###########################################################################
# axs[1].plot(x_data,Gauss2(x_data, popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]),'-')   #result line
axs[1].plot(x_data,Gauss1(x_data, popt[0],popt[1],popt[2]),'-')   #result line
axs[2].plot(x_data2,Gauss1(x_data2, popt2[0],popt2[1],popt2[2]),'-')   #result line
# axs[2].plot(x_data,Gauss2(x_data, popt2[0],popt2[1],popt2[2],popt2[3],popt2[4],popt2[5]),'-')
# axs[1].annotate(""+str(popt[5])+"", xy=(0.05,0.7), xycoords='axes fraction')
axs[1].annotate(""+str(popt[2])+"", xy=(0.1,0.9), xycoords='axes fraction')
axs[2].annotate(""+str(popt2[2])+"", xy=(0.1,0.9), xycoords='axes fraction')
axs[1].set_ylim(0.8,5000)

axs[1].set_ylim(0.8,np.max(hist)*2)
axs[2].set_ylim(0.8,np.max(hist2)*2)
# axs[2].set_ylim(0.8,20000)
# gain1 = popt[5]/(1.6*10**(-19)*10**12)
gain2 = popt[2]/(1.6*10**(-19)*10**12)
# print("gain is",gain1,'or',gain2)
print(gain2)

###########################################################################
axs[0].set_ylabel('ADC[counts]', fontsize=16)
# axs[0].set_xlabel('time[ns]')
axs[0].set_xlabel('bins', fontsize=16)
axs[1].set_xlabel('charge [pC]', fontsize=16)
# axs[2].set_xlabel('charge/bin [fC/s]', fontsize=16)
axs[2].set_xlabel('ADCpeak [counts]', fontsize=16)
axs[0].set_title('ch2(Dynode8)', fontsize=16)
# axs[1].set_title('charge', fontsize=16)
# axs[2].set_title('Peak Current', fontsize=16)
print(count,'waves',counts)
print(skipped,'skipped')
print('popt',popt)
print('pcov',pcov)
print('popt_current',popt2)
print('pcov_current',pcov2)
axs[0].tick_params(labelsize=18)
axs[1].tick_params(labelsize=18)
axs[2].tick_params(labelsize=18)

axs[2].grid()
plt.show()
