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
onedatatime = 15*10**(-9)   #[s]  15ns
count = 0
counts = 0
skipped = 0
conv = 2.05#2.65 #micto to wu

##### setting of histgram
bins_set =100
range_set = (0,1)
range_set2 = (0,100)
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
			
			if(len(ch1)>30 and len(ch1)<90 and max(ch1)>0):  ####### threshold ######
				Ch1 = np.asarray(ch1)
				count += 1
				# axs[0].plot(Ch1)
				if np.max(Ch1[5:10])<332:
					baseline = np.mean(Ch1[5:10])
					Ch1 = Ch1 - baseline
					# axs[0].plot(Ch1[:30])
				elif np.max(Ch1[10:15])<332:
					baseline = np.mean(Ch1[10:15])
					Ch1 = Ch1 - baseline
					# axs[0].plot(Ch1[:30])
				else:
					baseline = np.mean(Ch1[15:20])
					Ch1 = Ch1 - baseline
					# axs[0].plot(Ch1[:30])

				if np.min(Ch1)<-50:
					skipped+=1
				if np.max(Ch1)>5:    #hit
					if Ch1[2] > 1000:
						skipped+=1
					elif Ch1[4] > 1000:
						skipped+=1
					elif Ch1[29] > 1000:
						skipped+=1
					elif Ch1[28] > 1000:
						skipped+=1
					else:
						if counts <= 100:
							axs[0].plot(Ch1)
						max_ch1 = np.max(Ch1[15:25])
						max_Ch1 = np.append(max_Ch1,max_ch1)
						max1 = 21#ch1.index(max(ch1))
						Ch1_masked1=Ch1[max1-4 : max1+5]
						if len(Ch1_masked1) == 9:
							# axs[0].plot(Ch1_masked1)
							Ch1_masked = Ch1_masked1
						
						voltage = Ch1_masked * 0.44 * 10**(-3) #[V]  0.44[mV]
						current = voltage / 3600 * conv#[A]
						integrated_ch1 = np.sum(current)
						charge = integrated_ch1 * onedatatime *10**12  #[pC]
						ch1charge = np.append(ch1charge, charge)

						# adcsum = np.sum(Ch1_masked)/1000          ###########adc
						# ch1charge = np.append(ch1charge, adcsum)########adc
						counts += 1
						
				else:      #baseline
					Ch1_masked=Ch1[10-4 : 10+5]
					# axs[0].plot(Ch1_masked)

					voltage = Ch1_masked * 0.44 * 10**(-3) #[V]  0.44[mV]
					current = voltage / 3600 *conv#[A]
					integrated_ch1 = np.sum(current)
					charge = integrated_ch1 * onedatatime *10**12  #[pC]
					ch1charge = np.append(ch1charge, charge)

					# adcsum = np.sum(Ch1_masked)/1000         ############adc
					# ch1charge = np.append(ch1charge, adcsum)########adc
					counts += 1

		else:
			readstate=CH2_IN_PROGRESS
		continue

max_Ch1i = max_Ch1 #* 0.44 * 10**(-3) / (3600) * conv*10**12   #[pC/bin]
##################################added#######################################

hist, bins = np.histogram(ch1charge, bins = bins_set, range=range_set)
x_list = []
hist2, bins2 = np.histogram(max_Ch1i, bins = bins_set, range=range_set2)
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

limit = 10
limite = bins_set
ini1 = np.mean(ch1charge)  #initial value
ini2 = np.mean(max_Ch1i)

try:
	popt, pcov = curve_fit(Gauss1,x_data[limit:],y_data[limit:],p0=[0,1,ini1])
except RuntimeError:
	popt=[1,1,1]
	pcov= [[1, 1, 1],[1, 1, 1],[1, 1, 1]]
	print( 'Optimal parameters not found')
# popt, pcov = curve_fit(Gauss1,x_data[limit:limite],y_data[limit:limite])
try:
	popt2, pcov2 = curve_fit(Gauss1,x_data2[limit:],y_data2[limit:],p0=[0,1,ini2])
except RuntimeError:
	popt2=[1,1,1]
	pcov2= [[1, 1, 1],[1, 1, 1],[1, 1, 1]]
	print( 'Optimal parameters not found')
# axs[1].plot(x_data[limit],y_data[limit],'o',color='k')
#########################################################################

# cut=ch1charge[np.where(ch1charge<range_fit[1])]
# fit_charge = cut[np.where(cut>range_fit[0])]

axs[1].hist(ch1charge,bins = bins_set,log = True,range = range_set)
axs[2].hist(max_Ch1i,bins = bins_set,log = True,range = range_set2)
###########################################################################
# axs[1].plot(x_data,Gauss2(x_data, popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]),'-')   #result line
axs[1].plot(x_data,Gauss1(x_data, popt[0],popt[1],popt[2]),'-')   #result line
axs[2].plot(x_data2,Gauss1(x_data2, popt2[0],popt2[1],popt2[2]),'-')   #result line
# axs[1].annotate(""+str(popt[5])+"", xy=(0.05,0.7), xycoords='axes fraction')
axs[1].annotate(""+str(popt[2])+"", xy=(0.1,0.9), xycoords='axes fraction')
axs[2].annotate(""+str(popt2[2])+"", xy=(0.1,0.9), xycoords='axes fraction')
axs[1].set_ylim(0.8,np.max(hist)*2)
# axs[1].set_ylim(0.8,10000)
# axs[2].set_ylim(0.8,1000)
axs[2].set_ylim(0.8,np.max(hist2)*2)
# gain1 = popt[5]/(1.6*10**(-19)*10**12)
gain2 = popt[2]/(1.6*10**(-19)*10**12)
# print("gain is",gain1,'or',gain2)
print(gain2)

###########################################################################
axs[0].set_ylabel('ADC[counts]', fontsize=16)
# axs[0].set_xlabel('time[ns]')
axs[0].set_xlabel('bins', fontsize=16)
axs[1].set_xlabel('charge [pC]', fontsize=20)
# axs[1].set_xlabel('ADCpeak [counts]', fontsize=20)
# axs[2].set_xlabel('charge[pC]/bin ', fontsize=16)
axs[2].set_xlabel('ADCpeak [counts] ', fontsize=16)
axs[0].set_title('ch1(Anode)', fontsize=16)
# axs[1].set_title('charge', fontsize=20)
axs[2].set_title('Peak', fontsize=16)
print(count,'waves',counts)
print(skipped,"skipped")

print('popt',popt)
print('pcov',pcov)
print('popt_current',popt2)
print('pcov_current',pcov2)
axs[0].tick_params(labelsize=18)
axs[1].tick_params(labelsize=18)
axs[2].tick_params(labelsize=18)

axs[2].grid()
plt.show()
