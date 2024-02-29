import numpy as np
import matplotlib.pyplot as plt
import sys 
from scipy import integrate


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
counts = 0

fig, axs = plt.subplots(1,2, figsize=(11,3.7), constrained_layout=True)

for line in f:
	linestrip=line.rstrip()
	tokens=linestrip.split()
	if readstate==WAIT_FOR_HIT:
		if len(tokens)!=1:continue
		if tokens[0]!="V1Y": continue
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
			
			if(len(ch1)<500 and len(ch2)<500):
				if max(ch1) > -1:
					counts += 1
					if counts > 0 and counts < 100:
						axs[0].plot(ch1)
						axs[1].plot(ch2)
				
		else:
			readstate=CH2_IN_PROGRESS
		continue

print(counts, 'waveform')
axs[0].set_ylabel('ADC[counts]')
# axs[0].set_xlabel('time[nsec]')
# axs[0].set_ylim(324,340)
# axs[0].set_xlabel('bin')
axs[0].set_xlabel('time[nsec]')
axs[1].set_xlabel('time[nsec]')
axs[0].set_title('ch1(Anode)')
axs[1].set_title('ch2(Dy8)')
axs[0].set_xticks([ 0,6,12,18,24,30 ])
axs[0].set_xticklabels([ 0,100,200,300,400,500 ])
axs[1].set_xticks([ 0,6,12,18,24,30 ])
axs[1].set_xticklabels([ 0,100,200,300,400,500 ])
plt.show()
