import sys, time
import serial
import threading
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit  
import warnings
warnings.filterwarnings("ignore")
import sys 
from HDFwriterGen2DOMPMTinspection import HDFwriterGen2DOMPMTinspection

# collect wuBase data and save in file
# N.B. frequency and voltage settings are hard-coded in setup_lines[]
#
# usage: python take-spe-data.py port outfile runtime
#        port = serial port for UART connection to wuBase
#        outfile = name of output file
#        
# comms: start in ABR mode, 9600 baud
#        then switch to high baud rate for data taking
#        finally switch back to ABR / 9600 and confirm

if len(sys.argv)!=4:
	print("Usage:")
	print("python take-spe-data.py port outfile PMTID")
	sys.exit(1)

port = sys.argv[1]
outname = sys.argv[2]
PMTID    = sys.argv[3] #9795
runtime_text='20' # sys.argv[3]
wubaseID = -1
pdfname = outname + "_gainplots.pdf"

import subprocess
subprocess.run([r"/home/icecube603/Chiba/STM32Cube/bin/STM32_Programmer_CLI","-c",f"port={port}","br=230400","--start","0x08000000"])

baudrate=1181818   # highest exact value that works in test setup
#baudrate= 19200#*20  # highest exact value that works in test setup
#baudrate= 9600   # highest exact value that works in test setup

try:
	s=serial.Serial(port,baudrate,timeout=None)
except:
	print(f"Failed to open port \"{port}\", exiting.")
	sys.exit(1)

try:
	runtime=int(runtime_text)
	if runtime<=0: raise
except:
	print(f"Invalid run time \"{runtime_text}\", exiting.")
	sys.exit(1)

# define character sequence that denotes end of command response
# (error responses will be handled by timeout)
ok_bytes=bytearray("OK\n",'utf-8')
ok_bytes_len=len(ok_bytes)
okok_bytes=bytearray("OK\nOK\n",'utf-8')
okok_bytes_len=len(okok_bytes)

# buffer holding the last received part of response
# (data should be shuffled from the end to the beginning as needed so
# it doesn't run out of room; response_half_full should be used as a
# threshold when to do that shuffling depending on how much new data
# is available)
response_length=okok_bytes_len
response_half_full=response_length 
response_buffer=bytearray(response_half_full+response_length)

# indicator if send_and_receive is running
send_and_receive_running=False
# set flag to request stopping reception if send_and_receive is in a Thread
request_stop_reception=False
stop_timeout=10
# abort flag for send_and_receive
request_abort=False
# how many bytes have been received
send_and_receive_nbytes=0
# indicator that response to command started with ?
send_and_receive_command_error=False

# procedure for issuing command and collecting the response,
# writing each data chunk to output file while also watching
# for response terminator or timeout
# N.B. response terminator is both a character sequence and the
# stopping of further characters, i.e. an input data chunk will
# never include more characters after the terminator
def send_and_receive(sendline,timeout):
	global request_stop_reception
	global request_abort
	global send_and_receive_running
	global send_and_receive_nbytes
	global send_and_receive_command_error
	#
	request_stop_reception=False # clear error / abort flags at start
	request_abort=False
	send_and_receive_command_error=False
	#
	send_and_receive_running=True
	stop_reception_requested=False  # expect OK\nOK\n if stop requested
	response_n=0   # how many characters in response tail buffer are occupied
	send_and_receive_nbytes=0  # how many characters have been received in total
	reftime=time.time()
	reftime_nbytes=0  # how many bytes when reftime established
	# send the command

	if len(sendline)>0:
		# first read & save any characters remaining in serial input buffer, 
		# to avoid counting them as a response to the new command
		data=s.read(s.in_waiting or 1)
		ndata=len(data)

		if ndata>0: 
			f.write(data)
			send_and_receive_nbytes+=ndata
		# send the command
		linecr=sendline+'\n'
		print(str(time.ctime(time.time()))+"  "+sendline, flush=True)
		f.write(linecr.encode('utf-8'))
		s.write(linecr.encode('utf-8'))
	# collect the response and write to file, also watch for terminator
	while True:
		if request_stop_reception and not stop_reception_requested:
			s.write(b"OK\n")
			stop_reception_requested=True
		if request_abort:
			break
		data=s.read(s.in_waiting or 1)
		ndata=len(data)
		if ndata>0:
			if send_and_receive_nbytes==0: #check for ? at beginning of response
				if data[0]==b'?': command_response_error=True
			send_and_receive_nbytes+=ndata
			f.write(data)
			if ndata>response_length:
				response_buffer[:response_length]=data[-response_length:]
				response_n=response_length
			else:
				if response_n>response_half_full:
					nold=response_length-ndata
					response_buffer[:nold]=response_buffer[response_n-nold:response_n]
					response_buffer[nold:nold+ndata]=data
					response_n=response_length
				else:
					response_buffer[response_n:response_n+ndata]=data
					response_n+=ndata 
		else:
			if response_buffer[response_n-ok_bytes_len:response_n]==ok_bytes:
				if stop_reception_requested==False: break
				if response_buffer[response_n-okok_bytes_len:response_n]==okok_bytes: break
			if send_and_receive_nbytes>reftime_nbytes:
				reftime=time.time()
				reftime_nbytes=send_and_receive_nbytes
			else:
				if time.time()>reftime+timeout: break
	# log amount of collected data
	print(f"{time.ctime(time.time())}  Received {send_and_receive_nbytes} bytes", flush=True)
	send_and_receive_running=False
	return send_and_receive_nbytes


vlist = [80,82,84,86,88,90,92,94,96,98]
vlist = [76,78,80,82,84,86,88,90,92,94]
#vlist = [74,76,78,80,82,84,86,88,90,92]
#vlist = [85]

for i in range(len(vlist)):
	RETRY = 0
	outfile = outname + str(vlist[i]) + '.txt'
	try:
		f=open(outfile,'wb')
	except:
		print(f"Failed to open output file \"{outfile}\", exiting.")
		sys.exit(1)
	scanset = "quickscan" + ' ' + str(vlist[i])
	vset = "voltage" + ' ' + str(vlist[i])
	s.timeout=0  #timeout is implemented by hand, not in Serial

	# change baud rate from ABR to high value (latter should be exact value on wuBase)
	# N.B. have to send the command using a lower initial baud rate, in auto-baud mode
	s.baudrate=9600  # lower value used when sending the command to change to higher value
	time.sleep(1.)
	n=send_and_receive(f"Ubaud {baudrate}",1.0)
	s.baudrate=baudrate
	time.sleep(1.)
	print('setup')
	# send a series of setup lines
	# for each line, specify text, comms timeout and how long to sleep after command executed
	setup_lines=[("status",2.,0.1),
		("get_Uid",1.,0.5),
		("get_avg_temp",0.5,0.1),
		 (scanset,0.5,15),
		  ("get_freq",1.,0.1),
		(vset,3.,10.),
	     ("reportavg",0.5,0.2),
	     ("monstatus",0.5,0.5),	     
	     ("dac 1 650",0.5,0.2),              
	     ("fpgaload",0.5,0.5),
	     ("adcconfig",0.5,0.5),
	     ("fpgaload",1.,0.5),
	     ("fpgatrig 0",0.5,0.5),
	     ("fpgatrig 1",0.5,0.5),                          
	     ("flush_events",1.,1.)	  
		 ]
	
	for (sendline,timeout,sleeptime) in setup_lines:
		n=send_and_receive(sendline,timeout)
		if n==0:
			print(f"{time.ctime(time.time())}  No response, aborting. Retrying.")
			n=send_and_receive(sendline,timeout)
			if n==0:
				print(f"{time.ctime(time.time())}  No response, aborting.", flush=True)
				f.close()
				s.close()
				#sys.exit(1)
		time.sleep(sleeptime)


	# set up run start and duration
	runstart=time.time()
	print(f"{time.ctime(time.time())}  Start run, {runtime} seconds", flush=True)

	# start the reception in a thread
	request_stop_reception=False
	request_abort=False
	rx_thread=threading.Thread(target=send_and_receive,args=("send_batch -1 1",10.0))
	rx_thread.start()

	# wait for reception to end either because an interval greater than the timeout
	# didn't deliver any new data, or because the runtime is up
	while True:
		if send_and_receive_running==False:
			print(f"{time.ctime(time.time())}  Receive complete", flush=True)
			if send_and_receive_nbytes <= 10000:
				RETRY = 1
			break
		if time.time()>runstart+runtime and request_stop_reception==False:
			request_stop_reception=True
			print(f"{time.ctime(time.time())}  Request stop reception", flush=True)
		if time.time()>runstart+runtime+stop_timeout:
			request_abort=True
			print(f"{time.ctime(time.time())}  Request abort reception", flush=True)
		print(f"{time.ctime(time.time())}  (Progress: {send_and_receive_nbytes} bytes)", flush=True)
		time.sleep(1)

	# make sure the reception thread is really gone
	rx_thread.join(5)
	if rx_thread.is_alive():
		print(f"{time.ctime(time.time())}  Error: rx thread failed to complete", flush=True)

	if RETRY == 1:
		print(f'Restart the measurement of {vlist[i]} V')
		for (sendline,timeout,sleeptime) in setup_lines:
			n=send_and_receive(sendline,timeout)
			if n==0:
				n=send_and_receive(sendline,timeout)
				if n==0:
					print(f"{time.ctime(time.time())}  No response, aborting.", flush=True)
					f.close()
					s.close()
					sys.exit(1)
			time.sleep(sleeptime)


		# set up run start and duration
		runstart=time.time()
		print(f"{time.ctime(time.time())}  Start run, {runtime} seconds", flush=True)

		# start the reception in a thread
		request_stop_reception=False
		request_abort=False
		rx_thread=threading.Thread(target=send_and_receive,args=("send_batch -1 1",10.0))
		rx_thread.start()

		# wait for reception to end either because an interval greater than the timeout
		# didn't deliver any new data, or because the runtime is up
		while True:
			if send_and_receive_running==False:
				print(f"{time.ctime(time.time())}  Receive complete", flush=True)
				if send_and_receive_nbytes <= 10000:
					RETRY = 1
				break
			if time.time()>runstart+runtime and request_stop_reception==False:
				request_stop_reception=True
				print(f"{time.ctime(time.time())}  Request stop reception", flush=True)
			if time.time()>runstart+runtime+stop_timeout:
				request_abort=True
				print(f"{time.ctime(time.time())}  Request abort reception", flush=True)
			print(f"{time.ctime(time.time())}  (Progress: {send_and_receive_nbytes} bytes)", flush=True)
			time.sleep(1)

		# make sure the reception thread is really gone
		rx_thread.join(5)
		if rx_thread.is_alive():
			print(f"{time.ctime(time.time())}  Error: rx thread failed to complete", flush=True)

	n = send_and_receive("monstatus",0.5)
	# set voltage to zero and change the baud rate back to auto-baudrate mode (ABR)
	n=send_and_receive("voltage 0",3.0)
	#time.sleep(5.)
	n=send_and_receive("baud -1",1.0)
	s.baudrate=9600
	time.sleep(1.)

	# using new baudrate & ABR, verify comms still working
	n=send_and_receive("Ustatus",2.)


	f.close()
s.close()

############ save as hd5 ############
if(__name__ == "__main__"):

    MCUID    = 1234
    pwmfreq  = 12345
    temperature = 25
    PMTIDstr = "BB{0}".format(PMTID)
    userdict = {"Lasersetting":-1, "position":np.array([0,0,0]), "B-field":np.array([0,0,0])}


    # 2.0V/4096bits/1670ohm/60MHz //1.85is the correction value
    conversion_ch1 = 2.0*1.85/4096/1670.0/60e6*1e12 # pC/1 ADC sample calibrated by micro base
    conversion_ch2 = conversion_ch1/0.0115 # measured by Chiba
    yourname = "Chiba"
    description = "Data collected with no injection."


    #file    = sys.argv[1]
    #outfile = sys.argv[2]

    for i in range(len(vlist)):
        filename    = outname + str(vlist[i]) + '.txt'
        outfilename = outname + str(vlist[i]) + '.hd5'


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
        try:
            voltage10 = vol[0]
        except:
            voltage10 = 0
		
        # a class to save HDF file
        hdf = HDFwriterGen2DOMPMTinspection(Nsample_max)

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
#try:
subprocess.run(["python", "read_gain_pdf.py",outname,pdfname,PMTID])
#except:
	#print("failed to make a pdf file")
