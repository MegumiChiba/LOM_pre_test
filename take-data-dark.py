
import sys, time
import datetime
import serial
import threading

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

if len(sys.argv)!=7:
	print("Usage:")
	print("python take-spe-data.py port outfile waittime nfile PMTID sleeptime")
	sys.exit(1)

port=sys.argv[1]
outfile=sys.argv[2]
runtime_text=sys.argv[3]
Nfile = int(sys.argv[4])
PMTID = sys.argv[5]
sleeptime = sys.argv[6]

import subprocess
subprocess.run([r"/home/icecube603/Chiba/STM32Cube/bin/STM32_Programmer_CLI","-c",f"port={port}","br=230400","--start","0x08000000"])

# Setting
cvoltage=80.5						#sys.argv[4]
runtime_data = 60*5
sleep = float(sleeptime)    #[hours]

baudrate=1181818   # highest exact value that works in test setup
#baudrate= 19200#*20  # highest exact value that works in test setup
#baudrate= 9600   # highest exact value that works in test setup

try:
	s=serial.Serial(port,baudrate,timeout=None)
except:
	print(f"Failed to open port \"{port}\", exiting.")
	sys.exit(1)

#try:
#	filename = f'{outfile}{datetime.datetime.fromtimestamp(time.time())}.txt'
#	f=open(filename,'wb')
#except:
#	print(f"Failed to open output file \"{filename}\", exiting.")
#	sys.exit(1)

try:
	runtime=int(runtime_text)
	if runtime<=0: raise
except:
	print(f"Invalid run time \"{runtime_text}\", exiting.")
	sys.exit(1)

print('sleep',sleep,'hours')
for i in range (int(sleep*60)):
	time.sleep(60)
	print(i,'/',int(sleep*60))

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

quickscan = f"quickscan {cvoltage}"
voltage = f"voltage {cvoltage}"

s.timeout=0  #timeout is implemented by hand, not in Serial

# change baud rate from ABR to high value (latter should be exact value on wuBase)
# N.B. have to send the command using a lower initial baud rate, in auto-baud mode
s.baudrate=9600  # lower value used when sending the command to change to higher value

for n in range(int(Nfile)):
	try:
		#filename = f'{outfile}{time.ctime()}.txt'
		#filename = f'{outfile}{time.strftime("%H%M%S", time.strptime(time.ctime()))}.txt'
		filename = f'{outfile}{n}.txt'
		f=open(filename,'wb')
	except:
		print(f"Failed to open output file \"{filename}\", exiting.")
		sys.exit(1)

	tnow = time.strftime("%H%M%S", time.strptime(time.ctime()))
	f.write(bytearray(f"time {tnow}\n",'utf-8'))
	#f.write(time.strftime("%H%M%S", time.strptime(time.ctime())))
	n=send_and_receive(f"Ubaud {baudrate}",1.0)
	s.baudrate=baudrate
	time.sleep(1.)
	print('setup')
	# send a series of setup lines
	# for each line, specify text, comms timeout and how long to sleep after command executed
	setup_lines=[("status",2.,0.1),
			("get_Uid",1.,0.5),
			("get_avg_temp",0.5,0.1),
			 ("quickscanrange 95000 125000",0.5,0.1),
			(quickscan,0.5,15),
			("get_freq",1.,0.1),
			(voltage ,3.,10.),
				#("pwmfreq 110169",0.5,0.1),
			("reportavg",0.5,0.2),
			("monstatus",0.5,0.5),	     
			("dac 1 600",0.5,0.2),              
			("fpgaload",0.5,0.5),
			("adcconfig",0.5,0.5),
			("fpgaload",1.,0.5),
			("fpgatrig 0",0.5,0.5),
			#("fpgatrig 1",0.5,0.5),                          
			("flush_events",1.,1.)	  
			#("Uscan 100 100000 115000 -1",0.5,0.1),
			#("Uscanprint",0.5,0.1),
			]

	for (sendline,timeout,sleeptime) in setup_lines:
		while True:
			n=send_and_receive(sendline,timeout)
			if n > 2:
				time.sleep(sleeptime)
				break
			elif n<=2:
				print(f"{time.ctime(time.time())}  No response, aborting. Retrying.")
			time.sleep(sleeptime)


	# set up run start and duration
	runstart=time.time()
	print(f"{time.ctime(time.time())}  Start pre run, {runtime} seconds", flush=True)

	# start the reception in a thread
	request_stop_reception=False
	request_abort=False
	rx_thread=threading.Thread(target=send_and_receive,args=("send_batch -1 1",10.0))
	rx_thread.start()

	# wait for reception to end either because an interval greater than the timeout
	# didn't deliver any new data, or because the runtime is up
	while True:
		if send_and_receive_running==False:
			#print(f"{time.ctime(time.time())}  Receive complete", flush=True)
			#break
			dsize = 0
		if time.time()>runstart+runtime and request_stop_reception==False:
			request_stop_reception=True
			print(f"{time.ctime(time.time())}  Request stop reception", flush=True)
			break
		if time.time()>runstart+runtime+stop_timeout:
			request_abort=True
			print(f"{time.ctime(time.time())}  Request abort reception", flush=True)
		print(f"{time.ctime(time.time())}  (Progress: {send_and_receive_nbytes} bytes)", flush=True)
		time.sleep(1)

	# make sure the reception thread is really gone
	rx_thread.join(5)
	if rx_thread.is_alive():
		print(f"{time.ctime(time.time())}  Error: rx thread failed to complete", flush=True)

	### datataking mode
	time.sleep(1.)
	n=send_and_receive("monstatus",3.0)
	time.sleep(1.)
	print('setup')
	# send a series of setup lines
	# for each line, specify text, comms timeout and how long to sleep after command executed
	setup_lines=[("status",2.,0.1),
			("get_Uid",1.,0.5),
			("get_avg_temp",0.5,0.1),
			(quickscan,0.5,15),
			("get_freq",1.,0.1),
			(voltage ,3.,10.),
				#("pwmfreq 110169",0.5,0.1),
			("reportavg",0.5,0.2),
			("monstatus",0.5,0.5),	     
			("dac 1 600",0.5,0.2),              
			("fpgaload",0.5,0.5),
			("adcconfig",0.5,0.5),
			("fpgaload",1.,0.5),
			("fpgatrig 0",0.5,0.5),
			("fpgatrig 1",0.5,0.5),                          
			("flush_events",1.,1.)	  
			#("Uscan 100 100000 115000 -1",0.5,0.1),
			#("Uscanprint",0.5,0.1),
			]

	for (sendline,timeout,sleeptime) in setup_lines:
		while True:
			n=send_and_receive(sendline,timeout)
			if n > 2:
				time.sleep(sleeptime)
				break
			elif n<=2:
				print(f"{time.ctime(time.time())}  No response, aborting. Retrying.")
			time.sleep(sleeptime)
	# set up run start and duration
	
	runstart=time.time()
	print(f"{time.ctime(time.time())}  Start run, {runtime_data} seconds", flush=True)

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
			break
		if time.time()>runstart+runtime_data and request_stop_reception==False:
			request_stop_reception=True
			print(f"{time.ctime(time.time())}  Request stop reception", flush=True)
		if time.time()>runstart+runtime_data+stop_timeout:
			request_abort=True
			print(f"{time.ctime(time.time())}  Request abort reception", flush=True)
		print(f"{time.ctime(time.time())}  (Progress: {send_and_receive_nbytes} bytes)", flush=True)
		time.sleep(1)

	# make sure the reception thread is really gone
	rx_thread.join(5)
	if rx_thread.is_alive():
		print(f"{time.ctime(time.time())}  Error: rx thread failed to complete", flush=True)

	n=send_and_receive("monstatus",3.0)
	time.sleep(1.)
	if n < Nfile-1:
		f.close()
	else:
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

subprocess.run(["python","ascii2hdf_darkmonitor.py",outfile,sys.argv[4],PMTID])
subprocess.run(["python","read_pdf_darkmonitor.py",outfile,sys.argv[4]])
