# send email message to me for process notices


import smtplib
import mod_input
from email.mime.text import MIMEText

def SendMsg(Subject, MsgText):
	try:
		# form a message structure
		Msg = MIMEText(MsgText)
		
		# this is for me only, so hard-code the to and from addresses
		Msg['From'] = 'SITRUS'#'sol_claire@.asu.edu'
		Msg['To'] = mod_input.user_email
		
		# subject line
		Msg['Subject'] = Subject
		
		# send using the mail server on this machine
		Svr = smtplib.SMTP('localhost')
		Svr.sendmail(Msg['From'], Msg['To'], Msg.as_string())
		Svr.quit()
			
	except Exception as ex:
		
		# let me know if anything goes wrong
		print("EMail error", str(ex))

def runtime(time):
    # time = time #time.time() - start_time
    runtime = time / (60 * 60 * 24) # convert to days
    
    runtime_dd = int(runtime) # find the number of days
    time = time - (runtime_dd * (24 * 60 * 60)) # subtract days off of the total runtime in seconds
    
    runtime_hh = int(time / (60 * 60)) # find the total number of hours
    time = time - (runtime_hh * (60 * 60)) # subtract hours off the total runtime in seconds (left with just minutes and seconds)
    
    runtime_mm = int(time / 60) # find the total number of minutes
    time = time - (runtime_mm * 60) # subtract minutes off the total runtime in seconds (left with just seconds)
    
    runtime_ss = round(time,)
    
    if len(str(runtime_hh)) < 2:
        runtime_hh = f'0{runtime_hh}'
        
    if len(str(runtime_mm)) < 2:
        runtime_mm = f'0{runtime_mm}'
        
    if len(str(runtime_ss)) < 2:
        runtime_ss = f'0{runtime_ss}'

    return f'{runtime_dd}-{runtime_hh}:{runtime_mm}:{runtime_ss}'
