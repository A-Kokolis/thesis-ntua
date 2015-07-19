import shlex
from subprocess import Popen,PIPE
import sys
from time import sleep
from dvfs_config import *
import logging


logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s',datefmt='%d/%m/%Y %I:%M:%S %p',\
                    filename='infoli.log', filemode='a+',level=logging.DEBUG)

#sleep till the slack is reclaimed and then set the frequency to its default value
sleepTime=float(sys.argv[1])

print('sleepTime=',sleepTime)

sleep(sleepTime)

cmd = freqChangeCommand % (defaultFreq)

process=Popen(shlex.split(cmd),stdin=PIPE,stderr=PIPE,stdout=PIPE)

(out,err)=process.communicate()

logging.info('DVFS: frequency changed to %d',800000)

print('  ')
print('Frequency returned to default value')
print(' ')
