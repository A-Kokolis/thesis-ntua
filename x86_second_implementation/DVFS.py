import shlex
import os.path,time
from subprocess import PIPE,Popen
import logging
from dvfs_config import *

class dvfs(object):
    "executing dvfs changes"
    __name__="dvfs"
    def __init__(self,manager):
        self.manager=manager
        self.slack=0
        self.currentMultiplier=1
        self.determinedMultiplier=1
        self.slackReclaimingTime=0
        self.prevSlack=False
        self.freqProcThread=None

    def dvfsOperation(self,checkpoint):
        print("Inside dvfsOperation")
        #check whether we have reclaimed any execution slack so far
        if self.prevSlack:
            while self.freqProcThread.poll() is None:
                self.freqProcThread.kill()
            returnCode=self.freqProcThread.returncode
            #print('The returnCode is = ',returnCode)
            if (returnCode==0):
                print('Frequency proc has ended')
                self.slack=0
                self.currentMultiplier=defaultMultiplier
            else:
                print('I killed freq proc')

        totalckptOverhead=0

        #update slack adding the time overhead for restarting from a previous simStep
        #we do not calculate this overhead in updateSlack cause we also need the first rollBack
        rollbackTime=self.timeFromLastCkp(checkpoint)
        if rollbackTime>0:
            self.slack= self.slack - rollbackTime
        #print('rollbackTime=',rollbackTime)
        #print('slack with rollbackTime=',self.slack)

        #we also include the overhead due to checkpointing
        with open(ckptOverheadTimeFile) as f:
            for line in f:
                data=line.split()
                ckptOverhead=float(data[0])
	
	    mttf=sum(self.manager.mttf_values)/float(len(self.manager.mttf_values))
	    #print('DVFS: mttf=',mttf)	
        if not self.manager.mttr_values:
            print('DVFS: No TTR values yet')
            #print('ckptOverhead=',ckptOverhead)
            #print('ckptInterval=',self.manager.exec_list[-1])
            totalckptOverhead=(checkpoint/int(self.manager.exec_list[-1]))*ckptOverhead
            #print('totalckptOverhead=',totalckptOverhead)
            self.slack=self.slack-totalckptOverhead
            print('new slack=',self.slack)
            logging.info('DVFS: The slack is %f',self.slack)
        else:
            if (len(self.manager.checkpoints)>2):
                totalckptOverhead=((checkpoint-self.manager.checkpoints[-2])/int(self.manager.exec_list[-1]))*ckptOverhead
            #print('ckptOverhead=',ckptOverhead)
            #print('totalckptOverhead=',totalckptOverhead)
            self.slack=self.slack-totalckptOverhead
            #print(self.slack)
            self.updateSlack()
            print('slack=',self.slack)
            logging.info('DVFS: The slack is %f',self.slack)

        self.determineMultiplier(self.slack,mttf)
        print('determinedMultiplier=',self.determinedMultiplier)

        if (self.slack<0):
            if (self.determinedMultiplier>multiplierList[6]):
                multiplier=multiplierList[6]
                frequency=freqList[6]
            elif (self.determinedMultiplier>multiplierList[5]):
                multiplier=multiplierList[5]
                frequency=freqList[5]
            elif (self.determinedMultiplier>multiplierList[4]):
                multiplier=multiplierList[4]
                frequency=freqList[4]
            elif (self.determinedMultiplier>multiplierList[3]):
                multiplier=multiplierList[3]
                frequency=freqList[3]
            elif (self.determinedMultiplier>multiplierList[2]):
                multiplier=multiplierList[2]
                frequency=freqList[2]
            else:
                multiplier=multiplierList[1]
                frequency=freqList[1]

            if (multiplier!= self.currentMultiplier):
                self.currentMultiplier=multiplier
                print('New multiplier=',self.currentMultiplier)
                logging.info('DVFS: Multiplier changed to %f',self.currentMultiplier)
                self.changeFreqMul(frequency)

            if (self.currentMultiplier!=defaultMultiplier):
                self.slackReclaimingTime=self.determineSlackReclaimingTime(self.slack)
                print('slackReclaimingTime=',self.slackReclaimingTime)
                self.prevSlack=True
                self.startReclaimingProc(self.slackReclaimingTime)
        else:
            if (self.currentMultiplier!=defaultMultiplier):
                self.currentMultiplier=defaultMultiplier
                print('new multiplier=',self.currentMultiplier)
                logging.info('DVFS: Multiplier changed to %f',self.currentMultiplier)
                self.changeFreqMul(defaultFreq)

    def startReclaimingProc(self,reclaimTime):
       cmd=freqProcCommand % (reclaimTime)
       self.freqProcThread=Popen(shlex.split(cmd))

    def determineSlackReclaimingTime(self,slack):
        return (-slack)/(self.currentMultiplier-defaultMultiplier)

    def determineMultiplier(self,slack,mttf):
        self.determinedMultiplier=(-slack+defaultMultiplier*mttf)/mttf

    def changeFreqMul(self,frequency):
       print('Changing frequency to ',frequency)
       cmd=freqChangeCommand % (frequency)
       freqProcess=Popen(shlex.split(cmd),stdout=PIPE,stderr=PIPE,stdin=PIPE)
       (out,err)=freqProcess.communicate()
       logging.info('DVFS: frequency changed to %d',frequency)

    def updateSlack(self):
        #include the time overhead for the application to restore its data
        with open(applicationRestartTimeFile) as f:
            for line in f:
                data=line.split()
                appRestartTime=float(data[0])
                #print("  ")
                #print('applicationRestartTime=',data[0])
                #print("  ")

        self.slack=self.slack+(self.manager.mttf_values[-1]*self.currentMultiplier-self.manager.mttf_values[-1]*defaultMultiplier)-self.manager.mttr_values[-1]-appRestartTime

        #print('applicationRestartTime=',appRestartTime)


    def changeVoltage(self,divisor):
        print('Changing voltage to divisor ',divisor)
        cmd = VoltageChangePath % divisor
        voltChange=Popen(shlex.split(cmd),stdin=PIPE,stdout=PIPE)
        (output,err)=voltChange.communicate()
        print(output)

    def timeFromLastCkp(self,restartStep):
        file=ckpTimeFile % restartStep
        creationTime=os.path.getmtime(file)
        return (self.manager.timestamp-creationTime)

    def printTest(self,checkpoint):
        result=self.timeFromLastCkp(checkpoint)
        print('Rollback is ',result)
        with open("ckptRollback.txt","a+") as myfile:
            myfile.write("Rollback is "+str(result)+'\n')
            print(result)


