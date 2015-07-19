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
        self.boostMultiplier=1.3
        self.currentMultiplier=1
        self.slowDivisor=3
        self.boostDivisor=2
        self.r=0

    def printTest(self,checkpoint):
        result=self.timeFromLastCkp(checkpoint)
        print('Rollback is ',result)
        with open("ckptRollback.txt","a+") as myfile:
            myfile.write("Rollback is "+str(result)+'\n')
            print(result)


    def dvfsOperation(self,checkpoint):
        print("Inside dvfsOperation")
        rollbackTime=self.timeFromLastCkp(checkpoint)
        self.r=self.r+1
        #update slack adding the time overhead for restarting from a previous simStep
        #we do not calculate this overhead in updateSlack cause we also need the first rollBack
        if rollbackTime>0:
            self.slack= self.slack - rollbackTime
        #print('rollbackTime=',rollbackTime)
        print('slack with rollbackTime=',self.slack)

        #we also include the overhead due to checkpointing
        with open(ckptOverheadTimeFile) as f:
            for line in f:
                data=line.split()
                ckptOverhead=float(data[0])
        if not self.manager.mttr_values:
            print('DVFS: There is no time overhead to reclaim yet')
           # print('ckptOverhead=',ckptOverhead)
           # print('ckptInterval=',self.manager.exec_list[-1])
            totalckptOverhead=(checkpoint/int(self.manager.exec_list[-1]))*ckptOverhead
           # print('totalckptOverhead=',totalckptOverhead)
            self.slack=self.slack-totalckptOverhead
           # print('new slack=',self.slack)
            logging.info('DVFS: There is no time overhead to reclaim yet')
        else:
            if (len(self.manager.checkpoints)>1):
                totalckptOverhead=((checkpoint-self.manager.checkpoints[-2])/int(self.manager.exec_list[-1]))*ckptOverhead
            #print('ckptOverhead=',ckptOverhead)
            #print('totalckptOverhead=',totalckptOverhead)
            else:
                totalckptOverhead=0
            self.slack=self.slack-totalckptOverhead
            #print(self.slack)
            self.updateSlack()
            print('slack=',self.slack)
            print('self.r=',self.r)
            logging.info('DVFS: The slack is %f',self.slack)
            if self.r%2==0 :
                if (self.slack>=0):
                    print('We have positive slack...')
                    if self.currentMultiplier!=1:
                        print('We reduce frequency...')
                        self.changeVoltage(self.slowDivisor)
                        self.currentMultiplier=1
                        logging.info('DVFS :Voltage Changed. Current Multiplier is %f',self.currentMultiplier)
                else:
                    print('We have negative slack...')
                    if self.currentMultiplier!= self.boostMultiplier:
                        print('We increase frequency')
                        self.changeVoltage(self.boostDivisor)
                        self.currentMultiplier=self.boostMultiplier
                        logging.info('DVFS: Voltage Changed. Current Multiplier is %f',self.currentMultiplier)



    def updateSlack(self):
        #include the time overhead for the application to restore its data
        with open(applicationRestartTimeFile) as f:
            for line in f:
                data=line.split()
                appRestartTime=float(data[0])
                #print("  ")
                #print('applicationRestartTime=',data[0])
                #print("  ")

        if (self.currentMultiplier==self.boostMultiplier):
            self.slack=self.slack+(self.manager.mttf_values[-1]*self.boostMultiplier-self.manager.mttf_values[-1]*1)-self.manager.mttr_values[-1]-appRestartTime
        else:
            self.slack=self.slack-self.manager.mttr_values[-1]-appRestartTime

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
