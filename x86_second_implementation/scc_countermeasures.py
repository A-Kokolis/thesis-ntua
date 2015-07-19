import abc
import logging
from time import sleep, time
from subprocess import call, check_output
from monitors import corePinger
from config import sim_dump_location, devel
import infoli_diagnostics
import sys
import random


class countermeasure(object):
    ''' Countermeasure class '''
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def perform(self):
        return


''' defines an ordering among the different countermeasures, based on MTTR '''
countermeasure_enum = {
    'restartSimulation':0,
}



class restartSimulation(countermeasure):
    """ Restarts the simulation """
    __name__ = 'restartSimulation'

    def __init__(self, manager):
        self.manager = manager

    def perform(self):
        logging.info("performing the Restart Simulation countermeasure")
        print self.manager.checkpoints
        if any(isinstance(x, infoli_diagnostics.infoliOutputDivergence) for x in self.manager.failed_diagnostics()):   #infoli-specific
            # check if the SDC detection diagnostic has failed, and use the SDC checkpoint
            print sorted(self.manager.checkpoints)
            checkpoint = max(self.manager.checkpoints)
        else:
            checkpoint = max(self.manager.checkpoints)
        #randsleep=random.randint(2,4)
        #print('randsleep=',randsleep)
        #sleep(randsleep)
        print("The mttr_values are:",self.manager.mttr_values)
        print("Calling dvfs: ")
        self.manager.dvfs.dvfsOperation(checkpoint)

        print "Restarting from step" + str(checkpoint)
        logging.info("Restarting from step " + str(checkpoint))

        with self.manager.lock:
            self.manager.simrun([self.manager.restart_exec] + self.manager.exec_list[1:]) # use False as extra last argument to avoid piping stdout for diagnostics - useful for measurements
        logging.info("Restart Simulation countermeasure completed")
        return True

