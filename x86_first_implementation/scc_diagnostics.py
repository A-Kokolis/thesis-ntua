import logging
import abc
import sys
from threading import Lock
from monitors import monitor, stdoutMonitor,  fileReader, lineProcessor
from injectors import processExitInjector
from core_allocator import allocate_tasks


""" Diagnostics Interface """
class diagnostic(object):
    __metaclass__ = abc.ABCMeta
    #diagnostics must have an injector list if fault injection will be used

    def __init__(self):
        self.failed = False
        self.lock = Lock()

    def fail(self):
        with self.lock:
            if not self.failed:
                logging.error("%s diagnostic failed", self.__class__.__name__)
                self.failed = True
            if not self.manager.stopped:
                self.manager.stop()

    def reinit(self):
        self.failed = False
        self.reinitialize()

    def degrade(self):
        pass

    def completed(self):
        return True

    @abc.abstractmethod
    def reinitialize(self):
        return

    @abc.abstractmethod
    def countermeasure_procedure(self):
        return


class processExit(stdoutMonitor, diagnostic):
    def __init__(self, manager):
        self.manager = manager
        stdoutMonitor.__init__(self, self.manager.simulation)
        diagnostic.__init__(self)
        self.injectors = [processExitInjector(self)]

    def process_line(self, line):
            """ Checks for SCC FAILURE messages in the app's stdout.
                Returns true if the stdout line does not cause the diagnostic to fail
            """
            if line.find("FAILURE") != -1:
                core = line[23:29]  # core number
                if line[-12:-1] == "Interrupted": # ignore "Interrupted" messages
                    return True

                try:
                    error = int(line[-4:-1])
                except ValueError:
                    """ Ignore non-SCC messsages"""
                    print line  # DEBUG
                    return True

                if error != 255:
                    # ignore error code 255, as it occurs during manual killing of the process
                    logging.error('Core %s: Process failed with error value %d', core, error)
                    self.fail()
                    return False
            return True

    def reinitialize(self):
        self.switch_process(self.manager.simulation)

    def countermeasure_procedure(self):
        from scc_countermeasures import restartSimulation
        return [[restartSimulation(self.manager)]]

