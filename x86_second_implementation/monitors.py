from time import sleep
import os
import abc
import logging
import sys
from random import randrange
from threading import Thread
from subprocess import call, STDOUT
from Queue import Queue

''' Monitor Interface
    monitors are meant to be subclassed by platform-specific diagnostics.
    They handle thread-level operations on streams, files and job queues
'''
class monitor(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def wait(self):
        ''' Monitor objects must to be able to block until their work is completed '''
        return

''' Line processors can be used to '''
class lineProcessor(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def process_line(self, line):
        return

    @abc.abstractmethod
    def assert_line(self, line):
        return

    @abc.abstractmethod
    def break_condition(self, line):
        return False
''' stdout_monitor objects spawn a thread that waits on the stdout of
    a target process and scans its lines through the abstract method process_line.
'''
class stdoutMonitor(monitor):
    __metaclass__ = abc.ABCMeta

    def __init__(self, process):
        self.process = process
        self.kill_thread = False
        self._spawn_thread()

    def _spawn_thread(self):
        self.t = Thread(target=self.scan_stdout)
        self.t.daemon = True
        self.t.start()

    def wait(self):
        self.kill_thread = True
        sleep(1) # TODO: optional?

    def switch_process(self, process):
        self.kill_thread = False
        self.process = process
        self._spawn_thread()

    def scan_stdout(self):
        out = self.process.stdout
        try:
            for line in iter(out.readline, b''):
                valid = True

                line = line.decode(sys.stdout.encoding) # fix line encoding TODO:optional?
                if line != None:
                    valid = self.process_line(line)

                if self.kill_thread or not valid:
                    break
        except IOError:
            ''' handles the case where the stdout is closed while blocking on it'''
            pass

    @abc.abstractmethod
    def process_line(self, line):
        return



    # filename to follow, and line processor to be used for each line
    def __init__(self, filename, line_processor):
        self.filename = filename
        self.line_processor = line_processor
        self.fail = False
        done = False
        while not done:
           try:
               self.f = open(filename, 'r')
               done = True
           except IOError:
               done = False
               logging.error("%s could not be opened", filename)
               sleep(0.1)
        self._get_line_offsets()
        try:
            self.f.seek(self.line_offset[self.line_processor.simstep])
        except IndexError:
            self.f.seek(0)
        self.time = os.path.getmtime(self.filename)
        self._injectSDC = False
        self.temp_string = ""
        self.temp_saved = False
        self._spawn_follower()

    def _get_line_offsets(self):
        self.line_offset = []
        offset = 0
        for line in self.f:
            self.line_offset.append(offset)
            offset += len(line)
        if len(self.line_offset) == 0:
            self.line_offset.append(0)
        self.f.seek(0)

    def _spawn_follower(self):
        ''' Spawn a file follower thread '''
        t = Thread(target=self.follow)
        t.daemon = True
        t.start()

    def wait(self):
        ''' Close file and let the follower thread exit '''
        self.fail = True
        done = False
        while not done:
            try:
                self.f.close()
                done = True
            except IOError:
                sleep(0.3)
        self.temp_string = ""
        self.temp_saved = False

    def follow(self):
        if (self.time == os.path.getmtime(self.filename)):
            sleep(0.4)  #wait for the first edit of the file

        while True:
            if self.fail is True:
                break

            try:
                if (self.time != os.path.getmtime(self.filename)):
                    text = self.f.read()
                    self.time = os.path.getmtime(self.filename)
                    if len(text.splitlines()) == 0:
                        self.f.close()
                        self.f = open(self.filename, 'r')
                        self._get_line_offsets()
                        try:
                            self.f.seek(self.line_offset[self.line_processor.simstep])
                        except IndexError:
                            self.f.seek(0)
                    else:
                        self.process_linelist(text.splitlines())
            except (OSError, IOError, ValueError) as e:
                sleep(0.7)

    def tobits(self, s):
        return map(int, ''.join([bin(ord(i)).lstrip('0b').rjust(8,'0') for i in s]))

    def frombits(self, l):
        return "".join(chr(int("".join(map(str,l[i:i+8])),2)) for i in range(0,len(l),8))

    def injectSDC(self):
        self._injectSDC = True

    def process_linelist(self, lines):
        if len(lines) == 0 and self._injectSDC:
            print "first"
            self.line_processor.diagnostic.fail()

        for counter, line in enumerate(lines):
            if self.fail:
                return

            if len(line.split()) == 0:
                if self._injectSDC:
                    print "second"
                    self.line_processor.diagnostic.fail()
                    return
                else:
                    continue
            
            if self._injectSDC and not self.line_processor.break_condition(line):
                self._injectSDC = False
                oldline = line[:]
                linelist = self.tobits(line)
                bit_index = randrange(len(linelist)-2) #flip one bit at random
                linelist[0:4] = map(lambda x: x^1, linelist[0:4])
                line = self.frombits(linelist)
                print "corrupted line:"
                print line

            if self.temp_saved:
                if counter == 0 and \
                len(self.temp_string.split()) < self.line_processor.expected_length():
                    line = self.temp_string + line
                    self.temp_saved = False

            if self.fail:
                return

            if not self.line_processor.break_condition(line):
                try:
                    self.line_processor.assert_line(line)
                    self.line_processor.process_line(line)
                except AssertionError:
                    ''' the last line of the lines list should be merged with the
                        first line of the next read if the type does not match
                    '''
                    if len(line.strip().split()) >= self.line_processor.expected_length() and list(line.strip())[-1] != '-':
                        print line
                        self.line_processor.diagnostic.fail()
                        self.fail = True
                    elif counter == len(lines) - 1:
                        self.temp_string = line
                        self.temp_saved = True

