import subprocess
import multiprocessing
import time
import atexit
import os
import sys
#from __future__ import print_function

#Need to do more with thread safety
class Threads:
  """
  Class for managing threads, killing them on program exit.
  """
  def __init__(self):
    self.thread_classes = []
    self.total_allowed_threads = multiprocessing.cpu_count() - 1
    #atexit.register(self.kill_all); Disable default use of this for app-specific uses.

  def set_allowed_threads(self, n):
    self.total_allowed_threads = n

  def kill_all(self):

    for th in self.thread_classes:
      if th.sub_p and th.sub_p.poll() != 0:
        th.sub_p.kill()

      if th.multi_p:
        th.multi_p.terminate()

  def append(self, thread):
    self.thread_classes.append(thread)

  def __len__(self):
    return len(self.thread_classes)

  def n_alive(self):
    i = 0
    for th in self.thread_classes:

      if th.sub_p and th.sub_p.poll() == None:
        i+=1

      if th.multi_p and th.multi_p.exitcode == None:
        i+=1

    return i

  def is_alive(self, pid):
    th = self.thread_classes[pid]
    if th.sub_p and th.sub_p.poll() == None:
      return True
    elif th.multi_p and th.multi_p.exitcode == None:
      return True
    else:
      return False

  def get_exitcode(self, pid):
    th = self.thread_classes[pid]
    if th.sub_p:
      return th.sub_p.poll()
    if th.multi_p:
      return th.multi_p.exitcode

  def new_thread_allowed(self):
    if self.n_alive() == self.total_allowed_threads:
      return False
    else:
      return True

global threads
threads = Threads()


class Threader(object):
  """
  Class for starting 2 new threads.  One that runs  a system process and one that waits and prints info to std::out or whatever you currently have set as std::out.
  Use print interval to set the wait time between prints.
  Useful for GUI subprocessing.
  """
  def __init__(self, print_interval = 0):
    self.print_interval = print_interval
    self.sub_p = None; #A subprocess.Popen object
    self.multi_p = None; #A multiprocessing.Process object.

  def run_system_command(self, command):
    """
    Run a system command using Popen.  Prints out at end.  Probably should remove this.
    :param command:
    :return:
    """
    def start():
      ##You absolutely NEED to use PIPE for standard out or you will hang.
      p = subprocess.Popen(command, shell = True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
      wait_p = multiprocessing.Process(target=print_loop, args=(p, self.print_interval))
      wait_p.start()

      self.sub_p = p; self.multi_p = wait_p
      threads.append(self)
      print("Running Thread ID: "+repr(len(threads)+1)+" "+command)
      print("Total Running Threads: "+repr(threads.n_alive()))

      return p, wait_p


    if threads.new_thread_allowed():
      return start()

    else:
      print "Too many threads running.  Waiting to start next..."
      while not threads.new_thread_allowed():
        pass
      return start()

  def run_functions(self, functions):
    """
    Run a bunch of lambda functions together with multiprocessing
    :param functions:
    :return:
    """
    def run(fs):
      print "Running Thread ID: "+repr(len(threads)+1)
      print("Total Running Threads: "+repr(threads.n_alive()))

      for f in fs:
        f()

      print "Thread ID "+repr(len(threads)+1)+" done!"

    if threads.new_thread_allowed():
      p = multiprocessing.Process(target=run, args = (functions,))
    else:
      while not threads.new_thread_allowed():
        pass
      p = multiprocessing.Process(target=run, args = (functions,))

    p.start()
    self.multi_p = p

    threads.append(self)



def print_loop(p, print_interval = 0):

  #This all really does not work with the PIPE in subprocess.  Not sure how to get it to work.
  # The idea was to keep it printing stuff as things were running.  Now, it just prints at end.
  while p and not p.poll():
    if not print_interval == 0:
      time.sleep(print_interval)

    try:
      err, out = p.communicate()
      if err:
        print err
      if out:
        print out
    except ValueError:
      #Due to subprocess.Popen being fairly shitty.
      # If you don't set the returncode, the value error basically breaks the process, so returncode is NEVER set.
      # The valueerror comes from trying to use p.communicate on a process that is already done?  Perhaps this is done between the while and the communicate?
      # I don't know, but seems pretty fishy!!
      p.returncode = 0
      break

  print "Process done!"
  print repr(threads.n_alive())+" still running..."

def test_function(i, extra = ""):
  print extra+" "+repr(i)+"\n"
  time.sleep(3+i)
  print extra+" "+repr(i)+"done\n"

if __name__ == "__main__":
  functions = []
  functions.append(lambda: test_function(1, "one: "))
  functions.append(lambda: os.system("echo This works baby!; echo $PATH"))

  f2 = []
  f2.append(lambda: test_function(3, "three: "))

  f3 = [lambda: test_function(5, "five: ")]

  threader = Threader()
  threader.run_functions(functions)
  threader.run_functions(f2)
  threader.run_functions(f3)