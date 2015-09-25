import subprocess
import multiprocessing
import time

class Threader:
  """
  Class for starting 2 new threads.  One that runs  a system process and one that waits and prints info to std::out or whatever you currently have set as std::out.
  Use print interval to set the wait time between prints.
  Useful for GUI subprocessing.
  """
  def __init__(self, print_interval = 3):
    self.print_interval = print_interval

  def run_system_command(self, command):
    p = subprocess.Popen(command, shell = True)
    wait_p = multiprocessing.Process(target=print_loop, args=(p, self.print_interval))
    wait_p.start()
    return p, wait_p

def print_loop(self, p, print_interval = 3):
  while not p.poll():
    if not print_interval == 0:
      time.sleep(print_interval)

    err, out = p.communicate()
    if err:
      print err
    if out:
      print out
