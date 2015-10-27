
from multiprocessing import Process
import time

#Author: Jordan Willis



class JobDistributor():
    def __init__(self, jobs, limit=10, verbose=True):
        """
        Runs a set of processes using a given limit.
        When one is done, will start the next one untill all jobs from the list are complete.

        :param jobs: A list of multiprocessing.Subprocess objects.  Not yet started.
        :param limit:
        :param verbose:
        :return:
        """
        self.jobs = jobs
        self.limit = limit
        self.verbose = verbose
        self.sleep_time = 2

    def set_jobs_limit(self, limit):
        self.limit = limit

    def set_sleep_time(self, time_):
        self.sleep_time = time_

    def execute(self):
        job_counter = 0
        jobs_running = []
        while self.jobs:
            if self.verbose:
                print "jobs left -> {}".format(len(self.jobs))
            for i in self.jobs:
                if job_counter <= self.limit:
                    if self.verbose:
                        print "starting {}".format(i.name)
                    i.start()
                    job_counter += 1
                    self.jobs.remove(i)
                    jobs_running.append(i)
            for jobs_ in jobs_running:
                if jobs_.is_alive():
                    continue
                else:
                    if self.verbose:
                        print "print job {} is done, removing".format(jobs_.name)
                        jobs_running.remove(jobs_)
                        job_counter -= 1
            if self.verbose:
                print "sleeping distributor"
            time.sleep(self.sleep_time)
