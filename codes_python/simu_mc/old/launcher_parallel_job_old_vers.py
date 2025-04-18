import sys, os
from opengate_run import *

os.chdir("/sps/gdrmi2b/milewski/simulations/monte_carlo")


runJobs_opengate("ESRF_line_perfect_slits.py", jobs=1000)
