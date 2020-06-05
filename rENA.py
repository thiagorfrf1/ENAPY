# -*- coding: UTF-8 -*-
import time

from rpy2.robjects import r, pandas2ri, Formula
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
rENA = importr('rENA')


pandas2ri.activate()  # makes some conversions automatic

print("Keys")
pr = ro.r("print")
colMeans = ro.r("colMeans")
rescale = ro.r("scales::rescale")
rmax = ro.r("max")

def print_ena_functions():
    print(rENA._rpy2r.keys())


