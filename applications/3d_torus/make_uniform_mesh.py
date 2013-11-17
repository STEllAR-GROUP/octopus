#!/bin/bash

from numpy import zeros 
from optparse import OptionParser

args = OptionParser(usage="%prog [file] [floor]").parse_args()[1]

file = args[0]
floor = float(args[1])

f = open(file, 'r')

x = {}
y = {}
rho = {}
tracker = {}

master = {}

try:
    while True:
        line = f.next().split()

        x_here = float(line[0])
        y_here = float(line[1])
        rho_here = float(line[2])
        tracker_here = float(line[3])

        x[x_here] = 0
        y[y_here] = 0
        rho[rho_here] = 0
        tracker[tracker_here] = 0

        master[(x_here, y_here)] = (rho_here, tracker_here)

except StopIteration:
    pass

for i in sorted(x.iterkeys()):
    for j in sorted(y.iterkeys()):
        if (i, j) in master:
            print i, j, 0.0, master[(i, j)][0], master[(i, j)][1]
        else:
            print i, j, 0.0, floor, floor

    print 

