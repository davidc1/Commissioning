import sys
from ROOT import *
import matplotlib.pyplot as plt
import numpy

from root_numpy import root2array

# Use the file from command line:
if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

myfile = sys.argv[-1]
mytree = "waveformsub"

# Get the subtraction waveform for these blocks (one per plane):

# branches = ["run","subrun","event","stepSize",
#             "pedestal_0","rms_0","rmsCorrected_0",
#             "pedestal_1","rms_1","rmsCorrected_1",
#             "pedestal_2","rms_2","rmsCorrected_2"]

branches = ["run","subrun","event","stepSize",
            "pedestal_0","rms_0",
            "pedestal_1","rms_1",
            "pedestal_2","rms_2"]

branch_map = dict()
i = 0
for branch in branches:
    branch_map.update({branch: i})
    i += 1

data = root2array(myfile,mytree,branches=branches)

plots = []

plots.append(plt.figure(figsize=(20,10)))
plt.plot(data[0][branch_map["rms_0"]],label="U Wires RMS")
# plt.plot(data[0][branch_map["rmsCorrected_0"]],label="U Wires RMS corrected")
plt.legend()

plots.append(plt.figure(figsize=(20,10)))
plt.plot(data[0][branch_map["rms_1"]],label="V Wires RMS")
# plt.plot(data[0][branch_map["rmsCorrected_1"]],label="V Wires RMS corrected")
plt.legend()

plots.append(plt.figure(figsize=(20,10)))
plt.plot(data[0][branch_map["rms_2"]],label="Y Wires RMS")
# plt.plot(data[0][branch_map["rmsCorrected_2"]],label="Y Wires RMS corrected")
plt.legend()

# plots.append(plt.figure(figsize=(20,10)))
# plt.plot(data[0][branch_map["rms_2"]]/data[0][branch_map["rmsCorrected_2"]],label="Y Wire RMS ratio",color="r")
# plt.plot(data[0][branch_map["rms_0"]]/data[0][branch_map["rmsCorrected_0"]],label="U Wire RMS ratio",color="b")
# plt.plot(data[0][branch_map["rms_1"]]/data[0][branch_map["rmsCorrected_1"]],label="V Wire RMS ratio",color="g")
# plt.xlabel("Wire Number")
# plt.legend()


plots.append(plt.figure(figsize=(20,10)))
plt.plot(data[0][branch_map["pedestal_0"]],label="U Wire Pedestals")
plt.legend()


plt.show()