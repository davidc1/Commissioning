import sys
from ROOT import *
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy
import math

from root_numpy import root2array

def getBlocks(U = None,V=None,Y=None):
  if U != None:
    V = U + 7
    Y = U + 7
    return [U,V,Y]
  if V != None:
    U = V - 7
    Y = V
    return [U,V,Y]
  if Y != None:
    V = Y 
    U = Y - 7
    return [U,V,Y]  

def correlation(x,y):
  if len(x) != len(y):
    print "Error, can not determine correlation of vectors of different length"
    return None
  else:
    mean_x = numpy.mean(x)
    mean_y = numpy.mean(y)
    var_x  = numpy.std(x)
    var_y  = numpy.std(y)
  corr = 0
  if var_x == 0 or var_y == 0:
    return corr
  for i in xrange(len(x)):
    corr += (x[i] - mean_x)*(y[i] - mean_y)

  return corr / (var_x*var_y*len(x))

# Use the file from command line:
if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

myfile = sys.argv[-1]
mytree = "waveformsub"



nWires = [2400,2400,3456]

# First, get the step size of each plane:
stepsize = root2array(myfile,mytree,"stepSize")[0]
nsteps = [a/b for (a,b) in zip(nWires, stepsize)]
print nsteps


# Set up the branches needed for the correlation measurement:
branches = []
branches_map = dict()
j = 0
for i in xrange(len(stepsize)):
  nSteps = nWires[i]/stepsize[i]
  for step in xrange(nSteps):
    name = "subwaveform_"+str(i)+"_"+str(step)
    branches.append(name)
    branches_map.update({name:j})
    j += 1

# Also add on the rest of the data
tempList = ["run","subrun","event","stepSize",
            "pedestal_0","rms_0","rmsCorrected_0",
            "pedestal_1","rms_1","rmsCorrected_1",
            "pedestal_2","rms_2","rmsCorrected_2",
            "correlationMatrix"]



for item in tempList:
  branches.append(item)
  branches_map.update({item:j})
  j += 1




data = root2array(myfile,mytree,branches=branches)[0]


# Create a title for this plot:
title = "Correlated Noise, Run " + str(data[branches[branches_map['run']]])
title += ", Event " +  str(data[branches[branches_map['event']]])


# I guess copy the correlation matrix into a 2d array:
correlationMatrix = data[branches[branches_map["correlationMatrix"]]]
corrMatrix= numpy.ndarray((len(correlationMatrix),len(correlationMatrix)))

i = 0
for vec in correlationMatrix:
  j = 0
  for item in vec:
    if not math.isnan(item):
      corrMatrix[i][j] = item;
    else:
      corrMatrix[i][j] = 0;
    j += 1
  i += 1



# print numpy.amin(corrMatrix.all())


plt.rcParams.update({'font.size': 20})
fig = plt.figure(figsize=(20,10))
plt.suptitle(title,fontsize=30)

# ax2 = plt.subplot2grid((2,3), (1,0))
ax3 = plt.subplot2grid((2,3), (0,1), rowspan=2,colspan=2)

matpic = ax3.matshow(corrMatrix,origin='lower',cmap='RdBu')
plt.title("Subtracted Waveform Correlation Matrix")
ax3.get_xaxis().set_ticks([])
ax3.get_yaxis().set_ticks([])
plt.ylabel("Motherboard Pairs",fontsize=25)
# pos = fig.add_axes([0.95,0.15,0.05,0.7])
cbar = plt.colorbar(matpic)
cbar.set_clim(vmin=-1,vmax=1)

# Generate lines to divide U,V,Y:
divs = [a/b for a,b in zip(nWires,stepsize)]

plt.axhline(divs[0]-0.5)
plt.axvline(divs[0]-0.5)
plt.axhline(divs[0]+divs[1]-0.5)
plt.axvline(divs[0]+divs[1]-0.5)

# Compare these waveforms:
blocks = [2*10,2*17,2*10]

# Add arrows to show which wires we're using
arrow_length=8
ax3.arrow(blocks[0], blocks[0]-arrow_length-1.5, 0, arrow_length, head_width=0.5, head_length=1.0,color='b')
ax3.arrow(blocks[1]+nsteps[0], blocks[0]-arrow_length-1.5, 0, arrow_length, head_width=0.5, head_length=1.0,color='g')
ax3.arrow(blocks[2]+nsteps[0] + nsteps[1], blocks[0]-arrow_length-1.5, 0, arrow_length, head_width=0.5, head_length=1.0,color='r')
ax3.axhline(blocks[0],linestyle='--',color='black')
# matpic = ax.matshow(mat2_l3)
# plt.show()

plt.text(16,-8,"U Plane")
plt.text(16+nsteps[0],-8,"V Plane")
plt.text(8+nsteps[0] + nsteps[1],-8,"Y Plane")

# plt.tight_layout()


# Plot the wave forms:
ax1 = plt.subplot2grid((2,3), (0,0))

# Get the subtraction waveform for these blocks (one per plane):


wireStart = [a*b + 1 for a,b in zip(blocks,stepsize)]
wireEnd   = [a*b + b + 1 for a,b in zip(blocks,stepsize)]

waveform0_name = "subwaveform_0_"+str(blocks[0])
waveform1_name = "subwaveform_1_"+str(blocks[1])
waveform2_name = "subwaveform_2_"+str(blocks[2])




waveform0 = data[branches[branches_map[waveform0_name]]]
waveform1 = data[branches[branches_map[waveform1_name]]]
waveform2 = data[branches[branches_map[waveform2_name]]]

plt.plot(waveform0[0:500],label="U Wires " + str(wireStart[0]) + " to " + str(wireEnd[0]))
plt.plot(waveform1[0:500],label="V Wires " + str(wireStart[1]) + " to " + str(wireEnd[1]))
plt.plot(waveform2[0:500],label="Y Wires " + str(wireStart[2]) + " to " + str(wireEnd[2]))
plt.title("Sample of Correlated Waveforms")
plt.ylabel("[ADC]")
plt.legend(prop={'size':9})
plt.xlabel("Time Ticks [500 ns]")
plt.grid(True)



# Include the ffts:
ax2 = plt.subplot2grid((2,3), (1,0))


fft0 = numpy.fft.rfft(waveform0)
freqs0 = numpy.fft.rfftfreq(len(waveform0),0.5E-3)
fft1 = numpy.fft.rfft(waveform1)
freqs1 = numpy.fft.rfftfreq(len(waveform1),0.5E-3)
fft2 = numpy.fft.rfft(waveform2)
freqs2 = numpy.fft.rfftfreq(len(waveform2),0.5E-3)


plt.xlabel("kHz")
# plt.plot(freqs0,numpy.absolute(fft0),label="U Wires FFT",color="b")
plt.plot(freqs1,numpy.absolute(fft1),label="V Wires FFT",color="g")
plt.plot(freqs2,numpy.absolute(fft2),label="Y Wires FFT",color="r")
ax2.get_yaxis().set_ticks([])
plt.ylabel("Power")
plt.title("FFT of Correlated Waveforms")
plt.grid(True)
plt.legend(prop={'size':11})


plt.tight_layout()
fig.subplots_adjust(top=0.85)
fig.subplots_adjust(wspace=-0.2)

# Create a title for this plot:
save_title = "Correlated_Noise_Run_" + str(data[branches[branches_map['run']]])
save_title += "_Event_" +  str(data[branches[branches_map['event']]])


plt.savefig(save_title + ".png")
plt.savefig(save_title + ".pdf")


# plt.show()




# plt.legend(fontsize=20)


# b1 = plt.figure(figsize=(20,10))
# plt.xlabel("kHz")
# plt.plot(freqs0,numpy.absolute(fft0),label="U Wires FFT",color="b")
# plt.grid(True)
# plt.legend(fontsize=20)


# b2 = plt.figure(figsize=(20,10))
# plt.xlabel("kHz")
# plt.plot(freqs1,numpy.absolute(fft1),label="V Wires FFT",color="g")
# plt.grid(True)
# plt.legend(fontsize=20)


# b3 = plt.figure(figsize=(20,10))
# plt.xlabel("kHz")
# plt.plot(freqs2,numpy.absolute(fft2),label="Y Wires FFT",color="r")
# plt.grid(True)
# plt.legend(fontsize=20)







# plt.show()
# plt.savefig("mypng.png")

# # this next function applies the FFT to a numpy array
# fftvec = numpy.fft.rfft(wf_avg)
# # this next function correctly assigns the right scale to the x-axis
# # it takes the number of samples inputed to the FFT as the first argument
# # and the sampling time as the second argument.
# freqs  = numpy.fft.rfftfreq(len(wf_avg),0.5E-3)
# fig = plt.figure(figsize=(15,6))
# plt.plot(freqs,np.absolute(fftvec),color='r')
# plt.show()