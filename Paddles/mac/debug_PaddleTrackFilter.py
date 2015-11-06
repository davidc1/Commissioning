import sys

from ROOT import *
#gSystem.Load("libMy_Repo_scratch_area.so")

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE(s)\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from larlite import larlite as fmwk

# prepare viewer for GeoAlgo objects
from ROOT import geoalgo
from basictool.geoviewer import GeoViewer
import matplotlib.pyplot as plt

# load viewer
plt.ion()
viewer = GeoViewer()
viewer._use_box = False

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
for x in xrange(len(sys.argv)-1):
    my_proc.add_input_file(sys.argv[x+1])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)

# Specify output root file name
my_proc.set_ana_output_file("PaddleTrackFilter_output.root");

# Attach an analysis unit ... here we use a base class which does nothing.
# Replace with your analysis unit if you wish.
my_QRA = fmwk.PaddleTrackFilter()
#my_QRA.SetSaveHistos(True)
my_proc.add_process(my_QRA)

print
print  "Finished configuring ana_processor. Start event loop!"
print

while my_proc.process_event():

    print 'new event...'
    fidVol = my_QRA.getFiducialVol()
    mcsTop = my_QRA.getTopMuCS()
    mcsBot = my_QRA.getBottomMuCS()
    mcsTrj = my_QRA.getTrj()
    mcsTrjCon = my_QRA.getTrj_con()
    mcsTrjMuCS = my_QRA.getTrj_mucs()
    mcsHalfLine = my_QRA.getHalfLine()
    mcsLineSeg  = my_QRA.getLineseg()
    # add objects to viewer
    viewer.clear()
    viewer.add(fidVol,'FV.','k')
    viewer.add(mcsTop,'Top','r')
    viewer.add(mcsBot,'Bot','m')
    
    #for x in range(0,len(mcsTrj)):
    for i in xrange(len(mcsTrj)):
        viewer.add(mcsTrj[i],"Out","b")
        viewer.add(mcsLineSeg[i],'Prj ','y')
        #viewer.add(mcsHalfLine[i],"halfline","y")
    for j in xrange(len(mcsTrjCon)):
        viewer.add(mcsTrjCon[j],"In","g")
        print "hehe"
    for k in xrange(len(mcsTrjMuCS)):
        viewer.add(mcsTrjMuCS[k],"Tr_mucs","r")        
        print "STOP!!!!!!!!!I am a tagged MUON !!!!!!!!"
#    for h in xrange(len(mcsHalfLine)): 
        #print h
 #      viewer.add(mcsHalfLine[h],"halfline","y")       


    viewer.show()

    # go to next event
    usrinput = raw_input("Hit Enter: next evt  ||  q: exit viewer\n")
    if ( usrinput == "q" ):
        sys.exit(0)

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)
