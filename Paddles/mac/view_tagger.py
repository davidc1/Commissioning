import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE(s)\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from larlite import larlite as fmwk

# prepare viewer for GeoAlgo objects
#from ROOT import geoalgo
from basictool.geoviewer import GeoViewer,geoalgo
import matplotlib.pyplot as plt

# load viewer
plt.ion()
viewer = GeoViewer(width=4,height=8)
viewer._use_box = False

# Create ana_processor instance
my_proc = fmwk.ana_processor()

my_proc.set_verbosity(0)

# Config file
cfg=sys.argv[1]
if not cfg.endswith('.fcl'):
    print 'Config file needs to end with \'.fcl\' extension (sorry bad joke)'
    sys.exit(1)

# Set input root file
for x in xrange(len(sys.argv)-2):
    my_proc.add_input_file(sys.argv[x+2])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)

# Specify output root file name
my_proc.set_ana_output_file("PaddleTrackFilter_output.root");

# Attach an analysis unit ... here we use a base class which does nothing.
# Replace with your analysis unit if you wish.
myunit = fmwk.MuCSTagger()
myunit.configure(cfg)
my_proc.add_process(myunit)

print
print  "Finished configuring ana_processor. Start event loop!"
print
index=0
while my_proc.process_event():
    cos_score  = myunit.get_ctag_score()
    print 'cos_score is ',cos_score
    if ( cos_score != 0.5 and cos_score != 1.0): continue
        
    # add objects to viewer
    viewer.clear()
    print index
    index +=1
    print 'new event...'
    fidVol      = geoalgo.AABox(0,-116,0,253,116,1060)
    mcsTop      = myunit.upper_box()
    mcsBot      = myunit.lower_box()
    mcsupppt    = myunit.upper_pt()
    mcslowpt    = myunit.lower_pt()
    mcsTrj_v    = myunit.matched_trajectory()
    mcsDir_v    = myunit.matched_dir()

    viewer.add(fidVol,'FV.','k')
    viewer.add(mcsTop,'Top','r')
    viewer.add(mcsBot,'Bot','m')
    
    #for x in range(0,len(mcsTrj)):
    for i in xrange(len(mcsTrj_v)):
        print 'traj',i
        viewer.add(mcsTrj_v[i],"Traj","b")
        viewer.add(mcsupppt[i],'u','r')
        #viewer.add(mcslowpt[i],'l','m')   //enable when _hit_both_box == true
    for j in xrange(len(mcsDir_v)):
        print 'line',j
        viewer.add(mcsDir_v[j],"Dir","g")
    viewer.construct()
    viewer._ax.set_xlim(-150,300)
    viewer._ax.set_ylim(-120,350)
    viewer._ax.set_zlim(-10,1070)
    viewer._ax._axis3don = False
    viewer._fig.canvas.draw()#plt.show()

    # go to next event
    usrinput = raw_input("Hit Enter: next evt  ||  q: exit viewer\n")
    if ( usrinput == "q" ):
        sys.exit(0)

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)
