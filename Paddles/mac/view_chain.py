import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE(s)\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

import ROOT
from larlite import larlite as fmwk
from ROOT import flashana
from basictool.geoviewer import GeoViewer,geoalgo
import matplotlib.pyplot as plt
import numpy as np

plt.ion()
#viewer = GeoViewer(width=4,height=8)
#viewer._use_box = False

# Create ana_processor instance
my_proc = fmwk.ana_processor()

#my_proc.set_verbosity(0)

# Config file
cfg=sys.argv[1]
if not cfg.endswith('.fcl'):
    print 'Config file needs to end with \'.fcl\' extension (sorry bad joke)'
    sys.exit(1)

# Set input root file
for x in xrange(len(sys.argv)-2):
    my_proc.add_input_file(sys.argv[x+2])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kBOTH)

# Specify output root file name
#my_proc.set_ana_output_file("view_chain.root");

my_proc.set_output_file("view_chain.root")

# Attach an analysis unit ... here we use a base class which does nothing.
# Replace with your analysis unit if you wish.
tagger = fmwk.MuCSTagger()
tagger.configure(cfg)

ana = fmwk.MuCSOpStudy()
ana.configure(cfg)

mgr=ana.FlashMatchManager()
mgr.SetAlgo(flashana.NPtFilter())
mgr.SetAlgo(flashana.MaxNPEWindow())
mgr.SetAlgo(flashana.TimeCompatMatch())
mgr.Configure('flashmatch_mucs.fcl')

if not ana.track_producer() == tagger.producer():

    print 'Error: track producer is different between ana & tagger'
    sys.exit(1)

my_proc.add_process(tagger)
my_proc.add_process(ana)

my_proc.set_data_to_read(fmwk.data.kOpHit,ana.ophit_producer())
my_proc.set_data_to_read(fmwk.data.kTrack,ana.track_producer())
my_proc.set_data_to_read(fmwk.data.kCosmicTag,"MuCSTagger")
my_proc.set_data_to_read(fmwk.data.kOpFlash,ana.opflash_producer())

my_proc.set_data_to_write(fmwk.data.kOpHit,ana.ophit_producer())
my_proc.set_data_to_write(fmwk.data.kTrack,ana.track_producer())
my_proc.set_data_to_write(fmwk.data.kCosmicTag,"MuCSTagger")
my_proc.set_data_to_write(fmwk.data.kOpFlash,ana.opflash_producer())

print
print  "Finished configuring ana_processor. Start event loop!"
print

index=0
while my_proc.process_event():
    
    cos_score  = ana.get_ctag_score()
    print cos_score
    if ( cos_score != 0.5 and cos_score != 1.0): continue
    
    print 'event number:',index
    index +=1
    
    ophit_data = ana.get_ophit_flash().pe_v
    ophit_hypo = ana.get_ophit_hypo().pe_v
    ophit_veto = ana.get_veto()
    length = tagger.get_track_length()
    if (sum(ophit_data) == 0): continue
    ratio = float(sum(ophit_hypo))/float(sum(ophit_data))
    print ratio
    for x in range(0,31):
        if ophit_hypo[x]<7: ophit_hypo[x]=0
        if ophit_veto[x]:ophit_hypo[x]=0

    bins = np.arange(32)
    fig = plt.figure(figsize=(8,6))
    p1 = plt.plot(bins, ophit_data,  color = 'b',marker = '*')
    p2 = plt.plot(bins, ophit_hypo,  color = 'r',marker = 'o')
    plt.xlim(0,32)
    plt.text(25,10,'Hypo/OpReco:%0.3f\nTrack_Length:%0.1f'%(ratio,length),fontsize=15)
    plt.xlabel('PMT ID')
    plt.ylabel('PE')
    plt.legend( (p1[0], p2[0]), ('Data', 'Photon Library') )
    plt.grid()
    cos_score = 0.0
    usrinput = raw_input("Hit Enter: next evt  ||  q: exit viewer\n")
    if ( usrinput == "q"):
        sys.exit(0)
    # done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)
