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
my_proc.set_ana_output_file("mucs_tagger_ana.root");

my_proc.set_output_file("mucs_tagger.root")

# Attach an analysis unit ... here we use a base class which does nothing.
# Replace with your analysis unit if you wish.
tagger = fmwk.MuCSTagger.GetME()
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

# Let's run it.
my_proc.run();

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)
