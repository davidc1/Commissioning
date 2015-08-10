import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)


from ROOT import gSystem,TMath
from ROOT import larlite as fmwk
from ROOT import larutil

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
for x in xrange(len(sys.argv)-1):
    my_proc.add_input_file(sys.argv[x+1])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kBOTH)

# Specify analysis output root file name
my_proc.set_ana_output_file("");

# Specify data output root file name
my_proc.set_output_file("rawdigit_hits.root")

hitmaker = fmwk.RawHitMaker()
# set the number of ticks to use to measure baseline/rms
hitmaker.setNTicks(100)
# set the maximum RMS allowed for a channel to be 'good'
# if 'bad' hits will not be searched for on this channel
hitmaker.setRMSMax(10.0)
# set producer name for hits
hitmaker.setRawDigitProducer('daq')
# set threshold for a hit. hit window opened/closed when
# the baseline-subtracted ADCs are > than RMS * sigma cut
hitmaker.setSigmaCut(4.0)
# set minimum hit width
hitmaker.setMinHitWidth(2)
my_proc.add_process(hitmaker)

my_proc.set_data_to_write(fmwk.data.kHit,'rawhit')
#my_proc.set_data_to_write(fmwk.data.kRawDigit,'daq')



print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run()
#my_proc.run()

sys.exit()

