import sys,os

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from ROOT import gSystem
#gSystem.Load("libSimpleWFAna")
from ROOT import larlite as fmwk

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
run=int(sys.argv[1])
DATA_DIR="/home/david/uBooNE/data/WireBias/"
temp_files=[x for x in os.listdir(DATA_DIR) if x.find('LArLite-%03d-' % run)==0]
files={}
for f in temp_files:
    subrun=f.replace('LArLite-%03d-' % run,'')
    subrun=int(subrun.replace('.root',''))
    files[subrun]=f

subrun_start = files.keys()[0]
subrun_end   = files.keys()[-1]
if len(sys.argv)==4:
    subrun_start = int(sys.argv[2])
    subrun_end   = int(sys.argv[3])
    my_proc.set_ana_output_file("corr_run%03d_subrun%03d_%03d.root" % (run,subrun_start,subrun_end))
else:
    my_proc.set_ana_output_file("corr_run%03d.root" % run);

for x in files.keys():
    if x >= subrun_start and x<= subrun_end:
        my_proc.add_input_file('%s/%s' % (DATA_DIR,files[x]))

#my_proc.set_input_rootdir("scanner")
# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)

# Specify output root file name
my_proc.set_ana_output_file("noise_zigzag_run%03d.root" % run);

corr_unit = fmwk.ZigZag()
#corr_unit.setVerbose(True)
corr_unit.setZigZagMinLength(4)

contents = open('mapping_v01.txt','r')
crate_slot_map={}
for line in contents:

    if len(line.split()) < 4: continue
    (larch,crate,slot,ch) = [int(x) for x in line.split()]
    corr_unit.SetMap(larch,crate,slot,ch)
    if not crate in crate_slot_map.keys():
        crate_slot_map[crate]={}
    if not slot in crate_slot_map[crate].keys():
        crate_slot_map[crate][slot]={}
    if not ch in crate_slot_map[crate][slot].keys():
        crate_slot_map[crate][slot][ch]=larch
    

#reference = open('%s/noise_ch.txt' % os.environ['HOME'],'r').read().split('\n')
reference = open("zigzag_small.txt",'r')
for line in reference:
    chans = line.split(",")
    for l in chans:
        print l
        corr_unit.add(int(l))

my_proc.add_process(corr_unit)
#my_proc.add_process(fmwk.SimpleWFAna());

#my_proc.set_verbosity(fmwk.MSG.DEBUG)
print
print  "Finished configuring ana_processor. Start event loop!"
print

# Let's run it.
my_proc.run(5,1)

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)

