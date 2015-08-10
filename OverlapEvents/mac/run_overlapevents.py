import os,sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)


from ROOT import gSystem,TMath
from ROOT import larlite as fmwk
from ROOT import larutil
#from larlite import larlite as fmwk
#from larlite import larutil

# Set input root file
for x in xrange(len(sys.argv)-1):
    
    # Create ana_processor instance
    my_proc = fmwk.ana_processor()

    my_proc.add_input_file(sys.argv[x+1])

    fname = sys.argv[x+1]

    run    = int(fname[fname.find('_run_')+5:fname.find('_subrun_')])
    subrun = int(fname[fname.find('_subrun_')+8:fname.find('.root')])

    if (run <= 1236):
        continue

    print str(run)+' '+str(subrun)

    outfile = "rawdigit_laser_%05d_%04d.root"%(run,subrun)
    
    found = False
    ldir = os.listdir('./')
    for f in ldir:
        if (os.path.isfile(f) == True):
            if (f == outfile):
                found = True
                break
    if (found == True):
        continue

    print 'continuing...'
    # Specify IO mode
    my_proc.set_io_mode(fmwk.storage_manager.kBOTH)

    # Specify analysis output root file name
    #my_proc.set_ana_output_file("laser_finder.root");

    # Specify data output root file name
    my_proc.set_output_file(outfile)

    hitmaker = fmwk.OverlapEvents()

    my_proc.set_data_to_write(fmwk.data.kRawDigit,'laser')
    my_proc.add_process(hitmaker)

    print
    print  "Finished configuring ana_processor. Start event loop!"
    print

    my_proc.run()

sys.exit()

