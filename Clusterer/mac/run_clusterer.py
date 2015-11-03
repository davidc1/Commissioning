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



# Set input root file
for x in xrange(len(sys.argv)-1):

    # Create ana_processor instance
    my_proc = fmwk.ana_processor()

    fname = sys.argv[x+1]

    run = int(fname[fname.find('run_')+4:fname.find('run_')+11])
    subrun = int(fname[fname.find('subrun_')+7:fname.find('subrun_')+12])

    print 'run: %i, subrun: %i'%(run,subrun)

    my_proc.add_input_file(fname)

    # Specify IO mode
    my_proc.set_io_mode(fmwk.storage_manager.kBOTH)

    # Specify analysis output root file name
    my_proc.set_ana_output_file("");

    # Specify data output root file name
    my_proc.set_output_file("larlite_rawclusters_run_%05i_subrun_%05i.root"%(run,subrun))

    clusterer = fmwk.SimpleClusterer()
    #clusterer.setHitProducer('rawhit')
    #clusterer.setHitProducer('cchit')
    clusterer.setHitProducer('cchit')
    clusterer.setVerbose(True)
    
    my_proc.add_process(clusterer)
    
    my_proc.set_data_to_write(fmwk.data.kHit,'cchit')
    my_proc.set_data_to_write(fmwk.data.kCluster,'rawcluster')
    my_proc.set_data_to_write(fmwk.data.kAssociation,'rawcluster')



    print
    print  "Finished configuring ana_processor. Start event loop!"
    print
    
    my_proc.run()

sys.exit()

