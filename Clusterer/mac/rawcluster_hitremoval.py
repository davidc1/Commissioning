import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)


from larlite import larlite as fmwk

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
for x in xrange(len(sys.argv)-2):
    fname = sys.argv[x+1]
    my_proc.add_input_file(fname)
    
# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kBOTH)

# Specify analysis output root file name
my_proc.set_ana_output_file("");

# Specify data output root file name
my_proc.set_output_file(sys.argv[-1])

#hitproducer = 'shrhits2'
#hitproducer = 'gaushit'
hitproducer = 'shrhits'

clusterer = fmwk.SimpleClusterer()
clusterer.setHitProducer(hitproducer)
clusterer.setRadius(1.0)
clusterer.setCellSize(2.0)
clusterer.setMaxHitRMS(19)
clusterer.setUseVertex(True)
#clusterer.setHitProducer("gaushit")
clusterer.setVtxProducer("numuCC_vertex")
clusterer.setOutClusProducer("rawcluster")
clusterer.setVtxRadius(2.0)
clusterer.setVerbose(False)
clusterer.setMinTick(800)
clusterer.setMaxTick(5445)

my_proc.add_process(clusterer)

#my_proc.set_data_to_write(fmwk.data.kHit,hitproducer)
#my_proc.set_data_to_write(fmwk.data.kCluster,'rawcluster')
#my_proc.set_data_to_write(fmwk.data.kAssociation,'rawcluster')

print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run()

sys.exit()

