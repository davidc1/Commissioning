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
my_proc.set_ana_output_file("linearclusterlocalremoval.root");

# Specify data output root file name
my_proc.set_output_file(sys.argv[-1])

clusterproducer = 'shrcluster'

clusterer = fmwk.LinearClusterLocalRemoval()
clusterer.setClusterProducer(clusterproducer)
clusterer.setMaxLinearity( -4.5 )
clusterer.setMinNHits(50)
clusterer.setMaxRadius(7.0)
my_proc.add_process(clusterer)

#my_proc.set_data_to_write(fmwk.data.kHit,hitproducer)
#my_proc.set_data_to_write(fmwk.data.kCluster,'shrcluster')
#my_proc.set_data_to_write(fmwk.data.kAssociation,'shrcluster')
#my_proc.set_data_to_write(fmwk.data.kHit,"shrhits")

print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run(0,20)

sys.exit()

