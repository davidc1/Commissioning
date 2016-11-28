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
my_proc.set_ana_output_file("linearclusterremoval.root");

# Specify data output root file name
my_proc.set_output_file(sys.argv[-1])

clusterer = fmwk.PhotonClusterer()
clusterer.setClusterProducer('rawcluster')
clusterer.setNMin(2)
clusterer.setNMax(100)
clusterer.setQMin(60)
n_hits = [  3, 10, 20, 30, 50, 70, 100, 130, 180, 230]
lin_v  = [.02,.10,.13,.14,.15,.16, .17, .20, .20, .30]

for i,n in enumerate(n_hits):
    clusterer.setMaxLinearity( lin_v[i] )
    clusterer.setMinNHits( n )

my_proc.add_process(clusterer)

#my_proc.set_data_to_write(fmwk.data.kHit,hitproducer)
my_proc.set_data_to_write(fmwk.data.kCluster,'photon')
my_proc.set_data_to_write(fmwk.data.kAssociation,'photon')
#my_proc.set_data_to_write(fmwk.data.kHit,"shrhits")

print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run()

sys.exit()

