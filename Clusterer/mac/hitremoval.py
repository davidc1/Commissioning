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
my_proc.set_ana_output_file("hitremoval_ana.root");

# Specify data output root file name
my_proc.set_output_file(sys.argv[-1])

# prepare the various hit removal stages

_clusproducer   = "pandoraCosmic"
_vtxproducer    = "numuCC_vertex"
_hitproducer00  = "hit00"
_clusproducer01 = "cluster01"
_hitproducer02  = "hit02"
_clusproducer02 = "cluster02"
_hitproducer03  = "hit03"
_clusproducer03 = "cluster03"

def hitremoval00():

    algo = fmwk.LinearClusterSubsetRemoval()
    algo.setClusterProducer(_clusproducer)
    algo.setVertexProducer(_vtxproducer)
    algo.setOutHitProducer(_hitproducer00)

    n_hits = [  3, 10, 20, 30, 50, 70,100,130,180,230]
    lin_v  = [.02,.06,.10,.10,.10,.10,.10,.12,.20,.30]

    for i,n in enumerate(n_hits):
        algo.setMaxLinearity( lin_v[i] )
        algo.setMinNHits( n )

    return algo

def hitremoval01():

    algo = fmwk.SimpleClusterer()
    algo.setHitProducer(_hitproducer00)
    algo.setVtxProducer(_vtxproducer)
    algo.setOutClusProducer(_clusproducer01)
    algo.setRadius(1.0)
    algo.setCellSize(2.0)
    algo.setMaxHitRMS(19)
    algo.setUseVertex(True)
    algo.setVtxRadius(2.0)
    algo.setVerbose(False)
    algo.setMinTick(800)
    algo.setMaxTick(5445)

    return algo

def hitremoval02():

    algo = fmwk.LinearClusterRemoval()
    algo.setClusterProducer(_clusproducer01)
    algo.setOutHitProducer(_hitproducer02)
    algo.setOutClusterProducer(_clusproducer02)
    n_hits = [  3, 10, 20, 30, 50, 70,100,130,180,230]
    lin_v  = [.02,.10,.10,.13,.16,.20,.23,.30,.40,.45]
    for i,n in enumerate(n_hits):
        algo.setMaxLinearity( lin_v[i] )
        algo.setMinNHits( n )

    return algo

def hitremoval03():

    algo = fmwk.ClusterFilter()
    
    algo.setClusProducer(_clusproducer02)
    algo.setVtxProducer(_vtxproducer)
    algo.setOutHitProducer(_hitproducer03);
    algo.setOutClusterProducer(_clusproducer03);
    algo.setOutClusterAllProducer("clusterfilterall");
    algo.setMaxNHits(400)
    algo.setMaxArea(100*100)
    algo.setMaxDist(200)

    return algo

my_proc.add_process( hitremoval00() )
my_proc.add_process( hitremoval01() )
my_proc.add_process( hitremoval02() )
my_proc.add_process( hitremoval03() )

print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run()

sys.exit()

