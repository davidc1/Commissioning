import sys
from ROOT import gSystem
gSystem.Load("libCommissioning_MakeWFTTree")
from ROOT import sample

try:

    print "PyROOT recognized your class %s" % str(sample)

except NameError:

    print "Failed importing MakeWFTTree..."

sys.exit(0)

