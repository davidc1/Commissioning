import sys,os

for i in xrange(200):

    #fin = '/home/david/uboone/data/mcc7/bnb1pi0/larlite_mcc7_bnb_1pi0_%04i.root'%i
    fin = '/home/david/uboone/data/mcc7/v05_08/gamma/larlite_reco2d_%04i.root'%i

    if not (os.path.isfile(fin)):
        print '%s not a valid file'%fin 
        continue

    #fout = '/home/david/uboone/data/mcc7/bnb1pi0/larlite_rawclusters_%04i.root'%i
    fout = '/home/david/uboone/data/mcc7/v05_08/gamma/larlite_rawclusters_%04i.root'%i

    cmd = 'python clusterer.py %s %s'%(fin,fout)
    os.system(cmd)
