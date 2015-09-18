import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from ROOT import *


# Create ana_processor instance
my_proc = larlite.ana_processor()

# Set input root file
for x in xrange(len(sys.argv)-1):
    my_proc.add_input_file(sys.argv[x+1])


# Specify IO mode
my_proc.set_io_mode(larlite.storage_manager.kREAD)

# Specify analysis output root file name
# my_proc.set_ana_output_file("showerRecoUboone_ana.root")
# Specify data output root file name
# my_proc.set_output_file("noneOutput.root")



ana_unit=evd.DrawRawDigit()

ana_unit.SetCorrectData(True)
ana_unit.SetSaveData(True)
ana_unit.SetStepSizeByPlane(2*48,0)
ana_unit.SetStepSizeByPlane(2*48,1)
ana_unit.SetStepSizeByPlane(2*96,2)

my_proc.add_process(ana_unit)

# Add an ana unit to do the shower quality:


print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run(0,10)
# my_proc.process_event(2)