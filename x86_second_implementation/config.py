
# names of the diagnostics
diagnostics = ['processExit']

# Location of infoli files
sim_dump_location = '/home/apostolis/Desktop/mpi_programs/ckpt_infoli/output/'

# simrun path
simrun_path = '/home/apostolis/Desktop/mpi_programs/ckpt_infoli/'

# Infoli kill script path - RELATIVE TO RCCERUN
killfoli_path = '/home/apostolis/Desktop/mpi_programs/ckpt_infoli/killfoli.sh'

# injector input files
##here was benchmarkInjectorFile =sim_dump_location + 'injectors/benchmark.txt'
infoliInjectorFile = '/home/apostolis/Desktop/mpi_programs/ckpt_infoli/injectors/infolijector.txt'
coreFailureInjectorFile = '/home/apostolis/Desktop/mpi_programs/ckpt_infoli/injectors/corefailjector.txt'
coreShutdownInjectorFile = '/home/apostolis/Desktop/mpi_programs/ckpt_infoli/injectors/coreshutjector.txt'
processExitInjectorFile = '/home/apostolis/Desktop/mpi_programs/ckpt_infoli/injectors/procexitjector.txt'

# True if running on a development environment rather than the SCC
devel = True

# max number of elements to utilize for the MTTF estimation
moving_avg_N = 50

# if False, only DUE checkpoints will be used.
use_SDC_checkpoints = True

## also these lines where here in the source folder of SCC
# Checkpointing latency for the target application - in seconds
#latency = 0.23

#Checkpoint interval optimization precision - selected intervals will be divisible by:
#prec_interv = 100
