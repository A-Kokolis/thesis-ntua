# path to main.c program for voltage alterations. Needs the frequency devisor as a parameter
freqChangeCommand='cpufreq-set -r -f %d'

# path to ckpTime files that indicated when a checkpoint was taken. Needs as a parameter the simulation step we are restarting from
ckpTimeFile='/home/apostolis/Desktop/mpi_programs/ckpt_infoli/ckpTime/%d.txt'

#path to applicationRestartTimeFile which includes the restart time overhead of the application. (The time to restore its parameters)
applicationRestartTimeFile='/home/apostolis/Desktop/mpi_programs/ckpt_infoli/ckpTime/applicationRestartTime.txt'

#path to ckptOverheadTimeFile which includes the time overhead to store a checkpoint
ckptOverheadTimeFile='/home/apostolis/Desktop/mpi_programs/ckpt_infoli/ckpTime/ckpTime.txt'

#the default multiplier refers to the cpufrequency we wanted the app to run
defaultMultiplier=1
defaultFreq=1200000

#multiplierList depends on the default multiplier we select
multiplierList=[0.8,0.9,1,1.12,1.27,1.47,1.57]

#available frequencies list (check /sys/devices/system/cpu/cpu*/cpufreq/scaling_available_frequencies for your available frequencies)
freqList=[800000,1000000,1200000,1400000,1600000,1800000,2000000]



