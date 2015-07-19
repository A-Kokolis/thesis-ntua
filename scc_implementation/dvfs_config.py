# path to main.c program for voltage alterations. Needs the frequency devisor as a parameter
VoltageChangePath='/shared/apostolis/brain/main %d'

# path to ckpTime files that indicated when a checkpoint was taken. Needs as a parameter the simulation step we are restarting from
ckpTimeFile='/shared/apostolis/brain/ckpTime/%d.txt'

#path to applicationRestartTimeFile which includes the restart time overhead of the application. (The time to restore its parameters)
applicationRestartTimeFile='/shared/apostolis/brain/ckpTime/applicationRestartTime.txt'

#path to ckptOverheadTimeFile which includes the time overhead to store a checkpoint
ckptOverheadTimeFile='/shared/apostolis/brain/ckpTime/ckpTime.txt'



