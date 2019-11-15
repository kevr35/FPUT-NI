import numpy as np
import os


def SubmitJob(job_name, E, N, time, project = 'frgeeeph', run_time='48:00:00', mem_per_core='8', num_cores='8'):

   ###############################################
   # Create bash script used for the batch job
   ###############################################

   filename = job_name + '.sh'
   file = open(filename,"w+")

   file.write("#!/bin/bash -l \n \n" )

   file.write("# Set SCC Project \n")
   file.write("#$ -P "+project+" \n \n")

   file.write("# Specify hard time limit for the job. \n")
   file.write("#   The job will be aborted if it runs longer than this time. \n")
   file.write("#$ -l h_rt="+run_time+" \n \n")

   file.write("# Send an email when the job finishes or if it is aborted (by default no email is sent). \n")
   file.write("#$ -m ea \n \n")
    
   file.write("# Combine output and error files into a single file \n")
   file.write("#$ -j y \n \n")

   file.write("# Request a node with at least 4 GB of memory per core \n")
   file.write("#$ -l mem_per_core="+mem_per_core+"G \n \n")

   file.write("# Request a paralell environemtn with _ cores \n")
   file.write("#$ -pe omp "+num_cores+" \n \n")
   
   #file.write("# Request to run on Engineering resources only \n")
   #file.write("#$ -l buyin \n \n")

   #file.write("# Number of GPUs per CPU requestedz \n")
   #file.write("#$ -l gpus=0.25 \n \n")

   #file.write("# Minimum GPU Capability \n")
   #file.write("#$ -l gpu_c=3.5 \n \n")
    
   # ----- Assign Job Name -----
   file.write("# Give job a name \n")
   file.writelines(["#$ -N ", job_name, "\n \n"])

   file.write("# Specify the output file name \n")
   file.writelines(["#$ -o ", job_name, "\n \n"])

   # ----- Assign Variables -----

   file.writelines(["N=", str(N), "\n \n"])
   file.writelines(["e=", str(E), "\n \n"])
   file.writelines(["time=", str(time), "\n \n"])

   file.writelines(['echo -e \"$N\\n$e\\n$time\" |  ./pbeta \n'])

   file.close()


   ###################################################
   # Submit the batch job and delete the bash script
   ###################################################

   os.system('qsub '+filename)
   os.system('rm '+filename)

Sizes=[31]
Energies=[1]


for N in Sizes:
   for E in Energies:

      job_name = 'PBeta_EB='+str(round(E,3))+'_N=' + str(N)
      S = E* (N+1)
      
      if S < 10:
        time=(0.2024+0.0227*S)*(N+1)**3
      else:
        time = (0.5862 / S**(1/2))*(N+1)**3
      

      SubmitJob(job_name = job_name, E=E, N=N, time=time)