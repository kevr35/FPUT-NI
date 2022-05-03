import numpy as np
import os
import random

#504:00:00
def SubmitJob(job_name, E, N, time,k,spacing,theta,executeable,dest, project = 'frgeeeph', run_time='12:00:00', mem_per_core='1', num_cores='1'):

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
   
   file.write("# Specify the error file name \n")
   file.writelines(["#$ -e ",job_name, "_ERROR \n \n"])
   
   file.write("# Specify the output file name \n")
   file.writelines(["#$ -o ",dest,job_name, " \n \n"])
   
   
   file.write("#Make the environment actually parallel\n")
   file.write("export OMP_NUM_THREADS=$NSLOTS\n\n")
   

   

   # ----- Assign Variables -----

   file.writelines(["N=", str(N), "\n \n"])
   file.writelines(["e=", str(E), "\n \n"])
   file.writelines(["time=", str(time), "\n \n"])
   file.writelines(["m=", str(k), "\n \n"])
   file.writelines(["spacing=", str(spacing), "\n \n"])
   file.writelines(["theta=", str(theta), "\n \n"])
   
   
   file.writelines(['echo -e \"$N\\n$e\\n$time\\n$m\\n$spacing\\n$theta\" |  "/usr4/ugrad/kevr/AverageOverPhases/"',executeable, '\n\n'])
   
   
   file.close()


   ###################################################
   # Submit the batch job and delete the bash script
   ###################################################

   os.system('qsub '+filename)
   os.system('rm '+filename)

N=31
Es=np.arange(.04,.2,.0025)
k=1

for E in Es:
  destination="/projectnb2/frgeeeph/FPUT_files/AlphaPhaseAveraging/N="+str(N)+"./E="+str(E)+"/"
  os.mkdir(destination)
  spacing=5000
  jobs=100
  
  time=5*10**6
  
  theta=np.pi/2 #All P initially
  job_name = 'TodaEntropy_E='+str(round(E,6))+'_N=' + str(N)#+"_spacing="+str(spacing)#+'_k=' + str(k)
  SubmitJob(job_name = job_name, E=E, N=N,spacing=spacing, theta=theta, time=time, k=k,executeable='TodaPhase',dest=destination)
  
  for job in range(jobs):
  
    job_name = 'AlphaEntropy_E='+str(round(E,6))+'_N=' + str(N)+"_"+str(job)
  
    theta=2*np.pi*random.random()      
    SubmitJob(job_name = job_name, E=E, N=N,spacing=spacing,theta=theta, time=time, k=k,executeable='AlphaPhases',dest=destination)  
  
