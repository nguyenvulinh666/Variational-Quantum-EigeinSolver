#!/bin/bash

#specifies the interpreting shell for this job to be the Bash shell.
#$ -S /bin/bash

# The batch system should use the current directory as working directory.
#$ -cwd

# Name your job. Unless you use the -o and -e options, output will
# go to a unique file name.ojob_id for each job.
#$ -N Jobname

# Merge the standard error stream into the standard output stream
# instead of having two separate error and output streams.
#$ -j y

# Specify the hostname run
#$ -l h=compute-0-3

# Use the parallel environment "smp" with X threads.
#$ -pe smp 1

# Use reservation to stop starvation.
# To use parameter -R y, this will turn on the reservation of slots,
# i.e. you won't be jumped by processes requesting less slots.
#$ -R y

## <put commands to run your job here>
python run_VQE_modified_on_HPC_sample.py
