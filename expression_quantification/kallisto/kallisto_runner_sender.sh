#!/bin/bash

#$ -N cdi640
#$ -o /proj/omics4tb2/alomana/scratch/messages.cdi.kallisto.b640.o.txt
#$ -e /proj/omics4tb2/alomana/scratch/messages.cdi.kallisto.b640.e.txt
#$ -pe smp 64
#$ -S /bin/bash

cd /users/alomana
source .bash_profile 

echo "starting..."
cd /proj/omics4tb2/alomana/projects/cdi/deployment
time python kallisto_runner.640.py &> kallisto_runner_messages.b640.txt
echo "... run completed."
