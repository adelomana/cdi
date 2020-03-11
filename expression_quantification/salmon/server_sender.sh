#!/bin/bash

#$ -N salmon
#$ -o /proj/omics4tb2/alomana/scratch/messages.salmon.o.txt
#$ -e /proj/omics4tb2/alomana/scratch/messages.salmon.e.txt
#$ -pe smp 64
#$ -S /bin/bash

cd /users/alomana
source .bash_profile 

cd /proj/omics4tb2/alomana/projects/cdi/deployment
time python caller.py
