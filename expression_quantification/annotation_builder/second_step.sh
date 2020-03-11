#!/bin/bash

./gtf_to_fasta /proj/omics4tb2/alomana/projects/cdi/data/mmrf_annotation/Homo_sapiens.GRCh37.74.gtf.hs37d5.EGFRvIII.gtf /proj/omics4tb2/alomana/projects/cdi/data/mmrf_annotation/hs37d5_plusRibo_plusOncoViruses_plusERCC.fa /proj/omics4tb2/alomana/projects/cdi/data/mmrf_annotation/Homo_sapiens.GRCh37.74.gtf.hs37d5.EGFRvIII.fa

sed 's/^>[0-9]* />/g' /proj/omics4tb2/alomana/projects/cdi/data/mmrf_annotation/Homo_sapiens.GRCh37.74.gtf.hs37d5.EGFRvIII.fa > /proj/omics4tb2/alomana/projects/cdi/data/mmrf_annotation/Homo_sapiens_number_fix.GRCh37.74.gtf.hs37d5.EGFRvIII.fa
