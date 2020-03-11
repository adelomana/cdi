#!/bin/bash

ENST_TO_FILTER='/Volumes/omics4tb2/alomana/projects/cdi/data/mmrf_annotation/info/Homo_sapiens_GRCh37_74_ig_ENST_to_filter_out.tsv'
ENSG_TO_ENST_MAP='/Volumes/omics4tb2/alomana/projects/cdi/data/mmrf_annotation/info/Homo_sapiens_GRCh37_74_ENSG_to_ENST_map.txt'

# Change directory to salmon results
cd '/Volumes/omics4tb2/alomana/projects/cdi/results/salmon'
pwd

ALLSAMPLES="`echo \
				01_MM1S_0hr_R1 \
				02_MM1S_0hr_R2 \
				03_MM1S_0hr_R3 \
				04_MM1S_1hr_R1 \
				05_MM1S_1hr_R2 \
				06_MM1S_1hr_R3 \
				07_MM1S_4hr_R1 \
				08_MM1S_4hr_R2 \
				09_MM1S_4hr_R3 \
				10_MM1S_24hr_R1 \
				11_MM1S_24hr_R2 \
				12_MM1S_24hr_R3 \
				13_hFOB-MM1S_0hr_R1 \
				14_hFOB-MM1S_0hr_R2 \
				15_hFOB-MM1S_0hr_R3 \
				16_hFOB-MM1S_1hr_R1 \
				17_hFOB-MM1S_1hr_R2 \
				18_hFOB-MM1S_1hr_R3 \
				19_hFOB-MM1S_4hr_R1 \
				20_hFOB-MM1S_4hr_R2 \
				21_hFOB-MM1S_4hr_R3 \
				22_hFOB-MM1S_24hr_R1 \
				23_hFOB-MM1S_24hr_R2 \
				24_hFOB-MM1S_24hr_R3`"

for label in `echo $ALLSAMPLES` 
do
    SAMPLEDIR=/Volumes/omics4tb2/alomana/projects/cdi/results/salmon/${label}
    
    # Change directory to salmon results and rename files
    cd ${SAMPLEDIR}
    pwd
    cp quant.genes.sf ${label}_salmon_e74GTF_genes.sf
    cp quant.sf ${label}_salmon_e74GTF_transcripts.sf

    # Remove unwanted ENST
    awk -F '\t' 'FNR == NR {a[$1]=$1;next} !($1 in a) { OFS = "\t" ; print $0 }' ${ENST_TO_FILTER} ${label}_salmon_e74GTF_transcripts.sf > filtered_quant.sf

    # Calculate SUM of number of reads divided by Effective length for TPM calculation
    TPMSUM=`awk -F'\t' 'BEGIN { SUM = 0 } NR > 1 { SUM = SUM + ($5/$3)} END { print SUM }' filtered_quant.sf`
    echo $TPMSUM
    awk -F'\t' -v TPMSUM=${TPMSUM} 'NR == 1 { OFS = "\t" ; print $0 } ; NR > 1 { OFS = "\t" ; $4 = 1000000 * ($5 / $3 ) / TPMSUM ; print $0 } ' filtered_quant.sf > filtered_recalculated_quant.sf

    # Map each ENST to it's ENSG
    awk -F'\t' 'BEGIN{ OFS = "\t" ; print "Name","Gene","TPM" } FNR==NR{a[$2]=$1;next} ($1 in a) { OFS = "\t"; print $1,a[$1],$4 }' ${ENSG_TO_ENST_MAP} filtered_recalculated_quant.sf > filtered_recalculated_quant_with_ENSG.sf

    # Sum the transcript TPM's for each gene to create final filtered Gene level TPM file
    awk -F'\t' 'BEGIN{OFS = "\t" ; print "Name","TPM"} NR> 1 { a[$2]+=$3 } END { for (x in a) { OFS = "\t" ; print x,a[x] }}' filtered_recalculated_quant_with_ENSG.sf > ${label}_salmon_e74GTF_filtered_genes.sf

    # make final name of filtered transcript file
    cp filtered_recalculated_quant.sf ${label}_salmon_e74GTF_filtered_transcripts.sf
    
done
