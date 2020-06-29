###
### This script filters DEGs following the rules: (1) absolute log2 fold-change greater than one compared to initial conditions; (2) expression greater than 3 TPM; (3) relative standard error of the mean across replicates lower than one third; and (4) P < 0.05 and P-adjusted < 0.1 as called by DESeq2 [1] version 1.26.0.
###

import sys, numpy

###
### FUNCS
###

def metadata_reader():

    metadata = []
    with open(metadata_file, 'r') as f:
        next(f)
        for line in f:
            v = line.split('\t')
            sample = v[0]
            metadata.append(sample)
    
    return metadata

def read_DEGs(DEGs, comparison, trend): 

    '''
    Returns a dictionary of DEGs.
    '''

    working_file = DESeq2_folder + comparison + '_' + trend + '.tsv'
    
    with open(working_file, 'r') as f:
        next(f)
        for line in f:
            v = line.split('\t')
            
            ensembl = v[0]
            gene_name = v[1]
            biotype = v[2]
            description = v[3]
            basemean = float(v[4])
            logFC = float(v[5])
            pvalue = float(v[6])
            adjusted = float(v[7])

            info = (ensembl, gene_name, biotype, description, basemean, logFC, pvalue, adjusted)

            DEGs[comparison][trend].append(info)
                    
    print('{} {} \t DEGs: {}'.format(comparison, trend, len(DEGs[comparison][trend])))

    return DEGs

###
### MAIN
###

#
# 0. use-defined variables
#

DESeq2_folder = '/Users/alopez/projects_isb/cdi/results/deseq2/unfiltered/'
metadata_file = '/Users/alopez/projects_isb/cdi/data/metadata/metadata.txt'
expression_file = DESeq2_folder + 'DESeq2_TPM_values.tsv'
filtered_folder = '/Users/alopez/projects_isb/cdi/results/deseq2/filtered/'

comparisons = {}
comparisons['co_vs_mono_time_zero'] = [(13,14,15), (1,2,3)]

comparisons['mono_time_one_vs_time_zero'] = [(4,5,6), (1,2,3)]
comparisons['mono_time_four_vs_time_zero'] = [(7,8,9), (1,2,3)]
comparisons['mono_time_twentyfour_vs_time_zero'] = [(10,11,12), (1,2,3)]

comparisons['co_time_one_vs_time_zero'] = [(16,17,18), (13,14,15)]
comparisons['co_time_four_vs_time_zero'] = [(19,20,21), (13,14,15)]
comparisons['co_time_twentyfour_vs_time_zero'] = [(22,23,24), (13,14,15)]

trend_tags = ['up', 'down']

expression_threshold = 2
discrete_fc_threshold = 1
noise_threshold = 1/2

print('expression threshold: > {}'.format(expression_threshold))
print('discrete log2FC: > {}'.format(discrete_fc_threshold))
print('noise threshold: > {}'.format(noise_threshold))


#
# 1. read data
#

# 1.1. define the DEGs across experimental design
print('define DEGs')
DEGs = {}

for comparison in comparisons:
    DEGs[comparison] = {}
    for trend in trend_tags:
        DEGs[comparison][trend] = []
        DEGs = read_DEGs(DEGs, comparison, trend)

# 1.2. define metadata
print('define metadata')
sample_IDs = metadata_reader()
                
# 1.2. define expression
print('define expression')

expression = {}
for sample in sample_IDs:
    expression[sample] = {}
    
with open(expression_file, 'r') as f:
    next(f)
    for line in f:
        v = line.split('\t')
        
        gene_name = v[0]
        v = [float(element) for element in v[1:]]

        for i in range(len(v)):
            value = v[i]
            sampleID = sample_IDs[i]
            expression[sampleID][gene_name] = value

#
# 2. analysis
#

# 2.1. filter
print('filter DEGs')

for comparison in DEGs:
    for trend in trend_tags:                

        # define sample and reference sample labels
        sample_labels = [sample_IDs[index-1] for index in comparisons[comparison][0]]
        reference_labels = [sample_IDs[index-1] for index in comparisons[comparison][1]]
        print(sample_labels, reference_labels)

        # filters
        container = []
        before = len(DEGs[comparison][trend])
        for case in DEGs[comparison][trend]:
            including = True
            ensembl = case[0]

            # gather TPM expression
            ref = [expression[label][ensembl] for label in reference_labels]
            sam = [expression[label][ensembl] for label in sample_labels]

            # filter 1: identify low-expressed genes
            r = numpy.median(ref); s = numpy.median(sam)
            top = numpy.max([r, s])

            # filter 2: identify fold-changes using discrete values
            ###
            ###            [round(x, epsilon)/epsilon ] + 1 
            ###  FC = abs  -------------------------------- > 1
            ###            [round(y, epsilon)/epsilon ] + 1
            ###
            ###
            ###  epsilon = 1
            num = numpy.around(s) + 1
            den = numpy.around(r) + 1
            fc = num/den
            abs_log2FC = numpy.abs(numpy.log2(fc))

            # filter 3: identify noisy genes
            ref_int = numpy.around(ref) + 1
            sam_int = numpy.around(sam) + 1
            sem_ref = numpy.std(ref_int) / numpy.sqrt(len(ref_int))
            rsem_ref = sem_ref / numpy.mean(ref_int)
            sem_sam = numpy.std(sam_int) / numpy.sqrt(len(sam_int))
            rsem_sam = sem_sam / numpy.mean(sam_int)
            noise = numpy.max([rsem_ref, rsem_sam])

            # selection
            if abs_log2FC < discrete_fc_threshold:
                including = False
                info = 'WARNING: small change gene discarded. Expression changes from {:.3f} ({}) to {:.3f} ({}), resulting in abs_log2FC {:.3f}. {}, {}, {}'.format(r, den, s, num, abs_log2FC, case[0], case[1], case[3])
                print(info)

            if (including == True) and (top < expression_threshold):
                including = False
                info = 'WARNING: low-expression gene discarded. Expression changes from {:.3f} to {:.3f}. {}, {}, {}'.format(r, s, case[0], case[1], case[3])
                print(info)

            if (including == True) and (noise > noise_threshold):
                including = False
                info = 'WARNING: noisy gene. Ref: {}, RSEM {:.3f}. Sam: {}, RSEM {:.3f}. {}, {}, {}'.format(ref, rsem_ref, sam, rsem_sam, case[0], case[1], case[3])
                print(info)

            if including == True:
                content = list(case)
                content.append(r); content.append(s); content.append(abs_log2FC)
                container.append(content)

        # info about low-expression cases
        after = len(container)
        print('{} {} | DEGs final reduction from {} to {}'.format(comparison, trend, before, after))
        print('')

        # store reduced version
        DEGs[comparison][trend] = container

#
# 3. write results
#

# 3.1. store
for comparison in comparisons:
    for trend in trend_tags:
        storage = filtered_folder + '{}_{}_filtered.tsv'.format(comparison, trend)
        with open(storage, 'w') as f:
            f.write('ENSEMBL\tGene name\tBiotype\tDescription\tBase mean\tlog2FC\tP value\tAdjusted P-value\tReference expression (TPM)\tSample expression (TPM)\tDiscrete abs(log2FC)\n')
            for content in DEGs[comparison][trend]:
                line = ''
                for element in content:
                    if isinstance(element, str) == False:
                        if numpy.abs(element) < 0.05 and element != 0.:
                            sub = '{:.4E}'.format(element)
                        else:
                            sub = '{:.4f}'.format(element)
                    else:
                        sub = element
                    line=line+sub+'\t'
                line=line+'\n'
                f.write(line)
