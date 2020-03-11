import datetime, sys, os

def salmon_caller(label):

    printt('about to quantify {}'.format(label))

    sample_output_dir=results_dir+label
    executable='time salmon quant'
    flag1 = ' --index {} --libType A'.format(transcriptome_index)
    flag2 = ' --numBootstraps 100 '
    flag3 = '-1 ' + ' '.join([clean_fastq_dir+label+'_L00{}_R1_clean.fastq'.format(i+1) for i in range(4)])
    flag4 = ' -2 ' + ' '.join([clean_fastq_dir+label+'_L00{}_R2_clean.fastq'.format(i+1) for i in range(4)])
    flag5 = ' --threads {}'.format(threads)
    flag6 = ' --geneMap {}'.format(gtf_file)
    flag7 = ' --validateMappings --output {}'.format(sample_output_dir)
    command = executable + flag1 + flag2 + flag3 + flag4 + flag5 + flag6 + flag7
    

    print('')
    print(command)
    os.system(command)
    print('')

    return None

def printt(label):

    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(label)))

    return None

# 0. user defined variables
clean_fastq_dir = '/Volumes/omics4tb2/alomana/projects/cdi/data/clean_fastq/'
threads = 8
results_dir = '/Volumes/omics4tb2/alomana/projects/cdi/results/salmon_madrid/'
transcriptome_index = '/Volumes/omics4tb2/alomana/projects/cdi/data/mmrf_annotation/salmon.index'
gtf_file = '/Volumes/omics4tb2/alomana/projects/cdi/data/mmrf_annotation/Homo_sapiens.GRCh37.74.gtf.hs37d5.EGFRvIII.gtf'

# 1. recover labels
printt('recover labels...')

all_files=os.listdir(clean_fastq_dir)
labels=sorted(list(set(([element.split('_L')[0] for element in all_files]))))

# 2. call kallisto quant
if os.path.exists(results_dir) == False:
    os.mkdir(results_dir)

for label in labels:
    salmon_caller(label)
