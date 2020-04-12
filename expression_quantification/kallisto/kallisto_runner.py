###
### usage: time python kallisto_runner.py &> kallisto_runner_messages.txt
###

import sys,datetime,os

def kallisto_caller(label):

    printt('about to quantify {}'.format(label))

    sample_output_dir=results_dir+label
    executable='time kallisto quant'
    options=' -i {} -o {} --bias -t {} -b {} {} '.format(transcriptome_index,sample_output_dir,threads,boots,strand_flag)
    fastq_files=' '.join([clean_fastq_dir+label+'_L00{}_R1_clean.fastq'.format(i+1) + ' ' + clean_fastq_dir+label+'_L00{}_R2_clean.fastq'.format(i+1) for i in range(4)])
    command=executable+options+fastq_files

    print('')
    print(command)
    os.system(command)
    print('')

    return None

def printt(label):

    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(label)))

    return None

# 0. user defined variables
clean_fastq_dir='/Users/alomana/corona/cdi/data/clean_fastq/'
boots=100
threads=8
results_dir='/Users/alomana/corona/cdi/results/kallisto.{}/'.format(boots)
transcriptome_index='/Users/alomana/corona/cdi/data/homo_sapiens/transcriptome.idx'
strand_flag='--rf-stranded'

# 1. recover labels
printt('recover labels...')

all_files=os.listdir(clean_fastq_dir)
labels=sorted(list(set(([element.split('_L')[0] for element in all_files]))))

# 2. call kallisto quant
if os.path.exists(results_dir) == False:
    os.mkdir(results_dir)

for label in labels:
    kallisto_caller(label)
