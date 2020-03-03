###
### usage: time python read_cleaner.py &> messages.txt
###

import os,datetime,sys

def printt(label):

    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(label)))

    return None

def trimmomatic_caller(sample):

    executable='time java -jar {}trimmomatic-0.39.jar PE -threads {} -phred33 '.format(trimmomatic_path,number_threads)
    options=' ILLUMINACLIP:{}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'.format(adapter_file)

    input_files=os.listdir(raw_fastq_dir+samples[sample])
    for element in input_files:
        if 'R1_0' in element:
            input1=raw_fastq_dir+samples[sample]+'/'+element
        elif 'R2_0' in element:
            input2=raw_fastq_dir+samples[sample]+'/'+element
        else:
            printt('more files detected than expected. Exiting...')
            sys.exit()

    output1=clean_fastq_dir+sample+'_R1_clean.fastq'
    output2=clean_fastq_dir+sample+'_R2_clean.fastq'

    garbage1=clean_fastq_dir+sample+'_R1_garbage.fastq'
    garbage2=clean_fastq_dir+sample+'_R2_garbage.fastq'
        
    input_files=input1+' '+input2
    output_files=output1+' '+garbage1+' '+output2+' '+garbage2
    
    command=executable+input_files+' '+output_files+options

    printt('about to clean {}'.format(sample))
    print('')
    print(command)
    print('') 
    os.system(command)
    print('')

    return None

# 0. user defined variables
raw_fastq_dir='/Volumes/omics4tb2/Collaborations/CDI/FASTQ_Generation_2020-02-27_19_16_48Z-214660447/'
clean_fastq_dir='/Volumes/omics4tb2/alomana/projects/cdi/data/clean_fastq/'
trimmomatic_path='/Users/alomana/software/Trimmomatic-0.39/'
adapter_file=trimmomatic_path+'adapters/TruSeq3-PE-2.fa'
number_threads=8

# 1. recover samples
samples={}

folders=os.listdir(raw_fastq_dir)
folders.sort()

for folder in folders:
    sample_name=folder.split('-ds')[0]
    samples[sample_name]=folder

# 2. iterate Trimmomatic
for sample in samples:
    trimmomatic_caller(sample)
