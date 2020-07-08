import sys,os

# 0. user defined variables
kallisto_results='/Volumes/omics4tb2/alomana/projects/cdi/results/kallisto.640/'
expression_dir='/Volumes/omics4tb2/alomana/projects/cdi/results/expression.640/'

# 1. generating full expression matrix
print('generating expression matrix file...')

elements=os.listdir(kallisto_results)
kallisto_dirs=[element for element in elements if os.path.isdir(kallisto_results+element) == True]
kallisto_dirs.sort()

# 1.1. reading the genes
transcript_names=[]
one_file=kallisto_results+kallisto_dirs[0]+'/abundance.tsv'

with open(one_file,'r') as f:
    next(f)
    for line in f:
        vector=line.split('\t')
        name=vector[0]
        transcript_names.append(name)
transcript_names.sort()

# 1.2. reading expression
expression={}
for replicate in kallisto_dirs:
    expression[replicate]={}

    working_file=kallisto_results+replicate+'/abundance.tsv'
    with open(working_file,'r') as f:
        next(f)
        for line in f:
            vector=line.split('\t')
            name=vector[0]
            abundance=vector[-1].replace('\n','')
            expression[replicate][name]=abundance

# 1.3. write expression matrix
matrix_file=expression_dir+'expression_matrix.txt'
with open(matrix_file,'w') as f:
    
    for replicate in kallisto_dirs:
        f.write('\t{}'.format(replicate))
    f.write('\n')

    for name in transcript_names:
        f.write(name)
        for replicate in kallisto_dirs:
            f.write('\t{}'.format(expression[replicate][name]))
        f.write('\n')
