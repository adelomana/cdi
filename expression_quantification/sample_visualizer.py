import pandas

# 0. user-defined variables
expression_file = '/Volumes/omics4tb2/alomana/projects/cdi/results/expression.640/expression_matrix.txt'

# 1. read file
df = pandas.read_csv(expression_file,sep='\t',index_col=0)
print(df.head())
print(df.shape)

# 2. subset expression to those that have a min of > 10 TPMs
lower_values = df.min(axis=1)
sub = lower_values[lower_values > 10]

high = df[list(sub.index)]

print(high.head())
print(high.shape)
