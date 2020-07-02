###
### This script calculates the number of DEGs per regulon in order to have significance, based on the hypergeometric test
###

import scipy, scipy.stats
import matplotlib, matplotlib.pyplot

# We can then compute a probability of drawing X red marbles out of N from a jar containing n red marbles out of M in the following way:
#pval = hypergeom.sf(x-1, M, n, N)

# DEGs in a regulon is success
# x = 5       # number of drawn successes, aka the number of genes in a regulon that are DEGs
# M = 25,000  # population size, aka the number of genes in the genome
# n = 1900    # number of successes in the population, aka, the number of DEGs in the genome
# N = 10      # the sample size, aka the number of genes in the regulon.

###
# 0. user-defined variables
###
results_folder = '/Users/alopez/projects_isb/cdi/results/regulons/thresholds/'

M = 51505-1   # number of genes in the regerence genome
n = 1967      # number of DEGs in genome

###
# 1. analysis
###
membership_ranks = []
threshold_ranks = []

regulon_sizes = list(range(3, 200)) # 200 is above the our limit of regulon size
for regulon_size in regulon_sizes:
    membership_ranks.append(regulon_size)
    successes = list(range(1, regulon_size))
    for success in successes:
        N = regulon_size
        x = success
        pval = scipy.stats.hypergeom.sf(x-1, M, n, N)
        if (pval < 0.05) and (len(membership_ranks) != len(threshold_ranks)):
            print('in', x, N)
            threshold_ranks.append(x)
        print(x, N, pval)
    print()
matplotlib.pyplot.plot(membership_ranks, threshold_ranks, 'o-', color='black')
matplotlib.pyplot.savefig(results_folder + 'figure.pdf')

###
# 3. store the results
###
out_file = results_folder + 'results.txt'
with open(out_file, 'w') as f:
    f.write('# regulon size\tDEGs_pass\n')
    for i in range(len(membership_ranks)):
        f.write('{}\t{}\n'.format(membership_ranks[i], threshold_ranks[i]))
