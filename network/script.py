import os, sys

# 0.1. folders
data_dir = 'data/'
outdir_dir = 'results.MMRF_CoMMpass_IA14a_E74GTF_Salmon_Gene_TPM/'

if os.path.exists(outdir_dir) == False:
    os.mkdir(outdir_dir)

# 0.2. files
expression_file = '{}MMRF_CoMMpass_IA14a_E74GTF_Salmon_Gene_TPM.txt'.format(data_dir)
synonyms_file = '{}identifier_mappings.txt'.format(data_dir)
coexpression_results_json = '{}coexpression_results/coexpressionDictionary.json'.format(outdir_dir)
regulons_json = '{}mechanistic_inference_results/regulons.json'.format(outdir_dir)

# 1. coexpression
cmd = 'miner2-coexpr {} {} {}coexpression_results'.format(expression_file, synonyms_file, outdir_dir)

print()
print(cmd)
print()
os.system(cmd)

# 2. mechanistic inference
cmd = 'miner2-mechinf {} {} {} {}mechanistic_inference_results'.format(expression_file, synonyms_file, coexpression_results_json, outdir_dir)

print()
print(cmd)
print()
os.system(cmd)

# 3. network activity
cmd = 'miner2-bcmembers {} {} {} {}bicluster_membership_results'.format(expression_file, synonyms_file, regulons_json, outdir_dir)

print()
print(cmd)
print()
os.system(cmd)

# 4. subtypes
cmd = 'miner2-subtypes {} {} {} {}subtype_analysis_results'.format(expression_file, synonyms_file, regulons_json, outdir_dir)

print()
print(cmd)
print()
os.system(cmd)
