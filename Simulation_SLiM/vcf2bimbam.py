import sys,os
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=str, required=True)
parser.add_argument('--gen', type=str, required=True)
parser.add_argument('--sampleSize', type=str, required=True)
parser.add_argument('--trial', type=str, required=True)
parser.add_argument('--chrom', type=str, required=True)
parser.add_argument('--dir', type=str, required=True)

args = parser.parse_args()

# inputs
seed = args.seed
gen = args.gen
sampleSize = args.sampleSize
trial  = args.trial
chrom = args.chrom
dir_source = args.dir

print('chromosome '+chrom)

os.chdir(dir_source + '/slim_seed' + seed)

prefix = 'slim_seed'+seed+'_gen'+gen+'.sampleSize.'+sampleSize+'.trial.'+trial+'.chr.'+chrom

vcf_file = open(prefix+'.vcf',"r")

content = vcf_file.readlines()

# trim the header
l = -1
for c in content:
    l = l+1
    if ('#CHROM' in c):
        break

content = content[l:]

# indices for individuals from each populations
pop1 = np.array([0,1,2,3])
pop2 = np.array([4,5,6,7])
admx = np.array([8,9,10,11])

pops = [pop1,pop2,admx]
pops_name = ['pop1','pop2','admx']

splitted = []

i = 0
for line in content:
    line_new = line.strip().split('\t')

    # change 0,1 to ATCG
    if i > 0:
        ref = line_new[3]
        alt = line_new[4].split(',')
        alleles = [ref] + alt
        num_alt = len(alt)
        for i in np.arange(9,len(line_new)):
            line_new[i] = line_new[i].replace('|','')
            for j in range(num_alt+1):
                line_new[i] = line_new[i].replace(str(j),alleles[j])
        if len(set(''.join(line_new[9:]))) == 2:
            splitted.append(line_new)
    else:
        splitted.append(line_new)

    i = i + 1

splitted = np.array(splitted)

sample_names = splitted[0,9:]


# generate bimbam
output_final = np.array([])

n_sites = 0

for i in range(len(pops)):

    pop = pops[i]
    n_samples = np.shape(pop)[0]
    output_sample_names = sample_names[pop]
    output_genotypes = splitted[1:,pop+9]
    right_block = np.concatenate([[output_sample_names],output_genotypes])
    n_sites = np.shape(output_genotypes)[0]

    left_column = ['IND']
    for j in range(n_sites):
        left_column.append('rs'+str(j+1))
    left_column = np.array([left_column])
    output = np.concatenate((np.transpose(left_column),right_block),axis=1)
    output_list = []
    for j in range(np.shape(output)[0]):
        output_list.append(','.join(output[j,:].tolist()))
    output_list = '\n'.join(output_list)
    output = str(n_samples)+'\n'+str(n_sites)+'\n'+output_list

    # write the bimbamoutput
    pop_name = pops_name[i]
    output_filename = prefix + '.bimbam.'+pop_name+'.txt'
    text_file = open(output_filename, "w")
    text_file.write(output)
    text_file.close()

# write the position output
left_column = []
positions = splitted[1:,1]
for n in range(n_sites):
    left_column.append('rs'+str(n+1)+','+positions[n])
left_column = '\n'.join(left_column)
output_filename = prefix + '.bimbam.position.txt'
text_file = open(output_filename, "w")
text_file.write(left_column)
text_file.close()
