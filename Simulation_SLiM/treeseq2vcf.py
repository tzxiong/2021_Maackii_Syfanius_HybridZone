import sys,os
import numpy as np
import argparse
import pyslim, msprime, tskit

parser = argparse.ArgumentParser()
parser.add_argument('--dir', type=str, required=True)
parser.add_argument('--seed', type=str, required=True)
parser.add_argument('--gen', type=str, required=True)
parser.add_argument('--sampleSize', type=str, required=True)
parser.add_argument('--trialSize', type=str, required=True)

args = parser.parse_args()

# inputs
dir_source = args.dir
seed = args.seed
gen = args.gen
sampleSize = int(args.sampleSize)
trialSize=int(args.trialSize)

print("seed:", flush=True)
print(seed, flush=True)
print("gen:", flush=True)
print(gen, flush=True)

dir_source = dir_source +"/slim_seed"+seed

os.chdir(dir_source)

ts = pyslim.load("slim_seed"+seed+"_gen"+gen+".trees")

ts_mut = pyslim.SlimTreeSequence(msprime.sim_mutations(ts, rate=3e-7, keep=False))

alive = ts_mut.individuals_alive_at(0)

population_id = []
for i in alive:
    population_id.append(ts_mut.individual(i).metadata['subpopulation'])
population_id = np.array(population_id)

sp_size=sampleSize ## sample size

for trial in range(trialSize):

    print('trial:')
    print(trial+1)

    var_missing = True

    while var_missing:

        sp1_extreme_idx = np.random.choice(np.where(population_id==9)[0],sp_size, replace=False)
        sp2_extreme_idx = np.random.choice(np.where(population_id==10)[0],sp_size, replace=False)
        sp1_contact_idx = np.random.choice(np.where(population_id==1)[0],sp_size, replace=False)

        sp1_extreme = alive[sp1_extreme_idx].tolist()
        sp2_extreme = alive[sp2_extreme_idx].tolist()
        sp1_contact = alive[sp1_contact_idx].tolist()

        sp1_extreme_nodes = []
        for s in sp1_extreme:
            sp1_extreme_nodes.extend(ts_mut.individual(s).nodes)

        sp2_extreme_nodes = []
        for s in sp2_extreme:
            sp2_extreme_nodes.extend(ts_mut.individual(s).nodes)

        sp1_contact_nodes = []
        for s in sp1_contact:
            sp1_contact_nodes.extend(ts_mut.individual(s).nodes)

        ts_mut_subsampled = ts_mut.simplify(sp1_extreme_nodes
                                            + sp2_extreme_nodes
                                            + sp1_contact_nodes,filter_populations=False)

        test = False
        for v in ts_mut_subsampled.variants():
            test = (test | v.has_missing_data)

        if test:
            var_missing = True
        else:
            var_missing = False

        print(var_missing)

    output_population_order = [9,10,1]

    alive_sampled = ts_mut_subsampled.individuals_alive_at(0)
    population_sampled_id = []
    for i in alive_sampled:
        population_sampled_id.append(ts_mut_subsampled.individual(i).metadata['subpopulation'])

    output_individual_order = []
    for p in output_population_order:
        output_individual_order.extend(alive_sampled[np.where(np.array(population_sampled_id)==p)[0].tolist()])

    print(output_individual_order)  

    ## write samples to vcf
    ## first two groups of sp_size individuals are the extreme populations, and the last batch of sp_size
    ## individuals are from the contact zone
    names = []
    for i in range(len(output_individual_order)):
        names.append("SP"+str(ts_mut_subsampled.individual(output_individual_order[i]).metadata['subpopulation'])
                         +"ID"+str(i+1))
    print(names)
    prefix = "slim_seed"+seed+"_gen"+gen+".sampleSize."+str(sp_size)+".trial."+str(trial+1)

    with open(prefix+".vcf", "w") as vcf_file:
        ts_mut_subsampled.write_vcf(vcf_file,
                         individuals = output_individual_order,
                         individual_names = names)
