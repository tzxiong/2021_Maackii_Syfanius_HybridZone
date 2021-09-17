using PyCall
np = pyimport("numpy")
allel = pyimport("allel")

using LinearAlgebra, StatsBase, GLM, Statistics, Distributions, HypothesisTests
using DelimitedFiles
using SpecialMatrices
using FFTW
using ProgressMeter
using JLD, DataFrames,Query, CSV
using KernelDensity
using Interpolations,Combinatorics
using LaTeXStrings
using ArgParse


s_arg = ArgParseSettings()

@add_arg_table s_arg begin
    "--dir"
        arg_type = String
        required = true
    "--gen"
        arg_type = Int
        required = true
    "--elaiGen"
        arg_type = Int
        required = true
end

parsed_args = parse_args(ARGS, s_arg)

source_dir1 = parsed_args["dir"] # the folder where you can find "slim_seedXXX..." folders for each SLiM seed's simulations
gen = parsed_args["gen"] # the generation number of the ".tree" file from SLiM
elaiGen = parsed_args["elaiGen"] # the ELAI generation parameter


## read seeds from a seeds file
seeds = DataFrame(CSV.File(string(source_dir1,"/seeds.txt"); header=false))[:,1];


function get_fst_per_chr(dir,seed,gen,sampleSize,trial,chr,chr_L,n_segments,pop1,pop2;biallelic=true)

    vcf_data = allel.read_vcf(string(dir,"/slim_seed$(seed)/slim_seed$(seed)_gen$(gen).sampleSize.$(sampleSize).trial.$(trial).chr.$(chr).vcf"))
    gt_array = allel.GenotypeArray(vcf_data["calldata/GT"])
    pos = vcf_data["variants/POS"]

    chr_left = Int64(floor(pos[1]/chr_L))*chr_L + 1
    chr_right = chr_left + chr_L - 1

    segment_starts = Array(chr_left : chr_L / n_segments : chr_right)[1:n_segments]
    segment_ends = vcat(segment_starts[2:end].-1,[chr_right])

    ac_both = gt_array.subset(sel1=vcat(pop1,pop2)).count_alleles()

    if biallelic
        biallelic_pointer = ac_both.is_biallelic()
        ac_1 = gt_array.subset(sel0=biallelic_pointer,sel1=pop1).count_alleles()
        ac_2 = gt_array.subset(sel0=biallelic_pointer,sel1=pop2).count_alleles()
        pos = pos[biallelic_pointer]
    else
        ac_1 = gt_array.subset(sel1=pop1).count_alleles()
        ac_2 = gt_array.subset(sel1=pop2).count_alleles()
        pos = pos[:]
    end

    Fst = []

    pos_index = Array(1:size(pos)[1]) .- 1

    for (l1,l2) in zip(segment_starts,segment_ends)

        new_index = pos_index[(pos .>= l1) .& (pos .<= l2)]

        ac1_segment = ac_1.take(new_index)
        ac2_segment = ac_2.take(new_index)

        push!(Fst,allel.moving_hudson_fst(ac1_segment, ac2_segment,size(new_index)[1])[1])
    end

    return float.(Fst)
end

test_statistic(y,x) = cor(y,x) #* std(x) / sqrt(1/4)

function block_jacknife(y,x,blocksize)

    N_data = size(y)[1]
    repeats = Int64(N_data/blocksize)
    statistics = []
    for i in 1:repeats
        idx = vcat(Array(1:((i-1)*blocksize)),Array((i*blocksize+1):N_data))
        append!(statistics,test_statistic(y[idx],x[idx]))
    end

    return mean(statistics), sqrt(sum((statistics .- mean(statistics)).^2)*(size(statistics)[1]-1)/size(statistics)[1]), repeats

end

function get_entropy_fst_relationship(dir,seed,gen,elaiGen,sampleSize,trial,chr_L,n_segments,pop1,pop2;Entropy="Shannon")

    if !isfile(string(dir,"/slim_seed$(seed)","/slim_seed$(seed)_gen$(gen).sampleSize.$(sampleSize).trial.$(trial).entropy.elaiGen.$(elaiGen).n_segments.$(n_segments).FourierPrecision.1e-4.jld"))

        return NaN,NaN,NaN,NaN

    end

    # get fst
    Fst_all = Dict()

    for chr in 1:10
        fst_multi = get_fst_per_chr(dir,seed,gen,sampleSize,trial,chr,chr_L,n_segments,pop1,pop2;biallelic=false)
        Fst_all["$(chr),$(n_segments)"] = fst_multi
    end


    cd(string(dir,"/slim_seed$(seed)"))

    S = JLD.load(string(dir,"/slim_seed$(seed)","/slim_seed$(seed)_gen$(gen).sampleSize.$(sampleSize).trial.$(trial).entropy.elaiGen.$(elaiGen).n_segments.$(n_segments).FourierPrecision.1e-4.jld"))

    S = S["S_by_chr"]

    if Entropy=="Shannon"
        j = 1
    elseif Entropy =="Tsallis-2"
        j = 3
    elseif Entropy =="Renyi-2"
        j = 5
    end

    data_S_v = reshape([],0,2)
    data_S_f = reshape([],0,2)

    statistics_S_v = Dict()
    statistics_S_f = Dict()

    for chr in 1:10
        S_v = S["$chr"][j]
        S_f = S["$chr"][j+1]
        Fst = Fst_all["$(chr),$(n_segments)"]

        for s in 1:n_segments
            data_S_v = vcat(data_S_v,hcat(repeat([Fst[s]],50),S_v[:,s]))
            data_S_f = vcat(data_S_f,hcat(repeat([Fst[s]],50),S_f[:,s]))
            statistics_S_v["$(chr),$(s),mean"] = mean(S_v[:,s])
            statistics_S_v["$(chr),$(s),std"] = std(S_v[:,s])
            statistics_S_f["$(chr),$(s),mean"] = mean(S_f[:,s])
            statistics_S_f["$(chr),$(s),std"] = std(S_f[:,s])
        end

    end

    data_S_v = float.(data_S_v[.! isnan.(data_S_v[:,2].+data_S_v[:,1]),:])
    data_S_f = float.(data_S_f[.! isnan.(data_S_f[:,2].+data_S_f[:,1]),:])


    r = test_statistic(data_S_v[:,2],data_S_v[:,1])
    r_jk,ste_jk,n_blocks = block_jacknife(data_S_v[:,2],data_S_v[:,1],50)

    r_v = r
    r_bias_corr_v = n_blocks * r - (n_blocks - 1) * r_jk
    z_v = r < 0 ? abs(r)/ste_jk : 0

    r = test_statistic(data_S_f[:,2],data_S_f[:,1])
    r_jk,ste_jk,n_blocks = block_jacknife(data_S_f[:,2],data_S_f[:,1],50)

    r_f = r
    r_bias_corr_f = n_blocks * r - (n_blocks - 1) * r_jk
    z_f = r < 0 ? abs(r)/ste_jk : 0

    return r_v,r_f,z_v,z_f

end


## calculate the correlation between fst and entropy
dir = source_dir1



gen = gen
elaiGen_range = [elaiGen]
sampleSize = 4
trial = 1
chr_L = 1e6
n_segments_range = [1,10,20,50,100]
pop1 = [0,1,2,3]
pop2 = [4,5,6,7]

R_v_all = Dict()
R_f_all = Dict()
Z_v_all = Dict()
Z_f_all = Dict()

seed_range=Array(1:50)


for seed in seeds[seed_range]

    println(seed);flush(stdout)

for elaiGen in elaiGen_range

    for n_segments in n_segments_range

#         println([elaiGen,n_segments])

        R_v=[]
        R_f=[]
        Z_v=[]
        Z_f=[]

        i = 1

            for trial in [1]
                r_v,r_f,z_v,z_f= get_entropy_fst_relationship(dir,seed,gen,elaiGen,sampleSize,trial,chr_L,n_segments,pop1,pop2)
                push!(R_v,r_v)
                push!(R_f,r_f)
                push!(Z_v,z_v)
                push!(Z_f,z_f)
            end

        R_v_all["$(seed),$(elaiGen),$(n_segments)"]=R_v
        R_f_all["$(seed),$(elaiGen),$(n_segments)"]=R_f
        Z_v_all["$(seed),$(elaiGen),$(n_segments)"]=Z_v
        Z_f_all["$(seed),$(elaiGen),$(n_segments)"]=Z_f

    end

end
    
end


## save the pearson's correlation results

dir_output = "/OUTPUT_DIR"

for elaiGen in elaiGen_range

JLD.save(string(dir_output,"/S.vs.Fst.gen.$(gen).elaiGen.$(elaiGen).test_results.jld"),
    "R_v_all",
    R_v_all,
    "R_f_all",
    R_f_all,
    "Z_v_all",
    Z_v_all,
    "Z_f_all",
    Z_f_all)

end
