println("Julia starts!");flush(stdout)


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
    "--seed"
        arg_type = Int
        required = true
    "--gen"
        arg_type = Int
        required = true
    "--sampleSize"
        arg_type = Int
        required = true
    "--trial"
        arg_type = Int
        required = true
    "--nchrs"
        arg_type = Int
        required = true
    "--elaiGen"
        arg_type = Int
        required = true
    "--nsegments"
        arg_type = Int
        required = true
end

parsed_args = parse_args(ARGS, s_arg)

dir_source = parsed_args["dir"] # the folder where you can find "slim_seedXXX..." folders for each SLiM seed's simulations
seed = parsed_args["seed"] # SLiM seed
gen = parsed_args["gen"] # generation of SLiM output
sampleSize = parsed_args["sampleSize"] # sample size used in tskit subsampling
trial = parsed_args["trial"] # trial number in tskit
n_chrs = parsed_args["nchrs"] # number of chromosomes
elaiGen = parsed_args["elaiGen"] # generation used in ELAI for inference
n_segments = parsed_args["nsegments"] # number of segments per chromosome


function get_position(filename)
    
    ## input is the file .snpinfo.txt, with SNP information
    
    file = readdlm(filename)
    pos = file[2:end,2]
    
    return pos
        
end

function get_tracks(filename,pos)
    
    ## input is the file .ps21.txt, with ancestry information
    
    file = readdlm(filename)
    ANC = Dict()
    LEN = Dict()
    
    for i in 1:size(file)[1]
        
        a = file[i,1:2:end] ## ancestry of the first population, could be within [0,1]

        # transition test
        test = a[1:end-1] .!= a[2:end]
        # the SNP position to the left of the transition point
        left_SNP_pos = pos[1:end-1][test]
        # the SNP position to the right of the transition point
        right_SNP_pos = pos[2:end][test]
        # take the midpoint as the estimated transition point (can be fractional)
        transition_pos = (left_SNP_pos + right_SNP_pos)/2
    
        # store the ancestry (left SNP ancestry + final one on the right)
        anc = vcat(a[1:end-1][test],a[end]) # anc is the vector for ancestry
        
        lengths = vcat(transition_pos,pos[end]) - vcat(pos[1],transition_pos) # lengths is the vector for contiguous tracks
            
        ANC[i] = anc
        LEN[i] = lengths
        
    end  
    
    return ANC, LEN
                    
end

## using exp(i \theta) to represent bi-ancestric hybrid populations
complexify(anc) = exp.(im .* asin.(sqrt.(anc)))

## sampling along the ancestry track, and get discrete ancestries
function discretize_ancestry(anc,len,n)
    
    # anc: a vector of ancestries in interval [0,1]
    # len: lengths of each block, should be normalized!
    # n: a positive integer for the number of discrete blocks used for integration
    
    anc = complexify(anc)
    len = cumsum(len)
    X = Array(0:(len[end]/n):(len[end]))[1:n]
    
    discrete_anc = []
    
    for x in X
        index = findfirst(x .< len)
        append!(discrete_anc,anc[index])
    end
    
    return X, discrete_anc
    
end

# return the average of an arbitrary square matrix along each diagonal lines
function diag_average(M)
    # the average takes place from the upperright diag to the lowerleft diag
    
    n = size(M,1)
    m = []
    
    # uppertriangle
    for i in n:-1:1
        R = Array(0:n-i)
        append!(m,sum(collect(M[1+r,i+r] for r in R))/size(R)[1])
    end
    
    # lower triangle
    for j in 2:n
        R = Array(0:n-j)
        append!(m,sum(collect(M[j+r,1+r] for r in R))/size(R)[1])
    end
    
    return m
end

# re-partition an array of data
function re_partition(a,cumsum_l_old,cumsum_l_new)
    
    # a and cumsum_l_old must be of the same size
    # cumsum_l_old[end] == cumsum_l_new[end] must be true
    
    n = size(cumsum_l_new)[1]
    b = []
    
    for i in 1:n
        
        append!(b,a[findfirst(x->x>=cumsum_l_new[i],cumsum_l_old)])
    end
    
    return b
    
end

function to_fourier_spectral_entropy(a,cumsum_l,tolerance,entropy)
    # this is for transforming piece-wise contant functions on [0,L]
    
    cumsum_l = cumsum_l .* (1/cumsum_l[end]) # Fourier coefficients are invariant under time-rescaling
    x_right = cumsum_l
    x_left = vcat([0],cumsum_l[1:end-1])
    dl = x_right .- x_left
    
    total_energy = sum((abs.(a)).^2 .* dl)
    
    averaged = sum(a .* dl)
    
    energy_explained = abs(averaged)^2
    
    energy_spectrum=[energy_explained]
    
    n=1
    
#     prog = ProgressThresh(1e-4, "Minimizing:")
    
    while (total_energy - energy_explained)/total_energy > tolerance
        new_spectrum_mass_1 = sum(a.*(exp.(-im.*2 .*pi.*n.*x_right).- exp.(-im.*2 .*pi.*n.*x_left)))*im/(2*pi*n)
        new_spectrum_mass_2 = sum(a.*(exp.(-im.*2 .*pi.*(-n).*x_right).- exp.(-im.*2 .*pi.*(-n).*x_left)))*im/(2*pi*(-n))
        new_energy = abs(new_spectrum_mass_1)^2 + abs(new_spectrum_mass_2)^2
        energy_explained = energy_explained + new_energy
        
        n = n+1
        
        append!(energy_spectrum,new_energy)
        
#         ProgressMeter.update!(prog, (total_energy - energy_explained)/total_energy)
    end
    
    energy_spectrum = energy_spectrum ./ total_energy
    energy_spectrum_positive = energy_spectrum[energy_spectrum .> 0]
    S_shannon = -sum(energy_spectrum_positive .* log.(energy_spectrum_positive))
    S_tsallis2 =  1 - sum(energy_spectrum .^ 2)
    S_renyi2 =  -log(sum(energy_spectrum .^ 2))
    return S_shannon, S_tsallis2, S_renyi2
end

function to_vonNeumann_spectrum(anc,cumsum_l,n_samples)
    
    # anc and cumsum_l are dictionaries with entries corresponding to each samples
    
    Z = Complex.(zeros(n_samples,n_samples))
    for i in 1:n_samples
        for j in 1:n_samples
            
            l1 = cumsum_l[i]
            l2 = cumsum_l[j]
            l = sort(unique(vcat(l1,l2)))
            dl = l .- vcat([0],l[1:end-1])
                
            a1 = anc[i]
            a2 = anc[j]
                
            a1_new=re_partition(a1,l1,l)
            a2_new=re_partition(a2,l2,l)
                
            Z[i,j] = sum(a1_new .* conj.(a2_new) .* dl)/sum(dl)
            
        end
    end
    
    Z = Hermitian(Z)
    lambda = (eigen(Z).values)./n_samples
    
    return lambda
end

# ensemble analysis on ancestry estimation data from ELAI
function direct_spectrum_ensemble_ELAI(dir,seed,gen,sampleSize,trial,chr,elaiGen;n_segments=1)
    
    # go to the working directory of ELAI outputs, containing .snpinfo.txt and .ps21.txt
    cd(string(dir,"/slim_seed$(seed)/output"))
    
    # seed: the seed used for SLiM simulation
    # gen: the generation of SLiM output
    # sampleSize: the sampleSize of tskit re-sampling to generate sample treeseq
    # trial: the number of trial in tskit re-sampling to generate sample treeseq
    # chr: chromosome (1-10)
    # elaiGen: the generation used in ELAI, integer
    
    seed = string(seed)
    gen = string(gen)
    sampleSize = string(sampleSize)
    trial = string(trial)
    chr = string(chr)
    elaiGen = string(elaiGen)
    
    
     
    prefix = string("slim_seed$(seed)_gen$(gen).sampleSize.$(sampleSize).trial.$(trial).chr.$(chr).bimbam.admx.output.elaiGen.$(elaiGen).run.")
    pos = get_position(string(prefix,"1.snpinfo.txt"))
    
    S_vonNeumann_shannon = zeros(50,n_segments) # storing the von Neumann entropy
    S_vonNeumann_tsallis2 = zeros(50,n_segments)
    S_vonNeumann_renyi2 = zeros(50,n_segments)
    S_Fourier_shannon = zeros(50,n_segments) # storing the sample-averaged Fourier spectral entropy
    S_Fourier_tsallis2 = zeros(50,n_segments)
    S_Fourier_renyi2 = zeros(50,n_segments)
    
    # loading data
    for run in 1:50

        anc,len = get_tracks(string(prefix,string(run),".ps21.txt"),pos)
        
        n_samples = length(anc)
        
        windows = Array(0:1/n_segments:1)
        
        # normalize the lengths to a total length of 1 
        # complexify the signals
        cumsum_l_per_segment = Dict()
        anc_per_segment = Dict()
        for i in 1:n_samples
            cumsum_l_tmp = cumsum(len[i])./(cumsum(len[i])[end])
            len[i] = len[i]./sum(len[i])
            anc[i] = complexify(anc[i])
            for s in 1:n_segments
                w_left = windows[s]
                w_right = windows[s+1]
                idx_exterior = [findfirst(cumsum_l_tmp.>=w_right)]
                idx_interior = (cumsum_l_tmp.<w_right) .* (cumsum_l_tmp.>w_left)
                anc_per_segment["$i,$s"] = vcat(anc[i][idx_interior],anc[i][idx_exterior])
                cumsum_l_tmp_segment = vcat(cumsum_l_tmp[idx_interior],cumsum_l_tmp[idx_exterior])
                cumsum_l_per_segment["$i,$s"] = min.(cumsum_l_tmp_segment .- w_left,w_right-w_left)
            end
        end
        
        for s in 1:n_segments
            
            anc = Dict()
            cumsum_l = Dict()
            
            for i in 1:n_samples
                anc[i] = anc_per_segment["$i,$s"]
                cumsum_l[i] = cumsum_l_per_segment["$i,$s"]
            end
            
            # get von Neumann entropy
            lambda = to_vonNeumann_spectrum(anc,cumsum_l,n_samples)
            
            lambda = abs.(lambda)
            lambda_positive = lambda[lambda .> 0]
            S_shannon = -sum(lambda_positive .* log.(lambda_positive))
            S_tsallis2 =  1 - sum(lambda .^ 2)
            S_renyi2 =  -log(sum(lambda .^ 2))

            S_vonNeumann_shannon[run,s] = S_shannon
            S_vonNeumann_tsallis2[run,s] = S_tsallis2
            S_vonNeumann_renyi2[run,s] = S_renyi2
            
            # calculate the average Fourier spectral entropy
        
            S_shannon = 0
            S_tsallis2 = 0
            S_renyi2 = 0
        
            for i in 1:n_samples
                S_shannon_c, S_tsallis2_c, S_renyi2_c = to_fourier_spectral_entropy(anc[i],cumsum_l[i],1e-4,entropy)
                S_shannon = S_shannon + S_shannon_c
                S_tsallis2 = S_tsallis2 + S_tsallis2_c
                S_renyi2 = S_renyi2 + S_renyi2_c
            
            end
        
            S_shannon = S_shannon/n_samples
            S_renyi2 = S_renyi2/n_samples
            S_tsallis2 = S_tsallis2/n_samples
            
            S_Fourier_shannon[run,s] = S_shannon
            S_Fourier_tsallis2[run,s] = S_tsallis2
            S_Fourier_renyi2[run,s] = S_renyi2
        
        end
        
    end
    
    return [S_vonNeumann_shannon,S_Fourier_shannon,S_vonNeumann_tsallis2,S_Fourier_tsallis2,S_vonNeumann_renyi2,S_Fourier_renyi2]
    
end





## main part



entropy = "all"

S_by_chr = Dict()

Threads.@threads for chr in 1:n_chrs
    
    println("Working on chr $(chr)..");flush(stdout)
    
    S = direct_spectrum_ensemble_ELAI(dir_source,seed,gen,sampleSize,trial,chr,elaiGen;n_segments=n_segments);
    S_by_chr["$chr"] = S
    
    println("Chr $(chr) is done!");flush(stdout)
    
end

dir_output = string(dir_source,"/slim_seed$(seed)")

cd(dir_output)

save("slim_seed$(seed)_gen$(gen).sampleSize.$(sampleSize).trial.$(trial).entropy.elaiGen.$(elaiGen).n_segments.$(n_segments).FourierPrecision.1e-4.jld", "S_by_chr",S_by_chr)




