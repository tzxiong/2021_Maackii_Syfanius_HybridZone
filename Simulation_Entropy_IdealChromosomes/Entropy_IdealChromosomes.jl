using Plots; pyplot()
using Plots.Measures
# using StatsPlots; pyplot()
using PyCall

using LinearAlgebra, StatsBase, GLM, Statistics, Distributions, HypothesisTests
using DelimitedFiles
using SpecialMatrices
using FFTW
using ProgressMeter
using JLD, DataFrames,Query, CSV
using KernelDensity
using Interpolations,Combinatorics,Dierckx,LsqFit
using LaTeXStrings

function simulate_ancestry_unit_interval(λ;N=2,p=[0.5,0.5])
    
    # λ: the rate of break-point appearance
    # p: relative probability of each ancestry
    # N: total number of ploidy

    anc = Dict() # storing each haplotype's ancestry

    for n in 1:N
    
        breakpoints = sort(rand(Uniform(0,1),rand(Poisson(λ))))
        append!(breakpoints,1)
        categories = rand(Categorical(p),size(breakpoints)[1])
        anc["$(n),RightBoundary"] = breakpoints
        anc["$(n),Category"] = categories        
    end
    
    return anc
    
end

function simulate_ancestry_discrete_loci(λ,N_loci;N=2,p=[0.5,0.5])
    
    # λ: the probability of break-point appearance
    # p: relative probability of each ancestry
    # N: total number of ploidy

    anc = Dict() # storing each haplotype's ancestry

    for n in 1:N
        
        B = rand(Binomial(N_loci-1,λ),1)[1]
        breakpoints = sample(1:(N_loci-1),B)
        append!(breakpoints,N_loci)
        categories = rand(Categorical(p),size(breakpoints)[1])
        breakpoints = breakpoints ./ breakpoints[end]
        anc["$(n),RightBoundary"] = breakpoints
        anc["$(n),Category"] = categories        
    end
    
    return anc
    
end

function group_into_unphased_ancestry(anc;N=2,p_dims=2)
    
    boundaries = []
    for n in 1:N
        boundaries = vcat(boundaries,anc["$(n),RightBoundary"])
    end
    boundaries = sort(unique(boundaries)) 
    
    distributions = zeros(size(boundaries)[1],p_dims)
    
    i = 1
    for b in boundaries
        
        P = []
        
        for n in 1:N
            idx = findfirst(b .<= anc["$(n),RightBoundary"])
            append!(P,anc["$(n),Category"][idx])
        end
        
        for j in 1:p_dims
            distributions[i,j] = sum(P .=== j)/N
        end
        
        i = i+1    
    end
    
    return float.(boundaries),float.(distributions)
    
end

function plot_ancestry(boundaries,distributions)
    
    colors=[:blue,:red,:yellow,:brown]
    
    D = zeros(size(distributions))
    
    n_row = size(D)[1]
    for i in 1:n_row
        D[i,:] = cumsum(distributions[i,:])
    end
    
    n_col = size(D)[2]
    
    fig = plot(size=(500,100))
    
    for j in Array(n_col: -1:1)
    
        plot!(fig,vcat([0],boundaries),vcat([D[1,j]],D[:,j]),
            seriestype=:steppre,
            fillrange = 0,
            fillalpha=1,
#             xlims=(0,1),
#             ylims=(0,1),
            framestyle=:none,
            linecolor=colors[n_col-j+1],
            linewidth=0,
            fillcolor=colors[n_col-j+1],
            left_margin=0mm,
            label="")
    end
    
    return fig
end

function discretize_ancestry(boundaries,distributions;N_loci=5000)
    
    P = []
    X = (1:N_loci)./N_loci
    
    for x in X
        append!(P,distributions[findfirst(y->x<=y,boundaries)])
    end
    
    return P
    
end

complexify(A) = exp.(im .* pi ./2 .* A)

Lambda_b_phasor(A_complex) = begin
    C = Hermitian((transpose(A_complex)*conj.(A_complex)))
    λ = eigen(C).values
end

S_w_phasor(A_complex) = begin
    S_w_array = []
    for i in 1:size(A_complex)[2]
        z = A_complex[:,i]
        s = (abs.(fft(z))).^2
        S = [s[1]]
        s = s[2:end]
        while size(s)[1] >=2 
            append!(S,s[1]+s[end])
            s = s[2:(end-1)]
        end
        if size(s)[1] >0
            append!(S,s[1])
        end
        S = S ./ (sum(S))
        S = S[S .> 0]
        append!(S_w_array,-sum(S .* log.(S)))
    end
    
    mean(S_w_array)
end

S(λ) = begin
    l = λ[λ .> 0]
    l = l./(sum(l))
    sum(-l .* log.(l))
end
  
N_loci = 1000

N = 2
p = [2/4,2/4]

n_samples = 4

# Λ = 0:0.00001:0.001;S_dim = size(Λ)[1]
P = 0:0.01:1;S_dim = size(P)[1]

lambda = 0.001

ensemble_run = 1000

S_b_phasor_ensemble = reshape([],(S_dim,0))
# S_b_vector_ensemble = reshape([],(S_dim,0))
S_w_ensemble = reshape([],(S_dim,0))

@showprogress 1 for ens in 1:ensemble_run
    
S_b_phasor = []
# S_b_vector = []
S_w = []

# for lambda in Λ
for p0 in P

    A = reshape([],N_loci,0)
    
    p = [p0,1-p0]
    for i in 1:n_samples

        anc = simulate_ancestry_discrete_loci(lambda,N_loci;N=N,p=p)
        boundaries,distributions = group_into_unphased_ancestry(anc;N=N,p_dims=size(p)[1])
        a = discretize_ancestry(boundaries,distributions;N_loci=N_loci)
        A = hcat(A,a)

    end
    
    append!(S_w,S_w_phasor(complexify(A)))
    
    λ = reverse(Lambda_b_phasor(complexify(A)))
    λ = λ ./ (sum(λ))
    append!(S_b_phasor,S(λ))
    
end
        
S_b_phasor_ensemble = hcat(S_b_phasor_ensemble,S_b_phasor)
S_w_ensemble = hcat(S_w_ensemble,S_w)
    
end
