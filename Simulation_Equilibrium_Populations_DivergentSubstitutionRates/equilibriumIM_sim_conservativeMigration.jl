println("Launching Julia..");flush(stdout)

using LinearAlgebra, StatsBase, GLM, Statistics, Distributions, HypothesisTests
using ProgressMeter
using JLD, DataFrames,Query, CSV
using Interpolations,Combinatorics



# Rules on updating the variables
# use continuous-time approximation

function update_variables_equilibriumIM(variables,t,parameters)
    
    X1,X2,C1,C2,n = variables
    N1,N2,μ1,μ2,m,L = parameters
    
    # accumulate mutations
    pC1 = (μ1 * (1-X1) + μ2 * (X1)) * L
    pC2 = (μ1 * (1-X2) + μ2 * (X2)) * L
#     C1 = C1 + rand(Bernoulli(pC1))
#     C2 = C2 + rand(Bernoulli(pC2))
    
    # movement of lineages
    pX1 = m * (N2 * (1-X1) + N1 * (X1))/(N1 * (1-X1) + N2 * (X1))
    pX2 = m * (N2 * (1-X2) + N1 * (X2))/(N1 * (1-X2) + N2 * (X2))
#     X1 = pX1 * (1-X1) + (1-pX1) * X1
#     X2 = pX2 * (1-X2) + (1-pX2) * X2
    
    # determining if they coalesce
    
    pn = (X1==X2) ? 1/((1-X1) * N1 + X1 * N2) : 0
    
    pTotal = pC1 + pC2 + pX1 + pX2 + pn
    
#     append!(t,t[end]+rand(Exponential(pTotal)))
    
    i = rand(Categorical([pX1,pX2,pC1,pC2,pn]./pTotal))
    
    variables[i] = (variables[i] + 1) * (i===3) + (variables[i] + 1) * (i===4) + (1 - variables[i]) * (i===1) + (1 - variables[i]) * (i===2) + (1) * (i===5)
    
    return variables
    
end

function calculate_ensemble_rate_ratio_equilibriumIM(variables_init,parameters;ensemble_run = 1e3)
    
    X1,X2,C1,C2,n = variables_init
    
    C1_data = []
    C2_data = []

    for i in 1:ensemble_run
    
        variables = [X1,X2,C1,C2,n]
        
        t = [0.0]
    
        while variables[end] === 0
            variables = update_variables_equilibriumIM(variables,t,parameters)
        end
    
        append!(C1_data,variables[3])
        append!(C2_data,variables[4])
    
    end
    
#     ols = lm(@formula(Y ~ X+0),DataFrame(X=float.(C2_data), Y=float.(C1_data)))
#     slope = coef(ols)[1] # slope: number of mutations in lineages 2 divided by that in lineage 1
    
    slope = sum(C2_data)/sum(C1_data)
    C1_muts = mean(C1_data)
    C2_muts = mean(C2_data)
    
    return slope,C1_muts,C2_muts
     
end


# theoretical curve
r(fst,rmax)=(1+rmax+fst*(rmax-1))/(1+rmax-fst*(rmax-1))

function ensemble_simulation_equilibrium_IM(N1,N2,μ1,μ2,m_data,L;ensemble_run = 1e4)
    
    # default parameters
    # N1 = 1/6*1e4
    # N2 = 1/6*1e4
    # μ1 = 3e-9
    # μ2 = 3*3e-9
    # m_data = (0.1).^Array(1.4:0.25:5.4)
    # L = 1e4
    # ensemble_run = 1e4

    # initialize the variables
    # X1,X2 ∈ {1,-1} are variables determining the location of the two lineages
    X1 = 0
    X2 = 1
    # C1,C2 ∈ nonnegative intergers are counts of neutral substitutions of the two lineages
    C1 = 0
    C2 = 0
    # n = 0: not coalesced; n = 1: coalesced
    n = 0
    # t: time
    t = [0.0]

    slope_data = Float64[]
    D01_data = Float64[]
    D00_data = Float64[]
    D11_data = Float64[]

    for m in m_data

        println(m);flush(stdout)
        variables_init = [X1,X2,C1,C2,n]
        parameters = [N1,N2,μ1,μ2,m,L]
        slope1,C1_muts_1,C2_muts_1 = calculate_ensemble_rate_ratio_equilibriumIM([0,1,C1,C2,n],parameters;ensemble_run = ensemble_run)
        slope2,C1_muts_2,C2_muts_2 = calculate_ensemble_rate_ratio_equilibriumIM([1,1,C1,C2,n],parameters;ensemble_run = ensemble_run)
        slope3,C1_muts_3,C2_muts_3 = calculate_ensemble_rate_ratio_equilibriumIM([0,0,C1,C2,n],parameters;ensemble_run = ensemble_run)
        append!(slope_data,slope1)
        append!(D01_data,C1_muts_1 + C2_muts_1)
        append!(D11_data,C1_muts_2 + C2_muts_2)
        append!(D00_data,C1_muts_3 + C2_muts_3)

    end

    return slope_data,D01_data,D00_data,D11_data
    
end


slope_data_results = Dict()
D01_data_results = Dict()
D00_data_results = Dict()
D11_data_results = Dict()

N_scales = vcat(Array(0.1:0.2:0.9),Array(1:2:10))

μ_scales = [0.5,1,1.5,3,5]

N1 = 1/6*1e4
N2 = 1/6*1e4
μ1 = 3e-9
μ2 = 3e-9
L = 1e4
ensemble_run = 1e6

for i in N_scales
    
    for j in μ_scales
        
        println([i,j]);flush(stdout)

        slope_data,D01_data,D00_data,D11_data = ensemble_simulation_equilibrium_IM(N1*i,
                                                N2,
                                                μ1,
                                                j*μ2,
                                                (0.1).^Array(1.4:0.25:5.4),
                                                L,
                                                ensemble_run = ensemble_run)
        
        slope_data_results["$(i),$(j)"] = slope_data
        D01_data_results["$(i),$(j)"] = D01_data
        D00_data_results["$(i),$(j)"] = D00_data
        D11_data_results["$(i),$(j)"] = D11_data
        
    end
    
end

dir_output = "/n/home00/txiong/Research/2020_BarriersToGeneFlow/RateBias_Simulations/EquilibriumIM_VariousPopSize_VariousMutRates"

save(string(dir_output,"/results.conservativeMigration.jld"),
    "slope_data_results", slope_data_results,
    "D01_data_results", D01_data_results,
    "D00_data_results", D00_data_results,
    "D11_data_results", D11_data_results)
