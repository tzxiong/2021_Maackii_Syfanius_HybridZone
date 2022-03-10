# Rules on updating the variables
# use continuous-time approximation
function update_variables_equilibriumIM(variables,t,parameters)
    
    X1,X2,C1,C2,n = variables
    N1,N2,μ1,μ2,K1,K2,m,L,λ = parameters
    # K1: the number of demes in species 1
    # K2: (the number of demes in species 2) - 1
    # species 1: 1,2,3,4,...,K2
    # species 2: 0,-1,-2,-3,-4,...,-K1
    # N1: population size in species 1 / deme
    # N2: population size in species 2 / deme
    # μ1, μ2: substitution rates
    # m: migration rate between deme 1 and deme 0
    # λ: the ratio between migration rates within species and migration rates between species
    # the migration probability between conspecific demes are always 0.5
    
    # accumulate mutations
    pC1 = (μ2 * (X1>0) + μ1 * (X1<1)) * L
    pC2 = (μ2 * (X2>0) + μ1 * (X2<1)) * L
    
    # movement of lineages within species, +1(forward)
    pX1_w_f = ((X1>1) + (X1<0)) * λ * m
    pX2_w_f = ((X2>1) + (X2<0)) * λ * m
    
    # movement of lineages within species, -1(backward)
    pX1_w_b = ((X1>1) + (X1<0)) * λ * m
    pX2_w_b = ((X2>1) + (X2<0)) * λ * m
    
    # movement of lineages at the contact zone across species
    pX1_c_b = ((X1==1) + (X1==0)) * m
    pX2_c_b = ((X2==1) + (X2==0)) * m
    
    # movement of lineages at the contact zone back to its own species
    pX1_c_w = ((X1==1) & (K2>1) + (X1==0) & (K1<0)) * m * λ
    pX2_c_w = ((X2==1) & (K2>1) + (X2==0) & (K1<0)) * m * λ
    
    # determining if they coalesce
    
    pn = (X1==X2) ? 1/((X1>0) * N2 + (X1<=0) * N1) : 0
    
    rate_vec = [pC1,pC2,pX1_w_f,pX2_w_f,pX1_w_b,pX2_w_b,pX1_c_b,pX2_c_b,pX1_c_w,pX2_c_w,pn]
    
    pTotal = sum(rate_vec)
    
#     append!(t,t[end]+rand(Exponential(pTotal)))
    
    i = rand(Categorical(rate_vec./pTotal))
    
    if i===1
        variables[3] = variables[3] + 1
    elseif i===2
        variables[4] = variables[4] + 1
    elseif i===3
        variables[1] = (variables[1] < K2) ? (variables[1] + 1) : variables[1]
    elseif i===4
        variables[2] = (variables[2] < K2) ? (variables[2] + 1) : variables[2]
    elseif i===5
        variables[1] = (variables[1] > -K1) ? (variables[1] - 1) : variables[1]
    elseif i===6
        variables[2] = (variables[2] > -K1) ? (variables[2] - 1) : variables[2]
    elseif i===7
        variables[1] = 1 - variables[1]
    elseif i===8
        variables[2] = 1 - variables[2]
    elseif i===9
        variables[1] = variables[1] - (-1)^variables[1]
    elseif i===10
        variables[2] = variables[2] - (-1)^variables[2]
    elseif i===11
        variables[5] = 1
    end
    
    return variables
    
end

function calculate_ensemble_rate_ratio_equilibriumIM(variables_init,parameters;ensemble_run = 1e3)
    
    X1,X2,C1,C2,n = variables_init
    
    C1_data = []
    C2_data = []

    @showprogress 1 for i in 1:ensemble_run
    
        variables = [X1,X2,C1,C2,n]
        
        # println(variables)
        
        t = [0.0]
    
        while variables[end] === 0
            variables = update_variables_equilibriumIM(variables,t,parameters)
            # println(variables)
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



function ensemble_simulation_equilibrium_IM(N1,N2,μ1,μ2,K1,K2,m_data,L,λ;ensemble_run = 1e4)
    
    # default parameters
    # N1 = 1/6*1e4
    # N2 = 1/6*1e4
    # μ1 = 3e-9
    # μ2 = 3*3e-9
    # K1 = 10
    # K2 = 10
    # m_data = (0.1).^Array(1.4:0.25:5.4)
    # L = 1e4
    # λ = 2
    # ensemble_run = 1e4

    # initialize the variables
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

        println(m)
        parameters = [N1,N2,μ1,μ2,K1,K2,m,L,λ]
        slope1,C1_muts_1,C2_muts_1 = calculate_ensemble_rate_ratio_equilibriumIM([-K1,K2,C1,C2,n],parameters;ensemble_run = ensemble_run)
        slope2,C1_muts_2,C2_muts_2 = calculate_ensemble_rate_ratio_equilibriumIM([K2,K2,C1,C2,n],parameters;ensemble_run = ensemble_run)
        slope3,C1_muts_3,C2_muts_3 = calculate_ensemble_rate_ratio_equilibriumIM([-K1,-K1,C1,C2,n],parameters;ensemble_run = ensemble_run)
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

N_scales = [1]

μ_scales = [3]

N1 = 1/6*1e4
N2 = 1/6*1e4
μ1 = 3e-9
μ2 = 3e-9
K1 = 9 
K2 = 10
# m_data = (0.1).^Array(2:0.5:5)
m_data = (0.1).^Array(1.4:0.25:5.4)
L = 1e4
λ = 1

ensemble_run = 1e4

for i in N_scales
    
    for j in μ_scales
        
        println([i,j])

        slope_data,D01_data,D00_data,D11_data = ensemble_simulation_equilibrium_IM(N1*i,
                                                N2,
                                                μ1,
                                                j*μ2,
                                                K1,
                                                K2,
                                                m_data,
                                                L,
                                                λ,
                                                ensemble_run = ensemble_run)
        
        slope_data_results["$(i),$(j)"] = slope_data
        D01_data_results["$(i),$(j)"] = D01_data
        D00_data_results["$(i),$(j)"] = D00_data
        D11_data_results["$(i),$(j)"] = D11_data
        
    end
    
end
