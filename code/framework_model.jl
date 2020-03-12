using Pkg
Pkg.activate(".")

using Distributions
using Plots
using LinearAlgebra
using GaussianRandomFields
using Distances
using DataFrames
using CSV

reps = 1

@progress for rep in 1:reps

env_df = CSV.read("./landscape_data/env_$rep.csv")
disp_mat = CSV.read("./landscape_data/disp_mat_$rep.csv")
disp_mat = convert(Array, disp_mat)
disp_mat = disp_mat[:,2:end]
landscape = CSV.read("./landscape_data/landscape_$rep.csv")
time = CSV.read("./landscape_data/time_$rep.csv")

S = 50
dominants = trunc(Int,round(S*0.3))
M = size(disp_mat)[1]

env_mat2 = zeros(maximum(env_df.time), M)
for i in 1:maximum(env_df.time)
env_mat2[i,:] = env_df.env1[env_df.time .== i]
end

burn_in = time.burin_in[1]
generations = maximum(env_df.time) - burn_in
Tmax = burn_in + generations

#set up landscape
#plot(env_mat2, legend = false)

r = 5 #

z = reshape(repeat(rand(S), outer = [M]), S, M)'

a_Vect = exp.(range(log(1e-5),stop = log(1),length = 16))
a_Vect = a_Vect[a_Vect.<1]
σ_niche_V = exp.(range(log(0.001),stop = log(10),length = 13))
α_V = [0.5, 1.5] #[0, 0.5, 1, 1.5]

@progress for k in 1:(length(α_V)+2)
    #α = rand(Normal(0.07,0.01),S,S)#rand(LogNormal(0.001, 0.1), S,S)
if(k<=length(α_V))
    if(α_V[k] == 0)
        α = zeros(S,S)*1.0
        α_write = string(α_V[k])
    elseif(k <= length(α_V))
    α = rand(Uniform(0,α_V[k]),S,S)
    α_write = string(α_V[k])
end
elseif(k == (length(α_V)+1))
    α = ones(S,S)*1.0
    α_write = string("equal")
else
    α = rand(Uniform(0,1),S,S)
    α_hold = rand(Uniform(0,1),S,S)
    α[1:dominants,:] = rand(Uniform(1,1.5),dominants,S)
    α[LowerTriangular(α).>0] = α_hold[LowerTriangular(α).>0]
    α_write = string("patch_dynamics")
end
    α[diagind(α)] = repeat([1.0], outer = S)
    α = α*0.05

    @progress for i in 1:length(a_Vect)
        a = a_Vect[i]
        pd_disp = ones(M,S)*a
        pd_disp[:,1:dominants] = ones(M,dominants)*a*0.1
        @progress for j in 1:length(σ_niche_V)
            global Model_df

            σ_niche = σ_niche_V[j]

            N = rand(Poisson(0.5), M, S) * 1.0

            seedV = convert(Array{Int64}, collect(burn_in/(10):burn_in/(10):burn_in/2))
            sampV = convert(Array{Int64}, collect((3*burn_in):Tmax/(100):Tmax))
            N_save = N
            λ_save = zeros(M,S)
            env_save = z
            den_save = zeros(M,S)
            env_match_save = zeros(M,S)

            @progress for gen in 1:Tmax
    #global N
    #global N_save
    #global λ_save
    #global env_save
    #global env_match_save
    #global den_save
    if(any(y -> y == gen, seedV))
        N = N + rand(Poisson(0.5), M, S) * 1.0
    end


    x = reshape(repeat(env_mat2[gen,:], outer = [S]), M, S)

    env = exp.(-((x-z)/(2.0*σ_niche)).^2.0)

    #λ_v = r * N .* (1.0 ./ (1.0 .+ N * α)).*env + disp_mat * N .* a - N .* a

    density = N * α

    λ_v = r * N .* (1.0 ./ (1.0 .+ density)).*env

    λ_v[λ_v .< 0.0] .= 0.0

    N = [rand(Poisson(λ)) for λ in λ_v]

if k < (length(α_V)+2)
    emmigrants = [rand(Poisson(λ)) for λ in N * a]
    immigrants = [rand(Poisson(λ)) for λ in disp_mat * emmigrants]
else
    emmigrants = [rand(Poisson(λ)) for λ in N .* pd_disp]
    immigrants = [rand(Poisson(λ)) for λ in disp_mat * emmigrants * 0.5]
    N[rand(Binomial(1, 0.002), M, S).>0] .= 0
end
    N = N .- emmigrants .+ immigrants
    N[N .< 0.0] .= 0.0

    if(any(y -> y == gen, sampV))
        N_save = cat(dims = 3, N_save, N)
        λ_save = cat(dims = 3, λ_save, λ_v)
        env_save = cat(dims = 3, env_save, x)
        env_match_save = cat(dims = 3, env_match_save, env)
        den_save = cat(dims = 3, den_save, density)
    end
    N = N .* 1.0
end

N_save = N_save[:,:,2:end]
λ_save = λ_save[:,:,2:end]
env_save = env_save[:,:,2:end]
den_save = den_save[:,:,2:end]
env_match_save = env_match_save[:,:,2:end]

Model_df_1 = DataFrame(N = N_save[:],
lambda = λ_save[:],
density = den_save[:],
env_match = env_match_save[:],
env = env_save[:],
#x_coord = repeat(x_y_coords[:,1],outer = length(sampV)*S),
#y_coord = repeat(x_y_coords[:,2],outer = length(sampV)*S),
Species = repeat(1:S, inner = M, outer = length(sampV)),
Time = repeat(1:length(sampV), inner = S*M),
Patch = repeat(1:M, outer = length(sampV)*S),
z = repeat(z[1,:], inner = M, outer = length(sampV)),
dispersal = a,
sig_niche = σ_niche,
alpha = string(α_write))

Model_df_1 = Model_df_1[Model_df_1[:N].>0,:]

if a == a_Vect[1] && σ_niche == σ_niche_V[1] && k ==1
    Model_df = Model_df_1
else
    Model_df = [Model_df; Model_df_1]
end
end
end
end


CSV.write("./outputs/outputfile_$rep.csv", Model_df)
end
