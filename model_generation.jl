using DifferentialEquations
using Plots, StatsPlots
using Distributions, Random, FFTW
using Turing, ReverseDiff, Memoization
using LinearAlgebra:Diagonal
using Base.Threads
using FileIO, JLD2

#defining the models
Turing.setadbackend(:reversediff)
Turing.setrdcache(true)

# Ornstein-Uhlenbeck process (No noise)
@model ou(rn,T,delta_t) = begin
    ampl ~ Uniform(0.0,5.0)
    tau ~ Uniform(0.0,5.0)
    
    b = exp(-delta_t/tau)
    
    rn[1] ~ Normal(1,sqrt(ampl))
    
    for i=2:T
        rn[i] ~ Normal(rn[i-1]*b,sqrt(ampl*(1-b^2)))
    end
end

# Ornstein-Uhlenbeck process with added THERMAL Gaussian noise
@model oupn_thermal(rn,T,delta_t,::Type{R}=Vector{Float64}) where {R} = begin
    ampl ~ Uniform(0.0,2.0)
    tau ~ Uniform(0.1,2.0)
    noise_ampl ~ Uniform(0.0,0.5)
    
    b = exp(-delta_t/tau)
    r = R(undef, T)
    
    r[1] ~ Normal(0,sqrt(ampl))
    
    for i=2:T
        r[i] ~ Normal(r[i-1]*b,sqrt(ampl*(1-b^2)))
    end
    rn ~ MvNormal(r,noise_ampl)
end

# Ornstein-Uhlenbeck process with added MULTIPLICATIVE Gaussian noise (Mixed Model)
@model oupn_m(rn,T,delta_t,::Type{R}=Vector{Float64}) where {R} = begin
    ampl ~ Uniform(0.0,2.0)
    tau ~ Uniform(0.1,2.0)
    noise_ampl_m ~ Uniform(0.0,0.5)

    
    b = exp(-delta_t/tau)
    r = R(undef, T)
    
    r[1] ~ Normal(0,sqrt(ampl))
    
    for i=2:T
        r[i] ~ Normal(r[i-1]*b,sqrt(ampl*(1-b^2)))
    end
    rn ~ MvNormal(r,sqrt.(noise_ampl_m^2*abs.(r)))
end

# Ornstein-Uhlenbeck process with added MULTIPLICATIVE AND THERMAL Gaussian noise (Mixed Model)
@model oupn_mt(rn,T,delta_t,ratio,::Type{R}=Vector{Float64}) where {R} = begin
    ampl ~ Uniform(0.0,2.0)
    tau ~ Uniform(0.1,2.0)
    noise_ampl_t ~ Uniform(0.0,0.5)
    noise_ampl_m = noise_ampl_t*ratio

    
    b = exp(-delta_t/tau)
    r = R(undef, T)
    
    r[1] ~ Normal(0,sqrt(ampl))
    
    for i=2:T
        r[i] ~ Normal(r[i-1]*b,sqrt(ampl*(1-b^2)))
    end
    rn ~ MvNormal(r,sqrt.(noise_ampl_t^2 .+ noise_ampl_m^2*abs.(r)))
end

#note use double quotes for keyword argument (important in julia)
function data_generation(dt, Tnoise,ratio=1; model="t")
    μ = 0.0
    σ = sqrt(2)
    Θ = 1.0
    W = OrnsteinUhlenbeckProcess(Θ,μ,σ,0.0,1.0)
    prob = NoiseProblem(W,(0.0,7000.0))
    sol_pre = [solve(prob;dt=i).u for i in dt]
    sol = [i[1:1500] for i in sol_pre]

    noise_t = [rand.(Normal.(0,Tnoise),length(sol[i])) for i in 1:length(dt)]

    Mnoise=Tnoise

    noise_m = copy(sol)
    for i in eachindex(dt)
        #might be a bug ratio called twice on Mnoise
        noise_m[i] = [rand.(Normal.(0,ratio*Mnoise*sqrt.(abs.(j)))) for j in (sol[i])]
    end

    #combine base of keyword argument

    if model == "t"
        noisedata = noise_t
    elseif model == "m"
        noisedata = noise_m
    elseif model == "mt"
        noisedata = noise_t.+noise_m
    else
        throw(DomainError(model))
    end
    
    return noisedata .+ sol

end

# return multiplicative, thermal, and signal noise seperatly
function data_generation_sep(dt, Tnoise,ratio=1)
    μ = 0.0
    σ = sqrt(2)
    Θ = 1.0
    W = OrnsteinUhlenbeckProcess(Θ,μ,σ,0.0,1.0)
    prob = NoiseProblem(W,(0.0,1500.0))
    sol_pre = [solve(prob;dt=i).u for i in dt]
    sol = [i[1:1500] for i in sol_pre]

    noise_t = [rand.(Normal.(0,Tnoise),length(sol[i])) for i in 1:length(dt)]

    Mnoise=ratio*Tnoise

    noise_m = copy(sol)
    for i in eachindex(dt)
        noise_m[i] = [rand.(Normal.(0,Mnoise*sqrt.(abs.(j)))) for j in sol[i]]
    end

    return [noise_t, noise_m, sol]
end