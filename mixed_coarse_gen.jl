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

# Ornstein-Uhlenbeck process
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

function sample_model(dt, Tnoise, oupn, ratio)
    #initial genreated data
    μ = 0.0
    σ = sqrt(2)
    Θ = 1.0  #Θ is the inverse of the relaxation time
    #Tnoise = 0.2
    #dt = 0.05:0.05:2.0
    W = OrnsteinUhlenbeckProcess(Θ,μ,σ,0.0,1.0)
    prob = NoiseProblem(W,(0.0,500.0))
    sol_pre = [solve(prob;dt=i).u for i in dt]
    sol = [i[1:1500] for i in sol_pre]

    #generate  thermal noise (uncomment this block and comment out multiplicative block)
    noise_t = [rand.(Normal.(0,Tnoise),length(sol[i])) for i in 1:length(dt)]
    #noisedata50  = sol .+ noise

    #generate multiplicative noise
    for i in eachindex(dt)
       sol[i] = sol[i] + [rand.(Normal.(0,ratio*Tnoise*sqrt.(abs.(j)))) for j in (sol[i])]
    end

    noisedata50 = sol

    #combine the two
    noisedata50 = noisedata50 .+ noise

    @time oupndata = [sample(oupn(noisedata50[i],length(noisedata50[i]),dt[i],ratio), NUTS(0.65), 1000) for i in eachindex(dt)]
    #advi = ADVI(10, 1000)
    #@time oupndata = [vi(oupn(noisedata50[i],length(noisedata50[i]),dt[i]), advi) for i in eachindex(dt)]
    return oupndata

end

#note use double quotes for keyword argument (important in julia)
function data_generation(dt, Tnoise,ratio=1; model="t")
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


sample_rate = 0.05:0.05:1
num_cycles = 20
Tnoise = 0.2
ratio = 1
mean_total = [Threads.Atomic{Float64}(0) for i in sample_rate]
std_total = [Threads.Atomic{Float64}(0) for i in sample_rate]
data_total = Array{Chains}(undef, num_cycles, length(sample_rate))
@threads for j in 1:num_cycles
    noisedata = data_generation(sample_rate,0.2,model="mt")
    data_total[j,:] = @time [sample(oupn_mt(noisedata[i],length(noisedata[i]),sample_rate[i], 1), NUTS(0.65), 1000) for i in eachindex(sample_rate)]
end


FileIO.save("mixed_points1500_rate0.1:0.1:2_mcmc1000_20cycles.jld2","data",data_total)



# new formatting stuff
#=
mean_noise_t = mean.(map(x-> x[:noise_ampl_t], data))
std_noise_t = std.(map(x-> x[:noise_ampl_t], data))

mean_ampl = mean.(map(x-> x[:ampl], data))
mean_tau = mean.(map(x-> x[:tau], data))
std_ampl = std.(map(x-> x[:ampl], data))
std_tau = std.(map(x-> x[:tau], data))

plot(x_axis, data_parse(std_noise_t), label="Thermal noise", xlabel="Sample Rate", ylabel="STD", title="Mixed Gaussian STD")

plot!(x_axis, data_parse(std_tau), label="Tau")

plot!(x_axis, data_parse(std_ampl), label="Amplitude")

savefig("mixed.png")


plot(x_axis, data_parse(mean_noise_t), label="Thermal noise", xlabel="Sample Rate", ylabel="Mean", title="Mixed Gaussian Mean")
plot!(x_axis, data_parse(mean_tau), label="Tau")
plot!(x_axis, data_parse(mean_ampl), label="Amplitude")
true_mean_noise_t = [0.2 for i in x_axis]
plot!(x_axis, true_mean_noise_t, label="True Noise", line=(:dot, 2))
true_tau = [1.0 for i in x_axis]
plot!(x_axis, true_tau, label="True Tau and Amplitude", line=(:dot, 2))
=#




