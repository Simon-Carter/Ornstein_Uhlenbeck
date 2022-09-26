using DifferentialEquations
using Plots, StatsPlots
using Distributions, Random, FFTW
using Turing, ReverseDiff, Memoization, LinearAlgebra
using Base.Threads

#defining the models
Turing.setadbackend(:reversediff)
Turing.setrdcache(true)

# Ornstein-Uhlenbeck process
@model ou(rn,T,delta_t) = begin
    ampl ~ Uniform(0.0,5.0)
    tau ~ Uniform(0.0,5.0)
    
    b = exp(-delta_t/tau)
    
    rn[1] ~ Normal(0,sqrt(ampl))
    
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

# Ornstein-Uhlenbeck process with added MULTIPLICATIVE Gaussian noise
@model oupn_mult(rn,T,delta_t,::Type{R}=Vector{Float64}) where {R} = begin
    ampl ~ Uniform(0.0,2.0)
    tau ~ Uniform(0.1,2.0)
    noise_ampl ~ Uniform(0.0,0.5)
    
    b = exp(-delta_t/tau)
    r = R(undef, T)
    
    r[1] ~ Normal(0,sqrt(ampl))
    
    for i=2:T
        r[i] ~ Normal(r[i-1]*b,sqrt(ampl*(1-b^2)))
    end
    rn ~ MvNormal(r,noise_ampl*sqrt.(abs.(r)))
end


function sample_model(dt, Tnoise, oupn)
#initial genreated data
μ = 0.0
σ = sqrt(2)
Θ = 1.0  #Θ is the inverse of the relaxation time
#Tnoise = 0.2
#dt = 0.05:0.05:2.0
W = OrnsteinUhlenbeckProcess(Θ,μ,σ,0.0,1.0)
prob = NoiseProblem(W,(0.0,100.0))
sol = [solve(prob;dt=i).u for i in dt]

#generate  thermal noise (uncomment this block and comment out multiplicative block)
#noise = [rand.(Normal.(0,Tnoise),length(sol[i])) for i in 1:length(dt)]
#noisedata50  = sol .+ noise

#generate multiplicative noise

for i in eachindex(dt)
    sol[i] = sol[i] + [rand.(Normal.(0,Tnoise*sqrt.(abs.(j)))) for j in (sol[i])]
end
noisedata50 = sol

@time oupndata = [sample(oupn(noisedata50[i],length(noisedata50[i]),dt[i]), NUTS(0.65), 1000) for i in eachindex(dt)]
return oupndata

end

sample_rate = 0.1:0.1:0.1
mean_total = [Threads.Atomic{Float64}(0) for i in sample_rate]
std_total = [Threads.Atomic{Float64}(0) for i in sample_rate]
@threads for i in 1:10
    data = sample_model(sample_rate, 0.2, oupn_mult)
    for j in eachindex(mean_total)
        atomic_add!(mean_total[j], mean(data[j][:noise_ampl]))
        atomic_add!(std_total[j], std(data[j][:noise_ampl]))
    end
end

mean_total2 = [i[] for i in mean_total]/length(mean_total)
std_total2 = [i[] for i in std_total]/length(std_total)

plot(sample_rate, mean_total2, yerr = std_total2, title = "Thermal noise amplitude", xlabel = "sample rate in units of the relaxation time")

savefig("thermal.png")
