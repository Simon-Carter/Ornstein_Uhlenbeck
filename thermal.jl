using DifferentialEquations
using Plots, StatsPlots
using Distributions, Random, FFTW
using Turing, ReverseDiff, Memoization, LinearAlgebra

# Ornstein-Uhlenbeck process with added Gaussian noise
@model oupn(rn,T,delta_t,::Type{R}=Vector{Float64}) where {R} = begin
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

#initial genreated data
μ = 0.0
σ = sqrt(2)
Θ = 1.0
dt = 0.05:0.05:1.5
W = OrnsteinUhlenbeckProcess(Θ,μ,σ,0.0,1.0)
prob = NoiseProblem(W,(0.0,50.0))
sol = [solve(prob;dt=i).u for i in dt]

#generate noise
noisedata50 = [(sol[i] .+ [rand.(Normal.(0,0.2*sqrt.(abs.(sol[i][j])))) for j in 1:length(sol[i])]) for i in 1:30]
plot(sol[1])
plot!(noisedata50[1])

Turing.setadbackend(:reversediff)
Turing.setrdcache(true)

@time oupndata = [sample(oupn(noisedata50[i],length(noisedata50[i]),dt[i]), NUTS(0.65), 3000) for i in 1:30]

global mean_total = mean_total + [mean(oupndata[i][:noise_ampl]) for i in 1:30]
global std_total =  std_total + [std(oupndata[i][:noise_ampl]) for i in 1:30]
#end

