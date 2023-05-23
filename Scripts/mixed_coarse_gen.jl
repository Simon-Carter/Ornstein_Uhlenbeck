include("../model_generation.jl")

#generted data multiplicative and thermal noise for different sample rates

sample_rate = 0.1:0.1:2
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


# FileIO.save("../Data/mixed_points1500_rate0.1:0.1:2_mcmc1000_20cycles.jld2","data",data_total)