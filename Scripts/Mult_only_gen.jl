include("../model_generation.jl")

#generate data for a model with no thermal noise and only multiplicative noise

sample_rate = 0.05:0.05:1
num_cycles = 1
Tnoise = 0.2
ratio = 1
mean_total = [Threads.Atomic{Float64}(0) for i in sample_rate]
std_total = [Threads.Atomic{Float64}(0) for i in sample_rate]
data_total = Array{Chains}(undef, num_cycles, length(sample_rate))
@threads for j in 1:num_cycles
    noisedata = data_generation(sample_rate,0.2,model="m")
    data_total[j,:] = @time [sample(oupn_m(noisedata[i],length(noisedata[i]),sample_rate[i]), NUTS(0.65), 1000) for i in eachindex(sample_rate)]
end


FileIO.save("mult_0.05_4.jld2","data",data_total)