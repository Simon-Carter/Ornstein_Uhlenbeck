include("../model_generation.jl")


#generate data with a fixed sample rate but with various thermal to multiplicative noise ratios

sample_rate = [0.05]
ratio = 0.1:0.1:5
Tnoise = 0.2 ./(1 .+ ratio)
num_cycles = 1
mean_total = [Threads.Atomic{Float64}(0) for i in sample_rate]
std_total = [Threads.Atomic{Float64}(0) for i in sample_rate]
data_total = Array{Chains}(undef, length(ratio), length(sample_rate))
@threads for j in eachindex(ratio)
    noisedata = data_generation(sample_rate,Tnoise[j],ratio[j], model="mt")
    data_total[j,:] = @time [sample(oupn_mt(noisedata[i],length(noisedata[i]),sample_rate[i], ratio[j]), NUTS(0.65), 1000) for i in eachindex(sample_rate)]
end


    FileIO.save("mixed_ratio_0.05.jld2","data",data_total)



