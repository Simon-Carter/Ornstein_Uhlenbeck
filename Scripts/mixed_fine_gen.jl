include("../model_generation.jl")

#generates data for multiplicative and thermal model for a fine sample rate

sample_rate = 0.005:0.005:1
num_cycles = 1
Tnoise = 0.2
ratio = 1
mean_total = [Threads.Atomic{Float64}(0) for i in sample_rate]
std_total = [Threads.Atomic{Float64}(0) for i in sample_rate]
data_total = Array{Chains}(undef, num_cycles, length(sample_rate))
@threads for j in 1:num_cycles
    noisedata = data_generation(sample_rate,0.2,model="mt")
    data_total[j,:] = @time [sample(oupn_mt(noisedata[i],length(noisedata[i]),sample_rate[i], 1), NUTS(0.65), 1000) for i in eachindex(sample_rate)]
end


#FileIO.save("mixed_ratio_0.1:0.1:10_tnoise_0.2_numcycles_1.jld2","data",data_total)