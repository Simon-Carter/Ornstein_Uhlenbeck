include("../model_generation.jl")

sample_rate = 0.005:0.005:1
num_cycles = 2
Tnoise = 0.2
ratio = 1
mean_total = [Threads.Atomic{Float64}(0) for i in sample_rate]
std_total = [Threads.Atomic{Float64}(0) for i in sample_rate]
data_total = Array{Chains}(undef, num_cycles, length(sample_rate))
@threads for j in 1:num_cycles
    noisedata = data_generation(sample_rate,0.2,model="t")
    data_total[j,:] = @time [sample(oupn_thermal(noisedata[i],length(noisedata[i]),sample_rate[i]), NUTS(0.65), 1000) for i in eachindex(sample_rate)]
end


#FileIO.save("thermal.jld2","data",data_total)