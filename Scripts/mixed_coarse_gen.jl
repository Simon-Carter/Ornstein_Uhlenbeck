include("../model_generation.jl")

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


#FileIO.save("mixed_points1500_rate0.1:0.1:2_mcmc1000_20cycles.jld2","data",data_total)



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




