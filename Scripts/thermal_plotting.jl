using DifferentialEquations
using Plots, StatsPlots
using Distributions, Random, FFTW
using Turing, ReverseDiff, Memoization
using LinearAlgebra:Diagonal
using Base.Threads
using FileIO, JLD2
using LaTeXStrings

#parse and plot the thermal data generated from thermal_gen.jl

x_axis = 0.05:0.05:1

function data_parse(data)
    mean_data = vec(mean(data, dims=1))
    return mean_data
end

function gen_name(name, data)
    return ["$name$i" for i in 1:size(data, 2)]
end

#input the file that you saved in the gen script
data = FileIO.load("mult_0.05_4.jld2", "data")

mean_noise_t = mean.(map(x-> x[:noise_ampl_m], data))
std_noise_t = std.(map(x-> x[:noise_ampl_m], data))

mean_ampl = mean.(map(x-> x[:ampl], data))
mean_tau = mean.(map(x-> x[:tau], data))
std_ampl = std.(map(x-> x[:ampl], data))
std_tau = std.(map(x-> x[:tau], data))

col_data = hcat(x_axis, mean_noise_t', std_noise_t', mean_ampl', mean_tau', std_ampl', std_tau')

col = [gen_name("ratio_deltat_tau", x_axis), gen_name("mean_noise_thermal", mean_noise_t'),
gen_name("std_noise_thermal", std_noise_t'), gen_name("mean_amplitude", mean_ampl'),
gen_name("mean_tau", mean_tau'), gen_name("std_ampl", std_ampl'), gen_name("std_tau", std_tau')]

col = vcat(col...)

frame_data = DataFrame(col_data,:auto)

rename!(frame_data, col)

CSV.write("mult", frame_data)

#plotting stuff

plot(x_axis, data_parse(std_noise_t), label="Thermal noise", xlabel=L"\frac{Δt}{τ}", ylabel="STD", title="Thermal Gaussian STD")

#plot!(x_axis, data_parse(std_tau), label="Tau")

plot!(x_axis, data_parse(std_ampl), label="Amplitude")

savefig("thermal_std.png")

plot(x_axis, data_parse(mean_noise_t), label="Thermal noise", xlabel=L"\frac{Δt}{τ}", ylabel="Mean", title="Thermal Gaussian Mean")
#plot!(x_axis, data_parse(mean_tau), label="Tau")
plot!(x_axis, data_parse(mean_ampl), label="Amplitude")
true_mean_noise_t = [0.2 for i in x_axis]
plot!(x_axis, true_mean_noise_t, label="True Noise", line=(:dot, 2))
true_tau = [1.0 for i in x_axis]
plot!(x_axis, true_tau, label="True Amplitude", line=(:dot, 2))

savefig("thermal_mean.png")