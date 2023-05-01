using DifferentialEquations
using Plots, StatsPlots
using Distributions, Random, FFTW
using Turing, ReverseDiff, Memoization
using LinearAlgebra:Diagonal
using Base.Threads
using FileIO, JLD2
using LaTeXStrings
using DataFrames
using CSV

#thermal plotting

x_axis = 0.1:0.1:5
true_mean_noise_t = (map(x-> (0.2 ./(1 .+ x)), x_axis))

function data_parse(data)
    mean_data = vec(mean(data, dims=1))
    return mean_data
end

data = FileIO.load("mixed_ratio_0.05.jld2", "data")

mean_noise_t = mean.(map(x-> x[:noise_ampl_t], data))'
std_noise_t = std.(map(x-> x[:noise_ampl_t], data))'

mean_ampl = mean.(map(x-> x[:ampl], data))'
mean_tau = mean.(map(x-> x[:tau], data))'
std_ampl = std.(map(x-> x[:ampl], data))'
std_tau = std.(map(x-> x[:tau], data))'

data = Dict(:ratio_mulitplicative_thermal => x_axis,
:mean_noise_thermal => vec(mean_noise_t),
:std_noise_thermal => vec(std_noise_t),
:mean_amplitude => vec(mean_ampl),
:mean_tau => vec(mean_tau),
:std_ampl => vec(std_ampl),
:std_tau => vec(std_tau),
:true_mean_noise_t => true_mean_noise_t)

CSV_data = DataFrame(data)

CSV.write("mixed_ratio_0.05", CSV_data)

#= Plotting stuff

plot(x_axis, data_parse(std_noise_t), label="Thermal noise", xlabel="ratio of multiplicative to thermal noise", ylabel="STD", title="Mixed Gaussian STD")

#plot!(x_axis, data_parse(std_tau), label="Tau")

plot!(x_axis, data_parse(std_ampl), label="Amplitude")

savefig("mixed_ratio_std.png")

plot(x_axis, data_parse(mean_noise_t), label="Thermal noise", xlabel="ratio of multiplicative to thermal noise", ylabel="Mean", title="Mixed Gaussian Mean")
#plot!(x_axis, data_parse(mean_tau), label="Tau")
plot!(x_axis, data_parse(mean_ampl), label="Amplitude")
plot!(x_axis, true_mean_noise_t, label="True Noise", line=(:dot, 2))
true_tau = [1.0 for i in x_axis]
plot!(x_axis, true_tau, label="True Amplitude", line=(:dot, 2))

savefig("mixed_ratio_mean.png")

=#