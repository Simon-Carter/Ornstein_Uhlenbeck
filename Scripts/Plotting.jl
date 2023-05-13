# Julia plotting script for the HMC models and EM algorithm

using CSV
using DataFrames
using Plots
using LaTeXStrings

a = CSV.read("../Data/mixed_ratio", DataFrame)

b = CSV.read("../Data/mixed_ratio_1.0", DataFrame)

c = CSV.read("../Data/mixed_ratio_2.0", DataFrame)

d = CSV.read("../Data/mixed_ratio_4.0", DataFrame)

e = CSV.read("../Data/mixed_ratio_0.25", DataFrame)

f= CSV.read("../Data/mixed_ratio_0.05", DataFrame)

g= CSV.read("../Data/mult", DataFrame)

helmut_data = CSV.read("../Data/EMdata.csv", DataFrame)
mcmcmean_data = CSV.read("../Data/MCMCmean.csv", DataFrame)

default(plot_titlefont=font(22, "Computer Modern"),  titlefont=font(19, "Computer Modern"), guidefont=font(16, "Computer Modern"), legendfont=font(12, "Computer Modern"), fontfamily="Computer Modern",
foreground_color_legend=nothing, background_color_legend=nothing)

plot(a.ratio_mulitplicative_thermal, a.mean_amplitude)
plot!(b.ratio_mulitplicative_thermal, b.mean_amplitude)
plot!(c.ratio_mulitplicative_thermal, c.mean_amplitude)
plot!(d.ratio_mulitplicative_thermal, d.mean_amplitude)
plot(e.ratio_mulitplicative_thermal, e.mean_amplitude)
plot(f.ratio_mulitplicative_thermal, f.mean_amplitude)




time_m_sum = plot(a.ratio_mulitplicative_thermal[1:50], a.mean_amplitude[1:50], linecolor="grey50", thickness_scaling = 2, size=(1000,600), ribbon=a.std_ampl[1:50], color = "grey50", fillalpha=.5, 
label="Amplitude estimation", legend=:right, linewidth=2, ylim=(0,1.4))
ylabel!("Amplitude")
annotate!(1, 0.3, text("0.5𝜏", 26,"Computer Modern"))
xlabel!(L"\sigma_m / \sigma_t")
hline!([1], linestyle=:dash, linecolor="grey20", label="True amplitude", linewidth=2)
savefig(time_m_sum, "../Plots/ratio_0.5.png")

time_m_sum = plot(d.ratio_mulitplicative_thermal[1:50], d.mean_amplitude[1:50], linecolor="grey50", thickness_scaling = 2, size=(1000,600), ribbon=d.std_ampl, color = "grey50", fillalpha=.5, 
label="Amplitude estimation", legend=:right, linewidth=2, ylim=(0,1.4))
ylabel!("Amplitude")
xlabel!(L"\sigma_m / \sigma_t")
annotate!(0.75, 0.4, text("4𝜏", 26,"Computer Modern"))
hline!([1], linestyle=:dash, linecolor="grey20", label="True amplitude", linewidth=2)
savefig(time_m_sum, "../Plots/ratio_4.png")

time_m_sum = plot(e.ratio_mulitplicative_thermal, e.mean_amplitude, linecolor="grey50", thickness_scaling = 2, size=(1000,600), ribbon=e.std_ampl, color = "grey50", fillalpha=.5, 
label="Amplitude estimation", legend=:right, linewidth=2, ylim=(0,1.4))
ylabel!("Amplitude")
xlabel!(L"\sigma_m / \sigma_t")
annotate!(1, 0.5, text("0.2𝜏", 26,"Computer Modern"))
hline!([1], linestyle=:dash, linecolor="grey20", label="True amplitude", linewidth=2)
savefig(time_m_sum, "../Plots/ratio_0.2.png")

time_m_sum = plot(f.ratio_mulitplicative_thermal, f.mean_amplitude, linecolor="grey50", thickness_scaling = 2, size=(1000,600), ribbon=f.std_ampl, color = "grey50", fillalpha=.5, 
label="Amplitude estimation", legend=:topright, linewidth=2, ylim=(0,1.6))
ylabel!("Amplitude")
xlabel!(L"\sigma_m / \sigma_t")
annotate!(1, 0.5, text("0.05𝜏", 26,"Computer Modern"))
hline!([1], linestyle=:dash, linecolor="grey20", label="True amplitude", linewidth=2)
savefig(time_m_sum, "../Plots/ratio_0.05.png")


time_m_sum = plot(g.ratio_deltat_tau1, g.mean_amplitude1, linecolor="grey50", thickness_scaling = 2, size=(1000,600), ribbon=g.std_ampl1, color = "grey50", fillalpha=.5, 
label="Amplitude estimation", legend=:best, linewidth=2, ylim=(0,2.6))
ylabel!("Parameters")
xlabel!(L"\Delta t / \tau")
hline!([1], linestyle=:dash, linecolor="grey20", label="True amplitude and tau", linewidth=2)
plot!(g.ratio_deltat_tau1,g.mean_tau1, linecolor="grey80", linewidth=2, ribbon=g.std_tau1, color = "grey80",  label="Tau estimation",)
savefig(time_m_sum, "../Plots/mult_fail.png")

time_m_sum = plot(helmut_data.dt, helmut_data.A_mean, linecolor="grey50", thickness_scaling = 2, size=(1000,800), ribbon=helmut_data.dA_mean, color = "grey50", fillalpha=.5, 
label="Amplitude estimation", legend=:topleft, linewidth=2, ylim=(0,3.0), xlim=(0,3.0))
ylabel!("Parameter value")
xlabel!(L"\Delta t / \tau")
hline!([1], linestyle=:dash, linecolor="grey20", label="True Parameters", linewidth=2)
plot!(helmut_data.dt,helmut_data.tau_mean, linecolor="grey70", linewidth=2, ribbon=helmut_data.dtau_mean, color = "grey70", alpha=0.2, linealpha=0.2,  label="𝜏 estimation", foreground_color_legend = nothing)
savefig(time_m_sum, "../Plots/EM.png")

time_m_sum = plot(mcmcmean_data.dt, mcmcmean_data.A, linecolor="grey50", thickness_scaling = 2, size=(1000,800), ribbon=mcmcmean_data.dA, color = "grey50", fillalpha=.5, 
label="Amplitude estimation", legend=:topleft, linewidth=2, ylim=(0,3.0), xlim=(0,3.0))
ylabel!("Parameter value")
xlabel!(L"\Delta t / \tau")
hline!([1], linestyle=:dash, linecolor="grey20", label="True Parameters", linewidth=2)
plot!(mcmcmean_data.dt,mcmcmean_data.tau, linecolor="grey70", linewidth=2, ribbon=mcmcmean_data.dtau, color = "grey70", alpha=0.2, linealpha=0.2,  label="𝜏 estimation", foreground_color_legend = nothing)
savefig(time_m_sum, "../Plots/MCMC.png")