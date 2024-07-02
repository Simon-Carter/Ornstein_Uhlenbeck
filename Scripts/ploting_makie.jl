using LaTeXStrings
using DataFrames
using CSV
using Serialization
using CairoMakie
using Turing
using Statistics

include("./utils.jl")

#convert to dataframe to allow for  easy plotting
data = quick_open("../Data/true_value_redo.happy")

#convert to dataframe to allow for  easy plotting
df_data = DataFrame(quick_open("../Data/true_value_redo.happy")["prediction"])[:,3:1002]
mean_data=mean.(eachcol(df_data))



#= Depreciated, switched from plots.jl to makie
plot(quick_open("true_value_redo.happy")["noise_free"], linecolor="grey80", label="True amplitude", linewidth=2)
plot!(mean_data, linecolor="grey50", thickness_scaling = 2, size=(1000,600), color = "grey50", fillalpha=.5, 
label="Predicted Signal", legend=:topright, linewidth=2)
ylabel!("Amplitude")
xlabel!("Time")
savefig("dummy.png")
=#


#Set Default macki font
#Makie.set_theme!(fonts = (; regular = "sans roman serif"))


#experimentation on font size. For full sized graphs, 15 point font fits matches default,
    # Half size graphs go for 30
    # 1/3 go for 45pt etc

f = Figure()
ax = Axis(f[1, 1],
    xlabel = "Time (seconds)",
    xlabelsize=15,
    ylabel = "Amplitude",
    ylabelsize=15
    )
lines!(ax, 0:0.1:49.9, quick_open("../Data/true_value_redo.happy")["noise_free"][1:500], color="grey50", label = "True amplitude", linewidth=2)
lines!(ax, 0:0.1:49.9, mean_data[1:500], color="black", linewidth=2, label="Predicted signal")
axislegend(ax, position = :rt,  framevisible = false)
f
Makie.save("best_fit.png", f)




#plot helmut csv
t = collect(1:100)*0.1
data = CSV.read("../Data/EM.csv",  DataFrame)

eml_fig = Figure()
eml_ax = Axis(eml_fig[1, 1],
    xlabel = "Time (seconds)",
    xlabelsize=15,
    ylabel = "Amplitude",
    ylabelsize=15
    )

scatter!(eml_ax,t, data[!,:xn][1:100], color="grey90", label="OU + thermal", linewidth=2)
lines!(eml_ax,t, data[!,:x][1:100],linecolor="black", label="OU", linewidth=1.2)
eml_fig
plot!(t,data[1:100,1],ribbon = data[1:100,2], linecolor="grey40", fillcolor="grey50", label="EM prediction", linewidth=2)
plot!(xticks = (0:10:100, 0:10))
ylabel!("Amplitude")
xlabel!("Time (seconds)")


#plot equivalent in CairoMakie
data = CSV.read("EM.csv",  DataFrame)
f = Figure()
ax = Axis(f[1, 1],
    xlabel = L"\text{Time (seconds)}",
    xlabelsize=20,
    ylabel = L"\text{Amplitude}",
    ylabelsize=20
    )
band!(ax, 0:0.1:9.9, data[1:100,1] + data[1:100,2], data[1:100,1] - data[1:100,2],  color="grey60")
lines!(ax, 0:0.1:9.9, data[1:100,1],color="grey40", label="EM prediction", linewidth=3)
lines!(ax, 0:0.1:9.9, data[1:100,3], color = :black, label="OU", linewidth=3, linestyle=Linestyle([0.7, 1.4, 2.1, 2.8]))
scatter!(ax, 0:0.1:9.9, data[1:100,4], color="grey90", label="OU + thermal", strokecolor=:black, strokewidth=1)
axislegend(ax, position = (0.4, 1.0),  framevisible = false)
f
Makie.save("EM_fit.png", f)



# Julia plotting script for the HMC models and EM algorithm, used in paper

using CSV
using DataFrames
using LaTeXStrings

a = CSV.read("../Data/mixed_ratio", DataFrame)

b = CSV.read("../Data/mixed_ratio_1.0", DataFrame)

c = CSV.read("../Data/mixed_ratio_2.0", DataFrame)

d = CSV.read("../Data/mixed_ratio_4.0", DataFrame)

e = CSV.read("../Data/mixed_ratio_0.25", DataFrame)

fd= CSV.read("../Data/mixed_ratio_0.05", DataFrame)

g= CSV.read("../Data/mult", DataFrame)

helmut_data = CSV.read("../Data/EMdata.csv", DataFrame)
mcmcmean_data = CSV.read("../Data/MCMCmean.csv", DataFrame)


f = Figure()
ax = Axis(f[1, 1],
    xlabel = L"\Delta t /\tau",
    xlabelsize=30,
    ylabel = "Parameter value",
    ylabelsize=30,
    yticklabelsize  = 20,
    xticklabelsize  = 20
    )
band!(ax, helmut_data.dt,  helmut_data.A_mean - helmut_data.dA_mean, helmut_data.A_mean + helmut_data.dA_mean, color="grey70")
lines!(ax, helmut_data.dt, helmut_data.A_mean,color="black", linewidth=3, label="Amplitude estimate")

band!(ax, helmut_data.dt,helmut_data.tau_mean - helmut_data.dtau_mean, helmut_data.tau_mean + helmut_data.dtau_mean, color="grey70")
lines!(ax, helmut_data.dt,helmut_data.tau_mean ,color="grey70", linewidth=3, label="𝜏 estimate")

hlines!(ax, 1.0; linestyle=:dash, linewidth=3, color = "grey30", label="True value")

axislegend(ax, position = :lt,  framevisible = false, labelsize = 20)

xlims!(high = 3.5)
f
Makie.save("EM.png", f)





f = Figure()
ax = Axis(f[1, 1],
    xlabel = L"\Delta t/\tau",
    xlabelsize=30,
    ylabel = "Parameter value",
    ylabelsize=30,
    yticklabelsize  = 20,
    xticklabelsize  = 20
    )
band!(ax, mcmcmean_data.dt,  mcmcmean_data.A - mcmcmean_data.dA, mcmcmean_data.A + mcmcmean_data.dA, color="grey50", alpha=0.5)
band!(ax, mcmcmean_data.dt,mcmcmean_data.tau - mcmcmean_data.dtau, mcmcmean_data.tau + mcmcmean_data.dtau, color="grey70", alpha=0.5)

lines!(ax, mcmcmean_data.dt, mcmcmean_data.A,color="black", linewidth=3, label="Amplitude estimate")
lines!(ax, mcmcmean_data.dt,mcmcmean_data.tau ,color="grey40", linewidth=3, label="𝜏 estimate")

hlines!(ax, 1.0; linestyle=:dash, linewidth=3, color = "grey30", label="True value")

axislegend(ax, position = :lt,  framevisible = false, labelsize = 20)
f
Makie.save("MCMC.png", f)


# mult fail plot

f = Figure()
ax = Axis(f[1, 1],
    xlabel = L"\Delta t/\tau",
    xlabelsize=30,
    ylabel = "Parameter value",
    ylabelsize=30,
    yticklabelsize  = 20,
    xticklabelsize  = 20
    )
band!(ax, g.ratio_deltat_tau1, g.mean_amplitude1 - g.std_ampl1, g.mean_amplitude1 + g.std_ampl1, color="grey50", alpha=0.5)
band!(ax, g.ratio_deltat_tau1, g.mean_tau1 - g.std_tau1, g.mean_tau1 + g.std_tau1, color="grey80", alpha=0.5)

lines!(ax, g.ratio_deltat_tau1, g.mean_amplitude1,color="black", linewidth=3, label="Amplitude estimate")
lines!(ax, g.ratio_deltat_tau1,g.mean_tau1 ,color="grey50", linewidth=3, label="𝜏 estimate")

hlines!(ax, 1.0; linestyle=:dash, linewidth=3, color = "grey30", label="True value")

axislegend(ax, position = (0.25, 1.0),  framevisible = false, labelsize = 20)
f
Makie.save("mult_fail.png", f)




# The other stuff, kill me this took so long to do,

#0.5

f = Figure()
ax = Axis(f[1, 1],
    xlabel = L"\sigma_m / \sigma_t",
    xlabelsize=30,
    ylabel = "OU Amplitude estimate",
    ylabelsize=28,
    yticklabelsize  = 24,
    xticklabelsize  = 24
    )


band!(ax, a.ratio_mulitplicative_thermal[1:50], a.mean_amplitude[1:50]+a.std_ampl[1:50],  a.mean_amplitude[1:50]-a.std_ampl[1:50], color="grey50", alpha=0.5)

lines!(ax, a.ratio_mulitplicative_thermal[1:50], a.mean_amplitude[1:50], color="black", linewidth=3, label="OU Ampl estimate")

hlines!(ax, 1.0; linestyle=:dash, linewidth=3, color = "grey30", label="True value")

axislegend(ax, position = (0.8, 0.7),  framevisible = false, labelsize = 24)

text!(0.25, 0.2 ,text=L"\tau = 0.5", fontsize=36)
f
Makie.save("ratio_0.5.png", f)

#4


f = Figure()
ax = Axis(f[1, 1],
    xlabel = L"\sigma_m / \sigma_t",
    xlabelsize=30,
    ylabel = "OU Amplitude estimate",
    ylabelsize=28,
    yticklabelsize  = 24,
    xticklabelsize  = 24
    )


band!(ax, d.ratio_mulitplicative_thermal[1:50], d.mean_amplitude[1:50]+d.std_ampl[1:50],  d.mean_amplitude[1:50]-d.std_ampl[1:50], color="grey50", alpha=0.5)

lines!(ax, d.ratio_mulitplicative_thermal[1:50], d.mean_amplitude[1:50], color="black", linewidth=3, label="OU Amplitude estimate")

hlines!(ax, 1.0; linestyle=:dash, linewidth=3, color = "grey30", label="True value")

axislegend(ax, position = (0.8, 0.7),  framevisible = false, labelsize = 24)

text!(0.25, 0.2 ,text=L"\tau = 4.0", fontsize=36)
f
Makie.save("ratio_4.png", f)

#0.2


f = Figure()
ax = Axis(f[1, 1],
    xlabel = L"\sigma_m / \sigma_t",
    xlabelsize=30,
    ylabel = "OU Amplitude estimate",
    ylabelsize=28,
    yticklabelsize  = 24,
    xticklabelsize  = 24
    )


band!(ax, e.ratio_mulitplicative_thermal[1:50], e.mean_amplitude[1:50]+e.std_ampl[1:50],  e.mean_amplitude[1:50]-e.std_ampl[1:50], color="grey50", alpha=0.5)

lines!(ax, e.ratio_mulitplicative_thermal[1:50], e.mean_amplitude[1:50], color="black", linewidth=3, label="OU Ampl estimate")

hlines!(ax, 1.0; linestyle=:dash, linewidth=3, color = "grey30", label="True value")

axislegend(ax, position = (0.8, 0.7),  framevisible = false, labelsize = 24)

text!(0.25, 0.2 ,text=L"\tau = 0.2", fontsize=36)
f
Makie.save("ratio_0.2.png", f)


# 0.05


f = Figure()
ax = Axis(f[1, 1],
    xlabel = L"\sigma_m / \sigma_t",
    xlabelsize=30,
    ylabel = "OU Amplitude estimate}",
    ylabelsize = 28,
    yticklabelsize  = 24,
    xticklabelsize  = 24
    )


band!(ax, fd.ratio_mulitplicative_thermal[1:50], fd.mean_amplitude[1:50]+fd.std_ampl[1:50],  fd.mean_amplitude[1:50]-fd.std_ampl[1:50], color="grey50", alpha=0.5)

lines!(ax, fd.ratio_mulitplicative_thermal[1:50], fd.mean_amplitude[1:50], color="black", linewidth=3, label="OU Ampl estimate")

hlines!(ax, 1.0; linestyle=:dash, linewidth=3, color = "grey30", label="True value")

axislegend(ax, position = :rt,  framevisible = false, labelsize = 24)

text!(0.25, 0.2 ,text=L"\tau = 0.05", fontsize=36)
f
Makie.save("ratio_0.05.png", f)

nuts_data = CSV.read("../Data/nuts_chain.csv", DataFrame)
nuts_fig = Figure()
nuts_names = names(nuts_data)
nuts_labels = ["var OU", "var MN", "τ"]
for (i,name) in enumerate(nuts_names)
    ax = Axis(nuts_fig[i, 1]; ylabel=nuts_labels[i], ylabelsize=20)
    lines!(nuts_data[!,name],  color=:black)

#    hideydecorations!(ax; label=false)
    if i < length(nuts_names)
        hidexdecorations!(ax; grid=false)
    else
        ax.xlabel = "N"
    end
end
yspace = maximum(tight_yticklabel_spacing!, nuts_fig.content)
for ax in nuts_fig.content
    ax.yticklabelspace = yspace
end
nuts_fig

Makie.save("nuts_chains.png", nuts_fig)
