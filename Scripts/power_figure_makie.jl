### generates the test data for the spectra plots and the corresponding time series
using DSP, Statistics, Random, Distributions, DifferentialEquations, Measures
using LaTeXStrings
#pythonplot()


Makie.set_theme!(fonts = (; regular = "sans roman serif"))

#parameters
thermal_std = 2
ratio=1
dt = 0.01


#first calculate all the needed Data
μ = 0.0
σ = sqrt(2)
Θ = 1.0
W = OrnsteinUhlenbeckProcess(Θ,μ,σ,0.0,1.0)
prob = NoiseProblem(W,(0.0,1000.0))
sol = [solve(prob;dt=dt).u for i in 1:2000]

#calculate the ou power
ou_avg_per = welch_pgram.(sol, fs=(1/dt))
ou_avg_pow = mean(power.(ou_avg_per))

#calculate variations of ou power
sol2 = [i.^2 for i in sol]
sol_abs = [abs.(i) for i in sol]
ou_avg_per2=  welch_pgram.(sol2, fs=(1/dt))
ou_avg_pow2 = mean(power.(ou_avg_per2))

#calcualte thermal first and second order
noise_t = copy(sol)
noise_t2 = copy(sol)
for i in eachindex(sol)
noise_t[i] = ([rand.(Normal.(0,thermal_std)) for j in (sol[i])])
end

noise_t2 = [i.^2 for i in noise_t]

t_avg_per = welch_pgram.(noise_t, fs=(1/dt))
t_avg_pow = mean(power.(t_avg_per))
t_avg_per2 = welch_pgram.(noise_t2, fs=(1/dt))
t_avg_pow2 = mean(power.(t_avg_per2))

#multiplicative noise calculation, 1st and second order
noise_m = copy(sol)
noise_m2 = copy(sol)
for i in eachindex(sol)
    noise_m2[i] = ([rand.(Normal.(0,ratio*thermal_std*sqrt.(abs.(j)))).^2 for j in (sol[i])])
    noise_m[i] = ([rand.(Normal.(0,ratio*thermal_std*sqrt.(abs.(j)))) for j in (sol[i])])
end


mult_avg_per = welch_pgram.((noise_m), fs=(1/dt))
mult_avg_pow = mean(power.(mult_avg_per))

mult_avg_per2 = welch_pgram.((noise_m2), fs=(1/dt))
mult_avg_pow2 = mean(power.(mult_avg_per2))

#calcualte the sum of multiplicatice and ou


mult_sum = [(sol[i] .+ noise_m[i]) for i in eachindex(sol)]
mult_sum2 = [i.^2 for i in mult_sum]

mult_sum_per = welch_pgram.((mult_sum), fs=(1/dt))
mult_sum_pow = mean(power.(mult_sum_per))

mult_sum_per2 = welch_pgram.((mult_sum2), fs=(1/dt))
mult_sum_pow2 = mean(power.(mult_sum_per2))


#calcualte the sum of thermal and ou

t_sum = [(sol[i] .+ noise_t[i]) for i in eachindex(sol)]
t_sum2 = [i.^2 for i in t_sum]

t_sum_per = welch_pgram.((t_sum), fs=(1/dt))
t_sum_pow = mean(power.(t_sum_per))

t_sum_per2 = welch_pgram.((t_sum2), fs=(1/dt))
t_sum_pow2 = mean(power.(t_sum_per2))


# plotting

#regular ou Process
#plot equivalent in CairoMakie
f = Figure()
ax = Axis(f[1, 1],
    xlabel = L"\text{Time (seconds)}",
    xlabelsize=45,
    ylabel = L"\text{Measurement}",
    ylabelsize=45,
    yticklabelsize  = 30,
    xticklabelsize  = 30
    )
lines!(ax, 0:0.01:3, sol[1][1:301], color="black", linewidth=3)

Makie.save("OU_signal.png", f)


# power spectrum
f = Figure()
ax = Axis(f[1, 1],
    xlabel = L"\text{Frequency (hz)}",
    xlabelsize=45,
    ylabel = L"\text{Power}",
    ylabelsize=45,
    yticklabelsize  = 30,
    xticklabelsize  = 30
    )
lines!(ax, ou_avg_per[1].freq[2:200], ou_avg_pow[2:200], color="black", linewidth=3, label="OU")

Makie.save("OU_power.png", f)

#OU power 2nd

f = Figure()
ax = Axis(f[1, 1],
    xlabel = L"\text{Frequency (hz)}",
    xlabelsize=45,
    ylabel = L"Power$^2$",
    ylabelsize=45,
    yticklabelsize  = 30,
    xticklabelsize  = 30
    )
lines!(ax, ou_avg_per2[1].freq[2:200], ou_avg_pow2[2:200], color="black", linewidth=3, label="OU")

Makie.save("OU_power2.png", f)



########### OU + Thermal ######################

f = Figure()
ax = Axis(f[1, 1],
    xlabel = L"\text{Time (seconds)}",
    xlabelsize=45,
    ylabel = L"\text{Measurement}",
    ylabelsize=45,
    yticklabelsize  = 30,
    xticklabelsize  = 30
    )
lines!(ax, 0:0.01:3, t_sum[1][1:301], color="black", linewidth=3)

Makie.save("OU_thermal_time.png", f)

# Power spectrum

f = Figure()
ax = Axis(f[1, 1],
    xlabel = L"\text{Frequency (hz)}",
    xlabelsize=45,
    ylabel = L"\text{Power}",
    ylabelsize=45,
    yscale = log10,
    yticklabelsize  = 30,
    xticklabelsize  = 30
    )
lines!(ax, t_sum_per[1].freq[2:600], t_sum_pow[2:600], color="grey50", label="OU + thermal", linewidth=3)
lines!(ax, t_avg_per[1].freq[2:600], t_avg_pow[2:600], color="black", label="Thermal noise", linewidth=3)
axislegend(ax, position = :rt,  framevisible = false, labelsize=30)
Makie.save("OU_thermal_power.png", f)

# power spectrum second order

f = Figure()
ax = Axis(f[1, 1],
    xlabel = L"\text{Frequency (hz)}",
    xlabelsize=45,
    ylabel = L"Power$^2$",
    ylabelsize=45,
    yscale = log10,
    yticklabelsize  = 30,
    xticklabelsize  = 30
    )
lines!(ax, t_sum_per2[1].freq[2:600], t_sum_pow2[2:600], color="grey50", label="OU + thermal", linewidth=3)
lines!(ax, t_avg_per2[1].freq[2:600], t_avg_pow2[2:600], color="black", label="Thermal noise", linewidth=3)
axislegend(ax, position = :rt,  framevisible = false, labelsize=30)
Makie.save("OU_thermal_power2.png", f)




############################ Multiplicative ##################################



f = Figure()
ax = Axis(f[1, 1],
    xlabel = L"\text{Time (seconds)}",
    xlabelsize=45,
    ylabel = L"\text{Measurement}",
    ylabelsize=45,
    yticklabelsize  = 30,
    xticklabelsize  = 30
    )
lines!(ax, 0:0.01:3, mult_sum[1][1:301], color="black", linewidth=2)

Makie.save("OU_mult_time.png", f)

# Power spectrum

f = Figure()
ax = Axis(f[1, 1],
    xlabel = L"\text{Frequency (hz)}",
    xlabelsize=45,
    ylabel = L"\text{Power}",
    ylabelsize=45,
    yscale = log10,
    yticklabelsize  = 30,
    xticklabelsize  = 30
    )
lines!(ax, mult_sum_per[1].freq[2:600], mult_sum_pow[2:600], color="grey50", label="OU + multiplicative", linewidth=3)
lines!(ax, mult_avg_per[1].freq[2:600], mult_avg_pow[2:600], color="black", label="Multiplicative noise", linewidth=3)
axislegend(ax, position = :rt,  framevisible = false, labelsize=30)
Makie.save("OU_mult_power.png", f)

# power spectrum second order

f = Figure()
ax = Axis(f[1, 1],
    xlabel = L"\text{Frequency (hz)}",
    xlabelsize=45,
    ylabel = L"\text{Power}$^2$",
    ylabelsize=45,
    yscale = log10,
    yticklabelsize  = 30,
    xticklabelsize  = 30
    )
lines!(ax, mult_sum_per2[1].freq[2:800], mult_sum_pow2[2:800], color="grey50", label="OU + multiplicative", linewidth=3)
lines!(ax, mult_sum_per2[1].freq[2:800], mult_avg_pow2[2:800], color="black", label="Multiplicative noise", linewidth=3)
axislegend(ax, position = :rt,  framevisible = false, labelsize=30)
Makie.save("OU_mult_power2.png", f)
