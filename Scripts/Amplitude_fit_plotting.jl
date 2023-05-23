### generates the test data for the spectra plots and the corresponding time series
using DSP, Statistics, Random, Distributions, DifferentialEquations, Measures
using LaTeXStrings
using Plots
#pythonplot()


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


##test start

default(plot_titlefont=font(16, "Computer Modern"),  titlefont=font(14, "Computer Modern"), fontfamily="Computer Modern")

#test end




#Calculate the plots

#regular ou Process
time_ou = plot(0:0.01:3, sol[1][1:301], legend=false, label="OU", linecolor="grey50", linewidth=2)
title!("\nTime series")
ylabel!("Measurement")
xlabel!("Time, s")
power1_ou = plot(ou_avg_per[1].freq[2:200], ou_avg_pow[2:200], linecolor="grey50", legend=false, linewidth=3, ylim=(0,4))
ylabel!("Power")
xlabel!("Frequency, Hz")
title!("\nPower spectrum")
power2_ou = plot(ou_avg_per2[1].freq[2:200], ou_avg_pow2[2:200], linecolor="grey50",legend=false, linewidth=3, ylim=(0,4))
ylabel!(L"\textrm{Power^2}")
xlabel!("Frequency, Hz")
title!("\nPower spectrum (second order)")
ou_comb = plot(time_ou, power1_ou, power2_ou, size=(2400,600), margin=5mm, plot_title="Ornstein-Uhlenbeck Process", layout=(1,3),legend=false,
plot_titlevspan=0.1, thickness_scaling = 2)
savefig(ou_comb, "../Plots/ou_comb.png")

#combine thermal pure and thermal ou
time_t_sum = plot(0:0.01:3, t_sum[1][1:301], legend=false, label="OU", linecolor="grey50", linewidth=2)
title!("\nTime series")
ylabel!("Measurement")
xlabel!("Time, s")
power1_t_sum = plot(t_sum_per[1].freq[2:600], t_sum_pow[2:600], linecolor="grey50", yaxis=:log, ylim=(10^(-2.0),10^(1)), 
foreground_color_legend = nothing, background_color_legend = nothing, label="OU + thermal", linewidth=3)
plot!(t_avg_per[1].freq[2:600], t_avg_pow[2:600], linecolor="grey25", label="Thermal noise", linewidth=3)
ylabel!("Power")
xlabel!("Frequency, Hz")
title!("\nPower spectrum")
power2_t_sum = plot(t_sum_per2[1].freq[2:600], t_sum_pow2[2:600], linecolor="grey50", yaxis=:log, ylim=(10^(-1.0),10^(1)), 
foreground_color_legend = nothing, background_color_legend = nothing, label="OU + thermal", linewidth=3)
plot!(t_avg_per2[1].freq[2:600], t_avg_pow2[2:600], linecolor="grey25", label="Thermal noise", linewidth=3)
ylabel!(L"\textrm{Power^2}")
xlabel!("Frequency, Hz")
title!("\nPower spectrum (second order)")
t_sum_comb = plot(time_t_sum, power1_t_sum, power2_t_sum, size=(2400,600), margin=5mm, plot_title="Ornstein-Uhlenbeck Process with Thermal Noise", 
layout=(1,3), plot_titlevspan=0.1, thickness_scaling = 2)
savefig(t_sum_comb, "../Plots/2t_sum_comb.png")


#combine multiplicative pure and multiplicative ou
time_m_sum = plot(0:0.01:3, mult_sum[1][1:301], legend=false, label="OU", linecolor="grey50", linewidth=2)
title!( "\nTime series\n")
ylabel!("Measurement")
xlabel!("Time, s")
power1_m_sum = plot(mult_sum_per[1].freq[2:600], mult_sum_pow[2:600], linecolor="grey50", yaxis=:log, ylim=(10^(-2.0),10^(1)), 
label="OU + Multiplicative", foreground_color_legend = nothing, background_color_legend = nothing, linewidth=3)
plot!(mult_avg_per[1].freq[2:600], mult_avg_pow[2:600], linecolor="grey25", ylim=(10^(-2.0),10^(1)), label="Multiplicative noise", linewidth=3)
ylabel!("Power")
xlabel!("Frequency, Hz")
title!("\nPower spectrum\n")
power2_m_sum = plot(mult_sum_per2[1].freq[2:800], mult_sum_pow2[2:800], linecolor="grey50", yaxis=:log, ylim=(10^(-0.5),10^(1.5)), 
label="OU + Multiplicative", foreground_color_legend = nothing, background_color_legend = nothing, linewidth=3)
plot!(mult_sum_per2[1].freq[2:800], mult_avg_pow2[2:800], linecolor="grey25", yaxis=:log, label="Multiplicative noise", linewidth=3)
ylabel!(L"\textrm{Power^2}")
xlabel!("Frequency, Hz")
title!("\nPower spectrum (second order)\n", topmargin=5mm)
m_sum_comb = plot(time_m_sum, power1_m_sum, power2_m_sum, size=(2400,600), margin=5mm, layout=(1,3), 
plot_title="Ornstein-Uhlenbeck Process with Multiplicative Noise", plot_titlevspan=0.1, thickness_scaling = 2)
savefig(m_sum_comb, "../Plots/2m_sum_comb.png")