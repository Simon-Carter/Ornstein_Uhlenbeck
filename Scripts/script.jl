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

default(plot_titlefont=font(16, "Times Bold"), guidefont=font(12, "Times Bold"),  titlefont=font(14, "Times Bold"))

#test end

#Calculate the plots

#regular ou Process
time_ou = plot(0:0.01:3, sol[1][1:301], legend=false, label="OU", linecolor="grey50")
title!("\nTime series")
ylabel!("Measurement")
xlabel!("Time, s")
power1_ou = plot(ou_avg_per[1].freq[2:200], ou_avg_pow[2:200], linecolor="grey50", legend=false, linewidth=2, ylim=(0,4))
ylabel!("Power")
xlabel!("Frequency, Hz")
title!("\nPower spectrum")
power2_ou = plot(ou_avg_per2[1].freq[2:200], ou_avg_pow2[2:200], linecolor="grey50",legend=false, linewidth=2, ylim=(0,4))
ylabel!("Power^2")
xlabel!("Frequency, Hz")
title!("\nPower spectrum (second order)")
ou_comb = plot(time_ou, power1_ou, power2_ou, size=(2200,600), margin=5mm, plot_title="Ornstein-Uhlenbeck Process", layout=(1,3),legend=false,
plot_titlevspan=0.1, thickness_scaling = 2)
savefig(ou_comb, "../Plots/ou_comb.png")
#=
#Thermal noise
time_t = plot(0:0.01:50, noise_t[1][1:5001], legend=false, label="OU", linecolor="grey50")
title!("\nTime series")
ylabel!("Measurement")
xlabel!("Time, s")
power1_t = plot(t_avg_per[1].freq[2:2000], t_avg_pow[2:2000], linecolor="grey50", legend=false, ylim=(0,0.1))
ylabel!("Power")
xlabel!("Frequency, Hz")
title!("\nPower spectrum")
power2_t = plot(t_avg_per2[1].freq[2:2000], t_avg_pow2[2:2000], linecolor="grey50",legend=false, ylim=(0,1))
ylabel!("Power^2")
xlabel!("Frequency, Hz")
title!("\nPower spectrum (second order)")
t_comb = plot(time_t, power1_t, power2_t, size=(1800,600), margin=10mm, plot_title="Thermal Noise", layout=(1,3),legend=false)
savefig(t_comb, "t_comb.png")


#Multiplicative noise
time_m = plot(0:0.01:50, noise_m[1][1:5001], legend=false, label="OU", linecolor="grey50")
title!("\nTime series")
ylabel!("Measurement")
xlabel!("Time, s")
power1_m = plot(mult_avg_per[1].freq[2:200], mult_avg_pow[2:200], linecolor="grey50", legend=false, ylim=(0,0.1))
ylabel!("Power")
xlabel!("Frequency, Hz")
title!("\nPower spectrum")
power2_m = plot(mult_avg_per2[1].freq[2:600], mult_avg_pow2[2:600], linecolor="grey50",legend=false, yaxis=:log, ylim=(10^(-1.0),10^(1.5)))
ylabel!("Power^2")
xlabel!("Frequency, Hz")
title!("\nPower spectrum (second order)")
m_comb = plot(time_m, power1_m, power2_m, size=(1800,600), margin=10mm, plot_title="Multiplicative noise Process", layout=(1,3),legend=false)
savefig(m_comb, "m_comb.png")


#Multiplicative + ou
time_m_sum = plot(0:0.01:3, mult_sum[1][1:301], legend=false, label="OU", linecolor="grey50")
title!("\nTime series")
ylabel!("Measurement")
xlabel!("Time, s")
power1_m_sum = plot(mult_sum_per[1].freq[2:600], mult_sum_pow[2:600], linecolor="grey50", legend=false, yaxis=:log, ylim=(10^(-2.0),10^(1)))
ylabel!("Power")
xlabel!("Frequency, Hz")
title!("\nPower spectrum")
power2_m_sum = plot(mult_sum_per2[1].freq[2:800], mult_sum_pow2[2:800], linecolor="grey50",legend=false, yaxis=:log, ylim=(10^(0.0),10^(1.5)))
ylabel!("Power^2")
xlabel!("Frequency, Hz")
title!("\nPower spectrum (second order)")
m_sum_comb = plot(time_m_sum, power1_m_sum, power2_m_sum, size=(1800,600), margin=10mm, plot_title="Ornstein-Uhlenbeck Process with Multiplicative Noise", layout=(1,3),legend=false)
savefig(m_sum_comb, "m_sum_comb.png")

#Thermal + ou

time_t_sum = plot(0:0.01:3, t_sum[1][1:301], legend=false, label="OU", linecolor="grey50")
title!("\nTime series")
ylabel!("Measurement")
xlabel!("Time, s")
power1_t_sum = plot(t_sum_per[1].freq[2:600], t_sum_pow[2:600], linecolor="grey50", legend=false, yaxis=:log, ylim=(10^(-2.0),10^(1)))
ylabel!("Power")
xlabel!("Frequency, Hz")
title!("\nPower spectrum")
power2_t_sum = plot(t_sum_per2[1].freq[2:600], t_sum_pow2[2:600], linecolor="grey50",legend=false, yaxis=:log, ylim=(10^(-1.0),10^(1)))
ylabel!("Power^2")
xlabel!("Frequency, Hz")
title!("\nPower spectrum (second order)")
t_sum_comb = plot(time_t_sum, power1_t_sum, power2_t_sum, size=(1800,600), margin=10mm, plot_title="Ornstein-Uhlenbeck Process with Thermal Noise", layout=(1,3),legend=false)
savefig(t_sum_comb, "t_sum_comb.png")

=#
#combine thermal pure and thermal ou
time_t_sum = plot(0:0.01:3, t_sum[1][1:301], legend=false, label="OU", linecolor="grey50")
title!("\nTime series")
ylabel!("Measurement")
xlabel!("Time, s")
power1_t_sum = plot(t_sum_per[1].freq[2:600], t_sum_pow[2:600], linecolor="grey50", yaxis=:log, ylim=(10^(-2.0),10^(1)), 
foreground_color_legend = nothing, background_color_legend = nothing, label="OU + thermal", linewidth=2)
plot!(t_avg_per[1].freq[2:600], t_avg_pow[2:600], linecolor="grey25", label="Thermal noise")
ylabel!(L"\mathbf{Power}")
xlabel!(L"\mathbf{Frequency, Hz}")
title!("\nPower spectrum")
power2_t_sum = plot(t_sum_per2[1].freq[2:600], t_sum_pow2[2:600], linecolor="grey50", yaxis=:log, ylim=(10^(-1.0),10^(1)), 
foreground_color_legend = nothing, background_color_legend = nothing, label="OU + thermal", linewidth=2)
plot!(t_avg_per2[1].freq[2:600], t_avg_pow2[2:600], linecolor="grey25", label="Thermal noise")
ylabel!(L"\mathbf{Power^2}")
xlabel!(L"\mathbf{Frequency, Hz}")
title!("\nPower spectrum (second order)")
t_sum_comb = plot(time_t_sum, power1_t_sum, power2_t_sum, size=(2200,600), margin=5mm, plot_title="Ornstein-Uhlenbeck Process with Thermal Noise", 
layout=(1,3), plot_titlevspan=0.1, thickness_scaling = 1.5)
savefig(t_sum_comb, "../Plots/2t_sum_comb.png")

#combine multiplicative pure and multiplicative ou
time_m_sum = plot(0:0.01:3, mult_sum[1][1:301], legend=false, label="OU", linecolor="grey50")
title!( "Time series", weight="bold")
ylabel!(L"\mathbf{Measurement}")
xlabel!(L"\mathbf{Time, s}")
power1_m_sum = plot(mult_sum_per[1].freq[2:600], mult_sum_pow[2:600], linecolor="grey50", yaxis=:log, ylim=(10^(-2.0),10^(1)), 
label="OU + Multiplicative", foreground_color_legend = nothing, background_color_legend = nothing, linewidth=2)
plot!(mult_avg_per[1].freq[2:600], mult_avg_pow[2:600], linecolor="grey25", ylim=(10^(-2.0),10^(1)), label="Multiplicative noise", linewidth=2)
ylabel!(L"\mathbf{Power}")
xlabel!(L"\mathbf{Frequency, Hz}")
title!("Power spectrum")
power2_m_sum = plot(mult_sum_per2[1].freq[2:800], mult_sum_pow2[2:800], linecolor="grey50", yaxis=:log, ylim=(10^(-0.5),10^(1.5)), 
label="OU + Multiplicative", foreground_color_legend = nothing, background_color_legend = nothing, linewidth=2)
plot!(mult_sum_per2[1].freq[2:800], mult_avg_pow2[2:800], linecolor="grey25", yaxis=:log, label="Multiplicative noise", linewidth=2)
ylabel!(L"\mathbf{Power^2}")
xlabel!(L"\mathbf{Frequency, Hz}")
title!("Power spectrum (second order)")
m_sum_comb = plot(time_m_sum, power1_m_sum, power2_m_sum, size=(2200,600), margin=5mm, layout=(1,3), 
plot_title="Ornstein-Uhlenbeck Process with Multiplicative Noise", plot_titlevspan=0.1, thickness_scaling = 1.5)
savefig(m_sum_comb, "../Plots/2m_sum_comb.png")

#legacy
#=
power_2nd = plot(mult_sum_per[1].freq[2:200], mult_sum_pow[2:200], linecolor="black", label="OU with multiplicative noise", foreground_color_legend = nothing, background_color_legend = nothing)
plot!(mult_avg_per[1].freq[2:200], mult_avg_pow[2:200], linecolor="grey55", label="Multiplicative noise")
plot!(ou_avg_per2[1].freq[2:200], ou_avg_pow2[2:200], linecolor="grey80", label="OU")
ylabel!("Power")
xlabel!("Frequency, Hz")
title!("\nPower spectrum (second order)")

time_mult = plot(0:0.01:50, mult_sum[1][1:5001], legend=:topright, label="OU with multiplicative noise", linecolor="grey80", foreground_color_legend = nothing, background_color_legend = nothing)
plot!(0:0.01:50, sol[1][1:5001], legend=:topleft, label="OU", linecolor="grey50")
hline!([0], label="", color=:black)
title!("\nTime series")
ylabel!("Measurement")
xlabel!("Time, s")

out_mult = plot(time_mult, power_2nd, size=(1000,400), margin=10mm, plot_title="Ornstein-Uhlenbeck Process with Multiplicative Noise")
savefig(out_mult, "out_mult.png")
=#