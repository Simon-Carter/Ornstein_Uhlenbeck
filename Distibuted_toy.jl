# Import libraries.
@everywhere using Turing, 
using StatsPlots, Random

# Set the true probability of heads in a coin.
p_true = 0.5

# Iterate from having seen 0 observations to 100 observations.
Ns = 0:100

# Draw data from a Bernoulli distribution, i.e. draw heads or tails.
Random.seed!(12)
data = rand(Bernoulli(p_true), last(Ns))

# Declare our Turing model.
@everywhere @model function coinflip(y)
    # Our prior belief about the probability of heads in a coin.
    p ~ Beta(1, 1)

    # The number of observations.
    N = length(y)
    for n in 1:N
        # Heads or tails of a coin are drawn from a Bernoulli distribution.
        y[n] ~ Bernoulli(p)
    end
end

# Settings of the Hamiltonian Monte Carlo (HMC) sampler.
iterations = 1000
ϵ = 0.05
τ = 10

#setup multi-processing
const jobs = RemoteChannel(()->Channel{Int}(32));
const results = RemoteChannel(()->Channel{Chains}(32));

@everywhere function do_work(jobs, results, data, ϵ, τ, iterations) # define work function everywhere
    while true
        job_id = take!(jobs)
        chain = sample(coinflip(data), HMC(ϵ, τ), iterations)
        put!(results, chain)
    end
end

function make_jobs(n)
    for i in 1:n
        put!(jobs, i)
    end
end;

n = 5;

make_jobs(n)

for p in workers() # start tasks on the workers to process requests in parallel
    remote_do(do_work, p, jobs, results, data, ϵ, τ, iterations)
end

total_data = Vector{Chains}(undef, 5)

@elapsed while n > 0 # print out results
    result = take!(results)
    global total_data[n] =  result
    global n = n - 1
end