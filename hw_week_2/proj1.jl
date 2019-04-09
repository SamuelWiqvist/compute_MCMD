using PyPlot
using Statistics

# density function

function f(x)

    if x >= 0.08 && x <= 0.5
        γ = 1.3
    elseif x >= 0.5 && x <= 1
        γ = 2.2
    else
        γ = 2.7
    end

    return x^(-γ)

end

x = LinRange(0.08,120, 200)
density = zeros(length(x))

for i in 1:length(x); density[i] = f(x[i]); end

PyPlot.figure()
PyPlot.plot(x, density)


# hit-and-miss sampler
y_max = round(f(0.08))

function hitandmiss(N_samples)

    N_acc = 0

    samples = zeros(N_samples)

    for i = 1:N_samples

        generate = true

        while generate
            x_star = 0.08+(120-0.08)*rand()
            y_star = y_max*rand()

            if y_star <= f(x_star)
                generate = false
                samples[i] = x_star
            end

        end

    end

    return samples

end


# test

@time hitandmiss(100)

m = hitandmiss(1000)

PyPlot.figure()
h = PyPlot.plt[:hist](m,100)

# Compute prob for a supernova

N_samples = [100,300,1000]
prob_supernova = zeros(length(N_samples))

for i in 1:length(N_samples)
    prob_supernova[i] = length(findall(x -> x > 8, hitandmiss(N_samples[i])))/N_samples[i]
end

print(prob_supernova)

# plot prob supernova
N_samples = floor.(Int,LinRange(50,5000,200))
prob_supernova = zeros(length(N_samples))

for i in 1:length(N_samples)
    prob_supernova[i] = length(findall(x -> x > 8, hitandmiss(N_samples[i])))/N_samples[i]
end

PyPlot.figure()
PyPlot.plot(N_samples, prob_supernova)

# compute statistics

m = hitandmiss(5000)

mean(m)
median(m)
quantile(m, 0.25)

# How many stars in the cluster that we need for the sun

N_samples = floor.(Int,LinRange(100,10^4,100))
prob_supernova = zeros(length(N_samples))

for i in 1:length(N_samples)
    println(i)
    prob_supernova[i] = length(findall(x -> x > 25, hitandmiss(N_samples[i])))/N_samples[i]
end

PyPlot.figure()
PyPlot.plot(N_samples, prob_supernova)
PyPlot.plot(N_samples, zeros(length(N_samples)), "k")
