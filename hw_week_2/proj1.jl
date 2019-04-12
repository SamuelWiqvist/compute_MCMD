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

function hitandmiss(N_stars)

    N_acc = 0

    stars = zeros(N_stars)

    for i = 1:N_stars

        generate = true

        while generate
            x_star = 0.08+(120-0.08)*rand()
            y_star = y_max*rand()

            if y_star <= f(x_star)
                generate = false
                stars[i] = x_star
            end

        end

    end

    return stars

end


# test

@time hitandmiss(100)

m = hitandmiss(1000)

PyPlot.figure()
h = PyPlot.plt[:hist](m,100)

# Compute prob for a supernova

N_stars = [100,300,1000]
prob_supernova = zeros(length(N_stars))

for i in 1:length(N_stars)
    prob_supernova[i] = length(findall(x -> x > 8, hitandmiss(N_stars[i])))/N_stars[i]
end

print(prob_supernova)

# plot prob supernova
N_stars = floor.(Int,LinRange(50,5000,200))
prob_supernova = zeros(length(N_stars))

for i in 1:length(N_stars)
    prob_supernova[i] = length(findall(x -> x > 8, hitandmiss(N_stars[i])))/N_stars[i]
end

PyPlot.figure()
PyPlot.plot(N_stars, prob_supernova)

# compute statistics

N_clusters = 100
N_stars = 5000
prob_supernova = zeros(N_clusters)

for i in 1:N_clusters
    prob_supernova[i] = length(findall(x -> x > 8, hitandmiss(N_stars)))/N_stars
end

print(mean(prob_supernova))
print(median(prob_supernova))
print(quantile(prob_supernova, 0.25))

# How many stars in the cluster that we need for the sun

N_stars = floor.(Int,LinRange(100,10^4,100))
prob_supernova = zeros(length(N_stars))

for i in 1:length(N_stars)
    println(i)
    prob_supernova[i] = length(findall(x -> x > 25, hitandmiss(N_stars[i])))/N_stars[i]
end

PyPlot.figure()
PyPlot.plot(N_stars, prob_supernova)
PyPlot.plot(N_stars, zeros(length(N_stars)), "k")
