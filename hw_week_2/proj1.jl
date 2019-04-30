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
PyPlot.loglog(x, density)


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

N_clusters = 100
N_stars = [100,300,1000]
nbr_supernovae = zeros(N_clusters, length(N_stars))

# generate clusters and find number of supernovae in each cluster
for i in 1:length(N_stars)
    for j in 1:N_clusters
        nbr_supernovae[j,i] = length(findall(x -> x > 8, hitandmiss(N_stars[i])))
    end
end

# compute probabiltiy that a cluster contains at least one supernova
prob_supernovae = zeros(length(N_stars))

for i in 1:length(N_stars)
    prob_supernovae[i] =  length(findall(x -> x > 0, nbr_supernovae[:,i]))/N_clusters
end

print(prob_supernovae)

# plot prob supernova
N_stars = floor.(Int,LinRange(50,5000,50))
nbr_supernovae = zeros(N_clusters, length(N_stars))

# generate clusters and find number of supernovae in each cluster
for i in 1:length(N_stars)
    for j in 1:N_clusters
        nbr_supernovae[j,i] = length(findall(x -> x > 8, hitandmiss(N_stars[i])))
    end
end


# compute probabiltiy that a cluster contains at least one supernova
prob_supernovae = zeros(length(N_stars))

for i in 1:length(N_stars)
    prob_supernovae[i] =  length(findall(x -> x > 0, nbr_supernovae[:,i]))/N_clusters
end


PyPlot.figure()
PyPlot.plot(N_stars, prob_supernovae)

# compute statistics

N_clusters = 100
N_stars = 5000
nbr_supernovae = zeros(N_clusters)

for i in 1:N_clusters
    nbr_supernovae[i] = length(findall(x -> x > 8, hitandmiss(N_stars)))
end

print(mean(nbr_supernovae))
print(median(nbr_supernovae))
print(quantile(nbr_supernovae, 0.25))

# How many stars in the cluster that we need for the sun

N_clusters = 1000
nbr_stars_in_cluster = zeros(N_clusters)

for i in 1:N_clusters # loop over the number of clusters

    nbr_stars = 0
    generate_starts = true

    # generate stars in cluster i until we obtain one start  with mass > 25
    while generate_starts

        nbr_stars = nbr_stars + 1
        star_prop = hitandmiss(1)[1]

        if star_prop >= 25
            nbr_stars_in_cluster[i] = nbr_stars
            generate_starts = false
        end
    end

end

PyPlot.figure()
h = PyPlot.plt[:hist](nbr_stars_in_cluster,100)

println(mean(nbr_stars_in_cluster))
println(std(nbr_stars_in_cluster))
