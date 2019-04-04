# Warm-up exercises, week 1

using PyPlot
using KernelDensity


# warm-up

u = rand(10^5)

PyPlot.figure()
h = PyPlot.plt[:hist](u,100)
PyPlot.xlabel("x")
PyPlot.ylabel("Density")
PyPlot.savefig("fig/hw1/rs_warmup.eps", dpi=150)


# sample for f

N = 10^5
M = 1/100
samples = zeros(N)

f(x) = sin.(x).*sin.(x)./x.^2

for i = 1:N
    new_prop = true
    while new_prop
        x_star = 100*rand()
        u = rand()
        if u <= f(x_star)/(M*1/M)
            samples[i] = x_star
            new_prop = false
        end
    end
end

kde_approx = kde(samples,  boundary = (0,100))

PyPlot.figure()
PyPlot.plot(kde_approx.x[1:end-3],kde_approx.density[1:end-3], "--r")
PyPlot.plot(collect(kde_approx.x[1:end-3]), f(kde_approx.x[1:end-3]), "b")
PyPlot.xlabel("x")
PyPlot.ylabel("Density")
PyPlot.savefig("fig/hw1/rs_fig.eps", dpi=150)

# integration using hit-and-miss

I_true = 2/pi

N_try = floor.(Int,LinRange(100,100000, 50))

I_hitandmiss = zeros(length(N_try))

f(x) = cos(pi*x/2)

for h in 1:length(N_try)
    N_acc = 0
    for i in 1:N_try[h]
        x_star = rand()
        u = rand()
        if u <= f(x_star)
            N_acc = N_acc + 1
        end
    end
    I_hitandmiss[h] = N_acc/N_try[h]
end



PyPlot.figure()
PyPlot.plot(N_try, I_hitandmiss, "--*r")
PyPlot.plot(N_try, I_true*ones(length(N_try)), "b")



# integration using hit-and-miss with importance sampling

I_true = 2/pi

N_try = floor.(Int,LinRange(100,100000, 50))

I_hitandmiss = zeros(length(N_try))



f(x) = cos.(pi.*x./2)
g(x) = 1 .- x.^2


x = collect(LinRange(0,1, 50))

PyPlot.figure()
PyPlot.plot(x,f(x), "b")
PyPlot.plot(x,g(x), "r")

for h in 1:length(N_try)
    N_acc = 0
    for i in 1:N_try[h]

        new_prop = true

        while new_prop
            x_star = rand()
            u = rand()
            if u <= g(x_star)
                new_prop = false
            end
        end

        u = rand()
        if u <= f(x_star)
            N_acc = N_acc + 1
        end

    end

    I_hitandmiss[h] = N_acc/N_try[h]
end


PyPlot.figure()
PyPlot.plot(N_try, I_hitandmiss, "--*r")
PyPlot.plot(N_try, I_true*ones(length(N_try)), "b")
