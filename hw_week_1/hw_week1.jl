# Warm-up exercises, week 1

using PyPlot
using KernelDensity


# warm-up

u = rand(10^5)

PyPlot.figure()
h = PyPlot.plt[:hist](u,100)
PyPlot.xlabel("x")
PyPlot.ylabel("Freq.")
PyPlot.savefig("hw_week_1/fig/rs_warmup.eps", dpi=150)


# sample from f

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
PyPlot.savefig("hw_week_1/fig/rs_fig.eps", dpi=150)

# integration using hit-and-miss

I_true = 2/pi

N_try = floor.(Int,LinRange(100,10000, 50))

I_hitandmiss = zeros(length(N_try))

f(x) = cos.(pi.*x./2)

for h in 1:length(N_try)
    N_acc = 0
    for i in 1:N_try[h]
        x_star = rand()
        y = rand()
        if y <= f(x_star)
            N_acc = N_acc + 1
        end
    end
    I_hitandmiss[h] = N_acc/N_try[h]
end

# integration using hit-and-miss with importance sampling

I_importance = zeros(length(N_try))

g(x) = 1 .- x.^2

for h in 1:length(N_try)
    N_acc = 0
    samples = zeros(N_try[h])
    for i in 1:N_try[h]
        new_prop = true
        while new_prop
            x_star = rand()
            y = rand()
            if y <= g(x_star)
                new_prop = false
                if y <= f(x_star)
                    N_acc = N_acc + 1
                end
            end
        end
    end
    I_importance[h] = 2/3*N_acc/N_try[h]
end

PyPlot.figure()
PyPlot.plot(N_try, I_hitandmiss, "--*r")
PyPlot.plot(N_try, I_importance, "--*g")
PyPlot.plot(N_try, I_true*ones(length(N_try)), "b")
PyPlot.xlabel("Numbre of samples")
PyPlot.ylabel("Area")
PyPlot.savefig("hw_week_1/fig/rs_area_importance_and_hitandmiss.eps", dpi=150)
