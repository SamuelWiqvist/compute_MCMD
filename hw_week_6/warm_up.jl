using PyPlot
using Statistics

# Diff eq: dN(t)/dt = -f(t)N(t)
# f(t) = 1 + cos(t)
# N(=) = 1

# analytical solution

N(t) = exp(-(t+sin(t)))

t = 0:0.01:4

analytical_sol = zeros(length(t))

for i in 1:length(t); analytical_sol[i] = N(t[i]); end

PyPlot.figure()
PyPlot.plot(analytical_sol)


# veto algorithm

i = 0; t_old = 0; t_end = 4

f(t) = 1 + cos(t)

t_eval = zeros(10)

function veto(steps::Int)

    t_old = 0
    t_eval = zeros(steps)

    for j = 1:steps
        run = true
        while run
            t_new = (2*t_old-log(rand()))/2
            y = rand()*2
            if y > f(t_new)
                # continue
            else
                run = false
                t_old = t_new
                t_eval[j] = t_old
            end
        end
    end

    return t_eval
end

nbr_runs = 100
nbr_steps = 10

event_times = zeros(nbr_steps, nbr_runs)

veto(nbr_steps)

for i = 1:nbr_runs
    event_times[:,i] = veto(nbr_steps)
end


event_times


event_times_mean = mean(event_times, dims = 2)




PyPlot.figure()
PyPlot.plot(event_times_mean)

PyPlot.figure()
h = PyPlot.plt.hist(event_times[3,:],10)
