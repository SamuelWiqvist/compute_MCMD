using PyPlot
using Statistics
using KernelDensity
using Random
# Diff eq: dN(t)/dt = -f(t)N(t)
# f(t) = 1 + cos(t)
# N(0) = 1

# analytical solution

N(t) = exp(-t-sin(t))

t = 0:0.01:10

analytical_sol = zeros(length(t))

for i in 1:length(t); analytical_sol[i] = N(t[i]); end

PyPlot.figure()
PyPlot.plot(t, analytical_sol)
PyPlot.xlabel("t")
PyPlot.ylabel("N(t)")
PyPlot.savefig("analytical_sol.pdf")

# veto algorithm

f(t) = 1 + cos(t)

function veto(samples::Int)

    t_event = zeros(samples)

    for j = 1:samples
        run = true
        t_old = 0
        while run
            t_new = (2*t_old-log(rand()))/2
            y = rand()*2
            if y > f(t_new)
                t_old = t_new
            else
                run = false
                t_old = t_new
                t_event[j] = t_old
            end
        end
    end

    return t_event
end

Random.seed!(100)
event_times = veto(10^6)

PyPlot.figure()
h = PyPlot.plt[:hist](event_times,100,normed=1)
PyPlot.xlabel("t")
PyPlot.savefig("veto_event_times.pdf")

PyPlot.plot(t, analytical_sol)
