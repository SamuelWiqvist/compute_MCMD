using PyPlot

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

t_eval = veto(10)

PyPlot.figure()
PyPlot.plot(t_eval)

PyPlot.figure()
h = PyPlot.plt.hist(t_eval,10)
