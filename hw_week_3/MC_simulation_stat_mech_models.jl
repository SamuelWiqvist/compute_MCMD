# load packages

using PyPlot
using Printf
using LaTeXStrings

# TODO
# 1) fix sum where the first sum is over pairs of adjacent spins (every pair is counted once)
# 2) set parameters to resnoable avleus
# 3) rerun

# TODO
# 1) why do the results differ so much for different runs?

################################################################################
# Warmup: simulate the Ising model using the Metroplis algorithm
################################################################################

# generate start condiguration for stat S
dims = 100
S_start = zeros(dims,dims)
map!(x -> x = rand([1,-1]), S_start, S_start)

PyPlot.figure()
PyPlot.imshow(S_start, cmap="hot", interpolation="nearest")
PyPlot.savefig("hw_week_3/fig/ising_start_config.eps", format="eps", dpi=1000)

# the energy function for the ising model
function H_ising(S, J)

    energy = 0

    for i = 1:size(S,1)-1 # column
        for j = 1:size(S,2)-1 # row

            energy = energy + S[i,j]*S[i+1,j]
            energy = energy + S[i,j]*S[i,j+1]

         end
    end

    return -J*energy

end

# metroplis algorithm to simulate the ising model
function metroplis(S_start, iter, J, β)

    @printf "Starting Metoplis algorithm\n"

    nbr_stats = length(S_start) # set up
    S_configs = zeros(iter, size(S_start,1),size(S_start,2))
    a_vec = zeros(iter)
    energy_vec = zeros(iter)
    S_configs[1,:,:] = S_start

    H_old = H_ising(S_configs[1,:,:], J) # first iteration
    energy_vec[1] = H_old
    α_log = log(1)

    for i = 2:iter

        # print info
        if mod(i,10000) == 0
            @printf "-----------------------------------------------------------------\n"
            @printf "Iter: %0.f\n" i
            @printf "Percentage done: %.2f %%\n" i/iter*100
        end

        S_update = S_configs[i-1,:,:] # select state to flip at random
        S_flip = rand(1:nbr_stats) # flip state
        S_update[S_flip] = -1*S_update[S_flip] # set prop config

        H_new = H_ising(S_update,J) # compute energy for prop config

        # compute acc prob
        if H_new > H_old
            α_log = -β*(H_new-H_old)
        else
            α_log = log(1)
        end

        if log(rand()) < α_log # accapt new config
            S_configs[i,:,:] = S_update
            H_old = H_new
            a_vec[i] = 1
            energy_vec[i] = H_new
        else
            S_configs[i,:,:] = S_configs[i-1,:,:] # store old config
            energy_vec[i] = H_old

        end

    end

    return S_configs,energy_vec,a_vec

end

# set system parameters
J = 0.5 # interaction strength (we have the same interaction strength for all states)
β = 20 # tempering
iter = 200000 # nbr of MC interstions

# run Metroplis algorithm
S_configs, energy_vec, a_vec = @time metroplis(S_start, iter, J, β)

# compute avg acc prob
@printf "Avg. acc. prob.: %.2f %%\n" sum(a_vec)/iter*100

# plotting
PyPlot.figure(figsize=(8,5))
PyPlot.plot(1:iter, energy_vec)
PyPlot.xlabel("Iteration")
PyPlot.ylabel("Energy")
PyPlot.savefig("hw_week_3/fig/ising_energy_vs_iter.eps", format="eps", dpi=1000)

iter_plot = floor.(Int,LinRange(1,iter,100))

PyPlot.figure()
for i = iter_plot
    sleep(0.001)
    PyPlot.imshow(S_configs[i,:,:], cmap="hot", interpolation="nearest")
    PyPlot.xlabel(i)
end

PyPlot.figure(figsize=(8,12))

PyPlot.subplot(3,2,1)
PyPlot.imshow(S_configs[1,:,:], cmap="hot", interpolation="nearest")
PyPlot.xlabel(1)

PyPlot.subplot(3,2,2)
PyPlot.imshow(S_configs[1000,:,:], cmap="hot", interpolation="nearest")
PyPlot.xlabel(1000)

PyPlot.subplot(3,2,3)
PyPlot.imshow(S_configs[10000,:,:], cmap="hot", interpolation="nearest")
PyPlot.xlabel(10000)

PyPlot.subplot(3,2,4)
PyPlot.imshow(S_configs[100000,:,:], cmap="hot", interpolation="nearest")
PyPlot.xlabel(100000)

PyPlot.subplot(3,2,5)
PyPlot.imshow(S_configs[150000,:,:], cmap="hot", interpolation="nearest")
PyPlot.xlabel(150000)

PyPlot.subplot(3,2,6)
PyPlot.imshow(S_configs[200000,:,:], cmap="hot", interpolation="nearest")
PyPlot.xlabel(200000)
PyPlot.savefig("hw_week_3/fig/ising_config_conv.eps", format="eps", dpi=1000)

################################################################################
# Exercise: Computing the density of states for the Potts model in 2D with q = 10
################################################################################

# the energy function for the ising model
function H_potts(S, J)

    energy = 0

    for i = 1:size(S,1)-1 # column
        for j = 1:size(S,2)-1 # row
            
            energy = energy + δ(S[i,j], S[i+1,j])
            energy = energy + δ(S[i,j], S[i,j+1])

         end
    end

    return -J*energy

end

# help function for H_potts
function δ(s, s_star)

    if s == s_star
        return 1
    else
        return 0
    end

end

# Wang-Landau algorithm with one fixed f value
function wang_landau_one_iteration(S_start, iter_max, J, q, f, g_tilde, E)

    nbr_stats = length(S_start) # set up
    S_configs = zeros(iter_max, size(S_start,1),size(S_start,2))
    a_vec = zeros(iter_max)
    E_vec = zeros(iter_max)

    S_configs[1,:,:] = S_start # first iteration

    a_vec[1] = 1
    E_vec[1] = H_potts(S_start,J)

    # update g_tilde for first iteration
    g_update_idx = findlast(x -> x == E_vec[1], E)
    if typeof(g_update_idx) == Nothing
        g_update = 1
        append!(E,[E_vec[1]])
        append!(g_tilde,[f*g_update])
    else
        g_update = g_tilde[g_update_idx[1]]
        g_tilde[g_update_idx[1]] = f*g_update
    end


    α_log = log(1)

    min_energy = 0
    max_energy = 0

    for i = 2:iter_max

        # ordinary update
        S_update = S_configs[i-1,:,:] # select state to flip at random
        S_flip = rand(1:nbr_stats) # flip state
        S_update[S_flip] = rand(1:q) # set prop config

        H_new = H_potts(S_update,J) # compute energy for prop config

        # see if we have curent energy
        g_tilde_prop_idx = findlast(x -> x == H_new, E)

        if typeof(g_tilde_prop_idx) == Nothing
            g_tilde_prop = 1
        else
            g_tilde_prop = g_tilde[g_tilde_prop_idx[1]]
        end

        g_tilde_old_idx = findlast(x -> x == E_vec[i-1], E)
        g_tilde_old = g_tilde[g_tilde_old_idx[1]]

        α_log = log(g_tilde_old)-log(g_tilde_prop)

        if log(rand()) < min(0, α_log) # accapt new config
            S_configs[i,:,:] = S_update
            E_vec[i] = H_new
            a_vec[i] = 1
        else
            S_configs[i,:,:] = S_configs[i-1,:,:] # store old config
            E_vec[i] = E_vec[i-1]
        end

        # update energy function for current system
        g_update_idx = findlast(x -> x == E_vec[i], E)

        if typeof(g_update_idx) == Nothing
            g_update = 1
            append!(E,[E_vec[i]])
            append!(g_tilde,[f*g_update])
        else
            g_update = g_tilde[g_update_idx[1]]
            g_tilde[g_update_idx[1]] = f*g_update
        end

        if i == 1000

            min_energy = minimum(E_vec[1:i])
            max_energy = maximum(E_vec[1:i])


        end

        if mod(i,100) == 0 && i > 1000


            nbr_min = length(findall(x -> x <= min_energy, E_vec[1:i]))
            nbr_max = length(findall(x -> x >= max_energy, E_vec[1:i]))

            if nbr_min > 5 && nbr_max > 5
                return S_configs[1:i,:,:], a_vec[1:i], E_vec[1:i], i
            end


        end


    end

    return S_configs, a_vec, E_vec, iter_max

end


# Wang-Landau algorithm where the f value is decreased
function wang_landau(nbr_reps,iter, J, q, f, S_start, g_tilde, E)

    println("Starting Wang-Landau.")

    # full Wang-Landau algorithm
    f_save = zeros(nbr_reps)
    S_configs_last = zeros(iter, size(S_start,1),size(S_start,2))

    for i = 1:nbr_reps

        S_configs, a_vec, E_vec, iter_done = wang_landau_one_iteration(S_start, iter, J, q, f, g_tilde, E)

        f_save[i] = f

        # print info
        @printf "-----------------------------------------------------------------\n"
        @printf "Iter: %.0f\n" i
        @printf "f: %f\n" f
        @printf "MCMC iterations: %.0f\n" iter_done
        @printf "Acc. rate: %.2f %%\n" sum(a_vec)/iter_done*100

        #=
        PyPlot.figure()
        PyPlot.plot(E, g_tilde, "*-")
        PyPlot.xlabel("Energy")
        PyPlot.ylabel(L"\tilde{g}")
        PyPlot.savefig("hw_week_3/fig/potts_g_tilde_iteration_"*string(i)*".eps", format="eps", dpi=1000)
        =#

        # update t
        f = sqrt(f)

        #S_start = S_configs[end,:,:]
        map!(x -> x = rand(1:q), S_start, S_start)

        if i == nbr_reps
            S_configs_last = S_configs
        end

    end

    return f_save, S_configs_last

end


# run Wang-Landau algorithm

# algorithm settings
J = 0.5 # interaction strength (we have the same interaction strength for all states)
iter = 50000 # nbr of MC iterations
f = exp(1)
nbr_reps = 25 # such that exp(1)^((1/2)^25) \approx exp(10^(-8))
q = 5 # spins


# generate start condiguration for stat S
dims = 100
S_start = zeros(dims,dims)
map!(x -> x = rand(1:q), S_start, S_start)

# init g_tilde and E vectors
g_tilde = [1.]
E = [H_potts(S_start,J)]

PyPlot.figure()
PyPlot.imshow(S_start,cmap="hot", interpolation="nearest")
PyPlot.colorbar()
PyPlot.savefig("hw_week_3/fig/potts_start_config.eps", format="eps", dpi=1000)



f_save, S_configs_last = @time wang_landau(nbr_reps,iter, J, q, f, S_start, g_tilde, E)



PyPlot.figure()
PyPlot.plot(E, g_tilde, "*")
PyPlot.xlabel("Energy")
PyPlot.ylabel(L"\tilde{g}")
PyPlot.savefig("hw_week_3/fig/potts_g_tilde_final.eps", format="eps", dpi=1000)




PyPlot.figure()
PyPlot.plot(f_save)
PyPlot.xlabel("Iteration")
PyPlot.ylabel(L"f")
PyPlot.savefig("hw_week_3/fig/potts_f_vs_iter.eps", format="eps", dpi=1000)



T = 1

P_T = exp.(log.(g_tilde).*-E/T)

exp.(log.(g_tilde))

exp.(-E/T)



E


PyPlot.figure()
PyPlot.plot(E, P_T, "*")
PyPlot.xlabel("Energy")
PyPlot.ylabel(L"\tilde{g}")
PyPlot.savefig("hw_week_3/fig/potts_g_tilde_final.eps", format="eps", dpi=1000)



g_tilde

E


S_configs, a_vec, E_vec = wang_landau_one_iteration(S_start, iter, J, q, f, g_tilde, E)


PyPlot.figure()
nbr_in_bins, bins, plotobj = PyPlot.plt[:hist](E_vec,100)

diff_min_max_in_min = maximum(nbr_in_bins) - minimum(nbr_in_bins)



PyPlot.figure()
PyPlot.plot(E, g_tilde, "*")
PyPlot.xlabel("Energy")
PyPlot.ylabel(L"\tilde{g}")
PyPlot.savefig("hw_week_3/fig/potts_g_tilde_final.eps", format="eps", dpi=1000)




PyPlot.figure()
PyPlot.plot(f_save)
PyPlot.xlabel("Iteration")
PyPlot.ylabel(L"f")
PyPlot.savefig("hw_week_3/fig/potts_f_vs_iter.eps", format="eps", dpi=1000)


PyPlot.title("Scaling")


PyPlot.figure()

for i = 1:size(g_save,1)

    g_tilde = g_save[i,:,:]

    # plot g function
    energy = unique(g_tilde[2,:])
    g = zeros(length(energy))


    for i = 1:length(energy)

        idx = findlast(x -> x == energy[i], g_tilde[2,:])
        g[i] = g_tilde[1,idx]

    end


    g = g/sum(g)

    PyPlot.subplot(13,2,i)
    PyPlot.plot(energy, g, "*-")
    PyPlot.title(f_save[i])

end




for i = size(g_save,1)-4:size(g_save,1)

    g_tilde = g_save[i,:,:]

    # plot g function
    energy = unique(g_tilde[2,:])
    g = zeros(length(energy))


    for i = 1:length(energy)

        idx = findlast(x -> x == energy[i], g_tilde[2,:])
        g[i] = g_tilde[1,idx]

    end

    g = g/sum(g)


    PyPlot.figure()
    PyPlot.plot(energy, g, "*-")
    PyPlot.title(f_save[i])


end


iter_plot = floor.(Int,LinRange(1,iter,100))

PyPlot.figure()
for i = iter_plot
    sleep(0.001)
    PyPlot.imshow(S_configs_last[i,:,:], cmap="hot", interpolation="nearest")
    PyPlot.xlabel(i)
end



PyPlot.figure()
PyPlot.imshow(S_configs_last[1,:,:], cmap="hot", interpolation="nearest")

PyPlot.figure()
PyPlot.imshow(S_configs_last[1000,:,:], cmap="hot", interpolation="nearest")

PyPlot.figure()
PyPlot.imshow(S_configs_last[2000,:,:], cmap="hot", interpolation="nearest")

PyPlot.figure()
PyPlot.imshow(S_configs_last[5000,:,:], cmap="hot", interpolation="nearest")

PyPlot.figure()
PyPlot.imshow(S_configs_last[10000,:,:], cmap="hot", interpolation="nearest")


# code to test partial wang-landau
# run algorithm
#S_start = zeros(dims,dims)
#S_start = map!(x -> x = rand(1:q), S_start, S_start)
#S_configs, g_tilde, a_vec = @time wang_landau_one_iteration(S_start, iter, J, q, t)
