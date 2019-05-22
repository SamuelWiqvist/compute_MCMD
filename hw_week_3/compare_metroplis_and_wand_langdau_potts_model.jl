# load packages

using PyPlot
using Printf
using LaTeXStrings

# generate start condiguration for stat S
L = 10
N = L*L
E_min = -2*N

q = 10 # spins
J = 1 # interaction strength (we have the same interaction strength for all states)

S_start = ones(L,L)

map!(x -> x = rand(1:q), S_start, S_start)

PyPlot.figure()
PyPlot.imshow(S_start,cmap="hot", interpolation="nearest")
PyPlot.colorbar()

# the energy function for the Potts model
# code adaped from https://github.com/maxjkiss/q-state-potts-model
function H_potts(S, J)

    energy = 0
    L = size(S,1)

    for i in 0:L-1
        for j in 0:L-1

            energy = energy + δ(S[i+1,j+1], S[mod(i+1,L)+1,j+1])
            energy = energy + δ(S[i+1,j+1], S[i+1,mod(j+1,L)+1])

         end
    end

    return -J*energy

end

# help function for H_potts
δ(s, s_star) = s == s_star ? 1 : 0 # kronecker delta


################################################################################
# Metroplis algorithm
################################################################################



# metroplis algorithm to simulate the ising model
function metroplis(S_start, iter, J, β)

    @printf "Starting Metoplis algorithm\n"

    nbr_stats = length(S_start) # set up
    #S_configs = zeros(iter, size(S_start,1),size(S_start,2))
    a_vec = zeros(iter)
    energy_vec = zeros(iter)
    #S_configs[1,:,:] = S_start

    S_old = deepcopy(S_start) # first iteration


    E_old = H_potts(S_old, J) # first iteration
    energy_vec[1] = E_old
    α_log = log(1)

    for i in 2:iter

        # print info
        if mod(i,100000) == 0
            @printf "-----------------------------------------------------------------\n"
            @printf "Iter: %.0f\n" i
            @printf "Percentage done: %.2f %%\n" i/iter*100
            @printf "Current energy: %.2f\n" energy_vec[i-1]
        end

        # ordinary update
        S_update = deepcopy(S_old) # select site to flip at random

        s_flip = rand(1:nbr_stats) # flip ste
        S_update[s_flip] = rand(1:q) # set prop config

        E_new = H_potts(S_update,J) # compute energy for prop config

        # compute acc prob
        if E_new > E_old
            α_log = -β*(E_new-E_old)
        else
            α_log = log(1)
        end


        if log(rand()) < α_log # accapt new config
            #S_configs[i,:,:] = S_update
            S_old = deepcopy(S_update)
            E_old = E_new
            a_vec[i] = 1
            energy_vec[i] = E_new
        else
            #S_configs[i,:,:] = S_configs[i-1,:,:] # store old config
            energy_vec[i] = E_old

        end

    end

    return energy_vec,a_vec

end


# set system parameters
T = 0.7145 # run metroplis at T_c
β = 1/T # tempering
iter = 10^7 # nbr of MC interstions

# run Metroplis algorithm
energy_vec, a_vec = @time metroplis(S_start, iter, J, β)

# compute avg acc prob
@printf "Avg. acc. prob.: %.2f %%\n" sum(a_vec)/iter*100

# plotting
PyPlot.figure(figsize=(8,5))
PyPlot.plot(1:iter, energy_vec)
PyPlot.xlabel("Iteration")
PyPlot.ylabel("Energy")
#PyPlot.savefig("hw_week_3/fig/ising_energy_vs_iter.eps", format="eps", dpi=1000)

PyPlot.figure(figsize=(8,12))

PyPlot.subplot(3,2,1)
PyPlot.imshow(S_configs[1,:,:], cmap="hot", interpolation="nearest")
PyPlot.xlabel("Iteration 1")

PyPlot.subplot(3,2,2)
PyPlot.imshow(S_configs[1000,:,:], cmap="hot", interpolation="nearest")
PyPlot.xlabel("Iteration 1000")

PyPlot.subplot(3,2,3)
PyPlot.imshow(S_configs[10000,:,:], cmap="hot", interpolation="nearest")
PyPlot.xlabel("Iteration 10000")

PyPlot.subplot(3,2,4)
PyPlot.imshow(S_configs[100000,:,:], cmap="hot", interpolation="nearest")
PyPlot.xlabel("Iteration 100000")

PyPlot.subplot(3,2,5)
PyPlot.imshow(S_configs[150000,:,:], cmap="hot", interpolation="nearest")
PyPlot.xlabel("Iteration 150000")

PyPlot.subplot(3,2,6)
PyPlot.imshow(S_configs[200000,:,:], cmap="hot", interpolation="nearest")
PyPlot.xlabel("Iteration 200000")

#PyPlot.savefig("hw_week_3/fig/ising_config_conv.png", format="png")
