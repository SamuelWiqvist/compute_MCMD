# load packages

using PyPlot
using Printf
using LaTeXStrings

################################################################################
# Warmup: simulate the Ising model using the Metroplis algorithm
################################################################################

# generate start condiguration for stat S
dims = 100
S_start = zeros(dims,dims)
map!(x -> x = rand([1,-1]), S_start, S_start)

# the energy function for the ising model
function E_ising(S, J)

    energy = 0

    # sum over pairs of adjacent spins
    for i = 1:size(S,1)-1 # column
        for j = 1:size(S,2)-1 # row

            energy = energy + S[i,j]*S[i+1,j]
            energy = energy + S[i,j]*S[i,j+1]

         end
    end

    # last cloumn
    for i = 1:size(S,1)-1 # column

        j = size(S,2) # row
        energy = energy + S[i,j]*S[i+1,j]

    end


    # last row
    for j = 1:size(S,2)-1 # column

        i = size(S,1) # row
        energy = energy + S[i,j]*S[i,j+1]

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

    E_old = E_ising(S_configs[1,:,:], J) # first iteration
    energy_vec[1] = E_old
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

        E_new = E_ising(S_update,J) # compute energy for prop config

        # compute acc prob
        if E_new > E_old
            α_log = -β*(E_new-E_old)
        else
            α_log = log(1)
        end

        if log(rand()) < α_log # accapt new config
            S_configs[i,:,:] = S_update
            E_old = E_new
            a_vec[i] = 1
            energy_vec[i] = E_new
        else
            S_configs[i,:,:] = S_configs[i-1,:,:] # store old config
            energy_vec[i] = E_old

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
#PyPlot.savefig("hw_week_3/fig/ising_energy_vs_iter.eps", format="eps", dpi=1000)

PyPlot.figure(figsize=(8,12))

PyPlot.subplot(3,2,1)
PyPlot.imshow(S_configs[1,:,:], cmap="hot", interpolation="nearest")
PyPlot.xlabel("Iter. 1")

PyPlot.subplot(3,2,2)
PyPlot.imshow(S_configs[1000,:,:], cmap="hot", interpolation="nearest")
PyPlot.xlabel("Iter. 1000")

PyPlot.subplot(3,2,3)
PyPlot.imshow(S_configs[10000,:,:], cmap="hot", interpolation="nearest")
PyPlot.xlabel("Iter. 10000")

PyPlot.subplot(3,2,4)
PyPlot.imshow(S_configs[100000,:,:], cmap="hot", interpolation="nearest")
PyPlot.xlabel("Iter. 100000")

PyPlot.subplot(3,2,5)
PyPlot.imshow(S_configs[150000,:,:], cmap="hot", interpolation="nearest")
PyPlot.xlabel("Iter. 150000")

PyPlot.subplot(3,2,6)
PyPlot.imshow(S_configs[200000,:,:], cmap="hot", interpolation="nearest")
PyPlot.xlabel("Iter. 200000")

#PyPlot.savefig("hw_week_3/fig/ising_config_conv.png", format="png")
