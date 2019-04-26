# load packages

using PyPlot
using Printf
using LaTeXStrings

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


    # last cloumn
    for i = 1:size(S,1)-1 # column

        j = size(S,2) # row
        energy = energy + δ(S[i,j], S[i+1,j])

    end


    # last row
    for j = 1:size(S,2)-1 # column

        i = size(S,1) # row
        energy = energy + δ(S[i,j], S[i,j+1])

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

        E_new = H_potts(S_update,J) # compute energy for prop config

        # see if we have curent energy
        g_tilde_prop_idx = findlast(x -> x == E_new, E)

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
            E_vec[i] = E_new
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

        # update f
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
q = 10 # spins
J = -(q/2)*0.5 # interaction strength (we have the same interaction strength for all states)
iter = 50000 # nbr of MC iterations
f = exp(1)
nbr_reps = 25 # such that exp(1)^((1/2)^25) \approx exp(10^(-8))



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
#PyPlot.savefig("hw_week_3/fig/potts_start_config.eps", format="eps", dpi=1000)



f_save, S_configs_last = @time wang_landau(nbr_reps,iter, J, q, f, S_start, g_tilde, E)

PyPlot.figure()
PyPlot.plot(E, g_tilde, "*")
PyPlot.xlabel("Energy")
PyPlot.ylabel(L"\tilde{g}")
#PyPlot.savefig("hw_week_3/fig/potts_g_tilde_final.eps", format="eps", dpi=1000)

PyPlot.figure()
PyPlot.plot(f_save)
PyPlot.xlabel("Iteration")
PyPlot.ylabel(L"f")
#PyPlot.savefig("hw_week_3/fig/potts_f_vs_iter.eps", format="eps", dpi=1000)


T = 1000

P_T = exp.(log.(g_tilde).-E/T)

PyPlot.figure()
PyPlot.plot(E, P_T, "*")
