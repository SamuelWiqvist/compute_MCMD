# load packages

using PyPlot
using Printf
using LaTeXStrings

################################################################################
# Exercise: Computing the density of states for the Potts model in 2D with q = 10
################################################################################

# the energy function for the ising model
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
    idx_update = findfirst(x-> x == E_vec[1], E)
    g_tilde[idx_update] = f*g_tilde[idx_update]

    α_log = log(1)

    for i in 2:iter_max

        # ordinary update
        S_update = S_configs[i-1,:,:] # select site to flip at random

        # update 100 spins

        for j in 1:1
            spins = collect(1:q)
            s_flip = rand(1:nbr_stats) # flip site
            filter!(x->x≠S_update[s_flip],spins) # only select a new spin
            S_update[s_flip] = rand(spins) # set prop config
        end

        E_new = H_potts(S_update,J) # compute energy for prop config

        idx_old_energy = findfirst(x-> x == E_vec[i-1], E)
        g_tilde_old = g_tilde[idx_old_energy]

        idx_new_energy = findfirst(x-> x == E_new, E)
        g_tilde_prop = g_tilde[idx_new_energy]

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
        idx_update = findfirst(x -> x == E_vec[i], E)
        g_tilde[idx_update] = f*g_tilde[idx_update]

    end

    return S_configs, a_vec, E_vec, iter_max

end


# Wang-Landau algorithm where the f value is decreased
function wang_landau(nbr_reps,iter, J, q, f, S_start, g_tilde, E)

    println("Starting Wang-Landau.")

    # full Wang-Landau algorithm
    f_save = zeros(nbr_reps)
    S_configs_last = zeros(iter, size(S_start,1),size(S_start,2))

    for i in 1:nbr_reps

        S_start = ones(size(S_start))

        S_configs, a_vec, E_vec, iter_done = wang_landau_one_iteration(S_start, iter, J, q, f, g_tilde, E)

        f_save[i] = f

        # print info
        @printf "-----------------------------------------------------------------\n"
        @printf "Iter: %.0f\n" i
        @printf "f: %f\n" f
        @printf "MCMC iterations: %.0f\n" iter_done
        @printf "Acc. rate: %.2f %%\n" sum(a_vec)/iter_done*100
        @printf "E_min: %.0f\n" minimum(E_vec)
        @printf "E_max: %.0f\n" maximum(E_vec)

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


# test H_potts
L = 60
H_potts(zeros(L,L),1)
-2*L*L


# run Wang-Landau algorithm

# algorithm settings
q = 10 # spins
J = 1 # interaction strength (we have the same interaction strength for all states)
iter = 10000 # nbr of MC iterations
f = 1.1
nbr_reps = 18 # such that exp(1)^((1/2)^25) \approx exp(10^(-8))

# generate start condiguration for stat S
L = 60
N = L*L
E_min = -2*N
S_start = ones(L,L)
#map!(x -> x = rand(1:q), S_start, S_start)

E = collect(E_min:0)

g_tilde = ones(length(E))


PyPlot.figure()
PyPlot.imshow(S_start,cmap="hot", interpolation="nearest")
PyPlot.colorbar()
#PyPlot.savefig("hw_week_3/fig/potts_start_config.eps", format="eps", dpi=1000)



f_save, S_configs_last = @time wang_landau(nbr_reps,iter, J, q, f, S_start, g_tilde, E)


# init g_tilde and E vectors
E_final = LinRange(0,E_min, 100)
E_upper = E_final[1:end-1]
E_lower = E_final[2:end]

E_matrix = zeros(2,length(E_final)-1)
E_matrix[1,:] = E_upper
E_matrix[2,:] = E_lower

eval_point = zeros(length(E_final)-1)
g_tilde_final = ones(length(eval_point))
g_tilde_final_nbr_events = zeros(length(eval_point))

for i in 1:length(E_final)-1
    eval_point[i] = sum(E_matrix[:,i])/2
end



function findidx(E, E_matrix)

    for i = 1:size(E_matrix,2)

        if E <= round(E_matrix[1,i]) && E >= round(E_matrix[2,i])
            return i
        end
    end

end


for i in E
    println(i)
    idx = findidx(i, E_matrix)
    idx2 =  findfirst(x-> x == i, E)
    g_tilde_final[idx] = g_tilde_final[idx] + g_tilde[idx2]
    g_tilde_final_nbr_events[idx] = g_tilde_final_nbr_events[idx] + 1
end

g_tilde_final
g_tilde_final_nbr_events



g_tilde_final = g_tilde_final./g_tilde_final_nbr_events


# test findidx


minimum(g_tilde)
maximum(g_tilde)


E_vec = H_potts(S_start,J)

# update g_tilde for first iteration
idx_new_energy = findidx(E_vec, E_matrix)

S_start = ones(L,L)

f = 1.1


S_configs, a_vec, E_vec, iter_max = @time wang_landau_one_iteration(S_start, 10000, J, q, f, g_tilde, E)


sum(a_vec)/10000*100

E_vec

PyPlot.figure()
PyPlot.plot(E_vec)

minimum(E_vec)
maximum(E_vec)



minimum(g_tilde)
maximum(g_tilde)


E
minimum(E)


findall(x -> x != 0, g_tilde)



PyPlot.figure()
PyPlot.plot(eval_point, g_tilde_final, "*")
PyPlot.xlabel("Energy")
PyPlot.ylabel(L"\tilde{g}")

PyPlot.figure()
PyPlot.plot(eval_point, log.(g_tilde_final), "*")
PyPlot.xlabel("Energy")
PyPlot.ylabel(L"log(\tilde{g})")



PyPlot.figure()
PyPlot.semilogy(E, g_tilde, "*")
PyPlot.xlabel("Energy")
PyPlot.ylabel(L"\tilde{g}")


PyPlot.figure()
PyPlot.plot(E, log.(g_tilde), "*")
PyPlot.xlabel("Energy")
PyPlot.ylabel(L"log(\tilde{g})")

PyPlot.figure()
PyPlot.plot(E/N, P_T, "*")

#PyPlot.savefig("hw_week_3/fig/potts_g_tilde_final.eps", format="eps", dpi=1000)

PyPlot.figure()
PyPlot.plot(f_save)
PyPlot.xlabel("Iteration")
PyPlot.ylabel(L"f")
#PyPlot.savefig("hw_week_3/fig/potts_f_vs_iter.eps", format="eps", dpi=1000)


T = 100

P_T = exp.(log.(g_tilde_final).-eval_point/T)

PyPlot.figure()
PyPlot.plot(eval_point/N, P_T, "*")

PyPlot.figure()
PyPlot.semilogy(eval_point/N, P_T, "*")


# test energy function

L = 60

H_potts(zeros(L,L),1)

-2*L*L
