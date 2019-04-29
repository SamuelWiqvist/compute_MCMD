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
function wang_landau_one_iteration(S_start, iter_max, J, q, f, log_g_tilde, E_matrix)

    nbr_stats = length(S_start) # set up
    #S_configs = zeros(iter_max, size(S_start,1),size(S_start,2))
    a_vec = zeros(iter_max)
    E_vec = zeros(iter_max)

    S_old = deepcopy(S_start) # first iteration

    a_vec[1] = 1
    E_vec[1] = H_potts(S_start,J)
    E_visit_counter = zeros(length(log_g_tilde))

    # update log_g_tilde for first iteration
    idx_update = findidx(E_vec[1], E_matrix)
    log_g_tilde[idx_update] = log(f)+log_g_tilde[idx_update]
    E_visit_counter[idx_update] = 1

    α_log = log(1)

    for i in 2:iter_max

        # ordinary update
        S_update = deepcopy(S_old) # select site to flip at random

        # update 1 spin
        for j in 1:1
            s_flip = rand(1:nbr_stats) # flip site
            S_update[s_flip] = rand(1:q) # set prop config
        end


        idx_old_energy = findidx(E_vec[i-1], E_matrix) # compute log_g for old config
        log_g_tilde_old = log_g_tilde[idx_old_energy]

        E_new = H_potts(S_update,J) # compute g_tilde for prop config
        idx_new_energy = findidx(E_new, E_matrix)
        log_g_tilde_prop = log_g_tilde[idx_new_energy]

        α_log = log_g_tilde_old-log_g_tilde_prop

        if log(rand()) < min(0, α_log) # accapt new config
            S_old = deepcopy(S_update)
            E_vec[i] = E_new
            a_vec[i] = 1
        else
            E_vec[i] = E_vec[i-1]
        end

        # update log_g_tilde for current system
        idx_update = findidx(E_vec[i], E_matrix)
        log_g_tilde[idx_update] = log(f)+log_g_tilde[idx_update]
        E_visit_counter[idx_update] = E_visit_counter[idx_update] +1

    end

    return a_vec, E_vec, E_visit_counter, iter_max

end


# Wang-Landau algorithm where the f value is decreased
function wang_landau(nbr_reps,iter, J, q, f, S_start, log_g_tilde, E_matrix)

    println("Starting Wang-Landau.")

    # full Wang-Landau algorithm
    f_save = zeros(nbr_reps)

    for i in 1:nbr_reps

        map!(x -> x = rand(1:q), S_start, S_start)

        a_vec, E_vec, E_visit_counter, iter_done = wang_landau_one_iteration(S_start, iter, J, q, f, log_g_tilde, E_matrix)

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

    end

    return f_save

end


function findidx(E, E_matrix)

    for i = 1:size(E_matrix,2)

        if E <= round(E_matrix[1,i]) && E >= round(E_matrix[2,i])
            return i
        end
    end

end


# test H_potts
L = 10
H_potts(zeros(L,L),1)
-2*L*L


# run Wang-Landau algorithm

# algorithm settings
q = 10 # spins
J = 1 # interaction strength (we have the same interaction strength for all states)
iter = 10000000 # nbr of MC iterations
f = 2.7
nbr_reps = 25 # such that exp(1)^((1/2)^25) \approx exp(10^(-8))

# generate start condiguration for stat S
L = 10
N = L*L
E_min = -2*N
S_start = ones(L,L)
#map!(x -> x = rand(1:q), S_start, S_start)

# init g_tilde and E vectors
E = LinRange(0,E_min, 50)
E_upper = E[1:end-1]
E_lower = E[2:end]

E_matrix = zeros(2,length(E)-1)
E_matrix[1,:] = E_upper
E_matrix[2,:] = E_lower

eval_point = zeros(length(E)-1)
log_g_tilde = log.(ones(length(eval_point)))

for i in 1:length(E)-1
    eval_point[i] = sum(E_matrix[:,i])/2
end


# set start configuration
map!(x -> x = rand(1:q), S_start, S_start)

PyPlot.figure()
PyPlot.imshow(S_start,cmap="hot", interpolation="nearest")
PyPlot.colorbar()

f_save = @time wang_landau(nbr_reps,iter, J, q, f, S_start, log_g_tilde, E_matrix)

PyPlot.figure()
PyPlot.plot(f_save)
PyPlot.xlabel("Iteration")
PyPlot.ylabel(L"f")

PyPlot.figure()
PyPlot.plot(eval_point, log_g_tilde)
PyPlot.xlabel("Energy")
PyPlot.ylabel(L"log \tilde{g}")


# nomrmalize g
idx_norm = findidx(-2*N, E_matrix)

log_g_tilde_normalized = log_g_tilde .- log_g_tilde[idx_norm] .+ log(q)


PyPlot.figure()
PyPlot.plot(eval_point, log_g_tilde_normalized)
PyPlot.xlabel("Energy")
PyPlot.ylabel(L"log \tilde{g}")


T = 0.7145

# T_critical 0.7145
# 0.7 rigth mode
# 0.8 left mode

P_T = exp.(log_g_tilde_normalized.-eval_point/T)

maximum(P_T)

P_T_normalized = P_T ./ maximum(P_T)

maximum(P_T_normalized)


PyPlot.figure()
PyPlot.plot(eval_point/N, P_T_normalized)
PyPlot.ylabel(L"P_T(E)")
PyPlot.xlabel(L"E/N")

PyPlot.figure()
PyPlot.plot(eval_point,  P_T_normalized)
PyPlot.ylabel(L"P_T(E)")
PyPlot.xlabel(L"Energy")


# wang_landau one iteration

map!(x -> x = rand(1:q), S_start, S_start)

# only run first iteration of wang-landau
f = 2.7
log(f)
log_g_tilde = log.(ones(length(eval_point)))

a_vec, E_vec, E_visit_counter, iter_max = @time wang_landau_one_iteration(S_start, 10000000, J, q, f, log_g_tilde, E_matrix)

sum(a_vec)/10000000*100

PyPlot.figure()
PyPlot.plot(eval_point, E_visit_counter, "*-")
PyPlot.xlabel("Iteration")
PyPlot.ylabel("Energy")

PyPlot.figure()
PyPlot.plot(E_vec)
PyPlot.xlabel("Iteration")
PyPlot.ylabel("Energy")

PyPlot.figure()
PyPlot.plot(eval_point, log_g_tilde)
PyPlot.xlabel("Energy")
PyPlot.ylabel(L"log \tilde{g}")
