# load packages

using PyPlot
using Printf
using LaTeXStrings
using StatsBase

# ess function

function ess(x,k=1000)

    n = length(x)

    acf = autocor(x,1:k)

    return n/(1+2*sum(acf))

end


# generate start condiguration for stat S
L = 10
N = L*L
E_min = -2*N

q = 10 # spins test with q = 20
J = 1 # interaction strength (we have the same interaction strength for all states)

S_start = ones(L,L)

map!(x -> x = rand(1:q), S_start, S_start)

PyPlot.figure()
PyPlot.imshow(S_start,cmap="hot", interpolation="nearest")
PyPlot.colorbar()

T_c_exact = 1/log(1+sqrt(q))


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
# Metroplis
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
T = 0.588 #0.7145 # run metroplis at T_c
β = 1/T # tempering
iter = 2*10^6 # nbr of MC interstions

# run Metroplis algorithm
run_time_metroplis = @elapsed energy_vec, a_vec = metroplis(S_start, iter, J, β)

# compute avg acc prob
@printf "Avg. acc. prob.: %.2f %%\n" sum(a_vec)/iter*100

# plotting
PyPlot.figure(figsize=(8,5))
PyPlot.plot(1:iter, energy_vec)
PyPlot.xlabel("Iteration")
PyPlot.ylabel("Energy")


# calc ess
ess_metroplis = ess(energy_vec)
ess_per_sec_metroplis = ess_metroplis/run_time_metroplis



println(run_time_metroplis)
println(ess_metroplis)
println(ess_per_sec_metroplis)

################################################################################
# Wang-Landau
################################################################################

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
    iter_add = 500000

    for i in 1:nbr_reps

        map!(x -> x = rand(1:q), S_start, S_start)

        a_vec, E_vec, E_visit_counter, iter_done = wang_landau_one_iteration(S_start, iter, J, q, f, log_g_tilde, E_matrix)

        iter = iter + iter_add

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



# algorithm settings
iter = 2*10^6 # nbr of MC iterations
f = 2.7
nbr_reps = 25 # such that exp(1)^((1/2)^25) \approx exp(10^(-8))

# generate start condiguration for stat S
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


# only run first step
run_time_wl_one_step = @elapsed a_vec, E_vec, E_visit_counter, iter_done = @time wang_landau_one_iteration(S_start, iter, J, q, f, log_g_tilde, E_matrix)

# plotting
PyPlot.figure(figsize=(8,5))
PyPlot.plot(1:iter, E_vec)
PyPlot.xlabel("Iteration")
PyPlot.ylabel("Energy")

ess_wl_one_step = ess(E_vec)
ess_per_sec_wl_one_step = ess_wl_one_step/run_time_wl_one_step



println(run_time_wl_one_step)
println(ess_wl_one_step)
println(ess_per_sec_wl_one_step)

# run full w-l
f_save = @time wang_landau(nbr_reps,iter, J, q, f, S_start, log_g_tilde, E_matrix)


PyPlot.figure()
PyPlot.plot(f_save)
PyPlot.xlabel("Step")
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


# analytical value for T_c
T_c_exact = 1/log(1+sqrt(q)) # only holds for infinint system dimension


T = 0.588 #0.7138 #0.7145

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
