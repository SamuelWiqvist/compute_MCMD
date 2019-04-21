using PyPlot
using Printf


################################################################################
# Warmup: simulate the Ising model using the Metroplis algorithm
################################################################################

# generate start condiguration for stat S
dims = 100
S_start = zeros(dims,dims)
map!(x -> x = rand([1,-1]), S_start, S_start)

PyPlot.figure()
PyPlot.imshow(S_start, cmap="hot", interpolation="nearest")

# the energy function for the ising model
function H_ising(S, J)

    energy = 0

    for i = 1:size(S,1) # column
        for j = 1:size(S,2) # row

            # nearest neighbors for intearior states
            if i > 1 && j > 1 && i < size(S,1) && j < size(S,2)
                energy= energy + S[i,j]*S[i+1,j]
                energy= energy + S[i,j]*S[i-1,j]
                energy= energy + S[i,j]*S[i,j+1]
                energy= energy + S[i,j]*S[i,j-1]
            end

            # nearest neighbors for corner in first column
            if i == 1 && j == 1
                energy= energy + S[i,j]*S[i,j+1]
                energy= energy + S[i,j]*S[i+1,j]
            end

            if i == size(S,1) && j == 1
                energy= energy + S[i,j]*S[i,j+1]
                energy= energy + S[i,j]*S[i-1,j]
            end

            if i == 1 && j == size(S,2)
                energy= energy + S[i,j]*S[i,j-1]
                energy= energy + S[i,j]*S[i+1,j]
            end

            if i == size(S,1) && j == size(S,2)
                energy= energy + S[i,j]*S[i,j-1]
                energy= energy + S[i,j]*S[i-1,j]
            end

            # nearest neighbors for states in first column
            if i == 1 && j > 1 && j <= size(S,2)-1
                energy= energy + S[i,j]*S[i,j+1]
                energy= energy + S[i,j]*S[i+1,j]
                energy= energy + S[i,j]*S[i,j-1]
            end

            # nearest neighbors for states in last column
            if i == size(S,1) && j > 1 && j <= size(S,2)-1
                energy= energy + S[i,j]*S[i,j+1]
                energy= energy + S[i,j]*S[i-1,j]
                energy= energy + S[i,j]*S[i,j-1]
            end

            # nearest neighbors for states in first row
            if j == 1 && i > 1 && i <= size(S,1)-1
                energy= energy + S[i,j]*S[i,j+1]
                energy= energy + S[i,j]*S[i-1,j]
                energy= energy + S[i,j]*S[i+1,j]
            end

            # nearest neighbors for states in first last
            if j == size(S,2) && i > 1 && i <= size(S,1)-1
                energy= energy + S[i,j]*S[i,j-1]
                energy= energy + S[i,j]*S[i-1,j]
                energy= energy + S[i,j]*S[i+1,j]
            end
         end
    end

    return -J*energy

end

# metroplis algorithm to simulate the ising model
function metroplis(S_start, iter, J, β)

    nbr_stats = length(S_start) # set up
    S_configs = zeros(iter, size(S_start,1),size(S_start,2))
    a_vec = zeros(iter)
    energy_vec = zeros(iter)
    S_configs[1,:,:] = S_start

    H_old = H_ising(S_configs[1,:,:], J) # first iteration
    energy_vec[1] = H_old
    α_log = log(1)

    for i = 2:iter

        S_update = S_configs[i-1,:,:] # select state to flip at random
        S_flip = rand(1:nbr_stats) # flip state
        S_update[S_flip] = -1*S_update[S_flip] # set prop config

        H_new = H(S_update,J) # compute energy for prop config

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
sum(a_vec)/iter

# plotting
PyPlot.figure()
PyPlot.plot(1:iter, energy_vec)

iter_plot = floor.(Int,LinRange(1,iter,100))

PyPlot.figure()
for i = iter_plot
    sleep(0.001)
    PyPlot.imshow(S_configs[i,:,:], cmap="hot", interpolation="nearest")
    PyPlot.xlabel(i)
end

PyPlot.figure()
PyPlot.imshow(S_configs[1,:,:], cmap="hot", interpolation="nearest")

PyPlot.figure()
PyPlot.imshow(S_configs[1000,:,:], cmap="hot", interpolation="nearest")

PyPlot.figure()
PyPlot.imshow(S_configs[10000,:,:], cmap="hot", interpolation="nearest")

PyPlot.figure()
PyPlot.imshow(S_configs[100000,:,:], cmap="hot", interpolation="nearest")

PyPlot.figure()
PyPlot.imshow(S_configs[150000,:,:], cmap="hot", interpolation="nearest")

PyPlot.figure()
PyPlot.imshow(S_configs[200000,:,:], cmap="hot", interpolation="nearest")

################################################################################
# Potts model 2D q = 10
################################################################################


# generate start condiguration for stat S
dims = 100
q = 10
S_start = zeros(dims,dims)
map!(x -> x = rand(1:q), S_start, S_start)

PyPlot.figure()
PyPlot.imshow(S_start,cmap="hot", interpolation="nearest")
PyPlot.colorbar()

# the energy function for the ising model
function H_potts(S, J)

    energy = 0

    for i = 1:size(S,1) # column
        for j = 1:size(S,2) # row

            # nearest neighbors for intearior states
            if i > 1 && j > 1 && i < size(S,1) && j < size(S,2)
                energy= energy + δ(S[i,j],S[i+1,j])
                energy= energy + δ(S[i,j],S[i-1,j])
                energy= energy + δ(S[i,j],S[i,j+1])
                energy= energy + δ(S[i,j],S[i,j-1])
            end

            # nearest neighbors for corner in first column
            if i == 1 && j == 1
                energy= energy + δ(S[i,j],S[i,j+1])
                energy= energy + δ(S[i,j],S[i+1,j])
            end

            if i == size(S,1) && j == 1
                energy= energy + δ(S[i,j],S[i,j+1])
                energy= energy + δ(S[i,j],S[i-1,j])
            end

            if i == 1 && j == size(S,2)
                energy= energy + δ(S[i,j],S[i,j-1])
                energy= energy + δ(S[i,j],S[i+1,j])
            end

            if i == size(S,1) && j == size(S,2)
                energy= energy + δ(S[i,j],S[i,j-1])
                energy= energy + δ(S[i,j],S[i-1,j])
            end

            # nearest neighbors for states in first column
            if i == 1 && j > 1 && j <= size(S,2)-1
                energy= energy + δ(S[i,j],S[i,j+1])
                energy= energy + δ(S[i,j],S[i+1,j])
                energy= energy + δ(S[i,j],S[i,j-1])
            end

            # nearest neighbors for states in last column
            if i == size(S,1) && j > 1 && j <= size(S,2)-1
                energy= energy + δ(S[i,j],S[i,j+1])
                energy= energy + δ(S[i,j],S[i-1,j])
                energy= energy + δ(S[i,j],S[i,j-1])
            end

            # nearest neighbors for states in first row
            if j == 1 && i > 1 && i <= size(S,1)-1
                energy= energy + δ(S[i,j],S[i,j+1])
                energy= energy + δ(S[i,j],S[i-1,j])
                energy= energy + δ(S[i,j],S[i+1,j])
            end

            # nearest neighbors for states in first last
            if j == size(S,2) && i > 1 && i <= size(S,1)-1
                energy= energy + δ(S[i,j],S[i,j-1])
                energy= energy + δ(S[i,j],S[i-1,j])
                energy= energy + δ(S[i,j],S[i+1,j])
            end
         end
    end

    return -J*energy

end

function δ(s, s_star)

    if s == s_star
        return 1
    else
        return 0
    end

end


# Wang-Landau algorithm with one fixed t value
function wang_landau(S_start, iter, J, q, f)

    nbr_stats = length(S_start) # set up
    S_configs = zeros(iter, size(S_start,1),size(S_start,2))
    a_vec = zeros(iter)

    g_tilde = zeros(2,iter) #[g_tilde, E]
    g_tilde[1,:] = ones(iter)

    S_configs[1,:,:] = S_start
    H_old = H_potts(S_configs[1,:,:], J) # first iteration

    g_tilde[1,1] = 1
    g_tilde[2,1] = H_old
    α_log = log(1)

    for i = 2:iter

        # ordinary update
        S_update = S_configs[i-1,:,:] # select state to flip at random
        S_flip = rand(1:nbr_stats) # flip state
        S_update[S_flip] = rand(1:q) # set prop config

        H_new = H_potts(S_update,J) # compute energy for prop config

        # see if we have curent energy
        g_tilde_prop_idx = findlast(x -> x == H_new, g_tilde[2,:])

        if typeof(g_tilde_prop_idx) == Nothing
            g_tilde_prop = 1
        else
            g_tilde_prop = g_tilde[1,g_tilde_prop_idx[1]]
        end

        g_tilde_old = g_tilde[1,i-1]

        α_log = log(g_tilde_old)-log(g_tilde_prop)

        if log(rand()) < min(0, α_log) # accapt new config
            S_configs[i,:,:] = S_update
            g_tilde[2,i] = H_new
            a_vec[i] = 1
        else
            S_configs[i,:,:] = S_configs[i-1,:,:] # store old config
            g_tilde[2,i] = g_tilde[2,i-1]
        end

        # update energy function for current system
        g_update_idx = findlast(x -> x == g_tilde[2,i], g_tilde[2,1:i-1])

        if typeof(g_update_idx) == Nothing
            g_update = 1
        else
            g_update = g_tilde[1,g_update_idx[1]]
        end

        g_tilde[1,i] = f*g_update

    end

    return S_configs,g_tilde,a_vec

end


# Full Wang-Landau algorithm where the t value is decreased
function full_wang_landau(nbr_reps,iter, J, q, f, dims = 100)

    println("Starting Wang-Landau.")

    # full Wang-Landau algorithm
    S_start = zeros(dims,dims)
    g_save = zeros(nbr_reps,2,iter)
    f_save = zeros(nbr_reps)
    S_configs_last = zeros(iter, size(S_start,1),size(S_start,2))

    for i = 1:nbr_reps

        # run algorithm
        S_start = map!(x -> x = rand(1:q), S_start, S_start)

        S_configs, g_tilde, a_vec = wang_landau(S_start, iter, J, q, f)

        g_save[i,:,:] = g_tilde
        f_save[i] = f

        # print info
        @printf "-----------------------------------------------------------------\n"
        @printf "Iter: %f\n" i
        @printf "f: %f\n" f
        @printf "Acc. rate: %f %%\n" sum(a_vec)/iter*100

        # update t
        f = sqrt(f)

        if i == nbr_reps
            S_configs_last = S_configs
        end

    end

    return g_save, f_save, S_configs_last

end



J = 1 # interaction strength (we have the same interaction strength for all states)
iter = 10000 # nbr of MC interstions
f = 10
nbr_reps = 26 # such that 10^((1/2)^26) \approx exp(10^(-8))
q = 10

g_save, f_save, S_configs_last = @time full_wang_landau(nbr_reps, iter, J, q, f)


PyPlot.figure()
PyPlot.plot(f_save)
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
#S_configs, g_tilde, a_vec = @time wang_landau(S_start, iter, J, q, t)
