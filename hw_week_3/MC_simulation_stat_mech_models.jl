using PyPlot

# Warmup simulte the Ising model using the Metroplis algorighm

# generate start condiguration for stat S
dim = 100
S_start = zeros(dim,dim)
map!(x -> x = rand([1,-1]), S_start, S_start)

PyPlot.figure()
PyPlot.imshow(S_start,cmap="hot", interpolation="nearest")

# the energy function for the ising model
function H(S, J)

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

    nbr_stats = length(S_start)

    S_configs = zeros(iter, size(S_start,1),size(S_start,2))
    a_vec = zeros(iter)
    energy_vec = zeros(iter)
    S_configs[1,:,:] = S_start

    H_old = H(S_configs[1,:,:], J)
    energy_vec[1] = H_old

    α_log = log(1)

    for i = 2:iter

        # select state to flip at random
        S_update = S_configs[i-1,:,:]
        S_flip = rand(1:nbr_stats)
        S_update[S_flip] = -1*S_update[S_flip]

        H_new = H(S_update,J)

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
            S_configs[i,:,:] = S_configs[i-1,:,:]
            energy_vec[i] = H_old

        end

    end

    return S_configs,energy_vec,a_vec

end



# set system parameters
J = 0.5 # interaction strength (we have the same interaction strength for all states)
β = 20 # tempering
iter = 200000 # nbr of MC interstions


S_configs,energy_vec, a_vec = @time metroplis(S_start, iter, J, β)


sum(a_vec)/iter

PyPlot.figure()
PyPlot.plot(1:iter, energy_vec)

iter_plot = floor.(Int,LinRange(1,iter,100))

PyPlot.figure()
for i = iter_plot
    sleep(0.001)
    PyPlot.imshow(S_configs[i,:,:],cmap="hot", interpolation="nearest")
    PyPlot.xlabel(i)
end



PyPlot.figure()
PyPlot.imshow(S_configs[1,:,:],cmap="hot", interpolation="nearest")

PyPlot.figure()
PyPlot.imshow(S_configs[1000,:,:],cmap="hot", interpolation="nearest")

PyPlot.figure()
PyPlot.imshow(S_configs[10000,:,:],cmap="hot", interpolation="nearest")

PyPlot.figure()
PyPlot.imshow(S_configs[100000,:,:],cmap="hot", interpolation="nearest")

PyPlot.figure()
PyPlot.imshow(S_configs[150000,:,:],cmap="hot", interpolation="nearest")

PyPlot.figure()
PyPlot.imshow(S_configs[200000,:,:],cmap="hot", interpolation="nearest")
