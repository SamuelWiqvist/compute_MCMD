using PyPlot

# Warmup simulte the Ising model using the Metroplis algorighm


dim = 500

# generate start condiguration for stat S
S_start = zeros(dim,dim)
map!(x -> x = rand([1,-1]), S_start, S_start)

PyPlot.figure()
PyPlot.imshow(S_start,cmap="hot", interpolation="nearest")
PyPlot.colorbar()

# set system parameters
J = 1 # interaction strength (we have the same interaction strength for all states)
β = 10 # tempering

# the energy function for the ising model
function H(S, J)

    energy = 0


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


S_start


# test energy function
@time H(S_start, J)
