using PyPlot
using KernelDensity


include("data.jl")

PyPlot.figure()


'''
    generate_photon_path_lengths(μ::Real, N::Int)

Function to sample N photon path lengths for some attenuation coefficient μ.
'''
function generate_photon_path_lengths(μ::Real, N::Int)
    d = zeros(N)
    for i in 1:N; d[i] = -1/μ*log(rand()); end
    return d
end



PyPlot.figure()
for i = 1:size(data_h2o,1)
    μ = sum(data_h2o[i,end-1])
    N = 10^3
    d = generate_photon_path_lengths(μ, N)
    kde_est = kde(d)
    PyPlot.plot(kde_est.x, kde_est.density, label=string(data_h2o[i,1])*"KeV")
    PyPlot.title("H2O")
end
PyPlot.legend()

PyPlot.figure()
for i = 1:size(data_h2o,1)
    μ = sum(data_Al[i,end-1])
    N = 10^3
    d = generate_photon_path_lengths(μ, N)
    kde_est = kde(d)
    PyPlot.plot(kde_est.x, kde_est.density, label=string(data_Al[i,1])*"KeV")
    PyPlot.title("Al")
end
PyPlot.legend()



PyPlot.figure()
for i = 1:size(data_h2o,1)
    μ = sum(data_I[i,end-1])
    N = 10^3
    d = generate_photon_path_lengths(μ, N)
    kde_est = kde(d)
    PyPlot.plot(kde_est.x, kde_est.density, label=string(data_I[i,1])*"KeV")
    PyPlot.title("I")
end
PyPlot.legend()


PyPlot.figure()
for i = 1:size(data_h2o,1)
    μ = sum(data_Pb[i,end-1])
    N = 10^3
    d = generate_photon_path_lengths(μ, N)
    kde_est = kde(d)
    PyPlot.plot(kde_est.x, kde_est.density, label=string(data_Pb[i,1])*"KeV")
    PyPlot.title("Pb")
end
PyPlot.legend()


function compton_scattering(energy_in::Float64, material::String, N_samples::Int)


    # check inputs
    if !(energy_in in data_h2o[:,1]) # check energy level
        println("Energy level not allowed, allowed energy levelse are:")
        println(data_h2o[:,1])
    end

    if !(material in ["H2O"; "Al"; "I"; "Pb"]) # check material
        println("Material not allowed, allowed materials are:")
        println(["H2O"; "Al"; "I"; "Pb"])
    end

    # object parameters
    distance_to_object = 0.5
    x_length_object = 2
    z_length_object = 2
    y_length_object = 1


    hit_agels = zeros(2,N_samples)
    hit_possition = zeros(3,N_samples)

    nbr_samples = 0
    for i = 1:N_samples
        generate_new_sample = true
        while generate_new_sample
            generate_new_sample = true
            nbr_samples = nbr_samples + 1

            # sample the start angels from a point source at the origin
            θ = 2*rand() - 1 # sample polar angels
            ϕ = 2*π*rand() # sample azimuthal angel

            # compute position at object
            r_object = distance_to_object/(sin(θ)*sin(ϕ))
            x_object = r_object*sin(θ)*cos(ϕ)
            z_object = r_object*cos(θ)

            # check if we hit the object
            if x_object >= -x_length_object/2 && x_object <= x_length_object/2
                if z_object >= -z_length_object/2 && z_object <= z_length_object/2
                    generate_new_sample = false
                    hit_possition[:,i] = [x_object;distance_to_object;z_object]
                    hit_agels[:,i] = [θ;ϕ]
                end
            end
        end
    end

    println("Nbr samples used:")
    println(nbr_samples)

    # find attenuation coefficient μ
    idx_energy_levle = findfirst(x-> x == 1, data_h2o[:,1])

    if  material == "H2O"
        μ = data_h2o[idx_energy_levle,2]
    elseif material == "Al"
        μ = data_Al[idx_energy_levle,2]
    elseif material == "I"
        μ = data_I[idx_energy_levle,2]
    elseif material == "Pb"
        μ = data_Pb[idx_energy_levle,2]
    end

    # sample path-lengths
    d = generate_photon_path_lengths(μ, N_samples)


    # update possion for particels
    new_possition = zeros(3,N_samples)



    # check if a particel has left the object
    return d, hit_possition, hit_agels
end

data_h2o[:,1]

d, hit_possition, hit_agels = compton_scattering(0.05, "Al", 1)

findfirst(x-> x == 1, data_h2o[:,1])

hit_possition[:,1]

distance_to_object = 0.5
x_length_object = 2
z_length_object = 2
y_length_object = 1

θ = 2*rand() - 1 # sample polar angels
ϕ = 2*π*rand() # sample azimuthal angel
r = 1

r_object = distance_to_object/(sin(θ)*sin(ϕ))
x_object = r_object*sin(θ)*cos(ϕ)
z_object = r_object*cos(θ)

x_object >= -x_length_object/2 && x_object <= x_length_object/2
z_object >= -z_length_object/2 && z_object <= z_length_object/2

x_new(t) = t*x
