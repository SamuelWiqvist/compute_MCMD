using PyPlot
using KernelDensity
using LaTeXStrings

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


function compton_scattering(hv::Float64)


    # check inputs
    if !(hv in data_h2o[:,1]) # check energy level
        println("Energy level not allowed, allowed energy levelse are:")
        println(data_h2o[:,1])
    end

    # set parameters for compton scattering
    h = 1; m_0 = 1; c = 1

    λ = hv / (m_0*c^2)
    Q = (2*λ + 1)/(2*λ + 9)

    run_sampler = true
    nbr_sampels = 0
    hv_prime = 0
    θ = 0
    cosθ = 0
    while run_sampler

        nbr_sampels = nbr_sampels + 1

        R1 = rand()
        R2 = rand()
        R3 = rand()

        if R1 < Q
            ρ = 1 + 2*λ*R2
            if R3 > (4*(ρ-1))/ρ^2
                # generate new sample
            else
                cosθ = 1-2*R2
                θ = acos(cosθ)
                θ = θ*180/π # TODO fix this
                hv_prime = hv/ρ
                run_sampler = false
            end
        else
            ρ = (2*λ + 1)/(2*λ*R2 + 1)
            if R3 > ((1-(ρ-1)/λ)^2 + 1/ρ)/2
                # generate new sample
            else
                cosθ = 1-(ρ-1)/λ
                θ = acos(cosθ)
                θ = θ*180/π
                hv_prime = hv/ρ
                run_sampler = false
            end
        end
    end

    return [hv_prime; θ; cosθ; nbr_sampels]

end


N = 100000

energy_levels = data_h2o[:,1]


energy_levels


samples_compton_scattering = zeros(4,N,length(energy_levels))


for i in 1:length(energy_levels)
    for j = 1:N
        samples_compton_scattering[:,j,i] = compton_scattering(energy_levels[i])
    end
end

PyPlot.figure()
h = PyPlot.plt[:hist](samples_compton_scattering[4,:,7],100)


PyPlot.figure(figsize=(10,40))

energy_level = 0

for i in 1:2:2*length(energy_levels)
    global energy_level = energy_level + 1
    PyPlot.subplot(7,2,i)
    h = PyPlot.plt[:hist](samples_compton_scattering[1,:,energy_level],100)
    PyPlot.subplot(7,2,i+1)
    h = PyPlot.plt[:hist](samples_compton_scattering[2,:,energy_level],100)
end



energy_level = 7

PyPlot.figure()
PyPlot.subplot(1,2,1)
h = PyPlot.plt[:hist](samples_compton_scattering[1,:,energy_level],100)
PyPlot.subplot(1,2,2)
h = PyPlot.plt[:hist](samples_compton_scattering[3,:,energy_level],100)
