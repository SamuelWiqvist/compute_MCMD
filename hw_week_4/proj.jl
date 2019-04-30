using PyPlot
using KernelDensity

# task 1 (???)

# data for cross-sections for H2O, Al, I,

# Data for cross-sections for: water (H2O), aluminium (Al), Iodine (I) and Lead (Pb).
# The data is collected from: https://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html


# Photon    Coherent Incoher. Photoel. Tot. w/  Tot. wo/
# Energy    Scatter. Scatter. Absorb.  Coherent Coherent

data_h2o = [5.000E-02 1.936E-02 1.803E-01 2.725E-02 2.269E-01 2.076E-01;
            6.000E-02 1.392E-02 1.770E-01 1.493E-02 2.059E-01 1.920E-01;
            #8.000E-02 8.165E-03 1.697E-01 5.770E-03 1.837E-01 1.755E-01;
            1.000E-01 5.349E-03 1.626E-01 2.763E-03 1.707E-01 1.654E-01;
            #1.500E-01 2.442E-03 1.474E-01 7.308E-04 1.505E-01 1.481E-01;
            2.000E-01 1.388E-03 1.353E-01 2.887E-04 1.370E-01 1.356E-01;
            #3.000E-01 6.215E-04 1.179E-01 8.160E-05 1.186E-01 1.180E-01;
            4.000E-01 3.506E-04 1.058E-01 3.493E-05 1.061E-01 1.058E-01;
            #5.000E-01 2.247E-04 9.663E-02 1.883E-05 9.687E-02 9.665E-02;
            6.000E-01 1.561E-04 8.939E-02 1.173E-05 8.956E-02 8.940E-02;
            #8.000E-01 8.790E-05 7.856E-02 5.920E-06 7.866E-02 7.857E-02;
            1.000E+00 5.627E-05 7.066E-02 3.681E-06 7.072E-02 7.066E-02]

data_Al =  [5.000E-02 4.678E-02 1.496E-01 1.718E-01 3.681E-01 3.214E-01;
            6.000E-02 3.386E-02 1.483E-01 9.564E-02 2.778E-01 2.440E-01;
            #8.000E-02 2.005E-02 1.439E-01 3.783E-02 2.018E-01 1.817E-01;
            1.000E-01 1.324E-02 1.388E-01 1.840E-02 1.704E-01 1.572E-01;
            #1.500E-01 6.122E-03 1.267E-01 4.993E-03 1.378E-01 1.317E-01;
            2.000E-01 3.504E-03 1.168E-01 2.002E-03 1.223E-01 1.188E-01;
            #3.000E-01 1.580E-03 1.021E-01 5.743E-04 1.042E-01 1.026E-01;
            4.000E-01 8.934E-04 9.162E-02 2.480E-04 9.276E-02 9.187E-02;
            #5.000E-01 5.732E-04 8.374E-02 1.344E-04 8.445E-02 8.388E-02;
            6.000E-01 3.986E-04 7.754E-02 8.401E-05 7.802E-02 7.762E-02;
            #8.000E-01 2.245E-04 6.814E-02 4.252E-05 6.841E-02 6.818E-02;
            1.000E+00 1.438E-04 6.129E-02 2.643E-05 6.146E-02 6.132E-02]

data_I =   [5.000E-02 3.452E-01 1.093E-01 1.187E+01 1.232E+01 1.198E+01;
            6.000E-02 2.577E-01 1.110E-01 7.208E+00 7.577E+00 7.319E+00;
            #8.000E-02 1.584E-01 1.115E-01 3.240E+00 3.510E+00 3.352E+00;
            1.000E-01 1.073E-01 1.100E-01 1.725E+00 1.942E+00 1.835E+00;
            #1.500E-01 5.220E-02 1.036E-01 5.419E-01 6.978E-01 6.456E-01;
            2.000E-01 3.087E-02 9.709E-02 2.383E-01 3.663E-01 3.354E-01;
            #3.000E-01 1.440E-02 8.627E-02 7.650E-02 1.772E-01 1.628E-01;
            4.000E-01 8.295E-03 7.811E-02 3.525E-02 1.217E-01 1.134E-01;
            #5.000E-01 5.381E-03 7.175E-02 1.988E-02 9.701E-02 9.163E-02;
            6.000E-01 3.771E-03 6.663E-02 1.274E-02 8.313E-02 7.936E-02;
            #8.000E-01 2.143E-03 5.875E-02 6.596E-03 6.749E-02 6.534E-02;
            1.000E+00 1.379E-03 5.291E-02 4.120E-03 5.841E-02 5.703E-02]


data_Pb =  [5.000E-02 6.545E-01 9.478E-02 7.292E+00 8.042E+00 7.387E+00;
            6.000E-02 4.900E-01 9.734E-02 4.432E+00 5.020E+00 4.530E+00;
            #8.000E-02 3.078E-01 9.923E-02 2.012E+00 2.419E+00 2.112E+00;
            8.800E-02 2.632E-01 9.928E-02 1.547E+00 1.910E+00 1.647E+00;
            #8.800E-02 2.632E-01 9.928E-02 7.321E+00 7.684E+00 7.421E+00;
            1.000E-01 2.128E-01 9.894E-02 5.237E+00 5.549E+00 5.336E+00;
            #1.500E-01 1.049E-01 9.484E-02 1.815E+00 2.015E+00 1.910E+00;
            2.000E-01 6.260E-02 8.966E-02 8.464E-01 9.986E-01 9.360E-01;
            #3.000E-01 2.988E-02 8.036E-02 2.930E-01 4.032E-01 3.733E-01;
            4.000E-01 1.746E-02 7.310E-02 1.417E-01 2.323E-01 2.148E-01;
            #5.000E-01 1.143E-02 6.734E-02 8.257E-02 1.613E-01 1.499E-01;
            6.000E-01 8.060E-03 6.263E-02 5.406E-02 1.248E-01 1.167E-01;
            #8.000E-01 4.621E-03 5.537E-02 2.871E-02 8.870E-02 8.408E-02;
            1.000E+00 2.991E-03 4.993E-02 1.810E-02 7.102E-02 6.803E-02]

'''
    generate_photon_path_lengths(μ::Real, N::Int)

function to sample N photon path lengths for some attenuation coefficient μ.
'''
function generate_photon_path_lengths(μ::Real, N::Int)
    d = zeros(N)
    for i in 1:N; d[i] = -1/μ*log(rand()); end
    return d
end



PyPlot.figure()
for i = 1:size(data_h2o,1)
    μ = sum(data_h2o[i,2:end])
    N = 10^3
    d = generate_photon_path_lengths(μ, N)
    kde_est = kde(d)
    PyPlot.plot(kde_est.x, kde_est.density, label=string(data_h2o[i,1])*"KeV")
    PyPlot.title("H2O")
end
PyPlot.legend()

PyPlot.figure()
for i = 1:size(data_h2o,1)
    μ = sum(data_Al[i,2:end])
    N = 10^3
    d = generate_photon_path_lengths(μ, N)
    kde_est = kde(d)
    PyPlot.plot(kde_est.x, kde_est.density, label=string(data_Al[i,1])*"KeV")
    PyPlot.title("Al")
end
PyPlot.legend()



PyPlot.figure()
for i = 1:size(data_h2o,1)
    μ = sum(data_I[i,2:end])
    N = 10^3
    d = generate_photon_path_lengths(μ, N)
    kde_est = kde(d)
    PyPlot.plot(kde_est.x, kde_est.density, label=string(data_I[i,1])*"KeV")
    PyPlot.title("I")
end
PyPlot.legend()


PyPlot.figure()
for i = 1:size(data_h2o,1)
    μ = sum(data_Pb[i,2:end])
    N = 10^3
    d = generate_photon_path_lengths(μ, N)
    kde_est = kde(d)
    PyPlot.plot(kde_est.x, kde_est.density, label=string(data_Pb[i,1])*"KeV")
    PyPlot.title("Pb")
end
PyPlot.legend()
