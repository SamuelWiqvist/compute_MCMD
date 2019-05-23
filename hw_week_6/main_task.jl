using QuadGK

E_0 = 50.
p2max = E_0^2
α_s = 0.2

struct Particle{A,B}
    parent::A
    type::A
    E::B
    p::B
end


p = Particle("start", "q", E_0, E_0)


data = [p]

P_q_qg(z) = 4/3*(1+z^2)/(1-z)
P_g_gg(z) = 3*(1-z*(1-z))^2/(z*(1-z))


d[2]


data[1].type

d =
typeof(data)


Particle("q", E_0, E_0)


append!(d,[p])

P_q_qg(z) = 4/3*(1+z^2)/(1-z)
P_g_gg(z) = 3*(1-z*(1-z))^2/(z*(1-z))

P_q_qg(0)
P_q_qg(0.99)
10^-6

g_q_qg = quadgk(P_q_qg, 0, 1-10^-6, rtol=1e-3)[1]
g_g_qg = quadgk(P_g_gg, 10^-6, 1-10^-6, rtol=1e-3)[1]

function hitandmisss(f)

    y_max = f(1-10^-6)
    generate = true
    x_star = 0
    while generate

        x_star = 10^-6 + (1-10^-6-10^-6)*rand()
        y_star = y_max*rand()

        if y_star <= f(x_star)
            generate = false
        end
    end

    return x_star

end

hitandmisss(P_g_gg)



function f(p2, z, Ea, P)

    int_p2 = log(p2max) - log(p2)

    z_min = sqrt(p2)/Ea
    z_max = 1-sqrt(p2)/Ea

    println("p2")
    println(p2)
    println("z_min")

    println(z_min)
    println("z_max")

    println(z_max)

    println("Ea")
    println(Ea)

    int_z = quadgk(P, z_min, z_max, rtol=1e-3)[1]

    S = exp(-int_p2*α_s/π*int_z)

    return 1/p2*α_s/π*S

end

f(2499.98226035615, 0.865, 50, P_q_qg)

a = 1
function shower()

    steps = 0
    nbr_particles_old = 1

    g_q_qg = quadgk(P_q_qg, 0.0001, 1-0.0001, rtol=1e-3)[1]
    g_g_qg = quadgk(P_g_gg, 0.0001, 1-0.0001, rtol=1e-3)[1]

    p2_new = 0.
    z_new = 0.

    for j = 1:1

        steps = steps + 1
        nbr_particles_generate = 2*nbr_particles_old
        nbr_diff = nbr_particles_generate - nbr_particles_old

        println(nbr_particles_old)
        println(nbr_particles_generate)
        println(nbr_diff)

        # check if all particles are dead then end program

        for i in nbr_diff:div(nbr_diff+nbr_particles_generate-1,2)

            generate = true

            println(i)

            type_old = data[i].type
            Ea = data[i].E
            p2_old = data[i].p^2

            while generate

                if type_old == "q"
                    g = g_q_qg
                else
                    g = g_g_qg
                end

                p2_new = (g*p2_old + log(rand()))/g

                if sqrt(p2_new) <= 1
                    p2_new = -1
                    z_new = -1
                else

                    if type_old == "q"
                        z_new = hitandmisss(P_q_qg)
                        g_val = P_q_qg(z_new)
                        f_val = f(p2_new, z_new, Ea, P_q_qg)
                    else
                        z_new = hitandmisss(P_g_gg)
                        g_val = P_g_gg(z_new)
                        f_val = f(p2_new, z_new, Ea, P_g_gg)
                    end

                    if f_val/g_val <= rand()
                        generate = false
                    end
                end

            end

            println(p2_new)
            println(z_new)
            Eb = z_new*Ea
            Ec = (1-z_new)*Ea

            Pa = Particle(type_old, "q", Eb, sqrt(p2_new))
            Pb = Particle(type_old, "g", Ec, sqrt(p2_new))

            append!(data,[Pa])
            append!(data,[Pb])


        end

        nbr_particles_old = nbr_particles_generate

    end

end


p = Particle("start", "q", E_0, E_0)


data = [p]


println("...................")
shower()

data
