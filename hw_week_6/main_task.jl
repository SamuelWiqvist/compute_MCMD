using QuadGK

E_0 = 50.
pmax = E_0
p2max = pmax^2
α_s = 0.2

struct Decay{A,B}
    parent::A
    product::A
    E::B
    p::B
end

function I_a_bc(parent::String, z_min::Real, z_max::Real)

    if parent == "q"
        return α_s/(2*π)*(8/3)*log((1-z_min)/(1-z_max))
    else
        return  α_s/(2*π)*(z_max-z_min)
    end
end

function shower()

    steps = 0
    nbr_particles_old = 1

    p_new = 0.
    z_new = 0.

    #for j = 1:2

    simulate = true

    #while simulate

    for i = 1:3

        println("-----")
        steps = steps + 1
        nbr_particles_generate = 2*nbr_particles_old
        nbr_diff = nbr_particles_generate - nbr_particles_old
        nbr_dead = 0
        # check if all Decay s are dead then end program

        println(nbr_particles_generate)

        for i in nbr_diff:div(nbr_diff+nbr_particles_generate-1,2)

            parent = data[i].parent
            Ea = data[i].E
            p_old = data[i].p

            if p_old == 0

                Pa = Decay("dead", "dead", 0., 0.)
                Pb = Decay("dead", "dead", 0., 0.)

                append!(data,[Pa])
                append!(data,[Pb])

                nbr_dead =  nbr_dead + 1
            else

                #println(parent)
                #println(Ea)
                #println(p_old)

                sample = true

                while sample
                    # veto step in time, code adapted from https://github.com/kalaee/MCMD/tree/master/parton


                    z_min = 1/Ea
                    z_max = 1-1/Ea

                    #println("...................")

                    #println(z_min)
                    #println(z_max)


                    I = I_a_bc(parent, z_min, z_max)

                    #println(I)
                    p2_new = p_old^2*rand()^(1/I)
                    p_new = sqrt(p2_new)
                    #println(p_new)

                    if p_new <= 1
                        z_new = 0
                        p_new = 0
                        parent = "dead"
                    else

                        if parent == "q"
                            z_new = 1 - (1 - z_min)*((1-z_max)/(1-z_min))^(rand())
                        else
                            z_new = - z_min + rand()*(z_max-z_min)
                        end

                        p_old = p_new

                        #println(z_new)


                        if z_new <= z_max && z_new >= z_min
                            #println("test")
                            a_prob = (1 + z_new^2)/2
                            #println(a_prob)
                            if a_prob < rand()
                                sample = false
                            end
                        else
                            # contine veto alg.
                        end

                    end

                end

                Eb = z_new*Ea
                Ec = (1-z_new)*Ea

                if parent == "q"
                    Pa = Decay(parent, "q", Eb, p_new)
                    Pb = Decay(parent, "g", Ec, p_new)
                elseif parent == "g"
                    Pa = Decay(parent, "g", Eb, p_new)
                    Pb = Decay(parent, "g", Ec, p_new)
                else
                    Pa = Decay(parent, "dead", 0., 0.)
                    Pb = Decay(parent, "dead", 0., 0.)
                end


                append!(data,[Pa])
                append!(data,[Pb])

            end

        end

        nbr_particles_old = nbr_particles_generate

        println(nbr_dead)

        if nbr_dead == nbr_particles_generate
            simulate = false
        end

    end


end



p = Decay("q", "q", E_0, E_0)


data = [p]


println("...................")
shower()

data
