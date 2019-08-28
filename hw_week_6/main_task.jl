import Statistics.mean
import Statistics.std
using Printf

struct Decay{A,B}
    parent::A
    product::A
    E::B
    p::B
end

function I_a_bc(parent::String, z_min::Real, z_max::Real, α_s::Real = 0.2)

    if parent == "q"
        return α_s/(2*π)*(8/3)*log((1-z_min)/(1-z_max))
    else
        return  α_s/(2*π)*(z_max-z_min)
    end
end

function shower(data::Vector)

    nbr_particles_old = 1

    p_new = 0.
    z_new = 0.

    #for j = 1:2

    simulate = true

    steps = 0
    #while simulate

    while simulate

        steps = steps + 1

        nbr_particles_generate = 2*nbr_particles_old
        nbr_diff = nbr_particles_generate - nbr_particles_old
        nbr_dead = 0

        if steps == 1
            idx_old =  1:1
        else
            idx_old = nbr_diff:nbr_diff+nbr_particles_old-1
        end

        for i in idx_old

            parent = data[i].product
            Ea = data[i].E
            p_old = data[i].p

            if p_old == 0

                nbr_dead =  nbr_dead + 2

                p_new = 0
                z_new = 0
                parent = "dead"

            else


                sample = true

                while sample
                    # veto step in time, code adapted from https://github.com/kalaee/MCMD/tree/master/parton


                    z_min_tilde = 1/Ea # what should these values be???
                    z_max_tilde = 1-1/Ea


                    I = I_a_bc(parent, z_min_tilde, z_max_tilde)

                    println(z_min_tilde)
                    println(z_max_tilde)


                    if I < 0
                        println(parent)
                        println(I)
                        println(Ea)
                        println(z_min_tilde)
                        println(z_max_tilde)
                        sample = false
                        simulate = false
                    end

                    p2_new = p_old^2*rand()^(1/I)
                    p_new = sqrt(p2_new)


                    if p_new <= 1 #|| p_new == Inf
                        z_new = 0.
                        p_new = 0.
                        parent = "dead"
                        sample = false
                    else

                        if parent == "q"
                            z_new = 1 - (1 - z_min_tilde)*((1-z_max_tilde)/(1-z_min_tilde))^(rand())
                        else
                            z_new = - z_min_tilde + rand()*(z_max_tilde-z_min_tilde)
                        end

                        p_old = p_new

                        z_min = p_new/Ea # what should these values be???
                        z_max = 1-p_new/Ea


                        if z_new <= z_max && z_new >= z_min
                            a_prob = (1 + z_new^2)/2
                            if a_prob < rand()
                                sample = false
                            end
                        else
                            # contine veto alg.
                        end

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

        nbr_particles_old = nbr_particles_generate


        if nbr_dead == nbr_particles_generate
            simulate = false
        end

    end


end



E_0 = 100.
pmax = E_0
p2max = pmax^2

p = Decay("q", "q", E_0, E_0)
data = [p]

println("...................")
shower(data)

println(length(data))

data

# analyses for one run
nbr_g = 0
nbr_q_qg = 0

for i in data
    if i.product == "g"
        global nbr_g = nbr_g + 1
    end

    if  i.parent == "q" && i.product == "g"
        global nbr_q_qg = nbr_q_qg + 1
    end
end

println(nbr_g)
println(nbr_q_qg)

# mutiple runs

nbr_runs = 100000
stats = zeros(2,nbr_runs)

for i = 1:nbr_runs

    #println(i)
    E_0 = 400.
    pmax = E_0
    p2max = pmax^2

    p = Decay("q", "q", E_0, E_0)
    data = [p]

    shower(data)

    nbr_g = 0
    nbr_q_qg = 0

    for i in data
        if i.product == "g"
            nbr_g = nbr_g + 1
        end

        if  i.parent == "q" && i.product == "g"
            nbr_q_qg = nbr_q_qg + 1
        end
    end

    stats[:,i] = [nbr_g;nbr_q_qg]

end

@printf "----------------\n"
@printf "nbr_g (mena): %f\n" mean(stats[1,:])
@printf "nbr_g (std): %f\n" std(stats[1,:])

@printf "nbr_q_qg (mena): %f\n" mean(stats[2,:])
@printf "nbr_q_qg (std): %f\n" std(stats[2,:])
