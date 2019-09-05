import Statistics.mean
import Statistics.std
using Printf

struct Decay{A,B}
    parent::A
    product::A
    E::B
    p::B
end

P_q_qg(z::Real) = 4/3*(1+z^2)/(1-z)
P_g_gg(z::Real) = 3*(1-z*(1-z))^2/(z*(1-z))

function P_a_bc(parent::String, z::Real)

    if parent == "q"
        return P_q_qg(z)
    else
        return P_g_gg(z)
    end

end

function I_a_bc(parent::String, z_min::Real, z_max::Real, Ea::Real, α_s::Real = 0.2)

    if parent == "q"
        return α_s/(2*π)*(8/3)*log((1-z_min)/(1-z_max))
    else
        return  α_s/(2*π)*6*log( (z_max*(1-z_min) ) / ((z_min)*(1-z_max)) )
    end
end

function shower(data::Vector)

    nbr_particles_old = 1

    p_new = 0.
    z_new = 0.
    a_prob = 0.

    #for j = 1:2

    simulate = true

    steps = 0

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

                    z_min_tilde = 1/Ea
                    z_max_tilde = 1-1/Ea

                    I = I_a_bc(parent, z_min_tilde, z_max_tilde, Ea)

                    if I < 0
                        println(parent)
                        println(p_old)
                        println(I)
                        println(Ea)
                        println(z_min_tilde)
                        println(z_max_tilde)
                        sample = false
                        simulate = false
                    end

                    p2_new = p_old^2*rand()^(1/I)
                    p_new = sqrt(p2_new)
                    p_old = p_new


                    if p_new <= 1 #|| p_new == Inf
                        z_new = 0.
                        p_new = 0.
                        parent = "dead"
                        sample = false
                    elseif p_new > Ea/2
                        p_old = p_new
                        # prop to state for event
                    else

                        if parent == "q"
                            z_new = 1 - (1 - z_min_tilde)*((1-z_max_tilde)/(1-z_min_tilde))^(rand())
                        else
                            R = rand()
                            aR = ((z_max_tilde*(1-z_min_tilde) ) / ((z_min_tilde)*(1-z_max_tilde)))^R
                            z_new = aR/(1+aR)
                            #a = (z_min_tilde/(1-z_min_tilde))*(z_max_tilde/(1-z_max_tilde))^(R)*(z_min_tilde/(1-z_min_tilde))^(-R)
                            #z_new = 1/(1/a+1)
                        end

                        z_min = p_new/Ea # what should these values be???
                        z_max = 1-p_new/Ea

                        if z_new <= z_max && z_new >= z_min
                            if parent == "q"
                                a_prob = (1 + z_new^2)/2
                            else
                                a_prob = (1 - z_new*(1-z_new))^2/2
                            end
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
                Pb = Decay(parent, "q", Eb, p_new)
                Pc = Decay(parent, "g", Ec, p_new)
            elseif parent == "g"
                Pb = Decay(parent, "g", Eb, p_new)
                Pc = Decay(parent, "g", Ec, p_new)
            else
                Pb = Decay(parent, "dead", 0., 0.)
                Pc = Decay(parent, "dead", 0., 0.)
            end

            append!(data,[Pb])
            append!(data,[Pc])


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
@printf "nbr_q_qg (mena): %.2f\n" mean(stats[2,:])
@printf "nbr_q_qg (std): %.2f\n" std(stats[2,:])

@printf "nbr_g (mena): %.2f\n" mean(stats[1,:])
@printf "nbr_g (std): %.2f\n" std(stats[1,:])
