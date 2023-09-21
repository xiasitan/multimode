using QuantumOptics
using BenchmarkTools



samples = 1
evals = 100
cutoffs = [50:10:70;]
displacements = [1,2]

name = "ecd_timeevolution_master_1,2_test"

function f(N,α)

    # system parameter
    χ = 2*pi*1e6
    κ = 0

    # Operators
    bc = FockBasis(N-1)
    bq = SpinBasis(1//2)
    b = bq ⊗ bc

    sx = sigmax(bq)
    sy = sigmay(bq)
    sz = sigmaz(bq)

    Iq = identityoperator(bq)
  
    a = destroy(bc)
    at = create(bc)
    
    # Jump Operator

    J = [embed(b, 2, a*sqrt(κ))]
    
    # states
    g = spindown(bq)
    g_vac = g ⊗ fockstate(bc,0)
    
    # Hamiltonians
    H_dispersive= -χ/2 * sz ⊗ (at*a)

    function H_cav_drive(ϵ)
        # sign is chosen to match displace 
        # α = ϵ*T_max
        (conj(ϵ * 1im)a + ϵ * 1im*dagger(a)) 
    end
    
    function H_qubit_drive(δ, θ = 0)
        # θ = 0 -> rotate around x axis
        # θ = pi/2 -> rotate around y axis 
    
        # pi flip if ϵ*T_max = pi
        δ/2*(cos(θ)*sx + sin(θ)*sy)
    end

    function wait_evolution_loss(t, ρ, J)
        tout, ρt = timeevolution.master(t, ρ, H_dispersive, J)
        last(ρt)
    end

    function qubit_rotation_evolution_loss(δ, t, ρ, J, θ = 0)
        H_total = embed(b,1,H_qubit_drive(δ, θ)) + H_dispersive
        tout, ρt = timeevolution.master(t, ρ, H_total, J);
        last(ρt)
    end

    function displacement_evolution_loss(ϵ ,t , ρ, J)
        H_total = Iq ⊗ H_cav_drive(ϵ) + H_dispersive
        tout, ρt = timeevolution.master(t, ρ, H_total, J)
        last(ρt)
    end

    function ECD_evolution_loss(ϵ, t_displace, t_wait, t_pi, ρ0, J, θ = 0)
        # first displacement
        ρ1 = displacement_evolution_loss(ϵ, t_displace, ρ0, J)
        # wait
        ρ2 = wait_evolution_loss(t_wait, ρ1, J)
        # second displacement
        ρ3 = displacement_evolution_loss(-ϵ*cos(χ/2*maximum(t_wait)), t_displace, ρ2, J)
        # pi flip
        ρ4 = qubit_rotation_evolution_loss(π/maximum(t_pi), t_pi, ρ3, J, θ)
        # third displacement
        ρ5 = displacement_evolution_loss(-ϵ*cos(χ/2*maximum(t_wait)), t_displace, ρ4, J)
        # wait 
        ρ6 = wait_evolution_loss(t_wait, ρ5, J)
        # fourth displacement
        ρ7 = displacement_evolution_loss(ϵ*cos(χ*maximum(t_wait)), t_displace, ρ6, J)
        return ρ7
    end


    # time parameters
    t_wait = [0:0.01:3;].*0.05e-6
    t_pi = [0:0.1:3;].*1e-9
    t_displace= [0:1:3;].*1e-9;
    # calculate epsilon such that we implement an ecd(1)
    ϵ = 1/(2*maximum(t_displace)*sin(χ/2*maximum(t_wait)))
    #print("Displacement alpha for ECD(1) is ", ϵ*maximum(t_displace), "\n")

    ECD_evolution_loss(ϵ*α, t_displace, t_wait, t_pi, g_vac, J, π/2) 

end


println("Benchmarking: ", name)

# run benchmark and save data in nested dictionary = { alpha=1 : {N: dimension, {t = time} }}
nested_dict = Dict{Any, Dict{String, Vector{Float32}}}()
for α in displacements
    one_alpha_dict = Dict{String, Vector{Any}}()
    println()
    println("Max displacement: ", α)
    print("Cutoff: ")
    for N in cutoffs
        print(N," ")
        t = @belapsed f($N, $α) samples=samples evals=evals
        if haskey(one_alpha_dict, "N")
            # Key already exists, append the value to the existing list
            push!(one_alpha_dict["N"], N)
        else
            # Key doesn't exist, create a new key with a vector containing the value
            one_alpha_dict["N"] = [N]
        end

        if haskey(one_alpha_dict, "t")
            # Key already exists, append the value to the existing list
            push!(one_alpha_dict["t"], t)
        else
            # Key doesn't exist, create a new key with a vector containing the value
            one_alpha_dict["t"] = [t]
        end
        nested_dict["alpha=$α"] = one_alpha_dict
    end
end
println()

include("bmutils.jl")
bmutils.save(name, nested_dict)


