# A file that includes all the important system parameters

module Qcrew
using QuantumOptics
using Plots

########### System Parameters ###########
const N = 500
const χ = 1*2π*1e6

const χ_a = 1*2π*1e6
const χ_b = 1*2π*1e6

κ_cav = 1.0*2π*1e6 
κ_qubit = 0
γ_qubit = 0

bc = FockBasis(N)
bq = SpinBasis(1//2)

sx = sigmax(bq)
sy = sigmay(bq)
sz = sigmaz(bq)
sp = sigmap(bq)
sm = sigmam(bq)
g = spindown(bq)
e = spinup(bq)
Iq = identityoperator(bq)

a = destroy(bc)
at = create(bc)
b = destroy(bc)
bt = create(bc)
Ic = identityoperator(bc)

### states ### 
g = spindown(bq)
e = spinup(bq)
g_vac = g ⊗ fockstate(bc,0)
e_vac = e ⊗ fockstate(bc,0)
ge_vac = normalize(g+e) ⊗ fockstate(bc,0)

########### Jump Operators ###########
function photon_loss_cavity_jump_operator(κ = κ_cav,a = a)
    sqrt(κ)*a
end

function photon_loss_qubit_jump_operator(κ = κ_qubit,sm = sm)
    sqrt(κ)*sm
end

function dephasing_qubit_jump_operator(γ = γ_qubit,sz = sz)
    sqrt(κ)* 0.5*sqrt(γ)*sz
end


########### Ideal Gates ###########
function ECD_operator(β)
    sx ⊗ Ic * (sp ⊗ displace(bc, β/2) + sm ⊗ displace(bc, -β/2))
end

function U(u)
end

function V(V)
    
end
########### Hamiltonians ###########
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

########### 1-Mode ECD ###########
# One cavity mode..

function wait_evolution(t, state)
    tout, ψt = timeevolution.schroedinger(t, state, H_dispersive)
    last(ψt)
end

function wait_evolution_loss(t, ρ, J)
    tout, ρt = timeevolution.master(t, ρ, H_dispersive, J)
    last(ρt)
end


function qubit_rotation_evolution(δ, t,  state, θ = 0)
    H_total = H_qubit_drive(δ, θ) ⊗ Ic + H_dispersive
    tout, ψt = timeevolution.schroedinger(t, state, H_total);
    last(ψt)
end

function qubit_rotation_evolution_loss(δ, t, ρ, J, θ = 0)
    H_total = H_qubit_drive(δ, θ) ⊗ Ic + H_dispersive
    tout, ρt = timeevolution.master(t, ρ, H_total, J);
    last(ρt)
end

function displacement_evolution(ϵ ,t , state)
    H_total = Iq ⊗ H_cav_drive(ϵ) + H_dispersive
    tout, ψt = timeevolution.schroedinger(t, state, H_total)
    last(ψt)
end

function displacement_evolution_loss(ϵ ,t , ρ, J)
    H_total = Iq ⊗ H_cav_drive(ϵ) + H_dispersive
    tout, ρt = timeevolution.master(t, ρ, H_total, J)
    last(ρt)
end

function find_epsilon_1_for_ecd(ϵ_list, t_displace, t_wait, t_pi)
    #= 
    Function that calibrates an ECD 1 by finding the displacement ampltidue ϵ that implements an ECD(2) on vac and returns ϵ/2.
    If the given ϵ is too small, the ECD(2) will not be implemented and the condition n>0.99 will not be met. Then, the function will plot the average photon number in the cavity as a function of ϵ.
    
    =#
    abs_values = Vector{Float64}()
    for ϵ in ϵ_list
        state = ECD_evolution(ϵ, t_displace, t_wait, t_pi, g_vac)
        avg_photon_number = abs(expect(Iq⊗(at*a),state))
        push!(abs_values, avg_photon_number)
        if avg_photon_number > 0.99
            return ϵ/2
        end
    end
    plot(ϵ_list, abs_values)
end



    
function ECD_evolution(ϵ, t_displace, t_wait, t_pi, ψ0, θ = 0)
    # first displacement
    ψ1 = displacement_evolution(ϵ, t_displace, ψ0)
    # wait
    ψ2 = wait_evolution(t_wait, ψ1)
    # second displacement
    ψ3 = displacement_evolution(-ϵ*cos(χ/2*maximum(t_displace)), t_displace, ψ2)
    # pi flip
    ψ4 = qubit_rotation_evolution(π/maximum(t_pi), t_pi, ψ3, θ)
    # third displacement
    ψ5 = displacement_evolution(-ϵ*cos(χ/2*maximum(t_displace)), t_displace, ψ4)
    # wait 
    ψ6 = wait_evolution(t_wait, ψ5)
    # fourth displacement
    ψ7 = displacement_evolution(ϵ*cos(χ*maximum(t_displace)), t_displace, ψ6)
    return ψ7
end

function ECD_evolution_test(ϵ, t_displace, t_wait, t_pi, ψ0, θ = 0)
    # first displacement
    ψ1 = displacement_evolution(ϵ, t_displace, ψ0)
    # wait
    ψ2 = wait_evolution(t_wait, ψ1)
    # second displacement
    ψ3 = displacement_evolution(-ϵ*cos(χ/2*maximum(t_displace)), t_displace, ψ2)
    # pi flip
    ψ4 = qubit_rotation_evolution(π/maximum(t_pi), t_pi, ψ3, θ)
    # third displacement
    ψ5 = displacement_evolution(-ϵ*cos(χ/2*maximum(t_displace)), t_displace, ψ4)
    # wait 
    ψ6 = wait_evolution(t_wait, ψ5)
    # fourth displacement
    ψ7 = displacement_evolution(ϵ*cos(χ*maximum(t_displace)), t_displace, ψ6)
    return ψ1, ψ7
end

function ECD_evolution_loss(ϵ, t_displace, t_wait, t_pi, ρ0, J, θ = 0
    )
    # first displacement
    ρ1 = displacement_evolution_loss(ϵ, t_displace, ρ0, J)
    # wait
    ρ2 = wait_evolution_loss(t_wait, ρ1, J)
    # second displacement
    ρ3 = displacement_evolution_loss(-ϵ*cos(χ/2*maximum(t_displace)), t_displace, ρ2, J)
    # pi flip
    ρ4 = qubit_rotation_evolution_loss(π/maximum(t_pi), t_pi, ρ3, J, θ)
    # third displacement
    ρ5 = displacement_evolution_loss(-ϵ*cos(χ/2*maximum(t_displace)), t_displace, ρ4, J)
    # wait 
    ρ6 = wait_evolution_loss(t_wait, ρ5, J)
    # fourth displacement
    ρ7 = displacement_evolution_loss(ϵ*cos(χ*maximum(t_displace)), t_displace, ρ6, J)
    return ρ7
end

function ECD_evolution_loss_test(ϵ, t_displace, t_wait, t_pi, ρ0, J, θ = 0
    )
    # first displacement
    ρ1 = displacement_evolution_loss(ϵ, t_displace, ρ0, J)
    # wait
    ρ2 = wait_evolution_loss(t_wait, ρ1, J)
    # second displacement
    ρ3 = displacement_evolution_loss(-ϵ*cos(χ/2*maximum(t_displace)), t_displace, ρ2, J)
    # pi flip
    ρ4 = qubit_rotation_evolution_loss(π/maximum(t_pi), t_pi, ρ3, J, θ)
    # third displacement
    ρ5 = displacement_evolution_loss(-ϵ*cos(χ/2*maximum(t_displace)), t_displace, ρ4, J)
    # wait 
    ρ6 = wait_evolution_loss(t_wait, ρ5, J)
    # fourth displacement
    ρ7 = displacement_evolution_loss(ϵ*cos(χ*maximum(t_displace)), t_displace, ρ6, J)
    return ρ1, ρ7
end

function U_evolution_loss(u, ϵ_1, t_displace_1, t_wait_1, t_pi2, ρ0, J)
    # Note, pi2 and pi use same wait time, but amplitude is doubled 
    
    ρ1 = qubit_rotation_evolution_loss(π/2.0*1/maximum(t_pi2), t_pi2, ρ0, J, 3/2.0*π)
    
    ρ2 = ECD_evolution_loss(u*ϵ_1, t_displace_1, t_wait_1, t_pi2, ρ1, J,  1/2*π)

    ρ3 = qubit_rotation_evolution_loss(π/2.0*1/maximum(t_pi2), t_pi2, ρ2, J,  3/2*π)
 
end

function V_evolution_loss(u, ϵ_1, t_displace_1, t_wait_1, t_pi2, ρ0, J )
    
    ρ1 = qubit_rotation_evolution_loss(π/2*1/maximum(t_pi2), t_pi2, ρ0, J)

    ρ2 = ECD_evolution_loss(u*ϵ_1*1im, t_displace_1, t_wait_1, t_pi2, ρ1, J, π)

    ρ3 = qubit_rotation_evolution_loss(π/2*1/maximum(t_pi2), t_pi2, ρ2, J)
    
end

function UVs_evolution_loss(uvs, ϵ_1, t_displace_1, t_wait_1, t_pi2, ρ0, J)
    repetition = Int(length(uvs)/2)
    rhos = Vector{Any}()
    push!(rhos,ρ0)

    for i = 1:repetition
        print(uvs[i])
        push!(rhos,U_evolution_loss(uvs[i], ϵ_1, t_displace_1, t_wait_1, t_pi2, last(rhos), J) )
        print(uvs[i+1])
        push!(rhos ,V_evolution_loss(uvs[i+1], ϵ_1, t_displace_1, t_wait_1, t_pi2, last(rhos), J) )
    end 
    return rhos
end

########### 2-Mode ###########


g_vac_vac = g ⊗ fockstate(bc,0) ⊗ fockstate(bc,0)
e_vac_vac = e ⊗ fockstate(bc,0) ⊗ fockstate(bc,0)
ge_vac_vac = normalize(g+e) ⊗ fockstate(bc,0) ⊗ fockstate(bc,0)

H_dispersive_two_mode= -χ_a/2 * sz ⊗ (at*a) ⊗ Ic + -χ_b/2 * sz ⊗ Ic ⊗ (bt*b)


function wait_two_mode_evolution_loss(t, ρ, J)
    tout, ρt = timeevolution.master(t, ρ, H_dispersive_two_mode, J)
    last(ρt)
end

function qubit_rotation_two_mode_evolution_loss(δ, t, ρ, J, θ = 0)
    H_total = H_qubit_drive(δ, θ) ⊗ Ic ⊗ Ic + H_dispersive_two_mode
    tout, ρt = timeevolution.master(t, ρ, H_total, J);
    last(ρt)
end

function displacement_two_mode_evolution_loss(ϵ ,t , ρ, J, mode_index)
    # mode_index = 1 -> cavity a
    # mode_index = 2 -> cavity b
    # mode_index [1,2] -> both cavities
    if mode_index == 1
        H_total = Iq ⊗ H_cav_drive(ϵ) ⊗ Ic + H_dispersive_two_mode
    elseif mode_index == 2  
        H_total = Iq ⊗ Ic ⊗ H_cav_drive(ϵ) + H_dispersive_two_mode

    elseif mode_index == [1,2]
        if length(ϵ) != 2
            error("ϵ must be a vector of length 2")
        end
        H_total = Iq ⊗ H_cav_drive(ϵ[1]) ⊗ H_cav_drive(ϵ[2]) + H_dispersive_two_mode
    else
        error("mode_index must be 1, 2 or [1,2]")
    end

    tout, ρt = timeevolution.master(t, ρ, H_total, J)
    last(ρt)
end

function ECD_two_mode_evolution_loss(ϵ, t_displace, t_wait, t_pi, ρ0, J, θ = 0
    )
    # first displacement
    ρ1 = displacement_two_mode_evolution_loss(ϵ, t_displace, ρ0, J, [1,2])
    # wait
    ρ2 = wait_two_mode_evolution_loss(t_wait, ρ1, J)
    # second displacement
    ρ3 = displacement_two_mode_evolution_loss(-ϵ.*cos(χ/2*maximum(t_displace)), t_displace, ρ2, J, [1,2])
    # pi flip
    ρ4 =  qubit_rotation_two_mode_evolution_loss(π*1/maximum(t_pi), t_pi, ρ3, J,  θ)
    # third displacement
    ρ5 = displacement_two_mode_evolution_loss(-ϵ.*cos(χ/2*maximum(t_displace)), t_displace, ρ4, J, [1,2])
    # wait 
    ρ6 = wait_two_mode_evolution_loss(t_wait, ρ5, J)
    # fourth displacement
    ρ7 = displacement_two_mode_evolution_loss(ϵ.*cos(χ*maximum(t_displace)), t_displace, ρ6, J, [1,2])
    return ρ7
end

function U_two_mode_evolution_loss(u, ϵ_1, t_displace_1, t_wait_1, t_pi2, ρ0, J)
    # Note, pi2 and pi use same wait time, but amplitude is doubled 
    if length(u) != 2 || length(ϵ_1) != 2
        error("u and e_1 must be a vector of length 2")
    end
    
    ρ1 = qubit_rotation_two_mode_evolution_loss(π/2.0*1/maximum(t_pi2), t_pi2, ρ0, J, 3/2.0*π)
    
    ρ2 = ECD_evolution_two_mode_loss([u[1]*ϵ_1[1],u[2]*ϵ_1[2]], t_displace_1, t_wait_1, t_pi2, ρ1, J,  1/2*π)

    ρ3 = qubit_rotation_two_mode_evolution_loss(π/2.0*1/maximum(t_pi2), t_pi2, ρ2, J,  3/2*π)
 
end

function V_two_mode_evolution_loss(u, ϵ_1, t_displace_1, t_wait_1, t_pi2, ρ0, J )

    if length(u) != 2 || length(ϵ_1) != 2
        error("u and e_1 must be a vector of length 2")
    end
    
    ρ1 = qubit_rotation_two_mode_evolution_loss(π/2*1/maximum(t_pi2), t_pi2, ρ0, J)

    ρ2 = ECD_evolution_two_mode_loss([u[1]*ϵ_1[1],u[2]*ϵ_1[2]], t_displace_1, t_wait_1, t_pi2, ρ1, J, π)

    ρ3 = qubit_rotation_two_mode_evolution_loss(π/2*1/maximum(t_pi2), t_pi2, ρ2, J)
    
end

function UVs_two_mode_evolution_loss(uvs, ϵ_1, t_displace_1, t_wait_1, t_pi2, ρ0, J)
    repetition = Int(length(uvs)/2)
    rhos = Vector{Any}()
    push!(rhos,ρ0)

    for i = 1:repetition
        print(uvs[i])
        push!(rhos,U_evolution_loss(uvs[i], ϵ_1, t_displace_1, t_wait_1, t_pi2, last(rhos), J) )
        print(uvs[i+1])
        push!(rhos ,V_evolution_loss(uvs[i+1], ϵ_1, t_displace_1, t_wait_1, t_pi2, last(rhos), J) )
    end 
    return rhos
end


end





