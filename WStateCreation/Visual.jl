# import modules
using Plots
using QuantumOptics
using PyPlot



# Plotting Wigner Function
function quick_plot_wigner(state; 
    x = [-3:0.02:3;], 
    y = [-3:0.02:3;],
    cmap = :diverging_bwr_55_98_c37_n256,
    xlabel = "x", 
    ylabel = "y", 
    title = "Wigner function",
    size = (500, 500),
    g = sqrt(2), # such that α = Y + im X -> yes, kinda weird but thats how the displacement gate is defined in QuantumOptics
    ax = nothing,
    cbar=false)

    if ax isa typeof(nothing)
        ax = subplot()
    end
    
    W = wigner(state, x.*g, y.*g)
    
    cf = ax.contourf(x,y,W)
    ax.set_aspect("equal")
    ax.set_xlabel("Real")
    ax.set_ylabel("Im")
    ax.set_title("Wigner")
    if cbar
        colorbar(cf)
    end
end

# function quick_plot_wigner(state; 
#     x = [-3:0.02:3;], 
#     y = [-3:0.02:3;],
#     cmap = :diverging_bwr_55_98_c37_n256,
#     xlabel = "x", 
#     ylabel = "y", 
#     title = "Wigner function",
#     size = (500, 500),
#     g = sqrt(2), # such that α = Y + im X -> yes, kinda weird but thats how the displacement gate is defined in QuantumOptics
#     ax = nothing)
    
#     W = wigner(state, x.*g, y.*g)
#     broadcast(abs, W)
#     max = maximum(W) # used to scale the colorbar 
#     heatmap(x, y, W, size = size, clim = (-max,max),cmap = cmap, aspect_ratio=:equal, xlabel = xlabel, ylabel = ylabel, title = title, framestyle = :grid)
# end


#### CHAR FUNC #### 
function two_mode_char_func(state, α_list, β_list)
    # find use basis to 
    basis_a = basis(ptrace(state,1))
    basis_b = basis(ptrace(state,2))
    char_func_grid = Array{ComplexF64}(undef, length(α_list),length(β_list))

    for (i,α) in enumerate(α_list)
        for (j, β) in enumerate(β_list)
            char_func_grid[i,j] = expect(displace(basis_a,α) ⊗ displace(basis_b,β), state)
        end
    end

    # plot
    ax = subplot()
    cf = ax.contourf(α_list,β_list,real(char_func_grid))
    ax.set_aspect("equal")
    ax.set_xlabel("α")
    ax.set_ylabel("β")
    ax.set_title("Two-Mode Characteristic Function")
    colorbar(cf)
    gcf()
end



function char_func_point(basis ,state, x,y)
    expect(displace(basis, x+y*1im), state)
end 

function char_func(state, xvec = [-3:0.02:3;], yvec = [-3:0.02:3;], is_real = true)
    
    fb = basis(state)
    char_func_grid = char_func_point.(Ref(fb), Ref(state), xvec', yvec)

    # plot
    ax = subplot()
    if is_real
        cf = ax.contourf(xvec,yvec,real(char_func_grid))
    else 
        cf = ax.contourf(xvec,yvec,imag(char_func_grid))
    end

    ax.set_aspect("equal")
    ax.set_xlabel("α")
    ax.set_ylabel("β")
    ax.set_title("SingleMode Characteristic Function")
    colorbar(cf)
    gcf()
end

############ Retired ###############
# function char_func_point(basis ,state, x,y)
#     expect(displace(basis, x+y*1im), state)
# end 

# function characteristic_function(basis, state, xvec = [-3:0.02:3;], yvec = [-3:0.02:3;], is_real = true)
    
#     char_func_grid = char_func_point.(Ref(basis), Ref(state), xvec', yvec)

#     if is_real
#     heatmap(xvec,yvec,real(char_func_grid))
#     else
#     heatmap(xvec,yvec,imag(char_func_grid))
# end
    
# end

