# import modules
using Plots
using QuantumOptics

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
    )
    
    W = wigner(state, x.*g, y.*g)
    broadcast(abs, W)
    max = maximum(W) # used to scale the colorbar 
    heatmap(x, y, W, size = size, clim = (-max,max),cmap = cmap, aspect_ratio=:equal, xlabel = xlabel, ylabel = ylabel, title = title, framestyle = :grid)
end


function char_func_point(basis ,state, x,y)
    expect(displace(basis, x+y*1im), state)
end 

function characteristic_function(basis, state, xvec = [-3:0.02:3;], yvec = [-3:0.02:3;], is_real = true)
    
    char_func_grid = char_func_point.(Ref(basis), Ref(state), xvec', yvec)

    if is_real
    heatmap(xvec,yvec,real(char_func_grid))
    else
    heatmap(xvec,yvec,imag(char_func_grid))
end
    
end

