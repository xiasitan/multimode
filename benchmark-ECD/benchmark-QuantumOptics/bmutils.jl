module bmutils
using QuantumOptics
using JSON

# numerical thingy to set all subnormal numbers to zero
set_zero_subnormals(true)

function save(name, results)
    result_path = "C:/Users/jonat/desktop/Code/multimode/benchmark-ECD/results/julia-$name.json"
    f = open(result_path, "w")
    write(f, JSON.json(results))
    close(f)
end

end # module
