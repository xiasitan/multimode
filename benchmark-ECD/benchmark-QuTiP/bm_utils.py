import timeit
import json


result_path = "C:/Users/jonat/Desktop/Code/multimode/benchmark-ECD/results/results-QuTiP-{}.json"


def run_benchmark(f, *args, samples=5, evals=1):
    D = {"f": f, "args": args}
    t = timeit.repeat("f(*args)", globals=D, number=evals, repeat=samples)
    return min(t)/evals


def save(name, results):
    f = open(result_path.format(name), "w")
    json.dump(results, f)
    f.close()