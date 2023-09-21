import numpy as np
import qutip as qt
import timeit, functools
import json

samples = 10
evals = 1
cutoffs = range(50, 71, 10)
displacements = [1, 2]

name = "qutip_ecd_timeevolution_master_1,2_test"

def setup(N):
    options = qt.Options()
    options.nsteps = 1000000
    options.atol = 1e-8
    options.rtol = 1e-6
    return options


def f(N, alpha):
    # System parameters
    chi = 2 * np.pi * 1e6
    kappa = 0

    Iq = qt.qeye(2)
    sx = qt.sigmax()
    sy = qt.sigmay()
    sz = qt.sigmaz()

    Ic = qt.qeye(N)
    a = qt.destroy(N)
    at = a.dag()

    # Jump Operator
    J = [np.sqrt(kappa) * qt.tensor(Iq,a)]

    # States
    g = qt.fock(2, 0)  # Spin-down state
    g_vac = qt.tensor(g, qt.fock(N,0))  # g ⊗ |0⟩

    # Hamiltonians
    H_dispersive = -chi / 2 * qt.tensor(sz, at * a)

    def H_cav_drive(epsilon):
        return np.conj(epsilon * 1j)* a + epsilon * 1j * a.dag()


    def H_qubit_drive(delta, theta=0):
        return delta / 2 * (np.cos(theta) * sx + np.sin(theta) * sy)

    def wait_evolution_loss(t, rho, J):
        H_total = H_dispersive
        result = qt.mesolve(H_total, rho, t, c_ops=J)
        return result.states[-1]

    def qubit_rotation_evolution_loss(delta, t, rho, J, theta=0):
        H_total = qt.tensor(H_qubit_drive(delta, theta),Ic) + H_dispersive
        result = qt.mesolve(H_total, rho, t, c_ops=J)
        return result.states[-1]

    def displacement_evolution_loss(epsilon, t, rho, J):
        H_total = qt.tensor(Iq, H_cav_drive(epsilon)) + H_dispersive
        result = qt.mesolve(H_total, rho, t, c_ops=J)
        return result.states[-1]

    def ECD_evolution_loss(epsilon, t_displace, t_wait, t_pi, rho0, J, theta=0):
        rho1 = displacement_evolution_loss(epsilon, t_displace, rho0, J)
        rho2 = wait_evolution_loss(t_wait, rho1, J)
        rho3 = displacement_evolution_loss(-epsilon * np.cos(chi / 2 * max(t_wait)), t_displace, rho2, J)
        rho4 = qubit_rotation_evolution_loss(np.pi / max(t_pi), t_pi, rho3, J, theta)
        rho5 = displacement_evolution_loss(-epsilon * np.cos(chi / 2 * max(t_wait)), t_displace, rho4, J)
        rho6 = wait_evolution_loss(t_wait, rho5, J)
        rho7 = displacement_evolution_loss(epsilon * np.cos(chi * max(t_wait)), t_displace, rho6, J)
        return rho1, rho7

    # Time parameters
    t_wait = np.arange(0, 3 * 0.05e-6, 0.01 * 0.05e-6)
    t_pi = np.arange(0, 3 * 1e-9, 0.1 * 1e-9)
    t_displace = np.arange(0, 3 * 1e-9, 1 * 1e-9)

    # Calculate epsilon such that we implement an ecd(1)
    epsilon = 1 / (2 * max(t_displace) * np.sin(chi / 2 * max(t_wait)))

    return ECD_evolution_loss(epsilon * alpha, t_displace, t_wait, t_pi, g_vac, J, np.pi / 2)

print("Benchmarking:", name)
print("Cutoff: ", end="", flush=True)

# Run benchmark and save data in a nested dictionary = {alpha=1: {N: dimension, {t = time}}}

def run_benchmark(f, N, alpha, samples=5, evals=1):
    t = timeit.repeat(functools.partial(f, N, alpha),  number=evals, repeat=samples)
    return min(t)/evals

nested_dict = {}

for alpha in displacements:
    one_alpha_dict = {"N": [], "t": []}
    print()
    print("Max displacement:", alpha)
    print("Cutoff:", end=" ")
    
    for N in cutoffs:
        print(N, end=" ")
        t = run_benchmark(f,N,alpha, samples=samples, evals=evals)
        #t = timeit.repeat(functools.partial(f, N, alpha), number=evals, repeat=samples)
        one_alpha_dict["N"].append(N)
        one_alpha_dict["t"].append(t)
        nested_dict[f"alpha={alpha}"] = one_alpha_dict

    print()
# Save the data (you may need to adjust this part based on your specific requirements)
import json

Path = 'C:/Users/jonat/Desktop/Code/multimode/benchmark-ECD/results/'
full_path = Path + name +".json"
with open(full_path, "w") as json_file:
    json.dump(nested_dict, json_file)
    
print(t)
    
