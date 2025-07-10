import re
import math
import numpy as np
from collections import defaultdict
import pygad

def parse_thermo_file(filename):
    with open(filename, 'r') as f:
        lines = [line.strip('\n') for line in f.readlines()]
    species_data = {}
    i = 0
    while i < len(lines) - 3:
        header = lines[i]
        if len(header) < 80:
            i += 1
            continue
        name = header[:24].strip()
        name = name.split(' ')[0]
        try:
            T_low = float(header[45:55])
            T_mid = float(header[65:75])
            T_high = float(header[55:65])
        except:
            i += 1
            continue
        line2, line3, line4 = lines[i+1], lines[i+2], lines[i+3]
        try:
            coeffs = [
                float(line2[0:15]), float(line2[15:30]), float(line2[30:45]), float(line2[45:60]), float(line2[60:75]),
                float(line3[0:15]), float(line3[15:30]), float(line3[30:45]), float(line3[45:60]), float(line3[60:75]),
                float(line4[0:15]), float(line4[15:30]), float(line4[30:45]), float(line4[45:60])
            ]
        except:
            i += 4
            continue
        species_data[name] = {
            'T_low': T_low,
            'T_mid': T_mid,
            'T_high': T_high,
            'coeffs': coeffs
        }
        i += 4
    return species_data

def get_coefficients(thermo, specie, temperature):
    if temperature <= 1000:
        return thermo[specie]['coeffs'][0:7]
    else:
        return thermo[specie]['coeffs'][7:14]


# function that computes the del cp value using the formula delcp = cp of molecule - sum[cp of elements]*stoichiometry
def compute_del_cp(species_decomposition, specie, i):
    sum = 0
    print("for specie: ",specie)
    for element in species_decomposition[specie]:
        stoich = species_decomposition[specie][element]
        # get cp of element
        coefficient = get_coefficients(thermo, element, T)
        # store sum of cp of elements
        sum = sum + stoich * coefficient[i]
    
    # get cp for molecule
    coefficient = get_coefficients(thermo, specie, T)
    return coefficient[i]-sum

def parse_formula(formula):
    """
    Parses a chemical formula like CO2 into {'C': 1, 'O': 2}.
    """
    pattern = r'([A-Z][a-z]*)(\d*)'
    matches = re.findall(pattern, formula)
    atom_counts = {}
    for atom, count in matches:
        atom_counts[atom] = atom_counts.get(atom, 0) + (int(count) if count else 1)
    return atom_counts


def compute_stoich_and_initial_atoms(species_list, mole_list):
    element_set = set()
    species_atom_counts = []

    # Step 1: Parse each formula
    for formula in species_list:
        atom_dict = parse_formula(formula)
        species_atom_counts.append(atom_dict)
        element_set.update(atom_dict.keys())

    element_list = sorted(element_set)  # for consistent order

    # Step 2: Build stoich
    stoich = {}
    for element in element_list:
        stoich[element] = []
        for atom_dict in species_atom_counts:
            stoich[element].append(atom_dict.get(element, 0))

    # Step 3: Compute initial_atoms using mole_list
    initial_atoms = defaultdict(float)
    for atom_dict, moles in zip(species_atom_counts, mole_list):
        for elem, count in atom_dict.items():
            initial_atoms[elem] += moles * count

    return stoich, dict(initial_atoms)



## MAIN

#Step 1 Reading the thermo30.dat file into a variable of type dictionary called thermo
thermo = parse_thermo_file("thermo30.dat")
# print(type(thermo))
print(thermo)

#Step 1.1 Displaying available species
print("\n📘 Available species:")
thermo_keys = list(thermo.keys())
for i, sp in enumerate(thermo_keys):
    print(f"{i}: {sp}")

#Step 1.2 Taking input for species
sel = input("\nEnter species indices (comma-separated): ")
species_list = [thermo_keys[int(i.strip())] for i in sel.split(',')]
print(species_list)

#Step 3 Taking Temp and Pressure input from User
T = float(input("\nEnter temperature in K: "))
P = float(input("Enter pressure in bar: "))

#Step 3.1 Selecting the coefficients based on temperature

coefficients = {} # a dict to store the coefficients for selected species

# for all species, get their coefficients
for specie in species_list:
    coefficients[specie] = get_coefficients(thermo, specie, T)

#Step 5 Modified - Calculating enthalpy and entropy for each specie
enthalpy = {} # H (T) / RT
entropy = {}    # S(T) / R
gibbs = {} # G(T) / RT

for specie in species_list:
    enthalpy[specie] = (coefficients[specie][0]) + (coefficients[specie][1] * T/2) + (coefficients[specie][2]/3 * T * T) +  (coefficients[specie][3]/4 * T * T * T) + (coefficients[specie][4]/5 * T * T * T * T) + coefficients[specie][5]/T
    entropy[specie] = (coefficients[specie][0] * math.log(T)) + (coefficients[specie][1] * T) + (coefficients[specie][2]/2 * T * T) + (coefficients[specie][3]/3 * T * T * T) + (coefficients[specie][4]/4 * T * T * T * T) + (coefficients[specie][6])
    gibbs[specie] = enthalpy[specie] - entropy[specie]

print("Gibbs free energy: ")
print(gibbs)


#Step 6 # Using GA to minimize nG/RT

# taking initial moles user input
initial_moles = []
print("\nEnter initial moles of each selected species:")
for sp in species_list:
    initial_moles.append(float(input(f"{sp}: ")))

stoich , initial_atoms = compute_stoich_and_initial_atoms(species_list, initial_moles)
print(stoich)
print(initial_atoms)

# converting values of dictionary:gibbs into a list
G_RT = list(gibbs.values())

# function that computes nG/RT
def compute_nG_RT(n_i):
    n_out = np.sum(n_i)
    if n_out == 0:
        return 1e6  # avoid division by zero
    Y_i = n_i / n_out
    ln_term = np.log(np.clip(P * Y_i, 1e-12, None))  # avoid log(0)
    return n_out * np.sum(Y_i * (G_RT + ln_term))


# fitness function for Genetic Algo
def fitness_func(ga_instance, solution, _):
    n_i = np.array(solution)
    
    # Elemental balance penalty
    penalty = 0.0
    for element, coeffs in stoich.items():
        atoms = sum(n_i[i] * coeffs[i] for i in range(len(n_i)))
        penalty += abs(atoms - initial_atoms[element])

    gibbs = compute_nG_RT(n_i)
    return -gibbs - 1e6 * penalty  # PyGAD maximizes


num_species = len(species_list)
gene_space = [{'low': 0.0, 'high': 2.0} for _ in range(num_species)]

ga_instance = pygad.GA(
    num_generations=100,
    num_parents_mating=5,
    sol_per_pop=20,
    num_genes=num_species,
    gene_space=gene_space,
    fitness_func=fitness_func,
    mutation_percent_genes=20,
    crossover_type="single_point",
    mutation_type="random",
    stop_criteria=["saturate_20"]
)

ga_instance.run()

solution, solution_fitness, _ = ga_instance.best_solution()
print("Equilibrium mole distribution:", solution)
print("Minimized nG/RT:", -solution_fitness)


# 14 ,3,5,15



