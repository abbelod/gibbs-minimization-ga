# Chemical Equilibrium Minimization using Genetic Algorithm (PyGAD)

This project calculates the equilibrium composition of a chemical system by minimizing the **non-dimensional Gibbs free energy (nG/RT)** using a **Genetic Algorithm (GA)**. It uses NASA thermodynamic polynomial data and standard thermodynamic values extracted from a PDF.

---

## Features

-  Parses NASA `thermo30.dat` file for temperature-dependent Cp coefficients.
-  Reads standard Gibbs free energy (ΔGf) and enthalpy (ΔHf) values from a PDF file.
-  Calculates G(T)/RT for each species using NASA polynomials.
-  Accepts user input for selected species and initial mole amounts.
-  Enforces **elemental mass balance** using penalty-based constraints.
-  Minimizes total `nG/RT = n_out * ∑[Yi * (G/RT + ln(P * Yi))]` using [PyGAD](https://github.com/ahmedfgad/GeneticAlgorithmPython).
-  Displays equilibrium mole distribution and minimized Gibbs energy.

---

## Formula Used

The Gibbs free energy (non-dimensional) is minimized using:

nG/RT = n_out × Σ [ Yi × (Gi/RT + ln(P × Yi)) ]

Where:
- Yi = ni / n_out  → mole fraction of species i  
- Gi/RT is calculated using NASA polynomials  
- P is the pressure in bar
---

## Requirements

- Python 3.13
- `pygad`
- `numpy`

Install dependencies:

```bash
pip install pygad numpy
