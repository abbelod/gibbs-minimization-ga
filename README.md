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

\[
\frac{nG}{RT} = n_{\text{out}} \cdot \sum_i \left( Y_i \cdot \left(\frac{G_i}{RT} + \ln(P \cdot Y_i)\right) \right)
\]

Where:
- \( Y_i = \frac{n_i}{n_{\text{out}}} \) (mole fraction)
- \( G_i/RT \) is computed from NASA coefficients
- \( P \) is pressure in bar

---

## Requirements

- Python 3.13
- `pygad`
- `numpy`

Install dependencies:

```bash
pip install pygad numpy
