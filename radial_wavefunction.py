# radial wavefunction of hydrogen atom
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import physical_constants
from scipy.special import genlaguerre, factorial

# Constants
a0 = physical_constants['Bohr radius'][0]  # Bohr radius

def R_nl(n, l, r):
    """
    Radial wavefunction for a hydrogen atom.
    """
    rho = 2 * r / n / a0
    L = genlaguerre(n-l-1, 2*l+1)  # Generalized Laguerre polynomial
    norm = np.sqrt((2/n/a0)**3 * factorial(n-l-1) / 2/n/factorial(n+l))
    return norm * rho**l * np.exp(-rho/2) * L(rho)

# Ask for the principal quantum number
n = int(input("Enter the principal quantum number n: "))

# Generate an array of r values from 0 to 10*a0
r = np.linspace(0, 10*a0, 400)

# Create the plot
plt.figure(figsize=(6, 4))

# Define the legend names
legend_names = ['s', 'p', 'd', 'f'] + list(range(26))[4:]

# Calculate and plot the wavefunction for each orbital quantum number
for l in range(n):
    wavefunction = R_nl(n, l, r)
    plt.plot(r/a0, wavefunction, label=f'n={n}, l={legend_names[l]}')

plt.title(f'Radial wavefunctions R(r) for Hydrogen atom (n={n})')
plt.xlabel('Distance from nucleus (a0)')
plt.ylabel('Radial wavefunction R(r)')
plt.legend()
plt.grid(True)
plt.show()
