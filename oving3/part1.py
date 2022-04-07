import matplotlib.pyplot as plt
import matplotlib.widgets as wdgt
import numpy as np
from scipy import linalg
from scipy.constants import speed_of_light, pi, electron_volt, hbar, Planck

def finite_square_well(N, a, U):
    x_min = 0-(4*a/10)
    x_max = a+(4*a/10)
    dx = (x_max-x_min)/(N-1)

    xs = np.linspace(x_min, x_max, N)
    V = np.array([0 if x >= 0 and x < a else U for x in xs])
    H = np.zeros((N, N))
    for i in range(N):
        H[i][i] = (1/dx**2+V[i])
        if i > 0:
            H[i][i-1] = -1/(2*dx**2)
        if i < N-1:
            H[i][i+1] = -1/(2*dx**2)
    E, psi = linalg.eigh(H)
    return xs, V, E, psi

def find_info(E):
    # in here we assume SI units
    for state in range(6):
        E_n = E[state]
        # E = h*f => f = E/h
        f_n = E_n/Planck
        wavelength = "Ground state does not have"
        if state != 0:
            grnd_E = E[0]
            delta_E = E_n - grnd_E
            frequency = delta_E/Planck
            # lambda*f = c => lambda = c/f
            wavelength = speed_of_light/frequency
        print(f"{state = }\n{E_n = } J\n{f_n = } Hz\n{wavelength = } m\n")
    print(E[:6])


def main():
   # Constants not in scipy
    hartree = 4.3597447222071*10**-18 # Unit: Joules/hartree
    bohr = 5.29177210903*10**-11 # Unit: metres/bohr radii

    # Setting properties of well
    a = 0.5 # Ångstrom
    U = 40 # eV

    # Converting to hartree
    a = a * 10e-10/bohr # bohr radii
    U = U * electron_volt/hartree # hartree of U is Hartree = 4.3597447222071(85)×10−18 J
    
    N = 1000 

    xs, V, E, psi = finite_square_well(N, a, U)
    
    # Converting to SI units
    hartree = 4.3597447222071*10**-18
    bohr = 5.29177210903*10**-11
    E = E * hartree
    xs = xs * bohr
    V = V * hartree
    U = U * hartree
    psi = psi * hartree
    find_info(E)

    print(E[E<U])

    # Converting to nice plotting units
    E = E/electron_volt
    xs = xs/10e-10
    V = V/electron_volt
    psi = psi/electron_volt * 4

    # Plot 16 lowest eigenstates
    # fig, axs = plt.subplots(4, 4, sharex=True, sharey=True)
    # plt.xlabel("metres")
    # plt.ylabel("Joules")
    # for y in range(4):
    #     for x in range(4):
    #         i = y*4 + x
    #         axs[x, y].set_title(f"{i} state")
    #         axs[x, y].plot(xs, psi[:, i])
    #         axs[x, y].grid()
    # plt.subplots_adjust(wspace=0.4, hspace=0.4)
    # plt.show();

    #plot nicely
    pot, = plt.plot(xs, V)
    ax = plt.gca()
    legend_funcs = [pot]
    legend_texts = ["Potential"]
    for i in range(1, 6):
        x = i//3
        y = i%3
        color = next(ax._get_lines.prop_cycler)['color']
        line, = plt.plot(xs, E[i] + psi[:, i-1]/3, color=color)
        plt.plot(xs, E[i] + np.zeros_like(xs), color=color, linestyle="--")
        legend_funcs.append(line)
        legend_texts.append(f"State {i-1}")
    plt.legend(legend_funcs, legend_texts)
    plt.xlabel("Ångstrom")
    plt.ylabel("Electron volts")
    plt.show()

if __name__=="__main__":
    main()