import matplotlib.pyplot as plt
import matplotlib.widgets as wdgt
import numpy as np
from scipy import linalg
from scipy.constants import speed_of_light, Planck, Avogadro, electron_mass, proton_mass, electron_volt

def draw_interactive(xs, V, E, psi, N):
    fig, ax = plt.subplots()
    pot, = plt.plot(xs, V)
    # plt.ylim(-20e-14, 15e-13)
    schrod, = plt.plot(xs, E[0]+psi[:,0])
    plt.subplots_adjust(bottom=0.25)
    stateax = plt.axes([0.25, 0.15, 0.65, 0.03])
    state_number = wdgt.Slider(
        ax=stateax,
        label="Wave equation state",
        valmin=0,
        valmax=100,
        valstep=1,
        valinit=0
    )

    def update(_):
        state = state_number.val
        schrod.set_ydata(E[state]+psi[:,state])
        fig.canvas.draw_idle()

    state_number.on_changed(update)

    plt.show()

def find_info(E):
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


def part2():
    # Constants not in scipy
    hartree = 4.3597447222071*10**-18 # Unit: Joules/hartree
    bohr = 5.29177210903*10**-11 # Unit: metres/bohr radii

    # Constants from the material
    U_0 = 943*10**3/Avogadro/hartree # current unit: hartree
    r_0 = 1.08*10**-10/bohr # current unit: bohr radii
    alpha = 2.73*10**10*bohr # current unit: bohr radii
    m_nitrogen = 14 # Unit: proton mass
    m_reduced = m_nitrogen/2 #m*m/(m+m) 
    m = m_reduced*proton_mass/electron_mass # Converted to electron masses (hartree units)

    # Borders and resolution
    N = 2000
    x_min = 7.5*r_0/10 
    x_max = r_0*2.5
    dx = (x_max-x_min)/(N-1)

    xs = np.linspace(x_min, x_max, N)

    V = U_0 + U_0 *(np.exp(-2*alpha*(xs-r_0))-2*np.exp(-alpha*(xs-r_0))) # current unit: kJ/mol

    H = np.zeros((N, N))
    for i in range(N):
        H[i][i] = 1/(m*dx**2)+V[i]
        if i > 0:
            H[i][i-1] = -1/(m*2*dx**2)
        if i < N-1:
            H[i][i+1] = -1/(m*2*dx**2)
    E, psi = linalg.eigh(H)

     # Converting to SI units
    E = E * hartree
    xs = xs * bohr
    V = V * hartree
    psi = psi * hartree  # Done only for visual reasons

    find_info(E)

    # Converting to nice plotting units
    E = E/electron_volt
    xs = xs/10e-10
    V = V/electron_volt
    psi = psi/electron_volt/20

    # draw_interactive(xs, V, E, psi, N)

    # fig, axs = plt.subplots(2, 3)
    # pot, = axs[0, 0].plot(xs, V)
    # axs[0, 0].legend([pot], ["Potential"])
    # for i in range(1, 6):
    #     x = i//3
    #     y = i%3
    #     line, = axs[x, y].plot(xs, psi[:, i-1])
    #     axs[x, y].legend([line], [f"State {i-1}"])
    #     plt.xlabel("metres")
    #     plt.ylabel("Joules")
    # plt.show()

    pot, = plt.plot(xs, V)
    ax = plt.gca()
    legend_funcs = [pot]
    legend_texts = ["Potential"]
    for i in range(1, 6):
        x = i//3
        y = i%3
        color = next(ax._get_lines.prop_cycler)['color']
        line, = plt.plot(xs, E[i] + psi[:, i-1], color=color)
        line, = plt.plot(xs, E[i] + np.zeros_like(xs), color=color, linestyle="--")
        legend_funcs.append(line)
        legend_texts.append(f"State {i-1}")
    plt.legend(legend_funcs, legend_texts)
    plt.xlabel("Ångstrom")
    plt.ylabel("Electron volts")
    plt.show()



if __name__=="__main__":
    part2()