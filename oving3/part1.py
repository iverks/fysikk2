import matplotlib.pyplot as plt
import matplotlib.widgets as wdgt
import numpy as np
from scipy import linalg
from scipy.constants import speed_of_light, pi, electron_mass, hbar, Planck


def draw_interactive(xs, V, E, psi, N):
    fig, ax = plt.subplots()
    pot, = plt.plot(xs, V)
    schrod, = plt.plot(xs, E[0] + psi[:,0]*100)
    plt.subplots_adjust(bottom=0.25)

    ax_state = plt.axes([0.25, 0.15, 0.65, 0.03])
    state_number = wdgt.Slider(
        ax=ax_state,
        label="Wave equation state",
        valmin=0,
        valmax=100,
        valstep=1,
        valinit=0
    )

    ax_U = plt.axes([0.25, 0.1, 0.65, 0.03])
    U_slider = wdgt.Slider(
        ax = ax_U,
        label="Well wall height",
        valmin=0,
        valmax=500,
        valstep=1,
        valinit=200
    )
    ax_a = plt.axes([0.25, 0.05, 0.65, 0.03])
    a_slider = wdgt.Slider(
        ax=ax_a,
        label="Well width",
        valmin=0,
        valmax=100,
        valstep=1,
        valinit=1
    )

    def update_solution(_):
        U = U_slider.val
        a = a_slider.val
        state = state_number.val
        xs, V, E, psi = finite_square_well(N, a, U)
        pot.set_ydata(V)
        schrod.set_ydata(E[state]+psi[:,state]*100)
        fig.canvas.draw_idle()

    state_number.on_changed(update_solution)
    U_slider.on_changed(update_solution)
    a_slider.on_changed(update_solution)
    plt.show()

def finite_square_well(N, a, U):
    x_min = 0-0.5
    x_max = a+0.5
    dx = (x_max-x_min)/(N-1)

    xs = np.linspace(x_min, x_max, N)
    V = [0 if x >= 0 and x < a else U for x in xs]
    H = np.zeros((N, N))
    for i in range(N):
        H[i][i] = (1/dx**2+V[i])
        if i > 0:
            H[i][i-1] = -1/(2*dx**2)
        if i < N-1:
            H[i][i+1] = -1/(2*dx**2)
    E, psi = linalg.eigh(H)
    return xs, V, E, psi

def find_info(E, U_0):
    real_E = lambda E: (E + U_0)*2*electron_mass/hbar**2
    for state in range(6):
        E_n = real_E(E[state])
        # E = h*f => f = E/h
        f_n = E_n/Planck
        wavelength = "Ground state does not have"
        if state != 0:
            grnd_E = real_E(E[0])
            delta_E = E_n - grnd_E
            frequency = delta_E/Planck
            # lambda*f = c => lambda = c/f
            wavelength = speed_of_light/frequency
        print(f"{state = }\n{E_n = } J\n{f_n = } Hz\n{wavelength = } m\n")
    print(E[:6])


def main():
    a = 1
    U = 200 # J*kg/hbar^2
    N = 1000
    xs, V, E, psi = finite_square_well(N, a, U)
    
    find_info(E, 0)

    # draw_interactive(xs, V, E, psi, N)

    # Plot 16 lowest eigenstates
    # print(E[:16])
    # fig, axs = plt.subplots(4, 4)
    # for y in range(4):
    #     for x in range(4):
    #         i = y*4 + x
    #         axs[x, y].set_title(f"{i} state")
    #         axs[x, y].plot(xs, psi[:, i])
    #         axs[x, y].grid()
    # plt.show();

    # pot, = plt.plot(xs, V)
    # grnd, = plt.plot(xs, E[0] + psi[:,0])
    # first, = plt.plot(xs, E[1] + psi[:,1])
    # scnd, = plt.plot(xs, E[2] + psi[:,2])
    # plt.legend([pot, grnd, first, scnd], ["potential", "ground state", "first excited state", "second excited state"])
    # plt.show()

if __name__=="__main__":
    main()