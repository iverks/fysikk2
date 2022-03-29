import matplotlib.pyplot as plt
import matplotlib.widgets as wdgt
import numpy as np
from scipy import linalg
from scipy.constants import speed_of_light, pi, electron_mass, hbar, Planck

def draw_interactive(xs, V, E, psi, N):
    fig, ax = plt.subplots()
    pot, = plt.plot(xs, V)
    plt.ylim(-1500, 1000)
    schrod, = plt.plot(xs, E[0]+psi[:,0]*1000)
    plt.subplots_adjust(bottom=0.25)
    pickable_states = [i for i in range(10)] + [10*i for i in range(1, 10)] + [25*i for i in range(4, 8+1)]
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
        schrod.set_ydata(E[state]+psi[:,state]*1000)
        fig.canvas.draw_idle()

    state_number.on_changed(update)

    plt.show()
    
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


def part2():
    N = 2000
    x_min = 0.4
    x_max = 9
    dx = (x_max-x_min)/(N-1)

    xs = np.linspace(x_min, x_max, N)

    U_0 = 943 # *10**3 # J/mol
    r_0 = 1.08 # *10**-10 #m
    alpha = 2.73 # *10**10 # m^-1
    # U = lambda x: U_0*(np.exp(-2*alpha*(x-r_0))-2*np.exp(-alpha*(x-r_0)))
    # V = [U(x) for x in xs]
    V = U_0 + U_0*(np.exp(-2*alpha*(xs-r_0))-2*np.exp(-alpha*(xs-r_0)))

    H = np.zeros((N, N))
    for i in range(N):
        H[i][i] = (1/dx**2+V[i])
        if i > 0:
            H[i][i-1] = -1/(2*dx**2)
        if i < N-1:
            H[i][i+1] = -1/(2*dx**2)
    E, psi = linalg.eigh(H)

    find_info(E, 0)

    # draw_interactive(xs, V, E, psi, N)

    # fig, axs = plt.subplots(2, 3)
    # pot, = axs[0, 0].plot(xs, V)
    # axs[0, 0].legend([pot], ["Potential"])
    # for i in range(1, 6):
    #     x = i//3
    #     y = i%3
    #     line, = axs[x, y].plot(xs, psi[:, i])
    #     axs[x, y].legend([line], [f"State {i}"])
    # plt.show()


if __name__=="__main__":
    part2()