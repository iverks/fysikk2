import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import  ArtistAnimation

from IPython.display import HTML
FPS=30
Interval=50 #Delay between frames in milliseconds 


MINX = -10
MAXX = 10 # [m]
MINY = -10
MAXY = 10

x = np.linspace(MINX, MAXX, 500)

#system parameters
omega1=3*np.pi  
omega2=1.5*np.pi  
v = 1 #m/s

# Time-dependent surce amplitude
def A(t0):
     Amplitude=np.ones(t0.shape)*(1*t0/10)
     Amplitude=abs(np.ones(t0.shape)*np.cos(0.1*np.pi*t0))
     #Amplitude*=(Amplitude>0) 
     return Amplitude

def get_Lambda(omega, v):
    T1 = 2*np.pi/omega
    return v*T1

def y(x,t0,phi,omega,v):
    # calculated form given angular frequency and wave propagation velocity
    Lambda = get_Lambda(omega, v)
    k = 2*np.pi/Lambda
    y = A(t0-np.abs(x)/v)*np.sin(omega*t0 - np.abs(k*x) + phi) 
    return y  

fig, ax = plt.subplots()

xs = np.linspace(MINX, MAXX, 100)
ys = np.linspace(MINY, MAXY, 100).reshape(-1, 1)
xs, ys = np.meshgrid(xs, ys)
r = np.sqrt(xs**2 + ys**2)
plt.imshow(r)
plt.show

imgs = []
for frame_num in range(500):
    t=frame_num*Interval/1000 #time in [s]
    im = ax.imshow(y(r, t, 0, omega1, v), animated=True)
    imgs.append([im])

anim2 = ArtistAnimation(fig, imgs, interval=50, blit = True, repeat_delay=1000)
plt.close(anim2._fig)
plt.show()
HTML(anim2.to_jshtml())