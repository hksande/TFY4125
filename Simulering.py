import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from scipy.interpolate import CubicSpline

y_koordinater = np.asarray([0.251, 0.184, 0.171, 0.187, 0.149, 0.09,  0.077, 0.128])
x_koordinater = np.asarray([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])

# Kjører ring og ball.
# Masse ring: 0.0129 kg (metallring)
# Masse ball: 0.0237 kg
# Radius ytre ring: 0.0200 m
# Indre radius ring: 0.0180 m
# Radius ball: 0.0165 m

ytre_radius_ring = 0.020
indre_radius_ring = 0.018

c_ring = (1 + (indre_radius_ring**2) / ytre_radius_ring) / 2
c_ball = 2/5
g = 9.81
m_ring = 0.0129
m_ball = 0.0237

cs = CubicSpline(x_koordinater, y_koordinater, bc_type='natural')
xmin = 0.000
xmax = 1.401
dx = 0.001
x = np.arange(xmin, xmax, dx)
Nx = len(x)
y0 = 0.251
y = cs(x)
dy = cs(x, 1)
d2y = cs(x, 2)


def v(y0, y, g, c):
    return np.sqrt((2*g*(y0 - y))/(1 + c))


def v_horisontal(fart, beta):
    return fart * np.cos(beta)


def v_gjen_horisontal(fart_horisontal):
    arr = []
    for i in range(len(fart_horisontal) - 1):
        steg = 0.5 * (fart_horisontal[i] + fart_horisontal[i + 1])
        arr.append(steg)

    return np.asarray(arr)


def tid(fart_gjen_horisontal, dx):  # Finner strekning som funksjon av tiden.
    arr = dx / fart_gjen_horisontal
    t = 0
    tider = [0]
    for elem in arr:
        t += elem
        tider.append(t)
    return tider

def helningsvinkel(y_derivert):
    return np.arctan(y_derivert)


def krumning(dy, d2y):
    return d2y / ((1 + dy**2)**(3/2))


def sentripetal_a(g, y0, y, c, dy, d2y):
    return ((2*g*(y0-y))/(1+c))*((d2y)/(1+dy**2)**(3/2))


def normalkraft(M, g, beta, sentripetal_a):
    return M * (g * np.cos(beta) + sentripetal_a)


def friksjon(M, c, g, beta):
    return (c * M * g * np.sin(beta)) / (1 + c)


def friksjonskoeff(f, N):
    return abs(f/N)


def plot(xx, yy, y_min, y_max, navn, y_akse, x_akse='$x$ (m)'):
    plt.figure('y(x)', figsize=(12, 6))
    plt.tick_params(labelsize=14)
    plt.plot(xx, yy, '*', markersize=0.6, linestyle="solid")
    plt.title(navn, size=16)
    plt.xlabel(x_akse, fontsize=20)
    plt.ylabel(y_akse, fontsize=20)
    plt.ylim(y_min, y_max)
    plt.grid()
    plt.show()

def plotbegge(xx, yy, aa, bb, label1, label2, y_min, y_max, navn, y_akse, x_akse='$x$ (m)'):
    plt.figure('y(x)', figsize=(12, 6))
    plt.tick_params(labelsize=14)
    plt.plot(xx, yy, label='label1')
    plt.title(navn, size=16)
    plt.xlabel(x_akse, fontsize=20)
    plt.ylabel(y_akse, fontsize=20)
    plt.ylim(y_min, y_max)
    blue_patch = mpatches.Patch(color='blue', label=label1)
    red_patch = mpatches.Patch(color='red', label=label2)
    plt.legend(handles=[red_patch, blue_patch])
    plt.plot(aa, bb, label='label2')
    plt.xlabel(x_akse, fontsize=20)
    plt.ylabel(y_akse, fontsize=20)
    plt.ylim(y_min, y_max)
    plt.grid()
    #plt.savefig('figure_name')
    plt.show()


# Regn ut helningsvinkel og krumning.
beta = helningsvinkel(dy)
krumning = krumning(dy, d2y)

# Plot bane, helningsvinkel og krumning.
#plot(x, y, 0.0, 0.3, "Baneform", '$y(x)$ (m)')
#plot(x, np.rad2deg(beta), -22.0, 22.0, "Beta", "Grader")
#plot(x, krumning, -5.0, 3.0, "Krumning", "$K(x)$ (1/m)")

# Regn ut verdier for ring
ring_fart = v(y0, y, g, c_ring)
ring_friksjon = friksjon(m_ring, c_ring, g, beta)
ring_normalkraft = normalkraft(m_ring, g, beta, sentripetal_a(g, y0, y, c_ring, dy, d2y))
ring_friksjonskoeff = friksjonskoeff(ring_friksjon, ring_normalkraft)
ring_fart_horisontal = v_horisontal(ring_fart, beta)
ring_fart_gjen_horisontal = v_gjen_horisontal(ring_fart_horisontal)
print("Sluttfart ring: ", ring_fart_horisontal[-1])
ring_tid = tid(ring_fart_gjen_horisontal, dx)

# Regn ut verdier for ball
ball_fart = v(y0, y, g, c_ball)
ball_friksjon = friksjon(m_ball, c_ball, g, beta)
ball_normalkraft = normalkraft(m_ball, g, beta, sentripetal_a(g, y0, y, c_ball, dy, d2y))
ball_friksjonskoeff = friksjonskoeff(ball_friksjon, ball_normalkraft)
ball_fart_horisontal = v_horisontal(ball_fart, beta)
ball_fart_gjen_horisontal = v_gjen_horisontal(ball_fart_horisontal)
print("Sluttfart ball: ", ball_fart_horisontal[-1])
ball_tid = tid(ball_fart_gjen_horisontal, dx)

"""
# Plot verdier for ring
plot(x, ring_fart, 0.0, 2.0, "[Ring] Fart", '$v$ (m/s)')
plot(x, ring_normalkraft, 0.1, 0.4, "[Ring] Normalkraft", "$N$ (N)")
plot(x, ring_friksjonskoeff, -0.1, 0.25, "[Ring] |$f$/N|", "|$f$/N|")
plot(ring_tid, x, 0.0, 1.5, "[Ring] Strekning som funksjon av tid", '$x$ (m)', 't (s)')
plot(ring_tid, ring_fart, 0.0, 2.0, "[Ring] Fart som funksjon av tid", '$v$ (m/s)', 't (s)')

# Plot verdier for ball
plot(x, ball_fart, 0.0, 2.0, "[Kule] Fart", '$v$ (m/s)')
plot(x, ball_normalkraft, 0.1, 0.4, "[Kule] Normalkraft", "$N$ (N)")
plot(x, ball_friksjonskoeff, -0.1, 0.25, "[Kule] |$f$/N|", "|$f$/N|")
plot(ball_tid, x, 0.0, 1.5, "[Kule] Strekning som funksjon av tid", '$x$ (m)', 't (s)')
plot(ball_tid, ball_fart, 0.0, 2.0, "[Kule] Fart som funksjon av tid", '$v$ (m/s)', 't (s)')
"""

#plotbegge(x, ring_fart, x, ball_fart, "Ring", "Kule", 0.0, 2.0, "Fart", '$v(x)$ (m/s)')
#plotbegge(x, ring_normalkraft, x, ball_normalkraft, "Ring", "Kule", 0.1, 0.4, "Normalkraft", '$N$ (N)')
#plotbegge(x, ring_friksjonskoeff, x, ball_friksjonskoeff, "Ring", "Kule", -0.1, 0.25,"Forhold mellom friksjonskraften og normalkraften", '|$f$/N|')
plotbegge(ring_tid, x, ball_tid, x, "Ring", "Kule", 0.0, 1.5, "Strekning som funksjon av tid", '$x$ (m)', 't (s)')
plotbegge(ring_tid, ring_fart, ball_tid, ball_fart, "Ring", "Kule", 0.0, 2.0, "Fart som funksjon av tid",'$v$ (m/s)', 't (s)')
