import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib import animation
import scipy.linalg as la
plt.rcParams["figure.figsize"] = 6.4, 4.8

# Define figure and axis
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Define coordinates and connectivity for adenosine
molecular_coordinates = {
    "Atom": ['N(1)', 'C(2)', 'N(3)', 'C(4)', 'C(5)', 'C(6)', 'N(7)', 'C(8)', 'N(9)', 'N(10)', 'C(11)', 'C(12)', 'C(13)',
             'O(14)', 'C(15)', 'O(16)', 'O(17)', 'H(18)', 'H(19)', 'H(20)', 'C(21)', 'H(22)', 'P(23)', 'O(24)', 'O(25)',
             'O(26)', 'O(27)', 'H(28)', 'H(29)', 'H(30)', 'H(31)', 'H(32)', 'H(33)', 'H(34)', 'H(35)'],
    "x": [3.82, 2.498, 1.417, 1.688, 3.005, 4.081, 0.735, 1.697, 2.982, 5.269, -2.57, -1.256, -0.735, -1.229, -2.234,
          -1.443, -3.577, -0.564, -2.963, -1.074, -1.699, -3.136, -5.038, -4.985, -5.452, -1.392, -6.014, 2.314, 1.378,
          6.076, 5.449, -2.027, -2.471, -0.78, -1.066],
    "y": [-2.917, -3.397, -2.513, -1.155, -0.675, -1.544, 0.0, 1.165, 0.777, -1.108, 0.936, 0.252, 0.0, 1.066, 1.779,
          -0.949, 0.005, 0.848, 1.521, -0.989, 3.136, 1.957, 0.617, 2.093, 0.224, 3.906, 0.095, -4.482, 2.218, -1.78,
          -0.073, -1.478, 3.66, 2.992, 4.733],
    "z": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.121, 1.398, 0.0, -0.824, -0.089, 2.097, 0.83, 2.035,
          1.982, -0.381, 0.349, -0.716, 0.514, 0.607, -0.851, -0.781, 1.497, 0.0, 0.0, 0.0, 0.0, 1.58, 0.956, 0.961,
          -0.469]
}

molecular_coordinates_df = pd.DataFrame(molecular_coordinates)

molecular_connectivity = {
    "Atom Name": ['C(2)', 'N(3)', 'H(28)', 'C(4)', 'C(5)', 'N(7)', 'C(6)', 'N(9)', 'N(10)', 'C(8)', 'C(13)', 'H(29)',
                  'C(12)', 'O(14)', 'H(20)', 'C(11)', 'O(16)', 'H(18)', 'C(15)', 'O(17)', 'C(21)', 'H(19)', 'H(22)',
                  'P(23)', 'O(26)', 'H(33)', 'H(34)', 'O(24)', 'O(25)', 'O(27)', 'H(30)', 'H(31)', 'H(32)', 'H(35)',
                  'C(5)', 'N(9)', 'C(11)'],
    "Connectivity": ['N(1)', 'C(2)', 'C(2)', 'N(3)', 'C(4)', 'C(4)', 'N(1)', 'C(5)', 'C(6)', 'N(7)', 'N(7)', 'C(8)',
                     'C(13)', 'C(13)', 'C(13)', 'C(12)', 'C(12)', 'C(12)', 'O(14)', 'C(11)', 'C(15)', 'C(11)', 'C(15)',
                     'O(17)', 'C(21)', 'C(21)', 'C(21)', 'P(23)', 'P(23)', 'P(23)', 'N(10)', 'N(10)', 'O(16)', 'O(26)',
                     'C(6)', 'C(8)', 'C(15)']
}

molecular_connectivity_df = pd.DataFrame(molecular_connectivity)

# Animate
def updatefig(i):
    ax.clear()
    # Rotate molecular fragment around a bond
    n = 0
    while n <= len(molecular_coordinates_df["Atom"]) - 1:
        if molecular_coordinates_df["Atom"][n] in ['N(1)', 'C(2)', 'N(3)', 'C(4)', 'C(5)', 'C(6)', 'N(7)', 'C(8)',
                                                   'N(9)', 'N(10)', 'H(28)', 'H(29)', 'H(30)', 'H(31)']:
            theta = np.radians(i*5)
            A = np.array([[1, 0, 0], [0, np.cos(theta), -np.sin(theta)], [0, np.sin(theta), np.cos(theta)]])
            B = np.array([[float(molecular_coordinates_df["x"][n])], [float(molecular_coordinates_df["y"][n])],
                          [float(molecular_coordinates_df["z"][n])]])
            C = A @ B
            molecular_coordinates_df["x"][n] = C.item(0)
            molecular_coordinates_df["y"][n] = C.item(1)
            molecular_coordinates_df["z"][n] = C.item(2)
        n += 1

    # Assign colors
    n = 0
    list = []
    while n <= len(molecular_coordinates_df["Atom"]) - 1:
        if molecular_coordinates_df["Atom"][n][0] == "C":
            list.append('0.3')
        elif molecular_coordinates_df["Atom"][n][0] == "O":
            list.append('r')
        elif molecular_coordinates_df["Atom"][n][0] == "N":
            list.append('b')
        elif molecular_coordinates_df["Atom"][n][0] == "H":
            list.append('0.7')
        elif molecular_coordinates_df["Atom"][n][0] == "O":
            list.append('r')
        elif molecular_coordinates_df["Atom"][n][0] == "P":
            list.append('m')
        n += 1
    molecular_coordinates_df["color"] = list

    # Draw bonds
    xs = []
    ys = []
    zs = []
    n = 0
    m = 0
    while n <= len(molecular_connectivity_df) - 1:
        while m <= len(molecular_coordinates_df) - 1:
            if molecular_connectivity_df["Atom Name"][n] == molecular_coordinates_df["Atom"][m]:
                xs.append(molecular_coordinates_df["x"][m])
                ys.append(molecular_coordinates_df["y"][m])
                zs.append(molecular_coordinates_df["z"][m])
            elif molecular_connectivity_df["Connectivity"][n] == molecular_coordinates_df["Atom"][m]:
                xs.append(molecular_coordinates_df["x"][m])
                ys.append(molecular_coordinates_df["y"][m])
                zs.append(molecular_coordinates_df["z"][m])
            m += 1
        m = 0
        n += 1

    xs_list = [xs[x:x + 2] for x in range(0, len(xs), 2)]
    ys_list = [ys[y:y + 2] for y in range(0, len(ys), 2)]
    zs_list = [zs[z:z + 2] for z in range(0, len(zs), 2)]

    n = 0
    while n <= len(xs_list) - 1:
        line = art3d.Line3D(xs_list[n], ys_list[n], zs_list[n], color='k')
        ax.add_line(line)
        n += 1
    ax.scatter(molecular_coordinates_df["x"], molecular_coordinates_df["y"], molecular_coordinates_df["z"],
               color=molecular_coordinates_df["color"], s=40)
    ax.view_init(30,i*5)


anim = animation.FuncAnimation(fig, updatefig, frames=2000, interval=2000)

# Add bounding box to make axes equal
max_range = np.array([molecular_coordinates_df["x"].max() - molecular_coordinates_df["x"].min(),
                      molecular_coordinates_df["y"].max() - molecular_coordinates_df["y"].min(),
                      molecular_coordinates_df["z"].max() - molecular_coordinates_df["z"].min()]).max()
Xb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][0].flatten() + 0.5 * (
            molecular_coordinates_df["x"].max() - molecular_coordinates_df["x"].min())
Yb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][1].flatten() + 0.5 * (
            molecular_coordinates_df["y"].max() - molecular_coordinates_df["y"].min())
Zb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][2].flatten() + 0.5 * (
            molecular_coordinates_df["z"].max() - molecular_coordinates_df["z"].min())

for xb, yb, zb in zip(Xb, Yb, Zb):
    ax.plot([xb], [yb], [zb], 'w')

# Plot the figure
f = r"C:\Users\abate\Desktop\Python Programming\molecule model\rotating_adenosine.gif"
writergif = animation.PillowWriter(fps=35)
anim.save(f, writer=writergif)
# plt.show()








