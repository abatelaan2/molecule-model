import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d
plt.rcParams["figure.figsize"] = 6.4, 4.8

# Define figure and axis
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Import coordinate and connectivity data
molecular_coordinates_df = pd.read_excel(r"C:\Users\abate\Desktop\Python Programming\molecule model\Adenosine_Coordinates.xlsx", sheet_name=0)
molecular_connectivity_df = pd.read_excel(r"C:\Users\abate\Desktop\Python Programming\molecule model\Adenosine_Coordinates.xlsx", sheet_name=1)

# Assign colors
n=0
list=[]
while n<= len(molecular_coordinates_df["Atom"])-1:
    if molecular_coordinates_df["Atom"][n][0]=="C":
        list.append('0.3')
    elif molecular_coordinates_df["Atom"][n][0]=="O":
        list.append('r')
    elif molecular_coordinates_df["Atom"][n][0] == "N":
        list.append('b')
    elif molecular_coordinates_df["Atom"][n][0] == "H":
        list.append('0.7')
    elif molecular_coordinates_df["Atom"][n][0] == "O":
        list.append('r')
    elif molecular_coordinates_df["Atom"][n][0] == "P":
        list.append('m')
    n+=1
molecular_coordinates_df["color"]=list

# Draw bonds
xs=[]
ys=[]
zs=[]
n=0
m=0
while n<=len(molecular_connectivity_df)-1:
    while m<= len(molecular_coordinates_df)-1:
        if molecular_connectivity_df["Atom Name"][n]==molecular_coordinates_df["Atom"][m]:
            xs.append(molecular_coordinates_df["x"][m])
            ys.append(molecular_coordinates_df["y"][m])
            zs.append(molecular_coordinates_df["z"][m])
        elif molecular_connectivity_df["Connectivity"][n]==molecular_coordinates_df["Atom"][m]:
            xs.append(molecular_coordinates_df["x"][m])
            ys.append(molecular_coordinates_df["y"][m])
            zs.append(molecular_coordinates_df["z"][m])
        m+=1
    m=0
    n+=1

xs_list = [xs[x:x+2] for x in range(0, len(xs),2)]
ys_list = [ys[y:y+2] for y in range(0, len(ys),2)]
zs_list = [zs[z:z+2] for z in range(0, len(zs),2)]

n=0
while n<=len(xs_list)-1:
    line = art3d.Line3D(xs_list[n], ys_list[n], zs_list[n], color='k')
    ax.add_line(line)
    n+=1

# Add bounding box to make axes equal
max_range = np.array([molecular_coordinates_df["x"].max()-molecular_coordinates_df["x"].min(), molecular_coordinates_df["y"].max()-molecular_coordinates_df["y"].min(), molecular_coordinates_df["z"].max()-molecular_coordinates_df["z"].min()]).max()
Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(molecular_coordinates_df["x"].max()-molecular_coordinates_df["x"].min())
Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(molecular_coordinates_df["y"].max()-molecular_coordinates_df["y"].min())
Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(molecular_coordinates_df["z"].max()-molecular_coordinates_df["z"].min())

for xb, yb, zb in zip(Xb, Yb, Zb):
   ax.plot([xb], [yb], [zb], 'w')

# Plot the figure
ax.scatter(molecular_coordinates_df["x"], molecular_coordinates_df["y"], molecular_coordinates_df["z"], color=molecular_coordinates_df["color"], s=40)
plt.grid()
plt.show()
plt.savefig("molecule.png")



