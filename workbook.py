# -*- coding: utf-8 -*-
"""
Created on Thu May 12 09:28:16 2022

@author: ashley
"""

from diffpy.structure import Atom, Lattice, Structure
import numpy as np
from orix.crystal_map import CrystalMap
from orix.quaternion.rotation import Rotation

import tempfile

from diffpy.structure import Atom, Lattice, Structure
import matplotlib.pyplot as plt
import numpy as np

import maxflow
from orix import data, io, plot
from orix.crystal_map import CrystalMap, Phase, PhaseList
from orix.quaternion import Orientation, Rotation, symmetry
from orix.vector import Vector3d
plt.close('all')
path = r'C:\Users\ashle\Documents\GitHub\maxflow_for_matsci\Data\steel_ebsd.ang'

# Directly access *private* cache data path from module
# _target = data._fetcher.path / "sdss/sdss_ferrite_austenite.ang"

# Read each column from the file
euler1, euler2, euler3, x, y, iq, dp, phase_id, sem, fit  = np.loadtxt(path, unpack=True)

# Create a Rotation object from Euler angles
euler_angles = np.column_stack((euler1, euler2, euler3))
rotations = Rotation.from_euler(euler_angles)

# Create a property dictionary
properties = dict(iq=iq, dp=dp)

# Create unit cells of the phases
structures = [
    Structure(
        title="ferrite",
        atoms=[Atom("fe", [0] * 3)],
        lattice=Lattice(0.287, 0.287, 0.287, 90, 90, 90)
    ),
]
phase_list = PhaseList(
    names=["ferrite"],
    point_groups=["432"],
    structures=structures,
)

# Create a CrystalMap instance
xmap2 = CrystalMap(
    rotations=rotations,
    phase_id=phase_id,
    x=x,
    y=y,
    phase_list=phase_list,
    prop=properties,
)
xmap2.scan_unit = "um"

# print(xmap2)
# xmap2.plot(overlay='dp')
ipw = 1  #inplane weight--should be bet
ckey_m3m = plot.IPFColorKeyTSL(xmap2.phases["ferrite"].point_group, direction=Vector3d.zvector())
rgb_fe = ckey_m3m.orientation2color(xmap2["ferrite"].orientations)

fer_x = np.round(2*(xmap2['ferrite'].x))
fer_y = np.round(2*(xmap2['ferrite'].y))

g = maxflow.GraphFloat()
nodeids = g.add_grid_nodes((305, 305)) #int(np.sqrt(len(dp)))

for_network = fer_x, fer_y, rgb_fe[:,0], rgb_fe[:,1], rgb_fe[:,2]
for_network = np.asarray(for_network).T
# img = np.round(255*rgb_fe[:,0])
img = np.zeros((305,305))
for xx in range(len(for_network)):
    coordx, coordy = int(for_network[xx,0]), int(for_network[xx,1])
    img[coordx, coordy] = for_network[xx,2]
img = img.T
IDD = xmap2['ferrite'].id
g.add_grid_edges(nodeids, ipw)
g.add_grid_tedges(nodeids, img, 1-img)

bob = g.get_nx_graph()
connectivity = []
for (u,v,wt) in bob.edges.data('weight'):
    connectivity.append((u,v,wt))
connectivity = np.asarray(connectivity)
# what is difference between connectivity and adjacency matrix

g.maxflow()
sgm = g.get_grid_segments(nodeids)
img2 = np.int_(np.logical_not(sgm))
from matplotlib import pyplot as ppl
ppl.imshow(img2)
ppl.show()


