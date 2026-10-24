import sys
import numpy as np
from matplotlib import colormaps
from ovito.io import import_file
from ovito.vis import Viewport, OpenGLRenderer

src, out, nx = sys.argv[1], sys.argv[2], int(sys.argv[3])
tlo, thi = float(sys.argv[4]), float(sys.argv[5])
pipeline = import_file(src + "/dump*.dat", columns=["Particle Identifier", "Particle Type", "Position.X", "Position.Y", "Position.Z", "Velocity.X", "Velocity.Y", "Velocity.Z", "rho", "e", "Mass", "T"])
cmap = colormaps["inferno"]

def paint(frame, data):
    typ = np.asarray(data.particles["Particle Type"])
    T = np.asarray(data.particles["T"])
    c = cmap(np.clip((T - tlo) / (thi - tlo), 0, 1))[:, :3]
    c[typ == 3] = (0.6, 0.6, 0.6)
    data.particles_.create_property("Color")[...] = c
    tr = np.where(typ == 1, 0.55, 0.0)
    data.particles_.create_property("Transparency")[...] = tr
    data.cell_.vis.enabled = False

pipeline.modifiers.append(paint)
pipeline.add_to_scene()
pipeline.compute(0).particles.vis.radius = 0.5 / nx
vp = Viewport(type=Viewport.Type.Top)
vp.zoom_all()
r = OpenGLRenderer(antialiasing_level=3, order_independent_transparency=True)
n = pipeline.source.num_frames
for f in range(n):
    vp.render_image(filename=f"{out}/f{f:04d}.png", size=(600, 600), frame=f, renderer=r, background=(1, 1, 1))
print("frames", n)
