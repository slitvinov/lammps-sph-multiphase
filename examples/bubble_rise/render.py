import sys
import numpy as np
from ovito.io import import_file
from ovito.vis import Viewport, OpenGLRenderer

src, out, step, nx = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
pipeline = import_file(src + "/dump*.dat", columns=["Particle Identifier", "Particle Type", "Position.X", "Position.Y", "Position.Z", "Velocity.X", "Velocity.Y", "Velocity.Z", "rho", "e", "T"])

def paint(frame, data):
    liquid = np.asarray(data.particles["Particle Type"]) == 1
    c = data.particles_.create_property("Color")
    c[...] = np.where(liquid[:, None], (0.55, 0.55, 0.55), (0.0, 0.0, 0.5))
    t = data.particles_.create_property("Transparency")
    t[...] = np.where(liquid, 0.97, 0.0)

pipeline.modifiers.append(paint)
pipeline.add_to_scene()
pipeline.compute(0).particles.vis.radius = 0.5 / nx
vp = Viewport(type=Viewport.Type.Perspective, camera_pos=(1.6, -1.4, 0.9), camera_dir=(-1.1, 1.9, -0.4), fov=0.55)
vp.zoom_all()
r = OpenGLRenderer(antialiasing_level=3, order_independent_transparency=True)
n = pipeline.source.num_frames
fr = [int(a) for a in sys.argv[5].split(',')] if len(sys.argv) > 5 else range(0, n, step)
for k, f in enumerate(fr):
    vp.render_image(filename=f"{out}/f{k:04d}.png", size=(600, 600), frame=f, renderer=r, background=(1, 1, 1))
print("frames", n)
