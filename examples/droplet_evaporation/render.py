import sys
import numpy as np
from ovito.io import import_file
from ovito.vis import Viewport, TachyonRenderer

src, out, step = sys.argv[1], sys.argv[2], int(sys.argv[3])
pipeline = import_file(src + "/dump*.dat", columns=["Particle Identifier", "Particle Type", "Position.X", "Position.Y", "Position.Z", "Velocity.X", "Velocity.Y", "Velocity.Z", "rho", "e", "T"])

def paint(frame, data):
    t = np.asarray(data.particles["Particle Type"])
    c = data.particles_.create_property("Color")
    c[...] = np.where((t == 1)[:, None], (0.0, 0.0, 0.5), (0.78, 0.78, 0.78))
    data.cell_.vis.enabled = False

pipeline.modifiers.append(paint)
pipeline.add_to_scene()
pipeline.compute(0).particles.vis.radius = 0.5 / 40
vp = Viewport(type=Viewport.Type.Top)
vp.camera_pos = (0.5, 0.5, 1.0)
vp.fov = 0.5
r = TachyonRenderer(ambient_occlusion=False, shadows=False, direct_light_intensity=0.9)
n = pipeline.source.num_frames
for k, f in enumerate(range(0, n, step)):
    vp.render_image(filename=f"{out}/f{k:04d}.png", size=(600, 600), frame=f, renderer=r, background=(1, 1, 1))
print("frames", n, "rendered", len(range(0, n, step)))
