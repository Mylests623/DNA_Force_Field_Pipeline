from pymol import cmd
import numpy as np

def add_force_color_bar(name="force_scale", min_force=0, max_force=100, n_points=11, pos=[-10,0,0], radius=0.2):
    """
    Adds a linear dummy object with n_points colored along blue-white-red gradient.
    - name: PyMOL object name
    - min_force, max_force: numeric range corresponding to your heatmap
    - n_points: number of points in the bar
    - pos: [x,y,z] start position
    - radius: sphere radius for visibility
    """
    x0, y0, z0 = pos
    forces = np.linspace(min_force, max_force, n_points)

    for i, f in enumerate(forces):
        x = x0 + i * 0.5
        y, z = y0, z0
        cmd.pseudoatom(name, pos=[x, y, z], b=f, label=f"{f:.0f}")

    cmd.show("spheres", name)
    cmd.set("sphere_scale", radius, name)
    cmd.spectrum("b", "blue_white_red", name, minimum=min_force, maximum=max_force)
    cmd.rebuild()

cmd.extend("add_force_color_bar", add_force_color_bar)
print("Command 'add_force_color_bar' loaded. Example: add_force_color_bar(min_force=0,max_force=100)")
