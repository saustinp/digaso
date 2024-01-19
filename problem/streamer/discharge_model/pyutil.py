import numpy as np
import gmsh
import pyvista as pv

def plot(mesh, sol, pltflag):
    """
    Function for plotting thermal fin mesh and solution

    """

    # Concatenating adjacency matrices for the 5 subdomains into a single adjacency matrix for the whole domain
    theta = mesh['theta']
    coor = mesh['coor']

    # The face vector must be 1D, and each element starts with the number of vertices in the element (3 for triangles)
    faces = np.concatenate((3*np.ones((len(theta), 1)), theta), axis=1)
    faces = np.hstack(faces).astype(int)

    # Generate mesh
    mesh = pv.PolyData(coor, faces)

    # Generate plot object
    plot = pv.Plotter(window_size=[2000,1500])
    plot.set_background('white')
    # plot.add_title('2D Thermal Fin FEM',
    #                 color='black'
    #                 ) 

    # Assign mesh to plot object
    colorbar_args = dict(title='u',
                        vertical=0,
                        height=0.08,
                        width=0.6,
                        position_x=0.2,
                        position_y=0.05,
                        color='black'
                        )
    # colorbar_args = dict(title='Temperature',
    #                     vertical=0,
    #                     height=0.08,
    #                     width=0.6,
    #                     position_x=0.2,
    #                     position_y=0.05,
    #                     color='black'
    #                     )

    # Length of scalars vector determines if values are assigned to points (interpolated on element) or assigned to whole face
    plot.add_mesh(mesh, scalars = sol,
                        show_edges = True,
                        line_width=0.75,
                        label = 'u',
                        cmap='plasma',
                        scalar_bar_args=colorbar_args,
                        )

    if pltflag:
        # Render plot
        plot.show(cpos='xy')
        return plot
    else:
        return plot