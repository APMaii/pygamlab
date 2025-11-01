# gamvis.py
import numpy as np
from typing import List, Optional
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go
# ASE imports
try:
    import ase
    from ase import Atoms
    from ase.visualize import view
    HAS_ASE = True
except ImportError:
    HAS_ASE = False


# gamvis/visualizer.py
import numpy as np
from typing import List, Tuple
from ..Primatom import GAM_Atom , GAM_Bond , GAM_Molecule
#from gam_atom import GAM_Atom , GAM_Bond
#from .gamvis import MolecularVisualizer
from. gamvis_engine import GAMVisualizer
from PyQt5.QtWidgets import QApplication, QMessageBox
import sys
'''
try:
    import vtk
    HAS_vtk=True
except ImportError:
    HAS_vtk=False
'''
try:
    import pyvista as pv
    HAS_pv=True

except ImportError:
    HAS_pv=False




        
# Default colors for some common elements
ELEMENT_DEF_COLORS = {
    "H": "white",
    "C": "black",
    "O": "red",
    "N": "blue",
    "S": "yellow",
    "P": "orange",
    # Add more as needed
}



ELEMENT_SIZES = {
    "H": 5,
    "C": 15,
    "O": 12,
    "N": 12,
    "S": 18,
    "P": 18,
}


#---------
# Example color dictionary for common elements (extend as needed)
ATOM_COLORS = {
    "H": (1.0, 1.0, 1.0),   # white
    "C": (0.2, 0.2, 0.2),   # dark gray
    "N": (0.0, 0.0, 1.0),   # blue
    "O": (1.0, 0.0, 0.0),   # red
    "S": (1.0, 1.0, 0.0),   # yellow
    "P": (1.0, 0.6, 0.0),   # orange
}


ELEMENT_RADII = {
    "H": 1.2, "C": 1.7, "N": 1.55, "O": 1.52,
    "S": 1.8, "P": 1.8, "F": 1.47, "Cl": 1.75,
    "Br": 1.85, "I": 1.98, "Fe": 2.0, "Ca": 2.31,
    "Mg": 1.73
}


#=============================================
"""             2D    PLOT               """
#=============================================


def plot_2d(atoms: List["GAM_Atom"], 
            title: str = "",
            show_axis: bool = True,
            show_legend: bool = True,
            figsize: tuple = (6,6),
            xlim: Optional[tuple] = None,
            ylim: Optional[tuple] = None):
    """
    Plot a 2D projection of GAM_Atom objects using matplotlib.

    Args:
        atoms (List[GAM_Atom]): List of atom objects.
        title (str, optional): Plot title.
        show_axis (bool, optional): Whether to display axis.
        show_legend (bool, optional): Whether to show legend.
        figsize (tuple, optional): Size of the figure.
        xlim (tuple, optional): Limits for x-axis (min, max).
        ylim (tuple, optional): Limits for y-axis (min, max).
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Group atoms by element for coloring and legend
    elements = {}
    for atom in atoms:
        if atom.element not in elements:
            elements[atom.element] = {'x': [], 'y': []}
        elements[atom.element]['x'].append(atom.x)
        elements[atom.element]['y'].append(atom.y)
    
    # Plot atoms
    for element, coords in elements.items():
        color = ELEMENT_DEF_COLORS.get(element, "gray")  # default gray if unknown
        ax.scatter(coords['x'], coords['y'], color=color, label=element, s=100, edgecolors='k')
    
    ax.set_title(title)
    ax.axis('on' if show_axis else 'off')
    if show_legend:
        ax.legend()
    
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)
    
    plt.show()










#=============================================
"""             3D    PLOT               """
#=============================================


def plot_3d(atoms: List["GAM_Atom"],
            title: str = "",
            show_axis: bool = True,
            show_legend: bool = True,
            figsize: tuple = (8,8),
            elev: float = 30,
            azim: float = 45,
            xlim: Optional[tuple] = None,
            ylim: Optional[tuple] = None,
            zlim: Optional[tuple] = None):
    """
    Plot a professional 3D view of GAM_Atom objects using matplotlib.

    Args:
        atoms (List[GAM_Atom]): List of atom objects.
        title (str, optional): Plot title.
        show_axis (bool, optional): Show axis.
        show_legend (bool, optional): Show legend.
        figsize (tuple, optional): Figure size.
        elev (float, optional): Camera elevation angle.
        azim (float, optional): Camera azimuth angle.
        xlim, ylim, zlim (tuple, optional): Axis limits for zooming.
    """
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')

    # Group atoms by element
    elements = {}
    for atom in atoms:
        if atom.element not in elements:
            elements[atom.element] = {'x': [], 'y': [], 'z': []}
        elements[atom.element]['x'].append(atom.x)
        elements[atom.element]['y'].append(atom.y)
        elements[atom.element]['z'].append(atom.z)

    # Plot atoms
    for element, coords in elements.items():
        color = ELEMENT_DEF_COLORS.get(element, "gray")
        size = ELEMENT_SIZES.get(element, 100)
        ax.scatter(coords['x'], coords['y'], coords['z'], 
                   color=color, s=size, label=element, edgecolors='k', depthshade=True)

    ax.set_title(title, fontsize=16, weight='bold')
    ax.axis('on' if show_axis else 'off')

    if show_legend:
        ax.legend(title="Elements", fontsize=12)

    # Camera view
    ax.view_init(elev=elev, azim=azim)

    # Axis limits for zoom
    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)
    if zlim: ax.set_zlim(zlim)

    plt.show()
    
    
    






#=============================================
"""          2D    PLOT  PLOTLY            """
#=============================================

def plot_3d_plotly(atoms: List["GAM_Atom"], 
                   title: str = "3D Molecule"):
    """
    Professional interactive 3D molecular plot using Plotly.

    Args:
        atoms (List[GAM_Atom]): List of atom objects.
        title (str): Plot title.
    """
    fig = go.Figure()

    # Group atoms by element
    elements = {}
    for atom in atoms:
        if atom.element not in elements:
            elements[atom.element] = {'x': [], 'y': [], 'z': []}
        elements[atom.element]['x'].append(atom.x)
        elements[atom.element]['y'].append(atom.y)
        elements[atom.element]['z'].append(atom.z)

    # Add scatter for each element
    for element, coords in elements.items():
        fig.add_trace(go.Scatter3d(
            x=coords['x'],
            y=coords['y'],
            z=coords['z'],
            mode='markers',
            marker=dict(
                size=ELEMENT_SIZES.get(element, 10),
                color=ELEMENT_DEF_COLORS.get(element, "gray"),
                line=dict(width=1, color='black')
            ),
            name=element
        ))

    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z',
            aspectmode='data'  # ensures axes scaling matches coordinates
        ),
        margin=dict(l=0, r=0, b=0, t=30)
    )

    fig.show()
    




#=============================================
"""             3D ANIMATION               """
#=============================================


def animate_3d_vibration(atoms, 
                         title="Dynamic 3D Molecule", 
                         frames=200, 
                         vibration_strength=0.05):
    """
    Animate atoms with tiny vibrations in 3D using Matplotlib.
    
    Args:
        atoms: List of GAM_Atom objects
        title: Plot title
        frames: Number of animation frames
        vibration_strength: Maximum random displacement per axis
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(title)

    # Extract coordinates and elements
    coords = np.array([[atom.x, atom.y, atom.z] for atom in atoms])
    elements = [atom.element for atom in atoms]

    # Initialize scatter plot
    scatter = ax.scatter(coords[:,0], coords[:,1], coords[:,2],
                         s=[ELEMENT_SIZES.get(el, 100) for el in elements],
                         c=[ELEMENT_DEF_COLORS.get(el, "gray") for el in elements],
                         edgecolors='k')

    # Axis labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Equal aspect ratio
    max_range = np.array([coords[:,0].max()-coords[:,0].min(), 
                          coords[:,1].max()-coords[:,1].min(), 
                          coords[:,2].max()-coords[:,2].min()]).max() / 2.0

    mid_x = (coords[:,0].max()+coords[:,0].min()) * 0.5
    mid_y = (coords[:,1].max()+coords[:,1].min()) * 0.5
    mid_z = (coords[:,2].max()+coords[:,2].min()) * 0.5

    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    # Update function for animation
    def update(frame):
        vibration = (np.random.rand(*coords.shape) - 0.5) * 2 * vibration_strength
        new_coords = coords + vibration
        scatter._offsets3d = (new_coords[:,0], new_coords[:,1], new_coords[:,2])
        return scatter,

    # Keep a reference to the animation to prevent garbage collection
    ani = FuncAnimation(fig, update, frames=frames, interval=50, blit=False)
    plt.show()
    return ani  # <-- Important to keep the animation alive

    



#=============================================
"""             3D Pyvista               """
#=============================================

def pyvista_visualizer(atoms: List[GAM_Atom],
                     show_bonds: bool = True,
                     background_color: str = 'white',
                     window_size: Tuple[int, int] = (800, 600),
                     bond_threshold: float = 1.6,
                     atom_radius: float = 0.4):
    """
    Quick 3D visualization of atoms using PyVista.

    Parameters
    ----------
    atoms : list[GAM_Atom]
        List of atoms to visualize.
    show_bonds : bool
        Whether to display bonds between atoms (distance-based).
    background_color : str
        Background color of the visualization.
    window_size : tuple
        Size of the visualization window.
    bond_threshold : float
        Maximum distance to consider two atoms bonded.
    atom_radius : float
        Default radius for atom spheres.
    """
    if HAS_pv:

        plotter = pv.Plotter(window_size=window_size)
        plotter.set_background(background_color)

        positions = []
        for atom in atoms:
            pos = np.array([atom.x, atom.y, atom.z])
            positions.append(pos)

            # Get atom color
            color = ATOM_COLORS.get(atom.element, (0.5, 0.5, 0.5))  # default gray
            sphere = pv.Sphere(radius=atom_radius, center=pos)
            plotter.add_mesh(sphere, color=color, smooth_shading=True)

        positions = np.array(positions)

        # Add bonds (distance-based)
        if show_bonds:
            for i in range(len(atoms)):
                for j in range(i + 1, len(atoms)):
                    p1, p2 = positions[i], positions[j]
                    dist = np.linalg.norm(p2 - p1)
                    if dist <= bond_threshold:  # simple cutoff
                        direction = p2 - p1
                        direction /= np.linalg.norm(direction)

                        cylinder = pv.Cylinder(
                            center=(p1 + p2) / 2,
                            direction=direction,
                            radius=0.1,
                            height=dist
                        )
                        plotter.add_mesh(cylinder, color="gray", smooth_shading=True)

        plotter.camera_position = 'iso'
        plotter.show()
    else:
        print('You must install pyvista with "pip install pyvista"')
        return None


#=============================================
"""       EFFICIENT HIGH-QUALITY VIEWERS    """
#=============================================

def plot_3d_efficient_plotly(atoms: List["GAM_Atom"], 
                             title: str = "Efficient 3D Molecule",
                             camera_distance: float = 1.5,
                             show_bonds: bool = False,
                             bond_threshold: float = 1.6):
    """
    Efficient high-quality interactive 3D molecular viewer using Plotly.
    Optimized for performance with large molecules while maintaining visual quality.

    Args:
        atoms (List[GAM_Atom]): List of atom objects.
        title (str): Plot title.
        camera_distance (float): Camera distance multiplier for auto-zoom.
        show_bonds (bool): Whether to show bonds between atoms.
        bond_threshold (float): Maximum distance to consider atoms bonded.
    """
    fig = go.Figure()

    # Extract coordinates efficiently
    coords = np.array([[atom.x, atom.y, atom.z] for atom in atoms])
    elements = [atom.element for atom in atoms]
    
    # Group atoms by element for efficient rendering
    element_groups = {}
    for i, atom in enumerate(atoms):
        if atom.element not in element_groups:
            element_groups[atom.element] = {'indices': [], 'coords': []}
        element_groups[atom.element]['indices'].append(i)
        element_groups[atom.element]['coords'].append([atom.x, atom.y, atom.z])

    # Add scatter traces by element (more efficient than individual atoms)
    for element, data in element_groups.items():
        coords_array = np.array(data['coords'])
        
        # Create hover text with atom information
        hover_text = [f"Atom {data['indices'][i]+1}<br>Element: {element}<br>Position: ({coords_array[i,0]:.2f}, {coords_array[i,1]:.2f}, {coords_array[i,2]:.2f})" 
                     for i in range(len(coords_array))]
        
        fig.add_trace(go.Scatter3d(
            x=coords_array[:, 0],
            y=coords_array[:, 1],
            z=coords_array[:, 2],
            mode='markers',
            marker=dict(
                size=ELEMENT_SIZES.get(element, 10),
                color=ELEMENT_DEF_COLORS.get(element, "gray"),
                line=dict(width=1, color='black'),
                opacity=0.9
            ),
            name=element,
            hovertext=hover_text,
            hoverinfo='text'
        ))

    # Add bonds if requested (simplified for performance)
    if show_bonds and len(atoms) < 100:  # Only for smaller molecules
        bond_lines = []
        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                dist = np.linalg.norm(coords[j] - coords[i])
                if dist <= bond_threshold:
                    bond_lines.extend([coords[i], coords[j], [None, None, None]])
        
        if bond_lines:
            bond_array = np.array(bond_lines)
            fig.add_trace(go.Scatter3d(
                x=bond_array[:, 0],
                y=bond_array[:, 1],
                z=bond_array[:, 2],
                mode='lines',
                line=dict(color='gray', width=3),
                showlegend=False,
                hoverinfo='skip'
            ))

    # Calculate auto-zoom based on molecular extent
    if len(coords) > 0:
        center = coords.mean(axis=0)
        extent = np.max(coords, axis=0) - np.min(coords, axis=0)
        max_extent = np.max(extent)
        camera_eye = dict(x=camera_distance, y=camera_distance, z=camera_distance)
    else:
        camera_eye = dict(x=1.5, y=1.5, z=1.5)

    fig.update_layout(
        title=dict(text=title, x=0.5, font=dict(size=16)),
        scene=dict(
            xaxis_title='X (Å)',
            yaxis_title='Y (Å)',
            zaxis_title='Z (Å)',
            aspectmode='cube',  # Equal aspect ratio
            camera=dict(eye=camera_eye),
            bgcolor='rgba(240,240,240,0.1)'
        ),
        margin=dict(l=0, r=0, b=0, t=40),
        font=dict(size=12),
        showlegend=True,
        legend=dict(x=0.02, y=0.98),
        width=800,
        height=600
    )

    fig.show()


def plot_3d_lightweight_matplotlib(atoms: List["GAM_Atom"],
                                  title: str = "Lightweight 3D Molecule",
                                  style: str = "professional",
                                  figsize: tuple = (10, 8),
                                  dpi: int = 100,
                                  save_path: Optional[str] = None):
    """
    Lightweight but high-quality 3D molecular viewer using matplotlib.
    Optimized for fast rendering with professional appearance.

    Args:
        atoms (List[GAM_Atom]): List of atom objects.
        title (str): Plot title.
        style (str): Visual style ('professional', 'minimal', 'dark').
        figsize (tuple): Figure size in inches.
        dpi (int): Resolution for rendering.
        save_path (str, optional): Path to save the figure.
    """
    
    # Set style
    if style == "professional":
        plt.style.use('seaborn-v0_8-whitegrid')
        bg_color = 'white'
        grid_alpha = 0.3
    elif style == "minimal":
        plt.style.use('seaborn-v0_8-white')
        bg_color = 'white'
        grid_alpha = 0.1
    elif style == "dark":
        plt.style.use('dark_background')
        bg_color = 'black'
        grid_alpha = 0.2
    else:
        bg_color = 'white'
        grid_alpha = 0.3

    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor(bg_color)

    # Extract data efficiently
    coords = np.array([[atom.x, atom.y, atom.z] for atom in atoms])
    
    # Group by element for efficient plotting
    element_data = {}
    for atom in atoms:
        if atom.element not in element_data:
            element_data[atom.element] = {'coords': [], 'count': 0}
        element_data[atom.element]['coords'].append([atom.x, atom.y, atom.z])
        element_data[atom.element]['count'] += 1

    # Plot atoms by element with optimized settings
    for element, data in element_data.items():
        coords_array = np.array(data['coords'])
        color = ELEMENT_DEF_COLORS.get(element, "gray")
        size = ELEMENT_SIZES.get(element, 100)
        
        # Use efficient scatter plotting
        scatter = ax.scatter(coords_array[:, 0], coords_array[:, 1], coords_array[:, 2],
                           c=color, s=size, alpha=0.8, edgecolors='black', 
                           linewidths=0.5, label=f"{element} ({data['count']})")

    # Professional styling
    ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
    ax.set_xlabel('X (Å)', fontsize=12, labelpad=10)
    ax.set_ylabel('Y (Å)', fontsize=12, labelpad=10)
    ax.set_zlabel('Z (Å)', fontsize=12, labelpad=10)

    # Grid and appearance
    ax.grid(True, alpha=grid_alpha)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    
    # Make pane edges more subtle
    ax.xaxis.pane.set_edgecolor('gray')
    ax.yaxis.pane.set_edgecolor('gray')
    ax.zaxis.pane.set_edgecolor('gray')
    ax.xaxis.pane.set_alpha(0.1)
    ax.yaxis.pane.set_alpha(0.1)
    ax.zaxis.pane.set_alpha(0.1)

    # Equal aspect ratio and optimal viewing angle
    if len(coords) > 0:
        # Calculate center and extent
        center = coords.mean(axis=0)
        extent = coords.max(axis=0) - coords.min(axis=0)
        max_extent = extent.max()
        
        # Set equal limits around center
        margin = max_extent * 0.1  # 10% margin
        ax.set_xlim(center[0] - max_extent/2 - margin, center[0] + max_extent/2 + margin)
        ax.set_ylim(center[1] - max_extent/2 - margin, center[1] + max_extent/2 + margin)
        ax.set_zlim(center[2] - max_extent/2 - margin, center[2] + max_extent/2 + margin)

    # Optimal viewing angle
    ax.view_init(elev=20, azim=45)

    # Legend with optimized positioning
    legend = ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', 
                      frameon=True, fancybox=True, shadow=True)
    legend.get_frame().set_alpha(0.9)

    # Tight layout for better appearance
    plt.tight_layout()

    # Save if requested
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight', 
                   facecolor=bg_color, edgecolor='none')
        print(f"Figure saved to: {save_path}")

    plt.show()



#=============================================
"""             3D GAMVIS               """
#=============================================

#from gamlab.gamvisualizer import MolecularVisualizer


def gamvisualizer(atoms):
    """
    Main function to visualize molecules in an external window.
    
    Args:
        atoms: List of GAM_Atom objects or compatible atom objects with:
               - element: string (e.g., 'C', 'H', 'O', 'N')
               - x, y, z: float coordinates
    
    Usage:
        from professional_gamvis2 import gamvisualizer
        atoms = [GAM_Atom('C', 0, 0, 0), GAM_Atom('H', 1, 0, 0)]
        gamvisualizer(atoms)
    """
    print("=== PyGAMLab Visualizator Debug Info ===")
    print(f"Python version: {sys.version}")
    
    try:
        from PyQt5.QtCore import QT_VERSION_STR
        print(f"PyQt5 version: {QT_VERSION_STR}")
    except ImportError:
        print("PyQt5 version: Unknown")
    
    # Check if we're on a system that might need specific VTK setup
    import platform
    print(f"Platform: {platform.system()} {platform.release()}")
    
    # Create QApplication if it doesn't exist
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    
    # Set application style
    app.setStyle('Fusion')
    
    try:
        # Create and show the main window
        print("Creating main window...")
        # Don't load sample molecule if atoms are provided
        window = GAMVisualizer()
        
        # Load the provided atoms
        if atoms:
            print(f"Loading {len(atoms)} atoms...")
            window.load_atoms(atoms)
        else:
            print("No atoms provided, using sample molecule...")
        
        print("Showing window...")
        window.show()
        
        # Start the VTK interaction
        print("Initializing VTK interactor...")
        window.interactor.Initialize()
        print("Application ready!")
        
        # Start the event loop
        app.exec_()
        
    except Exception as e:
        print(f"ERROR: {str(e)}")
        import traceback
        traceback.print_exc()
        
        # Show error message
        if QApplication.instance():
            QMessageBox.critical(None, "Error", f"Failed to start application:\n{str(e)}")


#=============================================
"""              ASE INTEGRATION            """
#=============================================

def to_ase(atoms: List["GAM_Atom"]) -> Optional['ase.Atoms']:
    """
    Convert GAM_Atom objects to ASE Atoms object (if ASE is available).
    
    Args:
        atoms: List of GAM_Atom objects to convert
        
    Returns:
        ASE Atoms object or None if ASE not available
    """
    if not HAS_ASE:
        print("ASE not available. Install with: pip install ase")
        return None
    
    if not atoms:
        return Atoms()
    
    symbols = [atom.element for atom in atoms]
    positions = [[atom.x, atom.y, atom.z] for atom in atoms]
    
    atoms_obj = Atoms(symbols=symbols, positions=positions)
    
    return atoms_obj


def ase_visualizer(atoms: List["GAM_Atom"], 
                   viewer: str = "ase-gui",
                   **kwargs):
    """
    Visualize atoms using ASE's built-in visualizer.
    
    Args:
        atoms: List of GAM_Atom objects
        viewer: ASE viewer to use ('ase-gui', 'vmd', 'rasmol', etc.)
        **kwargs: Additional arguments passed to ase.visualize.view()
    """
    if not HAS_ASE:
        print("ASE not available. Install with: pip install ase")
        return None
    
    ase_atoms = to_ase(atoms)
    if ase_atoms is None:
        return None
    
    try:
        from ase.visualize import view
        view(ase_atoms, viewer=viewer, **kwargs)
        print(f"ASE visualization opened with {viewer} viewer")
    except Exception as e:
        print(f"Error opening ASE visualizer: {e}")
        print("Make sure the specified viewer is installed and available")


#=============================================
"""          MAIN VISUALIZATION FUNCTION    """
#=============================================

def Molecular_Visualizer(atoms: List["GAM_Atom"], 
                        format: str = "efficient_plotly",
                        **kwargs):
    """
    High-level molecular visualization interface supporting multiple backends and formats.

    This function provides a unified entry point to visualize molecular structures
    represented as a list of `GAM_Atom` objects. Users can choose from various
    rendering backends, ranging from lightweight static plots to high-end interactive
    3D visualizations with vibration animations.

    Parameters
    ----------
    atoms : List[GAM_Atom]
        A list of `GAM_Atom` instances representing the atoms of the molecule to visualize.
    format : str, optional
        Visualization backend/format to use. Default is `'efficient_plotly'`.
        Supported options include:
            - `'efficient_plotly'` : Interactive 3D visualization with Plotly (optimized for performance, recommended)
            - `'lightweight_matplotlib'` : Fast static 3D plot using Matplotlib
            - `'plotly'` : Standard interactive 3D visualization with Plotly
            - `'matplotlib'` : Standard 3D visualization using Matplotlib
            - `'2d'` : 2D projection of molecular structure
            - `'pyvista'` : Advanced 3D rendering using PyVista
            - `'ase'` : Visualize using ASE's built-in visualizer
            - `'gamvis'` : Professional VTK-based external 3D window
            - `'animation'` : Animated 3D visualization with vibration motion
    **kwargs : dict
        Additional keyword arguments specific to the chosen visualization backend.
        For example, camera angles, color maps, sphere sizes, or animation settings.

    Returns
    -------
    Optional[Any]
        Some visualization backends (e.g., `'animation'`) may return a handle
        or object representing the rendered scene; most backends render directly
        and return `None`.

    Notes
    -----
    - If `atoms` is empty or `None`, the function prints a message and exits.
    - Use `'efficient_plotly'` for interactive exploration in notebooks.
    - The `'gamvis'` format opens a dedicated VTK window with professional
      controls for rotation, lighting, and measurement.
    - All backend-specific options can be passed via `**kwargs` to customize the visualization.

    Examples
    --------
    >>> molecular_visualizer(atoms_list)
    >>> molecular_visualizer(atoms_list, format='pyvista', background='white')
    >>> molecular_visualizer(atoms_list, format='animation', vibration_amplitude=0.05)
    """
    if not atoms:
        print("No atoms provided for visualization")
        return
    
    print(f"Visualizing {len(atoms)} atoms using {format} format...")
    
    # Route to appropriate visualizer based on format
    if format == "efficient_plotly":
        plot_3d_efficient_plotly(atoms, **kwargs)
    elif format == "lightweight_matplotlib":
        plot_3d_lightweight_matplotlib(atoms, **kwargs)
    elif format == "plotly":
        plot_3d_plotly(atoms, **kwargs)
    elif format == "matplotlib":
        plot_3d(atoms, **kwargs)
    elif format == "2d":
        plot_2d(atoms, **kwargs)
    elif format == "pyvista":
        pyvista_visualizer(atoms, **kwargs)
    elif format == "ase":
        ase_visualizer(atoms, **kwargs)
    elif format == "gamvis":
        gamvisualizer(atoms)
    elif format == "animation":
        return animate_3d_vibration(atoms, **kwargs)
    else:
        print(f"Unknown format '{format}'. Available formats:")
        print("  - efficient_plotly (recommended)")
        print("  - lightweight_matplotlib")
        print("  - plotly")
        print("  - matplotlib")
        print("  - 2d")
        print("  - pyvista")
        print("  - ase")
        print("  - gamvis")
        print("  - animation")




