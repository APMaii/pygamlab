import sys
import numpy as np

# Export list for module usage

# Check and install dependencies
def check_dependencies():
    missing = []
    try:
        import PyQt5
    except ImportError:
        missing.append("PyQt5")
    
    try:
        import vtk
    except ImportError:
        missing.append("vtk")
        
    try:
        import numpy
    except ImportError:
        missing.append("numpy")
        
    if missing:
        print("Missing dependencies:")
        for dep in missing:
            print(f"  - {dep}")
        print("\nInstall with:")
        print(f"pip install {' '.join(missing)}")
        return False
    return True

if not check_dependencies():
    sys.exit(1)

from PyQt5.QtWidgets import (QApplication, QMainWindow, QVBoxLayout, QHBoxLayout, 
                             QWidget, QPushButton, QSlider, QLabel, QComboBox, 
                             QGroupBox, QFileDialog, QMessageBox, QGridLayout,
                             QCheckBox, QSpinBox)
from PyQt5.QtCore import QTimer, Qt
from PyQt5.QtGui import QFont, QPalette, QColor, QPixmap

import vtk

# Try different VTK Qt integration methods
VTK_QT_WIDGET = None
try:
    # Method 1: Modern VTK
    from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
    VTK_QT_WIDGET = QVTKRenderWindowInteractor
    print("Using vtkmodules.qt.QVTKRenderWindowInteractor")
except ImportError:
    try:
        # Method 2: Older VTK
        from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
        VTK_QT_WIDGET = QVTKRenderWindowInteractor
        print("Using vtk.qt.QVTKRenderWindowInteractor")
    except ImportError:
        try:
            # Method 3: Alternative import
            from vtkmodules.qt5.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
            VTK_QT_WIDGET = QVTKRenderWindowInteractor
            print("Using vtkmodules.qt5.QVTKRenderWindowInteractor")
        except ImportError:
            print("ERROR: Could not import VTK Qt widget!")
            print("Try installing with: pip install vtk[qt]")
            sys.exit(1)

# Mock GAM_Atom class for demonstration
'''
class GAM_Atom:
    def __init__(self, element, x, y, z):
        self.element = element
        self.coords = np.array([x, y, z])
        self.x, self.y, self.z = x, y, z
'''

# Element color scheme (RGB 0-1)
ELEMENT_COLORS = {
    "H": (1.0, 1.0, 1.0),     # White
    "C": (0.3, 0.3, 0.3),     # Dark gray
    "O": (1.0, 0.0, 0.0),     # Red
    "N": (0.0, 0.0, 1.0),     # Blue
    "S": (1.0, 1.0, 0.0),     # Yellow
    "P": (1.0, 0.5, 0.0),     # Orange
    "F": (0.0, 1.0, 0.0),     # Green
    "Cl": (0.0, 0.8, 0.0),    # Dark green
    "Br": (0.6, 0.2, 0.0),    # Brown
    "I": (0.4, 0.0, 0.8),     # Purple
    "Fe": (0.8, 0.5, 0.2),    # Iron brown
    "Ca": (0.8, 0.8, 0.8),    # Light gray
    "Mg": (0.5, 0.8, 0.5),    # Light green
}

# Element radii (van der Waals radii in Angstroms)
ELEMENT_RADII = {
    "H": 1.2, "C": 1.7, "N": 1.55, "O": 1.52,
    "S": 1.8, "P": 1.8, "F": 1.47, "Cl": 1.75,
    "Br": 1.85, "I": 1.98, "Fe": 2.0, "Ca": 2.31,
    "Mg": 1.73
}

class GAMVisualizer(QMainWindow):
    """
    Internal method to render all loaded atoms as 3D spheres in the VTK scene.

    This function handles the creation, coloring, and placement of 
    atom representations (spheres) in the VTK renderer. It also manages 
    the display of atom labels if enabled. 

    Responsibilities:
    - Clear any previously rendered spheres or labels from the scene.
    - Generate a VTK sphere for each atom, scaling the radius based on 
      element-specific covalent/atomic radii and the user-selected size factor.
    - Assign colors to spheres according to element type.
    - Add all sphere actors to the VTK renderer.
    - Render atom labels if the "Show Atom Labels" option is enabled.
    - Adjust the camera to fit the entire molecule in view.

    Notes
    -----
    - This method does not compute or render bonds between atoms; 
      bonding visualization should be handled by a separate function.
    - The function is automatically called when atoms are loaded or 
      when atom sizes, colors, or visibility options change.
    - Intended for internal use within `GAMVisualizer`; external code 
      should call `load_atoms()` to populate the visualizer and trigger rendering.

    Side Effects
    ------------
    - Modifies `self.sphere_actors` and `self.label_actors` lists.
    - Updates the VTK render window to reflect the new visualization.

    Example
    -------
    # Internal usage after loading atoms
    visualizer.load_atoms(atom_list)
    # This automatically calls render_molecules()
    """
    def __init__(self):
        super().__init__()
        self.atoms = []
        self.sphere_actors = []
        self.label_actors = []
        self.original_positions = []
        self.animation_timer = QTimer()
        self.animation_frame = 0
        self.vibration_amplitude = 0.05
        self.main_light = None
        self.axes_actor = None
        self.measurement_mode = False
        self.measurement_points = []
        
        self.init_ui()
        self.init_vtk()
        #if load_sample:
        #    self.setup_sample_molecule()
        
    def init_ui(self):
        self.setWindowTitle("PyGAMLab Visualizator - 3D Molecular Viewer")
        self.setGeometry(100, 100, 1400, 900)
        self.setStyleSheet("""
            QMainWindow { background-color: #2b2b2b; }
            QGroupBox { 
                font-weight: bold; 
                border: 2px solid #555; 
                border-radius: 5px; 
                margin-top: 1ex; 
                color: #ffffff;
            }
            QGroupBox::title { 
                subcontrol-origin: margin; 
                left: 10px; 
                padding: 0 5px 0 5px; 
            }
            QLabel { color: #ffffff; font-size: 11px; }
            QPushButton { 
                background-color: #404040; 
                color: white; 
                border: 1px solid #555; 
                padding: 5px; 
                border-radius: 3px; 
            }
            QPushButton:hover { background-color: #505050; }
            QPushButton:pressed { background-color: #606060; }
            QSlider::groove:horizontal {
                border: 1px solid #999999;
                height: 8px;
                background: qlineargradient(x1:0, y1:0, x2:0, y2:1, stop:0 #B1B1B1, stop:1 #c4c4c4);
                margin: 2px 0;
                border-radius: 4px;
            }
            QSlider::handle:horizontal {
                background: qlineargradient(x1:0, y1:0, x2:1, y2:1, stop:0 #b4b4b4, stop:1 #8f8f8f);
                border: 1px solid #5c5c5c;
                width: 18px;
                margin: -2px 0;
                border-radius: 3px;
            }
        """)
        
        # Create central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # Main layout
        main_layout = QHBoxLayout(central_widget)
        
        # Create VTK widget container for overlay
        vtk_container = QWidget()
        vtk_container_layout = QVBoxLayout(vtk_container)
        vtk_container_layout.setContentsMargins(0, 0, 0, 0)
        
        # VTK widget
        self.vtk_widget = VTK_QT_WIDGET(vtk_container)
        vtk_container_layout.addWidget(self.vtk_widget)
        
        # Create logo overlay
        self.logo_label = QLabel(vtk_container)
        try:
            pixmap = QPixmap("./pygamlab_gamvis.png")
            # Scale the logo 3 times bigger
            scaled_pixmap = pixmap.scaled(300, 150, Qt.KeepAspectRatio, Qt.SmoothTransformation)
            self.logo_label.setPixmap(scaled_pixmap)
        except Exception as e:
            # Fallback if image can't be loaded
            self.logo_label.setText("PyGAMLab")
            self.logo_label.setStyleSheet("""
                font-size: 14px; 
                font-weight: bold; 
                color: #333333;
                background-color: rgba(255, 255, 255, 200);
                padding: 5px;
                border-radius: 5px;
            """)
        
        # Style the logo with visible background (larger container)
        self.logo_label.setStyleSheet("""
            background-color: rgba(255, 255, 255, 220);
            border: 2px solid #cccccc;
            border-radius: 15px;
            padding: 20px;
            margin: 10px;
            min-width: 150px;
            min-height: 75px;
        """)
        
        # Position logo at top-left and ensure it stays visible
        self.logo_label.move(10, 10)
        self.logo_label.raise_()  # Bring to front
        self.logo_label.show()
        
        # Store reference to container for resize events
        self.vtk_container = vtk_container
        
        main_layout.addWidget(vtk_container, 3)
        
        # Control panel
        control_panel = self.create_control_panel()
        main_layout.addWidget(control_panel, 1)
        
    def create_control_panel(self):
        panel = QWidget()
        panel.setFixedWidth(300)
        layout = QVBoxLayout(panel)
        
        # File operations
        file_group = QGroupBox("File Operations")
        file_layout = QVBoxLayout(file_group)
        
        load_btn = QPushButton("Load Molecule")
        load_btn.clicked.connect(self.load_molecule)
        file_layout.addWidget(load_btn)
        
        save_btn = QPushButton("Save Snapshot")
        save_btn.clicked.connect(self.save_snapshot)
        file_layout.addWidget(save_btn)
        
        export_3d_btn = QPushButton("Export 3D Model")
        export_3d_btn.clicked.connect(self.export_3d_model)
        file_layout.addWidget(export_3d_btn)
        
        layout.addWidget(file_group)
        
        # Rendering controls
        render_group = QGroupBox("Rendering Controls")
        render_layout = QGridLayout(render_group)
        
        # Atom size slider
        render_layout.addWidget(QLabel("Atom Size:"), 0, 0)
        self.size_slider = QSlider(Qt.Horizontal)
        self.size_slider.setRange(10, 200)
        self.size_slider.setValue(100)
        self.size_slider.valueChanged.connect(self.update_atom_sizes)
        render_layout.addWidget(self.size_slider, 0, 1)
        self.size_label = QLabel("1.0x")
        render_layout.addWidget(self.size_label, 0, 2)
        
        # Background color
        render_layout.addWidget(QLabel("Background:"), 1, 0)
        self.bg_combo = QComboBox()
        self.bg_combo.addItems(["Black", "White", "Gray", "Blue", "Navy"])
        self.bg_combo.currentTextChanged.connect(self.change_background)
        render_layout.addWidget(self.bg_combo, 1, 1, 1, 2)
        
        # Render style
        render_layout.addWidget(QLabel("Style:"), 2, 0)
        self.style_combo = QComboBox()
        self.style_combo.addItems(["Spheres", "Points", "Wireframe"])
        self.style_combo.currentTextChanged.connect(self.change_render_style)
        render_layout.addWidget(self.style_combo, 2, 1, 1, 2)
        
        # Shadow control
        self.shadow_cb = QCheckBox("Enable Shadows")
        self.shadow_cb.setChecked(False)  # Default: no shadows
        self.shadow_cb.stateChanged.connect(self.toggle_shadows)
        render_layout.addWidget(self.shadow_cb, 3, 0, 1, 3)
        
        # Transparency control
        render_layout.addWidget(QLabel("Transparency:"), 4, 0)
        self.transparency_slider = QSlider(Qt.Horizontal)
        self.transparency_slider.setRange(0, 90)
        self.transparency_slider.setValue(0)
        self.transparency_slider.valueChanged.connect(self.update_transparency)
        render_layout.addWidget(self.transparency_slider, 4, 1)
        self.transparency_label = QLabel("0%")
        render_layout.addWidget(self.transparency_label, 4, 2)
        
        layout.addWidget(render_group)
        
        # Animation controls
        anim_group = QGroupBox("Animation Controls")
        anim_layout = QGridLayout(anim_group)
        
        # Animation toggle
        self.anim_btn = QPushButton("Start Animation")
        self.anim_btn.clicked.connect(self.toggle_animation)
        anim_layout.addWidget(self.anim_btn, 0, 0, 1, 3)
        
        # Vibration amplitude
        anim_layout.addWidget(QLabel("Vibration:"), 1, 0)
        self.vib_slider = QSlider(Qt.Horizontal)
        self.vib_slider.setRange(0, 100)
        self.vib_slider.setValue(20)
        self.vib_slider.valueChanged.connect(self.update_vibration)
        anim_layout.addWidget(self.vib_slider, 1, 1)
        self.vib_label = QLabel("0.05Å")
        anim_layout.addWidget(self.vib_label, 1, 2)
        
        # Animation speed
        anim_layout.addWidget(QLabel("Speed:"), 2, 0)
        self.speed_slider = QSlider(Qt.Horizontal)
        self.speed_slider.setRange(1, 100)
        self.speed_slider.setValue(50)
        self.speed_slider.valueChanged.connect(self.update_animation_speed)
        anim_layout.addWidget(self.speed_slider, 2, 1)
        self.speed_label = QLabel("50ms")
        anim_layout.addWidget(self.speed_label, 2, 2)
        
        layout.addWidget(anim_group)
        
        # Camera controls
        camera_group = QGroupBox("Camera Controls")
        camera_layout = QVBoxLayout(camera_group)
        
        reset_btn = QPushButton("Reset View")
        reset_btn.clicked.connect(self.reset_camera)
        camera_layout.addWidget(reset_btn)
        
        # Camera presets
        preset_layout = QHBoxLayout()
        preset_layout.addWidget(QLabel("Preset:"))
        self.preset_combo = QComboBox()
        self.preset_combo.addItems(["Front", "Back", "Top", "Bottom", "Left", "Right"])
        self.preset_combo.currentTextChanged.connect(self.set_camera_preset)
        preset_layout.addWidget(self.preset_combo)
        camera_layout.addLayout(preset_layout)
        
        # Auto-rotate
        self.auto_rotate_cb = QCheckBox("Auto Rotate")
        self.auto_rotate_cb.stateChanged.connect(self.toggle_auto_rotate)
        camera_layout.addWidget(self.auto_rotate_cb)
        
        layout.addWidget(camera_group)
        
        # Lighting controls
        lighting_group = QGroupBox("Lighting Controls")
        lighting_layout = QGridLayout(lighting_group)
        
        # Light intensity
        lighting_layout.addWidget(QLabel("Intensity:"), 0, 0)
        self.light_intensity_slider = QSlider(Qt.Horizontal)
        self.light_intensity_slider.setRange(10, 200)
        self.light_intensity_slider.setValue(100)
        self.light_intensity_slider.valueChanged.connect(self.update_light_intensity)
        lighting_layout.addWidget(self.light_intensity_slider, 0, 1)
        self.light_intensity_label = QLabel("1.0x")
        lighting_layout.addWidget(self.light_intensity_label, 0, 2)
        
        # Ambient lighting
        lighting_layout.addWidget(QLabel("Ambient:"), 1, 0)
        self.ambient_slider = QSlider(Qt.Horizontal)
        self.ambient_slider.setRange(0, 100)
        self.ambient_slider.setValue(100)  # Default to 1.0 (full ambient)
        self.ambient_slider.valueChanged.connect(self.update_ambient_light)
        lighting_layout.addWidget(self.ambient_slider, 1, 1)
        self.ambient_label = QLabel("1.0")
        lighting_layout.addWidget(self.ambient_label, 1, 2)
        
        layout.addWidget(lighting_group)
        
        # Display options
        display_group = QGroupBox("Display Options")
        display_layout = QVBoxLayout(display_group)
        
        # Atom labels
        self.labels_cb = QCheckBox("Show Atom Labels")
        self.labels_cb.setChecked(False)
        self.labels_cb.stateChanged.connect(self.toggle_atom_labels)
        display_layout.addWidget(self.labels_cb)
        
        # Axes
        self.axes_cb = QCheckBox("Show Coordinate Axes")
        self.axes_cb.setChecked(False)
        self.axes_cb.stateChanged.connect(self.toggle_axes)
        display_layout.addWidget(self.axes_cb)
        
        # Measurement mode
        self.measure_btn = QPushButton("Distance Measurement")
        self.measure_btn.setCheckable(True)
        self.measure_btn.clicked.connect(self.toggle_measurement_mode)
        display_layout.addWidget(self.measure_btn)
        
        layout.addWidget(display_group)
        
        # Information panel
        info_group = QGroupBox("Molecule Information")
        info_layout = QVBoxLayout(info_group)
        
        self.info_label = QLabel("No molecule loaded")
        self.info_label.setWordWrap(True)
        info_layout.addWidget(self.info_label)
        
        layout.addWidget(info_group)
        
        # Add stretch to push everything to top
        layout.addStretch()
        
        return panel
        
    def init_vtk(self):
        # Create renderer and render window
        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(0.0, 0.0, 0.0)  # Black background
        
        # Get render window
        render_window = self.vtk_widget.GetRenderWindow()
        render_window.AddRenderer(self.renderer)
        
        # Enable interaction
        self.interactor = self.vtk_widget.GetRenderWindow().GetInteractor()
        
        # Set interactor style for better camera control
        style = vtk.vtkInteractorStyleTrackballCamera()
        self.interactor.SetInteractorStyle(style)
        
        # Add lighting (no shadows by default)
        self.main_light = vtk.vtkLight()
        self.main_light.SetPosition(10, 10, 10)
        self.main_light.SetFocalPoint(0, 0, 0)
        self.main_light.SetIntensity(1.0)
        self.renderer.AddLight(self.main_light)
        
        # Disable shadows by default
        self.renderer.SetUseShadows(False)
        
        # Set ambient lighting to be minimal (let individual actors control their own ambient)
        self.renderer.SetAmbient(0.1, 0.1, 0.1)
        
        # Setup animation timer
        self.animation_timer.timeout.connect(self.animate_vibration)
    '''
    def setup_sample_molecule(self):
        """Create a sample molecule (caffeine-like structure) for demonstration"""
        sample_atoms = [
            GAM_Atom("C", 0.0, 0.0, 0.0),
            GAM_Atom("N", 1.4, 0.0, 0.0),
            GAM_Atom("C", 2.1, 1.2, 0.0),
            GAM_Atom("N", 1.4, 2.4, 0.0),
            GAM_Atom("C", 0.0, 2.4, 0.0),
            GAM_Atom("N", -0.7, 1.2, 0.0),
            GAM_Atom("O", 3.3, 1.2, 0.0),
            GAM_Atom("C", 2.1, -1.2, 0.0),
            GAM_Atom("H", -0.5, -0.9, 0.0),
            GAM_Atom("H", -0.5, 3.3, 0.0),
            GAM_Atom("H", 1.6, -2.1, 0.0),
            GAM_Atom("H", 2.6, -1.2, 0.9),
            GAM_Atom("H", 2.6, -1.2, -0.9),
        ]
        self.load_atoms(sample_atoms)
        
    '''
    
    def load_atoms(self, atoms):
        """Load atoms into the visualizer"""
        self.atoms = atoms
        self.original_positions = [np.array([atom.x, atom.y, atom.z]) for atom in atoms]
        self.render_molecules()
        self.update_info_panel()
        
    def render_molecules(self):
        #-----> here we must add the bonds between the atoms
        #we can have function for gettign teh bonds and tehn visualization
        """Render all atoms as spheres in the VTK scene"""
        # Clear existing actors
        for actor in self.sphere_actors:
            self.renderer.RemoveActor(actor)
        for actor in self.label_actors:
            self.renderer.RemoveActor(actor)
        self.sphere_actors.clear()
        self.label_actors.clear()
        
        for i, atom in enumerate(self.atoms):
            # Create sphere
            sphere = vtk.vtkSphereSource()
            radius = ELEMENT_RADII.get(atom.element, 1.5) * 0.3 * (self.size_slider.value() / 100.0)
            sphere.SetRadius(radius)
            sphere.SetPhiResolution(20)
            sphere.SetThetaResolution(20)
            
            # Create mapper
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(sphere.GetOutputPort())
            
            # Create actor
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            
            # Set position
            actor.SetPosition(atom.x, atom.y, atom.z)
            
            # Set color based on element
            color = ELEMENT_COLORS.get(atom.element, (0.8, 0.8, 0.8))
            actor.GetProperty().SetColor(color)
            actor.GetProperty().SetSpecular(0.3)
            actor.GetProperty().SetSpecularPower(20)
            actor.GetProperty().SetAmbient(1.0)  # Full ambient lighting by default
            
            # Add to renderer
            self.renderer.AddActor(actor)
            self.sphere_actors.append(actor)
            
        # Render labels if enabled
        if self.labels_cb.isChecked():
            self.render_labels()
            
        # Fit camera to show all atoms
        self.renderer.ResetCamera()
        self.vtk_widget.GetRenderWindow().Render()
        
    def update_atom_sizes(self, value):
        """Update atom sizes based on slider value"""
        scale = value / 100.0
        self.size_label.setText(f"{scale:.1f}x")
        
        for i, actor in enumerate(self.sphere_actors):
            atom = self.atoms[i]
            radius = ELEMENT_RADII.get(atom.element, 1.5) * 0.3 * scale
            
            # Get the sphere source and update radius
            mapper = actor.GetMapper()
            sphere_source = mapper.GetInputConnection(0, 0).GetProducer()
            sphere_source.SetRadius(radius)
            sphere_source.Update()
            
        self.vtk_widget.GetRenderWindow().Render()
        
    def change_background(self, color_name):
        """Change background color"""
        colors = {
            "Black": (0.0, 0.0, 0.0),
            "White": (1.0, 1.0, 1.0),
            "Gray": (0.5, 0.5, 0.5),
            "Blue": (0.2, 0.4, 0.8),
            "Navy": (0.0, 0.1, 0.3)
        }
        color = colors.get(color_name, (0.0, 0.0, 0.0))
        self.renderer.SetBackground(color)
        self.vtk_widget.GetRenderWindow().Render()
        
    def change_render_style(self, style):
        """Change rendering style"""
        for actor in self.sphere_actors:
            if style == "Spheres":
                actor.GetProperty().SetRepresentationToSurface()
            elif style == "Points":
                actor.GetProperty().SetRepresentationToPoints()
                actor.GetProperty().SetPointSize(5)
            elif style == "Wireframe":
                actor.GetProperty().SetRepresentationToWireframe()
        
        self.vtk_widget.GetRenderWindow().Render()
        
    def toggle_animation(self):
        """Start or stop vibration animation"""
        if self.animation_timer.isActive():
            self.animation_timer.stop()
            self.anim_btn.setText("Start Animation")
        else:
            self.animation_timer.start(50)  # 50ms interval
            self.anim_btn.setText("Stop Animation")
            
    def update_vibration(self, value):
        """Update vibration amplitude"""
        self.vibration_amplitude = value / 1000.0  # Convert to Angstroms
        self.vib_label.setText(f"{self.vibration_amplitude:.3f}Å")
        
    def update_animation_speed(self, value):
        """Update animation speed"""
        interval = 101 - value  # Reverse scale (higher value = faster)
        self.speed_label.setText(f"{interval}ms")
        if self.animation_timer.isActive():
            self.animation_timer.setInterval(interval)
            
    def animate_vibration(self):
        """Animate tiny vibrations for each atom"""
        self.animation_frame += 1
        
        for i, (actor, original_pos) in enumerate(zip(self.sphere_actors, self.original_positions)):
            # Generate small random vibrations
            vibration = np.random.normal(0, self.vibration_amplitude, 3)
            new_pos = original_pos + vibration
            actor.SetPosition(new_pos[0], new_pos[1], new_pos[2])
            
        self.vtk_widget.GetRenderWindow().Render()
        
    def reset_camera(self):
        """Reset camera to show all molecules"""
        self.renderer.ResetCamera()
        self.vtk_widget.GetRenderWindow().Render()
        
    def set_camera_preset(self, preset):
        """Set camera to predefined viewpoint"""
        camera = self.renderer.GetActiveCamera()
        
        if preset == "Front":
            camera.SetPosition(0, 0, 10)
            camera.SetViewUp(0, 1, 0)
        elif preset == "Back":
            camera.SetPosition(0, 0, -10)
            camera.SetViewUp(0, 1, 0)
        elif preset == "Top":
            camera.SetPosition(0, 10, 0)
            camera.SetViewUp(0, 0, -1)
        elif preset == "Bottom":
            camera.SetPosition(0, -10, 0)
            camera.SetViewUp(0, 0, 1)
        elif preset == "Left":
            camera.SetPosition(-10, 0, 0)
            camera.SetViewUp(0, 1, 0)
        elif preset == "Right":
            camera.SetPosition(10, 0, 0)
            camera.SetViewUp(0, 1, 0)
            
        camera.SetFocalPoint(0, 0, 0)
        self.renderer.ResetCamera()
        self.vtk_widget.GetRenderWindow().Render()
        
    def toggle_auto_rotate(self, state):
        """Toggle automatic camera rotation"""
        # Implementation would require a separate timer for rotation
        pass  # Placeholder for now
        
    def update_info_panel(self):
        """Update the information panel with molecule stats"""
        if not self.atoms:
            self.info_label.setText("No molecule loaded")
            return
            
        element_counts = {}
        for atom in self.atoms:
            element_counts[atom.element] = element_counts.get(atom.element, 0) + 1
            
        info_text = f"Total atoms: {len(self.atoms)}\n\n"
        info_text += "Element composition:\n"
        for element, count in sorted(element_counts.items()):
            info_text += f"{element}: {count}\n"
            
        self.info_label.setText(info_text)
        
    def load_molecule(self):
        """Load molecule from file (placeholder)"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Load Molecule", "", "XYZ files (*.xyz);;PDB files (*.pdb);;All files (*)"
        )
        if file_path:
            # Placeholder - would implement actual file parsing here
            QMessageBox.information(self, "Info", f"Would load molecule from:\n{file_path}")
            
    def save_snapshot(self):
        """Save current view as image"""
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Snapshot", "molecule_snapshot.png", "PNG files (*.png);;JPG files (*.jpg)"
        )
        if file_path:
            # Create window to image filter
            w2if = vtk.vtkWindowToImageFilter()
            w2if.SetInput(self.vtk_widget.GetRenderWindow())
            w2if.Update()
            
            # Write image
            writer = vtk.vtkPNGWriter()
            writer.SetFileName(file_path)
            writer.SetInputConnection(w2if.GetOutputPort())
            writer.Write()
            
            QMessageBox.information(self, "Success", f"Snapshot saved to:\n{file_path}")

    def export_3d_model(self):
        """Export 3D model for 3D printing or external use"""
        file_path, file_filter = QFileDialog.getSaveFileName(
            self, "Export 3D Model", "molecule_model.obj", 
            "OBJ files (*.obj);;STL files (*.stl);;PLY files (*.ply)"
        )
        if file_path:
            try:
                # Create a combined geometry of all visible actors
                append_filter = vtk.vtkAppendPolyData()
                
                # Add all sphere geometries
                for actor in self.sphere_actors:
                    mapper = actor.GetMapper()
                    sphere_source = mapper.GetInputConnection(0, 0).GetProducer()
                    
                    # Transform the geometry to world coordinates
                    transform = vtk.vtkTransform()
                    transform.SetMatrix(actor.GetMatrix())
                    
                    transform_filter = vtk.vtkTransformPolyDataFilter()
                    transform_filter.SetInputConnection(sphere_source.GetOutputPort())
                    transform_filter.SetTransform(transform)
                    transform_filter.Update()
                    
                    append_filter.AddInputConnection(transform_filter.GetOutputPort())
                
                append_filter.Update()
                
                # Choose writer based on file extension
                if file_path.lower().endswith('.obj'):
                    writer = vtk.vtkOBJWriter()
                elif file_path.lower().endswith('.stl'):
                    writer = vtk.vtkSTLWriter()
                elif file_path.lower().endswith('.ply'):
                    writer = vtk.vtkPLYWriter()
                else:
                    QMessageBox.warning(self, "Error", "Unsupported file format")
                    return
                
                writer.SetFileName(file_path)
                writer.SetInputConnection(append_filter.GetOutputPort())
                writer.Write()
                
                QMessageBox.information(self, "Success", f"3D model exported to:\n{file_path}")
                
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to export 3D model:\n{str(e)}")

    def toggle_shadows(self, state):
        """Toggle shadow rendering"""
        enable_shadows = state == Qt.Checked
        self.renderer.SetUseShadows(enable_shadows)
        
        # Update all sphere actors for shadow support
        for actor in self.sphere_actors:
            if enable_shadows:
                # Enable shadow mapping for better quality
                actor.GetProperty().SetLighting(True)
            else:
                actor.GetProperty().SetLighting(True)
        
        self.vtk_widget.GetRenderWindow().Render()
        
    def update_transparency(self, value):
        """Update atom transparency"""
        transparency = value / 100.0
        self.transparency_label.setText(f"{value}%")
        
        for actor in self.sphere_actors:
            actor.GetProperty().SetOpacity(1.0 - transparency)
            
        self.vtk_widget.GetRenderWindow().Render()
        
    def update_light_intensity(self, value):
        """Update main light intensity"""
        intensity = value / 100.0
        self.light_intensity_label.setText(f"{intensity:.1f}x")
        if self.main_light:
            self.main_light.SetIntensity(intensity)
        self.vtk_widget.GetRenderWindow().Render()
        
    def update_ambient_light(self, value):
        """Update ambient lighting"""
        ambient = value / 100.0
        self.ambient_label.setText(f"{ambient:.1f}")
        
        for actor in self.sphere_actors:
            actor.GetProperty().SetAmbient(ambient)
            
        self.vtk_widget.GetRenderWindow().Render()
                    
    def render_labels(self):
        """Render atom labels"""
        for i, atom in enumerate(self.atoms):
            # Create text source
            text_source = vtk.vtkVectorText()
            text_source.SetText(f"{atom.element}{i+1}")
            
            # Create mapper
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(text_source.GetOutputPort())
            
            # Create actor
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            
            # Position above the atom
            offset_radius = ELEMENT_RADII.get(atom.element, 1.5) * 0.3 * (self.size_slider.value() / 100.0)
            actor.SetPosition(atom.x, atom.y + offset_radius + 0.5, atom.z)
            actor.SetScale(0.3, 0.3, 0.3)
            
            # Set color (white for visibility)
            actor.GetProperty().SetColor(1.0, 1.0, 1.0)
            
            self.renderer.AddActor(actor)
            self.label_actors.append(actor)
            
    def toggle_atom_labels(self, state):
        """Toggle atom label visibility"""
        # Remove existing labels
        for actor in self.label_actors:
            self.renderer.RemoveActor(actor)
        self.label_actors.clear()
        
        # Re-render if enabled
        if state == Qt.Checked:
            self.render_labels()
            
        self.vtk_widget.GetRenderWindow().Render()
        
    def toggle_axes(self, state):
        """Toggle coordinate axes visibility"""
        if state == Qt.Checked:
            if self.axes_actor is None:
                # Create axes
                axes = vtk.vtkAxesActor()
                axes.SetTotalLength(2.0, 2.0, 2.0)
                axes.SetShaftType(0)  # Line shaft
                axes.SetTipType(0)    # Cone tip
                
                self.axes_actor = axes
                
            self.renderer.AddActor(self.axes_actor)
        else:
            if self.axes_actor:
                self.renderer.RemoveActor(self.axes_actor)
                
        self.vtk_widget.GetRenderWindow().Render()
        
    def toggle_measurement_mode(self, checked):
        """Toggle distance measurement mode"""
        self.measurement_mode = checked
        if checked:
            self.measure_btn.setText("Exit Measurement")
            self.measurement_points.clear()
            # Add click handler for measurements
            # This would require more complex VTK interaction setup
        else:
            self.measure_btn.setText("Distance Measurement")
            self.measurement_points.clear()
            
    def resizeEvent(self, event):
        """Handle window resize to keep logo positioned correctly"""
        super().resizeEvent(event)
        if hasattr(self, 'logo_label') and hasattr(self, 'vtk_container'):
            # Keep logo at top-left corner
            self.logo_label.move(10, 10)
            self.logo_label.raise_()
            
  






