#!/usr/bin/env python3

import os
import yaml
import json
import sys
import subprocess
import numpy as np
import meshio
import re

from scipy.spatial.transform import Rotation as R

from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QFileDialog, QLineEdit, QMessageBox, QCheckBox, QMainWindow, QFrame,
    QSplitter, QTextEdit, QGroupBox
)
from PyQt5.QtGui import QFont, QTextCharFormat, QColor, QSyntaxHighlighter
from PyQt5.QtCore import Qt, QRegExp

import pyvista as pv
from pyvistaqt import QtInteractor

from typing import Any, Dict, List, Union

def has_nvidia_gpu():
    try:
        subprocess.run(["nvidia-smi"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

def load_mesh(filename):
    try:
        mesh = meshio.read(filename)
        points = mesh.points[:, :3]
        for cell_block in mesh.cells:
            if cell_block.type == "triangle":
                triangles = cell_block.data
                break
        else:
            raise ValueError("No triangle cells found in the mesh.")
        return points, triangles
    except:
        # If the file does not exist, create a simple cylinder
        theta = np.linspace(0, 2*np.pi, 20)
        z = np.linspace(0, 1, 10)
        theta_grid, z_grid = np.meshgrid(theta, z)
        x = np.cos(theta_grid).flatten()
        y = np.sin(theta_grid).flatten()
        z = z_grid.flatten()
        points = np.column_stack([x, y, z])

        # Create simple triangles
        triangles = []
        for i in range(len(theta)-1):
            for j in range(len(z)-1):
                idx = i * len(z) + j
                triangles.append([idx, idx+1, idx+len(z)])
                triangles.append([idx+1, idx+len(z)+1, idx+len(z)])

        return points, np.array(triangles)

def normalize(v):
    norm = np.linalg.norm(v)
    return v / norm if norm != 0 else np.array([1.0, 0.0, 0.0])

def apply_sequential_quaternion_rotations(points, axes_angles):
    rotated_points = points.copy()
    rotation = R.identity()
    for axis, angle in axes_angles:
        axis = normalize(np.array(axis))
        r = R.from_rotvec(np.deg2rad(angle) * axis)
        rotation = r * rotation
    return rotation.apply(rotated_points)

def create_pyvista_mesh(points, triangles):
    faces = np.hstack([np.full((triangles.shape[0], 1), 3), triangles]).astype(np.int32)
    return pv.PolyData(points, faces)

def clean_yaml_content(content: str) -> str:
    lines = content.split('\n')
    cleaned_lines = [line.replace('\t', '  ') for line in lines]
    return '\n'.join(cleaned_lines)

def load_yaml_from_file(filename: str) -> Dict[str, Any]:
    try:
        with open(filename, 'r', encoding='utf-8') as file:
            content = file.read()
        cleaned_content = clean_yaml_content(content)
        config = yaml.safe_load(cleaned_content)
        return config
    except FileNotFoundError:
        print(f"❌ File '{filename}' not found")
        return None
    except Exception as e:
        print(f"❌ Error reading file: {e}")
        return None

def save_yaml_file(data: dict, filename: str):
    with open(filename, 'w') as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)

class YAMLSyntaxHighlighter(QSyntaxHighlighter):
    def __init__(self, parent=None):
        super().__init__(parent)
        
        # Format pour les clés (en noir)
        self.key_format = QTextCharFormat()
        self.key_format.setForeground(QColor("black"))
        self.key_format.setFontWeight(QFont.Bold)
        
        # Format pour les valeurs (en rouge)
        self.value_format = QTextCharFormat()
        self.value_format.setForeground(QColor("red"))
        
        # Format pour les commentaires (en vert)
        self.comment_format = QTextCharFormat()
        self.comment_format.setForeground(QColor("green"))
        
    def highlightBlock(self, text):
        # Commentaires
        comment_pattern = QRegExp(r'#.*')
        index = comment_pattern.indexIn(text)
        while index >= 0:
            length = comment_pattern.matchedLength()
            self.setFormat(index, length, self.comment_format)
            index = comment_pattern.indexIn(text, index + length)
        
        # Clés YAML (mot suivi de deux points)
        key_pattern = QRegExp(r'^(\s*)([^:\s]+)(\s*):')
        index = key_pattern.indexIn(text)
        while index >= 0:
            key_start = key_pattern.pos(2)
            key_length = len(key_pattern.cap(2))
            self.setFormat(key_start, key_length, self.key_format)
            index = key_pattern.indexIn(text, index + key_pattern.matchedLength())
        
        # Valeurs YAML (tout ce qui suit les deux points)
        value_pattern = QRegExp(r':\s*(.+)$')
        index = value_pattern.indexIn(text)
        while index >= 0:
            value_start = value_pattern.pos(1)
            value_length = len(value_pattern.cap(1))
            if value_start >= 0:  # Vérifier que la valeur existe
                self.setFormat(value_start, value_length, self.value_format)
            index = value_pattern.indexIn(text, index + value_pattern.matchedLength())

class CombinedInterface(QMainWindow):
    def __init__(self, yaml_file="settings.yml"):
        super().__init__()
        self.setWindowTitle("YAML Editor + 3D Visualization")
        
        # Fenêtre en plein écran
        self.showMaximized()

        # Variables
        self.file_path = yaml_file
        self.yaml_config = {}
        self.mesh_filename = None

        self.setup_ui()
        self.load_yaml_file()

    def setup_ui(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        # Main splitter
        main_splitter = QSplitter(Qt.Horizontal)
        central_widget.setLayout(QHBoxLayout())
        central_widget.layout().addWidget(main_splitter)

        # === Left Panel: YAML Editor ===
        left_panel = QFrame()
        left_panel.setFrameShape(QFrame.StyledPanel)
        left_panel.setMinimumWidth(500)
        left_layout = QVBoxLayout(left_panel)

                
        sol_layout = QVBoxLayout()
        self.sol_label = QLabel(".sol File:")
        self.sol_input = QLineEdit()
        self.sol_btn = QPushButton("Browse")
        self.sol_btn.clicked.connect(self.select_sol_file)

        sol_layout.addWidget(self.sol_label)
        sol_layout.addWidget(self.sol_input)
        sol_layout.addWidget(self.sol_btn)
        left_layout.addLayout(sol_layout)

        # Rotation parameters
        rotation_group = QGroupBox("Rotation Parameters")
        rotation_layout = QVBoxLayout(rotation_group)

        # Angle 1 and Axis 1
        angle1_layout = QVBoxLayout()
        angle1_label = QLabel("Angle 1 (°):")
        self.angle1_input = QLineEdit("90")
        self.angle1_input.textChanged.connect(self.update_mesh)
        angle1_layout.addWidget(angle1_label)
        angle1_layout.addWidget(self.angle1_input)
        rotation_layout.addLayout(angle1_layout)

        # Axis 1 layout (with horizontal arrangement for x, y, z)
        axe1_layout = QVBoxLayout()
        axe1_layout.addWidget(QLabel("Axis 1 (x y z):"))

        axis_input_layout = QHBoxLayout()
        self.axe1_inputs = [QLineEdit("1"), QLineEdit("0"), QLineEdit("0")]
        for inp in self.axe1_inputs:
            inp.textChanged.connect(self.update_mesh)
            axis_input_layout.addWidget(inp)

        axe1_layout.addLayout(axis_input_layout)
        rotation_layout.addLayout(axe1_layout)

        # Angle 2 and Axis 2
        angle2_layout = QVBoxLayout()
        angle2_label = QLabel("Angle 2 (°):")
        self.angle2_input = QLineEdit("0")
        self.angle2_input.textChanged.connect(self.update_mesh)
        angle2_layout.addWidget(angle2_label)
        angle2_layout.addWidget(self.angle2_input)
        rotation_layout.addLayout(angle2_layout)

        # Axis 2 layout (with horizontal arrangement for x, y, z)
        axe2_layout = QVBoxLayout()
        axe2_layout.addWidget(QLabel("Axis 2 (x y z):"))

        axis_input_layout = QHBoxLayout()
        self.axe2_inputs = [QLineEdit("1"), QLineEdit("0"), QLineEdit("0")]
        for inp in self.axe2_inputs:
            inp.textChanged.connect(self.update_mesh)
            axis_input_layout.addWidget(inp)

        axe2_layout.addLayout(axis_input_layout)
        rotation_layout.addLayout(axe2_layout)
        
        left_layout.addWidget(rotation_group)

        # Other parameters
        other_group = QGroupBox("Other Parameters")
        other_layout = QVBoxLayout(other_group)

        self.filled_checkbox = QCheckBox("Filled")
        self.filled_checkbox.setChecked(True)
        other_layout.addWidget(self.filled_checkbox)

        # Electrostatics
        elec_layout = QVBoxLayout()
        ce_label = QLabel("CE:")
        self.ce_input = QLineEdit("0")
        v_label = QLabel("V:")
        self.v_input = QLineEdit("0.0")
        elec_layout.addWidget(ce_label)
        elec_layout.addWidget(self.ce_input)
        elec_layout.addWidget(v_label)
        elec_layout.addWidget(self.v_input)
        other_layout.addLayout(elec_layout)

        # Detector
        det_layout = QVBoxLayout()
        zoom_label = QLabel("Zoom:")
        self.zoom_input = QLineEdit("1.0")
        meshsize_label = QLabel("Mesh Size:")
        self.meshsize_input = QLineEdit("1e-9")
        det_layout.addWidget(zoom_label)
        det_layout.addWidget(self.zoom_input)
        det_layout.addWidget(meshsize_label)
        det_layout.addWidget(self.meshsize_input)
        other_layout.addLayout(det_layout)

        left_layout.addWidget(other_group)

        # Action buttons
        button_layout = QHBoxLayout()

        self.save_yaml_btn = QPushButton("Save YAML")
        self.save_yaml_btn.clicked.connect(self.save_yaml_config)

        button_layout.addWidget(self.save_yaml_btn)
        left_layout.addLayout(button_layout)

        # Text area to display YAML - prend toute la hauteur restante
        yaml_label = QLabel("YAML Configuration:")
        left_layout.addWidget(yaml_label)
        
        self.yaml_display = QTextEdit()
        self.yaml_display.setReadOnly(True)
        
        # Appliquer le highlighter pour la coloration syntaxique
        self.highlighter = YAMLSyntaxHighlighter(self.yaml_display.document())
        
        # La zone de texte prend toute la hauteur restante
        left_layout.addWidget(self.yaml_display, 1)  # Le paramètre 1 fait que ce widget s'étend

        # === Right Panel: 3D Visualization ===
        right_panel = QFrame()
        right_panel.setFrameShape(QFrame.StyledPanel)
        right_layout = QVBoxLayout(right_panel)

        # Title
        title_label = QLabel("3D Visualization - Rotations")
        title_label.setFont(QFont("Arial", 14, QFont.Bold))
        title_label.setAlignment(Qt.AlignCenter)
        right_layout.addWidget(title_label)

        # 3D Visualizer
        self.plotter = QtInteractor()
        self.plotter.set_background("white")
        right_layout.addWidget(self.plotter.interactor)

        # Add panels to the splitter
        main_splitter.addWidget(left_panel)
        main_splitter.addWidget(right_panel)
        main_splitter.setSizes([500, 900])

    def select_yaml_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select a YAML file", "", "YAML files (*.yml *.yaml)")
        if path:
            self.file_path = path
            self.file_label.setText(os.path.basename(path))
            self.load_yaml_file()

    def select_sol_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select a .sol file", "", "SOL files (*.sol)")
        if path:
            self.sol_file_path = path
            self.sol_input.setText(path)

    def load_yaml_file(self):
        if self.file_path:
            self.yaml_config = load_yaml_from_file(self.file_path)
            if self.yaml_config:
                self.update_yaml_display()
                self.populate_fields_from_yaml()
                self.load_mesh_from_config()

    def load_mesh_from_config(self):
        if 'mesh' in self.yaml_config and 'filename' in self.yaml_config['mesh']:
            self.mesh_filename = self.yaml_config['mesh']['filename']
            try:
                self.original_points, self.triangles = load_mesh(self.mesh_filename)
                self.update_mesh()
            except Exception as e:
                print(f"Error loading mesh: {e}")

    def populate_fields_from_yaml(self):
        if not self.yaml_config:
            return

        # Fill fields from YAML
        if 'initial_magnetization' in self.yaml_config:
            self.sol_input.setText(str(self.yaml_config['initial_magnetization']))
            self.sol_input.setText('No sol file selected')
			
        if 'rotations' in self.yaml_config:
            rot = self.yaml_config['rotations']
            if 'angle1' in rot:
                self.angle1_input.setText(str(rot['angle1']))
            if 'axe1' in rot:
                for i, val in enumerate(rot['axe1']):
                    if i < len(self.axe1_inputs):
                        self.axe1_inputs[i].setText(str(val))
            if 'angle2' in rot:
                self.angle2_input.setText(str(rot['angle2']))
            if 'axe2' in rot:
                for i, val in enumerate(rot['axe2']):
                    if i < len(self.axe2_inputs):
                        self.axe2_inputs[i].setText(str(val))

        if 'filled' in self.yaml_config:
            self.filled_checkbox.setChecked(str(self.yaml_config['filled']).lower() == 'true')

        if 'electrostatics' in self.yaml_config:
            elec = self.yaml_config['electrostatics']
            if 'CE' in elec:
                self.ce_input.setText(str(elec['CE']))
            if 'V' in elec:
                self.v_input.setText(str(elec['V']))

        if 'detector' in self.yaml_config:
            det = self.yaml_config['detector']
            if 'zoom' in det:
                self.zoom_input.setText(str(det['zoom']))
            if 'meshSize' in det:
                self.meshsize_input.setText(str(det['meshSize']))

    def update_yaml_display(self):
        if self.yaml_config:
            yaml_text = yaml.dump(self.yaml_config, default_flow_style=False, sort_keys=False)
            self.yaml_display.setPlainText(yaml_text)

    def save_yaml_config(self):
        # Update YAML configuration with field values
        try:
            if not self.yaml_config:
                self.yaml_config = {}

            # Initial magnetization
 #           if self.sol_input.text():
 #               self.yaml_config['initial_magnetization'] = os.path.basename(self.sol_input.text())

            if not self.sol_file_path:
                raise ValueError("No .sol file selected")
                
            # Rotations
            if 'rotations' not in self.yaml_config:
                self.yaml_config['rotations'] = {}

            self.yaml_config['rotations']['angle1'] = float(self.angle1_input.text())
            self.yaml_config['rotations']['axe1'] = [float(inp.text()) for inp in self.axe1_inputs]
            self.yaml_config['rotations']['angle2'] = float(self.angle2_input.text())
            self.yaml_config['rotations']['axe2'] = [float(inp.text()) for inp in self.axe2_inputs]

            # Filled
            self.yaml_config['filled'] = str(self.filled_checkbox.isChecked()).lower()

            # Electrostatics
            if 'electrostatics' not in self.yaml_config:
                self.yaml_config['electrostatics'] = {}
            self.yaml_config['electrostatics']['CE'] = int(self.ce_input.text())
            self.yaml_config['electrostatics']['V'] = float(self.v_input.text())

            # Detector
            if 'detector' not in self.yaml_config:
                self.yaml_config['detector'] = {}
            self.yaml_config['detector']['zoom'] = float(self.zoom_input.text())
            self.yaml_config['detector']['meshSize'] = float(self.meshsize_input.text())

            self.update_yaml_display()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error updating: {e}")
            return

        # Save the file
        if self.file_path:
            base, ext = os.path.splitext(self.file_path)
            filename = f"{base}_ray{ext}"
        else:
            filename, _ = QFileDialog.getSaveFileName(self, "Save YAML file", "", "YAML files (*.yml *.yaml)")
            if not filename:
                return

        try:
            save_yaml_file(self.yaml_config, filename)
            QMessageBox.information(self, "Success", f"File saved: {filename}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error saving: {e}")

    def get_rotations_from_inputs(self):
        try:
            axes_angles = []

            # First rotation
            angle1 = float(self.angle1_input.text())
            axis1 = [float(inp.text()) for inp in self.axe1_inputs]
            axes_angles.append((axis1, angle1))

            # Second rotation
            angle2 = float(self.angle2_input.text())
            axis2 = [float(inp.text()) for inp in self.axe2_inputs]
            axes_angles.append((axis2, angle2))

            return axes_angles
        except ValueError:
            return [([1, 0, 0], 0), ([0, 1, 0], 0)]

    def update_mesh(self):
        if not hasattr(self, 'original_points') or not hasattr(self, 'triangles'):
            return

        try:
            axes_angles = self.get_rotations_from_inputs()
            rotated_points = apply_sequential_quaternion_rotations(self.original_points, axes_angles)
            mesh = create_pyvista_mesh(rotated_points, self.triangles)

            self.plotter.clear()
            self.plotter.add_mesh(mesh, color="lightblue", show_edges=False)
            self.plotter.add_axes()
            self.plotter.show_grid()

            center = mesh.center
            camera_position = [center[0], center[1], center[2] + 2]
            focal_point = center
            view_up = [0, 1, 0]
            self.plotter.camera_position = [camera_position, focal_point, view_up]
            self.plotter.reset_camera()
        except Exception as e:
            print(f"Error updating mesh: {e}")


if __name__ == "__main__":
    if has_nvidia_gpu():
        os.environ["__NV_PRIME_RENDER_OFFLOAD"] = "1"
        os.environ["__GLX_VENDOR_LIBRARY_NAME"] = "nvidia"
        print("NVIDIA detected, environment variables set.")
    else:
        print("No NVIDIA GPU detected.")

    yaml_file = sys.argv[1] if len(sys.argv) > 1 else "settings.yml"

    app = QApplication(sys.argv)
    window = CombinedInterface(yaml_file)
    window.show()
    sys.exit(app.exec_())
