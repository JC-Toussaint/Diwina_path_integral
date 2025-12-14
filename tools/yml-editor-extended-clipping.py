#!/usr/bin/env python3

import os
import yaml
import sys
import subprocess
import numpy as np
import meshio

from scipy.spatial.transform import Rotation as R

from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QFileDialog, QLineEdit, QMessageBox, QCheckBox, QMainWindow, QFrame,
    QSplitter, QTextEdit, QGroupBox, QDoubleSpinBox, QComboBox
)
from PyQt5.QtGui import QFont, QTextCharFormat, QColor, QSyntaxHighlighter
from PyQt5.QtCore import Qt, QRegExp, QLocale

import pyvista as pv
from pyvistaqt import QtInteractor

from typing import Any, Dict, List, Union

def has_nvidia_gpu():
    try:
        subprocess.run(["nvidia-smi"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

def load_sol_file(filename):
    """Load magnetization vectors from a .sol file."""
    magnetization = {}
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        idx = int(parts[0])
                        mx = float(parts[1])
                        my = float(parts[2])
                        mz = float(parts[3])
                        magnetization[idx] = np.array([mx, my, mz])
                    except (ValueError, IndexError):
                        continue
        return magnetization
    except Exception as e:
        raise RuntimeError(f"Error loading .sol file '{filename}': {e}")

def get_surface_nodes(mesh):
    """Extract surface nodes from the mesh."""
    surface_nodes = set()
    for cell_block in mesh.cells:
        if cell_block.type == "triangle":
            for tri in cell_block.data:
                surface_nodes.update(tri)
    return surface_nodes

def get_region_elements(mesh):
    """Return a dictionary mapping region name -> list of element descriptions."""
    if not hasattr(mesh, "cells"):
        raise ValueError("mesh has no attribute 'cells'")
    
    id_to_name = {}
    for name, arr in mesh.field_data.items():
        try:
            phys_id = int(np.asarray(arr).flatten()[0])
        except Exception:
            continue
        id_to_name[phys_id] = name
    
    phys_key = None
    if hasattr(mesh, "cell_data") and mesh.cell_data:
        for k in mesh.cell_data.keys():
            if k == "gmsh:physical":
                phys_key = k
                break
        if phys_key is None:
            for k in mesh.cell_data.keys():
                if "physical" in k.lower():
                    phys_key = k
                    break
    
    if phys_key is None:
        raise ValueError("No physical tags found in mesh.cell_data.")
    
    regions = {}
    global_idx = 0
    phys_tags_list = mesh.cell_data[phys_key]
    
    if len(phys_tags_list) != len(mesh.cells):
        raise ValueError("Mismatch between number of cell blocks and cell_data arrays.")
    
    for block_index, (cell_block, phys_tags) in enumerate(zip(mesh.cells, phys_tags_list)):
        cell_type = cell_block.type
        connectivity = np.asarray(cell_block.data, dtype=int)
        phys_tags = np.asarray(phys_tags, dtype=int).flatten()
        
        if connectivity.shape[0] != phys_tags.shape[0]:
            raise ValueError(f"Number of elements mismatch in block {block_index}.")
        
        for local_idx in range(connectivity.shape[0]):
            phys_id = int(phys_tags[local_idx])
            region_name = id_to_name.get(phys_id, f"phys_{phys_id}")
            elem = {
                "cell_type": cell_type,
                "block_index": block_index,
                "local_element_index": local_idx,
                "global_element_index": global_idx,
                "nodes": connectivity[local_idx].tolist()
            }
            regions.setdefault(region_name, []).append(elem)
            global_idx += 1
    
    return regions

def load_mesh(filename):
    """Load a mesh from a file using meshio."""
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"The file '{filename}' does not exist.")
    
    try:
        mesh = meshio.read(filename)
        for cell_block, cell_data in zip(mesh.cells, mesh.cell_data["gmsh:physical"]):
            cell_type = cell_block.type
            print(f"\nCell type: {cell_type}")
            for region_name, (region_id, dim) in mesh.field_data.items():
                element_indices = np.where(cell_data == region_id)[0]
                if len(element_indices) > 0:
                    print(f"  Region '{region_name}' (ID={region_id}, dim={dim}): {len(element_indices)} elements")
        
        regions = get_region_elements(mesh)
        points = mesh.points[:, :3]
        
        triangles = []
        for cell_block in mesh.cells:
            if cell_block.type == "triangle":
                triangles.append(cell_block.data)
        
        if not triangles:
            raise ValueError(f"No triangular cells found in mesh '{filename}'.")
        
        triangles = np.vstack(triangles)
        return mesh, points, triangles
    
    except Exception as e:
        raise RuntimeError(f"Error while loading mesh '{filename}': {e}")

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
        
        self.key_format = QTextCharFormat()
        self.key_format.setForeground(QColor("black"))
        self.key_format.setFontWeight(QFont.Bold)
        
        self.value_format = QTextCharFormat()
        self.value_format.setForeground(QColor("red"))
        
        self.comment_format = QTextCharFormat()
        self.comment_format.setForeground(QColor("green"))
        
    def highlightBlock(self, text):
        comment_pattern = QRegExp(r'#.*')
        index = comment_pattern.indexIn(text)
        while index >= 0:
            length = comment_pattern.matchedLength()
            self.setFormat(index, length, self.comment_format)
            index = comment_pattern.indexIn(text, index + length)
        
        key_pattern = QRegExp(r'^(\s*)([^:\s]+)(\s*):')
        index = key_pattern.indexIn(text)
        while index >= 0:
            key_start = key_pattern.pos(2)
            key_length = len(key_pattern.cap(2))
            self.setFormat(key_start, key_length, self.key_format)
            index = key_pattern.indexIn(text, index + key_pattern.matchedLength())
        
        value_pattern = QRegExp(r':\s*(.+)$')
        index = value_pattern.indexIn(text)
        while index >= 0:
            value_start = value_pattern.pos(1)
            value_length = len(value_pattern.cap(1))
            if value_start >= 0:
                self.setFormat(value_start, value_length, self.value_format)
            index = value_pattern.indexIn(text, index + value_pattern.matchedLength())


class CombinedInterface(QMainWindow):
    def __init__(self, yaml_file="settings.yml"):
        super().__init__()
        self.setWindowTitle("YAML Editor + 3D Visualization with Clipping")
        self.showMaximized()
        
        self.file_path = yaml_file
        self.yaml_config = {}
        self.mesh_filename = None
        self.sol_file_path = None
        self.magnetization_data = None
        self.mesh_object = None
        self.clip_enabled = False
        self.clip_mesh = None
        
        self.setup_ui()
        self.load_yaml_file()

    def setup_ui(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        main_splitter = QSplitter(Qt.Horizontal)
        central_widget.setLayout(QHBoxLayout())
        central_widget.layout().addWidget(main_splitter)
        
        # === Left Panel ===
        left_panel = QFrame()
        left_panel.setFrameShape(QFrame.StyledPanel)
        left_panel.setMinimumWidth(500)
        left_layout = QVBoxLayout(left_panel)
        
        # Sol file selection
        sol_layout = QVBoxLayout()
        self.sol_label = QLabel(".sol File:")
        self.sol_input = QLineEdit()
        self.sol_btn = QPushButton("Browse .sol")
        self.sol_btn.clicked.connect(self.select_sol_file)
        sol_layout.addWidget(self.sol_label)
        sol_layout.addWidget(self.sol_input)
        sol_layout.addWidget(self.sol_btn)
        left_layout.addLayout(sol_layout)
        
        # Vector display controls
        vector_group = QGroupBox("Magnetization Vector Display")
        vector_layout = QVBoxLayout(vector_group)

        self.show_vectors_checkbox = QCheckBox("Show Magnetization Vectors")
        self.show_vectors_checkbox.setChecked(False)
        self.show_vectors_checkbox.stateChanged.connect(self.update_mesh)
        vector_layout.addWidget(self.show_vectors_checkbox)

        self.show_volume_vectors_checkbox = QCheckBox("Show Volume Vectors (not just surface)")
        self.show_volume_vectors_checkbox.setChecked(False)
        self.show_volume_vectors_checkbox.stateChanged.connect(self.update_mesh)
        vector_layout.addWidget(self.show_volume_vectors_checkbox)

        scale_layout = QHBoxLayout()
        scale_layout.addWidget(QLabel("Vector Scale:"))
        self.vector_scale_spinbox = QDoubleSpinBox()
        self.vector_scale_spinbox.setRange(0.001, 100.0)
        self.vector_scale_spinbox.setValue(1.0)
        self.vector_scale_spinbox.setSingleStep(0.1)
        self.vector_scale_spinbox.setDecimals(3)
        self.vector_scale_spinbox.setLocale(QLocale(QLocale.English, QLocale.UnitedStates))
        self.vector_scale_spinbox.valueChanged.connect(self.on_vector_scale_changed)
        scale_layout.addWidget(self.vector_scale_spinbox)
        vector_layout.addLayout(scale_layout)

        left_layout.addWidget(vector_group)
        
        # === Clipping Controls ===
        clip_group = QGroupBox("Clipping Plane")
        clip_layout = QVBoxLayout(clip_group)
        
        self.clip_enabled_checkbox = QCheckBox("Enable Clipping")
        self.clip_enabled_checkbox.setChecked(False)
        self.clip_enabled_checkbox.stateChanged.connect(self.on_clip_enabled_changed)
        clip_layout.addWidget(self.clip_enabled_checkbox)
        
        normal_layout = QHBoxLayout()
        normal_layout.addWidget(QLabel("Normal:"))
        self.clip_normal_combo = QComboBox()
        self.clip_normal_combo.addItems(["+X", "-X", "+Y", "-Y", "+Z", "-Z", "Custom"])
        self.clip_normal_combo.currentTextChanged.connect(self.on_clip_normal_changed)
        normal_layout.addWidget(self.clip_normal_combo)
        clip_layout.addLayout(normal_layout)
        
        custom_normal_layout = QHBoxLayout()
        custom_normal_layout.addWidget(QLabel("Custom Normal (x y z):"))
        self.clip_normal_inputs = [QLineEdit("0"), QLineEdit("0"), QLineEdit("1")]
        for inp in self.clip_normal_inputs:
            inp.textChanged.connect(self.on_custom_normal_changed)
            custom_normal_layout.addWidget(inp)
        clip_layout.addLayout(custom_normal_layout)
        
        for inp in self.clip_normal_inputs:
            inp.setEnabled(False)
        
        position_layout = QVBoxLayout()
        position_layout.addWidget(QLabel("Position:"))
        self.clip_position_spinbox = QDoubleSpinBox()
        self.clip_position_spinbox.setRange(-1000.0, 1000.0)
        self.clip_position_spinbox.setValue(0.0)
        self.clip_position_spinbox.setSingleStep(0.01)
        self.clip_position_spinbox.setDecimals(3)
        self.clip_position_spinbox.setLocale(QLocale(QLocale.English, QLocale.UnitedStates))
        self.clip_position_spinbox.valueChanged.connect(self.on_clip_position_changed)
        position_layout.addWidget(self.clip_position_spinbox)
        clip_layout.addLayout(position_layout)
        
        self.clip_invert_checkbox = QCheckBox("Invert Clipping")
        self.clip_invert_checkbox.setChecked(False)
        self.clip_invert_checkbox.stateChanged.connect(self.update_mesh)
        clip_layout.addWidget(self.clip_invert_checkbox)
        left_layout.addWidget(clip_group)
        
        # Rotation parameters
        rotation_group = QGroupBox("Rotation Parameters")
        rotation_layout = QVBoxLayout(rotation_group)
        
        angle1_layout = QVBoxLayout()
        angle1_label = QLabel("Angle 1 (°):")
        self.angle1_input = QLineEdit("0")
        self.angle1_input.textChanged.connect(self.update_mesh)
        angle1_layout.addWidget(angle1_label)
        angle1_layout.addWidget(self.angle1_input)
        rotation_layout.addLayout(angle1_layout)
        
        axe1_layout = QVBoxLayout()
        axe1_layout.addWidget(QLabel("Axis 1 (x y z):"))
        axis_input_layout = QHBoxLayout()
        self.axe1_inputs = [QLineEdit("1"), QLineEdit("0"), QLineEdit("0")]
        for inp in self.axe1_inputs:
            inp.textChanged.connect(self.update_mesh)
            axis_input_layout.addWidget(inp)
        axe1_layout.addLayout(axis_input_layout)
        rotation_layout.addLayout(axe1_layout)
        
        angle2_layout = QVBoxLayout()
        angle2_label = QLabel("Angle 2 (°):")
        self.angle2_input = QLineEdit("0")
        self.angle2_input.textChanged.connect(self.update_mesh)
        angle2_layout.addWidget(angle2_label)
        angle2_layout.addWidget(self.angle2_input)
        rotation_layout.addLayout(angle2_layout)
        
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
        
        det_layout = QVBoxLayout()
        zoom_label = QLabel("Zoom:")
        self.zoom_input = QLineEdit("0.5")
        meshsize_label = QLabel("Mesh Size:")
        self.meshsize_input = QLineEdit("1e-9")
        det_layout.addWidget(zoom_label)
        det_layout.addWidget(self.zoom_input)
        det_layout.addWidget(meshsize_label)
        det_layout.addWidget(self.meshsize_input)
        other_layout.addLayout(det_layout)
        left_layout.addWidget(other_group)
        
        yaml_label = QLabel("YAML Configuration:")
        left_layout.addWidget(yaml_label)
        
        self.yaml_display = QTextEdit()
        self.yaml_display.setReadOnly(True)
        self.highlighter = YAMLSyntaxHighlighter(self.yaml_display.document())
        left_layout.addWidget(self.yaml_display, 1)
        
        # === Right Panel ===
        right_panel = QFrame()
        right_panel.setFrameShape(QFrame.StyledPanel)
        right_layout = QVBoxLayout(right_panel)
        
        button_layout = QHBoxLayout()
        self.save_yaml_btn = QPushButton("Save YAML")
        self.save_yaml_btn.clicked.connect(self.save_yaml_config)
        self.quit_btn = QPushButton("Quit")
        self.quit_btn.clicked.connect(self.quit_and_save)
        button_layout.addWidget(self.save_yaml_btn)
        button_layout.addWidget(self.quit_btn)
        button_layout.addStretch()
        right_layout.addLayout(button_layout)
        
        title_label = QLabel("3D Visualization - Rotations, Magnetization & Clipping")
        title_label.setFont(QFont("Arial", 14, QFont.Bold))
        title_label.setAlignment(Qt.AlignCenter)
        right_layout.addWidget(title_label)
        
        self.plotter = QtInteractor()
        self.plotter.set_background("white")
        right_layout.addWidget(self.plotter.interactor)
        
        main_splitter.addWidget(left_panel)
        main_splitter.addWidget(right_panel)
        main_splitter.setSizes([500, 900])

    def select_sol_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select a .sol file", "", "SOL files (*.sol)")
        if path:
            self.sol_file_path = path
            self.sol_input.setText(path)
            try:
                self.magnetization_data = load_sol_file(path)
                QMessageBox.information(self, "Success", 
                    f"Loaded magnetization data for {len(self.magnetization_data)} nodes")
                self.update_mesh()
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Error loading .sol file: {e}")
                self.magnetization_data = None

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
                self.mesh_object, self.original_points, self.triangles = load_mesh(self.mesh_filename)
                self.update_mesh()
            except Exception as e:
                print(f"Error loading mesh: {e}")

    def populate_fields_from_yaml(self):
        if not self.yaml_config:
            return
        
        if 'initial_magnetization' in self.yaml_config:
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
        try:
            if not self.yaml_config:
                self.yaml_config = {}
            
            if not hasattr(self, 'sol_file_path') or not self.sol_file_path:
                raise ValueError("No .sol file selected")
            
            if 'rotations' not in self.yaml_config:
                self.yaml_config['rotations'] = {}
            
            self.yaml_config['rotations']['angle1'] = float(self.angle1_input.text())
            self.yaml_config['rotations']['axe1'] = [float(inp.text()) for inp in self.axe1_inputs]
            self.yaml_config['rotations']['angle2'] = float(self.angle2_input.text())
            self.yaml_config['rotations']['axe2'] = [float(inp.text()) for inp in self.axe2_inputs]
            
            self.yaml_config['filled'] = str(self.filled_checkbox.isChecked()).lower()
            
            if 'electrostatics' not in self.yaml_config:
                self.yaml_config['electrostatics'] = {}
            self.yaml_config['electrostatics']['CE'] = int(self.ce_input.text())
            self.yaml_config['electrostatics']['V'] = float(self.v_input.text())
            
            if 'detector' not in self.yaml_config:
                self.yaml_config['detector'] = {}
            zoom = float(self.zoom_input.text())
            if zoom > 1.0:
                zoom = 1.0
                self.zoom_input.setText(str(zoom))
            self.yaml_config['detector']['zoom'] = zoom
            self.yaml_config['detector']['meshSize'] = float(self.meshsize_input.text())
            
            self.update_yaml_display()
        
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error updating: {e}")
            return
        
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
            angle1 = float(self.angle1_input.text())
            axis1 = [float(inp.text()) for inp in self.axe1_inputs]
            axes_angles.append((axis1, angle1))
            
            angle2 = float(self.angle2_input.text())
            axis2 = [float(inp.text()) for inp in self.axe2_inputs]
            axes_angles.append((axis2, angle2))
            
            return axes_angles
        except ValueError:
            return [([1, 0, 0], 0), ([0, 1, 0], 0)]

    def on_vector_scale_changed(self, value):
        if self.show_vectors_checkbox.isChecked():
            self.update_mesh()

    def on_clip_enabled_changed(self, state):
        self.clip_enabled = (state == Qt.Checked)
        self.update_mesh()

    def on_clip_normal_changed(self, text):
        is_custom = (text == "Custom")
        for inp in self.clip_normal_inputs:
            inp.setEnabled(is_custom)
        if not is_custom:
            self.update_mesh()

    def on_custom_normal_changed(self):
        if self.clip_normal_combo.currentText() == "Custom":
            self.update_mesh()

    def on_clip_position_changed(self, value):
        if self.clip_enabled:
            self.update_mesh()

    def get_clip_normal(self):
        preset = self.clip_normal_combo.currentText()
        normals = {
            "+X": [1, 0, 0],
            "-X": [-1, 0, 0],
            "+Y": [0, 1, 0],
            "-Y": [0, -1, 0],
            "+Z": [0, 0, 1],
            "-Z": [0, 0, -1]
        }
        if preset == "Custom":
            try:
                return normalize([float(inp.text()) for inp in self.clip_normal_inputs])
            except ValueError:
                return np.array([0, 0, 1])
        else:
            return np.array(normals.get(preset, [0, 0, 1]))

    def update_mesh(self):
        if not hasattr(self, 'original_points') or not hasattr(self, 'triangles'):
            return
        
        try:
            axes_angles = self.get_rotations_from_inputs()
            rotated_points = apply_sequential_quaternion_rotations(self.original_points, axes_angles)
            mesh = create_pyvista_mesh(rotated_points, self.triangles)
            
            self.plotter.clear()
            
            # Apply clipping if enabled
            if self.clip_enabled:
                clip_normal = self.get_clip_normal()
                clip_position = self.clip_position_spinbox.value()
                
                origin = mesh.center + clip_normal * clip_position
                
                clipped = mesh.clip(
                    normal=clip_normal,
                    origin=origin,
                    invert=self.clip_invert_checkbox.isChecked()
                )
                
                self.plotter.add_mesh(clipped, color="lightblue", show_edges=False)
                self.clip_mesh = clipped
                
                # Show the clipping plane
                plane_size = np.linalg.norm(mesh.bounds) * 0.5
                plane = pv.Plane(
                    center=origin,
                    direction=clip_normal,
                    i_size=plane_size,
                    j_size=plane_size
                )
                self.plotter.add_mesh(plane, color="yellow", opacity=0.3, show_edges=True)
            else:
                self.plotter.add_mesh(mesh, color="lightblue", show_edges=False)
                self.clip_mesh = None
            
            self.add_bounding_box(self.original_points, axes_angles)
            
            # Add magnetization vectors if available and enabled
            if (self.show_vectors_checkbox.isChecked() and 
                self.magnetization_data is not None and 
                self.mesh_object is not None):
                self.add_magnetization_vectors(axes_angles)
            
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

    def add_magnetization_vectors(self, axes_angles):
        """Add magnetization vectors to surface or volume nodes using optimized glyphs."""
        try:
            # Determine which nodes to show
            if self.show_volume_vectors_checkbox.isChecked():
                # Show all nodes that have magnetization data
                nodes_to_show = set(self.magnetization_data.keys())
                print(f"Showing volume vectors for all {len(nodes_to_show)} nodes with magnetization data")
            else:
                # Show only surface nodes
                nodes_to_show = get_surface_nodes(self.mesh_object)
                print(f"Showing surface vectors for {len(nodes_to_show)} surface nodes")
            
            positions = []
            vectors = []
            
            vector_scale = self.vector_scale_spinbox.value()
            print(f"Vector scale: {vector_scale}")
            
            # Collect all positions and vectors
            for node_idx in nodes_to_show:
                if node_idx in self.magnetization_data and node_idx < len(self.original_points):
                    node_pos = self.original_points[node_idx]
                    m_vec = self.magnetization_data[node_idx]
                    
                    positions.append(node_pos)
                    vectors.append(m_vec)
            
            if not positions:
                print("No magnetization data found for selected nodes!")
                return
            
            print(f"Found magnetization data for {len(positions)} nodes")
            
            # Convert to numpy arrays for batch rotation
            positions = np.array(positions)
            vectors = np.array(vectors)
            
            print(f"Position range: {positions.min(axis=0)} to {positions.max(axis=0)}")
            print(f"Vector magnitude range: {np.linalg.norm(vectors, axis=1).min():.3f} to {np.linalg.norm(vectors, axis=1).max():.3f}")
            
            # Apply rotations in batch
            rotated_positions = apply_sequential_quaternion_rotations(positions, axes_angles)
            rotated_vectors = apply_sequential_quaternion_rotations(vectors, axes_angles)
            
            # If clipping is enabled, filter vectors to only show those in the clipped region
            if self.clip_enabled and self.clip_mesh is not None:
                clip_normal = self.get_clip_normal()
                clip_position = self.clip_position_spinbox.value()
                
                # Get the original mesh center for calculating origin
                temp_mesh = create_pyvista_mesh(
                    apply_sequential_quaternion_rotations(self.original_points, axes_angles), 
                    self.triangles
                )
                origin = temp_mesh.center + clip_normal * clip_position
                
                # Calculate which points are on the visible side of the clip plane
                vectors_to_plane = rotated_positions - origin
                distances = np.dot(vectors_to_plane, clip_normal)
                
                # Keep points based on invert setting
                if self.clip_invert_checkbox.isChecked():
                    mask = distances < 0  # Keep points behind the plane
                else:
                    mask = distances >= 0  # Keep points in front of the plane
                
                # Filter positions and vectors
                rotated_positions = rotated_positions[mask]
                rotated_vectors = rotated_vectors[mask]
                
                print(f"After clipping: {len(rotated_positions)} vectors remaining")
                
                if len(rotated_positions) == 0:
                    print("No vectors in clipped region!")
                    return
            
            # Scale vectors for visualization
            scaled_vectors = rotated_vectors * vector_scale
            
            # Calculate colors based on z-component
            z_components = rotated_vectors[:, 2]
            
            up_count = np.sum(z_components > 0.5)
            down_count = np.sum(z_components < -0.5)
            inplane_count = len(z_components) - up_count - down_count
            
            print(f"Magnetization distribution - Up: {up_count}, Down: {down_count}, In-plane: {inplane_count}")
            
            # Create a PolyData object with points and vectors
            vector_poly = pv.PolyData(rotated_positions)
            vector_poly['vectors'] = scaled_vectors
            
            print(f"Creating glyphs with scale factor {vector_scale}")
            
            # Create arrows using glyphs
            arrows = vector_poly.glyph(
                orient='vectors',
                scale=False,
                factor=vector_scale,
                geom=pv.Arrow(tip_length=0.25, tip_radius=0.1, shaft_radius=0.05)
            )
            
            print(f"Created {arrows.n_cells} arrow cells, {arrows.n_points} points")
            
            # Add to plotter in black
            self.plotter.add_mesh(
                arrows, 
                color='black',
                opacity=0.8,
                show_scalar_bar=False
            )
            
            print("Magnetization vectors added successfully")
                    
        except Exception as e:
            print(f"Error adding magnetization vectors: {e}")
            import traceback
            traceback.print_exc()
            
    def add_bounding_box(self, original_points, axes_angles):
        """Add a bounding box around the mesh using 8 corner nodes and 12 edges."""
        xmin, ymin, zmin = np.min(original_points, axis=0)
        xmax, ymax, zmax = np.max(original_points, axis=0)
        
        corners = np.array([
            [xmin, ymin, zmin],  # 0
            [xmax, ymin, zmin],  # 1
            [xmax, ymax, zmin],  # 2
            [xmin, ymax, zmin],  # 3
            [xmin, ymin, zmax],  # 4
            [xmax, ymin, zmax],  # 5
            [xmax, ymax, zmax],  # 6
            [xmin, ymax, zmax],  # 7
        ])
        
        rotated_corners = apply_sequential_quaternion_rotations(corners, axes_angles)
        
        edges_x = [[0, 1], [3, 2], [4, 5]]
        edges_y = [[1, 2], [0, 3], [4, 7]]
        edges_z = [[0, 4], [1, 5], [3, 7]]
        edges   = [[7, 6], [5, 6], [2, 6]]
        
        for edge in edges_x:
            line_points = np.array([rotated_corners[edge[0]], rotated_corners[edge[1]]])
            line = pv.Line(line_points[0], line_points[1])
            self.plotter.add_mesh(line, color="red", line_width=4)
        
        for edge in edges_y:
            line_points = np.array([rotated_corners[edge[0]], rotated_corners[edge[1]]])
            line = pv.Line(line_points[0], line_points[1])
            self.plotter.add_mesh(line, color="green", line_width=4)
        
        for edge in edges_z:
            line_points = np.array([rotated_corners[edge[0]], rotated_corners[edge[1]]])
            line = pv.Line(line_points[0], line_points[1])
            self.plotter.add_mesh(line, color="blue", line_width=4)
        
        for edge in edges:
            line_points = np.array([rotated_corners[edge[0]], rotated_corners[edge[1]]])
            line = pv.Line(line_points[0], line_points[1])
            self.plotter.add_mesh(line, color="black", line_width=2)

    def quit_and_save(self):
        """Save YAML configuration before quitting"""
        try:
            self.save_yaml_config()
        except Exception as e:
            reply = QMessageBox.question(
                self, 
                'Error Saving', 
                f"Error saving configuration: {e}\n\nDo you still want to quit?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No
            )
            if reply == QMessageBox.No:
                return
        
        self.close()

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
