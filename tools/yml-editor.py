#!/usr/bin/env python3

import os
import yaml
import json
import sys

from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QLabel, QPushButton,
    QFileDialog, QHBoxLayout, QLineEdit, QMessageBox, QCheckBox
)
from PyQt5.QtGui import QFont

from typing import Any, Dict, List, Union

def print_tree(data: Any, prefix: str = "", is_last: bool = True, level: int = 0) -> None:
    """
    Displays the tree structure of a data structure in a hierarchical manner
    """
    if level == 0:
        print("📁 FeeLLGood Configuration")
        print("│")
    
    if isinstance(data, dict):
        items = list(data.items())
        for i, (key, value) in enumerate(items):
            is_last_item = (i == len(items) - 1)
            
            # Tree symbols
            if is_last_item:
                current_prefix = prefix + "└── "
                next_prefix = prefix + "    "
            else:
                current_prefix = prefix + "├── "
                next_prefix = prefix + "│   "
            
            # Display key with icon according to type
            if isinstance(value, dict):
                print(f"{current_prefix}📂 {key}/")
                print_tree(value, next_prefix, is_last_item, level + 1)
            elif isinstance(value, list):
                print(f"{current_prefix}📋 {key}[] ({len(value)} elements)")
                print_tree(value, next_prefix, is_last_item, level + 1)
            else:
                type_icon = "🔢" if isinstance(value, (int, float)) else "📝"
                print(f"{current_prefix}{type_icon} {key}: {repr(value)}")
    
    elif isinstance(data, list):
        for i, item in enumerate(data):
            is_last_item = (i == len(data) - 1)
            
            if is_last_item:
                current_prefix = prefix + "└── "
                next_prefix = prefix + "    "
            else:
                current_prefix = prefix + "├── "
                next_prefix = prefix + "│   "
            
            if isinstance(item, (dict, list)):
                print(f"{current_prefix}[{i}]")
                print_tree(item, next_prefix, is_last_item, level + 1)
            else:
                type_icon = "🔢" if isinstance(item, (int, float)) else "📝"
                print(f"{current_prefix}{type_icon} [{i}]: {repr(item)}")

def analyze_initial_magnetization(config: Dict[str, Any]) -> None:
    """
    Specific analysis of the initial_magnetization field
    """ 
    print("\n" + "="*80)
    print("⚠️ INITIAL_MAGNETIZATION FIELD")
    print("="*80)
    
    if 'initial_magnetization' not in config:
        print("❌ The 'initial_magnetization' field is not present in the configuration")
        return
    
    initial_mag = config['initial_magnetization']
    print(initial_mag)
				        

def clean_yaml_content(content: str) -> str:
    """
    Cleans YAML content by replacing tabs with spaces
    """
    # Replace tabs with spaces (2 spaces per tab)
    lines = content.split('\n')
    cleaned_lines = []
    
    for line in lines:
        # Replace tabs with 2 spaces
        cleaned_line = line.replace('\t', '  ')
        cleaned_lines.append(cleaned_line)
    
    return '\n'.join(cleaned_lines)

def load_yaml_from_file(filename: str = "paste.txt") -> Dict[str, Any]:
    """
    Loads the YAML file from the uploaded file
    """
    try:
        # Read file content
        with open(filename, 'r', encoding='utf-8') as file:
            content = file.read()
        
        # Clean content (replace tabs)
        cleaned_content = clean_yaml_content(content)
        
        # Parse YAML
        config = yaml.safe_load(cleaned_content)
        return config
        
    except FileNotFoundError:
        print(f"❌ File '{filename}' not found")
        return None
    except Exception as e:
        print(f"❌ Error reading file: {e}")
        return None

def analyze():
    """
    analyze function
    """            
    try:
        # Loading YAML file
        print(f"📂 Loading file {filename}...")
        config = load_yaml_from_file(filename)
        
        if config is None:
            print("❌ Unable to load configuration file")
            return
        
        print("✅ File loaded successfully!")
        print("\n🔍 FEELLGOOD CONFIGURATION FILE ANALYSIS")
        print("="*80)
        print(f"📊 Number of main sections: {len(config)}")
        print(f"🏷️  Sections: {', '.join(config.keys())}")
        print("\n")
        
        # Display complete tree structure
        print_tree(config)
        
        # Specific analysis of initial_magnetization field
        analyze_initial_magnetization(config)
        
        print("\n" + "="*80)
        print("✅ ANALYSIS COMPLETED")
        print("="*80)
        
    except yaml.YAMLError as e:
        print(f"❌ Error analyzing YAML: {e}")
        print("💡 Check YAML file indentation (use spaces, not tabs)")
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        
class SolSelector(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Sélection .sol et paramètres")
        self.file_path = None

        # Définir la police par défaut pour toute la fenêtre
        font = QFont("Arial", 12)  # Police Arial, taille 12
        self.setFont(font)

        layout = QVBoxLayout()

        # Fichier .sol
        self.file_label = QLabel("Aucun fichier sélectionné")
        self.select_file_btn = QPushButton("Sélectionner un fichier .sol")
        self.select_file_btn.clicked.connect(self.select_file)
        layout.addWidget(self.select_file_btn)
        layout.addWidget(self.file_label)

        # Paramètres existants
        layout.addWidget(QLabel("Angle 1 (degré)"))
        self.angle1_input = QLineEdit("90")
        layout.addWidget(self.angle1_input)

        layout.addWidget(QLabel("Axe 1 (x y z)"))
        axe1_layout = QHBoxLayout()
        self.axe1_inputs = [QLineEdit(str(val)) for val in [1, 0, 0]]
        for inp in self.axe1_inputs:
            axe1_layout.addWidget(inp)
        layout.addLayout(axe1_layout)

        layout.addWidget(QLabel("Angle 2 (degré)"))
        self.angle2_input = QLineEdit("0")
        layout.addWidget(self.angle2_input)

        layout.addWidget(QLabel("Axe 2 (x y z)"))
        axe2_layout = QHBoxLayout()
        self.axe2_inputs = [QLineEdit(str(val)) for val in [0, 1, 0]]
        for inp in self.axe2_inputs:
            axe2_layout.addWidget(inp)
        layout.addLayout(axe2_layout)

        # Nouveaux paramètres
        self.filled_checkbox = QCheckBox("filled")
        self.filled_checkbox.setChecked(True)
        layout.addWidget(self.filled_checkbox)

        layout.addWidget(QLabel("electrostatics.CE"))
        self.ce_input = QLineEdit("0")
        layout.addWidget(self.ce_input)

        layout.addWidget(QLabel("electrostatics.V"))
        self.v_input = QLineEdit("0.0")
        layout.addWidget(self.v_input)

        layout.addWidget(QLabel("detector.zoom"))
        self.zoom_input = QLineEdit("1.0")
        layout.addWidget(self.zoom_input)

        layout.addWidget(QLabel("detector.meshSize"))
        self.meshsize_input = QLineEdit("1e-9")
        layout.addWidget(self.meshsize_input)

        # Validation
        self.ok_btn = QPushButton("Valider")
        self.ok_btn.clicked.connect(self.validate)
        layout.addWidget(self.ok_btn)

        self.setLayout(layout)
        self.parameters = None

    def select_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Sélectionner un fichier .sol", "", "Fichiers .sol (*.sol)")
        if path:
            self.file_path = path
            self.file_label.setText(path)

    def validate(self):
        try:
            if not self.file_path:
                raise ValueError("Aucun fichier .sol sélectionné")

            angle1 = float(self.angle1_input.text())
            axe1 = [float(inp.text()) for inp in self.axe1_inputs]
            angle2 = float(self.angle2_input.text())
            axe2 = [float(inp.text()) for inp in self.axe2_inputs]

            filled = self.filled_checkbox.isChecked()
            ce = int(self.ce_input.text())
            v = float(self.v_input.text())
            zoom = float(self.zoom_input.text())
            mesh_size = float(self.meshsize_input.text())

            self.parameters = (
                self.file_path, angle1, axe1, angle2, axe2,
                filled, ce, v, zoom, mesh_size
            )
            self.close()
        except Exception as e:
            QMessageBox.critical(self, "Erreur", f"Erreur de saisie : {e}")

def get_sol_and_angles():
    app = QApplication(sys.argv)
    
    # MÉTHODE ALTERNATIVE : Définir la police pour toute l'application
    font = QFont("Arial", 12)
    app.setFont(font)
    
    window = SolSelector()
    window.show()
    app.exec_()
    return window.parameters

def update_settings(settings, parameters):
    settings['initial_magnetization'] = os.path.basename(parameters[0])

    settings['rotations']['angle1'] = parameters[1]
    settings['rotations']['axe1']   = parameters[2]
    settings['rotations']['angle2'] = parameters[3]
    settings['rotations']['axe2']   = parameters[4]

    settings['filled'] = str(parameters[5]).lower()
             
    settings['electrostatics']['CE'] =  parameters[6]
    settings['electrostatics']['V']  =  parameters[7]
          
    settings['detector']['zoom']     =  parameters[8]
    settings['detector']['meshSize'] =  parameters[9]
    
def save_yaml_file(data: dict, filename: str):
    """
    Sauvegarde une hiérarchie YAML (dictionnaire Python) dans un fichier.

    :param data: Dictionnaire représentant le contenu YAML
    :param filename: Chemin du fichier de sortie (.yml ou .yaml)
    """
    with open(filename, 'w') as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)
        
# --- Programme principal ---
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python dialog.py <setting.yml>")
        sys.exit(1)

    try:
        parameters = get_sol_and_angles()
        if parameters:
            print("Fichier :", parameters[0])
            print("Angle 1 :", parameters[1])
            print("Axe   1 :", parameters[2])
            print("Angle 2 :", parameters[3])
            print("Axe   2 :", parameters[4])
            print("Filled :", parameters[5])
            print("CE     :", parameters[6])
            print("V      :", parameters[7])
            print("Zoom   :", parameters[8])
            print("MeshSz :", parameters[9])

            filename = sys.argv[1]
            settings = load_yaml_from_file(filename)
        # Display complete tree structure
            print_tree(settings)
            print(settings)

            update_settings(settings, parameters)
        # Display complete tree structure
            print_tree(settings)
                
            base, ext = os.path.splitext(filename)
            filename = f"{base}_ray{ext}"  
            save_yaml_file(settings, filename)
                 
            print("Fichier mis à jour avec succès.")
        else:
            print("Opération annulée.")
    except Exception as e:
        print(f"Erreur : {e}")
        sys.exit(1)
