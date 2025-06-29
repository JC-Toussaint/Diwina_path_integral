#!/usr/bin/env python3
"""
FeeLLGood configuration file analyzer
Displays the complete YAML tree structure and analyzes the initial_magnetization field
"""

import yaml
import json
import sys
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

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python fg-yml <setting.yml>")
        sys.exit(1)

    try:
        filename = sys.argv[1]
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)
        
    analyze()
