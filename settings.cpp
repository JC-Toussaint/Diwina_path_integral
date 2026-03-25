#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>  // for sysconf(), gethostname()

#include "termcolor.h"
#include "tags.h"
#include "settings.h"

using namespace Nodes;

/***********************************************************************
 * Access to the default configuration embedded from the file
 * default-settings.yml.
 */

extern "C"
    {
    extern char _binary_default_settings_yml_start[];
    extern char _binary_default_settings_yml_end[];
    }

/***********************************************************************
 * Private helper functions.
 */

// Return the YAML document defining the defaults.
static std::string get_default_yaml()
    {
    return std::string(_binary_default_settings_yml_start,
                       _binary_default_settings_yml_end - _binary_default_settings_yml_start);
    }

// Bail out on errors.
static void error(const char *message)
    {
    std::cerr << "CONFIGURATION ERROR: " << message << "\n";
    exit(1);
    }

// Conditionally assign a variable if the node is defined and scalar.
template<typename T>
static bool assign(T &var, const YAML::Node &node)
    {
    if (node.IsScalar())
        {
        var = node.as<T>();
        return true;
        }
    return false;
    }

// Stringify a string: if it contains a newline, convert it to a multiline string in "literal style"
// (introduced by '|'). Otherwise leave it alone.
// Warning: this function assumes this is the value of a property at indentation level zero.
static const std::string str(std::string s)
    {
    // If it doesn't ave an eol, leave it alone.
    if (s.find('\n') == std::string::npos) return s;

    // Prepend "|\n".
    s.insert(0, "|\n");

    // Remove trailing eol.
    if (s.back() == '\n') s.resize(s.size() - 1);

    // Add indentation: replace "\n" with "\n  ".
    size_t pos = 0;
    while ((pos = s.find('\n', pos)) != std::string::npos)
        {
        s.replace(pos + 1, 0, "  ");
        pos += 3;
        }

    return s;
    }

/***********************************************************************
 * Public API.
 */

Settings::Settings()
    {
    verbose = 0;
    withTsv = true;
    read(YAML::Load(get_default_yaml()));  // load defaults
    }

void Settings::dumpDefaults() { std::cout << get_default_yaml(); }

void Settings::infos()
    {
    std::cout << "outputs:\n";
    std::cout << "  directory: " << r_path_output_dir << "\n";
    std::cout << "  file_basename: " << simName << "\n";
    std::cout << "mesh:\n";
    std::cout << "  filename: " << pbName << "\n";
    std::cout << "  length_unit: " << _scale << "\n";
    std::cout << "  volume_regions:\n";
    for (auto it = paramTetra.begin(); it != paramTetra.end(); ++it)
        {
        if (it->regName == "__default__")  // skip
            continue;
        std::cout << "    " << it->regName << ":\n";
        std::cout << "      Js: " << it->J << "\n";
        std::cout << "      lbp: " << it->lbp << "\n";
        std::cout << "      lbm: " << it->lbm << "\n";
        }
    std::cout << "  surface_regions:\n";
    for (auto it = paramFacette.begin(); it != paramFacette.end(); ++it)
        {
        if (it->regName == "__default__")  // skip
            continue;
        std::cout << "    " << it->regName << ":\n";
        }
    std::cout << "initial_magnetization: ";
    if (!sM.empty())
        std::cout << str(sM) << "\n";
    else if (restoreFileName.empty())
        std::cout << "[\"" << sMx << "\", \"" << sMy << "\", \"" << sMz << "\"]\n";
    else
        std::cout << restoreFileName << "\n";
    std::cout << "rotations:\n";
    std::cout << "  angle1: " << p_rot.angle1 << "\n";
    std::cout << "  axe1: [" << p_rot.axe1[0] << ", " << p_rot.axe1[1] << ", "
              << p_rot.axe1[2] << "]\n";
    std::cout << "  angle2: " << p_rot.angle2 << "\n";
    std::cout << "  axe2: [" << p_rot.axe2[0] << ", " << p_rot.axe2[1] << ", "
              << p_rot.axe2[2] << "]\n";
    std::cout << "filled: " << (filled ? "true" : "false") << "\n";
    std::cout << "electrostatics:\n";
    std::cout << "  CE: " << p_electrostatics.CE << "\n";
    std::cout << "  V: "  << p_electrostatics.V  << "\n";
    std::cout << "detector:\n";
    std::cout << "  zoom: " << p_detector.zoomFactor << "\n";
    std::cout << "  meshSize: " << p_detector.meshSize << "\n";
    }

void Settings::read(YAML::Node yaml)
    {
    YAML::Node outputs = yaml["outputs"];
    if (outputs)
        {
        if (!outputs.IsMap()) error("outputs should be a map.");
        if (assign(r_path_output_dir, outputs["directory"]))
            {
            // Normalize directory name.
            if (r_path_output_dir.empty())
                {
                r_path_output_dir = ".";
                }
            if (r_path_output_dir.length() > 1 && r_path_output_dir.back() == '/')
                {
                r_path_output_dir.pop_back();
                }
            }
        assign(simName, outputs["file_basename"]);
        }  // outputs

    YAML::Node mesh = yaml["mesh"];
    if (mesh)
        {
        if (!mesh.IsMap()) error("mesh should be a map.");
        assign(pbName, mesh["filename"]);
        if (assign(_scale, mesh["length_unit"]) && _scale <= 0)
            error("mesh.length_unit should be positive.");
        YAML::Node volumes = mesh["volume_regions"];
        if (volumes)
            {
            if (!volumes.IsMap()) error("mesh.volume_regions should be a map.");
            int default_idx = findTetraRegionIdx("__default__");
            for (auto it = volumes.begin(); it != volumes.end(); ++it)
                {
                std::string name = it->first.as<std::string>();
                YAML::Node volume = it->second;
                Tetra::prm p;
                if (default_idx >= 0) p = paramTetra[default_idx];
                p.regName = name;
                assign(p.J, volume["Js"]);
                assign(p.lbp, volume["lbp"]);
                assign(p.lbm, volume["lbm"]);
                if ((p.lbp==0) || (p.lbm==0))
                    {
				    std::cout << termcolor::bright_red << termcolor::blink << "Absorption coefficients lbp or lbm should be positive" << termcolor::reset << std::endl;
				    exit(1);
				    }
                paramTetra.push_back(p);
                }
            }  // mesh.volume_regions
        YAML::Node surfaces = mesh["surface_regions"];
        if (surfaces)
            {
            if (!surfaces.IsMap()) error("mesh.surface_regions should be a map.");
            int default_idx = findFacetteRegionIdx("__default__");
            for (auto it = surfaces.begin(); it != surfaces.end(); ++it)
                {
                std::string name = it->first.as<std::string>();
                YAML::Node surface = it->second;
                Facette::prm p;
                if (default_idx >= 0) p = paramFacette[default_idx];
                p.regName = name;
                paramFacette.push_back(p);
                }
            }  // mesh.surface_regions
        }      // mesh

    YAML::Node magnetization = yaml["initial_magnetization"];
    if (magnetization)
        {
        if (magnetization.IsScalar())
            {
            std::string s_mag = magnetization.as<std::string>();
            if (s_mag.find("function") == std::string::npos || s_mag.find('{') == std::string::npos)
                {
                restoreFileName = s_mag;
                }
            else
                {
                sM = s_mag;
                mag_parser.set_function(sM);
                }
            }
        else if (magnetization.IsSequence())
            {
            if (magnetization.size() != DIM)
                error("initial_magnetization should have three components.");
            sMx = magnetization[0].as<std::string>();
            sMy = magnetization[1].as<std::string>();
            sMz = magnetization[2].as<std::string>();
            mag_parser.set_expressions("x,y,z", sMx, sMy, sMz);
            }
        else
            {
            error("initial_magnetization should be a file name or a vector of expressions.");
            }
        }  // initial_magnetization

    YAML::Node rot = yaml["rotations"];
    if (rot)
       {
       if (rot["angle1"] && rot["axe1"]) {
          assign(p_rot.angle1, rot["angle1"]);
          if (rot["axe1"].size() != DIM) error("axe1 should have three components.");
          auto tmp1 = rot["axe1"].as<std::vector<double>>();
          p_rot.axe1 = Eigen::Vector3d(tmp1[0], tmp1[1], tmp1[2]);
	  p_rot.axe1.normalize();
          }
       if (rot["angle2"] && rot["axe2"]) {
          assign(p_rot.angle2, rot["angle2"]);
          if (rot["axe2"].size() != DIM) error("axe2 should have three components.");
          auto tmp2 = rot["axe2"].as<std::vector<double>>();
          p_rot.axe2 = Eigen::Vector3d(tmp2[0], tmp2[1], tmp2[2]);
	  p_rot.axe2.normalize();
          }
       }

    assign(filled, yaml["filled"]);

    YAML::Node electrostatics = yaml["electrostatics"];
    if (electrostatics)
       {
       assign(p_electrostatics.CE, electrostatics["CE"]);
       assign(p_electrostatics.V,  electrostatics["V"]);
       }
 
    YAML::Node detector = yaml["detector"];
    if (detector)
       {
       assign(p_detector.zoomFactor, detector["zoom"]);
       assign(p_detector.meshSize, detector["meshSize"]);
       }

    // outputs.file_basename defaults to base name of mesh.filename.
    if (simName.empty() && !pbName.empty())
        {
        simName = pbName;

        // Remove directory part.
        size_t pos = simName.rfind('/');
        if (pos != simName.npos) simName.erase(0, pos + 1);

        // Remove everything after last dot.
        pos = simName.rfind('.');
        if (pos != simName.npos) simName.erase(pos);
        }
    }

bool Settings::read(std::string filename)
    {
    std::ostringstream buffer;
    if (filename == "-")
        {
        buffer << std::cin.rdbuf();
        }
    else
        {
        std::ifstream input(filename, std::ios::in);
        buffer << input.rdbuf();
        input.close();
        }
    yaml_source = buffer.str();

    YAML::Node config = YAML::Load(yaml_source);
    if (config.IsNull()) return false;
    read(config);
    return true;
    }
