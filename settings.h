#ifndef feellgoodSettings_h
#define feellgoodSettings_h

/** \file feellgoodSettings.h
\brief many settings to give some parameters to the solver, boundary conditions for the problem, the
output file format wanted by the user. This is done mainly with the class Settings.
*/

#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>

#include "expression_parser.h"
#include "facette.h"
#include "rotation.h"
#include "electrostatics.h"
#include "detector.h"
#include "tetra.h"

/** \class Settings
container class to store many setting parameters, such as file names, parameters for the solver,
output file format. It also handles text user interation through terminal, and some parsing
functions.
*/
class Settings
    {
public:
    /** default constructor */
    Settings();

    /** print out the YAML document defining the default settings */
    static void dumpDefaults();

    /** some prints sent to terminal */
    void infos(void);

    /** read settings from a parsed YAML document */
    void read(YAML::Node);

    /** read settings from a YAML file. Returns true if a non-empty configuration is found. */
    bool read(std::string filename);

    /** getter for fileDisplayName */
    inline std::string getFileDisplayName(void) const { return fileDisplayName; }

    /** setter for fileDisplayName */
    inline void setFileDisplayName(std::string _s) { fileDisplayName = _s; }

    /** setter for .msh file name */
    inline void setPbName(std::string str) { pbName = str; }

    /** getter for problem file name */
    inline std::string getPbName(void) const { return pbName; }

    /** setter for .sol output file name */
    inline void setSimName(std::string str) { simName = str; }

    /** getter for output file name */
    inline std::string getSimName(void) const { return simName; }

    /** setter for geometrical scaling factor for physical coordinates of the mesh */
    inline void setScale(const double s) { _scale = s; }

    /** getter for geometrical scaling factor for physical coordinates of the mesh */
    inline double getScale(void) const { return _scale; }

    /** YAML source of these settings. **/
    std::string yaml_source;

    /** boolean flag to mention if you want output in txt tsv file format */
    bool withTsv;

    /** verbosity level, defaults to zero */
    int verbose;

    /** string for analytical definition of Mx */
    std::string sMx;

    /** string for analytical definition of My */
    std::string sMy;

    /** string for analytical definition of Mz */
    std::string sMz;

    /** string for a JavaScript function defining M, alternative to (sMx, sMy, sMz) */
    std::string sM;

    /** rotation parameters */
    ROT p_rot;
    
    /** whether a hollow system is filled by Cu or empty */
    bool filled;

    /** electrostatics parameters */
    ELECTROSTATICS p_electrostatics;

    /** detector parameters */
    DETECTOR p_detector;
    
    /** input file name for continuing a calculation (sol.in) */
    std::string restoreFileName;

    /** this vector contains the material parameters for all regions for all the tetrahedrons */
    std::vector<Tetra::prm> paramTetra;

    /** \return index of the region in volume region container  */
    inline int findTetraRegionIdx(const std::string name /**< [in] */) const
        {
        std::vector<Tetra::prm>::const_iterator result =
                std::find_if(paramTetra.begin(), paramTetra.end(),
                             [name](Tetra::prm const &p) { return (p.regName == name); });

        int idx(-2);
        if (result == paramTetra.end())
            {
            idx = -1;
            }
        else
            {
            idx = std::distance(paramTetra.begin(), result);
            }
        return idx;
        };

    /** this vector contains the material parameters for all regions for all the facettes */
    std::vector<Facette::prm> paramFacette;

    /** relative path for output files (to be implemented) */
    std::string r_path_output_dir;

    /** \return index of the region in surface region container  */
    inline int findFacetteRegionIdx(const std::string name /**< [in] */) const
        {
        std::vector<Facette::prm>::const_iterator result =
                std::find_if(paramFacette.begin(), paramFacette.end(),
                             [name](Facette::prm const &p) { return (p.regName == name); });
        int idx(-2);

        if (result == paramFacette.end())
            {
            idx = -1;
            }
        else
            {
            idx = std::distance(paramFacette.begin(), result);
            }
        return idx;
        };

    /** evaluation of the magnetization components through math expression, each component of the
     * magnetization is a function of (x,y,z). It is not safe to call this method simultaneously
     * from multiple threads.
     * \return unit vector
     */
    inline Eigen::Vector3d getMagnetization(const Eigen::Ref<Eigen::Vector3d> p) const
        {
        Eigen::Vector3d tmp = mag_parser.get_vector(p);
        tmp.normalize();
        return tmp;
        }

private:
    std::string fileDisplayName; /**< parameters file name : either a yaml file or standard input */
    double _scale;               /**< scaling factor from gmsh files to feellgood */
    std::string simName;         /**< simulation name */
    std::string pbName;          /**< mesh file, gmsh file format */
    VectorParser mag_parser;     /**< parser for the magnetization expressions */
    };

#endif /* feellgoodSettings_h */
