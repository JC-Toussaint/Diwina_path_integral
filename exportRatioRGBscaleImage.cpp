#include <cstdint>
#include <limits>
#include <vector>
#include <filesystem>
#include "termcolor.h"

#include <cmath>
#include <string>
#include <png.h>
#include <tinyxml2.h>
#include <cstring>

#include <stdexcept>
#include <sstream>
#include <cstring>

#include <unordered_map>
#include <functional>

#include "config.h"
#include "fem2d.h"
#include "write_png.h"

static Colormap jet_colors = {
    {0,    {0, 0, 1}},  // blue
    {0.25, {0, 1, 1}},  // cyan
    {0.5,  {0, 1, 0}},  // green
    {0.75, {1, 1, 0}},  // yellow
    {1,    {1, 0, 0}}   // red
};

/**
static Colormap read_colormap(const std::string& filename)
Description

This function loads a colormap from an XML file and returns it as a Colormap. A colormap is a sequence of points, where:

    Each point consists of a position (x) and its corresponding RGB color values (r, g, b).
    The XML file is expected to follow a specific structure with <Point> elements nested inside a <ColorMap> element, which is itself part of a <ColorMaps> root element.

The function parses the XML file, extracts the x, r, g, and b attributes of each <Point> element, and constructs the colormap.
Parameters

    filename: A string representing the path to the XML file containing the colormap definition.

Returns

A Colormap, which is defined as:

using Colormap = std::vector<std::pair<double, std::array<double, 3>>>;

Each element of the colormap contains:

    A double representing the position x in the range [0.0,1.0].
    A std::array<double, 3> containing the RGB components of the color, each in the range [0.0,1.0].

Behavior

    XML Parsing:
        The function uses tinyxml2 to parse the provided XML file.
        It checks for the presence of a <ColorMaps> root element and a <ColorMap> child element.
    Colormap Extraction:
        Inside the <ColorMap> element, the function iterates through all <Point> elements.
        It reads the x, r, g, and b attributes of each <Point> and adds them to the colormap as a key-value pair.
    Validation:
        If the XML file cannot be loaded or is malformed, an exception is thrown.
        If no valid points are found, the function throws an exception indicating an empty colormap.

Exceptions

    std::runtime_error if:
        The XML file cannot be opened.
        The file does not contain the required <ColorMaps> or <ColorMap> elements.
        No <Point> elements are found.
        Parsing errors occur for the attributes of any <Point>.
*/
static Colormap read_colormap(const std::string& filename) {
    // Vector to store the colormap points (x, [r, g, b]).
    Colormap colormap;

    // Load the XML document.
    tinyxml2::XMLDocument doc;
    if (doc.LoadFile(filename.c_str()) != tinyxml2::XML_SUCCESS) {
        throw std::runtime_error("Failed to load colormap file: " + filename);
    }

    // Find the root element.
    tinyxml2::XMLElement* root = doc.FirstChildElement("ColorMaps");
    if (!root) {
        throw std::runtime_error("Invalid colormap file: Missing <ColorMaps> root element.");
    }

    // Iterate over all <Colormap> elements.
    tinyxml2::XMLElement* colormapElement = root->FirstChildElement("ColorMap");
    if (!colormapElement) {
        throw std::runtime_error("Invalid colormap file: Missing <ColorMap> element.");
    }

    // Iterate over all <Point> elements inside the <Colormap>.
    tinyxml2::XMLElement* pointElement = colormapElement->FirstChildElement("Point");
    while (pointElement) {
        double x = 0.0, r = 0.0, g = 0.0, b = 0.0;

        // Parse attributes.
        pointElement->QueryDoubleAttribute("x", &x);
        pointElement->QueryDoubleAttribute("r", &r);
        pointElement->QueryDoubleAttribute("g", &g);
        pointElement->QueryDoubleAttribute("b", &b);

        // Add the point to the colormap.
        colormap.emplace_back(x, std::array<double, 3>{r, g, b});

        // Move to the next <Point>.
        pointElement = pointElement->NextSiblingElement("Point");
    }

    if (colormap.empty()) {
        throw std::runtime_error("No points found in the colormap.");
    }

    return colormap;
}

// Private method for managing export
int Fem2d::handleRGBexport(const Settings &settings, ExportType eType, const std::function<double(const Node2d&)>& valueExtractor) {
    const int NOD = node.size();
    const int N = static_cast<int>(sqrt(static_cast<double>(NOD)));
    if (NOD != N * N) {
        throw std::logic_error("Error: NOD is not equal to N squared.");
    }

    // Defining bounds
    double xmin = c[0] - l[0] / 2.;
    double xmax = c[0] + l[0] / 2.;
    double ymin = c[1] - l[1] / 2.;
    double ymax = c[1] + l[1] / 2.;

    // Allocating the pixel array
    std::vector<std::vector<double>> pixels(N, std::vector<double>(N));

    for (int np = 0; np < NOD; ++np) {
        Node2d& node_ = node[np];
        double x = node_.p[0];
        double y = node_.p[1];

        uint32_t i = static_cast<uint32_t>((N - 1) * (x - xmin) / (xmax - xmin) + 0.5);
        uint32_t j = static_cast<uint32_t>((N - 1) * (y - ymin) / (ymax - ymin) + 0.5);
        pixels[N-1-j][i] = valueExtractor(node_);
    }

    std::string filename, unit;

    switch (eType) {
   	case ExportType::CONTRAST:
		filename=settings.getSimName()+"_STXM_XMCD.png";
		unit = "dimless";
		break;
	case ExportType::MZ_INTEGRAL:
		filename=settings.getSimName()+"_MZ.png";
		unit = "A";
		break;
	case ExportType::PATH_LENGTH:
		filename=settings.getSimName()+"_PATH_LENGTH.png";
		unit = "m";
		break;
	case ExportType::HOLO_PHASE:
		filename=settings.getSimName()+"_HOLO_PHASE.png";
		unit = "rad";
		break;
	default:
		std::cout << termcolor::bright_red << termcolor::blink << "Error export type unknown" << termcolor::reset << std::endl;
		exit(1);
		}

    // Creating image with its metadata
    Image image{
        .width = N,
        .height = N,
        .pixels = pixels,
        .pixel_size = meshSize,
        .unit = unit,
        .metadata = {
            {"Software", "pathIntegral"},
            {"Version", pathIntegral_version},
            {"Simulation parameters", settings.yaml_source}
        }
    };

    // Saving image
    std::ofstream fout(filename, std::ios::binary);
    if (fout.fail()) {
        std::cerr << "Error: Unable to open the file " << filename << std::endl;
        return 1;
    }

    std::filesystem::path ColormapFile = "HoloScale.xml";

    Colormap colormap;
    if (std::filesystem::exists(ColormapFile)) {
       colormap = read_colormap(ColormapFile);
    } else {
        std::cout << "Colormap File " <<  ColormapFile << "is missing!"<< std::endl;
        std::cout << "replaced by jet colormap" << std::endl;
        colormap = jet_colors;
    }
    write_png_image(fout, image, colormap);
    
    fout.close();

    std::cout << "Image saved in " << filename << std::endl;
    return 0;
}

int Fem2d::exportRatioRGBscaleImage(const Settings &settings, ExportType eType) {
    // Map of extraction functions associated with each type
    static const std::unordered_map<ExportType, std::function<double(const Node2d&)>> valueExtractors = {
        {ExportType::CONTRAST,    [](const Node2d& node) { return node.contrast; }},
        {ExportType::MZ_INTEGRAL, [](const Node2d& node) { return node.Mz_integral; }},
        {ExportType::PATH_LENGTH, [](const Node2d& node) { return node.path_length; }},
        {ExportType::HOLO_PHASE, [this](const Node2d& node) { 
        	double phase = this->CE * this->V*node.path_length-CHARGE_ELECTRON/PLANCKS_HBAR*node.sol;
        return phase; }}
    };

    // Check if the type is supported
    auto it = valueExtractors.find(eType);
    if (it == valueExtractors.end()) {
        std::cerr << termcolor::bright_red << "Error: Export type unknown " << int(eType)<< termcolor::reset << std::endl;
        return 1;
    }

    // Apply the right extractor
    return handleRGBexport(settings, eType, it->second);
}

