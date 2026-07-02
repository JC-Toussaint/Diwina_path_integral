#include <cstdint>
#include <limits>
#include <vector>
#include <filesystem>
#include "termcolor.h"

#include <cmath>
#include <string>
#include <png.h>
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

// Read a colormap from a file.
// The file format is text, with four values per line: (x, r, g, b).
// All values should be between 0 and 1.
// Empty lines and lines starting with '#' are comments and are ignored.
static Colormap read_colormap(const std::string& filename) {
    std::ifstream f(filename);
    if (!f.is_open())
        throw std::ios_base::failure(filename);
    std::string line;
    Colormap colormap;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#')  // skip comments
            continue;
        std::istringstream input(line);
        double x, r, g, b;
        input >> x >> r >> g >> b;
        ColormapEntry entry = {x, {r, g, b}};
        colormap.push_back(entry);
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
		filename=settings.getSimName()+"_STXM_XMCD_RGB.png";
		unit = "dimless";
		break;
	case ExportType::MZ_INTEGRAL:
		filename=settings.getSimName()+"_MZ_RGB.png";
		unit = "A";
		break;
	case ExportType::PATH_LENGTH:
		filename=settings.getSimName()+"_PATH_LENGTH_RGB.png";
		unit = "m";
		break;
	case ExportType::HOLO_PHASE:
		filename=settings.getSimName()+"_HOLO_PHASE_RGB.png";
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
        .pixel_size = pixel_size,
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

    std::filesystem::path ColormapFile = "HoloScale.tsv";

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

