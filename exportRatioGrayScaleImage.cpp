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

// Private method for managing export
int Fem2d::handleExport(const Settings &settings, ExportType eType, const std::function<double(const Node2d&)>& valueExtractor) {
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

    write_png_image(fout, image);
    fout.close();

    std::cout << "Image saved in " << filename << std::endl;
    return 0;
}

int Fem2d::exportRatioGrayScaleImage(const Settings &settings, ExportType eType) {
    // Map of extraction functions associated with each type
    static const std::unordered_map<ExportType, std::function<double(const Node2d&)>> valueExtractors = {
        {ExportType::CONTRAST,    [](const Node2d& node) { return node.contrast; }},
        {ExportType::MZ_INTEGRAL, [](const Node2d& node) { return node.Mz_integral; }},
        {ExportType::PATH_LENGTH, [](const Node2d& node) { return node.path_length; }}
    };

    // Check if the type is supported
    auto it = valueExtractors.find(eType);
    if (it == valueExtractors.end()) {
        std::cerr << termcolor::bright_red << "Error: Export type unknown " << int(eType)<< termcolor::reset << std::endl;
        return 1;
    }

    // Apply the right extractor
    return handleExport(settings, eType, it->second);
}

