#ifndef WRITE_PNG_H
#define WRITE_PNG_H

#include <array>
#include <vector>
#include <string>
#include <fstream>

// A Colormap is a tabulated function from [0, 1] to [0, 1]^3.
// Entries should be sorted in ascendig order of the independent variable, which should be 0 (or 1)
// on the first (resp. last) entry.
using ColormapEntry = std::pair<double, std::array<double, 3>>;
using Colormap = std::vector<ColormapEntry>;

// Structure representing an image.
struct Image {
    using MetadataItem = std::pair<std::string, std::string>;

    int width;                                // Image width
    int height;                               // Image height
    std::vector<std::vector<double>> pixels;  // Pixel matrix (double)
    double pixel_size;                        // Physical size of a pixel
    std::string unit;                         // Pixel unit (e.g., "m")
    std::vector<MetadataItem> metadata;
};

// Function to write a PNG image.
// If the colormap is empty, write a 16-bit grayscale image,
// otherwise write an 8-bit palette-based image.
void write_png_image(std::ofstream& f, const Image& image, const Colormap &colormap = Colormap());

#endif // WRITE_PNG_H

