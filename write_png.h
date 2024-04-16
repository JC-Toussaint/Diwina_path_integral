#ifndef WRITE_PNG_H
#define WRITE_PNG_H

#include <vector>
#include <string>
#include <fstream>

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
void write_png_image(std::ofstream& f, const Image& image);

#endif // WRITE_PNG_H

