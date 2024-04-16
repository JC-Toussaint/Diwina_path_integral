#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <cassert>
#include <png.h>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <cstring>

#include "write_png.h"

static std::string number_to_string(double x, int precision = 6) {
    std::ostringstream stream;
    stream << std::setprecision(precision) << x;
    return stream.str();
}

// Helper function to find the minimum and maximum in a 2D vector.
std::pair<double, double> find_min_max(const std::vector<std::vector<double>>& matrix) {
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::lowest();

    for (const auto& row : matrix) {
        for (double val : row) {
            if (val < min) min = val;
            if (val > max) max = val;
        }
    }
    return {min, max};
}

void write_png_image(std::ofstream& f, const Image& image) {
    // Convenience variables.
    int width = image.width;
    int height = image.height;
    const std::vector<std::vector<double>>& pixels = image.pixels;

    // Grab current time.
    std::time_t now = std::time(nullptr);
    std::tm* local_time = std::localtime(&now);

    // Set up libpng data structures.
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    if (!png_ptr) {
        throw std::runtime_error("png_create_write_struct() failed");
    }

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
        png_destroy_write_struct(&png_ptr, nullptr);
        throw std::runtime_error("png_create_info_struct() failed");
    }

    if (setjmp(png_jmpbuf(png_ptr))) {
        png_destroy_write_struct(&png_ptr, &info_ptr);
        throw std::runtime_error("libpng error with longjmp()");
    }

    // Find min/max of the whole image.
    auto [min, max] = find_min_max(pixels);

    // Stringify min/max/range and pixel_size.
    std::string str_min = number_to_string(min);
    std::string str_max = number_to_string(max);
    std::string str_range = number_to_string(max - min);
    std::string pixel_size = number_to_string(image.pixel_size);

    // Stringify the creation time.
    //
    // Although the PNG specification suggests formatting the time as
    // per section 5.2.14 of RFC 1123, which is a variation on RFC 822,
    // we prefer the RFC 3339 format, which is a more modern standard
    // and easier to parse. Note that:
    //
    // - We write out local time. As the time zone is specified, UTC
    //   time can be easily deduced. This conveys more information than
    //   writing only UTC time.
    //
    // - We use a space to separate the date and time parts, for better
    //   readability. This is explicitly allowed by RFC 3339, but not by
    //   ISO 8601.
    //
    // - RFC 3339 requires the time zone to be formatted as +hh:mm or
    //   -hh:mm, but this is not supported by strftime(). When then
    //   use the conversion specification "%z", which outputs +hhmm or
    //   -hhmm, then manually insert the colon.
    //
    // Sample formatted time: "2024-12-09 12:34:56+01:00".
    
    // Stringify the creation time.
    std::ostringstream time_stream;
    time_stream << std::put_time(local_time, "%Y-%m-%d %H:%M:%S%z");
    std::string str_creation_time = time_stream.str();
    str_creation_time.insert(str_creation_time.size() - 2, 1, ':');

    // Initialize the image.
    //    png_set_write_fn:
    //    this function from libpng allows you to specify a custom mechanism to write data to a file or stream.
    //    Instead of using libpng’s default file I/O functions, you provide a custom callback to handle the output.

    // png_ptr:
    //    this is the main libpng structure for handling PNG files. It contains all state and configuration information for the current PNG operation.

    // &f:
    // this is a pointer to the std::ofstream object (passed as a reference here) that will be used to write the PNG data.
    // It is stored as the user-defined pointer (io_ptr) in the png_struct.

    // Lambda function ([](png_structp png_ptr, png_bytep data, png_size_t length) { ... }):
    // This lambda acts as the write callback.
    // Whenever libpng wants to write data, it calls this lambda with:
    //   png_ptr: the libpng structure.
    //   data: a pointer to the raw bytes that need to be written.
    //   length: the number of bytes to write.
    //   the lambda does the following:
    //       - uses png_get_io_ptr(png_ptr) to retrieve the io_ptr (which is &f, the std::ofstream object).
    //       - casts the io_ptr to a std::ofstream* so it can call methods on it.
    //       - writes the data buffer to the stream using write. The reinterpret_cast<const char*> ensures that
    //       the data pointer is interpreted as a const char*, which is what std::ofstream::write expects.

    // nullptr:
    // This is the flush callback, which is not used in this example. Passing nullptr means there’s no custom flush operation.
    
    // Initialize the image.
    png_set_write_fn(png_ptr, &f, [](png_structp png_ptr, png_bytep data, png_size_t length) {
        auto* stream = static_cast<std::ofstream*>(png_get_io_ptr(png_ptr));
        stream->write(reinterpret_cast<const char*>(data), length);
    }, nullptr);

    png_set_IHDR(png_ptr, info_ptr, width, height, 16, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

    // Add physical scale.
    png_set_sCAL(png_ptr, info_ptr, PNG_SCALE_METER, image.pixel_size, image.pixel_size);

    // Add pixel calibration data.
    const int pCAL_nparams = 2;
    char* pCAL_params[2] = { const_cast<char*>(str_min.c_str()), const_cast<char*>(str_range.c_str()) };
    png_set_pCAL(png_ptr, info_ptr, "Physical", 0, UINT16_MAX, PNG_EQUATION_LINEAR, pCAL_nparams,
                 image.unit.c_str(), pCAL_params);

    // Add creation time as tIME (image last-modification time).
    png_time creation_time;
    png_convert_from_time_t(&creation_time, now);
    png_set_tIME(png_ptr, info_ptr, &creation_time);

    // Add textual metadata. This is redundant with the sCAL and pCAL
    // chunks defined above. We duplicate the information because tEXt
    // chunks are easier to extract from the PNG file.
    std::vector<Image::MetadataItem> metadata = {
        {"Pixel width", pixel_size},
        {"Pixel height", pixel_size},
        {"Pixel unit", "m"},
        {"Value min", str_min},
        {"Value max", str_max},
        {"Value unit", image.unit},
        {"Creation Time", str_creation_time}
    };
    metadata.insert(metadata.end(), image.metadata.begin(), image.metadata.end());

    std::vector<png_text> comments;
    for (const auto& item : metadata) {
        const std::string& key = item.first;
        const std::string& value = item.second;
        if (value.empty())  // skip empty values
            continue;
        bool isASCII = std::all_of(value.begin(), value.end(), [](char c){
            return (c & 0x80) == 0;
        });
        comments.push_back({
            .compression = isASCII ? PNG_TEXT_COMPRESSION_NONE : PNG_ITXT_COMPRESSION_NONE,
            .key = const_cast<char*>(key.c_str()),
            .text = const_cast<char*>(value.c_str()),
            .text_length = isASCII ? value.size() : 0,
            .itxt_length = isASCII ? 0 : value.size(),
            .lang = nullptr,
            .lang_key = nullptr
        });
    }
    png_set_text(png_ptr, info_ptr, comments.data(), comments.size());

    // Prepare image data.
    std::vector<std::vector<uint16_t>> png_pixels(height, std::vector<uint16_t>(width));
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            png_pixels[i][j] = static_cast<uint16_t>(((pixels[i][j] - min) / (max - min) * UINT16_MAX) + 0.5);
        }
    }

    std::vector<png_bytep> row_pointers(height);
    for (int i = 0; i < height; ++i) {
        row_pointers[i] = reinterpret_cast<png_bytep>(png_pixels[i].data());
    }
    png_set_rows(png_ptr, info_ptr, row_pointers.data());

    // Write the file.
    png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_SWAP_ENDIAN, nullptr);

    // Cleanup.
    png_destroy_write_struct(&png_ptr, &info_ptr);
}
