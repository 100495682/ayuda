// imgsoa.cpp
#include "imgsoa.hpp"
#include "common/progargs.hpp"

#include <fstream>
#include <iostream>
#include <map>
#include <limits>
#include <cmath>
#include <algorithm>
#include <string>

namespace imgsoa {

  //COMPRESS

    Image::Image(dimensions dim_init, int max_color_init)
        : dim(dim_init), max_color(max_color_init),
          r(static_cast<std::vector<uint16_t>::size_type>(dim.height * dim.width)),
          g(static_cast<std::vector<uint16_t>::size_type>(dim.height * dim.width)),
          b(static_cast<std::vector<uint16_t>::size_type>(dim.height * dim.width)) {}

    std::map<std::tuple<uint16_t, uint16_t, uint16_t>, uint32_t> leerColores(Image& image) {
        std::map<std::tuple<uint16_t, uint16_t, uint16_t>, uint32_t> mapa_colores;

        for (std::size_t i = 0; i < image.r.size(); i++) {
            std::tuple<uint16_t, uint16_t, uint16_t> const color = std::make_tuple(image.r[i], image.g[i], image.b[i]);

            if (mapa_colores.find(color) == mapa_colores.end()) {
                mapa_colores[color] = 0;
                ++image.number_colors;
            }
        }

        uint32_t indice = 0;
        for (auto& [color, index] : mapa_colores) {
            mapa_colores[color] = indice;
            ++indice;
        }

        return mapa_colores;
    }

    void escribirColor(const std::tuple<uint16_t, uint16_t, uint16_t>& color, std::ofstream& file, int max_color) {
      if (max_color <= MAX_COLOR_1) {
        auto const red = static_cast<uint8_t>(std::get<0>(color));
        auto const green = static_cast<uint8_t>(std::get<1>(color));
        auto const blue = static_cast<uint8_t>(std::get<2>(color));

        binaryo::write_binary(file, red);
        binaryo::write_binary(file, green);
        binaryo::write_binary(file, blue);
      } else {
        auto const red = std::get<0>(color);
        auto const green = std::get<1>(color);
        auto const blue = std::get<2>(color);

        binaryo::write_binary(file, red);
        binaryo::write_binary(file, green);
        binaryo::write_binary(file, blue);
      }
    }

    void escribirTabla(Image &img,
                       std::map<std::tuple<uint16_t, uint16_t, uint16_t>, uint32_t> &listaColor,
                       std::string &archivo) {
        std::ofstream file(archivo, std::ios::binary | std::ios::app);

        if (!file.is_open()) {
            std::cerr << "Error al abrir el archivo " << archivo << "\n";
            return;
        }

        for (const auto& [color, index] : listaColor) {
            escribirColor(color, file, img.max_color);
        }
    }


    void escribirIndice(uint32_t indice, std::ofstream& file, uint32_t number_colors) {
        if (number_colors <= MAX_COLOR_1) {
            auto const indice_8bit = static_cast<uint8_t>(indice);
            binaryo::write_binary(file, indice_8bit);
        } else if (number_colors <= MAX_COLOR_2) {
            auto const indice_16bit = static_cast<uint16_t>(indice);
            binaryo::write_binary(file, indice_16bit);
        } else {
            binaryo::write_binary(file, indice);
        }
    }


    void escribirColores(Image &img,
                         std::map<std::tuple<uint16_t, uint16_t, uint16_t>, uint32_t> &listaColor,
                         std::string &archivo) {
        std::ofstream file(archivo, std::ios::binary | std::ios::app);

        if (!file.is_open()) {
            std::cerr << "Error al abrir el archivo " << archivo << "\n";
            return;
        }

        for (std::size_t i = 0; i < img.r.size(); i++) {
            std::tuple<uint16_t, uint16_t, uint16_t> const color = std::make_tuple(img.r[i], img.g[i], img.b[i]);
            const uint32_t indice = listaColor[color];
            escribirIndice(indice, file, img.number_colors);
        }
    }


    void compress(progargsCommon::parameters_files & params) {
        std::ifstream input(params.input_file , std::ios::binary);
        if (!input.is_open()) {
            std::cerr << "Error al abrir el archivo de entrada\n";
            return;
        }
        std::string magic_num;
        int width = 0;
        int height = 0;
        int max_color = 0;
        input >> magic_num >> width >> height >> max_color;
        input.ignore(1); // Ignorar el carácter de nueva línea después del encabezado

        dimensions const dim_init{.height = height, .width = width}; // Inicialización correcta
        Image image(dim_init, max_color);

        for (std::size_t i = 0; i < image.r.size(); i++) {
            if (image.max_color <= MAX_COLOR_1) {
                image.r[i] = binaryo::read_binary<uint8_t>(input);
                image.g[i] = binaryo::read_binary<uint8_t>(input);
                image.b[i] = binaryo::read_binary<uint8_t>(input);
            } else {
                image.r[i] = binaryo::read_binary<uint16_t>(input);
                image.g[i] = binaryo::read_binary<uint16_t>(input);
                image.b[i] = binaryo::read_binary<uint16_t>(input);
            }
        }
        input.close();
        std::map<std::tuple<uint16_t, uint16_t, uint16_t>, uint32_t> listaColor = leerColores(image);

        std::ofstream output(params.output_file, std::ios::binary);
        if (!output.is_open()) {
            std::cerr << "Error al abrir el archivo de salida\n";
            return;
        }
        output << "C6 " << image.dim.width << " " << image.dim.height << " " << image.max_color << " " << image.number_colors << "\n";
        output.close();
        escribirTabla(image, listaColor, params.output_file);
        escribirColores(image, listaColor, params.output_file);
        std::cout << "Compresión exitosa\n";
    }

  //RESIZE
    inline double interpolate(double start, double end, double factor) {
          return start + (factor * (end - start));
      }

      uint8_t bilinearInterpolate(
              const NeighboringPixels& neighbors,
              double x_factor,
              double y_factor) {
          double const top = interpolate(neighbors.top_left, neighbors.top_right, x_factor);
          double const bottom = interpolate(neighbors.bottom_left, neighbors.bottom_right, x_factor);
          return static_cast<uint8_t>(std::round(interpolate(top, bottom, y_factor)));
      }

      NeighboringPixels getNeighboringPixels(
              const ImageSOA& image,
              const PixelCoordinates& pixel_coords) {
          const std::vector<uint8_t>* channel_data = nullptr;

          switch (pixel_coords.channel) {
              case Channel::Red:
                  channel_data = &image.r;
                  break;
              case Channel::Green:
                  channel_data = &image.g;
                  break;
              case Channel::Blue:
                  channel_data = &image.b;
                  break;
              default:
                  throw std::invalid_argument("Invalid channel in getNeighboringPixels.");
          }

          const std::size_t width = image.width;
          const Coordinates& coords = pixel_coords.coords;

          NeighboringPixels neighbors{};
          neighbors.top_left = (*channel_data)[(coords.y_low * width) + coords.x_low];
          neighbors.top_right = (*channel_data)[(coords.y_low * width) + coords.x_high];
          neighbors.bottom_left = (*channel_data)[(coords.y_high * width) + coords.x_low];
          neighbors.bottom_right = (*channel_data)[(coords.y_high * width) + coords.x_high];

          return neighbors;
      }

      std::pair<double, double> calculateSourceCoordinates(
              std::size_t x_cord,
              std::size_t y_cord,
              const ImageSOA& original_image,
              const ImageSOA& resized_image) {
          double const src_x = (resized_image.width > 1)
                               ? (static_cast<double>(x_cord) * (static_cast<double>(original_image.width) - 1.0)) /
                                 (static_cast<double>(resized_image.width) - 1.0)
                               : 0.0;

          double const src_y = (resized_image.height > 1)
                               ? (static_cast<double>(y_cord) * (static_cast<double>(original_image.height) - 1.0)) /
                                 (static_cast<double>(resized_image.height) - 1.0)
                               : 0.0;

          return {src_x, src_y};
      }

      void resize(const ImageSOA& original_image, ImageSOA& resized_image) {
          // Resizing logic remains the same as previous resizeSOA implementation

          for (std::size_t y_cord = 0; y_cord < resized_image.height; ++y_cord) {
              for (std::size_t x_cord = 0; x_cord < resized_image.width; ++x_cord) {
                  auto [src_x, src_y] = calculateSourceCoordinates(
                          x_cord, y_cord, original_image, resized_image);

                  double const x_factor = src_x - std::floor(src_x);
                  double const y_factor = src_y - std::floor(src_y);

                  auto x_low = static_cast<std::size_t>(std::floor(src_x));
                  auto y_low = static_cast<std::size_t>(std::floor(src_y));
                  std::size_t const x_high = std::min(x_low + 1, original_image.width - 1);
                  std::size_t const y_high = std::min(y_low + 1, original_image.height - 1);

                  // Initialize Coordinates struct
                  Coordinates coords{};
                  coords.x_low = x_low;
                  coords.y_low = y_low;
                  coords.x_high = x_high;
                  coords.y_high = y_high;

                  // Initialize PixelCoordinates structs
                  PixelCoordinates const coords_r{coords, Channel::Red};
                  PixelCoordinates const coords_g{coords, Channel::Green};
                  PixelCoordinates const coords_b{coords, Channel::Blue};

                  NeighboringPixels const neighbors_r = getNeighboringPixels(original_image, coords_r);
                  NeighboringPixels const neighbors_g = getNeighboringPixels(original_image, coords_g);
                  NeighboringPixels const neighbors_b = getNeighboringPixels(original_image, coords_b);

                  std::size_t const index = (y_cord * resized_image.width) + x_cord;

                  resized_image.r[index] = bilinearInterpolate(neighbors_r, x_factor, y_factor);
                  resized_image.g[index] = bilinearInterpolate(neighbors_g, x_factor, y_factor);
                  resized_image.b[index] = bilinearInterpolate(neighbors_b, x_factor, y_factor);
              }
          }
      }

      void resize(const std::vector<std::string> &args){
          // Read the input image
          ImageSOA const original_image = readPPM(args[1]);

          // Create a resized image object
          ImageSOA resized_image(static_cast<std::size_t>(std::stoi(args[4])), static_cast<std::size_t>(std::stoi(args[progargsCommon::CINCO])));        // Perform the resizing operation
          resize(original_image, resized_image);

          // Write the resized image to the output file
          writePPM(args[2], resized_image);
      }

      ImageSOA readPPM(const std::string& filename) {
          std::ifstream file(filename, std::ios::binary);
          if (!file) {
              throw std::runtime_error("Could not open file: " + filename);
          }

          std::string format;
          std::size_t width = 0;
          std::size_t height = 0;
          int max_color = 0;

          file >> format;
          if (format != "P6") {
              throw std::runtime_error("Unsupported PPM format (expected P6)");
          }

          file >> width >> height >> max_color;
          file.ignore(1); // Skip the newline character

          if (max_color != MAX_COLOR_VALUE) {
              throw std::runtime_error("Unsupported max color value (expected 255)");
          }

          ImageSOA image(width, height);

          for (std::size_t i = 0; i < width * height; ++i) {
              image.r[i] = binaryo::read_binary<uint8_t>(file);
              image.g[i] = binaryo::read_binary<uint8_t>(file);
              image.b[i] = binaryo::read_binary<uint8_t>(file);
          }

          if (!file) {
              throw std::runtime_error("Error reading pixel data from file");
          }

          return image;
      }

      void writePPM(const std::string& filename, const ImageSOA& image) {
          std::ofstream file(filename, std::ios::binary);
          if (!file) {
              throw std::runtime_error("Could not open file for writing: " + filename);
          }

          file << "P6\n" << image.width << " " << image.height << "\n" << MAX_COLOR_VALUE << "\n";

          for (std::size_t i = 0; i < image.r.size(); ++i) {
              binaryo::write_binary(file, image.r[i]);
              binaryo::write_binary(file, image.g[i]);
              binaryo::write_binary(file, image.b[i]);
          }

          if (!file) {
              throw std::runtime_error("Error writing pixel data to file");
          }
      }


/*
    bool validateAndParseArgs(const std::vector<std::string>& args, progargsMaxLevel::parameters& params) {
    if (args.size() != 4) {
        std::cerr << "Usage: maxlevel <input_file> <output_file> <multiplier>\n";
        return false;
    }

    params.input_file = args[1];
    params.output_file = args[2];

    try {
        params.multiplier = std::stoi(args[3]);
      if (constexpr int max_multiplier = 65535;
          params.multiplier < 0 || params.multiplier > max_multiplier) {
            std::cerr << "Error: Multiplier must be between 0 and " << max_multiplier << ".\n";
            return false;
        }
    } catch (...) {
        std::cerr << "Error: Invalid multiplier value.\n";
        return false;
    }

    return true;
}

// Function to read the header and validate the PPM file
bool readAndValidateHeader(std::ifstream& input_file, PPMInfo& ppm_info) {
    std::string magic_number;
    input_file >> magic_number;
    constexpr int upper_limit_pixels = 65535;

    if (magic_number != "P6") {
        std::cerr << "Error: Unsupported PPM format. Expected P6 format.\n";
        return false;
    }

    input_file >> ppm_info.width >> ppm_info.height >> ppm_info.max_color_value;

    if (ppm_info.max_color_value > upper_limit_pixels) {
        std::cerr << "Error: Max color value exceeds supported range.\n";
        return false;
    }

    input_file.ignore(); // Skip newline
    return true;
}

// Function to read pixel data
bool readPixelData(std::ifstream& input_file, const PPMInfo& ppm_info, RGBData& rgb_data) {
    size_t const total_pixels =
        static_cast<size_t>(std::max(0, ppm_info.width)) * static_cast<size_t>(std::max(0, ppm_info.height));

    rgb_data.red.resize(total_pixels);
    rgb_data.green.resize(total_pixels);
    rgb_data.blue.resize(total_pixels);

    constexpr int og_pixel_max = 255;
    for (size_t i = 0; i < total_pixels; ++i) {
        if (ppm_info.max_color_value <= og_pixel_max) {
            rgb_data.red[i] = binaryo::read_binary<uint8_t>(input_file);
            rgb_data.green[i] = binaryo::read_binary<uint8_t>(input_file);
            rgb_data.blue[i] = binaryo::read_binary<uint8_t>(input_file);
        } else {
            constexpr int num_bits = 8;
            auto red = binaryo::read_binary<uint16_t>(input_file);
            auto green = binaryo::read_binary<uint16_t>(input_file);
            auto blue = binaryo::read_binary<uint16_t>(input_file);

            rgb_data.red[i] = (red >> num_bits) | (red << num_bits);
            rgb_data.green[i] = (green >> num_bits) | (green << num_bits);
            rgb_data.blue[i] = (blue >> num_bits) | (blue << num_bits);
        }
    }

    return true;
}

// Function to calculate new pixel intensities
void calculateNewIntensities(
    const RGBData& rgb_data,
    RGBData& new_rgb_data,
    int multiplier) {
    size_t const total_pixels = rgb_data.red.size();
    new_rgb_data.red.resize(total_pixels);
    new_rgb_data.green.resize(total_pixels);
    new_rgb_data.blue.resize(total_pixels);

    for (size_t i = 0; i < total_pixels; ++i) {
      constexpr int original_max = 255;
      new_rgb_data.red[i] = static_cast<uint16_t>((rgb_data.red[i] * multiplier) / original_max);
        new_rgb_data.green[i] = static_cast<uint16_t>((rgb_data.green[i] * multiplier) / original_max);
        new_rgb_data.blue[i] = static_cast<uint16_t>((rgb_data.blue[i] * multiplier) / original_max);
    }
}

void writePixelData(
    const progargsMaxLevel::parameters& params,
    const RGBData& rgb_data,
    const PPMInfo& ppm_info,
    int new_max_color_value) {

    // Ensure the directory for the output file exists
    std::filesystem::path const output_path(params.output_file);

    if (std::filesystem::path const directory = output_path.parent_path();
        !std::filesystem::exists(directory)) {
        std::cerr << "Directory " << directory << " does not exist. Creating it.\n";
        try {
            std::filesystem::create_directories(directory);
        } catch (const std::filesystem::filesystem_error& e) {
            std::cerr << "Error: Could not create directory " << directory << ": " << e.what() << "\n";
            return;
        }
    }

    // Open the output file in binary mode
    std::ofstream output(params.output_file, std::ios::binary);
    if (!output.is_open()) {
        std::cerr << "Error: Could not open file " << params.output_file << " for writing.\n";
        return;
    }

    // Write PPM header
    output << "P6\n" << ppm_info.width << " " << ppm_info.height << "\n" << new_max_color_value << "\n";

    // Write pixel data
    for (size_t i = 0; i < rgb_data.red.size(); ++i) {
        constexpr uint16_t max_pixel_value = 255;
        constexpr uint16_t min_pixel_value = 0;
        output.put(static_cast<char>(std::clamp(rgb_data.red[i], min_pixel_value, max_pixel_value)));
        output.put(static_cast<char>(std::clamp(rgb_data.green[i], min_pixel_value, max_pixel_value)));
        output.put(static_cast<char>(std::clamp(rgb_data.blue[i], min_pixel_value, max_pixel_value)));
    }

    output.close();
    std::cout << "Output successfully written to: " << params.output_file << "\n";
}

// Main maxlevel function
int maxlevel(const std::vector<std::string>& args) {
    progargsMaxLevel::parameters params;

    if (!validateAndParseArgs(args, params)) {
        return 1;
    }

    std::ifstream input_file(params.input_file, std::ios::binary);
    if (!input_file.is_open()) {
        std::cerr << "Error: Could not open file " << params.input_file << "\n";
        return 1;
    }

    PPMInfo ppm_info;
    if (!readAndValidateHeader(input_file, ppm_info)) {
        return 1;
    }

    RGBData rgb_data;
    if (!readPixelData(input_file, ppm_info, rgb_data)) {
        return 1;
    }

    input_file.close();

    RGBData new_rgb_data;
    calculateNewIntensities(rgb_data, new_rgb_data, params.multiplier);


    const int new_max_color_value = (params.multiplier > 255) ? 65535 : 255;
    writePixelData(params, new_rgb_data, ppm_info, new_max_color_value);

    std::cout << "Successfully wrote modified pixel data to " << params.output_file << "\n";
    return 0;
}
*/

/*
 *
int euclidean_distance_squared(const pixel_soa& pix1, const pixel_soa& pix2, size_t index) {
  int const dr1 = static_cast<int>(pix1.r[index]) - static_cast<int>(pix2.r[index]);
  int const dg2 = static_cast<int>(pix1.g[index]) - static_cast<int>(pix2.g[index]);
  int const db3 = static_cast<int>(pix1.b[index]) - static_cast<int>(pix2.b[index]);
  return (dr1 * dr1) + (dg2 * dg2) + (db3 * db3);
}

std::tuple<uint8_t, uint8_t, uint8_t> min_distance_vector(const pixel_soa& pix1, const pixel_soa& pix2) {
  int min_dist = std::numeric_limits<int>::max();
  uint8_t min_r = 0;
  uint8_t min_g = 0;
  uint8_t min_b = 0;

  for (size_t i = 0; i < pix2.r.size(); ++i) {
    int const dist = euclidean_distance_squared(pix1, pix2, i);
    if (dist == 1) {
      return {pix2.r[i], pix2.g[i], pix2.b[i]};
    }
    if (dist < min_dist) {
      min_dist = dist;
      min_r = pix2.r[i];
      min_g = pix2.g[i];
      min_b = pix2.b[i];
    }
  }

  return {min_r, min_g, min_b};
}

pixel_soa read_image(const std::string& filename) {
  std::ifstream file(filename, std::ios::binary);
  if (!file.is_open()) { std::cerr << "Error opening image file.\n"; }

  std::string format;
  unsigned int width = 0;
  unsigned int height = 0;
  unsigned int max_value = 0;
  file >> format >> width >> height >> max_value;
  file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  pixel_soa pixsoa;

  pixsoa.r.reserve(static_cast<int>(width * height));
  pixsoa.g.reserve(static_cast<int>(width * height));
  pixsoa.b.reserve(static_cast<int>(width * height));

  for (size_t i = 0; i < static_cast<int>(width * height); ++i) {
    pixsoa.r.push_back(static_cast<uint8_t>(file.get()));
    pixsoa.g.push_back(static_cast<uint8_t>(file.get()));
    pixsoa.b.push_back(static_cast<uint8_t>(file.get()));
  }

  return pixsoa;
}

std::map<std::tuple<uint8_t, uint8_t, uint8_t>, int> read_frequency_map(const std::string& filename) {
  std::ifstream file(filename, std::ios::binary);
  if (!file.is_open()) { std::cerr << "Error opening image file.\n"; }

  std::string format;
  unsigned int width = 0;
  unsigned int height = 0;
  unsigned int max_value = 0;
  file >> format >> width >> height >> max_value;
  file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  std::map<std::tuple<uint8_t, uint8_t, uint8_t>, int> freq_map;
  for (size_t i = 0; i < static_cast<int>(width * height); ++i) {
    auto red = static_cast<uint8_t>(file.get());
    auto green = static_cast<uint8_t>(file.get());
    auto blue = static_cast<uint8_t>(file.get());
    freq_map[{red, green, blue}]++;
  }

  return freq_map;
}

std::vector<std::tuple<uint8_t, uint8_t, uint8_t>> extract_sorted(const std::map<std::tuple<uint8_t, uint8_t, uint8_t>, int>& freq_map) {
  std::vector<std::pair<std::tuple<uint8_t, uint8_t, uint8_t>, int>> sorted_vec(freq_map.begin(), freq_map.end());
  std::ranges::sort(sorted_vec, [](auto const& aaa, auto const& bbb) { return aaa.second < bbb.second; });

  std::vector<std::tuple<uint8_t, uint8_t, uint8_t>> result;
  result.reserve(sorted_vec.size());
  for (const auto& [pixel, _] : sorted_vec) {
    result.push_back(pixel);
  }
  return result;
}

void choose_eliminated(size_t max_colors, const std::map<std::tuple<uint8_t, uint8_t, uint8_t>, int>& freq_map,
                       std::vector<std::tuple<uint8_t, uint8_t, uint8_t>>& to_eliminate,
                       std::vector<std::tuple<uint8_t, uint8_t, uint8_t>>& to_keep) {
  to_eliminate.reserve(max_colors);
  to_keep.reserve(freq_map.size());

  for (const auto& [pixel, _] : freq_map) {
    if (to_eliminate.size() < max_colors) {
      to_eliminate.push_back(pixel);
    } else {
      to_keep.push_back(pixel);
    }
  }
}

std::vector<std::pair<std::tuple<uint8_t, uint8_t, uint8_t>, std::tuple<uint8_t, uint8_t, uint8_t>>>
calculate_new_colors(const std::vector<std::tuple<uint8_t, uint8_t, uint8_t>>& to_eliminate,
                     const pixel_soa& pixsoa) {
  std::vector<std::pair<std::tuple<uint8_t, uint8_t, uint8_t>, std::tuple<uint8_t, uint8_t, uint8_t>>> pairs;
  pairs.reserve(to_eliminate.size());

  for (const auto& pixel : to_eliminate) {
    auto new_color = min_distance_vector(pixsoa, pixsoa);
    pairs.emplace_back(pixel, new_color);
  }
  return pairs;
}

void replace_pixels(pixel_soa& pixsoa, const std::vector<std::pair<std::tuple<uint8_t, uint8_t, uint8_t>, std::tuple<uint8_t, uint8_t, uint8_t>>>& pairs) {
  std::map<std::tuple<uint8_t, uint8_t, uint8_t>, std::tuple<uint8_t, uint8_t, uint8_t>> replacements;
  for (const auto& [from, to] : pairs) {
    replacements[from] = to;
  }

  for (size_t i = 0; i < pixsoa.r.size(); ++i) {
    auto ite = replacements.find({pixsoa.r[i], pixsoa.g[i], pixsoa.b[i]});
    if (ite != replacements.end()) {
      pixsoa.r[i] = std::get<0>(ite->second);
      pixsoa.g[i] = std::get<1>(ite->second);
      pixsoa.b[i] = std::get<2>(ite->second);
    }
  }
}
auto parametros(const std::string&  entrada) {
  std::ifstream file(entrada, std::ios::binary);
  if (!file.is_open()) {std::cerr << "Error al abrir la imagen.\n";}
  std::string formato;
  file >> formato;
  int ancho = 0;
  int alto = 0;
  int luz = 0;
  file >> ancho >> alto >> luz;
  std::vector<int> parametros{ancho, alto,luz};
  return parametros;
}
void save_image(const std::string& filename, const pixel_soa& pixsoa, const std::vector<int>& parametros) {
  std::ofstream file(filename, std::ios::binary);
  std::vector<int> const& param = parametros;
  int const width   = param[0];
  int const height     = param[1];
  int const max_value       = param[2];
  file << "P6\n" << width << " " << height << "\n" << max_value << "\n";

  for (size_t i = 0; i < pixsoa.r.size(); ++i) {
    file.put(static_cast<char>(pixsoa.r[i]));
    file.put(static_cast<char>(pixsoa.g[i]));
    file.put(static_cast<char>(pixsoa.b[i]));
  }
}

  void cutFreq(const std::vector<std::string> &args) {
  pixel_soa pixsoa = read_image(args[1]);
  images const archivos = {.input_file=args[1],.output_file=args[2]};;

  std::map<std::tuple<uint8_t, uint8_t, uint8_t>, int> const freq_map = read_frequency_map(args[1]);
  std::vector<std::tuple<uint8_t, uint8_t, uint8_t>> const sorted_pixels = extract_sorted(freq_map);
  std::vector<std::tuple<uint8_t, uint8_t, uint8_t>> to_eliminate;
  std::vector<std::tuple<uint8_t, uint8_t, uint8_t>> to_keep;
  choose_eliminated(std::stoi(args[4]), freq_map, to_eliminate, to_keep);
  std::vector<std::pair<std::tuple<uint8_t, uint8_t, uint8_t>, std::tuple<uint8_t, uint8_t, uint8_t>>> const pairs =
      calculate_new_colors(to_eliminate, pixsoa);
  replace_pixels(pixsoa, pairs);
  std::vector<int> const paramos = parametros(args[1]);
  save_image(args[1], pixsoa, paramos);
}
*/

} // namespace imgsoa
