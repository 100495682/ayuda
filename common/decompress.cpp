//
// Created by alejo on 14/11/24.
//
#include "decompress.hpp"
#include <fstream>
#include <iostream>
#include "common/progargs.hpp"
#include <limits>

namespace decompress{

  bool abrirFile(const std::string& filename, std::ifstream& file) {
    file.open(filename, std::ios::binary);
    if (!file.is_open()) {
      std::cerr << "Error al abrir el archivo de entrada: " << filename << "\n";
      return false;
    }
    return true;
  }

 bool readCPPMHeader(std::ifstream& file, Image& img, size_t& color_table_size) {
    std::string magic_number;
    file >> magic_number;
    if (magic_number != "C6") {
        std::cerr << "Invalid magic number: " << magic_number << "\n";
        return false;
    }

    file >> img.width >> img.height >> img.max_color >> color_table_size;
    // Ignorar el resto de la línea de la cabecera
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    return true;
}

    bool readColorTable(std::ifstream& file, int max_color, std::vector<Color>& color_table, size_t color_table_size) {
        color_table.resize(color_table_size);
        for (size_t i = 0; i < color_table_size; ++i) {
            if (max_color <= decompress::LONGITUD_CABECERA) { // Usar 255 para determinar si es 8 bits
                color_table[i].r = binaryo::read_binary<uint8_t>(file);
                color_table[i].g = binaryo::read_binary<uint8_t>(file);
                color_table[i].b = binaryo::read_binary<uint8_t>(file);
            } else {
                color_table[i].r = binaryo::read_binary<uint16_t>(file);
                color_table[i].g = binaryo::read_binary<uint16_t>(file);
                color_table[i].b = binaryo::read_binary<uint16_t>(file);
            }

            // Verificar si la lectura fue exitosa
            if (file.eof() || file.fail()) {
                std::cerr << "Error: Unexpected end of file or read failure while reading color table.\n";
                return false;
            }
        }
        return true;
    }

    bool readPixelIndices(std::ifstream& file, size_t num_pixels, std::vector<uint32_t>& pixel_indices, int index_size) {
        pixel_indices.resize(num_pixels);
        for (size_t i = 0; i < num_pixels; ++i) {
            if (index_size == 1) {
                pixel_indices[i] = binaryo::read_binary<uint8_t>(file);
            } else if (index_size == 2) {
                pixel_indices[i] = binaryo::read_binary<uint16_t>(file);
            } else if (index_size == 4) {
                pixel_indices[i] = binaryo::read_binary<uint32_t>(file);
            } else {
                std::cerr << "Error: Unsupported index size: " << index_size << "\n";
                return false;
            }

            // Verificar si la lectura fue exitosa
            if (file.eof() || file.fail()) {
                std::cerr << "Error: Unexpected end of file or read failure while reading pixel indices.\n";
                return false;
            }
        }
        return true;
}

    bool writePPM(const std::string& filename, const Image& img, const std::vector<Color>& color_table, const std::vector<uint32_t>& pixel_indices) {
        std::ofstream file(filename, std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "Could not open PPM file for writing: " << filename << "\n";
            return false;
        }

        file << "P6\n" << img.width << " " << img.height << "\n" << img.max_color << "\n";
        for (const auto& index : pixel_indices) {
            if (index >= color_table.size()) {
                std::cerr << "Error: Pixel index " << index << " out of range.\n";
                file.close();
                return false;
            }

          const auto& color = color_table[index];
          if (img.max_color <= decompress::LONGITUD_CABECERA) { // Usar 255 para 8 bits

            const auto red = static_cast<uint8_t>(color.r);
            const auto green = static_cast<uint8_t>(color.g);
            const auto blue = static_cast<uint8_t>(color.b);

            binaryo::write_binary<uint8_t>(file, red);
            binaryo::write_binary<uint8_t>(file, green);
            binaryo::write_binary<uint8_t>(file, blue);

          } else {

            binaryo::write_binary<uint16_t>(file, color.r);
            binaryo::write_binary<uint16_t>(file, color.g);
            binaryo::write_binary<uint16_t>(file, color.b);
          }

        }
        file.close();
        return true;
}

    bool initializeImage(const std::string& input_filename, Image& img, size_t& color_table_size, std::vector<Color>& color_table) {
    std::ifstream file;
    if (!abrirFile(input_filename, file)) {
        std::cerr << "Error: Could not open input file: " << input_filename << "\n";
        return false;
    }

    if (!readCPPMHeader(file, img, color_table_size)) {
        std::cerr << "Error: Could not read CPPM header.\n";
        file.close();
        return false;
    }

    if (!readColorTable(file, img.max_color, color_table, color_table_size)) {
        std::cerr << "Error: Could not read color table.\n";
        file.close();
        return false;
    }

    file.close();
    return true;
}

  // imgaos.cpp

  bool readPixelIndicesFromFile(const std::string& input_filename, size_t num_pixels, std::vector<uint32_t>& pixel_indices, int index_size) {
  std::ifstream file;
  if (!abrirFile(input_filename, file)) {
    std::cerr << "Error: Could not open input file: " << input_filename << "\n";
    return false;
  }

  // Move the pointer to the end of the header and color table
  Image img;
  size_t color_table_size = 0;
  std::vector<Color> color_table; // Remove const qualifier
  if (!readCPPMHeader(file, img, color_table_size)) {
    std::cerr << "Error: Could not read CPPM header.\n";
    file.close();
    return false;
  }

  if (!readColorTable(file, img.max_color, color_table, color_table_size)) {
    std::cerr << "Error: Could not read color table.\n";
    file.close();
    return false;
  }

  if (!readPixelIndices(file, num_pixels, pixel_indices, index_size)) {
    std::cerr << "Error: Could not read pixel indices.\n";
    file.close();
    return false;
  }

  file.close();
  return true;
}

    bool writeOutputFile(const std::string& filename, const Image& img, const std::vector<Color>& color_table, const std::vector<uint32_t>& pixel_indices) {
        std::ofstream file(filename, std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "Could not open PPM file for writing: " << filename << "\n";
            return false;
        }

        file << "P6\n" << img.width << " " << img.height << "\n" << img.max_color << "\n";
        for (const auto& index : pixel_indices) {
            if (index >= color_table.size()) {
                std::cerr << "Error: Pixel index " << index << " out of range.\n";
                file.close();
                return false;
            }
            const auto& color = color_table[index];
            if (img.max_color <= decompress::LONGITUD_CABECERA) { // Usar 255 para 8 bits
                const auto red = static_cast<uint8_t>(color.r);
                const auto green = static_cast<uint8_t>(color.g);
                const auto blue = static_cast<uint8_t>(color.b);
                binaryo::write_binary<uint8_t>(file, red);
                binaryo::write_binary<uint8_t>(file, green);
                binaryo::write_binary<uint8_t>(file, blue);
            } else {
                binaryo::write_binary<uint16_t>(file, color.r);
                binaryo::write_binary<uint16_t>(file, color.g);
                binaryo::write_binary<uint16_t>(file, color.b);
            }
        }
        if (file.fail()) {
            std::cerr << "Error: Failed to write to PPM file.\n";
            file.close();
            return false;
        }
        file.close();
        return true;
    }

  void decompress(const progargsCommon::parameters_files& params) {
      Image img;
      size_t color_table_size = 0;
      std::vector<Color> color_table;

      // Initialize the image by reading the input file
      if (!initializeImage(params.input_file, img, color_table_size, color_table)) {
        return;
      }

      const size_t num_pixels = static_cast<size_t>(img.width) * static_cast<size_t>(img.height);
      int index_size = 0;
      if (color_table_size <= MAX_COLORS_8BIT) {
        index_size = 1;
      } else if (color_table_size <= MAX_COLORS_16BIT) {
        index_size = 2;
      } else {
        index_size = 4;
      }

      std::vector<uint32_t> pixel_indices;
      if (!readPixelIndicesFromFile(params.input_file, num_pixels, pixel_indices, index_size)) {
        return;
      }

      // Check if any pixel index is out of range
      for (const auto& index : pixel_indices) {
        if (index >= color_table.size()) {
          std::cerr << "Error: Pixel index " << index << " out of range. Color table size: " << color_table.size() << "\n";
          return;
        }
      }

      // Write the output PPM file
      if (!writeOutputFile(params.output_file, img, color_table, pixel_indices)) {
        return;
      }

      std::cout << "Decompression completed successfully. Output file: " << params.output_file << "\n";
    }
  }