//
// Created by alejo on 14/11/24.
//

#ifndef DECOMPRESS_HPP
#define DECOMPRESS_HPP

#include <vector>
#include <iostream>
#include "common/binaryo.hpp"
#include "common/progargs.hpp"

namespace decompress {
  // Constantes para límites de colores
  constexpr int LONGITUD_CABECERA = 255;
  constexpr size_t MAX_COLORS_8BIT = 256;
  constexpr size_t MAX_COLORS_16BIT = 65536;
  constexpr size_t MAX_COLORS_32BIT = 4294967296;

  struct Color {
    uint16_t r;
    uint16_t g;
    uint16_t b;

    bool operator<(const Color& other) const;
  };

  struct Image {
    int width = 0;
    int height = 0;
    int max_color = 0;
    std::vector<Color> pixels;
  };

  struct CPPMData {
    Image img;
    std::vector<Color> color_table;
    std::vector<uint32_t> pixel_indices;
    int index_size;
  };

  bool abrirFile(const std::string& filename, std::ifstream& file);
  bool readCPPMHeader(std::ifstream& file, Image& img, size_t& color_table_size);
  bool readColorTable(std::ifstream& file, int max_color, std::vector<Color>& color_table, size_t color_table_size);
  bool readPixelIndices(std::ifstream& file, size_t num_pixels, std::vector<uint32_t>& pixel_indices, int index_size);
  bool writePPM(const std::string& filename, const Image& img, const std::vector<Color>& color_table, const std::vector<uint32_t>& pixel_indices);
  bool initializeImage(const std::string& input_filename, Image& img, size_t& color_table_size, std::vector<Color>& color_table);
  bool readPixelIndicesFromFile(const std::string& input_filename, size_t num_pixels, std::vector<uint32_t>& pixel_indices, int index_size);
  bool writeOutputFile(const std::string& filename, const Image& img, const std::vector<Color>& color_table, const std::vector<uint32_t>& pixel_indices);
  void decompress(const progargsCommon::parameters_files& params);
}

#endif //DECOMPRESS_HPP