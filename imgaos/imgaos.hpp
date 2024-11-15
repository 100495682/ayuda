#ifndef IMGAOS_HPP
#define IMGAOS_HPP

#include "common/progargs.hpp"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace imgaos {

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
  bool leerNumeroMagico(std::ifstream& file);
  bool leerHeader(std::ifstream& file, Image& img);
  bool leerPixels(std::ifstream& file, Image& img);
  bool buildColorTable(const Image& img, std::vector<Color>& color_table, std::map<Color, uint32_t>& color_map);
  bool determineIndexSize(const std::map<Color, uint32_t>& color_map, int& index_size);
  bool generatePixelIndices(const Image& img, const std::map<Color, uint32_t>& color_map, std::vector<uint32_t>& pixel_indices);
  void writeHeader(std::ofstream& file, const Image& picture, int color_table_size);
  void writeColors(std::ofstream& file, int max_color, const std::vector<Color>& color_table);
  void writePixelIndicesToFile(std::ofstream& file, const std::vector<uint32_t>& pixel_indices, int index_size);
  void loadImageAndPrepareTable(const std::string& input_file, Image& img, std::vector<Color>& color_table, std::map<Color, uint32_t>& color_map);
  void compress(progargsCommon::parameters_files & params);

  //RESIZE

  constexpr int MAX_COLOR_VALUE = 255;

    struct Red {
        uint16_t value;
        explicit Red(uint16_t value) : value(value) {}
    };

    struct Green {
        uint16_t value;
        explicit Green(uint16_t value) : value(value) {}
    };

    struct Blue {
        uint16_t value;
        explicit Blue(uint16_t value) : value(value) {}
    };

    struct Pixel {
        Red r;
        Green g;
        Blue b;
        Pixel() : r(0), g(0), b(0) {} // Default constructor
        Pixel(Red red, Green green, Blue blue) : r(red), g(green), b(blue) {}
    };

    struct ImageAOS {
        int width;
        int height;
        std::vector<Pixel> pixels;

        ImageAOS(int img_width, int img_height)
                : width(img_width), height(img_height),
                  pixels(static_cast<std::size_t>(img_width) * static_cast<std::size_t>(img_height)) {}

        Pixel &getPixel(int x_cord, int y_cord);
        [[nodiscard]] const Pixel &getPixel(int x_cord, int y_cord) const;
    };

    struct NeighboringPixels {
        uint16_t top_left;
        uint16_t top_right;
        uint16_t bottom_left;
        uint16_t bottom_right;
    };

    struct PixelCoordinates {
        int x_low;
        int y_low;
        int x_high;
        int y_high;
        char channel;
    };

    std::pair<double, double> calculateSourceCoordinates(int x_cord, int y_cord, const ImageAOS &original_image, const ImageAOS &resized_image);
    double interpolate(double start, double end, double factor);
    uint16_t bilinearInterpolate(const NeighboringPixels &neighbors, double x_factor, double y_factor);
    NeighboringPixels getNeighboringPixels(const ImageAOS &image, const PixelCoordinates &coords);

    ImageAOS readPPM(const std::string &filename);
    void writePPM(const std::string &filename, const ImageAOS &image);
    void resizeAOS(const ImageAOS &original_image, ImageAOS &resized_image);

    // Single-call resize function that uses the simplified read, process, and write approach
    void resize(const std::vector<std::string> &args);




  //CutFreq
  /*

  constexpr int SALTAR_LINEA = 255;

  struct pixel {
    uint8_t r;
    uint8_t g;
    uint8_t b;
    bool operator<(const pixel& other) const {
      return std::tie(r, g, b) < std::tie(other.r, other.g, other.b);
    }
    bool operator==(const pixel& other) const {
      return r == other.r && g == other.g && b == other.b;
    }
  };

  struct pares {
    pixel quitar;
    pixel poner;
  };

  struct pixel_soa {
    std::vector<uint8_t> r;
    std::vector<uint8_t> g;
    std::vector<uint8_t> b;
  };

  struct images {
    std::string entrada;
    std::string salida;
  };

  int euclidean_distance_squared(const pixel& pixel1, const pixel& pixel2);
  pares min_distance_vector(pixel const & pix1, std::vector<pixel> const & pix2);
  std::vector<pixel> leer_imagen(const std::string& entrada);
  std::map<pixel, int> recibir_mapa_frecuencias(std::string const & entrada);
  void sort_pixmap(std::map<pixel, int>& pixmap);
  std::vector<pixel> extract_sort(std::map<pixel, int> const& sorted_pixmap, std::vector<pixel> spixvex);
  void elegir_eliminados(size_t max_colors, const std::map<pixel, int>& pixmap, std::vector<pixel>& eliminar,
                         std::vector<pixel>& no_eliminar);
  void calculo_nuevos_colores(std::vector<pixel>& eliminar, const std::vector<pixel>& spixvex, std::vector<pares>& pareja);
  void cambiar_pixeles_imagen(std::vector<pixel>& pixvec, const std::vector<pares>& pareja,
                            std::map<pixel, pixel>& replacements);
  void guardar_imagen(const std::string& entrada, const std::string& salida, const std::vector<pixel>& pixvec);
  void cutFreq(const progargs::Parameters& params, size_t max_colors);
    */
}

#endif