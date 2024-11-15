// imgsoa.hpp
#ifndef IMGSOA_HPP
#define IMGSOA_HPP

#include <cstdint>
#include <vector>
#include <map>
#include <string>
#include "common/progargs.hpp"
#include "common/binaryo.hpp"

namespace imgsoa {

  constexpr int MAX_COLOR_1 = 255;
  constexpr int MAX_COLOR_2 = 65535;

  struct dimensions {
    int height;
    int width;
  };

  class Image {
    public:

        dimensions dim;
        int max_color;

        std::vector<uint16_t> r;
        std::vector<uint16_t> g;
        std::vector<uint16_t> b;

        uint32_t number_colors = 0;

        Image(dimensions dim_init, int max_color_init);
  };

  std::map<std::tuple<uint16_t, uint16_t, uint16_t>, uint32_t> leerColores(Image& image);
  void escribirTabla(Image &img,
                     std::map<std::tuple<uint16_t, uint16_t, uint16_t>, uint32_t> &listaColor,
                     std::string &archivo);
  void escribirColores(Image &img,
                     std::map<std::tuple<uint16_t, uint16_t, uint16_t>, uint32_t> &listaColor,
                     std::string &archivo);
  void escribirIndice(uint32_t indice, std::ofstream& file, uint32_t number_colors);
  void escribirColor(const std::tuple<uint16_t, uint16_t, uint16_t>& color, std::ofstream& file, int max_color);
  void compress(progargsCommon::parameters_files & params);
  //RESIZE

  constexpr int MAX_COLOR_VALUE = 255;

// Enum for color channels with specified underlying type
    enum class Channel : std::uint8_t {
        Red,
        Green,
        Blue
    };

// Structures

    struct ImageSOA {
        std::size_t width;
        std::size_t height;
        std::vector<uint8_t> r;
        std::vector<uint8_t> g;
        std::vector<uint8_t> b;

        explicit ImageSOA(std::size_t img_width = 0, std::size_t img_height = 0)
                : width(img_width),
                  height(img_height),
                  r(img_width * img_height),
                  g(img_width * img_height),
                  b(img_width * img_height) {}
    };

    struct Coordinates {
        std::size_t x_low = 0;
        std::size_t y_low = 0;
        std::size_t x_high = 0;
        std::size_t y_high = 0;
    };

    struct PixelCoordinates {
        Coordinates coords{};
        Channel channel = Channel::Red;

        PixelCoordinates() = default;
        PixelCoordinates(const Coordinates& coordinates, Channel channel)
                : coords(coordinates), channel(channel) {}
    };

    struct NeighboringPixels {
        uint8_t top_left = 0;
        uint8_t top_right = 0;
        uint8_t bottom_left = 0;
        uint8_t bottom_right = 0;
    };

// Function declarations

    std::pair<double, double> calculateSourceCoordinates(
            std::size_t x_cord,
            std::size_t y_cord,
            const ImageSOA& original_image,
            const ImageSOA& resized_image);

    uint8_t bilinearInterpolate(
            const NeighboringPixels& neighbors,
            double x_factor,
            double y_factor);

    NeighboringPixels getNeighboringPixels(
            const ImageSOA& image,
            const PixelCoordinates& pixel_coords);

// Function declarations for SOA operations
    ImageSOA readPPM(const std::string& filename);
    void writePPM(const std::string& filename, const ImageSOA& image);
    inline double interpolate(double start, double end, double factor);

// Existing resize function that resizes an image given input and output ImageSOA objects
    void resize(const ImageSOA& original_image, ImageSOA& resized_image);

// Updated resizeSOA function that handles file I/O and uses progargs::Parameters
  void resize(const std::vector<std::string> &args);

  //CUTFREQ
  /*

  struct pixel_soa {
    std::vector<uint8_t> r;
    std::vector<uint8_t> g;
    std::vector<uint8_t> b;
  };
  struct images {
    std::string input_file;
    std::string output_file;
  };
  int euclidean_distance_squared(const pixel_soa& pix1, const pixel_soa& pix2, size_t index);
  pixel_soa read_image(const std::string& filename);
  std::tuple<uint8_t, uint8_t, uint8_t> min_distance_vector(const pixel_soa& pix1, const pixel_soa& pix2);
  std::map<std::tuple<uint8_t, uint8_t, uint8_t>, int> read_frequency_map(const std::string& filename);
  std::vector<std::tuple<uint8_t, uint8_t, uint8_t>> extract_sorted(const std::map<std::tuple<uint8_t, uint8_t, uint8_t>, int>& freq_map);
  void choose_eliminated(size_t max_colors, const std::map<std::tuple<uint8_t, uint8_t, uint8_t>, int>& freq_map,
                         std::vector<std::tuple<uint8_t, uint8_t, uint8_t>>& to_eliminate,
                         std::vector<std::tuple<uint8_t, uint8_t, uint8_t>>& to_keep);
  std::vector<std::pair<std::tuple<uint8_t, uint8_t, uint8_t>, std::tuple<uint8_t, uint8_t, uint8_t>>>
  calculate_new_colors(const std::vector<std::tuple<uint8_t, uint8_t, uint8_t>>& to_eliminate,
                       const pixel_soa& pixsoa);

  void replace_pixels(pixel_soa& pixsoa, const std::vector<std::pair<std::tuple<uint8_t, uint8_t, uint8_t>, std::tuple<uint8_t, uint8_t, uint8_t>>>& pairs);
  auto parametros(const std::string&  entrada);
  void save_image(const std::string& filename, const pixel_soa& pixsoa, const std::vector<int>& parametros);

  void cutFreq(const std::vector<std::string> &args);
  */

  /*
  //CutFreq

  struct parameters_files {
    std::string input_file;
    std::string output_file;
    int max_colors;
  };

  //CUTFREQ

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

  struct images {
    std::string entrada;
    std::string salida;
  };

  struct pixel_soa {
    std::vector<uint8_t> r;
    std::vector<uint8_t> g;
    std::vector<uint8_t> b;
  };

  std::vector<std::tuple<uint8_t, uint8_t, uint8_t>> extract_sorted(const std::map<std::tuple<uint8_t, uint8_t, uint8_t>, int>& freq_map);
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
  std::map<std::tuple<uint8_t, uint8_t, uint8_t>, int> read_frequency_map(const std::string& filename);
  pixel_soa read_image(const std::string& filename);
  std::tuple<uint8_t, uint8_t, uint8_t> min_distance_vector(const pixel_soa& pix1, const pixel_soa& pix2);
  int euclidean_distance_squared(const pixel_soa& pix1, const pixel_soa& pix2, size_t index);
  void choose_eliminated(size_t max_colors, const std::map<std::tuple<uint8_t, uint8_t, uint8_t>, int>& freq_map,
                       std::vector<std::tuple<uint8_t, uint8_t, uint8_t>>& to_eliminate,
                       std::vector<std::tuple<uint8_t, uint8_t, uint8_t>>& to_keep);
  std::vector<std::pair<std::tuple<uint8_t, uint8_t, uint8_t>, std::tuple<uint8_t, uint8_t, uint8_t>>>
calculate_new_colors(const std::vector<std::tuple<uint8_t, uint8_t, uint8_t>>& to_eliminate,
                     const pixel_soa& pixsoa);
  void replace_pixels(pixel_soa& pixsoa, const std::vector<std::pair<std::tuple<uint8_t, uint8_t, uint8_t>, std::tuple<uint8_t, uint8_t, uint8_t>>>& pairs);
  auto parametros(const std::string&  entrada);
  void save_image(const std::string& filename, const pixel_soa& pixsoa, const std::vector<int>& parametros);
  */

}

#endif // IMGSOA_HPP
