#include "imgaos.hpp"

#include "common/binaryo.hpp"
#include "common/progargs.hpp"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cmath>

namespace imgaos {
  ///COMPRESS

  bool Color::operator<(const Color& other) const {
    if (r != other.r) { return r < other.r; }
    if (g != other.g) { return g < other.g; }
    return b < other.b;
}

  bool abrirFile(const std::string& filename, std::ifstream& file) {
      file.open(filename, std::ios::binary);
      if (!file.is_open()) {
          std::cerr << "Error al abrir el archivo de entrada: " << filename << "\n";
          return false;
      }
      return true;
  }

  bool leerNumeroMagico(std::ifstream& file) {
      std::string numero_magico;
      file >> numero_magico;
      if (numero_magico != "P6") {
          std::cerr << "Formato de archivo de entrada erróneo. Se esperaba 'P6', pero se encontró '" << numero_magico << "'.\n";
          return false;
      }
      return true;
  }

  bool leerHeader(std::ifstream& file, Image& img) {
      file >> img.width >> img.height >> img.max_color;
      file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      return !file.fail();
  }

  bool leerPixels(std::ifstream& file, Image& img) {
      const std::size_t num_pixels = static_cast<std::size_t>(img.width) * static_cast<std::size_t>(img.height);
      img.pixels.resize(num_pixels);

      if (img.max_color <= LONGITUD_CABECERA) {
          // Leer colores de 1 byte
          for (auto& color : img.pixels) {
              color.r = binaryo::read_binary<uint8_t>(file);
              color.g = binaryo::read_binary<uint8_t>(file);
              color.b = binaryo::read_binary<uint8_t>(file);
          }
      } else {
          // Leer colores de 2 bytes
          for (auto& color : img.pixels) {
              color.r = binaryo::read_binary<uint16_t>(file);
              color.g = binaryo::read_binary<uint16_t>(file);
              color.b = binaryo::read_binary<uint16_t>(file);
          }
      }

      return !file.fail();
  }

  bool buildColorTable(const Image& img, std::vector<Color>& color_table, std::map<Color, uint32_t>& color_map) {
      for (const auto& color : img.pixels) {
          if (color_map.find(color) == color_map.end()) {
              color_map[color] = static_cast<uint32_t>(color_map.size());
              color_table.push_back(color);
          }
      }
      return true;
  }

  bool determineIndexSize(const std::map<Color, uint32_t>& color_map, int& index_size) {
      const size_t num_colors = color_map.size();
      if (num_colors <= MAX_COLORS_8BIT) {
          index_size = 1;
      }
      else if (num_colors <= MAX_COLORS_16BIT) {
          index_size = 2;
      }
      else if (num_colors <= MAX_COLORS_32BIT) {
          index_size = 4;
      }
      else {
          std::cerr << "Número de colores excede el soporte máximo.\n";
          return false;
      }
      return true;
  }

  bool generatePixelIndices(const Image& img, const std::map<Color, uint32_t>& color_map, std::vector<uint32_t>& pixel_indices) {
      pixel_indices.reserve(img.pixels.size());
      for (const auto& color : img.pixels) {
          auto iterator = color_map.find(color);
          if (iterator != color_map.end()) {
              pixel_indices.push_back(iterator->second);
          } else {
              std::cerr << "Error: Píxel no encontrado en el mapa de colores.\n";
              return false;
          }
      }
      return true;
  }

  void writeHeader(std::ofstream& file, const Image& picture, int color_table_size) {
      file << "C6 " << picture.width << " " << picture.height << " " << picture.max_color << " " << color_table_size << "\n";
  }

  void writeColors(std::ofstream& file, int max_color, const std::vector<Color>& color_table) {
      constexpr int limit1 = 255;

      for (const auto& color : color_table) {
          if (max_color <= limit1) {
              binaryo::write_binary<uint8_t>(file, static_cast<uint8_t>(color.r));
              binaryo::write_binary<uint8_t>(file, static_cast<uint8_t>(color.g));
              binaryo::write_binary<uint8_t>(file, static_cast<uint8_t>(color.b));
          } else {
              binaryo::write_binary<uint16_t>(file, color.r);
              binaryo::write_binary<uint16_t>(file, color.g);
              binaryo::write_binary<uint16_t>(file, color.b);
          }
      }
  }

  void writePixelIndicesToFile(std::ofstream& file, const std::vector<uint32_t>& pixel_indices, int index_size) {
      for (const auto& index : pixel_indices) {
          if (index_size == 1) {
              binaryo::write_binary<uint8_t>(file, static_cast<uint8_t>(index));
          } else if (index_size == 2) {
              binaryo::write_binary<uint16_t>(file, static_cast<uint16_t>(index));
          } else {
              binaryo::write_binary<uint32_t>(file, index);
          }
      }
  }

    void loadImageAndPrepareTable(const std::string& input_file, Image& img, std::vector<Color>& color_table, std::map<Color, uint32_t>& color_map) {
  std::ifstream file;
  if (!abrirFile(input_file, file)) {
    std::cerr << "Error al abrir el archivo de entrada.\n";
    return;
  }
  if (!leerNumeroMagico(file) || !leerHeader(file, img) || !leerPixels(file, img)) {
    file.close();
    std::cerr << "Error al leer el archivo de entrada.\n";
    return;
  }
  file.close();

  for (const auto& color : img.pixels) {
    if (color_map.find(color) == color_map.end()) {
      color_map[color] = static_cast<uint32_t>(color_map.size());
      color_table.push_back(color);
    }
  }
}

    void compress(progargsCommon::parameters_files & params) {
      Image img;
      std::vector<Color> color_table;
      std::map<Color, uint32_t> color_map;

      loadImageAndPrepareTable(params.input_file, img, color_table, color_map);

      int index_size = 0;
      if (!determineIndexSize(color_map, index_size)) {
        std::cerr << "Error: No se pudo determinar el tamaño del índice.\n";
        return;
      }

      std::vector<uint32_t> pixel_indices;
      if (!generatePixelIndices(img, color_map, pixel_indices)) {
        std::cerr << "Error: No se pudo generar los índices de píxeles.\n";
        return;
      }

      std::ofstream output_file(params.output_file, std::ios::binary);
      if (!output_file) {
        std::cerr << "Error al abrir el archivo de salida.\n";
        return;
      }

      writeHeader(output_file, img, static_cast<int>(color_table.size()));
      writeColors(output_file, img.max_color, color_table);
      writePixelIndicesToFile(output_file, pixel_indices, index_size);
      output_file.close();
      std::cout << "Compresión completada exitosamente.\n";
}

  //RESIZE

  Pixel &ImageAOS::getPixel(int x_cord, int y_cord) {
    return pixels[(static_cast<std::size_t>(y_cord) * static_cast<std::size_t>(width)) + static_cast<std::size_t>(x_cord)];
  }

  const Pixel &ImageAOS::getPixel(int x_cord, int y_cord) const {
    return pixels[(static_cast<std::size_t>(y_cord) * static_cast<std::size_t>(width)) + static_cast<std::size_t>(x_cord)];
  }

  ImageAOS readPPM(const std::string &filename) {
    std::ifstream file(filename, std::ios::binary);

    std::string format;
    int width = 0;
    int height = 0;
    int max_color = 0;
    file >> format;
    file >> width >> height >> max_color;
    file.ignore(1);

    ImageAOS image(width, height);
    for (std::size_t i = 0; i < static_cast<std::size_t>(width) * static_cast<std::size_t>(height); ++i) {
      const auto red = binaryo::read_binary<uint8_t>(file);
      const auto green = binaryo::read_binary<uint8_t>(file);
      const auto blue = binaryo::read_binary<uint8_t>(file);
      image.pixels[i] = Pixel(Red(red), Green(green), Blue(blue));
    }

    return image;
  }

  void writePPM(const std::string &filename, const ImageAOS &image) {
    std::ofstream file(filename, std::ios::binary);

    file << "P6\n" << image.width << " " << image.height << "\n" << MAX_COLOR_VALUE << "\n";
    for (const auto &pixel : image.pixels) {
      binaryo::write_binary(file, static_cast<uint8_t>(pixel.r.value));
      binaryo::write_binary(file, static_cast<uint8_t>(pixel.g.value));
      binaryo::write_binary(file, static_cast<uint8_t>(pixel.b.value));
    }
  }

  std::pair<double, double> calculateSourceCoordinates(int x_cord, int y_cord, const ImageAOS &original_image, const ImageAOS &resized_image) {
    const double src_x = x_cord * ((original_image.width - 1) / static_cast<double>(resized_image.width - 1));
    const double src_y = y_cord * ((original_image.height - 1) / static_cast<double>(resized_image.height - 1));
    return {src_x, src_y};
  }

  double interpolate(double start, double end, double factor) {
    return start + (factor * (end - start));
  }

  uint16_t bilinearInterpolate(const NeighboringPixels &neighbors, double x_factor, double y_factor) {
    const double top = interpolate(neighbors.top_left, neighbors.top_right, x_factor);
    const double bottom = interpolate(neighbors.bottom_left, neighbors.bottom_right, x_factor);
    return static_cast<uint16_t>(std::round(interpolate(top, bottom, y_factor)));
  }

  NeighboringPixels getNeighboringPixels(const ImageAOS &image, const PixelCoordinates &coords) {
    NeighboringPixels neighbors = {};

    switch (coords.channel) {
      case 'r':
        neighbors.top_left = image.getPixel(coords.x_low, coords.y_low).r.value;
      neighbors.top_right = image.getPixel(coords.x_high, coords.y_low).r.value;
      neighbors.bottom_left = image.getPixel(coords.x_low, coords.y_high).r.value;
      neighbors.bottom_right = image.getPixel(coords.x_high, coords.y_high).r.value;
      break;
      case 'g':
        neighbors.top_left = image.getPixel(coords.x_low, coords.y_low).g.value;
      neighbors.top_right = image.getPixel(coords.x_high, coords.y_low).g.value;
      neighbors.bottom_left = image.getPixel(coords.x_low, coords.y_high).g.value;
      neighbors.bottom_right = image.getPixel(coords.x_high, coords.y_high).g.value;
      break;
      case 'b':
        neighbors.top_left = image.getPixel(coords.x_low, coords.y_low).b.value;
      neighbors.top_right = image.getPixel(coords.x_high, coords.y_low).b.value;
      neighbors.bottom_left = image.getPixel(coords.x_low, coords.y_high).b.value;
      neighbors.bottom_right = image.getPixel(coords.x_high, coords.y_high).b.value;
      break;
      default:
        throw std::invalid_argument("Unexpected channel value in getNeighboringPixels.");
    }

    return neighbors;
  }

  void resizeAOS(const ImageAOS &original_image, ImageAOS &resized_image) {
    resized_image = ImageAOS(resized_image.width, resized_image.height);
    for (int y_cord = 0; y_cord < resized_image.height; ++y_cord) {
      for (int x_cord = 0; x_cord < resized_image.width; ++x_cord) {
        auto [src_x, src_y] = calculateSourceCoordinates(x_cord, y_cord, original_image, resized_image);

        const int x_low = static_cast<int>(std::floor(src_x));
        const int y_low = static_cast<int>(std::floor(src_y));
        const int x_high = std::min(x_low + 1, original_image.width - 1);
        const int y_high = std::min(y_low + 1, original_image.height - 1);

        const double x_factor = src_x - x_low;
        const double y_factor = src_y - y_low;

        const PixelCoordinates red_coords = {.x_low = x_low, .y_low = y_low, .x_high = x_high, .y_high = y_high, .channel = 'r'};
        const PixelCoordinates green_coords = {.x_low = x_low, .y_low = y_low, .x_high = x_high, .y_high = y_high, .channel = 'g'};
        const PixelCoordinates blue_coords = {.x_low = x_low, .y_low = y_low, .x_high = x_high, .y_high = y_high, .channel = 'b'};

        const NeighboringPixels red_neighbors = getNeighboringPixels(original_image, red_coords);
        const NeighboringPixels green_neighbors = getNeighboringPixels(original_image, green_coords);
        const NeighboringPixels blue_neighbors = getNeighboringPixels(original_image, blue_coords);

        resized_image.getPixel(x_cord, y_cord) = Pixel(
                Red(bilinearInterpolate(red_neighbors, x_factor, y_factor)),
                Green(bilinearInterpolate(green_neighbors, x_factor, y_factor)),
                Blue(bilinearInterpolate(blue_neighbors, x_factor, y_factor))
        );
      }
    }
  }

  void resize(const std::vector<std::string> &args) {
    const ImageAOS original_image = readPPM(args[1]);
    ImageAOS resized_image(std::stoi(args[4]), std::stoi(args[progargsCommon::CINCO]));
    resizeAOS(original_image, resized_image);
    writePPM(args[2], resized_image);
  }

  //CutFreq
/*

  int euclidean_distance_squared(const pixel& pixel1, const pixel& pixel2) {
    return  ((pixel1.r - pixel2.r) * (pixel1.r - pixel2.r))
          + ((pixel1.g - pixel2.g) * (pixel1.g - pixel2.g))
          + ((pixel1.b - pixel2.b) * (pixel1.b - pixel2.b));
  }

  pares min_distance_vector(pixel const & pix1, std::vector<pixel> const & pix2) {
    int min_dist = std::numeric_limits<int>::max();
    pixel min_p  = pix2[0];
    for (auto const & pix : pix2) {
      int const dist = euclidean_distance_squared(pix1, pix);
      if (dist == 1) { return {.quitar = pix1, .poner = min_p}; }
      if (dist < min_dist) {
        min_dist = dist;
        min_p    = pix;
      }
    }
    return {.quitar = pix1, .poner = min_p};
  }

  auto leer_imagen(const std::string& entrada) {
    std::ifstream file(entrada, std::ios::binary);
    if (!file.is_open()) { std::cerr << "Error al abrir la imagen.\n"; }
    std::string formato;
    file >> formato;
    unsigned int ancho = 0;
    unsigned int alto  = 0;
    unsigned int luz   = 0;
    file >> ancho >> alto >> luz;
    file.ignore(SALTAR_LINEA, '\n');
    std::vector<pixel> pixvec(static_cast<std::vector<int>::size_type>(ancho * alto));
    for (auto & pix : pixvec) {
      pix.r = file.get();
      pix.g = file.get();
      pix.b = file.get();
    }
    return pixvec;
  }

  auto recibir_mapa_frecuencias(std::string const & entrada) {
    std::ifstream file(entrada, std::ios::binary);
    if (!file.is_open()) { std::cerr << "Error al abrir la imagen.\n"; }
    std::string formato;
    file >> formato;
    unsigned int ancho = 0;
    unsigned int alto  = 0;
    unsigned int luz   = 0;
    file >> ancho >> alto >> luz;
    file.ignore(SALTAR_LINEA, '\n');
    std::map<pixel, int> pixmap;
    std::vector<pixel> const pixvec(static_cast<std::vector<int>::size_type>(ancho * alto));
    for (const auto & pix : pixvec) {

      pixmap[pix]++;

    }
    return pixmap;
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
  void sort_pixmap(std::map<pixel, int> pixmap) {
    std::vector<std::pair<pixel, int>> sorted_pixmap(pixmap.begin(), pixmap.end());
    std::ranges::sort(sorted_pixmap.begin(), sorted_pixmap.end(), [](auto const & aaa, auto const & bbb) {
      return aaa.second < bbb.second;
    });
  }

  std::vector<pixel> extract_sort(std::map<pixel, int> const& sorted_pixmap, std::vector<pixel> spixvex) {
    for (auto const& pix : sorted_pixmap) {
      spixvex.push_back(pix.first);
    };
    return spixvex;
  }

  void elegir_eliminados(long unsigned max_colors, const std::map<pixel, int>& pixmap, std::vector<pixel> eliminar,
                         std::vector<pixel> no_eliminar) {
    for (auto const & ver : pixmap) {
      if (eliminar.size() < max_colors) {
        eliminar.push_back(ver.first);  //
      } else {
        no_eliminar.push_back(ver.first);
      }
    }
  }

  void calculo_nuevos_colores(std::vector<pixel>& eliminar, const std::vector<pixel>& spixvex, std::vector<pares> & pareja) {
    pareja.reserve(eliminar.size());
    for (auto const & elim_px : eliminar) { pareja.push_back(min_distance_vector(elim_px, spixvex)); }
  }

  void cambiar_pixeles_imagen(std::vector<pixel> pixvec, const std::vector<pares>& pareja,
                            std::map<pixel, pixel> replacements) {
    for (auto const & [quitar, poner] : pareja) { replacements[quitar] = poner; }
    auto px_end = pixvec.end();  // Iterador al final de pixvec
    for (auto px_iter = pixvec.begin(); px_iter != px_end; ++px_iter) {
      auto ita = replacements.find(*px_iter);  // Buscar el pixel actual en el mapa de reemplazos
      if (ita != replacements.end()) {
        *px_iter = ita->second;  // Si se encuentra, reemplazar el pixel
      }
    }
  }

  void guardar_imagen(const images& archivos, const std::vector<pixel>& pixvec) {
    std::ofstream imgout(archivos.entrada, std::ios::binary);
    std::string const formato = "P6";
    std::vector<int> const param = parametros(archivos.entrada);
    int const ancho   = param[0];
    int const alto     = param[1];
    int const luz       = param[2];
    imgout << formato << "\n" << ancho << " " << alto << "\n" << luz << "\n";

    for (auto const & pix : pixvec) {
      imgout.put(static_cast<char>(pix.r));
      imgout.put(static_cast<char>(pix.g));
      imgout.put(static_cast<char>(pix.b));
    }
  }

  void cutFreq(const progargsCompress::parameters_files &params, long unsigned max_colors) {

    // CREAMOS VECTOR DE PIXELES Y MAPA DE FRECUENCIAS
    images const archivos = {.entrada=params.input_file,.salida=params.output_file};
    std::vector<pixel> const pixvec = leer_imagen(params.input_file); // Entrada leída y devuelve vector pixvec
    std::map<pixel, int> const pixmap =recibir_mapa_frecuencias(params.input_file); // Devuelve el mapa de frecuencia de vectores
    std::map<pixel,int>const& sorted_pixmap = pixmap;
    sort_pixmap(sorted_pixmap); // Ordenamos el mapa de recuencias

    // CREAMOS VECTORES PARA OPERAR COLORES A ELIMINAR Y A GUARDAR
    std::vector<pixel> eliminar;
    std::vector<pixel> no_eliminar;
    eliminar.reserve(max_colors);
    no_eliminar.reserve(sorted_pixmap.size());
    std::vector<pixel> const spixvex(pixmap.size());
    extract_sort(sorted_pixmap, spixvex); // Este vector tiene todos los colores únicos
    elegir_eliminados(max_colors, pixmap, eliminar, no_eliminar);

    // CALCULAMOS LAS PAREJAS DE COLORES A ELIMINAR Y REEMPLAZAMOS
    std::vector<pares> pareja;
    calculo_nuevos_colores(eliminar, spixvex, pareja);
    std::map<pixel, pixel> const replacements;
    cambiar_pixeles_imagen(pixvec, pareja, replacements);

    // OUTPUT IMAGEN
    guardar_imagen(archivos, pixvec);
  }
  */
}