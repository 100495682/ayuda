#include "imgaos/imgaos.hpp"

#include <fstream>
#include <gtest/gtest.h>

/*
  TEST(ImgaosTest, AbrirFile) {
    std::ifstream file;
    const bool result = imgaos::abrirFile("test_image.ppm", file);
    EXPECT_TRUE(result);
    EXPECT_TRUE(file.is_open());
  }

  TEST(ImgaosTest, LeerNumeroMagico) {
    std::ifstream file("test_image.ppm");
    ASSERT_TRUE(file.is_open());
    const bool result = imgaos::leerNumeroMagico(file);
    EXPECT_TRUE(result);
  }

  TEST(ImgaosTest, LeerHeader) {
    std::ifstream file("test_image.ppm");
    ASSERT_TRUE(file.is_open());
    imgaos::Image img;
    const bool result = imgaos::leerHeader(file, img);
    EXPECT_TRUE(result);
    EXPECT_GT(img.width, 0);
    EXPECT_GT(img.height, 0);
  }

  TEST(ImgaosTest, LeerPixels) {
    std::ifstream file("test_image.ppm");
    ASSERT_TRUE(file.is_open());
    imgaos::Image img;
    ASSERT_TRUE(imgaos::leerHeader(file, img));
    const bool result = imgaos::leerPixels(file, img);
    EXPECT_TRUE(result);
    EXPECT_EQ(img.pixels.size(), img.width * img.height);
  }

  TEST(ImgaosTest, BuildColorTable) {
    imgaos::Image img;
    img.width = 2;
    img.height = 2;
    img.pixels = {{255, 0, 0}, {0, 255, 0}, {0, 0, 255}, {255, 255, 255}};
    std::vector<imgaos::Color> color_table;
    std::map<imgaos::Color, uint32_t> color_map;
    const bool result = imgaos::buildColorTable(img, color_table, color_map);
    EXPECT_TRUE(result);
    EXPECT_EQ(color_table.size(), 4);
  }

  TEST(ImgaosTest, DetermineIndexSize) {
    std::map<imgaos::Color, uint32_t> color_map = {{{255, 0, 0}, 0}, {{0, 255, 0}, 1}, {{0, 0, 255}, 2}, {{255, 255, 255}, 3}};
    int index_size;
    const bool result = imgaos::determineIndexSize(color_map, index_size);
    EXPECT_TRUE(result);
    EXPECT_EQ(index_size, 2);
  }

  TEST(ImgaosTest, GeneratePixelIndices) {
    imgaos::Image img;
    img.width = 2;
    img.height = 2;
    img.pixels = {{255, 0, 0}, {0, 255, 0}, {0, 0, 255}, {255, 255, 255}};
    std::map<imgaos::Color, uint32_t> color_map = {{{255, 0, 0}, 0}, {{0, 255, 0}, 1}, {{0, 0, 255}, 2}, {{255, 255, 255}, 3}};
    std::vector<uint32_t> pixel_indices;
    const bool result = imgaos::generatePixelIndices(img, color_map, pixel_indices);
    EXPECT_TRUE(result);
    EXPECT_EQ(pixel_indices.size(), 4);
  }

  TEST(ImgaosTest, WriteHeader) {
    std::ofstream file("output.img");
    ASSERT_TRUE(file.is_open());
    imgaos::Image img;
    img.width = 2;
    img.height = 2;
    img.max_color = 255;
    imgaos::writeHeader(file, img, 4);
    file.close();
    std::ifstream input_file("output.img");
    EXPECT_TRUE(input_file.is_open());
  }

  TEST(ImgaosTest, WriteColors) {
    std::ofstream file("output.img");
    ASSERT_TRUE(file.is_open());
    std::vector<imgaos::Color> color_table = {{255, 0, 0}, {0, 255, 0}, {0, 0, 255}, {255, 255, 255}};
    imgaos::writeColors(file, 255, color_table);
    file.close();
    std::ifstream input_file("output.img");
    EXPECT_TRUE(input_file.is_open());
  }

  TEST(ImgaosTest, WritePixelIndicesToFile) {
    std::ofstream file("output.img");
    ASSERT_TRUE(file.is_open());
    std::vector<uint32_t> pixel_indices = {0, 1, 2, 3};
    imgaos::writePixelIndicesToFile(file, pixel_indices, 2);
    file.close();
    std::ifstream input_file("output.img");
    EXPECT_TRUE(input_file.is_open());
  }

  TEST(ImgaosTest, LoadImageAndPrepareTable) {
    imgaos::Image img;
    std::vector<imgaos::Color> color_table;
    std::map<imgaos::Color, uint32_t> color_map;
    imgaos::loadImageAndPrepareTable("test_image.ppm", img, color_table, color_map);
    EXPECT_GT(img.width, 0);
    EXPECT_GT(img.height, 0);
    EXPECT_GT(color_table.size(), 0);
    EXPECT_GT(color_map.size(), 0);
  }

  TEST(ImgaosTest, Compress) {
    // Prepare input parameters
    progargs::parameters_files params;
    params.input_file = "test_image.ppm";
    params.output_file = "compressed_image.img";

    // Call the compress function
    imgaos::compress(params);

    // Verify the output file is created
    std::ifstream output_file(params.output_file);
    EXPECT_TRUE(output_file.is_open());

    // Additional checks can be added here to verify the contents of the output file
  }
*/
