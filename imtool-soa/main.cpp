#include "common/decompress.hpp"
#include "common/progargs.hpp"
#include "imgaos/imgaos.hpp"
#include "imgsoa/imgsoa.hpp"

#include <cstring>
#include <gsl/gsl>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char* argv[]) {
  if (!progargsCommon::comprobarArgc(argc)) {
    return -1;
  }

  gsl::span const args_view{argv, gsl::narrow<std::size_t>(argc)};
  std::vector<std::string> const arguments(args_view.begin(), args_view.end());

  progargsCommon::parameters_files params;
  params.input_file = arguments[1];
  params.output_file = arguments[2];

  if (!progargsCommon::comprobarArgumentos(arguments)) {
    std::cerr << "Invalid arguments.\n";
    return 1;
  }
  if (strcmp(arguments[3].c_str(), "resize") == 0) {
    imgsoa::resize(arguments);
  }else if (arguments[3] == "compress") {
    imgsoa::compress(params);
  }else{
    std::cerr << "Unsupported function: " << arguments[3] << "\n";
    return 1;
  }

  return 0;
}