//
// Created by alejo on 10/11/24.
//

#include "common/info.hpp"

#include <fstream>
#include <iostream>
#include <string>

namespace common {
  void info(const progargsCommon::parameters_files& params) {
    std::ifstream imagen(params.input_file, std::ios::binary);
    if (imagen.fail()) {
      std::cerr << "Error opening file " << params.input_file << '\n';
      return;
    }
    std::string identification;
    int width     = 0;
    int height    = 0;
    int max_color = 0;
    imagen >> identification >> width >> height >> max_color;

    std::cout << "Input: "<< params.input_file <<'\n'<<"Output: "<<params.output_file<<'\n'<<"Operation: info"<<"\n";
    std::cout <<"Image size: "<<width<<"x"<<height<<"\n"<<"Max level: "<<max_color;
  }
}