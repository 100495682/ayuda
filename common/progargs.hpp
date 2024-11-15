#ifndef COMMON_PROGARGS_HPP
#define COMMON_PROGARGS_HPP

#include <vector>
#include <string>

namespace progargsCommon {

  constexpr int SEIS = 6;
 constexpr int CINCO = 5;
  constexpr int MAX_COLOR_VALUE_16BIT = 65535;

  struct parameters_files {
    std::string input_file;
    std::string output_file;
  };

  bool comprobarArgc(int &argc);
  bool comprobarArgumentos(std::vector<std::string> const &args);
  bool pert(const std::string& operation);
}

#endif // COMMON_PROGARGS_HPP