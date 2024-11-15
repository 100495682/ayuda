//
// Created by alejo on 11/11/24.
//

#include <iostream>

#ifndef BINARYO_HPP
#define BINARYO_HPP

namespace binaryo{
  template <typename T>
T   read_binary(std::istream & input) {
    T value;
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
    input.read(reinterpret_cast<char*>(&value), sizeof(value));
    return value;
  }

  template <typename T>
  void write_binary(std::ostream & output, T const & value) {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
    output.write(reinterpret_cast<char const*>(&value), sizeof(value));
  }
}
#endif //BINARYO_HPP
