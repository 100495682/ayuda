#include <cstring>
#include <iostream>
#include <string>
#include <vector>
#include "progargs.hpp"
#include <algorithm>
#include <ranges>


namespace progargsCommon {

    bool comprobarArgc(int &argc) {
        return argc >= 4;
    }

    bool comprobarArgumentos(std::vector<std::string> const &args) {
        bool const contained = pert(args[3]);
        if (contained) {
            if (strcmp(args[3].c_str(), "compress") == 0) {
                if (args.size() != 4) {
                    std::cerr << "Error: Wrong number of arguments; Needed 4 parameters\nParameters introduced: "
                              << args.size() << "\n";
                    return false;
                }
            }
            if (strcmp(args[3].c_str(), "resize") == 0) {
                if (args.size() != SEIS) {
                    std::cerr << "Error: Wrong number of arguments; Needed 6 parameters\nParameters introduced: "
                              << args.size() << "\n";
                    return false;
                }
            }
            if (strcmp(args[3].c_str(), "maxLevel") == 0) {
                if (args.size() != CINCO) {
                    std::cerr << "Error: Wrong number of arguments; Needed 5 parameters\nParameters introduced: "
                              << args.size() << "\n";
                    return false;
                }
            }
            if (strcmp(args[3].c_str(), "cutfreq") == 0) {
                if (args.size() != CINCO) {
                    std::cerr << "Error: Wrong number of arguments; Needed 5 parameters\nParameters introduced: "
                              << args.size() << "\n";
                    return false;
                }
            }
        } else {
            std::cerr << "Wrong operation introduced; " << args[3] << " is not a valid operation";
            return false;
        }
        return true;
    }


  bool pert(const std::string& operation) {
      std::vector<std::string> const operations = {"compress", "resize", "cutfreq", "maxLevel", "info"};
      return std::ranges::any_of(operations, [&operation](const std::string& element) {
          return strcmp(element.c_str(), operation.c_str()) == 0;
      });
    }

} // namespace progargsCommon