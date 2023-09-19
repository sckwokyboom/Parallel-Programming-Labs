#ifndef MAIN_CPP_MATRIXERROR_HPP
#define MAIN_CPP_MATRIXERROR_HPP

//#include <exception>
#include <stdexcept>
#include <string>

class matrix_error : public std::runtime_error {
public:
  explicit matrix_error(const std::string &msg) : std::runtime_error(msg) {}
};

#endif //MAIN_CPP_MATRIXERROR_HPP
