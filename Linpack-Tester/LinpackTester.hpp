#ifndef LINPACKTEST_LINPACKTESTER_HPP
#define LINPACKTEST_LINPACKTESTER_HPP

#include "Matrix.hpp"

class LinpackTester {
public:
  explicit LinpackTester() = default;

  virtual ~LinpackTester() = default;

  virtual double invokeTest() = 0;
};

#endif //LINPACKTEST_LINPACKTESTER_HPP
