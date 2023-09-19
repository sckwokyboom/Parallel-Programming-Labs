#include <gtest/gtest.h>
#include "Matrix.hpp"

TEST(Matrix, CreateTest) {
  Matrix A(50, 10, 'n');
  Matrix B(50, 50, 't');
  Matrix C(50, 1, 'r');
  Matrix D(50, 2, 'r');
  Matrix E(50, 3, 'e');
  Matrix F(50, 100, 'h');
}

TEST(Matrix, AriphmeticTest) {
  Matrix A(50, 50, 'e');
  Matrix B(A);
  A += A;
  A -= B;
  Matrix C = A + A - A * A;
  C *= A;
  for (int i = 0; i < 50; ++i) {
    for (int j = 0; j < 50; ++j) {
      ASSERT_EQ(A[i * 50 + j], C[i * 50 + j]);
    }
  }
  Matrix vector(50, 1, 'r');
  Matrix vectorMultiple = A * vector;
  for (int i = 0; i < 50; ++i) {
    ASSERT_EQ(vector[i], vectorMultiple[i]);
  }
}