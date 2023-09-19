#include <cmath>
#include <mpi.h>
#include <iostream>

#include "matrix_error.hpp"
#include "Matrix.hpp"


void Matrix::createArrayOfChunkSizes() {
  mass_count_s = new int[num_of_processes_]();
  int countOfChunkRows = height_ / num_of_processes_;
  for (int i = 0; i < num_of_processes_; ++i) {
    mass_count_s[i] = countOfChunkRows * width_;
  }
  for (int i = 0; i < height_ % num_of_processes_; ++i) {
    mass_count_s[i] += width_;
  }
}

void Matrix::createArrayOfDisps() {
  mass_disp_s = new int[num_of_processes_]();
  int countOfChunkRows = height_ / num_of_processes_;
  for (int i = 0; i < num_of_processes_; ++i) {
    mass_disp_s[i] = countOfChunkRows * width_ * i;
  }
  for (int i = 1; i < height_ % num_of_processes_; ++i) {
    mass_disp_s[i] += width_;
  }
  for (int i = height_ % num_of_processes_; i < num_of_processes_; ++i) {
    mass_disp_s[i] += width_ * (height_ % num_of_processes_);
  }
}

Matrix::Matrix(const int &height, const int &width, const char &flag = '-') : height_(height), width_(width) {
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_of_process_);
  MPI_Comm_size(MPI_COMM_WORLD, &num_of_processes_);

  if (rank_of_process_ == 0) {
    array_ = new double[height_ * width_]();
    if (flag == 'n') {
      fillFullMatrixWithNum(0);
    }
    if (flag == 'e') {
      fillMatrixDiagonalWithNum(1);
    }
    if (flag == 'r') {
      randomFillFullMatrix();
    }
    if (flag == 't') {
      fillFullMatrixWithNum(1.0);
      fillMatrixDiagonalWithNum(2.0);
    }
    if (flag == 'h') {
      fillMatrixDiagonalWithNum(height_);
    }
  }
}

Matrix::~Matrix() {
  delete[] mass_disp_s;
  delete[] mass_count_s;
  delete[] chunk_;
  if (rank_of_process_ == 0 || wasBroadCasting) {
    delete[] array_;
  }
}

Matrix::Matrix(const Matrix &other) : height_(other.height_), width_(other.width_),
                                      count_of_chunk_rows_(other.count_of_chunk_rows_),
                                      isSeparated(other.isSeparated),
                                      isGathered(other.isGathered),
                                      rank_of_process_(other.rank_of_process_),
                                      num_of_processes_(other.num_of_processes_),
                                      wasBroadCasting(other.wasBroadCasting) {
  if (rank_of_process_ == 0 || wasBroadCasting) {
    array_ = new double[height_ * width_];
    for (int i = 0; i < height_; ++i) {
      for (int j = 0; j < width_; ++j) {
        array_[i * width_ + j] = other.array_[i * width_ + j];
      }
    }
  }

  if (isSeparated) {
    chunk_ = new double[count_of_chunk_rows_ * width_];
    for (int i = 0; i < count_of_chunk_rows_; ++i) {
      for (int j = 0; j < width_; ++j) {
        chunk_[i * width_ + j] = other.chunk_[i * width_ + j];
      }
    }
    mass_disp_s = new int[num_of_processes_];
    mass_count_s = new int[num_of_processes_];
    for (int i = 0; i < num_of_processes_; ++i) {
      mass_disp_s[i] = other.mass_disp_s[i];
      mass_count_s[i] = other.mass_count_s[i];
    }
  }
}

Matrix::Matrix(const Matrix &&other) noexcept: height_(other.height_), width_(other.width_),
                                               count_of_chunk_rows_(other.count_of_chunk_rows_),
                                               isSeparated(other.isSeparated),
                                               isGathered(other.isGathered),
                                               rank_of_process_(other.rank_of_process_),
                                               num_of_processes_(other.num_of_processes_),
                                               wasBroadCasting(other.wasBroadCasting) {
  if (rank_of_process_ == 0) {
    array_ = other.array_;
  }

  if (isSeparated) {
    chunk_ = other.chunk_;
    mass_count_s = other.mass_count_s;
    mass_disp_s = other.mass_disp_s;
  }
}

Matrix &Matrix::operator=(Matrix other) {
  std::swap(other.array_, array_);
  std::swap(other.chunk_, chunk_);
  std::swap(other.mass_disp_s, mass_disp_s);
  std::swap(other.mass_count_s, mass_count_s);
  width_ = other.width_;
  height_ = other.height_;
  isSeparated = other.isSeparated;
  isGathered = other.isGathered;
  rank_of_process_ = other.rank_of_process_;
  num_of_processes_ = other.num_of_processes_;
  count_of_chunk_rows_ = other.count_of_chunk_rows_;
  wasBroadCasting = other.wasBroadCasting;
  return *this;
}

Matrix &Matrix::operator+=(const Matrix &other) {
  if (width_ != other.width_ || height_ != other.height_) {
    throw matrix_error("Invalid operation +=, wrong matrix dimensions.");
  }

  if (isSeparated && other.isSeparated) {
    for (int i = 0; i < count_of_chunk_rows_; ++i) {
      for (int j = 0; j < width_; ++j) {
        chunk_[i * width_ + j] += other.chunk_[i * width_ + j];
      }
    }
  }

  if (!isSeparated && !other.isSeparated) {
    for (int i = 0; i < height_; ++i) {
      for (int j = 0; j < width_; ++j) {
        array_[i * width_ + j] += other.array_[i * width_ + j];
      }
    }
  }
  isGathered = false;
  return *this;
}

Matrix &Matrix::operator*=(const Matrix &other) {
  if (width_ != other.height_) {
    throw matrix_error("Invalid operation *=, wrong matrix dimensions.");
  }

  if (isSeparated && other.isSeparated) {
    throw matrix_error("Invalid multiple, the both matrices are separated.");
  }

  if (other.array_ == nullptr || !other.wasBroadCasting) {
    throw matrix_error("Invalid multiple, need to broad cast matrix.");
  }

  Matrix temp(height_, other.width_, 'n');
  temp.scatterMatrixData();

  if (isSeparated && !other.isSeparated) {
    for (int i = 0; i < temp.count_of_chunk_rows_; ++i) {
      double *c = temp.chunk_ + i * other.width_;

      for (int k = 0; k < other.height_; ++k) {
        const double *b = other.array_ + k * other.width_;
        double a = chunk_[i * width_ + k];
        for (int j = 0; j < other.width_; ++j)
          c[j] += a * b[j];
      }
    }
  }

  if (!isSeparated && !other.isSeparated) {
    for (int i = 0; i < height_; ++i) {
      double *c = temp.array_ + i * other.width_;

      for (int k = 0; k < other.height_; ++k) {
        const double *b = other.array_ + k * other.width_;
        double a = array_[i * other.height_ + k];
        for (int j = 0; j < other.width_; ++j)
          c[j] += a * b[j];
      }
    }
  }
  isGathered = false;
  temp.gatherMatrixData();
  return (*this = temp);
}

Matrix &Matrix::operator*=(const double &scalar) {
  if (isSeparated) {
    for (int i = 0; i < count_of_chunk_rows_; ++i) {
      for (int j = 0; j < width_; ++j) {
        chunk_[i * width_ + j] *= scalar;
      }
    }
  }

  if (!isSeparated) {
    for (int i = 0; i < height_; ++i) {
      for (int j = 0; j < width_; ++j) {
        array_[i * width_ + j] *= scalar;
      }
    }
  }
  isGathered = false;
  return *this;
}

Matrix &Matrix::operator-=(const Matrix &other) {
  if (width_ != other.width_ || height_ != other.height_) {
    throw matrix_error("Invalid operation -=, wrong matrix dimensions.");
  }

  if (isSeparated && other.isSeparated) {
    for (int i = 0; i < count_of_chunk_rows_; ++i) {
      for (int j = 0; j < width_; ++j) {
        chunk_[i * width_ + j] -= other.chunk_[i * width_ + j];
      }
    }
  }

  if (wasBroadCasting && !other.wasBroadCasting) {
    for (int i = other.mass_disp_s[other.rank_of_process_] / width_;
         i < other.mass_count_s[other.rank_of_process_] / width_ +
             other.mass_disp_s[other.rank_of_process_] / width_; ++i) {
      for (int j = 0; j < width_; ++j) {
        array_[i * width_ + j] -= other.chunk_[(i - other.mass_disp_s[other.rank_of_process_] / width_) * width_ + j];
      }
    }
  }

  if (!wasBroadCasting && other.wasBroadCasting) {
    for (int i = mass_disp_s[rank_of_process_] / width_; i < mass_count_s[rank_of_process_] / width_ +
                                                             mass_disp_s[rank_of_process_] / width_; ++i) {
      for (int j = 0; j < width_; ++j) {
        chunk_[(i - mass_disp_s[rank_of_process_] / width_) * width_ + j] -= other.array_[i * width_ + j];
      }
    }
  }


  if (!isSeparated && !other.isSeparated && other.wasBroadCasting && wasBroadCasting) {
    for (int i = 0; i < height_; ++i) {
      for (int j = 0; j < width_; ++j) {
        array_[i * width_ + j] -= other.array_[i * width_ + j];
      }
    }
  }
  isGathered = false;
  return *this;
}

double &Matrix::operator[](const int &index) {
  return chunk_[index];
}

double Matrix::operator[](const int &index) const {
  return chunk_[index];
}

void Matrix::printFullMatrix() {
  if (rank_of_process_ == 0) {
    if (!isGathered) {
      std::cout << "Please, gather matrix" << std::endl;
      return;
    }
    for (size_t i = 0; i < height_; ++i) {
      for (size_t j = 0; j < width_; ++j) {
        std::cout << array_[i * width_ + j] << ' ';
      }
      std::cout << std::endl;
    }
  }
}

void Matrix::printChunk() {
  for (size_t i = 0; i < count_of_chunk_rows_; ++i) {
    for (size_t j = 0; j < width_; ++j) {
      std::cout << chunk_[i * width_ + j] << ' ';
    }
    std::cout << std::endl;
  }
}

double Matrix::calculateNormOfVector() const {
  if (rank_of_process_ == 0 || wasBroadCasting) {
    if (width_ != 1) {
      return -1;
    }
    double result = 0;
    for (int i = 0; i < height_; ++i) {
      result += array_[i] * array_[i];
    }
    result = sqrt(result);
    return result;
  }
  return -1;
}

Matrix operator+(const Matrix &m1, const Matrix &m2) {
  if (m1.width_ != m2.width_ || m1.height_ != m2.height_) {
    throw matrix_error("Invalid operation +, wrong matrix dimensions.");
  }
  Matrix temp(m1);
  temp += m2;

  return temp;
}

Matrix operator*(const Matrix &m1, const Matrix &m2) {
  if (m1.width_ != m2.height_) {
    throw matrix_error("Invalid operation *, wrong matrix dimensions.");
  }

  Matrix temp = Matrix(m1);
  temp *= m2;
  return temp;
}

Matrix operator*(const double &scalar, const Matrix &matrix) {
  Matrix temp = matrix;
  temp *= scalar;
  return temp;
}

Matrix operator-(const Matrix &m1, const Matrix &m2) {
  if (m1.width_ != m2.width_ || m1.height_ != m2.height_) {
    throw matrix_error("Invalid operation -, wrong matrix dimensions.");
  }
  Matrix temp(m1);
  temp -= m2;
  return temp;
}

void Matrix::fillMatrixDiagonalWithNum(const double &num) {
  if (rank_of_process_ == 0 || wasBroadCasting) {
    for (long long i = 0; i < height_; ++i) {
      for (long long j = 0; j < width_; ++j) {
        if (i == j) {
          array_[i * width_ + j] = num;
        }
      }
    }
  }
}

void Matrix::randomFillFullMatrix() {
  if (rank_of_process_ == 0 || wasBroadCasting) {
    for (long long i = 0; i < height_; ++i) {
      for (long long j = 0; j < width_; ++j) {
        array_[i * width_ + j] = (double) random();
      }
    }
  }
}

void Matrix::fillFullMatrixWithNum(const double &num) {
  if (rank_of_process_ == 0 || wasBroadCasting) {
    for (long long i = 0; i < height_; ++i) {
      for (long long j = 0; j < width_; ++j) {
        array_[i * width_ + j] = (double) num;
      }
    }
  }
}

void Matrix::randomFillMatrixDiagonal() {
  if (rank_of_process_ == 0 || wasBroadCasting) {
    for (long long i = 0; i < height_; ++i) {
      for (long long j = 0; j < width_; ++j) {
        if (i == j) {
          array_[i * width_ + j] = (double) random();
        }
      }
    }
  }
}

void Matrix::broadCastMatrixData() {
  if (wasBroadCasting) {
    return;
  }
  if (rank_of_process_ != 0) {
    array_ = new double[height_ * width_];
  }
  MPI_Bcast(array_, height_ * width_, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  wasBroadCasting = true;
}

void Matrix::scatterMatrixData() {
  isSeparated = true;
  isGathered = false;
  count_of_chunk_rows_ =
          height_ / num_of_processes_ +
          (rank_of_process_ < (height_ % num_of_processes_) ? 1 : 0);
  chunk_ = new double[count_of_chunk_rows_ * width_]();
  createArrayOfChunkSizes();
  createArrayOfDisps();
  if (rank_of_process_ == 0) {
    if (MPI_SUCCESS ==
        MPI_Scatterv(array_, mass_count_s, mass_disp_s, MPI_DOUBLE, chunk_, count_of_chunk_rows_ * width_, MPI_DOUBLE,
                     0,
                     MPI_COMM_WORLD)) {
    } else {
    }
  }

  if (rank_of_process_ != 0) {
    if (MPI_SUCCESS ==
        MPI_Scatterv(nullptr, nullptr, nullptr, MPI_DOUBLE, chunk_, count_of_chunk_rows_ * width_, MPI_DOUBLE, 0,
                     MPI_COMM_WORLD)) {
    } else {
    }
  }

}

void Matrix::gatherMatrixData() {
  if (MPI_SUCCESS !=
      MPI_Gatherv(chunk_, count_of_chunk_rows_ * width_, MPI_DOUBLE, array_, mass_count_s, mass_disp_s, MPI_DOUBLE, 0,
                  MPI_COMM_WORLD)) {
    std::cout << "Sad!" << std::endl;
  }
  isGathered = true;
}
