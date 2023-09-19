#ifndef LINPACKTEST_MATRIX_HPP
#define LINPACKTEST_MATRIX_HPP


class Matrix {
public:
  double *array_ = nullptr;
  int height_;
  int width_;
  double *chunk_ = nullptr;
  int *mass_count_s = nullptr;
  int *mass_disp_s = nullptr;
  int count_of_chunk_rows_ = 0;
  int rank_of_process_ = 0;
  int num_of_processes_ = 0;
  bool isSeparated = false;
  bool isGathered = true;
  bool wasBroadCasting = false;

  void createArrayOfChunkSizes();

  void createArrayOfDisps();

public:
  Matrix() = default;

  explicit Matrix(const int &height, const int &width, const char &flag);

  ~Matrix();

  Matrix(const Matrix &other);

  Matrix(const Matrix &&other) noexcept;

  void broadCastMatrixData();

  void scatterMatrixData();

  void gatherMatrixData();

  void randomFillMatrixDiagonal();

  void fillMatrixDiagonalWithNum(const double &num);

  void fillFullMatrixWithNum(const double &num);

  void randomFillFullMatrix();

  Matrix &operator=(Matrix other);

  Matrix &operator+=(const Matrix &other);

  Matrix &operator*=(const Matrix &other);

  Matrix &operator*=(const double &scalar);

  Matrix &operator-=(const Matrix &other);

  double &operator[](const int &index);

  double operator[](const int &index) const;

  double A_1();

  double A_inf();

  void printFullMatrix();

  void printChunk();

  double calculateNormOfVector() const;

  friend Matrix operator+(const Matrix &m1, const Matrix &m2);

  friend Matrix operator*(const Matrix &m1, const Matrix &m2);

  friend Matrix operator*(const double &scalar, const Matrix &matrix);

  friend Matrix operator-(const Matrix &m1, const Matrix &m2);
};

Matrix operator+(const Matrix &m1, const Matrix &m2);

Matrix operator*(const Matrix &m1, const Matrix &m2);

Matrix operator*(const double &scalar, const Matrix &matrix);

Matrix operator-(const Matrix &m1, const Matrix &m2);

#endif //LINPACKTEST_MATRIX_HPP
