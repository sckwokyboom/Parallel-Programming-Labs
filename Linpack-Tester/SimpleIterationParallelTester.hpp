#ifndef LINPACKTESTER_HPP_SIMPLEITERATIONPARALLELTESTER_HPP
#define LINPACKTESTER_HPP_SIMPLEITERATIONPARALLELTESTER_HPP

#include <mpi.h>
//#include <chrono>
#include "LinpackTester.hpp"

class SimpleIterationParallelTester : public LinpackTester {
private:
  const int matrixHeight_;
  long long floatOperationsCounter = 0;
  const double epsilon_;
  const double tau_;
  Matrix A;
  Matrix x;
  Matrix b;
  double normVector_b = 0;
  int rank_ = 0;
  char *isNeedToTerminate = nullptr;
  double flopsResult = 0;
  double timeInSec = 0;
  static const int countOfFloatOpsToCalcSqrt = 6;
  const int countOfFloatOpsToCalcNorm;
  const int countOfFloatOpsToApproximateAlgValue;
public:
  explicit SimpleIterationParallelTester(const int &matrixSize = 700,
                                         const double &epsilon = 0.1,
                                         const double &tau = 0.0001) :
          matrixHeight_(matrixSize),
          epsilon_(epsilon),
          tau_(tau),
          A(Matrix(matrixSize, matrixSize, 't')),
          x(Matrix(matrixSize, 1, 'n')),
          b(Matrix(matrixSize, 1, 'h')),
          countOfFloatOpsToCalcNorm(2 * matrixSize),
          countOfFloatOpsToApproximateAlgValue(2 * matrixSize * matrixSize + 2 * matrixSize) {

    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    A.scatterMatrixData();
    b.broadCastMatrixData();
    x.broadCastMatrixData();
    normVector_b = b.calculateNormOfVector();
    floatOperationsCounter += countOfFloatOpsToCalcNorm + countOfFloatOpsToCalcSqrt;
    isNeedToTerminate = new char(0);
  }

  ~SimpleIterationParallelTester() override {
    delete isNeedToTerminate;
  }

  double invokeTest() override {
    // The unique number of the current process in the MPI_COMM_WORLD communicator.


    // The total number of processes in the MPI_COMM_WORLD communicator.
    int commSize;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    double start;
    while (!isTerminationCriterion()) {
      start = MPI_Wtime();
      x -= (tau_ * ((A * x) - b));
      timeInSec += (MPI_Wtime() - start);
      floatOperationsCounter += countOfFloatOpsToApproximateAlgValue;
    }

    flopsResult = (double) floatOperationsCounter / timeInSec;
    floatOperationsCounter = 0;
    timeInSec = 0;
    return flopsResult;
  }

  double getTimeResultInSec() {
    return timeInSec;
  }

private:
  bool isTerminationCriterion() {
    const long countOfFloatOpsToCalcTempMatrix = 2 * matrixHeight_ * matrixHeight_;
    floatOperationsCounter += countOfFloatOpsToCalcTempMatrix;

    double start = MPI_Wtime();
    Matrix temp = A * x - b;
    timeInSec += (MPI_Wtime() - start);
    start = MPI_Wtime();
    temp.gatherMatrixData();

    if (rank_ == 0) {
      double terminationCriterionNumber = (temp.calculateNormOfVector()) / normVector_b;
//      std::cout << terminationCriterionNumber << std::endl;
      floatOperationsCounter += countOfFloatOpsToCalcNorm + countOfFloatOpsToCalcSqrt + 1;
      if (terminationCriterionNumber < epsilon_) {
        *isNeedToTerminate = 1;
      }
    }
    timeInSec += (MPI_Wtime() - start);
    MPI_Bcast(isNeedToTerminate, 1, MPI_CHAR, 0, MPI_COMM_WORLD);
    return *isNeedToTerminate;
  }
};

#endif //LINPACKTESTER_HPP_SIMPLEITERATIONPARALLELTESTER_HPP
