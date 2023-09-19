#include "LinpackResult.h"
#include <mpi.h>

LinpackResult::LinpackResult(double flopsResult) : resultInFlops_(flopsResult) {
  MPI_Comm_rank(MPI_COMM_WORLD, &rankOfProcess);
  MPI_Comm_size(MPI_COMM_WORLD, &numOfProcess);
}

std::string LinpackResult::getResultInFLOPS() {
  if (!displayOnlyOnMainProcess) {
    result_ = std::to_string(resultInFlops_) + " FLOPS";
    return result_;
  }


  if (rankOfProcess == 0) {
    return std::to_string(generalResultInFlops_) + " FLOPS";
  }
  return "";
}

std::string LinpackResult::getResultInKFLOPS() {
  if (!displayOnlyOnMainProcess) {
    result_ = std::to_string(resultInFlops_ / 1000.0) + " KFLOPS";
    return result_;
  }
  if (rankOfProcess == 0) {
    return std::to_string(generalResultInFlops_ / 1000.0) + " KFLOPS";
  }
  return "";
}

std::string LinpackResult::getResultInMFLOPS() {
  if (!displayOnlyOnMainProcess) {
    result_ = std::to_string(resultInFlops_ / 1000000.0) + " MFLOPS";
    return result_;
  }
  if (rankOfProcess == 0) {
    return std::to_string(generalResultInFlops_ / 1000000.0) + " MFLOPS";
  }
  return "";
}

std::string LinpackResult::getResultInGFLOPS() {
  if (!displayOnlyOnMainProcess) {
    result_ = std::to_string(resultInFlops_ / 1000000000.0) + " GFLOPS";
    return result_;
  }
  if (rankOfProcess == 0) {
    return std::to_string(generalResultInFlops_/ 1000000000.0) + " GFLOPS";
  }
  return "";
}

std::string LinpackResult::getResultInTFLOPS() {
  if (!displayOnlyOnMainProcess) {
    result_ = std::to_string(resultInFlops_ / 1000000000000.0) + " TFLOPS";
    return result_;
  }
  if (rankOfProcess == 0) {
    return std::to_string(generalResultInFlops_ / 1000000000000.0) + " TFLOPS";
  }
  return "";
}

std::string LinpackResult::getResultInPFLOPS() {
  if (!displayOnlyOnMainProcess) {
    result_ = std::to_string(resultInFlops_ / 1000000000000000.0) + " PFLOPS";
    return result_;
  }
  if (rankOfProcess == 0) {
    return std::to_string(generalResultInFlops_ / 1000000000000000.0) + " PFLOPS";
  }
  return "";
}

std::string LinpackResult::getResultInEFLOPS() {
  if (!displayOnlyOnMainProcess) {
    result_ = std::to_string(resultInFlops_ / 1000000000000000000.0) + " EFLOPS";
    return result_;
  }
  if (rankOfProcess == 0) {
    return std::to_string(generalResultInFlops_ / 1000000000000000000.0) + " EFLOPS";
  }
  return "";
}

void LinpackResult::gatherData() {
  double *resultOnProcess = new double[1]();
  resultOnProcess[0] = resultInFlops_;
  double *generalResult = nullptr;
  if (rankOfProcess == 0) {
    generalResult = new double[numOfProcess]();
  }
  MPI_Gather(resultOnProcess, 1, MPI_DOUBLE, generalResult, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (rankOfProcess == 0) {
    for (int i = 0; i < numOfProcess; ++i) {
      generalResultInFlops_ += generalResult[i];
    }
    delete[] generalResult;
  }
  delete[] resultOnProcess;
}

void LinpackResult::setGatheredResult() {
  displayOnlyOnMainProcess = true;
}
