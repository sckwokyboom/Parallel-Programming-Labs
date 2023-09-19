#include <iostream>
#include "SimpleIterationParallelTester.hpp"
#include "mpi.h"
#include "LinpackResult.h"

int main() {
  int rank;
  MPI_Init(nullptr, nullptr);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  SimpleIterationParallelTester tester;
  LinpackResult result = LinpackResult(tester.invokeTest());
  std::cout << result.getResultInGFLOPS() << " in process " << rank << std::endl;
  result.gatherData();
  result.setGatheredResult();
  std::cout << std::endl;
  if (rank == 0) {
    std::cout << "General result: " << result.getResultInGFLOPS() << std::endl;
  }
  MPI_Finalize();
  return 0;
}