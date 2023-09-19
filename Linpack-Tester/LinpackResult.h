#ifndef MAIN_CPP_LINPACKRESULT_H
#define MAIN_CPP_LINPACKRESULT_H

#include <string>

class LinpackResult {
private:
  std::string result_;
  double resultInFlops_ = 0;
  double generalResultInFlops_ = 0;
  bool displayOnlyOnMainProcess = false;
  int rankOfProcess = 0;
  int numOfProcess = 0;
public:
  explicit LinpackResult(double flopsResult);

  std::string getResultInFLOPS();

  std::string getResultInKFLOPS();

  std::string getResultInMFLOPS();

  std::string getResultInGFLOPS();

  std::string getResultInTFLOPS();

  std::string getResultInPFLOPS();

  std::string getResultInEFLOPS();

  void gatherData();

  void setGatheredResult();

};


#endif //MAIN_CPP_LINPACKRESULT_H
