#ifndef __UNPACKER2D__
#define __UNPACKER2D__

#include <TObject.h>
#include <map>
#include <string>

#include "Unpacker2.h"

class Unpacker2D : public Unpacker2
{
public:
  void UnpackSingleStep(std::string inputFile, std::string inputPath, std::string outputPath, std::string configFile, int numberOfEvents,
                        int refChannelOffset, std::string TDCcalibFile);

protected:
  void ParseConfigFile();
  void DistributeEventsSingleStep();
  void BuildEvent(EventIII* e, std::map<UInt_t, std::vector<UInt_t>>* m, std::map<UInt_t, UInt_t>* refTimes);
  bool loadTDCcalibFile(const char* calibFile);

  const double kCoarseConstant = 2.3529;
  const double kFineConstant = 0.0201612903;
  const int kModularOffset = 2100;
  const int kNumberOfChannels = 105;
  const int kNumberOfModules = 24;
  const int kReferenceChannel = 104;
};

#endif
