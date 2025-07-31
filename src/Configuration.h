#ifndef _Configuration_h_
#define _Configuration_h_

#include <string>
#include <vector>
#include <TObject.h>

namespace galmc {
  class Configuration : public TObject {
  public:
    Configuration() = default;
    Configuration(const int argc, const char** argv);
  public:
    std::string fOutFilename;
    std::string fMagneticFieldConfig;
    int fSeedMagneticField;
    int fSeedPropagation;
    int fSeedInjection;
    unsigned int fNparticles;
    std::vector<double> fTimes;
    double fXinjection;
    double fYinjection;
    double fZinjection;
    double fRmaxInjection;
    double fRbiasInjection;
    double fAlphaBiasInjection;
    double fH0injection;
    double fDiffusionConstant;
    double fDiffusionDelta;
    double fLgRmin;
    double fLgRmax;
    int fPdgId;
    int fNthreads;
    // observer
    double fObserverX;
    double fObserverY;
    double fObserverZ;
    double fObserverR2;
    double fObserverH = -1;
    // grid for precalculations
    double fXmin;
    double fXmax;
    double fYmin;
    double fYmax;
    double fZmin;
    double fZmax;
    unsigned int fNx;
    unsigned int fNy;
    unsigned int fNz;
    bool fPreCalcGas;

    // step sizes
    double fMaxStep;
    // fudge factors
    double fRandBFudge2;
    // nuclear interactions
    int fSigmaInelOpt;
    int fSigmaFragOpt;
    // maximum diffusion volume
    double fDiffZmax;
    double fDiffR2max;
    // outer scale of turbulence
    double fLmax;
    // threshold between diffusion and tracking
    double fTrackingThresh;
    ClassDefNV(Configuration, 5);
  };

}

#endif
