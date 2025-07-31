#ifndef _GalMC_h_
#define _GalMC_h_

#include "Particle.h"
#include <TRandom3.h>

class TFile;
class TTree;

namespace galmc {

  class Configuration;
  class Particle;
  class Propagator;

  class GalMC {
  public:
    GalMC(const Configuration& c);
    ~GalMC();
    void Run();

  private:
    bool AtObserver(const Particle& p, const int iSnapShot);
    Particle CreateParticle();
    void StoreParticle(const Particle& p, const int iSnapShot);
    const Configuration& fConfig;
    TFile* fOutfile;
    const Particle* fParticlePtr;
    const Propagator* fPropagator;
    std::vector<TTree*> fSnapshots;
    std::vector<unsigned int> fObserverParticles;
    TRandom3 fRandom;
  };

}

#endif
