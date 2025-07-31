#include "GalMC.h"
#include "Configuration.h"
#include "Propagator.h"
#include "Gas.h"
#include "PDGParticleIds.h"
#include "Nuclei.h"

#include <utl/Units.h>
#include <utl/DebugOutput.h>
#include <utl/MathConstants.h>

#include <TRandom3.h>
#include <TFile.h>
#include <TTree.h>

#include <iostream>

#ifdef _WITH_OPENMP_
#include <omp.h>
#endif

using namespace std;
using namespace utl;

namespace galmc {

  GalMC::GalMC(const Configuration& c) :
    fConfig(c),
    fRandom(fConfig.fSeedInjection)
  {
#ifdef _WITH_OPENMP_
    omp_set_num_threads(fConfig.fNthreads);
#endif

    fOutfile = new TFile(fConfig.fOutFilename.c_str(), "RECREATE");
    if (!fOutfile || fOutfile->IsZombie())
      FATAL("cannot open " + fConfig.fOutFilename);
    for (unsigned int i = 0; i < fConfig.fTimes.size(); ++i) {
      ostringstream treeName;
      treeName << "snapShot" << i;
      ostringstream treeTitle;
      treeTitle << "t = " << fConfig.fTimes[i] / megayear << " Myr";
      auto tree = new TTree(treeName.str().c_str(), treeTitle.str().c_str());
      Particle tmp(0, 0, 0, TVector3());
      fParticlePtr = &tmp;
      tree->Branch(Particle::Class_Name(), Particle::Class_Name(),
                   &fParticlePtr);
      fSnapshots.push_back(tree);
      fObserverParticles.push_back(0);
    }
    fPropagator = new Propagator(fConfig);
    Gas::Init(fConfig);
    Nuclei::Init(Nuclei::ESigmaInel(fConfig.fSigmaInelOpt),
                 Nuclei::ESigmaFrag(fConfig.fSigmaFragOpt));
  }

  GalMC::~GalMC()
  {
    fOutfile->cd();
    fConfig.Write();
    fOutfile->Write();
    fOutfile->Close();
    delete fOutfile;
    delete fPropagator;
  }

  void
  GalMC::Run()
  {
    int nPart = 0;
    const unsigned int n = fConfig.fNparticles;
    #pragma omp parallel for
    for (unsigned int i = 0; i < n; ++i) {
      vector<Particle> stack;
      stack.push_back(CreateParticle());
      while (!stack.empty()) {
        auto particle = stack.back();
        stack.pop_back();
        for (unsigned int i = 0; i < fConfig.fTimes.size(); ++i) {
          const double time = fConfig.fTimes[i];
          if (time >= particle.GetTime()) {
            fPropagator->Run(particle, stack, time);
            if (particle.GetTime() >= time) {
              if (AtObserver(particle, i)) {
                #pragma omp critical
                {
                  particle.IncrementSaveCounter();
                  StoreParticle(particle, i);
                }
              }
            }
            else
              break;
          }
        }
      }
      if (nPart && nPart%100000 == 0) {
        #pragma omp critical
        cout << " --> propagating particle # " << nPart
             << " (" << int(nPart / double(fConfig.fNparticles)*1000)/10.
             << "%)\n";
        for (unsigned int iT = 0; iT < fConfig.fTimes.size(); ++iT) {
          #pragma omp critical
          cout << "      t = " << fConfig.fTimes[iT] / megayear
               << " Myr, nTot = " << fSnapshots[iT]->GetEntries()
               << ", nObs = " << fObserverParticles[iT] << endl;
        }
      }
      #pragma omp atomic
      ++nPart;
    }
  }

  namespace {
    template<typename T>
    TVector3
    UniformInDisk(const double radius, const double z, T& rand)
    {
      const double r = sqrt(rand.Uniform()) * radius;
      const double phi = rand.Uniform(0, kTwoPi);
      const double x = r*sin(phi);
      const double y = r*cos(phi);
      return TVector3(x, y, z);
    }
  }

  Particle
  GalMC::CreateParticle()
  {
    static int particleCount = 0;
    const double h = fRandom.Exp(fConfig.fH0injection);
    const double z = fRandom.Uniform() < 0.5 ? -h : h;
    double w = 1;
    TVector3 pos;
    if (fConfig.fRbiasInjection < 0)
      pos = UniformInDisk(fConfig.fRmaxInjection, z, fRandom);
    else {
      const double r0 = fConfig.fRbiasInjection;
      const double R = fConfig.fRmaxInjection;
      const double alpha = fConfig.fAlphaBiasInjection;
      const double beta =
        pow(R/r0, 2) * (1 - alpha)/alpha / (2*(1-(1+R/r0)*exp(-R/r0)));
      double r = 0;
      if (fRandom.Uniform() < alpha) {
        pos = UniformInDisk(R, z, fRandom);
        r = sqrt(pow(pos.X()-fConfig.fXinjection, 2) +
                 pow(pos.Y()-fConfig.fYinjection, 2));
      }
      else {
        double rGal = 2*fConfig.fRmaxInjection;
        while (rGal > fConfig.fRmaxInjection) {
          /*
            dn/dr ~ r * exp(-r/r0)
            random number according to
            Sec. 17.6.1 of Walck's Handbook of statistical distributions
          */
          r = -log(fRandom.Uniform()*fRandom.Uniform())*r0;
          const double phi = fRandom.Uniform(0, kTwoPi);
          const double x = r*sin(phi);
          const double y = r*cos(phi);
          pos = TVector3(fConfig.fXinjection,
                         fConfig.fYinjection,
                         fConfig.fZinjection) +
                                          TVector3(x, y, z);
          rGal = sqrt(pow(pos.X(), 2) + pow(pos.Y(), 2));
        }
      }
      w = 1 + beta * exp(-r/r0);
    }

    const double t = fRandom.Uniform(0, fConfig.fTimes.back());
    const double Z  = ParticleConst::GetNucleusCharge(fConfig.fPdgId);
    const double lgR = fRandom.Uniform(fConfig.fLgRmin, fConfig.fLgRmax);
    const double R = pow(10, lgR);
    const double p = R * Z * eplus;
    const auto& properties =
      Nuclei::GetInstance().GetNuclearProperties(fConfig.fPdgId);
    const double m = properties.fMassAmu * kAtomicMass;
    const double mc2 = m * kSpeedOfLight2;
    const double m2c4 = pow(mc2, 2);
    const double E = sqrt(p*p*kSpeedOfLight2 + m2c4);
    const double Ek = E - mc2;
    Particle particle(fConfig.fPdgId, Ek, t, pos, 1./w, particleCount);
    ++particleCount;
    return particle;
  }

  bool
  GalMC::AtObserver(const Particle& particle, const int iTime)
  {
    const auto& pos = particle.GetPosition();
    if (abs(pos.Z() - fConfig.fObserverZ) < fConfig.fObserverH) {
      const double r2 =
        pow(pos.X() - fConfig.fObserverX, 2) +
        pow(pos.Y() - fConfig.fObserverY, 2);
      if (r2 < fConfig.fObserverR2) {
        ++fObserverParticles[iTime];
        return true;
      }
    }
    return false;
  }

  void
  GalMC::StoreParticle(const Particle& p, const int i)
  {
    Particle pp(p);
    pp.ConvertUnits();
    fParticlePtr = &pp;
    fSnapshots[i]->Fill();
  }

}
