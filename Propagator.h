#ifndef _Propagator_h_
#define _Propagator_h_

#include <ruqi/MagneticField.h>
#include <TRandom3.h>
#include "Diffusion.h"
#include <tuple>

namespace galmc {

  class Configuration;
  class Particle;
  class Diffusion;

  class Propagator {
  public:
    Propagator(const Configuration& c);
    /// propagate particle to time t and fill secondaries to stack
    void Run(Particle& p, std::vector<Particle>& stack,
             const double t) const;
    /// calculate "safe distance" following coherent field given
    /// relative tolerance on diffusion coefficient
    double GetSafeDistance(const Particle& p, const double relTol,
                           const double maxDist = 10 * utl::kpc,
                           const unsigned int nStep = 100) const;

    const Diffusion& GetDiffusion() const { return fDiffusion; }
  private:
    /// tangent unit vector along coherent magnetic field
    TVector3 GetBcohTangent(const TVector3& pos) const;
    /// RK4 step along coherent field
    /// (returns delta of x1 = x0 + delta)
    TVector3 RK4Step(const TVector3 x0, const double dx) const;
    /// diffusion coefficient and turbulence level (Dpar, Dperp, eta)
    std::tuple<double, double, double> GetD(const TVector3& pos, const Particle& p) const;
    /// safe distance anisotropic case
    double SafeDistAniso(const Particle& p,
                         const double Dpar, const double Dperp,
                         const double relTol, const double maxDist,
                         const double nStep, const int nRefineMax = 1) const;
    /// safe distance isotropic case
    double SafeDistIso(const Particle& p,
                       const double Dpar, const double Dperp,
                       const double relTol, const double maxDist,
                       const double nStep, const int nRefineMax = 1) const;
    /// check if inside diffusive volume
    bool InDiffusionVolume(const Particle& p) const;

    const Configuration& fConfig;
    const ruqi::MagneticField fMagneticField;
    const Diffusion fDiffusion;
    mutable TRandom3 fRandom;
  };

}

#endif
