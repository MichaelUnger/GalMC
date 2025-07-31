#ifndef _Diffusion_h_
#define _Diffusion_h_

#include <cmath>
#include <iostream>

#include "Particle.h"
#include <utl/Units.h>
#include <utl/PhysicalConstants.h>
#include <utl/DebugOutput.h>

namespace galmc {

  namespace {
    inline
    std::pair<double, double>
    CummingsPD1(const Particle& particle,
                const double D0)
    {
      const double rhoRef = 4.886 * utl::GeV / utl::eplus / utl::kSpeedOfLight;
      const double rho = particle.GetPByQ();
      const double delta = rho < rhoRef ? -0.631 : 0.570;
      const double beta = particle.GetBeta();
      const double D = beta*D0*pow(rho/rhoRef, delta);
      return std::pair<double, double>(D, D);
    }

    inline
    std::pair<double, double>
    SimpleDiffusion(const Particle& particle,
                    const double D0,
                    const double delta)
    {
      const double rhoRef = 10*utl::GeV/ utl::eplus / utl::kSpeedOfLight;
      const double rho = particle.GetPByQ();
      const double beta = particle.GetBeta();
      const double D = beta*D0*pow(rho/rhoRef, delta);
      return std::pair<double, double>(D, D);
    }
  }

  class Diffusion {
  public:
    // isotropic diffusion for Diso > 0
    Diffusion(const double Lmax, const double Diso = -1,
              const double deltaIso = 0.5) :
      fLmax(Lmax), fIsoD(Diso), fDeltaIso(deltaIso) { }

    // returns parallel and perpendicular diffusion coefficient
    std::pair<double, double>
    GetCoefficient(const Particle& particle,
                   const double B02,
                   const double b02)
      const
    {
      if (fIsoD > 0)
        return SimpleDiffusion(particle, fIsoD, fDeltaIso);
      const double B2 = (b02 + B02);
      const double eta = std::max(1e-8, b02 / B2);
      const double rL = particle.GetLarmorRadius(sqrt(B2));
      const double v = particle.GetVelocity();
      const double rho = rL / fLmax;
      return GetCoefficient(rho, eta, rL, v);
    }

    // returns parallel and perpendicular diffusion coefficient
    std::pair<double, double>
    GetCoefficient(const double rho,
                   const double eta,
                   const double rL,
                   const double v)
      const
    {
      if (fIsoD > 0)
        FATAL("not supported");
      const double Dpar = DParByRv(rho, eta) * rL * v;
      const double Dperp = pow(SqrtDPerpByDpar(rho, eta), 2) * Dpar;
      return std::pair<double, double>(Dpar, Dperp);
    }

  private:

    // from test particle simulation
    double
    SqrtDPerpByDpar(const double rho, const double eta)
      const
    {
      const double B0b02 = (1 - eta)/eta;
      const double d00 = 6/49.;
      const double d01 = 0.93;
      const double d0 = pow((B0b02+d00)/d00, -d01);
      const double delta = (1-exp(-5*(1-eta)))/10.;
      const double x0 = 7/25.+11/400.*pow(B0b02, -0.69);
      const double DoverRc =
        d0 * pow(rho/1e-4, delta) / (1 + pow(rho/x0, 1.6));
      return DoverRc;
    }

    // from test particle simulation
    double
    DParByRv(const double rho, const double eta)
      const
    {
      const double d0 = 10 + 176 * (1 - eta) / eta;
      const double d1 = 1/2.*eta;
      const double d2 = 25/4. + 7 *(1-eta) / eta;
      const double alpha0 = -2/3.;
      const double alpha2 = 5/4.;
      const double DoverLc =
        d0 * pow(rho/1e-4, alpha0) +
        d1 +
        d2 * pow(rho, alpha2);
      return DoverLc;
    }
    const double fLmax;
    const double fIsoD;
    const double fDeltaIso;
  };
}

#endif
