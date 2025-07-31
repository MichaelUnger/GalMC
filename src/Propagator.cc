#include "Propagator.h"
#include "Configuration.h"
#include "Particle.h"
#include "PDGParticleIds.h"
#include "Gas.h"
#include "Nuclei.h"
#include "TrackingUtilities.h"

#include <utl/Units.h>
#include <utl/PhysicalConstants.h>

#include <TRandom3.h>

#include <iostream>
#include <sstream>
#include <limits>

using namespace std;
using namespace ruqi;
using namespace utl;

#undef _PROPDEBUG_


namespace galmc {

  Propagator::Propagator(const Configuration& c) :
    fConfig(c),
    fMagneticField(c.fMagneticFieldConfig),
    fDiffusion(fConfig.fLmax, fConfig.fDiffusionConstant,
               fConfig.fDiffusionDelta),
    fRandom(fConfig.fSeedPropagation)
  {
    const Particle p(ParticleConst::eProton, 10*GeV, 0, TVector3());
    double z = 0;
    while (z < 15* kpc) {
      const auto bFields = fMagneticField.GetFields(-8.5*kpc, 0, z);
      const double b2 = std::get<MagneticField::eRandomVariance>(bFields) /
        fConfig.fRandBFudge2;
      const auto& Bvec = std::get<MagneticField::eCoherent>(bFields);
      const double B2 = Bvec.SquaredLength();

      const TVector3 B(Bvec.x, Bvec.y, Bvec.z);
      const auto uB = B.Unit();
      const auto u1 = uB.Orthogonal();
      const auto u2 = uB.Cross(u1);

      const auto D = fDiffusion.GetCoefficient(p, B2, b2);
      const TVector3 Dvec = uB*D.first + u1*D.second + u2*D.second;
      cout << z / kpc << " eta = " << b2 / (b2 + B2) << " D =["
           << D.first / (cm2/s) <<", "
           << D.second / (cm2/s) << "] "
           << "Bxy " << sqrt(Bvec.x*Bvec.x + Bvec.y*Bvec.y)/microgauss
           << "Bz " << Bvec.z / microgauss
           << "b " << sqrt(b2) / microgauss
           << endl;
      z+= 1*kpc;
    }
  }

  TVector3
  Propagator::RK4Step(const TVector3 x0, const double delta)
    const
  {
    const TVector3 k1 = delta * GetBcohTangent(x0);
    const TVector3 k2 = delta * GetBcohTangent(x0 + k1 * 0.5);
    const TVector3 k3 = delta * GetBcohTangent(x0 + k2 * 0.5);
    const TVector3 k4 = delta * GetBcohTangent(x0 + k3);
    return 1/6. * (k1 + 2*k2 + 2*k3 + k4);
  }


  double
  Propagator::SafeDistAniso(const Particle& p,
                            const double Dpar, const double Dperp,
                            const double relTol, const double maxDist,
                            const double nStep, const int nRefineMax)
    const
  {
    const auto& pos = p.GetPosition();
    double delta = maxDist / nStep;

    auto pos1 = pos;
    auto pos2 = pos;

    int nRefine = 0;
    double d = delta;
    while (d < maxDist) {

      const auto delta1 = RK4Step(pos1, delta);
      const auto nextPos1 = pos1 + delta1;
      const auto D1 = GetD(nextPos1, p);
      const auto delta2 = RK4Step(pos2, -delta);
      const auto nextPos2 = pos2 + delta2;
      const auto D2 = GetD(nextPos2, p);
      const double rDiff1par = std::abs(get<0>(D1) - Dpar) / Dpar;
      const double rDiff1perp = std::abs(get<1>(D1) - Dperp) / Dperp;
      const double rDiff2par = std::abs(get<0>(D2) - Dpar) / Dpar;
      const double rDiff2perp = std::abs(get<1>(D2) - Dperp) / Dperp;

      /*
        cout << delta/kpc << " [" << pos1.X()/kpc << " " << pos1.Y()/kpc
             << " " << pos1.Z()/kpc << "] " << (pos1-pos).Mag()/kpc << " "
             << rDiff1par << " " << rDiff1perp << " " << Dpar << endl;
        cout << delta/kpc << " [" << pos2.X()/kpc << " " << pos2.Y()/kpc
             << " " << pos2.Z()/kpc << "] " << (pos2-pos).Mag()/kpc << " "
             << rDiff2par << " " << rDiff2perp << "||" << endl;
      */
      if (rDiff1par > relTol || rDiff1perp > relTol ||
          rDiff2par > relTol || rDiff2perp > relTol) {
        if (nRefine >= nRefineMax) {
          return
            std::max(delta,
                     sqrt(std::max((pos1-pos).Mag2(), (pos2-pos).Mag2())));
        }
        else {
          ++nRefine;
          d -= delta;
          delta /= 10;
        }
      }
      else {
        pos1 = nextPos1;
        pos2 = nextPos2;
      }
      d += delta;
    }
    // WARN("reached maxDist");
    return maxDist;
  }

  double
  Propagator::SafeDistIso(const Particle& p,
                          const double Dpar, const double Dperp,
                          const double relTol, const double maxDist,
                          const double nStep, const int nRefineMax)
    const
  {
    const double v = 1/sqrt(3);
    const vector<TVector3> dirs =
      {
       // a) six sides
       TVector3(1, 0, 0),
       TVector3(0, 1, 0),
       TVector3(0, 0, 1),
       TVector3(-1, 0, 0),
       TVector3(0, -1, 0),
       TVector3(0, 0, -1),
       // b) eight corners
       TVector3(-v, -v, -v),
       TVector3(v, -v, -v),
       TVector3(-v, -v, v),
       TVector3(-v, v, -v),
       TVector3(v, v, v),
       TVector3(-v, v, v),
       TVector3(v, v, -v),
       TVector3(v, -v, v)
      };

    const auto pos = p.GetPosition();

    double delta = maxDist / nStep;
    int nRefine = 0;
    double d = delta;
    while (d < maxDist) {

      for (unsigned int i = 0; i < dirs.size(); ++i) {
        const auto& dir = dirs[i];
        const auto nextPos = pos + d*dir;
        const auto D = GetD(nextPos, p);
        const double rDiffpar = std::abs(get<0>(D) - Dpar) / Dpar;
        const double rDiffperp = std::abs(get<1>(D) - Dperp) / Dperp;
        /*
        cout << i << " "
             << delta/kpc << " [" << nextPos.X()/kpc << " " << nextPos.Y()/kpc
             << " " << nextPos.Z()/kpc << "] " << (nextPos-pos).Mag()/kpc << " "
             << rDiffpar << " " << rDiffperp << endl;
        */
        if (rDiffpar > relTol || rDiffperp > relTol) {
          if (nRefine >= nRefineMax)
            return d;
          else {
            ++nRefine;
            d -= delta;
            delta /= 10;
          }
          break;
        }
      }
      d += delta;
    }

    // WARN("reached maxDist");

    return d;
  }

  double
  Propagator::GetSafeDistance(const Particle& p,
                              const double relTol,
                              const double maxDist,
                              const unsigned int nStep)
    const
  {
    const auto pos = p.GetPosition();
    const auto D = GetD(pos, p);
    const double Dpar = get<0>(D);
    const double Dperp = get<1>(D);
    const double eta = get<2>(D);

    if (!Dpar || !Dperp) {
      WARN("zero D at [" +
           to_string(pos.X()/kpc) + ", " +
           to_string(pos.Y()/kpc) + ", " +
           to_string(pos.Z()/kpc) + "] kpc ");
      return -1;
    }

    // follow field line if anisotropic
    if (eta < 0.9)
      return SafeDistAniso(p, Dpar, Dperp, relTol, maxDist, nStep);
    else // near isotropic case
      return SafeDistIso(p, Dpar, Dperp, relTol, maxDist, nStep);
  }


  TVector3
  Propagator::GetBcohTangent(const TVector3& pos)
    const
  {
    const auto bFields =
      fMagneticField.GetFields(pos.X(), pos.Y(), pos.Z());
    const Vector3 bFieldCoh = std::get<MagneticField::eCoherent>(bFields);
    const double mag = bFieldCoh.Length();
    if (mag)
      return TVector3(bFieldCoh.x/mag, bFieldCoh.y/mag, bFieldCoh.z/mag);
    else {
      WARN("zero field at" +
           to_string(pos.X()/kpc) +
           to_string(pos.Y()/kpc) +
           to_string(pos.Z()/kpc));
      return TVector3(0, 0, 0);
    }
  }


  void
  SetPropParameters(const Particle& p, double& sigmaP,
                    double& sigmaHe, double& lDec)
  {
    const auto& nuclei = Nuclei::GetInstance();
    const auto pdgId = p.GetId();
    const auto& prop = nuclei.GetNuclearProperties(pdgId);
    const double lgEkPerNGeV = log10(p.GetKineticEnergy() / p.GetA() / GeV);
    sigmaP = nuclei.GetSigmaInel(pdgId, ParticleConst::eProton, lgEkPerNGeV);
    sigmaHe = nuclei.GetSigmaInel(pdgId, ParticleConst::eHe4, lgEkPerNGeV);
    lDec = prop.IsStable() ? numeric_limits<double>::max() :
      p.GetGamma() * prop.fBetaDecayHalfLife * kSpeedOfLight;
  }

  bool
  DoReaction(Particle& p, const bool decay, const double heProb)
  {
    const auto& nuclei = Nuclei::GetInstance();
    if (decay) {
      const auto& prop = nuclei.GetNuclearProperties(p.GetId());
      const int idDec = prop.fBetaDecayPdgCode;
      /*
        cout << " decay! ["
        << ParticleConst::GetName(p.GetId()) << " --> "
        << ParticleConst::GetName(idDec) << "] X=" << endl;
      */
      if (idDec <= 0) {
        cerr << " wtf decay? " << idDec;
        p.UpdateParticleId(ParticleConst::eProtonNucleus);
        return false;
      }
      p.UpdateParticleId(idDec);
    }
    else {
      const double lgEkPerNuclGeV =
        log10(p.GetKineticEnergy() / p.GetA() / GeV);

      int fragmentId;
      do {
        fragmentId =
          nuclei.RandomFragment(p.GetId(), lgEkPerNuclGeV, heProb);
        if (fragmentId == 0 ||
            ParticleConst::GetNucleusMassNumber(fragmentId) < 3) {
          /*
            cout << " breakup !! ["
            << ParticleConst::GetName(p.GetId()) << "] X="
            << p.GetX() / (g/cm2)
            << endl;
          */
          p.UpdateParticleId(ParticleConst::eProtonNucleus);
          return false;
        }
      }
      while (fragmentId &&
             ParticleConst::GetNucleusMassNumber(fragmentId) < 3);
      /*
        cout << " interaction!! ["
        << ParticleConst::GetName(p.GetId()) << " --> "
        << ParticleConst::GetName(fragmentId) << "] X="
        << p.GetX() / (g/cm2)
        << " t = " << p.GetTime() / megayear << endl;
      */
      p.UpdateParticleId(fragmentId);
    }
    return true;
  }

  tuple<double, double, double>
  Propagator::GetD(const TVector3& pos, const Particle& p)
    const
  {
     const auto bFields =
       fMagneticField.GetFields(pos.X(), pos.Y(), pos.Z());
     const double b02 = std::get<MagneticField::eRandomVariance>(bFields) /
       fConfig.fRandBFudge2;
     const Vector3 bFieldCoh = std::get<MagneticField::eCoherent>(bFields);
     const double B02 = bFieldCoh.SquaredLength();
     if (B02 > 0 || b02 > 0) {
       const double eta = b02 / (B02 + b02);
       const auto D = fDiffusion.GetCoefficient(p, B02, b02);
       return make_tuple(D.first, D.second, eta);
     }
     else
       return make_tuple(0., 0., 0.);
  }

  inline
  void
  PrintParticle(const Particle& p, const double n)
  {
    cout << "p=["
         << p.GetInternalId() << ", "
         << p.GetPosition().X() / kpc << ", "
         << p.GetPosition().Y() / kpc << ", "
         << p.GetPosition().Z() / kpc << ", "
         << p.GetTime() / megayear << ", "
         << p.GetX() / (g/cm2) << ", "
         << n / (1/cm3)
         << "]" << endl;
  }

  inline
  void
  PrintStep(const Particle& p, const int nStep, const double B02, const double b02,
            const double lMax, const double relTime)
  {
    cout << " step "
         << " " << p.GetInternalId()
         << " " << nStep
         << " " << p.GetTime() /megayear
         << " " << p.GetPosition().X()/kpc
         << " " << p.GetPosition().Y()/kpc
         << " " << p.GetPosition().Z()/kpc
         << " " << p.GetDirection().X()
         << " " << p.GetDirection().Y()
         << " " << p.GetDirection().Z()
         << " " << p.GetX() / (g/cm2)
         << " " << sqrt(B02) / microgauss
         << " " << sqrt(b02) / microgauss
         << " " << p.GetLarmorRadius(sqrt(B02+b02))/kpc
         << " " << b02 / (b02 + B02)
         << " " << lMax/kpc
         << " " << relTime
         << endl;
  }

  bool
  Propagator::InDiffusionVolume(const Particle& p)
    const
  {
    const auto& pos = p.GetPosition();
    const double r2 = pow(pos.X(), 2) + pow(pos.Y(), 2);
    return std::abs(pos.Z()) < fConfig.fDiffZmax && r2 < fConfig.fDiffR2max;
  }

  void
  Propagator::Run(Particle& p,
                  std::vector<Particle>& /*particleStack*/,
                  const double time)
    const
  {

#ifdef _PROPDEBUG_
    cout << " -------------- " << endl;
    PrintParticle(p, 0);
#endif

    // sanity check
    if (!InDiffusionVolume(p)) {
      ostringstream warn;
      warn << " particle injected outside diffusion volume at"
           << " " << p.GetPosition().X() / kpc
           << " " << p.GetPosition().Y() / kpc
           << " " << p.GetPosition().Z() / kpc
           << " --> skip! ";
      WARN(warn.str());
      return;
    }

    const double minField2 = pow(1e-6*microgauss, 2);

    const auto& gas = Gas::GetInstance();

    double sigmaP, sigmaHe, lDec;
    SetPropParameters(p, sigmaP, sigmaHe, lDec);

    const bool isotropic = fConfig.fDiffusionConstant > 0;

    const double accuracy = 1e-4;
    const double rL = p.GetLarmorRadius(6*microgauss);
    const double t0 = (2*kPi*rL)/kSpeedOfLight;
    const double dt = t0 / 10;
    const double minStepSize = t0/1e6;
    const double lMin = kAU;
    const double lMax = fConfig.fLmax;
    const unsigned int nWavesPerDecade = 25;
    const unsigned int nWaves = log10(lMax/lMin) * nWavesPerDecade;
    MyRK5ErrorScaling errScaling(rL*accuracy, accuracy);
    MagField mf(fMagneticField, lMin, lMax, nWaves, fConfig.fRandBFudge2);


    double xTrack = 0;
    double lTrack = 0;
    int nStepTrack = 0;
    int nStep = 0;
    while (p.GetTime() < time) {

      #warning disentangle tracking from stepping?
      const double trackLength = fConfig.fMaxStep;
#warning fixed 2.5 sigma (99%)
      const double relTol = 0.1;
      const double maxStep = GetSafeDistance(p, relTol, trackLength*10, 100)/2.575;
      const double maxStep2 = pow(maxStep, 2);

      const auto lastPos = p.GetPosition();
      const auto& pos = p.GetPosition();
      const auto& dir = p.GetDirection();

      const auto bFields =
        fMagneticField.GetFields(pos.X(), pos.Y(), pos.Z());
      const double b02 = std::get<MagneticField::eRandomVariance>(bFields) /
        fConfig.fRandBFudge2;
      const Vector3 bFieldCoh = std::get<MagneticField::eCoherent>(bFields);
      const double B02 = bFieldCoh.SquaredLength();
      if (b02 < minField2 && B02 < minField2) {
#ifdef _PROPDEBUG_
        cout << " zero-field region " << " "
             << pos.X() / kpc << " " << pos.Y() / kpc
             << " " << pos.Z() / kpc << endl;
#endif
            return;
      }
      const auto n = gas.GetNumberDensities(pos);

      const auto D = fDiffusion.GetCoefficient(p, B02, b02);
      const double Dpar = D.first;
      const double Dperp = D.second;

      // step proposal
      const double dTend = time - p.GetTime();
      const double dTdiff =
        std::min(maxStep2 / std::max(Dpar, Dperp) / 2, dTend);
      const double v = p.GetVelocity();

#ifdef _PROPDEBUG_
      if (!p.GetTime())
        PrintStep(p, nStep, B02, b02, lMax, dTdiff / (lMax/v));
#endif
      const bool doTracking = dTdiff / (lMax/v) < fConfig.fTrackingThresh;
      //const bool doTracking = sqrt(2*Dpar*dTdiff)/lMax <  fConfig.fTrackingThresh;
      if (doTracking) {
        // start at random direction if |dir| = 0
        // (i.e. at first step or when returning from diffusive regime)
        if (!p.GetDirection().Mag2()) {
          p.GetDirection() = RandomDirection<TVector3>(fRandom);
        }
        try {
          typedef double Vector6[6];
          typedef ChargeInMagneticFieldODE<MagField> ODE;
          const double q = p.GetZ()*eplus;
          ODE chargeMotionODE(q, p.GetEk(), mf);
          RK5ODEIntegrator<ODE> integrator(chargeMotionODE);

          Vector6 xu = {pos.X(), pos.Y(), pos.Z(), dir.X(), dir.Y(), dir.Z()};

          double t = 0;
          const double dT = trackLength / kSpeedOfLight;
          auto it = integrator.AdaptiveBegin(0, dt/10, xu, 1., errScaling);
          it.SetMinStepSize(minStepSize);
          while (++it && t < dT) {
            if (!it) {
              cerr << " tracking failed !" << endl;
              break;
            }
            const double deltaT = it.GetX() - t;
            t = it.GetX();
            const auto result = it.GetY();
            for (unsigned int i = 0; i < 6; ++i)
              xu[i] = result[i];
            const bool renorm = true;
            if (renorm) {
              const Vector3 uNew = Vector3(xu[3], xu[4], xu[5]).Norm();
              xu[3] = uNew.x;
              xu[4] = uNew.y;
              xu[5] = uNew.z;
            }
            p.GetPosition() = TVector3(xu[0], xu[1], xu[2]);
            p.GetDirection() = TVector3(xu[3], xu[4], xu[5]);
            p.GetTime() += deltaT;
            const double dL = deltaT * kSpeedOfLight;
            p.GetPathLength() += dL;
            const double dX =
              dL*(n.first * kProtonMass + n.second * kHelium4Mass);
            p.GetX() += dX;
            xTrack += dX;
            lTrack += dL;
            ++nStep;
            ++nStepTrack;
          }
#ifdef _PROPDEBUG_
          PrintStep(p, nStep, B02, b02, lMax, dTdiff / (lMax/v));
#endif
        }
        catch (std::runtime_error& err) {
          ostringstream errMsg;
          errMsg << err.what() << " at "
                 << pos.X() << " " << pos.Y() << " " << pos.Z();
          ERROR(errMsg.str());
          p.SetSaveCounter(-999);
          return;
        }
      }
      else { // diffusion
        p.GetDirection() = TVector3();
        TVector3 dCorr;
        try {
          const double dLdiff = v * dTdiff;

          // check interaction or decay
          const double denom = (n.first * sigmaP + n.second * sigmaHe);
          const double lambdaInt = denom > 0 ?
            1 / denom : numeric_limits<double>::max();
          const double dLint = fRandom.Exp(lambdaInt);
          const double dLdec = fRandom.Exp(lDec);

          double dL = dLdiff;
          double dT = dTdiff;
          bool fullyFragmented = false;
          if (dLint < dLdiff || dLdec < dLdiff) {
            dL = std::min(dLint, dLdec);
            dT = dL / v;
            const double heProb = lambdaInt * (n.second * sigmaHe);
            fullyFragmented = !DoReaction(p, dLdec < dLint, heProb);
            SetPropParameters(p, sigmaP, sigmaHe, lDec);
          }

          TVector3 delta;
          const double stepPar = sqrt(2*Dpar*dT);
          #warning use turbulence level instead of minField2 (think extragal!)
          if (B02 < minField2 || isotropic) {
            delta.SetXYZ(fRandom.Gaus(0, stepPar),
                         fRandom.Gaus(0, stepPar),
                         fRandom.Gaus(0, stepPar));
            p.GetPosition() += delta;
            p.GetDirection() = delta.Unit();
          }
          else {
            // derivatives at start position
            const auto uB = GetBcohTangent(p.GetPosition());
            const auto u1 = uB.Orthogonal().Unit();
            const auto u2 = uB.Cross(u1);
            const double eps = stepPar / 100;
            const double dDpar =
              get<0>(GetD(p.GetPosition() + uB*eps, p)) - Dpar;
            const double dDperp1 =
              get<1>(GetD(p.GetPosition() + u1*eps, p)) - Dperp;
            const double dDperp2 =
              get<1>(GetD(p.GetPosition() + u2*eps, p)) - Dperp;
            const double dDparDz = dDpar / eps;
            const double dDparDx = dDperp1 / eps;
            const double dDparDy = dDperp2 / eps;
            dCorr = (dDparDz * uB + dDparDx*u1 + dDparDy*u2) * dT;

            //  propagate parallel to magnetic field
            const double xiB = fRandom.Gaus(0, stepPar);
            const double magStep = 20*pc;
            const int nS = std::max(1, int(std::abs(xiB) / magStep));
            const double dStep = xiB / nS;

            int iS = 0;
            while (iS < nS) {
              const auto delta = RK4Step(p.GetPosition(), dStep);
              p.GetPosition() += delta;
              p.GetDirection() = delta.Unit();
              if (!InDiffusionVolume(p))
                break;
              ++iS;
            }

            if (iS < nS) {
              const double corr = sqrt(iS/double(nS));
              dT *= corr;
              dL *= corr;
            }
            else {
              // apply transverse diffusion
              const double stepPerp = sqrt(2*Dperp*dT);
              const auto uB = GetBcohTangent(p.GetPosition());

              const auto u1 = uB.Orthogonal().Unit();
              const double xi1 = fRandom.Gaus(0, stepPerp);
              p.GetPosition() += xi1 * u1;

              const auto u2 = uB.Cross(u1);
              const double xi2 = fRandom.Gaus(0, stepPerp) ;
              p.GetPosition() += xi2 * u2;
            }
          }

          p.SetCurrentDiff(Dpar, Dperp);
          p.GetTime() += dT;
          p.GetPathLength() += dL;
          const double dX =
            dL*(n.first * kProtonMass + n.second * kHelium4Mass);
          p.GetX() += dX;
#ifdef _PROPDEBUG_
          PrintStep(p, 0, B02, b02, lMax, dTdiff / (lMax/v));
#endif
          ++nStep;
          if (fullyFragmented)
            return;
        }
        catch (std::runtime_error& err) {
          ostringstream errMsg;
          errMsg << err.what() << " at "
                 << pos.X() << " " << pos.Y() << " " << pos.Z();
          ERROR(errMsg.str());
          throw errMsg.str();
        }
      }
#ifdef _PROPDEBUG_
      PrintParticle(p, n.first);
#endif

      if (!InDiffusionVolume(p)) {
#ifdef _PROPDEBUG_
        cout << " left Galaxy at ["
             << p.GetPosition().X() / kpc << ", "
             << p.GetPosition().Y() / kpc << ", "
             << p.GetPosition().Z() / kpc << "]"
             << endl;
#endif
        return;
      }
    }
    return;
  }
}
