#ifndef _galmc_TrackingUtilities_h_
#define _galmc_TrackingUtilities_h_

#include <utl/RK5ODEIntegrator.h>
#include <utl/Vector3.h>
#include <ruqi/GJ94Random.h>

// auxillary code for tracking
namespace galmc {

  class MyRK5ErrorScaling {
  public:
    MyRK5ErrorScaling(const double accX, const double accU) :
      fAccX(accX), fAccU(accU) {}

    template<class VectorType>
    void
    operator()(const unsigned int n, const double /*dx*/,
               const VectorType& /*y*/, const VectorType& /*dYdX*/,
               VectorType& errorScaling)
      const
    {
      if (n != 6)
        FATAL("n!=6");

      for (unsigned int i = 0; i < 3; ++i)
        errorScaling[i] = fAccX;
      for (unsigned int i = 3; i < 6; ++i)
        errorScaling[i] = fAccU;

    }
  private:
    const double fAccX;
    const double fAccU;
  };

  template<class MagneticField>
  class ChargeInMagneticFieldODE {
  public:
    ChargeInMagneticFieldODE(const double charge,
                             const double energy,
                             const MagneticField& mag)
      : fQcE(charge * pow(utl::kSpeedOfLight, 2) / energy), fMagField(mag) { }

    /// calculate derivatives
    template<typename Vector6>
    bool
    operator()(const double /*x*/, const Vector6& y, Vector6& dYdX)
      const
    {
      const utl::Vector3& b = fMagField(y[0], y[1], y[2]);
      dYdX[0] = utl::kSpeedOfLight * y[3];  // u_x
      dYdX[1] = utl::kSpeedOfLight * y[4];  // u_y
      dYdX[2] = utl::kSpeedOfLight * y[5];  // u_z
      dYdX[3] = fQcE * (y[4]*b.z - y[5]*b.y);
      dYdX[4] = fQcE * (y[5]*b.x - y[3]*b.z);
      dYdX[5] = fQcE * (y[3]*b.y - y[4]*b.x);
      return true;
    }

    /// number of equations
    operator unsigned int() const
    { return 6; }

  private:
    const double fQcE;
    const MagneticField& fMagField;
  };


  class MagField {

  public:
    MagField(const ruqi::MagneticField& bfield,
             const double lMin, const double lMax,
             const unsigned int nWaves,
             const double randFudge2) :
      fMagneticField(bfield),
      fGJ94(lMin, lMax, nWaves, 1000*lMax, 4, 5./3.),
      fRandBFudge2(randFudge2)
    {
      /*
      cout << " dMax of GJ94 is " << fGJ94.GetDistMaxFast() / pc
           << " pc " << endl;
      cout << " nWaves is " << nWaves << endl;
      */
    }

    utl::Vector3
    operator()(const double x,
               const double y,
               const double z)
      const
    {
      const auto bFields = fMagneticField.GetFields(x, y, z);
      const double b2 =
        std::get<ruqi::MagneticField::eRandomVariance>(bFields) / fRandBFudge2;
      const auto& B = std::get<ruqi::MagneticField::eCoherent>(bFields);
      const auto bRand = fGJ94.GetDeltaBFast(utl::Vector3(x, y, z))*sqrt(b2);
      return bRand + B;
    }

  private:
    const ruqi::MagneticField& fMagneticField;
    const ruqi::GJ94Random fGJ94;
    const double fRandBFudge2;
  };

  template<typename Vector, typename RandGen>
  Vector
  RandomDirection(RandGen& rand)
  {
    const double phi = rand.Uniform(0, utl::kTwoPi);
    const double cosTheta = rand.Uniform(-1, 1);
    const double theta = acos(cosTheta);
    const double sinTheta = sin(theta);
    return Vector(sinTheta*cos(phi), sinTheta*sin(phi), cosTheta).Unit();
  }

}
#endif
