#ifndef _Gas_h_
#define _Gas_h_

class TVector3;
class TH3D;
namespace DRAGON {
  class Galaxy;
}

#include <utl/Singleton.h>
#include <utility>
#include <memory>

namespace galmc {

  class Configuration;

  class Gas : public utl::Singleton<Gas> {
  public:
    static void Init(const Configuration& c);
    static Gas& GetInstance();

    static bool IsInitialized() { return fIsInitialized; }

    double GetH2Density(const TVector3& pos) const;
    double GetHIDensity(const TVector3& pos) const;
    double GetHIIDensity(const TVector3& pos) const;
    // number density of <proton, helium>
    std::pair<double, double> GetNumberDensities(const TVector3& pos) const;

    ~Gas();
  private:
    static bool fIsInitialized;
    static TH3D* fH2;
    static TH3D* fHI;
    static TH3D* fHII;
    static std::shared_ptr<DRAGON::Galaxy> fHIModel;
    static std::shared_ptr<DRAGON::Galaxy> fHIIModel;
    static std::shared_ptr<DRAGON::Galaxy> fH2Model;
    static std::shared_ptr<DRAGON::Galaxy> fXCOModel;
  };

}

#endif
