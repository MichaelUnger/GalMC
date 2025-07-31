#include "Gas.h"
#include "Configuration.h"

#include <d2mw.h>
#include <galaxy.h>

#include <utl/Units.h>
#include <utl/PhysicalConstants.h>
#include <utl/DebugOutput.h>

#include <iostream>
#include <sstream>

#include <TVector3.h>
#include <TH3D.h>
#include <TROOT.h>
#include <TGraph.h>

using namespace std;

namespace galmc {

  bool Gas::fIsInitialized = false;
  TH3D* Gas::fH2 = nullptr;
  TH3D* Gas::fHII = nullptr;
  TH3D* Gas::fHI = nullptr;
  std::shared_ptr<DRAGON::Galaxy> Gas::fHIModel  = nullptr;
  std::shared_ptr<DRAGON::Galaxy> Gas::fHIIModel  = nullptr;
  std::shared_ptr<DRAGON::Galaxy> Gas::fH2Model  = nullptr;
  std::shared_ptr<DRAGON::Galaxy> Gas::fXCOModel  = nullptr;

  inline
  double
  Interpolate(TH3D& h, const TVector3& pos)
  {
    const double x = pos.X();
    const double y = pos.Y();
    const double z = pos.Z();
    if (x > h.GetXaxis()->GetBinCenter(1) &&
        x < h.GetXaxis()->GetBinCenter(h.GetNbinsX()) &&
        y > h.GetYaxis()->GetBinCenter(1) &&
        y < h.GetYaxis()->GetBinCenter(h.GetNbinsY()) &&
        z > h.GetZaxis()->GetBinCenter(1) &&
        z < h.GetZaxis()->GetBinCenter(h.GetNbinsZ()))
      return h.Interpolate(x, y, z);
    else
      return 0;
  }

  double
  GetKachel(const TVector3& p)
  {
    static vector<double> ri = {0., 1., 2., 4., 6., 8.5, 10., 11., 12.};
    static vector<double> Ni = {10.,3.,1.5,0.7,0.5,0.3,0.1,0.02,0.001};
    static TGraph gKachel(ri.size(), &ri.front(), &Ni.front());
    const double r = sqrt(pow(p.X(), 2) + pow(p.Y(), 2));
    const double zKachel = exp(-pow(p.Z()/utl::kpc/0.21, 2));
    const double nKachel = gKachel.Eval(r/utl::kpc);
    return std::max(1e-4, zKachel*nKachel) *1/utl::cm3;
  }

  std::pair<double, double>
  Gas::GetNumberDensities(const TVector3& p)
    const
  {
    const double heliumFraction = 0.1;
    if (fH2) {
      const double nH2 = GetH2Density(p);
      const double nHI = GetHIDensity(p);
      const double nHII = GetHIIDensity(p);
      const double nP = 2*nH2 + nHI + nHII;
      return std::pair<double, double>(nP, heliumFraction*nP);
    }
    else {
      const Vector3d pos(p.X()/utl::kpc*kpc,
                         p.Y()/utl::kpc*kpc,
                         p.Z()/utl::kpc*kpc);
      const double dragonToRUQI = (1/utl::cm3)/(1/cm3);
      const double nHI = fHIModel->get(pos)*dragonToRUQI;
#warning AAAAAAAAAAAAA
      const double nHII = 0;//fHIIModel->get(pos)*dragonToRUQI;
      const double nH2 = 0;//fXCOModel->get(pos)*fH2Model->get(pos)*dragonToRUQI;

#warning AAAAAAAAAAAAAAAAA
      return std::pair<double, double>(GetKachel(p), 0);

      const double nP = 2*nH2 + nHI + nHII;
      return std::pair<double, double>(nP, heliumFraction*nP);
    }
  }

  Gas&
  Gas::GetInstance()
  {
    auto& gas = Singleton<Gas>::GetInstance();
    if (!gas.IsInitialized())
      FATAL("before you call GetInstance, "
            "you must call Init()");
    return gas;
  }

  void
  Gas::Init(const Configuration& config)
  {
    DRAGON::D2MW mw;
    mw.set_HI("Nakanishi2003");
    //mw.set_HI("Ferriere2007");
    mw.set_H2("Bronfmann1988");
    mw.set_XCO("Evoli2012");
    mw.set_HII("Cordes1991");
    //mw.set_HII("Ferriere2007");

    fHIModel = mw.create_HI();
    fH2Model = mw.create_H2();
    fHIIModel = mw.create_HII();
    fXCOModel = mw.create_XCO();

    auto HI = fHIModel;
    auto H2 = fH2Model;
    auto HII = fHIIModel;
    auto XCO = fXCOModel;

    const double dragonToRUQI = (1/utl::cm3)/(1/cm3);
    vector<double> rVec = {-7.5*kpc, -8.5*kpc, -9.5*kpc};
    for (const auto r : rVec) {
      double sumI = 0;
      double sum2 = 0;
      double sumII = 0;
      const double dz = 0.025*kpc;
      for (double z = -20 * kpc; z <= 20 * kpc; z += dz) {
        #warning AAAAAAAAAAA
        auto tmp = TVector3(r/kpc*utl::kpc, 0, z/kpc*utl::kpc);
        sumI += GetKachel(tmp)/dragonToRUQI;
        //  auto pos = Vector3d(r, 0, z);
      //        sumI += HI->get(pos);
       //        sum2 += XCO->get(pos) * H2->get(pos);
        //        sumII += HII->get(pos);
      }
      const double delta = dz/kpc*utl::kpc;
      sumI *= (dragonToRUQI*utl::kProtonMass)*delta;
      sumII *= (dragonToRUQI*utl::kProtonMass)*delta;
      sum2 *= (dragonToRUQI*utl::kProtonMass)*delta;
      cout << " column density r = " << r/kpc << " kpc is ("
           << sumI / (utl::kSolarMass/utl::pc2) << " + "
           << 2*sum2 / (utl::kSolarMass/utl::pc2) << " + "
           << sumII / (utl::kSolarMass/utl::pc2) << ") = "
           << (sumI + 2*sum2 + sumII) / (utl::kSolarMass/utl::pc2)
           << " Msun/pc2 " << endl;
      /*
      cout << " column density is ("
         << sumI / (utl::mg/utl::cm2) << " + "
           << 2*sum2 / (utl::mg/utl::cm2) << " + "
           << sumII / (utl::mg/utl::cm2) << ") = "
           << (sumI + 2*sum2 + sumII) / (utl::mg/utl::cm2)
           << " mg/cm2 " << endl;
      cout << " column density is ("
           << sumI / (utl::kProtonMass/utl::cm2) << " + "
           << 2*sum2 / (utl::kProtonMass/utl::cm2) << " + "
           << sumII / (utl::kProtonMass/utl::cm2) << ") = "
           << (sumI + 2*sum2 + sumII) / (utl::kProtonMass/utl::cm2)
           << " H/cm2 " << endl;
      */
    }

    if (config.fPreCalcGas) {
      const auto currDir = gDirectory;
      gROOT->cd();
      ostringstream s;
      s << "3*(" << config.fNx << "*" << config.fNy << "*" << config.fNz
        << ") = "
        << 3*config.fNx*config.fNy*config.fNz*sizeof(double)/1024/1024 << " MB";
      INFO("allocating gas maps: " + s.str());
      fH2 = new TH3D("H2", "",
                     config.fNx, config.fXmin, config.fXmax,
                     config.fNy, config.fYmin, config.fYmax,
                     config.fNz, config.fZmin, config.fZmax);
      fHII = new TH3D("HII", "",
                      config.fNx, config.fXmin, config.fXmax,
                      config.fNy, config.fYmin, config.fYmax,
                      config.fNz, config.fZmin, config.fZmax);
      fHI = new TH3D("HI", "",
                     config.fNx, config.fXmin, config.fXmax,
                     config.fNy, config.fYmin, config.fYmax,
                     config.fNz, config.fZmin, config.fZmax);
      INFO("filling gas maps");
      for (unsigned int i = 0; i < config.fNx; ++i) {
        const double x = fH2->GetXaxis()->GetBinCenter(i+1);
        for (unsigned int j = 0; j < config.fNy; ++j) {
          const double y = fH2->GetYaxis()->GetBinCenter(j+1);
          for (unsigned int k = 0; k < config.fNz; ++k) {
            const double z = fH2->GetZaxis()->GetBinCenter(k+1);
            Vector3d pos(x/utl::kpc*kpc, y/utl::kpc*kpc, z/utl::kpc*kpc);
            fHI->SetBinContent(i+1, j+1, k+1, HI->get(pos)*dragonToRUQI);
            fHII->SetBinContent(i+1, j+1, k+1, HII->get(pos)*dragonToRUQI);
            fH2->SetBinContent(i+1, j+1, k+1,
                               XCO->get(pos)*H2->get(pos)*dragonToRUQI);
          }
        }
      }
      INFO("done filling gas maps");

      currDir->cd();
    }
    fIsInitialized = true;
  }

  Gas::~Gas()
  {
    delete fH2;
    delete fHI;
    delete fHII;
  }

  double
  Gas::GetH2Density(const TVector3& p)
    const
  {
    if (!fH2) {
      const Vector3d pos(p.X()/utl::kpc*kpc,
                         p.Y()/utl::kpc*kpc,
                         p.Z()/utl::kpc*kpc);
      const double dragonToRUQI = (1/utl::cm3)/(1/cm3);
      const double nH2 = fXCOModel->get(pos)*fH2Model->get(pos)*dragonToRUQI;
      return nH2;
    }
    else
      return Interpolate(*fH2, p);
  }

  double
  Gas::GetHIDensity(const TVector3& p)
    const
  {
    if (!fHI) {
      const Vector3d pos(p.X()/utl::kpc*kpc,
                         p.Y()/utl::kpc*kpc,
                         p.Z()/utl::kpc*kpc);
      const double dragonToRUQI = (1/utl::cm3)/(1/cm3);
      const double nHI = fHIModel->get(pos)*dragonToRUQI;
      return nHI;
    }
    else
      return Interpolate(*fHI, p);
  }

  double
  Gas::GetHIIDensity(const TVector3& p)
    const
  {
    if (!fHII) {
      const Vector3d pos(p.X()/utl::kpc*kpc,
                         p.Y()/utl::kpc*kpc,
                         p.Z()/utl::kpc*kpc);
      const double dragonToRUQI = (1/utl::cm3)/(1/cm3);
      const double nHII = fHIIModel->get(pos)*dragonToRUQI;
      return nHII;
    }
    else
      return Interpolate(*fHII, p);
  }

}
