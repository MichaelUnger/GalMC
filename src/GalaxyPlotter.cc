#include "GalaxyPlotter.h"
#include "Diffusion.h"
#include "Particle.h"
#include "Configuration.h"
#include "PDGParticleIds.h"
#include "Gas.h"
#include "Nuclei.h"

#include <ruqi/MagneticField.h>
#include <utl/DebugOutput.h>
#include <utl/Units.h>

#include <TFile.h>
#include <TH3D.h>
#include <TH2D.h>
#include <TRandom3.h>

#include <boost/program_options.hpp>

using namespace std;
using namespace utl;
using namespace ruqi;

namespace po = boost::program_options;

namespace galmc {
  GalaxyPlotter::GalaxyPlotter(const int argc, const char** argv)
  {
    string outputFilename;
    string configFilename;
    po::options_description options("Options");
    options.add_options()
      ("help,h", "write this message")
      ("output,o", po::value<string>(&outputFilename)->
       default_value("galPlot.root"), "output filename.")
      ("config,c", po::value<string>(&configFilename)->
       default_value("JF12b.xml"), "magnetic field configuration filename.")
      ("xMin", po::value<double>(&fXmin)->
       default_value(-20), "minimum x [kpc].")
      ("xMax", po::value<double>(&fXmax)->
       default_value(20), "maximum x [kpc].")
      ("yMin", po::value<double>(&fYmin)->
       default_value(-20), "minimum y [kpc].")
      ("yMax", po::value<double>(&fYmax)->
       default_value(20), "maximum y [kpc].")
      ("zMin", po::value<double>(&fZmin)->
       default_value(-10.5), "minimum z [kpc].")
      ("zMax", po::value<double>(&fZmax)->
       default_value(10.5), "maximum z [kpc].")
      ("nX", po::value<unsigned int>(&fNx)->
       default_value(100), "number of x bins.")
      ("nY", po::value<unsigned int>(&fNy)->
       default_value(100), "number of y bins.")
      ("nZ", po::value<unsigned int>(&fNz)->
       default_value(101), "number of z bins.")
      ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, options), vm);

    if (vm.count("help")) {
      cerr << options << flush;
      exit(EXIT_SUCCESS);
    }

    try {
      po::notify(vm);
    }
    catch (const exception& e) {
      ERROR("\n\t" + string(e.what()));
      cerr << options << flush;
      exit(EXIT_FAILURE);
    }

    MagneticField mf(configFilename);

    TFile outfile(outputFilename.c_str(), "RECREATE");

    FillXS();
    FillDandB(mf);
    FillGas();

    outfile.Write();
    outfile.Close();

  }

  void
  GalaxyPlotter::FillXS()
  {
    INFO("filling XS....");
    const double lgE = 0.999;
    Nuclei::Init(Nuclei::eWellisch96, Nuclei::eOpt12LiFudge);
    const auto& nuclei = Nuclei::GetInstance();

    for (int isHe = 0; isHe < 2; ++isHe) {

      const auto& csInel = isHe ?
        nuclei.GetSigmaInelHeliumTarget() :
        nuclei.GetSigmaInelProtonTarget();
      const double A2 = isHe ? 4 : 1;

      for (const auto iter : csInel) {
        const auto& gxs = iter.second;
        TGraph g;
        const double A1 = ParticleConst::GetNucleusMassNumber(iter.first);
        const double R0 = 1.47*fermi;
        const double delta = 1.12;
        const double sigma0 =
          kPi * pow(R0, 2) * pow(pow(A1, 1/3.) + pow(A2, 1/3.) - delta, 2);

        for (int i = 0; i < gxs.GetN(); ++i) {
          g.SetPoint(i, *(gxs.GetX()+i), *(gxs.GetY()+i) / sigma0);
        }
        g.GetXaxis()->SetTitle("lg((E/n)/GeV)");
        g.GetYaxis()->SetTitle("#sigma/#sigma_{0}");
        g.SetName((string("xs") + (isHe ? "He" : "P") +
                   to_string(iter.first)).c_str());
        g.Write();
      }
    }


    const auto& frag = nuclei.GetCSProtonTarget();
    TH2D hBreakUpProbP("hBreakUpProbP",
                       ";number of neutrons;number of protons;P(any frag>He)",
                       30, 0.5, 30.5, 30, 0.5, 30.5);

    TH2D hFragN("hFragN", ";number of neutrons;number of protons;"
                "number of fragmentation channels",
                30, 0.5, 30.5, 30, 0.5, 30.5);
    TH2D hFragB("hFragB", ";number of neutrons;number of protons;"
                "#sigma(A+p#rightarrow B + X) / mb",
                30, 0.5, 30.5, 30, 0.5, 30.5);
    TH2D hFragBe("hFragBe", ";number of neutrons;number of protons;"
                 "#sigma(A+p#rightarrow Be + X) / mb",
                 30, 0.5, 30.5, 30, 0.5, 30.5);
    TH2D hFragLi("hFragLi", ";number of neutrons;number of protons;"
                 "#sigma(A+p#rightarrow Li + X) / mb",
                 30, 0.5, 30.5, 30, 0.5, 30.5);

    for (int i = 0; i < hFragN.GetNbinsX(); ++i) {
      const int nN = hFragN.GetXaxis()->GetBinCenter(i+1);
      for (int j = 0; j < hFragN.GetNbinsY(); ++j) {
        const int nP = hFragN.GetYaxis()->GetBinCenter(j+1);
        const int pdgId = ParticleConst::GetNucleusId(nP, nP+nN);
        if (frag.count(pdgId)) {
          const auto& xs = frag.at(pdgId);
          hFragN.SetBinContent(i+1, j+1, xs.size());
          double sum = 0;
          for (const auto iter : xs) {
            const int id = iter.first;
            if (ParticleConst::GetNucleusCharge(id) > 2)
              sum += iter.second.Eval(lgE);
            if (ParticleConst::GetNucleusCharge(id) == 5) {
              const double sigma = iter.second.Eval(lgE);
              const double c = hFragB.GetBinContent(i+1, j+1);
              hFragB.SetBinContent(i+1, j+1, c + sigma/millibarn);
            }
            if (ParticleConst::GetNucleusCharge(id) == 4) {
              const double sigma = iter.second.Eval(lgE);
              const double c = hFragBe.GetBinContent(i+1, j+1);
              hFragBe.SetBinContent(i+1, j+1, c + sigma/millibarn);
            }
            if (ParticleConst::GetNucleusCharge(id) == 3) {
              const double sigma = iter.second.Eval(lgE);
              const double c = hFragLi.GetBinContent(i+1, j+1);
              hFragLi.SetBinContent(i+1, j+1, c + sigma/millibarn);
            }
          }
          const double sig0 = nuclei.GetSigmaInelProtonTarget().at(pdgId).Eval(lgE);
          hBreakUpProbP.SetBinContent(i+1, j+1, sum/sig0);
        }
      }
    }
    hBreakUpProbP.Write();
    hFragB.Write();
    hFragBe.Write();
    hFragLi.Write();
    hFragN.Write();
  }

  void
  GalaxyPlotter::FillGas()
  {
    INFO("filling gas....");
    string options = "galaxyPlotter";
    const char* argv = options.c_str();
    Gas::Init(Configuration(1, &argv));

    TH3D hIHist("hI", ";x/kpc;y/kpc;z/kpc",
                  fNx, fXmin, fXmax,
                  fNy, fYmin, fYmax,
                  fNz, fZmin, fZmax);
    TH3D hIIHist("hII", ";x/kpc;y/kpc;z/kpc",
                  fNx, fXmin, fXmax,
                  fNy, fYmin, fYmax,
                  fNz, fZmin, fZmax);
    TH3D h2Hist("h2", ";x/kpc;y/kpc;z/kpc",
                  fNx, fXmin, fXmax,
                  fNy, fYmin, fYmax,
                  fNz, fZmin, fZmax);
    const auto& gas = Gas::GetInstance();
    for (unsigned int iX = 0; iX < fNx; ++iX) {
      const double x = hIHist.GetXaxis()->GetBinCenter(iX+1) * kpc;
      for (unsigned int iY = 0; iY < fNy; ++iY) {
        const double y = hIHist.GetYaxis()->GetBinCenter(iY+1) * kpc;
        for (unsigned int iZ = 0; iZ < fNz; ++iZ) {
          const double z = hIHist.GetZaxis()->GetBinCenter(iZ+1) * kpc;
          const TVector3 pos(x, y, z);
          h2Hist.SetBinContent(iX+1, iY+1, iZ+1,
                               gas.GetH2Density(pos)*2/(1/cm3));
          hIHist.SetBinContent(iX+1, iY+1, iZ+1,
                               gas.GetHIDensity(pos)/(1/cm3));
          hIIHist.SetBinContent(iX+1, iY+1, iZ+1,
                                gas.GetHIIDensity(pos)/(1/cm3));
        }
      }
    }
    hIHist.Write();
    hIIHist.Write();
    h2Hist.Write();

    // 1D plots
    vector<double> z, nZ, nZ2, nZI, nZII;
    for (int iZ = 0; iZ < 500; ++iZ) {
      const TVector3 pos(-8.5*kpc, 0, iZ/100.*kpc);
      const auto rho = gas.GetNumberDensities(pos);
      z.push_back(pos.Z()/kpc);
      nZ.push_back(rho.first/(1/cm3));
      nZ2.push_back(2*gas.GetH2Density(pos)/(1/cm3));
      nZI.push_back(gas.GetHIDensity(pos)/(1/cm3));
      nZII.push_back(gas.GetHIIDensity(pos)/(1/cm3));
    }
    TGraph gZ(z.size(), &z.front(), &nZ.front());
    gZ.SetName("nVsZ");
    gZ.SetTitle(";z/kpc;n/(1/cm^{3})");
    gZ.Write();
    TGraph gZI(z.size(), &z.front(), &nZI.front());
    gZI.SetName("nIVsZ");
    gZI.SetTitle(";z/kpc;n(HI)/(1/cm^{3})");
    gZI.Write();
    TGraph gZII(z.size(), &z.front(), &nZII.front());
    gZII.SetName("nIIVsZ");
    gZII.SetTitle(";z/kpc;n(HII)/(1/cm^{3})");
    gZII.Write();
    TGraph gZ2(z.size(), &z.front(), &nZ2.front());
    gZ2.SetName("n2VsZ");
    gZ2.SetTitle(";z/kpc;n(H2)/(1/cm^{3})");
    gZ2.Write();

    vector<double> r, nR, nR2, nRI, nRII;
    for (int iR = 0; iR < 200; ++iR) {
      const TVector3 pos(iR/10.*kpc, 0, 0);
      const auto rho = gas.GetNumberDensities(pos);
      r.push_back(pos.X()/kpc);
      nR.push_back(rho.first/(1/cm3));
      nR2.push_back(2*gas.GetH2Density(pos)/(1/cm3));
      nRI.push_back(gas.GetHIDensity(pos)/(1/cm3));
      nRII.push_back(gas.GetHIIDensity(pos)/(1/cm3));
    }
    TGraph gR(r.size(), &r.front(), &nR.front());
    gR.SetName("nVsR");
    gR.SetTitle(";r/kpc;n/(g/cm^{3})");
    gR.Write();
    TGraph gRI(r.size(), &r.front(), &nRI.front());
    gRI.SetName("nIVsR");
    gRI.SetTitle(";r/kpc;n(HI)/(g/cm^{3})");
    gRI.Write();
    TGraph gRII(r.size(), &r.front(), &nRII.front());
    gRII.SetName("nIIVsR");
    gRII.SetTitle(";r/kpc;n(HII)/(g/cm^{3})");
    gRII.Write();
    TGraph gR2(r.size(), &r.front(), &nR2.front());
    gR2.SetName("n2VsR");
    gR2.SetTitle(";r/kpc;n(H2)/(g/cm^{3})");
    gR2.Write();
  }

  void
  GalaxyPlotter::FillDandB(const MagneticField& mf)
  {
    TH3D bXYHist("bXY", ";x/kpc;y/kpc;z/kpc",
                 fNx, fXmin, fXmax,
                 fNy, fYmin, fYmax,
                 fNz, fZmin, fZmax);
    TH3D bZHist("bZ", ";x/kpc;y/kpc;z/kpc",
                fNx, fXmin, fXmax,
                fNy, fYmin, fYmax,
                fNz, fZmin, fZmax);
    TH3D bRand("bRand", ";x/kpc;y/kpc;z/kpc",
               fNx, fXmin, fXmax,
               fNy, fYmin, fYmax,
               fNz, fZmin, fZmax);
    TH3D eta("eta", ";x/kpc;y/kpc;z/kpc",
             fNx, fXmin, fXmax,
             fNy, fYmin, fYmax,
             fNz, fZmin, fZmax);
    TH3D dPerpHist("dPerp", ";x/kpc;y/kpc;z/kpc",
                   fNx, fXmin, fXmax,
                   fNy, fYmin, fYmax,
                   fNz, fZmin, fZmax);
    TH3D dParHist("dPar", ";x/kpc;y/kpc;z/kpc",
                  fNx, fXmin, fXmax,
                  fNy, fYmin, fYmax,
                  fNz, fZmin, fZmax);
    TH3D dZHist("dZ", ";x/kpc;y/kpc;z/kpc",
                fNx, fXmin, fXmax,
                fNy, fYmin, fYmax,
                fNz, fZmin, fZmax);
    TH3D dZMCHist("dZMC", ";x/kpc;y/kpc;z/kpc",
                  fNx, fXmin, fXmax,
                  fNy, fYmin, fYmax,
                  fNz, fZmin, fZmax);

    INFO("filling D and B....");

    const double minField2 = pow(1e-6*microgauss, 2);
    TRandom3 rand(0);

    const Diffusion diff(100*pc);
    for (unsigned int iX = 0; iX < fNx; ++iX) {
      const double x = dPerpHist.GetXaxis()->GetBinCenter(iX+1) * kpc;
      for (unsigned int iY = 0; iY < fNy; ++iY) {
        const double y = dPerpHist.GetYaxis()->GetBinCenter(iY+1) * kpc;
        for (unsigned int iZ = 0; iZ < fNz; ++iZ) {
          const double z = dPerpHist.GetZaxis()->GetBinCenter(iZ+1) * kpc;
          const TVector3 pos(x, y, z);
          const auto bFields = mf.GetFields(x, y, z);
          const double b2 = std::get<MagneticField::eRandomVariance>(bFields);
          const auto& Bvec = std::get<MagneticField::eCoherent>(bFields);
          const double B2 = Bvec.SquaredLength();
          if (b2 < minField2 && B2 < minField2)
            continue;
          const double bXY = sqrt(pow(Bvec.x, 2) + pow(Bvec.y, 2));
          bXYHist.SetBinContent(iX+1, iY+1, iZ+1, bXY / microgauss);
          bZHist.SetBinContent(iX+1, iY+1, iZ+1, Bvec.z / microgauss);
          bRand.SetBinContent(iX+1, iY+1, iZ+1, sqrt(b2) / microgauss);
          eta.SetBinContent(iX+1, iY+1, iZ+1, b2 / (b2 + B2));
          const TVector3 B(Bvec.x, Bvec.y, Bvec.z);
          const Particle p(ParticleConst::eCarbon12, 20*6*GeV, 0, pos);
          const auto D = diff.GetCoefficient(p, B2, b2);
          dPerpHist.SetBinContent(iX+1, iY+1, iZ+1, D.second/(1e28*cm2/s));
          dParHist.SetBinContent(iX+1, iY+1, iZ+1, D.first/(1e28*cm2/s));
          if (B2 < minField2) {
            dZHist.SetBinContent(iX+1, iY+1, iZ+1, D.second/(1e28*cm2/s));
            dZMCHist.SetBinContent(iX+1, iY+1, iZ+1, D.second/(1e28*cm2/s));
          }
          else {
            const auto uB = B.Unit();
            const auto u2 = uB.Orthogonal();
            const auto u1 = uB.Cross(u2);
            // const auto Dvec = uB*D.first;// + u1*D.second + u2*D.second;
            const double step = sqrt(2*D.first);
            const double stepPerp = sqrt(2*D.second);
            double sumZ = 0;
            const int n = 1;
            for (int i = 0; i < n; ++i) {
              const auto delta =
                uB*rand.Gaus(0, step) +
                u1*rand.Gaus(0, stepPerp) +
                u2*rand.Gaus(0, stepPerp);
              sumZ += pow(delta.Z(), 2) / 2;
            }
            dZMCHist.SetBinContent(iX+1, iY+1, iZ+1, sumZ/n/(1e28*cm2/s));
            const double ct = abs(uB.Z());
            const double ct2 = ct*ct;
            const double st2 = 1 - ct2;
            const double Dz = D.first * ct2 + D.second * st2;
            dZHist.SetBinContent(iX+1, iY+1, iZ+1, Dz/(1e28*cm2/s));
          }
        }
      }
    }
    bXYHist.Write();
    bZHist.Write();
    bRand.Write();
    eta.Write();
    dPerpHist.Write();
    dParHist.Write();
    dZHist.Write();
    dZMCHist.Write();
  }
}
