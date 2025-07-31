#include "Nuclei.h"
#include "PDGParticleIds.h"

#include <utl/PhysicalConstants.h>
#include <utl/Units.h>
#include <utl/DebugOutput.h>

#include <TRandom3.h>

#include <iostream>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <string>

using namespace utl;
using namespace std;

namespace galmc {

  Nuclei::ESigmaInel Nuclei::fSigmaInelOption = Nuclei::eUninitialized;
  Nuclei::ESigmaFrag Nuclei::fSigmaFragOption = Nuclei::eOpt12;

  namespace {
    inline
    double
    SafeEval(const double xVal, const TGraph& g)
    {
      const int n = g.GetN();
      if (n <= 0)
        return 0;
      const double xMin = *g.GetX();
      const double xMax = *(g.GetX() + n - 1);
      const double xx = std::min(std::max(xMin, xVal), xMax);
      return g.Eval(xx);
    }
  }

  void
  Nuclei::Init(const ESigmaInel optInel, const ESigmaFrag optFrag)
  {
    if (fSigmaInelOption != eUninitialized)
      FATAL("reinitialization not allowed");
    fSigmaInelOption = optInel;
    fSigmaFragOption = optFrag;
  }

  const Nuclei&
  Nuclei::GetInstance()
  {
    if (fSigmaInelOption == eUninitialized)
      FATAL("need to call Init() before GetInstance()");
    static Nuclei db;
    return db;
  }


  Nuclei::Nuclei()
  {
    const string datadir = string(GMCPATH) + "/DataFiles";
    ReadNuclearProperties(datadir);
    ReadInelasticCrossSection(datadir);
    ReadFragmentationCrossSection(datadir);
  }


  void
  Nuclei::ReadFragmentationCrossSection(const string& datadir)
  {
    const double fudgeFactorLi = 1.15;
    const string opt = fSigmaFragOption <= 12 ? to_string(12) :
      to_string(fSigmaFragOption);
    const string filename = datadir + "/sigProdGALPROP17_OPT" +
      opt + ".dat";
    INFO("fragmentation option: " + to_string(fSigmaFragOption));
    ifstream fragfile(filename);
    if (!fragfile.good())
      FATAL("Can not opening fragmentation file");

    string line;
    const string comment = "#";
    while (getline(fragfile, line)) {
      if (line[0] == '#' || line.empty())
        continue;
      string idIn, idOut, dummy;
      istringstream data(line);
      data >> idIn >> dummy >> idOut;
      if (!data.eof() && !data.good())
        FATAL("error 1 reading nuclear data " + line);
      if (dummy != "->")
        FATAL("error 2 reading nuclear data " + line);
      const int pdgIn = ParticleConst::GetNucleusId(idIn);
      const int pdgOut = ParticleConst::GetNucleusId(idOut);
      vector<double> sigmaP;
      vector<double> sigmaHe;
      vector<double> lgE;
      for (unsigned int i = 0; i < 41; ++i) {
        getline(fragfile, line);
        if (!fragfile.good())
          FATAL("unexpected end of file");
        istringstream data(line);
        double sP, sHe;
        data >> sP >> sHe;
        if (!data.eof() && !data.good())
          FATAL("error 3 reading nuclear data " + line);
        const double fudge =
          fSigmaFragOption == eOpt12LiFudge &&
          ParticleConst::GetNucleusCharge(pdgOut) == 3 &&
          (ParticleConst::GetNucleusCharge(pdgIn) == 6 ||
           ParticleConst::GetNucleusCharge(pdgIn) == 8)
          ? fudgeFactorLi : 1;
        sigmaP.push_back(sP*millibarn*fudge);
        sigmaHe.push_back(sHe*millibarn*fudge);
        const double log10E = -2+i/40.*3;
        lgE.push_back(log10E);
      }
      fFragmentationCrossSectionProtonTarget[pdgIn][pdgOut] =
        TGraph(sigmaP.size(), &lgE.front(), &sigmaP.front());
      fFragmentationCrossSectionHeliumTarget[pdgIn][pdgOut] =
        TGraph(sigmaP.size(), &lgE.front(), &sigmaHe.front());
      fFragmentationCrossSectionProtonTarget[pdgIn][pdgOut].SetBit(TGraph::kIsSortedX);
      fFragmentationCrossSectionHeliumTarget[pdgIn][pdgOut].SetBit(TGraph::kIsSortedX);
    }
  }

  int
  Nuclei::RandomFragment(const int pdgId, const double lgEkPerNuclGeV,
                         const double heProb)
    const
  {
    const bool isHe = gRandom->Uniform() < heProb;
    const auto& fragCrossSections =
      isHe ? fFragmentationCrossSectionHeliumTarget :
      fFragmentationCrossSectionProtonTarget;
    const auto& csFrags = fragCrossSections.at(pdgId);
    vector<double> cumuSigma;
    vector<int> pdgIds;
    double sum = 0;
    for (const auto& cs : csFrags) {
      // p, d, t, 3He, 4He considered ablation products
      if (ParticleConst::GetNucleusCharge(cs.first) > 2) {
        const double s = SafeEval(lgEkPerNuclGeV, cs.second);
        sum += s;
        cumuSigma.push_back(sum);
        pdgIds.push_back(cs.first);
      }
    }

    if (sum <= 0) {
      // ERROR("sum is zero");
      return 0;
    }

    const double sigmaTot = isHe ?
      GetSigmaInel(pdgId, ParticleConst::eHe4, lgEkPerNuclGeV) :
      GetSigmaInel(pdgId, ParticleConst::eProton, lgEkPerNuclGeV);

    double norm = sum;
    if (sum < sigmaTot) {
      cumuSigma.push_back(sigmaTot);
      pdgIds.push_back(1000000010);
      norm = sigmaTot;
    }

    const double r = gRandom->Uniform() * norm;
    auto it = upper_bound(cumuSigma.begin(), cumuSigma.end(), r);
    if (it == cumuSigma.end())
      FATAL("index error");
    const int index = it - cumuSigma.begin();
    return pdgIds[index];
    /*
    for (unsigned int i = 0; i < cumuSigma.size(); ++i) {
      if (cumuSigma[i] >= r) {
        return pdgIds[i];
      }
    }
    return 0;
    */
  }


  double
  Nuclei::AlphaProtonSigmaRatio(const int pdgId)
    const
  {
    const double A = ParticleConst::GetNucleusMassNumber(pdgId);
    // Ferrando et al, Phys.Rev.C 37 (1988) 1490
    return 2.1 / pow(A, 0.055);

    /*
      // Ferrando:
      TF1* a = new TF1("a", "2.1/pow(x,0.055)",4,56);
      // Gaisser:
      TF1* c = new TF1("c", "pow(1.59+pow(x,1/3.)-1.12, 2)/"
      "pow(1+pow(x,1/3.)-1.12, 2)",4,56);
      c->SetLineColor(kBlue);
      a->Draw();
      c->Draw("SAME");
     */

  }

  void
  Nuclei::ReadNuclearProperties(const string& datadir)
  {
    INFO("unstable nuclei:");
    ifstream nucdatafile(datadir + "/crcharts_Zmax30_ghost84.dat");
    string line;
    const string comment = "#";
    while (getline(nucdatafile, line)) {
      if (line[0] == '#' || line.empty())
        continue;
      istringstream data(line);
      double amu, r, t, tErr;
      int A, Z;
      string name, dec, dummy, tunit, decName;
      data >> amu >> Z >> A >> name >> dec >> r >> dummy >> t >> tErr
           >> tunit >> decName;

      if (!data.good())
        FATAL("error reading nuclear data");
      if (t > 0) {
        if (tunit == "s")
          t *= second;
        else if (tunit == "d")
          t *= day;
        else if (tunit == "yr")
          t *= year;
        else if (tunit == "kyr")
          t *= 1e3*year;
        else if (tunit == "Myr")
          t *= megayear;
        else
          FATAL("unknown unit");
      }
      if (Z > 0) {
        const int pdgId = ParticleConst::GetNucleusId(Z, A);
        const int daughterId = t > 0 ? ParticleConst::GetNucleusId(decName) : 0;
        fNuclei.emplace(pdgId, NuclearProperties(amu, A, Z, t, daughterId));
        if (t > 0)
          cout << setw(5) << name << " -->" << setw(5) << decName
               << " tau/Myr = " << t /megayear << endl;
      }
    }
    INFO("read properties of " + to_string(fNuclei.size()) + " nuclei");
  }

  void
  Nuclei::ReadInelasticCrossSection(const string& datadir)
  {
    string xsModel;
    switch (fSigmaInelOption) {
    case eLetaw83:
      xsModel = "Letaw83";
      break;
    case eWellisch96:
      xsModel = "Wellish96";
      break;
    case eBarashenkov01:
      xsModel = "Barashenkov01";
      break;
    case eCRN6:
      xsModel = "CRN6";
      break;
    default:
      FATAL("unknown model");
    }
    INFO("inel. xs model: " + xsModel);

    string line;
    const string comment = "#";

    if (fSigmaInelOption == eCRN6) {
      ifstream xsfile(datadir + "/SigTotCRN6.txt");
      while (getline(xsfile, line)) {
        if (line[0] == '#' || line.empty())
          continue;
        istringstream header(line);
        string name;
        int Z, A;
        header >> name >> Z >> A;
        const int pdgId = ParticleConst::GetNucleusId(Z, A);
        vector<double> lgE;
        vector<double> sigmaP;
        vector<double> sigmaHe;
        for (unsigned int i = 0; i < 97; ++i) {
          getline(xsfile, line);
          if (!xsfile.good())
            FATAL("unexpected end of file");
          istringstream data(line);
          double lE, sP, sHe;
          data >> lE >> sP >> sHe;
          if (data.fail())
            FATAL("error reading XS_inel data");
          sigmaP.push_back(sP*millibarn);
          sigmaHe.push_back(sHe*millibarn);
          lgE.push_back(lE);
        }

        fInelasticCrossSectionP[pdgId] =
          TGraph(sigmaP.size(), &lgE.front(), &sigmaP.front());
        fInelasticCrossSectionP[pdgId].SetBit(TGraph::kIsSortedX);
        fInelasticCrossSectionHe[pdgId] =
          TGraph(sigmaHe.size(), &lgE.front(), &sigmaHe.front());
        fInelasticCrossSectionHe[pdgId].SetBit(TGraph::kIsSortedX);
      }
    }
    else {
      ifstream xsfile(datadir + "/SigTotGAL_" +
                      to_string(fSigmaInelOption) + ".txt");
      while (getline(xsfile, line)) {
        if (line[0] == '#' || line.empty())
          continue;
        const string name = line;
        const int pdgId = ParticleConst::GetNucleusId(name);
        vector<double> lgE;
        vector<double> sigmaP;
        vector<double> sigmaHe;
        for (unsigned int i = 0; i < 41; ++i) {
          getline(xsfile, line);
          if (!xsfile.good())
            FATAL("unexpected end of file");
          istringstream data(line);
          double s;
          data >> s;
          if (data.fail())
            FATAL("error reading XS_inel data");
          const double E = 0.01*pow(1000, i/40.);
          sigmaP.push_back(s*millibarn);
          const double heFac = AlphaProtonSigmaRatio(pdgId);
          sigmaHe.push_back(s*millibarn*heFac);
          lgE.push_back(log10(E));
        }
        fInelasticCrossSectionP[pdgId] =
          TGraph(sigmaP.size(), &lgE.front(), &sigmaP.front());
        fInelasticCrossSectionP[pdgId].SetBit(TGraph::kIsSortedX);
        fInelasticCrossSectionHe[pdgId] =
          TGraph(sigmaHe.size(), &lgE.front(), &sigmaHe.front());
        fInelasticCrossSectionHe[pdgId].SetBit(TGraph::kIsSortedX);
      }
      INFO("read sigma_inel of " + to_string(fInelasticCrossSectionP.size()) +
           " nuclei");
    }
  }

  /*
      const string outfilename(datadir + "/nuclei.txt");
      ofstream out(outfilename);
      out << name << " " << ParticleConst::GetNucleusCharge(pdgId) << " "
          <<  ParticleConst::GetNucleusMassNumber(pdgId) << endl;
      out.close();
      INFO("wrote particle list to " + outfilename);
  */

  double
  Nuclei::GetSigmaInel(const int projectileId,
                       const int targetId,
                       const double lgEkPerNuclGeV)
    const
  {
    if (targetId != ParticleConst::eProton &&
        targetId != ParticleConst::eProtonNucleus &&
        targetId != ParticleConst::eHe4)
      FATAL("target not implemented");
    const auto& cs = targetId == ParticleConst::eHe4 ?
      fInelasticCrossSectionHe :
      fInelasticCrossSectionP;
    const auto iter = cs.find(projectileId);
    if (iter == cs.end())
      FATAL("unknown particle " + to_string(projectileId));
    return SafeEval(lgEkPerNuclGeV, iter->second);
  }


  Nuclei::~Nuclei()
  {}

  bool
  Nuclei::HasNuclearProperties(const int pdgCode)
    const
  {
    return fNuclei.count(pdgCode);
  }


  const NuclearProperties&
  Nuclei::GetNuclearProperties(const int pdgCode)
    const
  {
    const auto iter = fNuclei.find(pdgCode);
    if (iter == fNuclei.end())
      FATAL("unknown pdgCode " + std::to_string(pdgCode));
    return iter->second;
  }

}
