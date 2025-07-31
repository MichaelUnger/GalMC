#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE GalMC
#include <boost/test/unit_test.hpp>

#include "Particle.h"
#include "PDGParticleIds.h"
#include "Diffusion.h"
#include "Nuclei.h"

#include <utl/Units.h>
#include <utl/Verbosity.h>

#include <iostream>

using namespace std;
using namespace utl;
using namespace galmc;

struct Setup {
  Setup() {
    Verbosity::SetLevel(Verbosity::eSilent);
    Nuclei::Init(Nuclei::eBarashenkov01, Nuclei::eOpt12);
  }
};

BOOST_GLOBAL_FIXTURE(Setup);

BOOST_AUTO_TEST_SUITE(ParticleTest)
BOOST_AUTO_TEST_CASE(Kinematics)
{
  // particle physics rule of thumb: rL/m = 3.3 m * (p/GeV/c)/(z B/T)
  const double mpc2 = 1.672621923e-27*kg * kSpeedOfLight2;
  const double p = 1*GeV/kSpeedOfLight;
  const double E = sqrt(pow(p*kSpeedOfLight, 2) + pow(mpc2, 2));
  const double Ek = E - mpc2;
  Particle p1(ParticleConst::eProton, Ek, 0, TVector3());
  const double rL1 = p1.GetLarmorRadius(1*tesla);
  BOOST_CHECK_CLOSE(rL1, 3.3*m, 0.1/percent);

  // astro rule of thumb: rL/kpc = 1.1 kpc * (E/EeV)/(z B/muG)
  Particle p2(ParticleConst::eProton, 1*EeV, 0, TVector3());
  const double rL2 = p2.GetLarmorRadius(1*microgauss);
  BOOST_CHECK_CLOSE(rL2, 1.1*kpc, 0.1/percent);

  // velocity
  BOOST_CHECK_CLOSE(p1.GetBeta(), 0.729256, 1e-6/percent);
  BOOST_CHECK_CLOSE(p2.GetBeta(), 1, 1e-6/percent);

}

BOOST_AUTO_TEST_CASE(NuclearData)
{
  const auto& nuclei = Nuclei::GetInstance();
  vector<int> Z = {1, 14, 26};
  vector<int> A = {2, 28, 56};
  vector<double> lgE = {-2, -2+9/40., 2};
  vector<double> sigma = {1.470000e+02*millibarn,
                          7.337035e+02*millibarn,
                          7.472042e+02*millibarn};
  for (unsigned int i = 0; i < Z.size(); ++i) {
    const int id = ParticleConst::GetNucleusId(Z[i], A[i]);
    BOOST_CHECK_CLOSE(nuclei.GetSigmaInel(id, ParticleConst::eProton, lgE[i]),
                      sigma[i], 1e-6/percent);
  }

  const auto& feProp =
    nuclei.GetNuclearProperties(ParticleConst::GetNucleusId(26, 56));
  BOOST_CHECK_EQUAL(feProp.IsStable(), 1);
  const auto& beProp =
    nuclei.GetNuclearProperties(ParticleConst::GetNucleusId(4, 10));
  BOOST_CHECK_EQUAL(beProp.IsStable(), 0);
  BOOST_CHECK_CLOSE(beProp.fBetaDecayHalfLife, 1.387*megayear, 1e-6/percent);

  const bool printout = false;

  for (unsigned int isHelium = 0; isHelium < 2; ++isHelium) {
    const vector<int> pdgIds =
      { ParticleConst::GetNucleusId("12C"),
        ParticleConst::GetNucleusId("10B"),
        ParticleConst::GetNucleusId("28Si"),
        ParticleConst::GetNucleusId("56Fe") };
    for (auto pdgId : pdgIds) {
      const double lgEnergy = 1;
      const double sigmaTot = isHelium ?
        nuclei.GetSigmaInel(pdgId, ParticleConst::eHe4, lgEnergy) :
        nuclei.GetSigmaInel(pdgId, ParticleConst::eProton, lgEnergy);

      const auto& fragments =
        isHelium ? nuclei.GetCSHeliumTarget().at(pdgId) :
        nuclei.GetCSProtonTarget().at(pdgId);
      if (true) {
        ostringstream info;
        info << ParticleConst::GetName(pdgId) << "+" << (isHelium?"He":"p");
        cout << "---- " << info.str() << endl;
      }
      map<int, unsigned int> fragmentCounts;
      const int N = 100000;
      for (unsigned int i = 0; i < N; ++i) {
        const int fragId = nuclei.RandomFragment(pdgId, lgEnergy, isHelium);
        ++fragmentCounts[fragId];
      }
      for (const auto f : fragmentCounts) {
        const int fragId = f.first;
        if (ParticleConst::GetNucleusCharge(fragId) > 2) {
          const double expected =
            fragments.at(f.first).Eval(lgEnergy)/sigmaTot;
          const double observed = f.second / double(N);
          const double expError = sqrt((expected - pow(expected, 2))/N);
          if (printout) {
            cout << fragId << " " << expected << " "
                 << observed << "+/-" << expError << " "
                 << (observed - expected) / expError << endl;
          }
          BOOST_CHECK( std::abs((observed - expected) / expError) < 3);
        }
      }
    }
  }

}

BOOST_AUTO_TEST_CASE(Basics)
{
  Particle p1(ParticleConst::eProton, 1*GeV, 0, TVector3());
  const double protonMass = 1.672621923e-27*kg;
  BOOST_CHECK_CLOSE(p1.GetMass(), protonMass, 1e-6/percent);
  BOOST_CHECK_EQUAL(p1.GetZ(), 1);
  BOOST_CHECK_EQUAL(p1.GetA(), 1);
  Particle p2(ParticleConst::eCarbon12, 1*GeV, 0, TVector3());
  const double c12mass = 12 * kAtomicMass;
  BOOST_CHECK_CLOSE(p2.GetMass(), c12mass, 1e-6/percent);
  BOOST_CHECK_EQUAL(p2.GetZ(), 6);
  BOOST_CHECK_EQUAL(p2.GetA(), 12);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(DiffusionTest)
BOOST_AUTO_TEST_CASE(Coefficient)
{
  Particle p1(ParticleConst::eProton, 1e15*eV, 0, TVector3());
  const double b2 = pow(1*microgauss, 2);
  const double B2 = pow(4*microgauss, 2);
  const Diffusion diff(25*pc);
  const auto D = diff.GetCoefficient(p1, B2, b2);
  BOOST_CHECK_CLOSE(D.first, 3.09238e+30 * (cm2/s), 1e-6/percent);
  BOOST_CHECK_CLOSE(D.second, 8.79281e+26 * (cm2/s), 1e-6/percent);
}

BOOST_AUTO_TEST_CASE(IsotropicDiffusion)
{
  Particle p1(ParticleConst::eCarbon12, 5*12*GeV, 0, TVector3());
  const Diffusion diff(25*pc, 10*cm2/s, 0.5);
  const auto D = diff.GetCoefficient(p1, 0, 0);
  BOOST_CHECK_EQUAL(D.first, D.second);
  cout << p1.GetPByQ() / (GeV/kSpeedOfLight) << endl;
  BOOST_CHECK_CLOSE(D.first/(cm2/s), 10.68964640, 1e-6/percent);
}
BOOST_AUTO_TEST_SUITE_END()

//BOOST_AUTO_TEST_SUITE(GasTest)
//BOOST_AUTO_TEST_CASE(Density)
//{
//}
//BOOST_AUTO_TEST_SUITE_END()
