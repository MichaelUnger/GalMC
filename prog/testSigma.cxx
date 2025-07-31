#include "Nuclei.h"
#include <utl/Units.h>
#include "PDGParticleIds.h"
#include "Configuration.h"

#include <iostream>
#include <cmath>

using namespace utl;
using namespace std;
using namespace galmc;

int
main(const int argc, const char** argv)
{
  const Configuration config(argc, argv);
  Nuclei::Init(Nuclei::ESigmaInel(config.fSigmaInelOpt),
               Nuclei::ESigmaFrag(config.fSigmaFragOpt));

  const Nuclei& nuclei = Nuclei::GetInstance();

  const vector<int> projectileIds =
    { ParticleConst::GetNucleusId(6, 12), ParticleConst::GetNucleusId(8, 16) };
  const vector<vector<int>> fragments =
    {
     {ParticleConst::GetNucleusId(3, 6), ParticleConst::GetNucleusId(3, 7)},
     {ParticleConst::GetNucleusId(4, 7), ParticleConst::GetNucleusId(4, 9),
      ParticleConst::GetNucleusId(4, 10)},
     {ParticleConst::GetNucleusId(5, 10), ParticleConst::GetNucleusId(5, 11)}
    };

  const double EkPerN = 10;
  const double lgE = log10(EkPerN);
  for (const auto frag : fragments) {
    for (const auto f : frag) {
      cout << "sigmaInel(" << ParticleConst::GetName(f) << "+p) = "
           << nuclei.GetSigmaInel(f, ParticleConst::eProton, lgE) / millibarn
           << " mb " << endl;
    }
  }

  for (const auto proj : projectileIds) {
    const auto& fmap = nuclei.GetCSProtonTarget().at(proj);
    for (const auto frag : fragments) {
      cout << "sigma(" << ParticleConst::GetName(proj) << "+p->"
           << ParticleConst::GetElementName(frag.front()) << ") = (";
      double sum = 0;
      for (const auto f : frag) {
        double val = fmap.at(f).Eval(lgE);
        cout << val / millibarn << "+";
        sum += val;
      }
      cout << "\b) mb = " << sum/millibarn << " mb" << endl;
    }
  }

  return 1;
}
