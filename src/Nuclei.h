#ifndef _utl_Nuclei_h_
#define _utl_Nuclei_h_

#include "NuclearProperties.h"
#include <map>
#include <TGraph.h>

namespace galmc {
  class Nuclei {
  public:
    enum ESigmaFrag {
      eOpt12LiFudge = 10,
      eOpt12 = 12,
      eOpt22 = 22
    };

    enum ESigmaInel {
      eUninitialized = -1,
      eLetaw83 = 0,
      eWellisch96 = 1,
      eBarashenkov01 = 2,
      eCRN6 = 3,
      eFirstSigmaInel = eLetaw83,
      eLastSigmaInel = eCRN6
    };

  public:
    static void Init(const ESigmaInel optInel, const ESigmaFrag optFrag);
    static const Nuclei& GetInstance();

    /// Test whether a nucleus is in the database
    bool HasNuclearProperties(const int pdgCode) const;
    const NuclearProperties& GetNuclearProperties(const int pdgCode) const;
    double GetSigmaInel(const int projectileId,
                        const int targetId,
                        const double lgEkPerNuclGeV) const;
    const std::map<int, TGraph>& GetSigmaInelProtonTarget() const
    { return fInelasticCrossSectionP; }
    const std::map<int, TGraph>& GetSigmaInelHeliumTarget() const
    { return fInelasticCrossSectionHe; }

    /// dice a random fragment
    int RandomFragment(const int pdgId, const double lgEkPerNucGeV,
                       const double heProb) const;

    using FragmentationMap = std::map<int, std::map<int, TGraph>>;
    const FragmentationMap& GetCSProtonTarget() const
    { return fFragmentationCrossSectionProtonTarget; }
    const FragmentationMap& GetCSHeliumTarget() const
    { return fFragmentationCrossSectionHeliumTarget; }


  protected:
    Nuclei();
    ~Nuclei();

  private:
    double AlphaProtonSigmaRatio(const int pdgId) const;
    void ReadFragmentationCrossSection(const std::string& datadir);
    void ReadInelasticCrossSection(const std::string& datadir);
    void ReadNuclearProperties(const std::string& datadir);
    std::map<int, NuclearProperties> fNuclei;
    std::map<int, TGraph> fInelasticCrossSectionP;
    std::map<int, TGraph> fInelasticCrossSectionHe;
    FragmentationMap fFragmentationCrossSectionProtonTarget;
    FragmentationMap fFragmentationCrossSectionHeliumTarget;
    static ESigmaInel fSigmaInelOption;
    static ESigmaFrag fSigmaFragOption;
  };
}

#endif
