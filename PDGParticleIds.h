#ifndef _PDGParticleIds_h_
#define _PDGParticleIds_h_

#include <string>

namespace galmc {

  namespace ParticleConst {

    const int kNucleusOffset = 1000000000;

    enum EType {
      eUndefined = 0,
      eElectron = 11, ePositron = -11,
      eNuElectron = 12, eAntiNuElectron = -12,
      eMuon = 13, eAntiMuon = -13,
      eNuMuon = 14, eAntiNuMuon = -14,
      eTau = 15, eAntiTau = -15,
      eNuTau = 16, eAntiNuTau = -16,
      ePhoton = 22,
      eRhoDiffZero = 110,
      ePiZero = 111,
      eRhoZero = 113,
      eKaonZeroLong = 130,
      ePiPlus = 211, ePiMinus = -211,
      eRhoPlus = 213, eRhoMinus = -213,
      eEta = 221, eOmega = 223,
      eKaonZeroShort = 310, eKaonZeroStar = 313,
      eKaonPlus = 321, eKaonMinus = -321,
      eNeutron = 2112, eAntiNeutron = -2112,
      eProton = 2212, eAntiProton = -2212,
      eDeltaPlusPlus = 2224, eDeltaPlus = 2214,
      eDeltaZero = 2114, eDeltaMinus = 1114,
      eLambda = 3122, eAntiLambda = -3122,
      eLambdac = 4122,
      // Selected nuclei.
      eProtonNucleus =   1000010010,
      eDeuteron =   1000010020,
      eTriton =     1000010030,
      eHelium3 =    1000020030,
      eAlpha =      1000020040,
      eHelium4 =    eAlpha,
      eBeryllium7 = 1000040070,
      eBeryllium9 = 1000040090,
      eCarbon12   = 1000060120,
      eArgon40 =    1000180400,
      eScandium45 = 1000200450,
      eIron56 =       1000260560,
      eLead208 =       1000822080,
      eHe4 = eAlpha,
      eBe7 = eBeryllium7,
      eBe9 = eBeryllium9,
      eAr40 = eArgon40,
      eSc45 = eScandium45,
      eFe56 = eIron56,
      ePb208 = eLead208
    };

    bool IsNucleus(int pid);
    int GetNucleusCharge(int pid);
    unsigned int GetNucleusMassNumber(int pid);
    int GetNucleusId(unsigned int Z, unsigned int A);
    // 10Be, 1H etc
    int GetNucleusId(const std::string& name);
    std::string GetElementName(const int id);
    std::string  GetName(const int id);

  }
}
#endif
