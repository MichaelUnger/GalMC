#ifndef _utl_NuclearProperties_h_
#define _utl_NuclearProperties_h_

namespace galmc {
  class NuclearProperties {
  public:
    NuclearProperties(const double amu = 0, const int Z = 0, const int A = 0,
                      const double tau = 0, const int pdg = 0) :
      fMassAmu(amu), fZ(Z), fA(A), fBetaDecayHalfLife(tau),  fBetaDecayPdgCode(pdg) {}
    bool IsStable() const { return fBetaDecayHalfLife <= 0; }
    const double fMassAmu;
    const int fZ;
    const int fA;
    const double fBetaDecayHalfLife;
    const int fBetaDecayPdgCode;
  };
}

#endif
