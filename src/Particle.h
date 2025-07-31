#ifndef _Particle_h_
#define _Particle_h_
#include <Rtypes.h>

#include <TVector3.h>

#include <utl/PhysicalConstants.h>

namespace galmc {

  class
  Particle {

  public:
    Particle() = default;
    Particle(const int pdgId, const double Ek,
             const double t, const TVector3& pos, const double weight = 1,
             const int id = 0);

    void UpdateEnergy(const double Ek);
    void UpdateParticleId(const int pdgId);

    void SetWeight(const double w) { fWeight = w; }
    double GetWeight() const { return fWeight; }

    int GetInjectedId() const { return fInjectedPdgId; }
    TVector3 GetInjectedPosition() const { return fInjectedPosition; }
    double GetInjectedEk() const { return fInjectedEk; }
    double GetInjectedPByQ() const;
    double GetInjectedTime() const { return fInjectedTime; }

    const TVector3& GetPosition() const { return fPosition; }
    const TVector3& GetDirection() const { return fDirection; }
    void SetPosition(const TVector3& p) { fPosition = p; }
    void SetDirection(const TVector3& d) { fDirection = d; }

    double GetMomentum() const { return fP; }
    double GetEk() const { return fEk; }
    double GetEnergy() const { return fE; }
    double GetGamma() const { return fE / (fMass * utl::kSpeedOfLight2); }
    double GetKineticEnergy() const { return fEk; }
    double GetBeta() const { return fP / fE * utl::kSpeedOfLight; }
    double GetVelocity() const { return fP / fE * utl::kSpeedOfLight2; }
    double GetPByQ() const { return fPByQ; }
    double GetLarmorRadius(const double B) const { return fPByQ/B; }
    int GetId() const { return fPdgId; }
    int GetParentId() const { return fParentId; }
    int GetZ() const { return fZ; }
    int GetA() const { return fA; }
    double GetMass() const { return fMass; }


    double& GetX() { return fX; }
    double GetX() const { return fX; }
    double& GetTime() { return fTime; }
    double GetTime() const { return fTime; }
    TVector3& GetPosition() { return fPosition; }
    TVector3& GetDirection() { return fDirection; }
    double& GetPathLength() { return fPathLength; }
    double GetPathLength() const { return fPathLength; }

    void SetCurrentDiff(const double par, const double perp)
    { fDpar = par; fDperp = perp; }
    std::pair<double, double> GetCurrentDiff()
    { return std::pair<double, double>(fDpar, fDperp); }

    void SetInternalId(const int id) { fId = id; }
    int GetInternalId() const { return fId; }

    void IncrementSaveCounter() { ++fSaveCounter; }
    void SetSaveCounter(const int c) { fSaveCounter = c; }

    void SetCPUTime(const double t) { fCPUTime = t; }

    // convert internal to "human-readable" units and vice versa
    void ConvertUnits(const bool internalToExternal = true);

  private:
    void SetMassAndCharge();

    // injected properties
    int fInjectedPdgId = 0;
    double fInjectedTime = 0;
    double fInjectedEk = 0;
    TVector3 fInjectedPosition;
    double fWeight = 1;

    // internal serial number of particle during this run
    int fId = 0;
    // how often the particle state was saved
    int fSaveCounter = 0;

    // parent particle (0 for primary)
    int fParentId = 0;
    // PDG particle id
    int fPdgId = 0;
    // current time
    double fTime = 0;
    // current position
    TVector3 fPosition;
    // current direction
    TVector3 fDirection;
    // kinetic energy
    double fEk = 0;

    // derived quantities
    double fMc2 = 0;
    double fMass = 0;
    int fZ = 0;
    int fA = 0;
    double fE = 0;
    double fP = 0;
    double fPByQ = 0;
    bool fIsNucleus = false;

    // current diffusion coefficients
    double fDpar = 0;
    double fDperp = 0;

    // current traversed matter density
    double fX = 0;
    // current traversed path length
    double fPathLength = 0;

    // computing time used for this particle (only useful for nThreads=1)
    double fCPUTime;

    // switch from internal to human-readable units
    bool fIsInternalUnits = true;

    ClassDefNV(Particle, 9);
  };

}

#endif
