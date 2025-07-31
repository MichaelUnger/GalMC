#include "Particle.h"
#include "DatabasePDG.h"
#include "PDGParticleIds.h"
#include "Nuclei.h"

#include <utl/Units.h>
#include <utl/DebugOutput.h>
#include <TParticlePDG.h>

ClassImp(galmc::Particle)

namespace galmc {

  Particle::Particle(const int pdgId, const double Ek,
                     const double t, const TVector3& pos,
                     const double weight, const int id) :
    fInjectedPdgId(pdgId),
    fInjectedTime(t),
    fInjectedEk(Ek),
    fInjectedPosition(pos),
    fWeight(weight),
    fId(id),
    fPdgId(pdgId),
    fTime(t),
    fPosition(pos)
  {
    SetMassAndCharge();
    UpdateEnergy(Ek);
  }

  void
  Particle::SetMassAndCharge()
  {
    // set mass, charge etc
    fIsNucleus = ParticleConst::IsNucleus(fPdgId);
    if (fIsNucleus) {
      const auto& nuclei = Nuclei::GetInstance();
      const auto& properties = nuclei.GetNuclearProperties(fPdgId);
      fMass =  properties.fMassAmu * utl::kAtomicMass;
      fZ = ParticleConst::GetNucleusCharge(fPdgId);
      fA = ParticleConst::GetNucleusMassNumber(fPdgId);
    }
    else {
      const auto& pdgParticle = DatabasePDG::GetInstance().GetParticle(fPdgId);
      fMass =  pdgParticle.Mass() * utl::GeV/utl::kSpeedOfLight2;
      fZ = pdgParticle.Charge() / 3;
      if (fPdgId == ParticleConst::eProton ||
          fPdgId == ParticleConst::eAntiProton ||
          fPdgId == ParticleConst::eNeutron ||
          fPdgId == ParticleConst::eAntiNeutron)
        fA = 1;
      else
        fA = 0;
    }
  }

  void
  Particle::UpdateEnergy(const double Ek)
  {
    fEk = Ek;
    const double mc2 = fMass * utl::kSpeedOfLight2;
    fE = Ek + mc2;
    fP = sqrt(fE*fE - mc2*mc2) / utl::kSpeedOfLight;
    fPByQ = fP / (std::abs(fZ)*utl::eplus);
  }

  double
  Particle::GetInjectedPByQ()
    const
  {
    if (!fIsInternalUnits)
      FATAL("not in internal units");
    const double Ek = fInjectedEk;
    const auto& nuclei = Nuclei::GetInstance();
    const auto& properties = nuclei.GetNuclearProperties(fInjectedPdgId);
    const double mass =  properties.fMassAmu * utl::kAtomicMass;
    const double mc2 = mass * utl::kSpeedOfLight2;
    const double E = Ek + mc2;
    const double p = sqrt(E*E - mc2*mc2) / utl::kSpeedOfLight;
    const double Z = ParticleConst::GetNucleusCharge(fInjectedPdgId);
    const double pByQ = p / (Z*utl::eplus);
    return pByQ;
  }

  void
  Particle::UpdateParticleId(const int id)
  {
    const double oldA = fA;
    fParentId = fPdgId;
    fPdgId = id;
    SetMassAndCharge();
    const double newA = fA;
    if (oldA < newA)
      FATAL("oldA < newA???");
#warning Ek or E??
    const double Ek = fEk / oldA * newA;
    UpdateEnergy(Ek);
  }

   void
  Particle::ConvertUnits(const bool internalToExternal)
  {
    using namespace utl;
    const double energyConv = internalToExternal ? 1 / GeV : GeV;
    const double GeVc = GeV/utl::kSpeedOfLight;
    const double momentumConv = internalToExternal ? 1 / GeVc : GeVc;
    const double lengthConv = internalToExternal ? 1 / kpc : kpc;
    const double timeConv = internalToExternal ? 1 / megayear : megayear;
    const double cpuTimeConv = internalToExternal ? 1 / second : second;
    const double diffConv = internalToExternal ? 1 / (cm2/s) : (cm2/s);
    const double xConf = internalToExternal ? 1 / (g/cm2) : (g/cm2);
    if (internalToExternal) {
      if (!fIsInternalUnits)
        return;
      fIsInternalUnits = false;
    }
    else {
      if (fIsInternalUnits)
        return;
      fIsInternalUnits = true;
    }

    fInjectedPosition *= lengthConv;
    fInjectedTime *= timeConv;
    fInjectedEk *= energyConv;
    fPosition  *= lengthConv;
    fTime *= timeConv;
    fDpar *= diffConv;
    fDperp *= diffConv;
    fX *= xConf;

    fEk *= energyConv;
    fE *= energyConv;
    fP *= momentumConv;

    fPathLength *= lengthConv;
    fPByQ *= momentumConv;

    fCPUTime *= cpuTimeConv;

  }

}
