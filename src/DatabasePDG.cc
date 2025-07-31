#include "DatabasePDG.h"
#include "PDGParticleIds.h"

#include <utl/PhysicalConstants.h>
#include <utl/Units.h>
#include <utl/DebugOutput.h>

#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include <iostream>
#include <sstream>

using namespace utl;

namespace galmc {

  DatabasePDG&
  DatabasePDG::GetInstance()
  {
    static DatabasePDG db;
    return db;
  }


  DatabasePDG::DatabasePDG() :
    fDB(TDatabasePDG::Instance())
  {
  }


  DatabasePDG::~DatabasePDG()
  {}


  bool
  DatabasePDG::HasParticle(const int pdgCode)
    const
  {
    TParticlePDG* p = fDB->GetParticle(pdgCode);
    return p != 0;
  }


  TParticlePDG&
  DatabasePDG::GetParticle(const int pdgCode)
    const
  {
    TParticlePDG* p = fDB->GetParticle(pdgCode);

    if (!p) {
      FATAL("unknown particle");
    }
    return *p;
  }
}
