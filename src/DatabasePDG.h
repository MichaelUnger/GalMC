#ifndef _utl_DatabasePDG_h_
#define _utl_DatabasePDG_h_

class TDatabasePDG;
class TParticlePDG;

namespace galmc {
  class DatabasePDG {
  public:
    static DatabasePDG& GetInstance();

    /// Test whether a particle is in the database
    bool HasParticle(const int pdgCode) const;
    /// Returns particle from database or placeholder if none is found
    TParticlePDG& GetParticle(const int pdgCode) const;

  protected:
    DatabasePDG();
    ~DatabasePDG();
    TDatabasePDG* const fDB;
  };
}

#endif
