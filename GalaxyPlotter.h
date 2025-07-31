#ifndef _GalaxyPlotter_h_
#define _GalaxyPlotter_h_

namespace ruqi {
  class MagneticField;
}

namespace galmc {
  class GalaxyPlotter {
  public:
    GalaxyPlotter(const int argc, const char** argv);
  private:
    void FillDandB(const ruqi::MagneticField&);
    void FillGas();
    void FillXS();
    double fXmin;
    double fXmax;
    double fYmin;
    double fYmax;
    double fZmin;
    double fZmax;
    unsigned int fNx;
    unsigned int fNy;
    unsigned int fNz;
  };
}

#endif
