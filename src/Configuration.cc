#include "Configuration.h"
#include "PDGParticleIds.h"
#include "Nuclei.h"

#include <utl/DebugOutput.h>
#include <utl/Units.h>
#include <utl/PhysicalConstants.h>

#include <boost/program_options.hpp>

#include <sstream>
#include <iostream>
#include <iomanip>

using namespace std;
namespace po = boost::program_options;

namespace galmc {


  inline
  string
  PrintConfig(const po::variables_map& vm) {
    ostringstream ret;
    for (const auto& p : vm) {
      ret << "" << setw(10) <<  p.first;
      ret << " = ";

      bool is_char;
      try {
        boost::any_cast<const char *>(p.second.value());
        is_char = true;
      }
      catch (const boost::bad_any_cast &) {
        is_char = false;
      }
      bool is_str;
      try {
        boost::any_cast<string>(p.second.value());
        is_str = true;
      } catch (const boost::bad_any_cast &) {
        is_str = false;
      }

      if (((boost::any)p.second.value()).type() == typeid(int))
        ret << vm[p.first].as<int>();
      else if (((boost::any)p.second.value()).type() == typeid(unsigned int))
        ret << vm[p.first].as<unsigned int>();
      else if (((boost::any)p.second.value()).type() == typeid(bool))
        ret << vm[p.first].as<bool>();
      else if (((boost::any)p.second.value()).type() == typeid(double))
        ret << vm[p.first].as<double>();
      else if (is_char)
        ret << vm[p.first].as<const char * >();
      else if (is_str) {
        string temp = vm[p.first].as<string>();
        if (temp.size()) {
          ret << temp;
        } else {
          ret << "[empty string]";
        }
      }
      else { // Assumes that the only remainder is vector<double>
        try {
          const auto v = vm[p.first].as<vector<double> >();
          ret << "[";
          for (unsigned int i = 0; i < v.size(); ++i) {
            ret << v[i];
            if (i + 1 < v.size())
              ret << ",";
          }
          ret << "]";
        } catch (const boost::bad_any_cast &) {
          ostringstream warn;
          warn << "UnknownType("
               << ((boost::any)p.second.value()).type().name()
               << ")";
          WARN(warn.str());
        }
      }
      if (((boost::any)p.second.value()).empty())
        ret << " (empty)";
      if (vm[p.first].defaulted() || p.second.defaulted())
        ret << " (default)";
      ret << "\n";
    }
    const string s = ret.str();
    if (!s.empty())
      return s.substr(0, s.size()-1); // strip last "\n"
    else
      return s;
  }

  Configuration::Configuration(const int argc, const char** argv)
  {
    po::options_description options("Options");
    options.add_options()
      ("help,h", "write this message")
      //------------------ file names
      ("config,c", po::value<string>()->default_value("galmc.config"),
       "config file name.")
      ("output,o", po::value<string>(&fOutFilename)->
       default_value("galmc.root"), "output filename.")
      ("gmf-config", po::value<string>(&fMagneticFieldConfig)->
       default_value("JF12b.xml"), "magnetic field configuration filename.")
      //------------------ general run options
      ("seedI", po::value<int>(&fSeedInjection)->default_value(0),
       "random seed for injection.")
      ("seedP", po::value<int>(&fSeedPropagation)->default_value(0),
       "random seed for propagation.")
      ("seedM", po::value<int>(&fSeedMagneticField)->default_value(0),
       "random seed for field realization.")
      ("nthreads", po::value<int>(&fNthreads)
       ->default_value(1), "number of threads.")
      ("npart,n", po::value<unsigned int>(&fNparticles)->default_value(1),
       "number of particles.")
      ("Dconst", po::value<double>(&fDiffusionConstant)->default_value(0),
       "hom.&iso. diffusion if != 0 [cm^2/s]")
      ("delta", po::value<double>(&fDiffusionDelta)->default_value(0),
       "D(R) = Dconst*(R/10GV)^delta")
      //----------------- tracking options
      ("maxStep", po::value<double>(&fMaxStep)->default_value(0.1),
       "maximum step size [kpc].")
      //------------------ particle distribution and type
      ("xInj", po::value<double>(&fXinjection)->default_value(0),
       "injection center point x [kpc].")
      ("yInj", po::value<double>(&fYinjection)->default_value(0),
       "injection center point y [kpc].")
      ("zInj", po::value<double>(&fZinjection)->default_value(0),
       "injection center point z [kpc].")
      ("rMaxInj", po::value<double>(&fRmaxInjection)->default_value(16),
       "max injection radius [kpc].")
      ("rBias", po::value<double>(&fRbiasInjection)->default_value(-1),
       "injection biasing radius [kpc] (unused if < 0).")
      ("alphaBias", po::value<double>(&fAlphaBiasInjection)->default_value(1),
       "fraction of uniformly injectied particles for biasing")
      ("hInj", po::value<double>(&fH0injection)->default_value(0.20),
       "injection scale height [kpc].")
      ("Rmin", po::value<double>(&fLgRmin)->default_value(30),
       "minimum rigidity [GV/c].")
      ("Rmax", po::value<double>(&fLgRmax)->default_value(30),
       "maximum rigidity [GV/c].")
      ("idInj,i", po::value<int>(&fPdgId)
       ->default_value(ParticleConst::eFe56), "PDG id.")
      //------------------ snapshot times
      ("times,t", po::value<vector<double>>(&fTimes)->required(),
       "snap shot times [Myr].")
      //------------------ observer
      ("observerX", po::value<double>(&fObserverX)
       ->default_value(-8.5), "observer X position [kpc]")
      ("observerY", po::value<double>(&fObserverY)
       ->default_value(0), "observer X position [kpc]")
      ("observerZ", po::value<double>(&fObserverZ)
       ->default_value(0), "observer X position [kpc]")
      ("observerR", po::value<double>(&fObserverR2)
       ->default_value(1), "observer radius [kpc]")
      ("observerH", po::value<double>(&fObserverH)
       ->default_value(0.5), "observer height [kpc]")
      //------------------ grid for pre-calculation (gas etc)
      ("preCalcGas", po::value<bool>(&fPreCalcGas)->
       default_value(0), "precalculate gas densities in grid.")
      ("xMin", po::value<double>(&fXmin)->
       default_value(-30), "minimum x [kpc].")
      ("xMax", po::value<double>(&fXmax)->
       default_value(30), "maximum x [kpc].")
      ("yMin", po::value<double>(&fYmin)->
       default_value(-30), "minimum y [kpc].")
      ("yMax", po::value<double>(&fYmax)->
       default_value(30), "maximum y [kpc].")
      ("zMin", po::value<double>(&fZmin)->
       default_value(-10.5), "minimum z [kpc].")
      ("zMax", po::value<double>(&fZmax)->
       default_value(10.5), "maximum z [kpc].")
      ("nX", po::value<unsigned int>(&fNx)->
       default_value(100), "number of x bins.")
      ("nY", po::value<unsigned int>(&fNy)->
       default_value(100), "number of y bins.")
      ("nZ", po::value<unsigned int>(&fNz)->
       default_value(101), "number of z bins.")
      //------------------ misc
      ("bFudge", po::value<double>(&fRandBFudge2)->
       default_value(1), "fudge factor: b --> b/f")
      //------------------ hadronic stuff
      ("sigmaInel", po::value<int>(&fSigmaInelOpt)->
       default_value(Nuclei::eBarashenkov01),
       "0/1/2/3 = Letaw83/Wellisch96/Barashenkov01/CRN6")
      ("sigmaFrag", po::value<int>(&fSigmaFragOpt)->
       default_value(Nuclei::eOpt12),
       "10/12/22 = fudgeLi12/opt12/op22")
      //------------------ diffusion volume
      ("diffZmax", po::value<double>(&fDiffZmax)->
       default_value(1000), "maximum diffusion height [kpc].")
      ("diffRmax", po::value<double>(&fDiffR2max)->
       default_value(1000), "maximum diffusion radius [kpc].")
      ("lMax", po::value<double>(&fLmax)->
       default_value(50), "outer scale of turbulence [pc].")
      ("trackingThresh", po::value<double>(&fTrackingThresh)->
       default_value(-1), "threshold between tracking and diffusion.")
      ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, options), vm);

    if (vm.count("help")) {
      cerr << options << flush;
      exit(EXIT_SUCCESS);
    }

    if (vm.count("config") > 0) {
      ifstream configStream(vm["config"].as<string>().c_str());
      po::store(po::parse_config_file(configStream, options),
                vm);
    }
    try {
      po::notify(vm);
    }
    catch (const exception& e) {
      ERROR("\n\t" + string(e.what()));
      cerr << options << flush;
      exit(EXIT_FAILURE);
    }
    INFO("\n" + PrintConfig(vm));
    for (auto& t : fTimes)
      t *= utl::megayear;
    fXinjection *= utl::kpc;
    fYinjection *= utl::kpc;
    fZinjection *= utl::kpc;
    fRmaxInjection *= utl::kpc;
    fRbiasInjection *= utl::kpc;
    fH0injection *= utl::kpc;
    fObserverX *= utl::kpc;
    fObserverY *= utl::kpc;
    fObserverZ *= utl::kpc;
    fLgRmin *= utl::GeV / utl::kSpeedOfLight / utl::eplus;
    fLgRmax *= utl::GeV / utl::kSpeedOfLight / utl::eplus;
    fLgRmin = log10(fLgRmin);
    fLgRmax = log10(fLgRmax);
    fObserverR2 *= utl::kpc;
    fObserverR2 = pow(fObserverR2, 2);
    fObserverH *= utl::kpc;
    fDiffusionConstant *= (utl::cm2 / utl::s);
    fMaxStep *= utl::kpc;
    fRandBFudge2 = pow(fRandBFudge2, 2);
    if (fSigmaInelOpt < Nuclei::eFirstSigmaInel ||
        fSigmaInelOpt > Nuclei::eLastSigmaInel)
      FATAL("invalid sigma_inel option");
    if (fSigmaFragOpt != Nuclei::eOpt12 &&
        fSigmaFragOpt != Nuclei::eOpt22 &&
        fSigmaFragOpt != Nuclei::eOpt12LiFudge)
      FATAL("invalid sigma_frag option");
    fDiffZmax  *= utl::kpc;
    fDiffR2max = pow(fDiffR2max*utl::kpc, 2);
    fLmax *= utl::pc;

    if (fRbiasInjection > 0) {
      const double r0 = fRbiasInjection;
      const double R = fRmaxInjection;
      const double alpha = fAlphaBiasInjection;
      const double beta =
        pow(R/r0, 2) * (1 - alpha)/alpha / (2*(1-(1+R/r0)*exp(-R/r0)));
      ostringstream info;
      info << " weights between 1 and " << 1+beta;
      INFO(info.str());
    }
  }
}
