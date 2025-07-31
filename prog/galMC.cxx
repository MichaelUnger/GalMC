#include "Configuration.h"
#include "GalMC.h"

using namespace galmc;

int
main(const int argc, const char** argv)
{
  Configuration c(argc, argv);

  GalMC galMC(c);
  galMC.Run();

  return 0;
}
