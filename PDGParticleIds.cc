#include "PDGParticleIds.h"

#include <utl/DebugOutput.h>

#include <map>
#include <iostream>

using namespace std;

namespace galmc {

  namespace ParticleConst {

    bool
    IsNucleus(int pid)
    {
      if (pid > kNucleusOffset)
        return true;
      else
        return false;
    }

    int
    GetNucleusCharge(int pid)
    {
      if (!IsNucleus(pid)) {
        if (pid == 1000010010)
          return 1;
        else if (pid == 1000000010)
          return 0;
        else
          throw std::runtime_error("The particle is not a nucleus." + to_string(pid));
      }
      else
        return (pid / 10000) % 1000;
    }

    unsigned int
    GetNucleusMassNumber(int pid)
    {
      if (!IsNucleus(pid)) {
        if (pid == 1000010010 || pid == 1000000010)
          return 1;
        else
          throw std::runtime_error("The particle is not a nucleus." + to_string(pid));
      }
      else
        return (pid / 10) % 1000;
    }

    int
    GetNucleusId(unsigned int Z, unsigned int A)
    {
      return kNucleusOffset + Z * 10000 + A * 10;
    }


    namespace {
      const map<string, int> nameToCharge =
        {
         {"H", 1},{"He", 2},{"Li", 3},{"Be", 4},{"B", 5},{"C", 6},{"N", 7},
         {"O", 8},{"F", 9},{"Ne", 10},{"Na", 11},{"Mg", 12},{"Al", 13},
         {"Si", 14}, {"P", 15},{"S", 16},{"Cl", 17},{"Ar", 18},{"K", 19},
         {"Ca", 20},{"Sc", 21}, {"Ti", 22},{"V", 23},{"Cr", 24},{"Mn", 25},
         {"Fe", 26},{"Co", 27},{"Ni", 28},{"Cu", 29},{"Zn", 30}
        };
      const map<int, string> chargeToName =
        {
         {1, "H"}, {2, "He"}, {3, "Li"}, {4, "Be"}, {5, "B"}, {6, "C"},
         {7, "N"}, {8, "O"}, {9, "F"}, {10, "Ne"}, {11, "Na"}, {12, "Mg"},
         {13, "Al"}, {14, "Si"}, {15, "P"}, {16, "S"}, {17, "Cl"}, {18, "Ar"},
         {19, "K"}, {20, "Ca"}, {21, "Sc"}, {22, "Ti"}, {23, "V"}, {24, "Cr"},
         {25, "Mn"}, {26, "Fe"}, {27, "Co"}, {28, "Ni"}, {29, "Cu"}, {30, "Zn"}
        };
    }

    int
    GetNucleusId(const string& name)
    {
      auto elementNameStart = name.find_first_not_of("0123456789");

      if (elementNameStart < 1 || elementNameStart == std::string::npos)
        FATAL("cannot decipher " + name);
      const string massString = name.substr(0, elementNameStart);
      const int A = stoi(massString);
      const string elementString = name.substr(elementNameStart, name.size() - elementNameStart);
      const auto iter = nameToCharge.find(elementString);
      if (iter != nameToCharge.end())
        return GetNucleusId(iter->second, A);
      else
        FATAL("unknown element " + elementString);
    }

    string
    GetElementName(const int id)
    {
      const int Z = GetNucleusCharge(id);
      return chargeToName.at(Z);
    }

    string
    GetName(const int id)
    {
      return to_string(GetNucleusMassNumber(id)) + GetElementName(id);
    }
  }
}
