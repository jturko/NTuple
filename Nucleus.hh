#ifndef __NUCLEUS_HH
#define __NUCLEUS_HH

#include <iostream>
#include <math.h>
#include <iomanip>
#include <stdlib.h>
#include <fstream>
#include <string.h>
#include <string>

class Nucleus{
 public:
  Nucleus();
  Nucleus(const char*);
  Nucleus(int, int, double, const char*);
  Nucleus(int Z, int N, const char*);
  void SetZ(int);
  void SetN(int);
  void SetMassExcess(double);  
  void SetMass(double);  
  void SetMass();  
  void SetSymbol(const char*);  
  int GetZfromSymbol(const char*);  
  int GetZ();
  int GetN();
  int GetA();
  double GetMassExcess();
  double GetMass();
  double GetRadius();
  const char* GetSymbol();
 private:
  int fZ;
  int fN;
  double fMass;
  double fMassExcess;
  std::string fSymbol;
};
#endif
