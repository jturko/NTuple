#ifndef __HIT_HH
#define __HIT_HH

#include <iostream>
#include <fstream>
#include <string>

#include "TROOT.h"
#include "TMath.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TFile.h"

#include "ParticleMC.hh"
#include "TistarSettings.hh"

#include "Settings.hh"

#ifndef PI
#define PI                       (TMath::Pi())
#endif

enum Direction { kForward, kBackward, kUndefined };

class HitSim {
public:
  HitSim(Settings * settings);
 
  void Clear();
  void SetFirstDeltaE(ParticleMC& firstbarrel, Direction direction);
  void SetSecondDeltaE(ParticleMC& secondbarrel, Direction direction);
  void SetPad(ParticleMC& pad);

  TVector3 FirstPosition(bool doubleSidedFirstLayer, bool smear);
  TVector3 SecondPosition(bool smear);

  double GetFirstDeltaEEnergy(bool verbose = false);
  double GetSecondDeltaEEnergy(bool verbose = false);
  double GetPadEnergy();

private:
  Settings*     fSett;
  TistarSettings* fTRexSett;
  TRandom* fRand;
  ParticleMC* fFirstDeltaE;
  ParticleMC* fSecondDeltaE;
  ParticleMC* fPad;
  Direction fFirstDirection;
  Direction fSecondDirection;

  // variables to hold results
  TVector3 fFirstPosition;
  TVector3 fSecondPosition;
  double fFirstEnergy;
  double fSecondEnergy;
};
#endif
