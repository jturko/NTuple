#ifndef __PARTICLE_HH
#define __PARTICLE_HH

#include <iostream>
#include <vector>
#include "TObject.h"

#include "TVector3.h"
#include "TLorentzVector.h"

class Particle : public TObject {
public:
  Particle(){};
  void Clear(){
    fType = 0;
    fTime = 0;
    fVertexZ = 0;            ///////////
    fEnergy = 0;
    fRecEnergy = 0;
    fDeltaE = 0;
    fERest = 0;
    fMult = -1;
    fPosition.SetXYZ(0,0,0);
    fDetected.SetPxPyPzE(0.,0.,0.,0.);
    fReconstructed.SetPxPyPzE(0.,0.,0.,0.);
    fEjectile.SetPxPyPzE(0.,0.,0.,0.);
    fRecoil.SetPxPyPzE(0.,0.,0.,0.);
    fDetector = -1;
    fRings.clear();
    fStrips.clear();
    fStripPos.clear();
    /*
    fEjectileE =0;
    fEjectileTh =0;
    fRecoilE =0;
    fRecoilTh =0;
    */
  };
  //set
  void SetType(int type){
    fType = type;
    //0 unident, 1 e-, 2 p, 3 d, 4 t, 6 alpha, 7 3He
    switch(fType%10){
    case 2:
      fMass = 938783.; // proton in keV/c^2 
      break;
    case 3:
      fMass = 1876122.; // deuteron in keV/c^2 
      break;
    case 4:
      fMass = 2809432.; // triton in keV/c^2 
      break;
    case 6:
      fMass = 3727379.109; // alpha in keV/c^2 
      break;
    case 7:
      fMass = 2809413.34; // 3He nucleus in keV/c^2 
      break;
    default:
      fMass=938783.; // proton in keV/c^2 
      break;
    }
    
  };
  void SetDetector(int det){
    fDetector = det;
  }; 
  void SetRings(std::vector<int> rings){
    fRings = rings;
  }
  void SetStrips(std::vector<int> strips){
    fStrips = strips;
  }
  void SetStripPos(std::vector<double> stripPos){
    fStripPos = stripPos;
  }
  void SetMultiplicity(int mult){
    fMult = mult;
  }; 
  void SetEnergy(double energy){
    fEnergy = energy;
  }; 
  void SetDeltaE(double energy){
    fDeltaE = energy;
  }; 
  void SetERest(double energy){
    fERest = energy;
  }; 
  void SetRecEnergy(double energy){
    fRecEnergy = energy;
  }; 
  void SetTime(double time){
    fTime = time;
  };
  void SetVertexZ(double vertexZ){
    fVertexZ = vertexZ;
  };
  void SetPosition(double x, double y, double z){
    fPosition.SetXYZ(x,y,z);
  };
  void SetPosition(TVector3 pos){
    fPosition = pos;
  };
  void SetDetected(){
    fDetected.SetVect(fPosition);
    fDetected.SetE(fMass+fEnergy);
    if(fPosition.Mag()>0)
      fDetected.SetRho( sqrt( (fMass+fEnergy)*(fMass+fEnergy) - fMass*fMass ) );
  }
  void SetReconstructed(){
    fReconstructed.SetVect(fPosition);
    fReconstructed.SetE(fMass+fRecEnergy);
    if(fPosition.Mag()>0)
      fReconstructed.SetRho( sqrt( (fMass+fRecEnergy)*(fMass+fRecEnergy) - fMass*fMass ) );
  }
  void SetEjectile(TLorentzVector lor){
    fEjectile = lor;
  }
  void SetRecoil(TLorentzVector lor){
    fRecoil = lor;
  }
  /*
  void SetEjectile(double energy, double angle){
    fEjectileE = energy;
    fEjectileTh = angle;
  };
  void SetRecoil(double energy, double angle){
    fRecoilE = energy;
    fRecoilTh = angle;
  };
  */
  //get
  int GetType(){
    return fType;
  };
  int GetDetector(){
    return fDetector;
  }; 
  std::vector<int> GetRings(){
    return fRings;
  }
  std::vector<int> GetStrips(){
    return fStrips;
  }
  std::vector<double> GetStripPos(){
    return fStripPos;
  }
  int GetMultiplicity(){
    return fMult;
  }; 
  double GetEnergy(){
    return fEnergy;
  };  
  double GetDeltaE(){
    return fDeltaE;
  };  
  double GetERest(){
    return fERest;
  };  
  double GetRecEnergy(){
    return fRecEnergy;
  };  
  double GetTime(){
    return fTime;
  };
  double GetVertexZ(){
    return fVertexZ;
  };
  TVector3 GetPosition(){
    return fPosition;
  };
  TLorentzVector GetDetected(){
    return fDetected;
  }
  TLorentzVector GetReconstructed(){
    return fReconstructed;
  }
  /*
  double GetEjectileE(){
    return fEjectileE;
  }
  double GetEjectileTh(){
    return fEjectileTh;
  }
  double GetRecoilE(){
    return fRecoilE;
  }
  double GetRecoilTh(){
    return fRecoilTh;
  }
  */
  TLorentzVector GetEjectile(){
    return fEjectile;
  }
  TLorentzVector GetRecoil(){
    return fRecoil;
  }

protected:
  int fType;
  int fDetector;//0 -3  forward, 10-13 backward, 20-23 cd
  std::vector<int> fRings;
  std::vector<int> fStrips;
  std::vector<double> fStripPos;
  int fMult;//0 -3  forward, 10-13 backward, 20-23 cd
  double fEnergy;
  double fRecEnergy;
  double fDeltaE;
  double fERest;
  TVector3 fPosition;
  double fTime;
  double fVertexZ;
  double fMass;
  TLorentzVector fDetected;
  TLorentzVector fReconstructed;
  TLorentzVector fEjectile;
  TLorentzVector fRecoil;
  /*
  double fEjectileE;
  double fEjectileTh;
  double fRecoilE;
  double fRecoilTh;
  */
  ClassDef(Particle, 3);
};

#endif
