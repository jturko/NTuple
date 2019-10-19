/************************************************************************
 * Kinematic class for elastic scattering, Coulex and transfer reactions
 ***********************************************************************/

#ifndef __KINEMATICS_HH
#define __KINEMATICS_HH

#include <iostream>
#include <fstream>
#include <string>

#include "TMath.h"
#include "TSpline.h"
#include "TGraph.h"
#include "TLorentzVector.h"

#include "Nucleus.hh"

#ifndef PI
#define PI                       (TMath::Pi())
#endif

class Kinematics {
public:
  Kinematics();
  Kinematics(Nucleus* projectile, Nucleus* target, double ebeam);
  Kinematics(Nucleus* projectile, Nucleus* target, Nucleus* recoil, Nucleus* ejectile, double ebeam, double ex3);

  void Initial();
  void FinalCm();
  void Final(double angle, int part, bool upper = false);
  void SetEBeam (double ebeam) { fEBeam = ebeam; Initial(); };
  //void SetAngles(double angle, int part);
  void SetAngles(double angle, int part, bool upper=false);
  //get
  TSpline3* Evslab(double thmin, double thmax, double size, int part = 2, bool upper = false);
  TSpline3* Evscm(double thmin, double thmax, double size, int part = 2);
  TSpline3* labvscm(double thmin, double thmax, double size, int part = 2);
  TSpline3* cmvslab(double thmin, double thmax, double size, int part = 2);
  double GetQValue(){
    return fQValue;
  }
  double GetCmEnergy(double ebeam);
  double GetCmEnergy();
  double GetBeamEnergy(double LabAngle, double LabEnergy);
  double GetMaxAngle(double vcm);
  double GetMaxAngle(int part);
  double NormalkinEnergy();
  bool CheckMaxAngle(double angle, int part);
  double GetExcEnergy(TLorentzVector recoil, bool verbose = false);
  double ELab(double angle_lab, int part);
  double GetElab(int i){
    return fE[i];
  }
  double GetM(int i){
    return fM[i];
  }
  double GetTlab(int i){
    //    cout << " GetTlab start, return fT[i] " << fT[i] << "\ti" << i << endl;
    return fT[i];
  }
  double GetEcm(int i){
    return fEcm[i];
  }
  double GetTcm(int i){
    return fTcm[i];
  }
  double GetThetalab(int i){
    //    cout << " GetThetalab start " << endl;
    if(fTheta[i]<1e-5){
      //  cout << " GetThetalab stop, return 000000000000 " << endl;
      return 0;
    }
    else{
      //    cout << " GetThetalab stop, return fTheta[i] " << fTheta[i] << "\ti" << i << endl;
      return fTheta[i];
    }
  }
  double GetThetacm(int i){
    if(fThetacm[i]<1e-5)
      return 0;
    else
      return fThetacm[i];
  }
  double GetBetacm(){
    return fBeta_cm;
  }
  double GetGammacm(){
    return fGamma_cm;
  }
  double GetBetacm(int i){
    return fBetacm[i];
  }
  double GetVcm(int i){
    return fVcm[i];
  }
  double GetV(int i){
    return fV[i];
  } 
  TSpline3* Ruthvscm(double thmin, double thmax, double size); 
  TSpline3* Ruthvslab(double thmin, double thmax, double size, int part); 
  double Angle_lab2cm(double vcm, double angle_lab);
  double Angle_lab2cminverse(double vcm, double angle_lab, bool upper);
  double Angle_cm2lab(double vcm, double angle_cm);
  double Sigma_cm2lab(double angle_cm, double sigma_cm);
  double Sigma_lab2cm(double angle_cm, double sigma_lab);
  void Transform2cm(double &angle, double &sigma);
  void Transform2cm(double &angle, double &errangle, double &sigma, double &errsigma);
  void AngleErr_lab2cm(double angle, double &err);
  void SigmaErr_lab2cm(double angle, double err, double &sigma, double &errsigma);
  //double Sigma_cm2labnew(double vcm, double angle_cm, double sigma_cm);
  double Rutherford(double angle_cm);

    
private:

  double fTCm_i;
  double fTCm_f;

  Nucleus* fParticle[4];
  double fM[4];
  double fQValue;
  double fEBeam;
  double fT[4];
  double fE[4];
  double fP[4];
  double fV[4];
  double fTheta[4];

  double fTcm[4];
  double fEcm[4];
  double fPcm[4];
  double fVcm[4];
  double fBetacm[4];
  double fThetacm[4];

  double fBeta_cm;
  double fGamma_cm;

  double Pcm_em(double, double);
  double P_tm(double, double);
  double E_tm(double, double);
  double T_em(double, double);
  double betacm_tm(double, double);
  double V_pe(double, double);
  double E_final(int);
  double T_final(int);
};
#endif
