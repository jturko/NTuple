#ifndef __CONVERTER_HH
#define __CONVERTER_HH

#include <vector>
#include <utility>
#include <map>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
//#include "TH3I.h"
#include "THnSparse.h"
#include "TVector3.h"

#include "Settings.hh"
#include "Griffin.hh"

#include "Particle.hh"
#include "ParticleMC.hh"
#include "Kinematics.hh"

class Converter {
public:
    Converter(std::vector<std::string>&, const std::string&, Settings*);
    ~Converter();

    bool Run();

private:
    bool AboveThreshold(double, int);
    bool InsideTimeWindow();
    bool DescantNeutronDiscrimination();

    // GRIFFIN
    void CheckGriffinCrystalAddback();
    void SupressGriffin();
    void SupressGriffinByNeighbouringAncillaryBgos();
    void SupressGriffinBySceptar();
    void AddbackGriffin();
    void AddbackGriffinNeighbour();
    void AddbackGriffinNeighbourVector();
    // LaBr
    void CheckLaBrDetectorAddback();
    void SupressLaBr();
    void SupressLaBrByNeighbouringGriffinShields();
    void AddbackLaBr();
    // EightPi
    void CheckEightPiDetectorAddback();
    void SupressEightPi();
    void AddbackEightPi();
    // Ancillary BGO
    void CheckAncillaryBgoCrystalAddback();
    void AddbackAncillaryBgo();
    // SCEPTAR
    void CheckSceptarDetectorAddback();
    void AddbackSceptar();
    // DESCANT
    void CheckDescantDetectorAddback();
    void AddbackDescant();
    // Paces
    void CheckPacesDetectorAddback();
    void AddbackPaces();
    
    // TI-STAR
    int CalculateTistarStripNumber(int layerNb, TVector3 particlePos, TVector3 stripPos, TVector3 stripDim);
    int CalculateTistarRingNumber (int layerNb, TVector3 particlePos, TVector3 stripPos, TVector3 stripDim);
    void FillTistarVectors();
    void ClearTistarVectors();
    void FillTistarParticleMCs();
    void ClearTistarParticleMCs();
    void CreateTistarHistograms(Kinematics*);

    void PrintStatistics();

    TH1F* Get1DHistogram(std::string, std::string, int, double, double);
    TH1F* Get1DHistogram(std::string, std::string);
    TH2F* Get2DHistogram(std::string, std::string, int, double, double, int, double, double);
    TH2F* Get2DHistogram(std::string, std::string);
    //TH3I* Get3DHistogram(std::string, std::string);
    THnSparseF* GetNDHistogram(std::string, std::string);

    void FillHistDetector1DGamma(TH1F* hist1D, std::vector<Detector>* detector, std::string hist_name, std::string hist_dir);
    void FillHistDetector2DGammaGamma(TH2F* hist2D, std::vector<Detector>* detector, std::string hist_name, std::string hist_dir);

    void FillHistDetector1DGammaNR(TH1F* hist1D, std::vector<Detector>* detector, std::string hist_name, std::string hist_dir);
    void FillHistDetector2DGammaGammaNR(TH2F* hist2D, std::vector<Detector>* detector, std::string hist_name, std::string hist_dir);

    void FillHist2DGriffinSceptarHitPattern(TH2F* hist2D, std::vector<Detector>* detector1, std::vector<Detector>* detector2, std::string hist_name, std::string hist_dir);
    void FillHist2DGriffinHitPattern(TH2F* hist2D, std::vector<Detector>* detector, std::string hist_name, std::string hist_dir);

    TVector3 GriffinCrystalCenterPosition(int cry, int det);
    bool AreGriffinCrystalCenterPositionsWithinVectorLength(int cry1, int det1, int cry2, int det2);

    double transX(double x, double y, double z, double theta, double phi);
    double transY(double x, double y, double z, double theta, double phi);
    double transZ(double x, double y, double z, double theta, double phi);

    Settings* fSettings;
    TChain fChain;
    TFile* fOutput;
    TTree fTree;
    TRandom3 fRandom;

    Int_t LaBrGriffinNeighbours_det[8][3];
    Int_t LaBrGriffinNeighbours_cry[8][3];

    Int_t GriffinAncillaryBgoNeighbours_det[16][2];
    Int_t GriffinAncillaryBgoNeighbours_cry[16][2];

    Int_t GriffinSceptarSuppressors_det[16][4];

    Int_t GriffinNeighbours_counted[16];
    Int_t GriffinNeighbours_det[16][4];

    double GriffinDetCoords[16][5];

    double GriffinCryMap[64][64];
    double GriffinCryMapCombos[52][2];

    double GriffinDetMap[16][16];
    double GriffinDetMapCombos[7][2];

    TVector3 GriffinCrystalCenterVectors[64];

    //branches of input tree/chain
    Int_t fEventNumber;
    Int_t fTrackID;
    Int_t fParentID;
    Int_t fStepNumber;
    Int_t fParticleType;
    Int_t fProcessType;
    Int_t fSystemID;
    Int_t fCryNumber;
    Int_t fDetNumber;
    Double_t fDepEnergy;
    Double_t fPosx;
    Double_t fPosy;
    Double_t fPosz;
    Double_t fTime;
    Int_t    fTargetZ;
    Double_t fTargetA;

    double LightOutput(double E, std::vector<double> & coeff) {
        return ( coeff[0]*E-coeff[1]*(1.-TMath::Exp(-1.0*coeff[2]*TMath::Power(E,coeff[3]))) );
    }

    bool   fSceptarHit;

    //branches of output tree
    // GRIFFIN
    std::vector<Detector>* fGriffinCrystal;
    std::vector<Detector>* fGriffinDetector;
    std::vector<Detector>* fGriffinNeighbour;
    std::vector<Detector>* fGriffinNeighbourVector;
    std::vector<Detector>* fGriffinArray;
    std::vector<Detector>* fGriffinBgo;
    std::vector<Detector>* fGriffinBgoBack;

    // LaBr
    std::vector<Detector>* fLaBrArray;
    std::vector<Detector>* fLaBrDetector;

    // EightPi
    std::vector<Detector>* fEightPiArray;
    std::vector<Detector>* fEightPiDetector;
    std::vector<Detector>* fEightPiBgoDetector;

    // Ancillary BGO
    std::vector<Detector>* fAncillaryBgoArray;
    std::vector<Detector>* fAncillaryBgoCrystal;
    std::vector<Detector>* fAncillaryBgoDetector;

    // Sceptar
    std::vector<Detector>* fSceptarArray;
    std::vector<Detector>* fSceptarDetector;

    // Descant
    std::vector<Detector>* fDescantArray;
    std::vector<Detector>* fDescantBlueDetector;
    std::vector<Detector>* fDescantGreenDetector;
    std::vector<Detector>* fDescantRedDetector;
    std::vector<Detector>* fDescantWhiteDetector;
    std::vector<Detector>* fDescantYellowDetector;

    // Testcan
    std::vector<Detector>* fTestcanDetector;

    // Paces
    std::vector<Detector>* fPacesArray;
    std::vector<Detector>* fPacesDetector;

    // TI-STAR
    std::vector<Detector>* fTISTARArray;
    std::vector<Detector>* fTISTARLayer1;
    std::vector<Detector>* fTISTARLayer2;
    std::vector<Detector>* fTISTARLayer3;
    
    //histograms
    std::map<std::string,TList*> fHistograms;
    
    // from the TRex-derived generator tree for TI-STAR
    TChain fTISTARGenChain; 
    
    Double_t fTISTARGenReactionBeamEnergy;
    Double_t fTISTARGenReactionBeamEnergyCM;
    Double_t fTISTARGenReactionX;
    Double_t fTISTARGenReactionY;
    Double_t fTISTARGenReactionZ;
    Double_t fTISTARGenRecoilTheta;
    Double_t fTISTARGenRecoilPhi;
    Double_t fTISTARGenRecoilEnergy;
    Int_t    fTISTARGenReaction;
    
    std::vector<Particle> * fTISTARParticleVector;
    std::vector<ParticleMC>* fTISTARFirstDeltaE[4];
    std::vector<ParticleMC>* fTISTARSecondDeltaE[2];
    std::vector<ParticleMC>* fTISTARPad[2];

    // same trackers as in the TRex/TI-STAR sensitive detector
    // as we loop over the ti-star hits, we will fill these, 
    // then once the last hit/entry of the ntuple has been processed
    // for that given event, we will fill the ParticleMC class and do
    // the TI-STAR analysis
    std::vector<int>        fTISTARFirstLayerStripNb[4];
    std::vector<double>     fTISTARFirstLayerStripEnergy[4];
    std::vector<int>        fTISTARFirstLayerStripA[4];
    std::vector<int>        fTISTARFirstLayerStripZ[4];
    std::vector<int>        fTISTARFirstLayerStripTrackID[4];
    std::vector<double>     fTISTARFirstLayerStripTime[4];
    std::vector<TVector3>   fTISTARFirstLayerStripPos[4]; // leila
    std::vector<int>        fTISTARFirstLayerStripStopped[4];
    std::vector<int>        fTISTARFirstLayerRingNb[4];
    std::vector<double>     fTISTARFirstLayerRingEnergy[4];
    std::vector<int>        fTISTARFirstLayerRingA[4];
    std::vector<int>        fTISTARFirstLayerRingZ[4];
    std::vector<int>        fTISTARFirstLayerRingTrackID[4];
    std::vector<double>     fTISTARFirstLayerRingTime[4];
    std::vector<int>        fTISTARFirstLayerRingStopped[4];
    
    std::vector<int>        fTISTARSecondLayerStripNb[2];
    std::vector<double>     fTISTARSecondLayerStripEnergy[2];
    std::vector<int>        fTISTARSecondLayerStripA[2];
    std::vector<int>        fTISTARSecondLayerStripZ[2];
    std::vector<int>        fTISTARSecondLayerStripTrackID[2];
    std::vector<double>     fTISTARSecondLayerStripTime[2];
    std::vector<TVector3>   fTISTARSecondLayerStripPos[2]; // leila
    std::vector<int>        fTISTARSecondLayerStripStopped[2];
    std::vector<int>        fTISTARSecondLayerRingNb[2];
    std::vector<double>     fTISTARSecondLayerRingEnergy[2];
    std::vector<int>        fTISTARSecondLayerRingA[2];
    std::vector<int>        fTISTARSecondLayerRingZ[2];
    std::vector<int>        fTISTARSecondLayerRingTrackID[2];
    std::vector<double>     fTISTARSecondLayerRingTime[2];
    std::vector<int>        fTISTARSecondLayerRingStopped[2];

};

#endif
