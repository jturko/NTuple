#ifndef __SETTINGS_HH
#define __SETTINGS_HH

#include <string>
#include <map>
#include <vector>

#include "TF1.h"

class Settings {
public:
    Settings(std::string, int);
    ~Settings(){};

    std::string NtupleName() {
        return fNtupleName;
    }

    int VerbosityLevel() {
        return fVerbosityLevel;
    }

    int BufferSize() {
        return fBufferSize;
    }

    int SortNumberOfEvents() {
        return fSortNumberOfEvents;
    }

    bool WriteTree() {
        return fWriteTree;
    }

    bool Write2DHist() {
        return fWrite2DHist;
    }

    //bool Write3DHist() {
    //    return fWrite3DHist;
    //}

    bool WriteNDHist() {
        return fWriteNDHist;
    }
    
    bool Write2DSGGHist() {
        return fWrite2DSGGHist;
    }

    bool WriteGriffinAddbackVector() {
        return fWriteGriffinAddbackVector;
    }

    double GriffinAddbackVectorLengthmm() {
        return fGriffinAddbackVectorLengthmm;
    }

    double GriffinAddbackVectorDepthmm() {
        return fGriffinAddbackVectorDepthmm;
    }

    double GriffinAddbackVectorCrystalFaceDistancemm() {
        return fGriffinAddbackVectorCrystalFaceDistancemm;
    }

    double Resolution(int systemID, int detectorID, int crystalID, double en) {
        if(fResolution.find(systemID) != fResolution.end()) {
            return fResolution[systemID].at(detectorID).at(crystalID).Eval(en);
        }
        return 0.;
    }
    double Threshold(int systemID, int detectorID, int crystalID) {
        if(fThreshold.find(systemID) != fThreshold.end()) {
            return fThreshold[systemID].at(detectorID).at(crystalID);
        }
        return 0.001;
    }
    double ThresholdWidth(int systemID, int detectorID, int crystalID) {
        if(fThresholdWidth.find(systemID) != fThresholdWidth.end()) {
            return fThresholdWidth[systemID].at(detectorID).at(crystalID);
        }
        return 0.;
    }
    double TimeWindow(int systemID, int detectorID, int crystalID) {
        if(fTimeWindow.find(systemID) != fTimeWindow.end()) {
            return fTimeWindow[systemID].at(detectorID).at(crystalID);
        }
        return 0.;
    }

    int NofBins(std::string directoryName) {
        if(fNofBins.find(directoryName) != fNofBins.end()) {
            return fNofBins[directoryName];
        }
        return 10000;
    }
    double RangeLow(std::string directoryName) {
        if(fRangeLow.find(directoryName) != fRangeLow.end()) {
            return fRangeLow[directoryName];
        }
        return 0.5;
    }
    double RangeHigh(std::string directoryName) {
        if(fRangeHigh.find(directoryName) != fRangeHigh.end()) {
            return fRangeHigh[directoryName];
        }
        return 10000.5;
    }

    enum Ion { kDeuteron, kCarbon, kProton, kAlpha, kElectron, kNeutron, kOther, kMax };
    double Quenching(Ion scintIon){
        if(scintIon < 0 || scintIon >= kMax) { std::cout << "error : scintIon outside of the enum range!" << std::endl; return -1.; }
        else return fQuenching[scintIon];
    }
    
    //double ProtonCoeff(int i) { return fProtonCoeff[i]; }
    //double DeuteronCoeff(int i) { return fDeuteronCoeff[i]; }
    //double CarbonCoeff(int i) { return fCarbonCoeff[i]; }
    //double BeCoeff(int i) { return fBeCoeff[i]; }
    //double BCoeff(int i) { return fBCoeff[i]; }
    //double AlphaCoeff(int i) { return fAlphaCoeff[i]; }

    std::vector<double> & ProtonCoeff() { return fProtonCoeff; }
    std::vector<double> & DeuteronCoeff() { return fDeuteronCoeff; }
    std::vector<double> & CarbonCoeff() { return fCarbonCoeff; }
    std::vector<double> & BeCoeff() { return fBeCoeff; }
    std::vector<double> & BCoeff() { return fBCoeff; }
    std::vector<double> & AlphaCoeff() { return fAlphaCoeff; }

private:
    std::string fNtupleName;

    int fVerbosityLevel;
    int fBufferSize;
    int fSortNumberOfEvents;

    bool fWriteTree;
    bool fWrite2DHist;
    //bool fWrite3DHist;
    bool fWriteNDHist;
    bool fWrite2DSGGHist;
    bool fWriteGriffinAddbackVector;

    double fGriffinAddbackVectorLengthmm;
    double fGriffinAddbackVectorDepthmm;
    double fGriffinAddbackVectorCrystalFaceDistancemm;

    std::map<int,std::vector<std::vector<TF1> > > fResolution;
    std::map<int,std::vector<std::vector<double> > > fThreshold;
    std::map<int,std::vector<std::vector<double> > > fThresholdWidth;
    std::map<int,std::vector<std::vector<double> > > fTimeWindow;

    std::map<std::string,int> fNofBins;
    std::map<std::string,double> fRangeLow;
    std::map<std::string,double> fRangeHigh;
    
    std::vector<double> fQuenching;

    std::vector<double> fProtonCoeff;
    std::vector<double> fDeuteronCoeff;
    std::vector<double> fCarbonCoeff;
    std::vector<double> fBeCoeff;
    std::vector<double> fBCoeff;
    std::vector<double> fAlphaCoeff;
    
};

#endif
