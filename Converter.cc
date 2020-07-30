#include "Converter.hh"

#include <iostream>
#include <iomanip>

#include "TMath.h"

#include "Utilities.hh"
#include "LightYield.hh"

#include "TistarSettings.hh"
#include "ParticleMC.hh"

#include "HitSim.hh"
#include "Compound.hh"
#include "Kinematics.hh"
#include "Reconstruction.hh"

#include "TSpline.h"

Converter::Converter(std::vector<std::string>& inputFileNames, const std::string& outputFileName, Settings* settings)
    : fSettings(settings) {
    TistarSettings * trex_settings = NULL;
    //create TChain to read in all input files
    for(auto fileName = inputFileNames.begin(); fileName != inputFileNames.end(); ++fileName) {
        if(!FileExists(*fileName)) {
            std::cerr<<"Failed to find file '"<<*fileName<<"', skipping it!"<<std::endl;
            continue;
        }
        //add sub-directory and tree name to file name
        //fileName->append(fSettings->NtupleName());
        //fChain.Add(fileName->c_str());
        std::string ntuple_name = *fileName + fSettings->NtupleName();
        fChain.Add(ntuple_name.c_str());
        
        std::string generator_ntuple_name = *fileName + fSettings->TISTARGenNtupleName();
        fTISTARGenChain.Add(generator_ntuple_name.c_str());

        fRandom.SetSeed(1);
        if(!trex_settings) {
            trex_settings = static_cast<TistarSettings*>(fChain.GetFile()->Get("settings"));
            trex_settings->Print();
            fSettings->SetTistarSettings(trex_settings);
        }
    }
 
    fTISTARGenGammaEnergy = 0;
    fTISTARGenGammaTheta = 0;
    fTISTARGenGammaPhi = 0;

    for(int i=0; i<4; i++) fTISTARFirstDeltaE[i] = new std::vector<ParticleMC>;
    for(int i=0; i<2; i++) fTISTARSecondDeltaE[i] = new std::vector<ParticleMC>;
    for(int i=0; i<2; i++) fTISTARPad[i] = new std::vector<ParticleMC>;

    // isn't f->GetListOfKeys()->Contains("graph") what you are looking for?

    //----------------------------------------------------------------------------------------------------
    // The following assumes detector 1 is placed is position 1, and detector 2 is placed in position 2, etc.
    // But that need not be the case!

    // Set LaBr3 Griffin Neighbours

    // LaBr3 detector 0 has three neighbouring Griffin detectors
    // Griffin detectors 0, 3, and 4. The "crystal" numbers for the side and extension suppressors for
    // those three detectors are 2, 1, and 4, respectively.
    LaBrGriffinNeighbours_det[0][0] = 0;
    LaBrGriffinNeighbours_det[0][1] = 3;
    LaBrGriffinNeighbours_det[0][2] = 4;
    LaBrGriffinNeighbours_cry[0][0] = 2;
    LaBrGriffinNeighbours_cry[0][1] = 0;
    LaBrGriffinNeighbours_cry[0][2] = 3;
    // next LaBr3 detector...
    LaBrGriffinNeighbours_det[1][0] = 0;
    LaBrGriffinNeighbours_det[1][1] = 1;
    LaBrGriffinNeighbours_det[1][2] = 2;
    LaBrGriffinNeighbours_cry[1][0] = 0;
    LaBrGriffinNeighbours_cry[1][1] = 2;
    LaBrGriffinNeighbours_cry[1][2] = 3;
    LaBrGriffinNeighbours_det[2][0] = 1;
    LaBrGriffinNeighbours_det[2][1] = 2;
    LaBrGriffinNeighbours_det[2][2] = 8;
    LaBrGriffinNeighbours_cry[2][0] = 0;
    LaBrGriffinNeighbours_cry[2][1] = 2;
    LaBrGriffinNeighbours_cry[2][2] = 3;
    LaBrGriffinNeighbours_det[3][0] = 2;
    LaBrGriffinNeighbours_det[3][1] = 3;
    LaBrGriffinNeighbours_det[3][2] = 10;
    LaBrGriffinNeighbours_cry[3][0] = 0;
    LaBrGriffinNeighbours_cry[3][1] = 2;
    LaBrGriffinNeighbours_cry[3][2] = 3;
    LaBrGriffinNeighbours_det[4][0] = 4;
    LaBrGriffinNeighbours_det[4][1] = 12;
    LaBrGriffinNeighbours_det[4][2] = 15;
    LaBrGriffinNeighbours_cry[4][0] = 1;
    LaBrGriffinNeighbours_cry[4][1] = 2;
    LaBrGriffinNeighbours_cry[4][2] = 0;
    LaBrGriffinNeighbours_det[5][0] = 6;
    LaBrGriffinNeighbours_det[5][1] = 12;
    LaBrGriffinNeighbours_det[5][2] = 13;
    LaBrGriffinNeighbours_cry[5][0] = 1;
    LaBrGriffinNeighbours_cry[5][1] = 0;
    LaBrGriffinNeighbours_cry[5][2] = 2;
    LaBrGriffinNeighbours_det[6][0] = 8;
    LaBrGriffinNeighbours_det[6][1] = 13;
    LaBrGriffinNeighbours_det[6][2] = 14;
    LaBrGriffinNeighbours_cry[6][0] = 1;
    LaBrGriffinNeighbours_cry[6][1] = 0;
    LaBrGriffinNeighbours_cry[6][2] = 2;
    LaBrGriffinNeighbours_det[7][0] = 10;
    LaBrGriffinNeighbours_det[7][1] = 14;
    LaBrGriffinNeighbours_det[7][2] = 15;
    LaBrGriffinNeighbours_cry[7][0] = 1;
    LaBrGriffinNeighbours_cry[7][1] = 0;
    LaBrGriffinNeighbours_cry[7][2] = 2;

    // GRIFFIN detector 1 (0 index) has two neighbouring ancillary BGO detectors
    // ancillary BGO detectors 0 and 1. The "crystal" or "segmentation" numbers for those BGOs are 1 and 2 respectively
    GriffinAncillaryBgoNeighbours_det[0][0] = 0;
    GriffinAncillaryBgoNeighbours_det[0][1] = 1;
    GriffinAncillaryBgoNeighbours_cry[0][0] = 1;
    GriffinAncillaryBgoNeighbours_cry[0][1] = 2;
    // Next GRIFFIN detector
    GriffinAncillaryBgoNeighbours_det[1][0] = 1;
    GriffinAncillaryBgoNeighbours_det[1][1] = 2;
    GriffinAncillaryBgoNeighbours_cry[1][0] = 1;
    GriffinAncillaryBgoNeighbours_cry[1][1] = 2;
    GriffinAncillaryBgoNeighbours_det[2][0] = 2;
    GriffinAncillaryBgoNeighbours_det[2][1] = 3;
    GriffinAncillaryBgoNeighbours_cry[2][0] = 1;
    GriffinAncillaryBgoNeighbours_cry[2][1] = 2;
    GriffinAncillaryBgoNeighbours_det[3][0] = 3;
    GriffinAncillaryBgoNeighbours_det[3][1] = 0;
    GriffinAncillaryBgoNeighbours_cry[3][0] = 1;
    GriffinAncillaryBgoNeighbours_cry[3][1] = 2;
    GriffinAncillaryBgoNeighbours_det[4][0] = 0;
    GriffinAncillaryBgoNeighbours_det[4][1] = 4;
    GriffinAncillaryBgoNeighbours_cry[4][0] = 1;
    GriffinAncillaryBgoNeighbours_cry[4][1] = 2;
    GriffinAncillaryBgoNeighbours_det[5][0] = 9999; // no ancillary bgo neighbours
    GriffinAncillaryBgoNeighbours_det[5][1] = 9999;
    GriffinAncillaryBgoNeighbours_cry[5][0] = 9999;
    GriffinAncillaryBgoNeighbours_cry[5][1] = 9999;
    GriffinAncillaryBgoNeighbours_det[6][0] = 1;
    GriffinAncillaryBgoNeighbours_det[6][1] = 5;
    GriffinAncillaryBgoNeighbours_cry[6][0] = 0;
    GriffinAncillaryBgoNeighbours_cry[6][1] = 0;
    GriffinAncillaryBgoNeighbours_det[7][0] = 9999;
    GriffinAncillaryBgoNeighbours_det[7][1] = 9999;
    GriffinAncillaryBgoNeighbours_cry[7][0] = 9999;
    GriffinAncillaryBgoNeighbours_cry[7][1] = 9999;
    GriffinAncillaryBgoNeighbours_det[8][0] = 2;
    GriffinAncillaryBgoNeighbours_det[8][1] = 6;
    GriffinAncillaryBgoNeighbours_cry[8][0] = 0;
    GriffinAncillaryBgoNeighbours_cry[8][1] = 0;
    GriffinAncillaryBgoNeighbours_det[9][0] = 9999;
    GriffinAncillaryBgoNeighbours_det[9][1] = 9999;
    GriffinAncillaryBgoNeighbours_cry[9][0] = 9999;
    GriffinAncillaryBgoNeighbours_cry[9][1] = 9999;
    GriffinAncillaryBgoNeighbours_det[10][0] = 3;
    GriffinAncillaryBgoNeighbours_det[10][1] = 7;
    GriffinAncillaryBgoNeighbours_cry[10][0] = 0;
    GriffinAncillaryBgoNeighbours_cry[10][1] = 0;
    GriffinAncillaryBgoNeighbours_det[11][0] = 9999;
    GriffinAncillaryBgoNeighbours_det[11][1] = 9999;
    GriffinAncillaryBgoNeighbours_cry[11][0] = 9999;
    GriffinAncillaryBgoNeighbours_cry[11][1] = 9999;
    GriffinAncillaryBgoNeighbours_det[12][0] = 4;
    GriffinAncillaryBgoNeighbours_det[12][1] = 5;
    GriffinAncillaryBgoNeighbours_cry[12][0] = 2;
    GriffinAncillaryBgoNeighbours_cry[12][1] = 1;
    GriffinAncillaryBgoNeighbours_det[13][0] = 5;
    GriffinAncillaryBgoNeighbours_det[13][1] = 6;
    GriffinAncillaryBgoNeighbours_cry[13][0] = 2;
    GriffinAncillaryBgoNeighbours_cry[13][1] = 1;
    GriffinAncillaryBgoNeighbours_det[14][0] = 6;
    GriffinAncillaryBgoNeighbours_det[14][1] = 7;
    GriffinAncillaryBgoNeighbours_cry[14][0] = 2;
    GriffinAncillaryBgoNeighbours_cry[14][1] = 1;
    GriffinAncillaryBgoNeighbours_det[15][0] = 7;
    GriffinAncillaryBgoNeighbours_det[15][1] = 4;
    GriffinAncillaryBgoNeighbours_cry[15][0] = 2;
    GriffinAncillaryBgoNeighbours_cry[15][1] = 1;

    // GRIFFIN detector 1 (0 index) has two sceptar suppressors
    // SCEPTAR detectors 0 and 1. The detector numbers for those paddles are 1 and 2 respectively
    GriffinSceptarSuppressors_det[0][0] = 0;
    GriffinSceptarSuppressors_det[0][1] = 1;
    GriffinSceptarSuppressors_det[0][2] = 9999;
    GriffinSceptarSuppressors_det[0][3] = 9999;
    // Next GRIFFIN detector
    GriffinSceptarSuppressors_det[1][0] = 2;
    GriffinSceptarSuppressors_det[1][1] = 9999;
    GriffinSceptarSuppressors_det[1][2] = 9999;
    GriffinSceptarSuppressors_det[1][3] = 9999;
    GriffinSceptarSuppressors_det[2][0] = 3;
    GriffinSceptarSuppressors_det[2][1] = 9999;
    GriffinSceptarSuppressors_det[2][2] = 9999;
    GriffinSceptarSuppressors_det[2][3] = 9999;
    GriffinSceptarSuppressors_det[3][0] = 0;
    GriffinSceptarSuppressors_det[3][1] = 4;
    GriffinSceptarSuppressors_det[3][2] = 9999;
    GriffinSceptarSuppressors_det[3][3] = 9999;
    GriffinSceptarSuppressors_det[4][0] = 5;
    GriffinSceptarSuppressors_det[4][1] = 10;
    GriffinSceptarSuppressors_det[4][2] = 9999;
    GriffinSceptarSuppressors_det[4][3] = 9999;
    GriffinSceptarSuppressors_det[5][0] = 5;
    GriffinSceptarSuppressors_det[5][1] = 10;
    GriffinSceptarSuppressors_det[5][2] = 6;
    GriffinSceptarSuppressors_det[5][3] = 11;
    GriffinSceptarSuppressors_det[6][0] = 6;
    GriffinSceptarSuppressors_det[6][1] = 11;
    GriffinSceptarSuppressors_det[6][2] = 7;
    GriffinSceptarSuppressors_det[6][3] = 12;
    GriffinSceptarSuppressors_det[7][0] = 7;
    GriffinSceptarSuppressors_det[7][1] = 12;
    GriffinSceptarSuppressors_det[7][2] = 9999;
    GriffinSceptarSuppressors_det[7][3] = 9999;
    GriffinSceptarSuppressors_det[8][0] = 7;
    GriffinSceptarSuppressors_det[8][1] = 12;
    GriffinSceptarSuppressors_det[8][2] = 8;
    GriffinSceptarSuppressors_det[8][3] = 13;
    GriffinSceptarSuppressors_det[9][0] = 8;
    GriffinSceptarSuppressors_det[9][1] = 13;
    GriffinSceptarSuppressors_det[9][2] = 9999;
    GriffinSceptarSuppressors_det[9][3] = 9999;
    GriffinSceptarSuppressors_det[10][0] = 9;
    GriffinSceptarSuppressors_det[10][1] = 14;
    GriffinSceptarSuppressors_det[10][2] = 9999;
    GriffinSceptarSuppressors_det[10][3] = 9999;
    GriffinSceptarSuppressors_det[11][0] = 9;
    GriffinSceptarSuppressors_det[11][1] = 14;
    GriffinSceptarSuppressors_det[11][2] = 5;
    GriffinSceptarSuppressors_det[11][3] = 10;
    GriffinSceptarSuppressors_det[12][0] = 16;
    GriffinSceptarSuppressors_det[12][1] = 15;
    GriffinSceptarSuppressors_det[12][2] = 9999;
    GriffinSceptarSuppressors_det[12][3] = 9999;
    GriffinSceptarSuppressors_det[13][0] = 17;
    GriffinSceptarSuppressors_det[13][1] = 9999;
    GriffinSceptarSuppressors_det[13][2] = 9999;
    GriffinSceptarSuppressors_det[13][3] = 9999;
    GriffinSceptarSuppressors_det[14][0] = 18;
    GriffinSceptarSuppressors_det[14][1] = 19;
    GriffinSceptarSuppressors_det[14][2] = 9999;
    GriffinSceptarSuppressors_det[14][3] = 9999;
    GriffinSceptarSuppressors_det[15][0] = 18;
    GriffinSceptarSuppressors_det[15][1] = 17;
    GriffinSceptarSuppressors_det[15][2] = 9999;
    GriffinSceptarSuppressors_det[15][3] = 9999;

    GriffinNeighbours_det[0][0] = 5;
    GriffinNeighbours_det[0][1] = 9999;
    GriffinNeighbours_det[0][2] = 9999;
    GriffinNeighbours_det[0][3] = 9999;
    // Next detector
    GriffinNeighbours_det[1][0] = 7;
    GriffinNeighbours_det[1][1] = 9999;
    GriffinNeighbours_det[1][2] = 9999;
    GriffinNeighbours_det[1][3] = 9999;
    GriffinNeighbours_det[2][0] = 9;
    GriffinNeighbours_det[2][1] = 9999;
    GriffinNeighbours_det[2][2] = 9999;
    GriffinNeighbours_det[2][3] = 9999;
    GriffinNeighbours_det[3][0] = 11;
    GriffinNeighbours_det[3][1] = 9999;
    GriffinNeighbours_det[3][2] = 9999;
    GriffinNeighbours_det[3][3] = 9999;
    GriffinNeighbours_det[4][0] = 5;
    GriffinNeighbours_det[4][1] = 11;
    GriffinNeighbours_det[4][2] = 9999;
    GriffinNeighbours_det[4][3] = 9999;
    GriffinNeighbours_det[5][0] = 0;
    GriffinNeighbours_det[5][1] = 12;
    GriffinNeighbours_det[5][2] = 4;
    GriffinNeighbours_det[5][3] = 6;
    GriffinNeighbours_det[6][0] = 5;
    GriffinNeighbours_det[6][1] = 7;
    GriffinNeighbours_det[6][2] = 9999;
    GriffinNeighbours_det[6][3] = 9999;
    GriffinNeighbours_det[7][0] = 1;
    GriffinNeighbours_det[7][1] = 13;
    GriffinNeighbours_det[7][2] = 6;
    GriffinNeighbours_det[7][3] = 8;
    GriffinNeighbours_det[8][0] = 7;
    GriffinNeighbours_det[8][1] = 9;
    GriffinNeighbours_det[8][2] = 9999;
    GriffinNeighbours_det[8][3] = 9999;
    GriffinNeighbours_det[9][0] = 2;
    GriffinNeighbours_det[9][1] = 14;
    GriffinNeighbours_det[9][2] = 8;
    GriffinNeighbours_det[9][3] = 10;
    GriffinNeighbours_det[10][0] = 9;
    GriffinNeighbours_det[10][1] = 11;
    GriffinNeighbours_det[10][2] = 9999;
    GriffinNeighbours_det[10][3] = 9999;
    GriffinNeighbours_det[11][0] = 3;
    GriffinNeighbours_det[11][1] = 15;
    GriffinNeighbours_det[11][2] = 4;
    GriffinNeighbours_det[11][3] = 10;
    GriffinNeighbours_det[12][0] = 5;
    GriffinNeighbours_det[12][1] = 9999;
    GriffinNeighbours_det[12][2] = 9999;
    GriffinNeighbours_det[12][3] = 9999;
    GriffinNeighbours_det[13][0] = 7;
    GriffinNeighbours_det[13][1] = 9999;
    GriffinNeighbours_det[13][2] = 9999;
    GriffinNeighbours_det[13][3] = 9999;
    GriffinNeighbours_det[14][0] = 9;
    GriffinNeighbours_det[14][1] = 9999;
    GriffinNeighbours_det[14][2] = 9999;
    GriffinNeighbours_det[14][3] = 9999;
    GriffinNeighbours_det[15][0] = 11;
    GriffinNeighbours_det[15][1] = 9999;
    GriffinNeighbours_det[15][2] = 9999;
    GriffinNeighbours_det[15][3] = 9999;

    /////////////////////////////////////////////////////////////////////
    // Coords for GRIFFIN
    // Note that the GRIFFIN lampshade angles are rotated by 45 degrees with respect to those of TIGRESS.
    // Modified coords for TIGRESS are below!
    /////////////////////////////////////////////////////////////////////
    double thisGriffinDetCoords[16][5];
    // theta
    thisGriffinDetCoords[0][0] 	= 45.0;
    thisGriffinDetCoords[1][0] 	= 45.0;
    thisGriffinDetCoords[2][0] 	= 45.0;
    thisGriffinDetCoords[3][0] 	= 45.0;
    thisGriffinDetCoords[4][0] 	= 90.0;
    thisGriffinDetCoords[5][0] 	= 90.0;
    thisGriffinDetCoords[6][0] 	= 90.0;
    thisGriffinDetCoords[7][0] 	= 90.0;
    thisGriffinDetCoords[8][0] 	= 90.0;
    thisGriffinDetCoords[9][0] 	= 90.0;
    thisGriffinDetCoords[10][0] 	= 90.0;
    thisGriffinDetCoords[11][0] 	= 90.0;
    thisGriffinDetCoords[12][0] 	= 135.0;
    thisGriffinDetCoords[13][0] 	= 135.0;
    thisGriffinDetCoords[14][0] 	= 135.0;
    thisGriffinDetCoords[15][0] 	= 135.0;
    // phi
    thisGriffinDetCoords[0][1] 	= 67.5;
    thisGriffinDetCoords[1][1] 	= 157.5;
    thisGriffinDetCoords[2][1] 	= 247.5;
    thisGriffinDetCoords[3][1] 	= 337.5;
    thisGriffinDetCoords[4][1] 	= 22.5;
    thisGriffinDetCoords[5][1] 	= 67.5;
    thisGriffinDetCoords[6][1] 	= 112.5;
    thisGriffinDetCoords[7][1] 	= 157.5;
    thisGriffinDetCoords[8][1] 	= 202.5;
    thisGriffinDetCoords[9][1] 	= 247.5;
    thisGriffinDetCoords[10][1] 	= 292.5;
    thisGriffinDetCoords[11][1] 	= 337.5;
    thisGriffinDetCoords[12][1] 	= 67.5;
    thisGriffinDetCoords[13][1] 	= 157.5;
    thisGriffinDetCoords[14][1] 	= 247.5;
    thisGriffinDetCoords[15][1] 	= 337.5;
    // yaw (alpha)
    thisGriffinDetCoords[0][2] 	= 0.0;
    thisGriffinDetCoords[1][2] 	= 0.0;
    thisGriffinDetCoords[2][2] 	= 0.0;
    thisGriffinDetCoords[3][2] 	= 0.0;
    thisGriffinDetCoords[4][2] 	= 0.0;
    thisGriffinDetCoords[5][2] 	= 0.0;
    thisGriffinDetCoords[6][2] 	= 0.0;
    thisGriffinDetCoords[7][2] 	= 0.0;
    thisGriffinDetCoords[8][2] 	= 0.0;
    thisGriffinDetCoords[9][2] 	= 0.0;
    thisGriffinDetCoords[10][2] 	= 0.0;
    thisGriffinDetCoords[11][2] 	= 0.0;
    thisGriffinDetCoords[12][2] 	= 0.0;
    thisGriffinDetCoords[13][2] 	= 0.0;
    thisGriffinDetCoords[14][2] 	= 0.0;
    thisGriffinDetCoords[15][2] 	= 0.0;
    // pitch (beta)
    thisGriffinDetCoords[0][3] 	= -45.0;
    thisGriffinDetCoords[1][3] 	= -45.0;
    thisGriffinDetCoords[2][3] 	= -45.0;
    thisGriffinDetCoords[3][3] 	= -45.0;
    thisGriffinDetCoords[4][3] 	= 0.0;
    thisGriffinDetCoords[5][3] 	= 0.0;
    thisGriffinDetCoords[6][3] 	= 0.0;
    thisGriffinDetCoords[7][3] 	= 0.0;
    thisGriffinDetCoords[8][3] 	= 0.0;
    thisGriffinDetCoords[9][3] 	= 0.0;
    thisGriffinDetCoords[10][3] 	= 0.0;
    thisGriffinDetCoords[11][3] 	= 0.0;
    thisGriffinDetCoords[12][3] 	= 45.0;
    thisGriffinDetCoords[13][3] 	= 45.0;
    thisGriffinDetCoords[14][3] 	= 45.0;
    thisGriffinDetCoords[15][3] 	= 45.0;
    // roll (gamma)
    thisGriffinDetCoords[0][4] 	= 67.5;
    thisGriffinDetCoords[1][4] 	= 157.5;
    thisGriffinDetCoords[2][4] 	= 247.5;
    thisGriffinDetCoords[3][4] 	= 337.5;
    thisGriffinDetCoords[4][4] 	= 22.5;
    thisGriffinDetCoords[5][4] 	= 67.5;
    thisGriffinDetCoords[6][4] 	= 112.5;
    thisGriffinDetCoords[7][4] 	= 157.5;
    thisGriffinDetCoords[8][4] 	= 202.5;
    thisGriffinDetCoords[9][4] 	= 247.5;
    thisGriffinDetCoords[10][4] 	= 292.5;
    thisGriffinDetCoords[11][4] 	= 337.5;
    thisGriffinDetCoords[12][4] 	= 67.5;
    thisGriffinDetCoords[13][4] 	= 157.5;
    thisGriffinDetCoords[14][4] 	= 247.5;
    thisGriffinDetCoords[15][4] 	= 337.5;
    memcpy(GriffinDetCoords, thisGriffinDetCoords, sizeof(GriffinDetCoords));


    // Detector Method
    double thisGriffinDetMap[16][16] = {
        {0,60,90,60,60,45,60,90,120,135,120,90,90,120,180,120},
        {60,0,60,90,120,90,60,45,60,90,120,135,120,90,120,180},
        {90,60,0,60,120,135,120,90,60,45,60,90,180,120,90,120},
        {60,90,60,0,60,90,120,135,120,90,60,45,120,180,120,90},
        {60,120,120,60,0,45,90,135,180,135,90,45,60,120,120,60},
        {45,90,135,90,45,0,45,90,135,180,135,90,45,90,135,90},
        {60,60,120,120,90,45,0,45,90,135,180,135,60,60,120,120},
        {90,45,90,135,135,90,45,0,45,90,135,180,90,45,90,135},
        {120,60,60,120,180,135,90,45,0,45,90,135,120,60,60,120},
        {135,90,45,90,135,180,135,90,45,0,45,90,135,90,45,90},
        {120,120,60,60,90,135,180,135,90,45,0,45,120,120,60,60},
        {90,135,90,45,45,90,135,180,135,90,45,0,90,135,90,45},
        {90,120,180,120,60,45,60,90,120,135,120,90,0,60,90,60},
        {120,90,120,180,120,90,60,45,60,90,120,135,60,0,60,90},
        {180,120,90,120,120,135,120,90,60,45,60,90,90,60,0,60},
        {120,180,120,90,60,90,120,135,120,90,60,45,60,90,60,0}
    };
    memcpy(GriffinDetMap, thisGriffinDetMap, sizeof(GriffinDetMap));

    double thisGriffinDetMapCombos[7][2] = {
        {0, 16},
        {45, 32},
        {60, 48},
        {90, 64},
        {120, 48},
        {135, 32},
        {180, 16}
    };
    memcpy(GriffinDetMapCombos, thisGriffinDetMapCombos, sizeof(GriffinDetMapCombos));

    double thisGriffinCryMap[64][64] = {
        {0.0000, 19.131, 27.184, 19.131, 49.631, 60.157, 46.607, 33.166, 72.817, 91.582, 88.418, 69.473, 49.631, 65.195, 76.694, 62.720, 60.157, 76.694, 86.721, 71.054, 44.341, 63.403, 66.891, 48.703, 53.690, 71.054, 65.195, 46.607, 78.429, 93.836, 82.965, 67.049, 103.31, 119.84, 108.95, 93.279, 116.60, 135.66, 131.30, 113.11, 108.95, 126.31, 133.39, 114.81, 86.164, 101.57, 112.95, 97.035, 88.418, 107.18, 110.53, 91.582, 114.81, 130.37, 117.28, 103.31, 160.87, 180.00, 160.87, 152.82, 119.84, 130.37, 146.83, 133.39},
        {19.131, 0.0000, 19.131, 27.184, 65.195, 71.054, 53.690, 46.607, 91.582, 110.53, 107.18, 88.418, 60.157, 71.054, 86.721, 76.694, 49.631, 62.720, 76.694, 65.195, 25.235, 44.341, 48.703, 31.860, 46.607, 60.157, 49.631, 33.166, 82.965, 93.836, 78.429, 67.049, 117.28, 130.37, 114.81, 103.31, 135.66, 154.77, 148.14, 131.30, 119.84, 133.39, 146.83, 130.37, 86.164, 97.035, 112.95, 101.57, 69.473, 88.418, 91.582, 72.817, 108.95, 119.84, 103.31, 93.279, 180.00, 160.87, 152.82, 160.87, 108.95, 114.81, 133.39, 126.31},
        {27.184, 19.131, 0.0000, 19.131, 76.694, 86.721, 71.054, 60.157, 88.418, 107.18, 110.53, 91.582, 46.607, 53.690, 71.054, 65.195, 33.166, 49.631, 60.157, 46.607, 31.860, 48.703, 44.341, 25.235, 65.195, 76.694, 62.720, 49.631, 101.57, 112.95, 97.035, 86.164, 130.37, 146.83, 133.39, 119.84, 131.30, 148.14, 154.77, 135.66, 103.31, 114.81, 130.37, 117.28, 67.049, 78.429, 93.836, 82.965, 72.817, 91.582, 88.418, 69.473, 126.31, 133.39, 114.81, 108.95, 160.87, 152.82, 160.87, 180.00, 93.279, 103.31, 119.84, 108.95},
        {19.131, 27.184, 19.131, 0.0000, 62.720, 76.694, 65.195, 49.631, 69.473, 88.418, 91.582, 72.817, 33.166, 46.607, 60.157, 49.631, 46.607, 65.195, 71.054, 53.690, 48.703, 66.891, 63.403, 44.341, 71.054, 86.721, 76.694, 60.157, 97.035, 112.95, 101.57, 86.164, 114.81, 133.39, 126.31, 108.95, 113.11, 131.30, 135.66, 116.60, 93.279, 108.95, 119.84, 103.31, 67.049, 82.965, 93.836, 78.429, 91.582, 110.53, 107.18, 88.418, 133.39, 146.83, 130.37, 119.84, 152.82, 160.87, 180.00, 160.87, 103.31, 117.28, 130.37, 114.81},
        {49.631, 65.195, 76.694, 62.720, 0.0000, 19.131, 27.184, 19.131, 49.631, 60.157, 46.607, 33.166, 72.817, 91.582, 88.418, 69.473, 108.95, 126.31, 133.39, 114.81, 86.164, 101.57, 112.95, 97.035, 60.157, 76.694, 86.721, 71.054, 44.341, 63.403, 66.891, 48.703, 53.690, 71.054, 65.195, 46.607, 78.429, 93.836, 82.965, 67.049, 103.31, 119.84, 108.95, 93.279, 116.60, 135.66, 131.30, 113.11, 119.84, 130.37, 146.83, 133.39, 88.418, 107.18, 110.53, 91.582, 114.81, 130.37, 117.28, 103.31, 160.87, 180.00, 160.87, 152.82},
        {60.157, 71.054, 86.721, 76.694, 19.131, 0.0000, 19.131, 27.184, 65.195, 71.054, 53.690, 46.607, 91.582, 110.53, 107.18, 88.418, 119.84, 133.39, 146.83, 130.37, 86.164, 97.035, 112.95, 101.57, 49.631, 62.720, 76.694, 65.195, 25.235, 44.341, 48.703, 31.860, 46.607, 60.157, 49.631, 33.166, 82.965, 93.836, 78.429, 67.049, 117.28, 130.37, 114.81, 103.31, 135.66, 154.77, 148.14, 131.30, 108.95, 114.81, 133.39, 126.31, 69.473, 88.418, 91.582, 72.817, 108.95, 119.84, 103.31, 93.279, 180.00, 160.87, 152.82, 160.87},
        {46.607, 53.690, 71.054, 65.195, 27.184, 19.131, 0.0000, 19.131, 76.694, 86.721, 71.054, 60.157, 88.418, 107.18, 110.53, 91.582, 103.31, 114.81, 130.37, 117.28, 67.049, 78.429, 93.836, 82.965, 33.166, 49.631, 60.157, 46.607, 31.860, 48.703, 44.341, 25.235, 65.195, 76.694, 62.720, 49.631, 101.57, 112.95, 97.035, 86.164, 130.37, 146.83, 133.39, 119.84, 131.30, 148.14, 154.77, 135.66, 93.279, 103.31, 119.84, 108.95, 72.817, 91.582, 88.418, 69.473, 126.31, 133.39, 114.81, 108.95, 160.87, 152.82, 160.87, 180},
        {33.166, 46.607, 60.157, 49.631, 19.131, 27.184, 19.131, 0.0000, 62.720, 76.694, 65.195, 49.631, 69.473, 88.418, 91.582, 72.817, 93.279, 108.95, 119.84, 103.31, 67.049, 82.965, 93.836, 78.429, 46.607, 65.195, 71.054, 53.690, 48.703, 66.891, 63.403, 44.341, 71.054, 86.721, 76.694, 60.157, 97.035, 112.95, 101.57, 86.164, 114.81, 133.39, 126.31, 108.95, 113.11, 131.30, 135.66, 116.60, 103.31, 117.28, 130.37, 114.81, 91.582, 110.53, 107.18, 88.418, 133.39, 146.83, 130.37, 119.84, 152.82, 160.87, 180.00, 160.87},
        {72.817, 91.582, 88.418, 69.473, 49.631, 65.195, 76.694, 62.720, 0.0000, 19.131, 27.184, 19.131, 49.631, 60.157, 46.607, 33.166, 103.31, 119.84, 108.95, 93.279, 116.60, 135.66, 131.30, 113.11, 108.95, 126.31, 133.39, 114.81, 86.164, 101.57, 112.95, 97.035, 60.157, 76.694, 86.721, 71.054, 44.341, 63.403, 66.891, 48.703, 53.690, 71.054, 65.195, 46.607, 78.429, 93.836, 82.965, 67.049, 160.87, 180.00, 160.87, 152.82, 119.84, 130.37, 146.83, 133.39, 88.418, 107.18, 110.53, 91.582, 114.81, 130.37, 117.28, 103.31},
        {91.582, 110.53, 107.18, 88.418, 60.157, 71.054, 86.721, 76.694, 19.131, 0.0000, 19.131, 27.184, 65.195, 71.054, 53.690, 46.607, 117.28, 130.37, 114.81, 103.31, 135.66, 154.77, 148.14, 131.30, 119.84, 133.39, 146.83, 130.37, 86.164, 97.035, 112.95, 101.57, 49.631, 62.720, 76.694, 65.195, 25.235, 44.341, 48.703, 31.860, 46.607, 60.157, 49.631, 33.166, 82.965, 93.836, 78.429, 67.049, 180.00, 160.87, 152.82, 160.87, 108.95, 114.81, 133.39, 126.31, 69.473, 88.418, 91.582, 72.817, 108.95, 119.84, 103.31, 93.279},
        {88.418, 107.18, 110.53, 91.582, 46.607, 53.690, 71.054, 65.195, 27.184, 19.131, 0.0000, 19.131, 76.694, 86.721, 71.054, 60.157, 130.37, 146.83, 133.39, 119.84, 131.30, 148.14, 154.77, 135.66, 103.31, 114.81, 130.37, 117.28, 67.049, 78.429, 93.836, 82.965, 33.166, 49.631, 60.157, 46.607, 31.860, 48.703, 44.341, 25.235, 65.195, 76.694, 62.720, 49.631, 101.57, 112.95, 97.035, 86.164, 160.87, 152.82, 160.87, 180.00, 93.279, 103.31, 119.84, 108.95, 72.817, 91.582, 88.418, 69.473, 126.31, 133.39, 114.81, 108.95},
        {69.473, 88.418, 91.582, 72.817, 33.166, 46.607, 60.157, 49.631, 19.131, 27.184, 19.131, 0.0000, 62.720, 76.694, 65.195, 49.631, 114.81, 133.39, 126.31, 108.95, 113.11, 131.30, 135.66, 116.60, 93.279, 108.95, 119.84, 103.31, 67.049, 82.965, 93.836, 78.429, 46.607, 65.195, 71.054, 53.690, 48.703, 66.891, 63.403, 44.341, 71.054, 86.721, 76.694, 60.157, 97.035, 112.95, 101.57, 86.164, 152.82, 160.87, 180.00, 160.87, 103.31, 117.28, 130.37, 114.81, 91.582, 110.53, 107.18, 88.418, 133.39, 146.83, 130.37, 119.84},
        {49.631, 60.157, 46.607, 33.166, 72.817, 91.582, 88.418, 69.473, 49.631, 65.195, 76.694, 62.720, 0.0000, 19.131, 27.184, 19.131, 53.690, 71.054, 65.195, 46.607, 78.429, 93.836, 82.965, 67.049, 103.31, 119.84, 108.95, 93.279, 116.60, 135.66, 131.30, 113.11, 108.95, 126.31, 133.39, 114.81, 86.164, 101.57, 112.95, 97.035, 60.157, 76.694, 86.721, 71.054, 44.341, 63.403, 66.891, 48.703, 114.81, 130.37, 117.28, 103.31, 160.87, 180.00, 160.87, 152.82, 119.84, 130.37, 146.83, 133.39, 88.418, 107.18, 110.53, 91.582},
        {65.195, 71.054, 53.690, 46.607, 91.582, 110.53, 107.18, 88.418, 60.157, 71.054, 86.721, 76.694, 19.131, 0.0000, 19.131, 27.184, 46.607, 60.157, 49.631, 33.166, 82.965, 93.836, 78.429, 67.049, 117.28, 130.37, 114.81, 103.31, 135.66, 154.77, 148.14, 131.30, 119.84, 133.39, 146.83, 130.37, 86.164, 97.035, 112.95, 101.57, 49.631, 62.720, 76.694, 65.195, 25.235, 44.341, 48.703, 31.860, 108.95, 119.84, 103.31, 93.279, 180.00, 160.87, 152.82, 160.87, 108.95, 114.81, 133.39, 126.31, 69.473, 88.418, 91.582, 72.817},
        {76.694, 86.721, 71.054, 60.157, 88.418, 107.18, 110.53, 91.582, 46.607, 53.690, 71.054, 65.195, 27.184, 19.131, 0.0000, 19.131, 65.195, 76.694, 62.720, 49.631, 101.57, 112.95, 97.035, 86.164, 130.37, 146.83, 133.39, 119.84, 131.30, 148.14, 154.77, 135.66, 103.31, 114.81, 130.37, 117.28, 67.049, 78.429, 93.836, 82.965, 33.166, 49.631, 60.157, 46.607, 31.860, 48.703, 44.341, 25.235, 126.31, 133.39, 114.81, 108.95, 160.87, 152.82, 160.87, 180.00, 93.279, 103.31, 119.84, 108.95, 72.817, 91.582, 88.418, 69.473},
        {62.720, 76.694, 65.195, 49.631, 69.473, 88.418, 91.582, 72.817, 33.166, 46.607, 60.157, 49.631, 19.131, 27.184, 19.131, 0.0000, 71.054, 86.721, 76.694, 60.157, 97.035, 112.95, 101.57, 86.164, 114.81, 133.39, 126.31, 108.95, 113.11, 131.30, 135.66, 116.60, 93.279, 108.95, 119.84, 103.31, 67.049, 82.965, 93.836, 78.429, 46.607, 65.195, 71.054, 53.690, 48.703, 66.891, 63.403, 44.341, 133.39, 146.83, 130.37, 119.84, 152.82, 160.87, 180.00, 160.87, 103.31, 117.28, 130.37, 114.81, 91.582, 110.53, 107.18, 88.418},
        {60.157, 49.631, 33.166, 46.607, 108.95, 119.84, 103.31, 93.279, 103.31, 117.28, 130.37, 114.81, 53.690, 46.607, 65.195, 71.054, 0.0000, 19.131, 27.184, 19.131, 44.341, 48.703, 31.860, 25.235, 88.418, 91.582, 72.817, 69.473, 131.30, 135.66, 116.60, 113.11, 160.87, 180.00, 160.87, 152.82, 131.30, 135.66, 154.77, 148.14, 88.418, 91.582, 110.53, 107.18, 44.341, 48.703, 66.891, 63.403, 62.720, 76.694, 65.195, 49.631, 133.39, 126.31, 108.95, 114.81, 130.37, 119.84, 133.39, 146.83, 60.157, 71.054, 86.721, 76.694},
        {76.694, 62.720, 49.631, 65.195, 126.31, 133.39, 114.81, 108.95, 119.84, 130.37, 146.83, 133.39, 71.054, 60.157, 76.694, 86.721, 19.131, 0.0000, 19.131, 27.184, 48.703, 44.341, 25.235, 31.860, 91.582, 88.418, 69.473, 72.817, 135.66, 131.30, 113.11, 116.60, 180.00, 160.87, 152.82, 160.87, 135.66, 131.30, 148.14, 154.77, 91.582, 88.418, 107.18, 110.53, 48.703, 44.341, 63.403, 66.891, 49.631, 60.157, 46.607, 33.166, 119.84, 108.95, 93.279, 103.31, 117.28, 103.31, 114.81, 130.37, 46.607, 53.690, 71.054, 65.195},
        {86.721, 76.694, 60.157, 71.054, 133.39, 146.83, 130.37, 119.84, 108.95, 114.81, 133.39, 126.31, 65.195, 49.631, 62.720, 76.694, 27.184, 19.131, 0.0000, 19.131, 66.891, 63.403, 44.341, 48.703, 110.53, 107.18, 88.418, 91.582, 154.77, 148.14, 131.30, 135.66, 160.87, 152.82, 160.87, 180.00, 116.60, 113.11, 131.30, 135.66, 72.817, 69.473, 88.418, 91.582, 31.860, 25.235, 44.341, 48.703, 65.195, 71.054, 53.690, 46.607, 130.37, 114.81, 103.31, 117.28, 103.31, 93.279, 108.95, 119.84, 33.166, 46.607, 60.157, 49.631},
        {71.054, 65.195, 46.607, 53.690, 114.81, 130.37, 117.28, 103.31, 93.279, 103.31, 119.84, 108.95, 46.607, 33.166, 49.631, 60.157, 19.131, 27.184, 19.131, 0.0000, 63.403, 66.891, 48.703, 44.341, 107.18, 110.53, 91.582, 88.418, 148.14, 154.77, 135.66, 131.30, 152.82, 160.87, 180.00, 160.87, 113.11, 116.60, 135.66, 131.30, 69.473, 72.817, 91.582, 88.418, 25.235, 31.860, 48.703, 44.341, 76.694, 86.721, 71.054, 60.157, 146.83, 133.39, 119.84, 130.37, 114.81, 108.95, 126.31, 133.39, 49.631, 65.195, 76.694, 62.720},
        {44.341, 25.235, 31.860, 48.703, 86.164, 86.164, 67.049, 67.049, 116.60, 135.66, 131.30, 113.11, 78.429, 82.965, 101.57, 97.035, 44.341, 48.703, 66.891, 63.403, 0.0000, 19.131, 27.184, 19.131, 44.341, 48.703, 31.860, 25.235, 88.418, 91.582, 72.817, 69.473, 131.30, 135.66, 116.60, 113.11, 160.87, 180.00, 160.87, 152.82, 131.30, 135.66, 154.77, 148.14, 88.418, 91.582, 110.53, 107.18, 44.341, 63.403, 66.891, 48.703, 97.035, 101.57, 82.965, 78.429, 154.77, 135.66, 131.30, 148.14, 93.836, 93.836, 112.95, 112.95},
        {63.403, 44.341, 48.703, 66.891, 101.57, 97.035, 78.429, 82.965, 135.66, 154.77, 148.14, 131.30, 93.836, 93.836, 112.95, 112.95, 48.703, 44.341, 63.403, 66.891, 19.131, 0.0000, 19.131, 27.184, 48.703, 44.341, 25.235, 31.860, 91.582, 88.418, 69.473, 72.817, 135.66, 131.30, 113.11, 116.60, 180.00, 160.87, 152.82, 160.87, 135.66, 131.30, 148.14, 154.77, 91.582, 88.418, 107.18, 110.53, 25.235, 44.341, 48.703, 31.860, 86.164, 86.164, 67.049, 67.049, 135.66, 116.60, 113.11, 131.30, 82.965, 78.429, 97.035, 101.57},
        {66.891, 48.703, 44.341, 63.403, 112.95, 112.95, 93.836, 93.836, 131.30, 148.14, 154.77, 135.66, 82.965, 78.429, 97.035, 101.57, 31.860, 25.235, 44.341, 48.703, 27.184, 19.131, 0.0000, 19.131, 66.891, 63.403, 44.341, 48.703, 110.53, 107.18, 88.418, 91.582, 154.77, 148.14, 131.30, 135.66, 160.87, 152.82, 160.87, 180.00, 116.60, 113.11, 131.30, 135.66, 72.817, 69.473, 88.418, 91.582, 31.860, 48.703, 44.341, 25.235, 101.57, 97.035, 78.429, 82.965, 131.30, 113.11, 116.60, 135.66, 67.049, 67.049, 86.164, 86.164},
        {48.703, 31.860, 25.235, 44.341, 97.035, 101.57, 82.965, 78.429, 113.11, 131.30, 135.66, 116.60, 67.049, 67.049, 86.164, 86.164, 25.235, 31.860, 48.703, 44.341, 19.131, 27.184, 19.131, 0.0000, 63.403, 66.891, 48.703, 44.341, 107.18, 110.53, 91.582, 88.418, 148.14, 154.77, 135.66, 131.30, 152.82, 160.87, 180.00, 160.87, 113.11, 116.60, 135.66, 131.30, 69.473, 72.817, 91.582, 88.418, 48.703, 66.891, 63.403, 44.341, 112.95, 112.95, 93.836, 93.836, 148.14, 131.30, 135.66, 154.77, 78.429, 82.965, 101.57, 97.035},
        {53.690, 46.607, 65.195, 71.054, 60.157, 49.631, 33.166, 46.607, 108.95, 119.84, 103.31, 93.279, 103.31, 117.28, 130.37, 114.81, 88.418, 91.582, 110.53, 107.18, 44.341, 48.703, 66.891, 63.403, 0.0000, 19.131, 27.184, 19.131, 44.341, 48.703, 31.860, 25.235, 88.418, 91.582, 72.817, 69.473, 131.30, 135.66, 116.60, 113.11, 160.87, 180.00, 160.87, 152.82, 131.30, 135.66, 154.77, 148.14, 60.157, 71.054, 86.721, 76.694, 62.720, 76.694, 65.195, 49.631, 133.39, 126.31, 108.95, 114.81, 130.37, 119.84, 133.39, 146.83},
        {71.054, 60.157, 76.694, 86.721, 76.694, 62.720, 49.631, 65.195, 126.31, 133.39, 114.81, 108.95, 119.84, 130.37, 146.83, 133.39, 91.582, 88.418, 107.18, 110.53, 48.703, 44.341, 63.403, 66.891, 19.131, 0.0000, 19.131, 27.184, 48.703, 44.341, 25.235, 31.860, 91.582, 88.418, 69.473, 72.817, 135.66, 131.30, 113.11, 116.60, 180.00, 160.87, 152.82, 160.87, 135.66, 131.30, 148.14, 154.77, 46.607, 53.690, 71.054, 65.195, 49.631, 60.157, 46.607, 33.166, 119.84, 108.95, 93.279, 103.31, 117.28, 103.31, 114.81, 130.37},
        {65.195, 49.631, 62.720, 76.694, 86.721, 76.694, 60.157, 71.054, 133.39, 146.83, 130.37, 119.84, 108.95, 114.81, 133.39, 126.31, 72.817, 69.473, 88.418, 91.582, 31.860, 25.235, 44.341, 48.703, 27.184, 19.131, 0.0000, 19.131, 66.891, 63.403, 44.341, 48.703, 110.53, 107.18, 88.418, 91.582, 154.77, 148.14, 131.30, 135.66, 160.87, 152.82, 160.87, 180.00, 116.60, 113.11, 131.30, 135.66, 33.166, 46.607, 60.157, 49.631, 65.195, 71.054, 53.690, 46.607, 130.37, 114.81, 103.31, 117.28, 103.31, 93.279, 108.95, 119.84},
        {46.607, 33.166, 49.631, 60.157, 71.054, 65.195, 46.607, 53.690, 114.81, 130.37, 117.28, 103.31, 93.279, 103.31, 119.84, 108.95, 69.473, 72.817, 91.582, 88.418, 25.235, 31.860, 48.703, 44.341, 19.131, 27.184, 19.131, 0.0000, 63.403, 66.891, 48.703, 44.341, 107.18, 110.53, 91.582, 88.418, 148.14, 154.77, 135.66, 131.30, 152.82, 160.87, 180.00, 160.87, 113.11, 116.60, 135.66, 131.30, 49.631, 65.195, 76.694, 62.720, 76.694, 86.721, 71.054, 60.157, 146.83, 133.39, 119.84, 130.37, 114.81, 108.95, 126.31, 133.39},
        {78.429, 82.965, 101.57, 97.035, 44.341, 25.235, 31.860, 48.703, 86.164, 86.164, 67.049, 67.049, 116.60, 135.66, 131.30, 113.11, 131.30, 135.66, 154.77, 148.14, 88.418, 91.582, 110.53, 107.18, 44.341, 48.703, 66.891, 63.403, 0.0000, 19.131, 27.184, 19.131, 44.341, 48.703, 31.860, 25.235, 88.418, 91.582, 72.817, 69.473, 131.30, 135.66, 116.60, 113.11, 160.87, 180.00, 160.87, 152.82, 93.836, 93.836, 112.95, 112.95, 44.341, 63.403, 66.891, 48.703, 97.035, 101.57, 82.965, 78.429, 154.77, 135.66, 131.30, 148.14},
        {93.836, 93.836, 112.95, 112.95, 63.403, 44.341, 48.703, 66.891, 101.57, 97.035, 78.429, 82.965, 135.66, 154.77, 148.14, 131.30, 135.66, 131.30, 148.14, 154.77, 91.582, 88.418, 107.18, 110.53, 48.703, 44.341, 63.403, 66.891, 19.131, 0.0000, 19.131, 27.184, 48.703, 44.341, 25.235, 31.860, 91.582, 88.418, 69.473, 72.817, 135.66, 131.30, 113.11, 116.60, 180.00, 160.87, 152.82, 160.87, 82.965, 78.429, 97.035, 101.57, 25.235, 44.341, 48.703, 31.860, 86.164, 86.164, 67.049, 67.049, 135.66, 116.60, 113.11, 131.30},
        {82.965, 78.429, 97.035, 101.57, 66.891, 48.703, 44.341, 63.403, 112.95, 112.95, 93.836, 93.836, 131.30, 148.14, 154.77, 135.66, 116.60, 113.11, 131.30, 135.66, 72.817, 69.473, 88.418, 91.582, 31.860, 25.235, 44.341, 48.703, 27.184, 19.131, 0.0000, 19.131, 66.891, 63.403, 44.341, 48.703, 110.53, 107.18, 88.418, 91.582, 154.77, 148.14, 131.30, 135.66, 160.87, 152.82, 160.87, 180.00, 67.049, 67.049, 86.164, 86.164, 31.860, 48.703, 44.341, 25.235, 101.57, 97.035, 78.429, 82.965, 131.30, 113.11, 116.60, 135.66},
        {67.049, 67.049, 86.164, 86.164, 48.703, 31.860, 25.235, 44.341, 97.035, 101.57, 82.965, 78.429, 113.11, 131.30, 135.66, 116.60, 113.11, 116.60, 135.66, 131.30, 69.473, 72.817, 91.582, 88.418, 25.235, 31.860, 48.703, 44.341, 19.131, 27.184, 19.131, 0.0000, 63.403, 66.891, 48.703, 44.341, 107.18, 110.53, 91.582, 88.418, 148.14, 154.77, 135.66, 131.30, 152.82, 160.87, 180.00, 160.87, 78.429, 82.965, 101.57, 97.035, 48.703, 66.891, 63.403, 44.341, 112.95, 112.95, 93.836, 93.836, 148.14, 131.30, 135.66, 154.77},
        {103.31, 117.28, 130.37, 114.81, 53.690, 46.607, 65.195, 71.054, 60.157, 49.631, 33.166, 46.607, 108.95, 119.84, 103.31, 93.279, 160.87, 180.00, 160.87, 152.82, 131.30, 135.66, 154.77, 148.14, 88.418, 91.582, 110.53, 107.18, 44.341, 48.703, 66.891, 63.403, 0.0000, 19.131, 27.184, 19.131, 44.341, 48.703, 31.860, 25.235, 88.418, 91.582, 72.817, 69.473, 131.30, 135.66, 116.60, 113.11, 130.37, 119.84, 133.39, 146.83, 60.157, 71.054, 86.721, 76.694, 62.720, 76.694, 65.195, 49.631, 133.39, 126.31, 108.95, 114.81},
        {119.84, 130.37, 146.83, 133.39, 71.054, 60.157, 76.694, 86.721, 76.694, 62.720, 49.631, 65.195, 126.31, 133.39, 114.81, 108.95, 180.00, 160.87, 152.82, 160.87, 135.66, 131.30, 148.14, 154.77, 91.582, 88.418, 107.18, 110.53, 48.703, 44.341, 63.403, 66.891, 19.131, 0.0000, 19.131, 27.184, 48.703, 44.341, 25.235, 31.860, 91.582, 88.418, 69.473, 72.817, 135.66, 131.30, 113.11, 116.60, 117.28, 103.31, 114.81, 130.37, 46.607, 53.690, 71.054, 65.195, 49.631, 60.157, 46.607, 33.166, 119.84, 108.95, 93.279, 103.31},
        {108.95, 114.81, 133.39, 126.31, 65.195, 49.631, 62.720, 76.694, 86.721, 76.694, 60.157, 71.054, 133.39, 146.83, 130.37, 119.84, 160.87, 152.82, 160.87, 180.00, 116.60, 113.11, 131.30, 135.66, 72.817, 69.473, 88.418, 91.582, 31.860, 25.235, 44.341, 48.703, 27.184, 19.131, 0.0000, 19.131, 66.891, 63.403, 44.341, 48.703, 110.53, 107.18, 88.418, 91.582, 154.77, 148.14, 131.30, 135.66, 103.31, 93.279, 108.95, 119.84, 33.166, 46.607, 60.157, 49.631, 65.195, 71.054, 53.690, 46.607, 130.37, 114.81, 103.31, 117.28},
        {93.279, 103.31, 119.84, 108.95, 46.607, 33.166, 49.631, 60.157, 71.054, 65.195, 46.607, 53.690, 114.81, 130.37, 117.28, 103.31, 152.82, 160.87, 180.00, 160.87, 113.11, 116.60, 135.66, 131.30, 69.473, 72.817, 91.582, 88.418, 25.235, 31.860, 48.703, 44.341, 19.131, 27.184, 19.131, 0.0000, 63.403, 66.891, 48.703, 44.341, 107.18, 110.53, 91.582, 88.418, 148.14, 154.77, 135.66, 131.30, 114.81, 108.95, 126.31, 133.39, 49.631, 65.195, 76.694, 62.720, 76.694, 86.721, 71.054, 60.157, 146.83, 133.39, 119.84, 130.37},
        {116.60, 135.66, 131.30, 113.11, 78.429, 82.965, 101.57, 97.035, 44.341, 25.235, 31.860, 48.703, 86.164, 86.164, 67.049, 67.049, 131.30, 135.66, 116.60, 113.11, 160.87, 180.00, 160.87, 152.82, 131.30, 135.66, 154.77, 148.14, 88.418, 91.582, 110.53, 107.18, 44.341, 48.703, 66.891, 63.403, 0.0000, 19.131, 27.184, 19.131, 44.341, 48.703, 31.860, 25.235, 88.418, 91.582, 72.817, 69.473, 154.77, 135.66, 131.30, 148.14, 93.836, 93.836, 112.95, 112.95, 44.341, 63.403, 66.891, 48.703, 97.035, 101.57, 82.965, 78.429},
        {135.66, 154.77, 148.14, 131.30, 93.836, 93.836, 112.95, 112.95, 63.403, 44.341, 48.703, 66.891, 101.57, 97.035, 78.429, 82.965, 135.66, 131.30, 113.11, 116.60, 180.00, 160.87, 152.82, 160.87, 135.66, 131.30, 148.14, 154.77, 91.582, 88.418, 107.18, 110.53, 48.703, 44.341, 63.403, 66.891, 19.131, 0.0000, 19.131, 27.184, 48.703, 44.341, 25.235, 31.860, 91.582, 88.418, 69.473, 72.817, 135.66, 116.60, 113.11, 131.30, 82.965, 78.429, 97.035, 101.57, 25.235, 44.341, 48.703, 31.860, 86.164, 86.164, 67.049, 67.049},
        {131.30, 148.14, 154.77, 135.66, 82.965, 78.429, 97.035, 101.57, 66.891, 48.703, 44.341, 63.403, 112.95, 112.95, 93.836, 93.836, 154.77, 148.14, 131.30, 135.66, 160.87, 152.82, 160.87, 180.00, 116.60, 113.11, 131.30, 135.66, 72.817, 69.473, 88.418, 91.582, 31.860, 25.235, 44.341, 48.703, 27.184, 19.131, 0.0000, 19.131, 66.891, 63.403, 44.341, 48.703, 110.53, 107.18, 88.418, 91.582, 131.30, 113.11, 116.60, 135.66, 67.049, 67.049, 86.164, 86.164, 31.860, 48.703, 44.341, 25.235, 101.57, 97.035, 78.429, 82.965},
        {113.11, 131.30, 135.66, 116.60, 67.049, 67.049, 86.164, 86.164, 48.703, 31.860, 25.235, 44.341, 97.035, 101.57, 82.965, 78.429, 148.14, 154.77, 135.66, 131.30, 152.82, 160.87, 180.00, 160.87, 113.11, 116.60, 135.66, 131.30, 69.473, 72.817, 91.582, 88.418, 25.235, 31.860, 48.703, 44.341, 19.131, 27.184, 19.131, 0.0000, 63.403, 66.891, 48.703, 44.341, 107.18, 110.53, 91.582, 88.418, 148.14, 131.30, 135.66, 154.77, 78.429, 82.965, 101.57, 97.035, 48.703, 66.891, 63.403, 44.341, 112.95, 112.95, 93.836, 93.836},
        {108.95, 119.84, 103.31, 93.279, 103.31, 117.28, 130.37, 114.81, 53.690, 46.607, 65.195, 71.054, 60.157, 49.631, 33.166, 46.607, 88.418, 91.582, 72.817, 69.473, 131.30, 135.66, 116.60, 113.11, 160.87, 180.00, 160.87, 152.82, 131.30, 135.66, 154.77, 148.14, 88.418, 91.582, 110.53, 107.18, 44.341, 48.703, 66.891, 63.403, 0.0000, 19.131, 27.184, 19.131, 44.341, 48.703, 31.860, 25.235, 133.39, 126.31, 108.95, 114.81, 130.37, 119.84, 133.39, 146.83, 60.157, 71.054, 86.721, 76.694, 62.720, 76.694, 65.195, 49.631},
        {126.31, 133.39, 114.81, 108.95, 119.84, 130.37, 146.83, 133.39, 71.054, 60.157, 76.694, 86.721, 76.694, 62.720, 49.631, 65.195, 91.582, 88.418, 69.473, 72.817, 135.66, 131.30, 113.11, 116.60, 180.00, 160.87, 152.82, 160.87, 135.66, 131.30, 148.14, 154.77, 91.582, 88.418, 107.18, 110.53, 48.703, 44.341, 63.403, 66.891, 19.131, 0.0000, 19.131, 27.184, 48.703, 44.341, 25.235, 31.860, 119.84, 108.95, 93.279, 103.31, 117.28, 103.31, 114.81, 130.37, 46.607, 53.690, 71.054, 65.195, 49.631, 60.157, 46.607, 33.166},
        {133.39, 146.83, 130.37, 119.84, 108.95, 114.81, 133.39, 126.31, 65.195, 49.631, 62.720, 76.694, 86.721, 76.694, 60.157, 71.054, 110.53, 107.18, 88.418, 91.582, 154.77, 148.14, 131.30, 135.66, 160.87, 152.82, 160.87, 180.00, 116.60, 113.11, 131.30, 135.66, 72.817, 69.473, 88.418, 91.582, 31.860, 25.235, 44.341, 48.703, 27.184, 19.131, 0.0000, 19.131, 66.891, 63.403, 44.341, 48.703, 130.37, 114.81, 103.31, 117.28, 103.31, 93.279, 108.95, 119.84, 33.166, 46.607, 60.157, 49.631, 65.195, 71.054, 53.690, 46.607},
        {114.81, 130.37, 117.28, 103.31, 93.279, 103.31, 119.84, 108.95, 46.607, 33.166, 49.631, 60.157, 71.054, 65.195, 46.607, 53.690, 107.18, 110.53, 91.582, 88.418, 148.14, 154.77, 135.66, 131.30, 152.82, 160.87, 180.00, 160.87, 113.11, 116.60, 135.66, 131.30, 69.473, 72.817, 91.582, 88.418, 25.235, 31.860, 48.703, 44.341, 19.131, 27.184, 19.131, 0.0000, 63.403, 66.891, 48.703, 44.341, 146.83, 133.39, 119.84, 130.37, 114.81, 108.95, 126.31, 133.39, 49.631, 65.195, 76.694, 62.720, 76.694, 86.721, 71.054, 60.157},
        {86.164, 86.164, 67.049, 67.049, 116.60, 135.66, 131.30, 113.11, 78.429, 82.965, 101.57, 97.035, 44.341, 25.235, 31.860, 48.703, 44.341, 48.703, 31.860, 25.235, 88.418, 91.582, 72.817, 69.473, 131.30, 135.66, 116.60, 113.11, 160.87, 180.00, 160.87, 152.82, 131.30, 135.66, 154.77, 148.14, 88.418, 91.582, 110.53, 107.18, 44.341, 48.703, 66.891, 63.403, 0.0000, 19.131, 27.184, 19.131, 97.035, 101.57, 82.965, 78.429, 154.77, 135.66, 131.30, 148.14, 93.836, 93.836, 112.95, 112.95, 44.341, 63.403, 66.891, 48.703},
        {101.57, 97.035, 78.429, 82.965, 135.66, 154.77, 148.14, 131.30, 93.836, 93.836, 112.95, 112.95, 63.403, 44.341, 48.703, 66.891, 48.703, 44.341, 25.235, 31.860, 91.582, 88.418, 69.473, 72.817, 135.66, 131.30, 113.11, 116.60, 180.00, 160.87, 152.82, 160.87, 135.66, 131.30, 148.14, 154.77, 91.582, 88.418, 107.18, 110.53, 48.703, 44.341, 63.403, 66.891, 19.131, 0.0000, 19.131, 27.184, 86.164, 86.164, 67.049, 67.049, 135.66, 116.60, 113.11, 131.30, 82.965, 78.429, 97.035, 101.57, 25.235, 44.341, 48.703, 31.860},
        {112.95, 112.95, 93.836, 93.836, 131.30, 148.14, 154.77, 135.66, 82.965, 78.429, 97.035, 101.57, 66.891, 48.703, 44.341, 63.403, 66.891, 63.403, 44.341, 48.703, 110.53, 107.18, 88.418, 91.582, 154.77, 148.14, 131.30, 135.66, 160.87, 152.82, 160.87, 180.00, 116.60, 113.11, 131.30, 135.66, 72.817, 69.473, 88.418, 91.582, 31.860, 25.235, 44.341, 48.703, 27.184, 19.131, 0.0000, 19.131, 101.57, 97.035, 78.429, 82.965, 131.30, 113.11, 116.60, 135.66, 67.049, 67.049, 86.164, 86.164, 31.860, 48.703, 44.341, 25.235},
        {97.035, 101.57, 82.965, 78.429, 113.11, 131.30, 135.66, 116.60, 67.049, 67.049, 86.164, 86.164, 48.703, 31.860, 25.235, 44.341, 63.403, 66.891, 48.703, 44.341, 107.18, 110.53, 91.582, 88.418, 148.14, 154.77, 135.66, 131.30, 152.82, 160.87, 180.00, 160.87, 113.11, 116.60, 135.66, 131.30, 69.473, 72.817, 91.582, 88.418, 25.235, 31.860, 48.703, 44.341, 19.131, 27.184, 19.131, 0.0000, 112.95, 112.95, 93.836, 93.836, 148.14, 131.30, 135.66, 154.77, 78.429, 82.965, 101.57, 97.035, 48.703, 66.891, 63.403, 44.341},
        {88.418, 69.473, 72.817, 91.582, 119.84, 108.95, 93.279, 103.31, 160.87, 180.00, 160.87, 152.82, 114.81, 108.95, 126.31, 133.39, 62.720, 49.631, 65.195, 76.694, 44.341, 25.235, 31.860, 48.703, 60.157, 46.607, 33.166, 49.631, 93.836, 82.965, 67.049, 78.429, 130.37, 117.28, 103.31, 114.81, 154.77, 135.66, 131.30, 148.14, 133.39, 119.84, 130.37, 146.83, 97.035, 86.164, 101.57, 112.95, 0.0000, 19.131, 27.184, 19.131, 71.054, 65.195, 46.607, 53.690, 110.53, 91.582, 88.418, 107.18, 71.054, 60.157, 76.694, 86.721},
        {107.18, 88.418, 91.582, 110.53, 130.37, 114.81, 103.31, 117.28, 180.00, 160.87, 152.82, 160.87, 130.37, 119.84, 133.39, 146.83, 76.694, 60.157, 71.054, 86.721, 63.403, 44.341, 48.703, 66.891, 71.054, 53.690, 46.607, 65.195, 93.836, 78.429, 67.049, 82.965, 119.84, 103.31, 93.279, 108.95, 135.66, 116.60, 113.11, 131.30, 126.31, 108.95, 114.81, 133.39, 101.57, 86.164, 97.035, 112.95, 19.131, 0.0000, 19.131, 27.184, 60.157, 49.631, 33.166, 46.607, 91.582, 72.817, 69.473, 88.418, 65.195, 49.631, 62.720, 76.694},
        {110.53, 91.582, 88.418, 107.18, 146.83, 133.39, 119.84, 130.37, 160.87, 152.82, 160.87, 180.00, 117.28, 103.31, 114.81, 130.37, 65.195, 46.607, 53.690, 71.054, 66.891, 48.703, 44.341, 63.403, 86.721, 71.054, 60.157, 76.694, 112.95, 97.035, 86.164, 101.57, 133.39, 114.81, 108.95, 126.31, 131.30, 113.11, 116.60, 135.66, 108.95, 93.279, 103.31, 119.84, 82.965, 67.049, 78.429, 93.836, 27.184, 19.131, 0.0000, 19.131, 76.694, 62.720, 49.631, 65.195, 88.418, 69.473, 72.817, 91.582, 46.607, 33.166, 49.631, 60.157},
        {91.582, 72.817, 69.473, 88.418, 133.39, 126.31, 108.95, 114.81, 152.82, 160.87, 180.00, 160.87, 103.31, 93.279, 108.95, 119.84, 49.631, 33.166, 46.607, 60.157, 48.703, 31.860, 25.235, 44.341, 76.694, 65.195, 49.631, 62.720, 112.95, 101.57, 86.164, 97.035, 146.83, 130.37, 119.84, 133.39, 148.14, 131.30, 135.66, 154.77, 114.81, 103.31, 117.28, 130.37, 78.429, 67.049, 82.965, 93.836, 19.131, 27.184, 19.131, 0.0000, 86.721, 76.694, 60.157, 71.054, 107.18, 88.418, 91.582, 110.53, 53.690, 46.607, 65.195, 71.054},
        {114.81, 108.95, 126.31, 133.39, 88.418, 69.473, 72.817, 91.582, 119.84, 108.95, 93.279, 103.31, 160.87, 180.00, 160.87, 152.82, 133.39, 119.84, 130.37, 146.83, 97.035, 86.164, 101.57, 112.95, 62.720, 49.631, 65.195, 76.694, 44.341, 25.235, 31.860, 48.703, 60.157, 46.607, 33.166, 49.631, 93.836, 82.965, 67.049, 78.429, 130.37, 117.28, 103.31, 114.81, 154.77, 135.66, 131.30, 148.14, 71.054, 60.157, 76.694, 86.721, 0.0000, 19.131, 27.184, 19.131, 71.054, 65.195, 46.607, 53.690, 110.53, 91.582, 88.418, 107.18},
        {130.37, 119.84, 133.39, 146.83, 107.18, 88.418, 91.582, 110.53, 130.37, 114.81, 103.31, 117.28, 180.00, 160.87, 152.82, 160.87, 126.31, 108.95, 114.81, 133.39, 101.57, 86.164, 97.035, 112.95, 76.694, 60.157, 71.054, 86.721, 63.403, 44.341, 48.703, 66.891, 71.054, 53.690, 46.607, 65.195, 93.836, 78.429, 67.049, 82.965, 119.84, 103.31, 93.279, 108.95, 135.66, 116.60, 113.11, 131.30, 65.195, 49.631, 62.720, 76.694, 19.131, 0.0000, 19.131, 27.184, 60.157, 49.631, 33.166, 46.607, 91.582, 72.817, 69.473, 88.418},
        {117.28, 103.31, 114.81, 130.37, 110.53, 91.582, 88.418, 107.18, 146.83, 133.39, 119.84, 130.37, 160.87, 152.82, 160.87, 180.00, 108.95, 93.279, 103.31, 119.84, 82.965, 67.049, 78.429, 93.836, 65.195, 46.607, 53.690, 71.054, 66.891, 48.703, 44.341, 63.403, 86.721, 71.054, 60.157, 76.694, 112.95, 97.035, 86.164, 101.57, 133.39, 114.81, 108.95, 126.31, 131.30, 113.11, 116.60, 135.66, 46.607, 33.166, 49.631, 60.157, 27.184, 19.131, 0.0000, 19.131, 76.694, 62.720, 49.631, 65.195, 88.418, 69.473, 72.817, 91.582},
        {103.31, 93.279, 108.95, 119.84, 91.582, 72.817, 69.473, 88.418, 133.39, 126.31, 108.95, 114.81, 152.82, 160.87, 180.00, 160.87, 114.81, 103.31, 117.28, 130.37, 78.429, 67.049, 82.965, 93.836, 49.631, 33.166, 46.607, 60.157, 48.703, 31.860, 25.235, 44.341, 76.694, 65.195, 49.631, 62.720, 112.95, 101.57, 86.164, 97.035, 146.83, 130.37, 119.84, 133.39, 148.14, 131.30, 135.66, 154.77, 53.690, 46.607, 65.195, 71.054, 19.131, 27.184, 19.131, 0.0000, 86.721, 76.694, 60.157, 71.054, 107.18, 88.418, 91.582, 110.53},
        {160.87, 180.00, 160.87, 152.82, 114.81, 108.95, 126.31, 133.39, 88.418, 69.473, 72.817, 91.582, 119.84, 108.95, 93.279, 103.31, 130.37, 117.28, 103.31, 114.81, 154.77, 135.66, 131.30, 148.14, 133.39, 119.84, 130.37, 146.83, 97.035, 86.164, 101.57, 112.95, 62.720, 49.631, 65.195, 76.694, 44.341, 25.235, 31.860, 48.703, 60.157, 46.607, 33.166, 49.631, 93.836, 82.965, 67.049, 78.429, 110.53, 91.582, 88.418, 107.18, 71.054, 60.157, 76.694, 86.721, 0.0000, 19.131, 27.184, 19.131, 71.054, 65.195, 46.607, 53.690},
        {180.00, 160.87, 152.82, 160.87, 130.37, 119.84, 133.39, 146.83, 107.18, 88.418, 91.582, 110.53, 130.37, 114.81, 103.31, 117.28, 119.84, 103.31, 93.279, 108.95, 135.66, 116.60, 113.11, 131.30, 126.31, 108.95, 114.81, 133.39, 101.57, 86.164, 97.035, 112.95, 76.694, 60.157, 71.054, 86.721, 63.403, 44.341, 48.703, 66.891, 71.054, 53.690, 46.607, 65.195, 93.836, 78.429, 67.049, 82.965, 91.582, 72.817, 69.473, 88.418, 65.195, 49.631, 62.720, 76.694, 19.131, 0.0000, 19.131, 27.184, 60.157, 49.631, 33.166, 46.607},
        {160.87, 152.82, 160.87, 180.00, 117.28, 103.31, 114.81, 130.37, 110.53, 91.582, 88.418, 107.18, 146.83, 133.39, 119.84, 130.37, 133.39, 114.81, 108.95, 126.31, 131.30, 113.11, 116.60, 135.66, 108.95, 93.279, 103.31, 119.84, 82.965, 67.049, 78.429, 93.836, 65.195, 46.607, 53.690, 71.054, 66.891, 48.703, 44.341, 63.403, 86.721, 71.054, 60.157, 76.694, 112.95, 97.035, 86.164, 101.57, 88.418, 69.473, 72.817, 91.582, 46.607, 33.166, 49.631, 60.157, 27.184, 19.131, 0.0000, 19.131, 76.694, 62.720, 49.631, 65.195},
        {152.82, 160.87, 180.00, 160.87, 103.31, 93.279, 108.95, 119.84, 91.582, 72.817, 69.473, 88.418, 133.39, 126.31, 108.95, 114.81, 146.83, 130.37, 119.84, 133.39, 148.14, 131.30, 135.66, 154.77, 114.81, 103.31, 117.28, 130.37, 78.429, 67.049, 82.965, 93.836, 49.631, 33.166, 46.607, 60.157, 48.703, 31.860, 25.235, 44.341, 76.694, 65.195, 49.631, 62.720, 112.95, 101.57, 86.164, 97.035, 107.18, 88.418, 91.582, 110.53, 53.690, 46.607, 65.195, 71.054, 19.131, 27.184, 19.131, 0.0000, 86.721, 76.694, 60.157, 71.054},
        {119.84, 108.95, 93.279, 103.31, 160.87, 180.00, 160.87, 152.82, 114.81, 108.95, 126.31, 133.39, 88.418, 69.473, 72.817, 91.582, 60.157, 46.607, 33.166, 49.631, 93.836, 82.965, 67.049, 78.429, 130.37, 117.28, 103.31, 114.81, 154.77, 135.66, 131.30, 148.14, 133.39, 119.84, 130.37, 146.83, 97.035, 86.164, 101.57, 112.95, 62.720, 49.631, 65.195, 76.694, 44.341, 25.235, 31.860, 48.703, 71.054, 65.195, 46.607, 53.690, 110.53, 91.582, 88.418, 107.18, 71.054, 60.157, 76.694, 86.721, 0.0000, 19.131, 27.184, 19.131},
        {130.37, 114.81, 103.31, 117.28, 180.00, 160.87, 152.82, 160.87, 130.37, 119.84, 133.39, 146.83, 107.18, 88.418, 91.582, 110.53, 71.054, 53.690, 46.607, 65.195, 93.836, 78.429, 67.049, 82.965, 119.84, 103.31, 93.279, 108.95, 135.66, 116.60, 113.11, 131.30, 126.31, 108.95, 114.81, 133.39, 101.57, 86.164, 97.035, 112.95, 76.694, 60.157, 71.054, 86.721, 63.403, 44.341, 48.703, 66.891, 60.157, 49.631, 33.166, 46.607, 91.582, 72.817, 69.473, 88.418, 65.195, 49.631, 62.720, 76.694, 19.131, 0.0000, 19.131, 27.184},
        {146.83, 133.39, 119.84, 130.37, 160.87, 152.82, 160.87, 180.00, 117.28, 103.31, 114.81, 130.37, 110.53, 91.582, 88.418, 107.18, 86.721, 71.054, 60.157, 76.694, 112.95, 97.035, 86.164, 101.57, 133.39, 114.81, 108.95, 126.31, 131.30, 113.11, 116.60, 135.66, 108.95, 93.279, 103.31, 119.84, 82.965, 67.049, 78.429, 93.836, 65.195, 46.607, 53.690, 71.054, 66.891, 48.703, 44.341, 63.403, 76.694, 62.720, 49.631, 65.195, 88.418, 69.473, 72.817, 91.582, 46.607, 33.166, 49.631, 60.157, 27.184, 19.131, 0.0000, 19.131},
        {133.39, 126.31, 108.95, 114.81, 152.82, 160.87, 180.00, 160.87, 103.31, 93.279, 108.95, 119.84, 91.582, 72.817, 69.473, 88.418, 76.694, 65.195, 49.631, 62.720, 112.95, 101.57, 86.164, 97.035, 146.83, 130.37, 119.84, 133.39, 148.14, 131.30, 135.66, 154.77, 114.81, 103.31, 117.28, 130.37, 78.429, 67.049, 82.965, 93.836, 49.631, 33.166, 46.607, 60.157, 48.703, 31.860, 25.235, 44.341, 86.721, 76.694, 60.157, 71.054, 107.18, 88.418, 91.582, 110.53, 53.690, 46.607, 65.195, 71.054, 19.131, 27.184, 19.131, 0}
    };
    memcpy(GriffinCryMap, thisGriffinCryMap, sizeof(GriffinCryMap));

    double thisGriffinCryMapCombos[52][2] = {
        {0.0000, 64},
        {19.131, 128},
        {25.235, 64},
        {27.184, 64},
        {31.860, 64},
        {33.166, 48},
        {44.341, 128},
        {46.607, 96},
        {48.703, 128},
        {49.631, 96},
        {53.690, 48},
        {60.157, 96},
        {62.720, 48},
        {63.403, 64},
        {65.195, 96},
        {66.891, 64},
        {67.049, 64},
        {69.473, 64},
        {71.054, 96},
        {72.817, 64},
        {76.694, 96},
        {78.429, 64},
        {82.965, 64},
        {86.164, 64},
        {86.721, 48},
        {88.418, 128},
        {91.582, 128},
        {93.279, 48},
        {93.836, 64},
        {97.035, 64},
        {101.57, 64},
        {103.31, 96},
        {107.18, 64},
        {108.95, 96},
        {110.53, 64},
        {112.95, 64},
        {113.11, 64},
        {114.81, 96},
        {116.60, 64},
        {117.28, 48},
        {119.84, 96},
        {126.31, 48},
        {130.37, 96},
        {131.30, 128},
        {133.39, 96},
        {135.66, 128},
        {146.83, 48},
        {148.14, 64},
        {152.82, 64},
        {154.77, 64},
        {160.87, 128},
        {180.00, 64}
    };
    memcpy(GriffinCryMapCombos, thisGriffinCryMapCombos, sizeof(GriffinCryMapCombos));

    //----------------------------------------------------------------------------------------------------

    if(fSettings->WriteGriffinAddbackVector()) {
        for(int i = 0; i < 4; i++) {
            for(int j = 0; j < 16; j++) {
                GriffinCrystalCenterVectors[i+(j*4)] = GriffinCrystalCenterPosition(i,j);
                //std::cout << "cry = " << i << " det = " << j << " : x = " << GriffinCrystalCenterVectors[i+(j*4)].X() << " mm - y = " << GriffinCrystalCenterVectors[i+(j*4)].Y() << " mm - z = " << GriffinCrystalCenterVectors[i+(j*4)].Z() << " mm" << std::endl;
            }
        }
    }

    fSceptarHit = false;

    //add branches to input chain
    fChain.SetBranchAddress("eventNumber", &fEventNumber);
    fChain.SetBranchAddress("trackID", &fTrackID);
    fChain.SetBranchAddress("parentID", &fParentID);
    fChain.SetBranchAddress("stepNumber", &fStepNumber);
    fChain.SetBranchAddress("particleType", &fParticleType);
    fChain.SetBranchAddress("processType", &fProcessType);
    fChain.SetBranchAddress("systemID", &fSystemID);
    fChain.SetBranchAddress("detNumber", &fDetNumber);
    fChain.SetBranchAddress("cryNumber", &fCryNumber);
    fChain.SetBranchAddress("depEnergy", &fDepEnergy);
    fChain.SetBranchAddress("posx", &fPosx);
    fChain.SetBranchAddress("posy", &fPosy);
    fChain.SetBranchAddress("posz", &fPosz);
    fChain.SetBranchAddress("time", &fTime);
    fChain.SetBranchAddress("targetZ", &fTargetZ);
    fChain.SetBranchAddress("targetA", &fTargetA);

    // add branches from the TRex derived generators
    // treeGen
    fTISTARGenChain.SetBranchAddress("reactionEnergy",      &fTISTARGenReactionBeamEnergy);
    fTISTARGenChain.SetBranchAddress("reactionEnergyCM",    &fTISTARGenReactionBeamEnergyCM);
    fTISTARGenChain.SetBranchAddress("reactionX",           &fTISTARGenReactionX);
    fTISTARGenChain.SetBranchAddress("reactionY",           &fTISTARGenReactionY);
    fTISTARGenChain.SetBranchAddress("reactionZ",           &fTISTARGenReactionZ);
    fTISTARGenChain.SetBranchAddress("recoilTheta",         &fTISTARGenRecoilTheta);
    fTISTARGenChain.SetBranchAddress("recoilPhi",           &fTISTARGenRecoilPhi);
    fTISTARGenChain.SetBranchAddress("recoilEnergy",        &fTISTARGenRecoilEnergy);
    fTISTARGenChain.SetBranchAddress("ejectileTheta",       &fTISTARGenEjectileTheta);
    fTISTARGenChain.SetBranchAddress("ejectilePhi",         &fTISTARGenEjectilePhi);
    fTISTARGenChain.SetBranchAddress("ejectileEnergy",      &fTISTARGenEjectileEnergy);
    fTISTARGenChain.SetBranchAddress("reaction",            &fTISTARGenReaction);
    
    TBranch * branchGammaEnergy = 0;
    TBranch * branchGammaTheta = 0;
    TBranch * branchGammaPhi = 0;
    fTISTARGenChain.SetBranchAddress("gammaEnergy", &fTISTARGenGammaEnergy, &branchGammaEnergy);
    fTISTARGenChain.SetBranchAddress("gammaTheta",  &fTISTARGenGammaTheta,  &branchGammaTheta);
    fTISTARGenChain.SetBranchAddress("gammaPhi",    &fTISTARGenGammaPhi,    &branchGammaPhi);

    //create output file
    fOutput = new TFile(outputFileName.c_str(),"recreate");
    if(!fOutput->IsOpen()) {
        std::cerr<<"Failed to open file '"<<outputFileName<<"', check permissions on directory and disk space!"<<std::endl;
        throw;
    }

    //set tree to belong to output file
    if(fSettings->WriteTree())
        fTree.SetDirectory(fOutput);

    //create branches for output tree
    // GRIFFIN
    fGriffinCrystal       = new std::vector<Detector>;
    fGriffinDetector      = new std::vector<Detector>;
    fGriffinNeighbour     = new std::vector<Detector>;
    fGriffinNeighbourVector = new std::vector<Detector>;
    fGriffinArray         = new std::vector<Detector>;
    fGriffinBgo           = new std::vector<Detector>;
    fGriffinBgoBack       = new std::vector<Detector>;
    fTree.Branch("GriffinCrystal",&fGriffinCrystal, fSettings->BufferSize());
    fTree.Branch("GriffinDetector",&fGriffinDetector, fSettings->BufferSize());
    fTree.Branch("GriffinNeighbour",&fGriffinNeighbour, fSettings->BufferSize());
    fTree.Branch("GriffinNeighbourVector",&fGriffinNeighbourVector, fSettings->BufferSize());
    fTree.Branch("GriffinArray",&fGriffinArray, fSettings->BufferSize());
    fTree.Branch("GriffinBgo",&fGriffinBgo, fSettings->BufferSize());
    fTree.Branch("GriffinBgo",&fGriffinBgoBack, fSettings->BufferSize());

    // LaBr
    fLaBrArray            = new std::vector<Detector>;
    fLaBrDetector         = new std::vector<Detector>;
    fTree.Branch("LaBrArray",&fLaBrArray, fSettings->BufferSize());
    fTree.Branch("LaBrDetector",&fLaBrDetector, fSettings->BufferSize());

    // EightPi
    fEightPiArray            = new std::vector<Detector>;
    fEightPiDetector         = new std::vector<Detector>;
    fEightPiBgoDetector         = new std::vector<Detector>;
    fTree.Branch("EightPiArray",&fEightPiArray, fSettings->BufferSize());
    fTree.Branch("EightPiDetector",&fEightPiDetector, fSettings->BufferSize());
    fTree.Branch("EightPiBgoDetector",&fEightPiBgoDetector, fSettings->BufferSize());

    // Ancillary Detector
    fAncillaryBgoCrystal  = new std::vector<Detector>;
    fAncillaryBgoDetector = new std::vector<Detector>;
    fAncillaryBgoArray    = new std::vector<Detector>;
    fTree.Branch("AncillaryBgoCrystal",&fAncillaryBgoCrystal, fSettings->BufferSize());
    fTree.Branch("AncillaryBgoDetector",&fAncillaryBgoDetector, fSettings->BufferSize());
    fTree.Branch("AncillaryBgoArray",&fAncillaryBgoArray, fSettings->BufferSize());

    // SCEPTAR
    fSceptarArray         = new std::vector<Detector>;
    fSceptarDetector      = new std::vector<Detector>;
    fTree.Branch("SceptarArray",&fSceptarArray, fSettings->BufferSize());
    fTree.Branch("SceptarDetector",&fSceptarDetector, fSettings->BufferSize());

    // DESCANT
    fDescantArray            = new std::vector<Detector>;
    fTree.Branch("DescantArray",&fDescantArray, fSettings->BufferSize());
    fDescantBlueDetector = new std::vector<Detector>;
    fDescantGreenDetector = new std::vector<Detector>;
    fDescantRedDetector = new std::vector<Detector>;
    fDescantWhiteDetector = new std::vector<Detector>;
    fDescantYellowDetector = new std::vector<Detector>;
    fTree.Branch("DescantBlueDetector",&fDescantBlueDetector, fSettings->BufferSize());
    fTree.Branch("DescantGreenDetector",&fDescantGreenDetector, fSettings->BufferSize());
    fTree.Branch("DescantRedDetector",&fDescantRedDetector, fSettings->BufferSize());
    fTree.Branch("DescantWhiteDetector",&fDescantWhiteDetector, fSettings->BufferSize());
    fTree.Branch("DescantYellowDetector",&fDescantYellowDetector, fSettings->BufferSize());

    // Testcan
    fTestcanDetector         = new std::vector<Detector>;
    fTree.Branch("TestcanDetector", &fTestcanDetector, fSettings->BufferSize());

    // PACES
    fPacesArray         = new std::vector<Detector>;
    fPacesDetector      = new std::vector<Detector>;
    fTree.Branch("PacesArray",&fPacesArray, fSettings->BufferSize());
    fTree.Branch("PacesDetector",&fPacesDetector, fSettings->BufferSize());

    // TI-STAR
    fTISTARArray        = new std::vector<Detector>;
    fTISTARLayer1       = new std::vector<Detector>;
    fTISTARLayer2       = new std::vector<Detector>;
    fTISTARLayer3       = new std::vector<Detector>;
    fTree.Branch("TISTARArray",         &fTISTARArray,       fSettings->BufferSize());     
    fTree.Branch("TISTARLayer1",        &fTISTARLayer1,      fSettings->BufferSize());     
    fTree.Branch("TISTARLayer2",        &fTISTARLayer2,      fSettings->BufferSize());     
    fTree.Branch("TISTARLayer3",        &fTISTARLayer3,      fSettings->BufferSize());     
    
    fTISTARParticleVector = new std::vector<Particle>;
    fTree.Branch("TISTARParticleVector", &fTISTARParticleVector, fSettings->BufferSize());

}

Converter::~Converter() {
    if(fOutput->IsOpen()) {
        if(fSettings->WriteTree())
            fTree.Write("tree");
        for(auto list = fHistograms.begin(); list != fHistograms.end(); ++list) {
            fOutput->mkdir(list->first.c_str());
            fOutput->cd(list->first.c_str());
            list->second->Write();
        }
        fOutput->Close();
    }
}

bool Converter::Run() {
    int status;
    int eventNumber = 0;
    int trackID = 0;
    // for 3d gamma-gamma correlations
    int cry1 = 0;
    int cry2 = 0;
    int det1 = 0;
    int det2 = 0;
    int index = 0;
    double cry1energy = 0;
    double cry2energy = 0;
    double det1energy = 0;
    double det2energy = 0;
    double angle = 0;
    double norm = 0;

    int descantArrayHits;
    int descantBlueHits;
    int descantGreenHits;
    int descantRedHits;
    int descantWhiteHits;
    int descantYellowHits;

    int testcanHits;   

//    double buffer1 = 0;
//    double buffer2 = 0;

    double smearedEnergy;
    std::map<int,int> belowThreshold;
    std::map<int,int> outsideTimeWindow;

    TH1F* hist1D = NULL;
    TH2F* hist2D = NULL;
    //TH3I* hist3D;
    THnSparseF* histND = NULL;
    long int nEntries = fChain.GetEntries();

    //  const char* charbuffer;
    //  std::string stringbuffer;
    //  std::stringstream ss;

// ======================================================================================
    
    bool dontSmear = true;
    bool doubleSidedFirstLayer = false;    

    TistarSettings * sett = fSettings->GetTistarSettings();
    
    std::cout<<"beam N = "<<sett->GetProjectileA()-sett->GetProjectileZ()<<", Z = "<<sett->GetProjectileZ()
        << " on target N = "<<sett->GetTargetA()-sett->GetTargetZ()<<", Z = "<<sett->GetTargetZ()<<std::endl;

    double beamEnergy = sett->GetBeamEnergy(); //initial beam energy (total) in MeV

    bool isSolid = (sett->GetGasTargetLength() == 0.);
    if (fSettings->VerbosityLevel()) {
        if(isSolid) std::cout<<"Using a solid target!"<<std::endl;
        else        std::cout<<"Using a gas taget!"<<std::endl;
    }

    TTree* rectr = new TTree("rectr","reconstructed events");
    std::vector<Particle>* ParticleBranch = new std::vector<Particle>;
    rectr->Branch("Particle",&ParticleBranch);

    HitSim * hit = new HitSim(fSettings);

    std::string massfile = sett->GetMassFile();
    if(fSettings->VerbosityLevel()) std::cout<<"Massfile = "<<massfile<<std::endl;

    if(fSettings->VerbosityLevel()) std::cout<<"creating compound \""<<sett->GetTargetMaterialName().c_str()<<"\""<<std::endl;
    Compound* targetMat = new Compound(sett->GetTargetMaterialName().c_str());
    std::cout << "Target Material from Settings File: " << sett->GetTargetMaterialName().c_str() <<std::endl;
    if(!isSolid) targetMat->SetDensity(targetMat->GetDensity()*sett->GetTargetPressure()/6.24151e+08); // 4He gas density in Geant Material Lib.
    if(fSettings->VerbosityLevel()) std::cout<<"creating compound \"MY\""<<std::endl;
    Compound* foilMat = new Compound("MY");
    //if(fSettings->VerbosityLevel()) std::cout<<"creating compound \""<<sett->GetVacuumChamberGas().c_str()<<"\""<<std::endl;
    //Compound* chamberGasMat = new Compound(sett->GetVacuumChamberGas().c_str());
    if(fSettings->VerbosityLevel()) std::cout<<"creating compound \""<<"helium"<<"\""<<std::endl;
    Compound* chamberGasMat = new Compound("helium");
    chamberGasMat->SetDensity(chamberGasMat->GetDensity()*sett->GetVacuumChamberGasPressure()/6.24151e+08); // internal geant4 units are such that 1 bar = 6.24151e+08 whatevers *** changed by Leila
    std::cout<<"chamberGasMat density: "<<chamberGasMat->GetDensity()<<std::endl;
    Compound* layerMat = new Compound("silicon");

    // 1n transfer
    Nucleus* projectile = new Nucleus(sett->GetProjectileZ(), sett->GetProjectileA() - sett->GetProjectileZ(),   massfile.c_str());
    Nucleus* target     = new Nucleus(sett->GetTargetZ(),     sett->GetTargetA()-sett->GetTargetZ(),             massfile.c_str());
    Nucleus* ejectile   = new Nucleus(sett->GetProjectileZ(), sett->GetProjectileA()-sett->GetProjectileZ() + 1, massfile.c_str());
    Nucleus* recoil     = new Nucleus(sett->GetTargetZ(),     sett->GetTargetA()-sett->GetTargetZ() - 1,         massfile.c_str());

    std::cout<<"projectile "<<projectile->GetSymbol()<<" ("<<projectile->GetA()<<", "<<projectile->GetZ()<<"; "<<projectile->GetMass()<<")"<<std::endl;
    std::cout<<"target "<<target->GetSymbol()<<" ("<<target->GetA()<<", "<<target->GetZ()<<"; "<<target->GetMass()<<")"<<std::endl;
    std::cout<<"ejectile "<<ejectile->GetSymbol()<<" ("<<ejectile->GetA()<<", "<<ejectile->GetZ()<<"; "<<ejectile->GetMass()<<")"<<std::endl;
    std::cout<<"recoil "<<recoil->GetSymbol()<<" ("<<recoil->GetA()<<", "<<recoil->GetZ()<<"; "<<recoil->GetMass()<<")"<<std::endl;

    Kinematics* transferP = new Kinematics(projectile, target, recoil, ejectile, beamEnergy, 0.); //reaction.GetGroundStateTransferP();
    Reconstruction* beamTarget = new Reconstruction(projectile, targetMat);

    double beamEnergyRec;      //beam energy at reaction, reconstructed 
    double recoilEnergyRec;    //recoil energy, reconstructed
    double recoilEnergyRecdE;
    double recoilEnergyRecErest;
    double recoilThetaRec;
    double recoilPhiRec;

    double targetThickness = sett->GetTargetThicknessMgPerCm2();
    double targetLength    = sett->GetTargetPhysicalLength();
    double targetForwardZ  =   targetLength/2.;
    double targetBackwardZ = - targetLength/2.;
    if(fSettings->VerbosityLevel()) { 
        std::cout <<"Target Thickness from Input File: "<< targetThickness<<" mg/cm2"<<std::endl;
        std::cout <<"Target ForwardZ from Input File: "<< targetForwardZ<<" mm"<<std::endl;
        std::cout <<"Target BackwardZ from Input File: "<< targetBackwardZ<<" mm"<<std::endl;
        std::cout <<"Target Length from Input File: "<< targetLength<<" mm"<<std::endl;
    }
    TSpline3* energyInTarget = beamTarget->Thickness2EnergyAfter(beamEnergy, targetThickness, targetThickness/1000., true);
    energyInTarget->Write("energyInTarget");

    // variables for recoil energy loss reconstruction (assuming max. recoil energy of 100 MeV!!!)
    Reconstruction* recoilTarget = new Reconstruction(recoil, targetMat);
    Reconstruction* recoilFoil = new Reconstruction(recoil, foilMat);
    Reconstruction* recoilLayer = new Reconstruction(recoil, layerMat);
    Reconstruction* recoilChamberGas = new Reconstruction(recoil, chamberGasMat);
    TSpline3* recoilTargetRange  = recoilTarget->Energy2Range(100., 0.1, !isSolid);
    TSpline3* recoilTargetEnergy = recoilTarget->Range2Energy(100., 0.1, !isSolid);
    TSpline3* recoilFoilRange  = recoilFoil->Energy2Range(100., 0.1, false);
    TSpline3* recoilFoilEnergy = recoilFoil->Range2Energy(100., 0.1, false);
    TSpline3* recoilLayerRange  = recoilLayer->Energy2Range(100., 0.1, false);
    TSpline3* recoilLayerEnergy = recoilLayer->Range2Energy(100., 0.1, false);
    TSpline3* recoilChamberGasRange  = recoilChamberGas->Energy2Range(100., 0.1, true);
    TSpline3* recoilChamberGasEnergy = recoilChamberGas->Range2Energy(100., 0.1, true); 
    
    double foilDistance = sett->GetTargetDiameter()/2.;

    double firstLayerDistance =     sett->GetLayerPositionVector()[0][0].x(); 
    double secondLayerDistance =    sett->GetLayerPositionVector()[1][0].x();
    double padDistance =            sett->GetLayerPositionVector()[2][0].x();

    double targetWidthMgCm2 = foilDistance*100.*targetMat->GetDensity();// convert from mm to cm, and from there to mg/cm2 using density
    double foilThicknessMgCm2 = sett->GetTargetMylarThickness()*100.*1.39;// convert from um to cm, and from there to mg/cm2 using the solid mylar density 1.39 g/cm3.
    //double firstGasLayerThicknessMgCm2 = 300.*chamberGasMat->GetDensity(); // distance between foil and 1. layer hard-coded to 3 mm *** original
    double firstGasLayerThicknessMgCm2 = (firstLayerDistance - foilDistance)*100.*chamberGasMat->GetDensity(); // distance between foil and 1. layer, foil distance is not hard-coded to 3 mm any more *** changed by Leila, chamber gas density = gas 4he density 0.164e-3 g/cm3
    double secondGasLayerThicknessMgCm2 = (secondLayerDistance - firstLayerDistance)*100.*chamberGasMat->GetDensity(); // distance 1. and 2. layer in cm and converted to mg/cm2 using density
    double thirdGasLayerThicknessMgCm2 = (padDistance - secondLayerDistance)*100.*chamberGasMat->GetDensity(); // distance 2. layer and pad in cm and converted to mg/cm2 using density
    double secondLayerThicknessMgCm2 = sett->GetLayerDimensionVector()[1][0].x()*100.*2.328; // density Silicon = 2.328 g/cm3
    double firstLayerThicknessMgCm2 =  sett->GetLayerDimensionVector()[0][0].x()*100.*2.328; // density Silicon = 2.328 g/cm3

    std::cout<<"\n foil distance "<<foilDistance<<" mm => target width = "<<targetWidthMgCm2<<" mg/cm^2"<<std::endl;
    std::cout<<"first layer distance "<<firstLayerDistance<<" mm and second layer: "<<secondLayerDistance<<" mm"<<std::endl;
    std::cout<<"foil thickness "<<sett->GetTargetMylarThickness()<<" mm => foilThicknessMgCm2 "<<foilThicknessMgCm2<<" mg/cm^2"<<std::endl;
    std::cout<<"1. layer thickness "<<sett->GetLayerDimensionVector()[0][0].x()<<" mm => 1. layerThicknessMgCm2 "<<firstLayerThicknessMgCm2<<" mg/cm^2"<<std::endl;
    std::cout<<"2. layer thickness "<<sett->GetLayerDimensionVector()[1][0].x()<<" mm => 2. layerThicknessMgCm2 "<<secondLayerThicknessMgCm2<<" mg/cm^2"<<std::endl;
    //std::cout<<"distance foil - 1. layer 3. mm => 1. gas layer thickness = "<<firstGasLayerThicknessMgCm2<<" mg/cm^2"<<std::endl; //  original
    std::cout<<"distance foil - 1. layer "<<firstLayerDistance - foilDistance<<" mm => 1. gas layer thickness = "<<firstGasLayerThicknessMgCm2<<" mg/cm^2"<<std::endl; // Leila
    std::cout<<"distance 1. layer - 2. layer "<<secondLayerDistance - firstLayerDistance<<" mm and chambergasdensity "<<chamberGasMat->GetDensity()<<" => 2. gas layer thickness = "<<secondGasLayerThicknessMgCm2<<" mg/cm^2"<<std::endl;
    std::cout<<"distance 2. layer - pad "<<padDistance - secondLayerDistance<<" mm => 3. gas layer thickness = "<<thirdGasLayerThicknessMgCm2<<" mg/cm^2"<<std::endl;
  
    Particle part;

    // create and save energy vs. theta-lab splines for reaction at front/middle/back of target
    //std::cout<<"beam energy at front/middle/back of target: "<<beamEnergy/sett->GetProjectileA()<<"/";
    std::cout<<"beam energy at front/middle/back of target: "<<beamEnergy<<"/";
    //transferP->SetEBeam(beamEnergy/sett->GetProjectileA());
    transferP->SetEBeam(beamEnergy);
    TSpline3* front = transferP->Evslab(0., 180., 1.);
    front->Write("RecoilEVsThetaLabFront");
    //std::cout<<beamTarget->EnergyAfter(beamEnergy, -3, true)/sett->GetProjectileA()<<"/";
    std::cout<<energyInTarget->Eval(targetThickness/2.)/1000.<<"/";
    //transferP->SetEBeam(beamTarget->EnergyAfter(beamEnergy, -3, true)/sett->GetProjectileA());
    transferP->SetEBeam(energyInTarget->Eval(targetThickness/2.)/1000.);
    TSpline3* middle = transferP->Evslab(0., 180., 1.);
    middle->Write("RecoilEVsThetaLabMiddle");
    //std::cout<<beamTarget->EnergyAfter(beamEnergy, -3, true)/sett->GetProjectileA()<<std::endl;
    std::cout<<energyInTarget->Eval(targetThickness)/1000.<<std::endl;
    //transferP->SetEBeam(beamTarget->EnergyAfter(beamEnergy, -3, true)/sett->GetProjectileA());
    transferP->SetEBeam(energyInTarget->Eval(targetThickness)/1000.);

    TSpline3* back = transferP->Evslab(0., 180., 1.);
    back->Write("RecoilEVsThetaLabBack");

    double dE1Eloss=0;
    double dE2Eloss=0;
    double dE1ElossRange=0;
    double dE2ElossRange=0;
    double dE1MeasuredCorr=0;
    double dE1MeasMinRec;
    double dE2MeasMinRec;

    UInt_t nofLevels = fTISTARGenChain.GetMaximum("reaction")+1;

    TRandom3 rndm;
    
    CreateTistarHistograms(transferP);

    for(int i = 0; i < nEntries; ++i) {
        status = fChain.GetEntry(i);
        if(status == -1) {
            std::cerr<<"Error occured, couldn't read entry "<<i<<" from tree "<<fChain.GetName()<<" in file "<<fChain.GetFile()->GetName()<<std::endl;
            continue;
        } else if(status == 0) {
            std::cerr<<"Error occured, entry "<<i<<" in tree "<<fChain.GetName()<<" in file "<<fChain.GetFile()->GetName()<<" doesn't exist"<<std::endl;
            return false;
        }        

        // we write every generated event (not just the ones that create hits), so we must get the entry corresponding to the correct event
        status = fTISTARGenChain.GetEntry(fEventNumber); 
        if(status == -1) {
            std::cerr<<"Error occured, couldn't read entry "<<fEventNumber<<" from tree "<<fTISTARGenChain.GetName()<<" in file "<<
                       fTISTARGenChain.GetFile()->GetName()<<std::endl;
            continue;
        } else if(status == 0) {
            std::cerr<<"Error occured, entry "<<fEventNumber<<" in tree "<<fTISTARGenChain.GetName()<<" in file "<<
                       fTISTARGenChain.GetFile()->GetName()<<" doesn't exist"<<std::endl;
            return false;
        }

        //if this entry is from the next event, we fill the tree with everything we've collected so far (after SupressGriffinion) and reset the vector(s)
        if((fEventNumber != eventNumber) && ((fSettings->SortNumberOfEvents()==0)||(fSettings->SortNumberOfEvents()>=eventNumber))) {
            if(fSettings->VerbosityLevel() > 1) {
                std::cout<<"-> Last hit of the event processed, adding to histograms..."<<std::endl;
            }
            // This method checks that all the "crystal" hits are unique, that is, they have different crystal and detector IDs.
            // If they have the same crystal and detector IDs, then we sum the energies together.
            // Normally the Geant4 simulation would sum energy deposits on the same volume, but if we ran the code in "step mode",
            // or if we merged two ntuples together, this would not be true. This method checks and will do what "hit mode" in Geant4 normally does for us!
            // This could slow things down, if that's the case we can put a flag in the settings such that these checks are not done everytime.
            CheckGriffinCrystalAddback();
            CheckLaBrDetectorAddback();
            CheckEightPiDetectorAddback();
            CheckAncillaryBgoCrystalAddback();
            CheckSceptarDetectorAddback();
            CheckDescantDetectorAddback();
            CheckPacesDetectorAddback();


            for(int j = 0; j < 16; j++) {
                GriffinNeighbours_counted[j] = 0;
            }

            // number of descant hits
            descantBlueHits = fDescantBlueDetector->size();
            descantGreenHits = fDescantGreenDetector->size();
            descantRedHits = fDescantRedDetector->size();
            descantWhiteHits = fDescantWhiteDetector->size();
            descantYellowHits = fDescantYellowDetector->size();
            descantArrayHits = descantBlueHits + descantGreenHits + descantRedHits + descantWhiteHits + descantYellowHits;

            // number of Testcan hits
            testcanHits = fTestcanDetector->size();

            //---------------------------------------------------------------------
            // Unsuppressed GRIFFIN
            //---------------------------------------------------------------------
            AddbackGriffin();
            if(fSettings->WriteGriffinAddbackVector())
                AddbackGriffinNeighbourVector();

            //statistics histograms
            hist1D = Get1DHistogram("GriffinCrystalMultiplicityUnsup","Statistics");
            hist1D->Fill(fGriffinCrystal->size());
            hist1D = Get1DHistogram("GriffinBgoMultiplicityUnsup","Statistics");
            hist1D->Fill(fGriffinBgo->size());
            hist1D = Get1DHistogram("GriffinDetectorMultiplicityUnsup","Statistics");
            hist1D->Fill(fGriffinDetector->size());
            hist1D = Get1DHistogram("GriffinCrystalHitPattern","Statistics");
            for(size_t firstDet = 0; firstDet < fGriffinCrystal->size(); ++firstDet) {
                hist1D->Fill((4*fGriffinCrystal->at(firstDet).DetectorId())+fGriffinCrystal->at(firstDet).CrystalId());
            }
            hist1D = Get1DHistogram("GriffinDetectorHitPattern","Statistics");
            for(size_t firstDet = 0; firstDet < fGriffinDetector->size(); ++firstDet) {
                hist1D->Fill((fGriffinDetector->at(firstDet).DetectorId()));
            }
            hist1D = Get1DHistogram("SceptarDetectorHitPattern","Statistics");
            for(size_t firstDet = 0; firstDet < fSceptarDetector->size(); ++firstDet) {
                hist1D->Fill((fSceptarDetector->at(firstDet).DetectorId()));
            }
            hist1D = Get1DHistogram("DescantArrayMultiplicity","Statistics");
            hist1D->Fill(fDescantBlueDetector->size()+fDescantGreenDetector->size()+fDescantRedDetector->size()+fDescantWhiteDetector->size()+fDescantYellowDetector->size());
            hist1D = Get1DHistogram("DescantBlueMultiplicity","Statistics");
            hist1D->Fill(fDescantBlueDetector->size());
            hist1D = Get1DHistogram("DescantGreenMultiplicity","Statistics");
            hist1D->Fill(fDescantGreenDetector->size());
            hist1D = Get1DHistogram("DescantRedMultiplicity","Statistics");
            hist1D->Fill(fDescantRedDetector->size());
            hist1D = Get1DHistogram("DescantWhiteMultiplicity","Statistics");
            hist1D->Fill(fDescantWhiteDetector->size());
            hist1D = Get1DHistogram("DescantYellowMultiplicity","Statistics");
            hist1D->Fill(fDescantYellowDetector->size());
            hist1D = Get1DHistogram("DescantBlueHitPattern","Statistics");
            for(size_t firstDet = 0; firstDet < fDescantBlueDetector->size(); ++firstDet) {
                hist1D->Fill((fDescantBlueDetector->at(firstDet).DetectorId()));
            }
            hist1D = Get1DHistogram("DescantGreenHitPattern","Statistics");
            for(size_t firstDet = 0; firstDet < fDescantGreenDetector->size(); ++firstDet) {
                hist1D->Fill((fDescantGreenDetector->at(firstDet).DetectorId()));
            }
            hist1D = Get1DHistogram("DescantRedHitPattern","Statistics");
            for(size_t firstDet = 0; firstDet < fDescantRedDetector->size(); ++firstDet) {
                hist1D->Fill((fDescantRedDetector->at(firstDet).DetectorId()));
            }
            hist1D = Get1DHistogram("DescantWhiteHitPattern","Statistics");
            for(size_t firstDet = 0; firstDet < fDescantWhiteDetector->size(); ++firstDet) {
                hist1D->Fill((fDescantWhiteDetector->at(firstDet).DetectorId()));
            }
            hist1D = Get1DHistogram("DescantYellowHitPattern","Statistics");
            for(size_t firstDet = 0; firstDet < fDescantYellowDetector->size(); ++firstDet) {
                hist1D->Fill((fDescantYellowDetector->at(firstDet).DetectorId()));
            }
            hist1D = Get1DHistogram("TestcanMultiplicity","Statistics");
            hist1D->Fill(fTestcanDetector->size());
            for(size_t firstDet = 0; firstDet < fTestcanDetector->size(); ++firstDet) {
                hist1D->Fill((fTestcanDetector->at(firstDet).DetectorId()));
            }
            hist1D = Get1DHistogram("TISTARMultiplicity","Statistics");
            hist1D->Fill(fTISTARArray->size());
            hist1D = Get1DHistogram("TISTARHitPattern","Statistics");
            for(size_t firstDet = 0; firstDet < fTISTARArray->size(); ++firstDet) {
                hist1D->Fill((fTISTARArray->at(firstDet).DetectorId()));
            }            
    

            // GRIFFIN Crystal
            FillHistDetector1DGamma(hist1D, fGriffinCrystal, "griffin_crystal_unsup_edep_cry", "Griffin1D");
            FillHistDetector1DGammaNR(hist1D, fGriffinCrystal, "griffin_crystal_unsup_edep_cry_nr", "0RES_Griffin1D");

            FillHistDetector2DGammaGamma(hist2D, fGriffinCrystal, "griffin_crystal_unsup_edep_cry_matrix", "Griffin2D");
            FillHistDetector2DGammaGammaNR(hist2D, fGriffinCrystal, "griffin_crystal_unsup_edep_cry_matrix_nr", "0RES_Griffin2D");

            // GRIFFIN Detector / Clover
            FillHistDetector1DGamma(hist1D, fGriffinDetector, "griffin_crystal_unsup_edep", "Griffin1D");
            FillHistDetector1DGammaNR(hist1D, fGriffinDetector, "griffin_crystal_unsup_edep_nr", "0RES_Griffin1D");

            if(fSceptarHit) {
                FillHist2DGriffinSceptarHitPattern(hist2D, fGriffinDetector, fSceptarDetector, "griffin_crystal_sceptar_hit_pattern","Griffin2D");

                FillHistDetector1DGamma(hist1D, fGriffinDetector, "griffin_crystal_unsup_sceptar_coin_edep", "Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinDetector, "griffin_crystal_unsup_sceptar_coin_edep_nr", "0RES_Griffin1D");
                FillHistDetector1DGamma(hist1D, fGriffinCrystal, "griffin_crystal_unsup_sceptar_coin_edep_cry", "Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinCrystal, "griffin_crystal_unsup_sceptar_coin_edep_cry_nr", "0RES_Griffin1D");

                if(fSettings->Write2DSGGHist()) {
                    FillHistDetector2DGammaGamma(hist2D, fGriffinDetector, "griffin_crystal_unsup_sceptar_coin_edep_matrix","Griffin2D");
                    FillHistDetector2DGammaGammaNR(hist2D, fGriffinDetector, "griffin_crystal_unsup_sceptar_coin_edep_matrix_nr","0RES_Griffin2D");
                    FillHistDetector2DGammaGamma(hist2D, fGriffinCrystal, "griffin_crystal_unsup_sceptar_coin_edep_cry_matrix","Griffin2D");
                    FillHistDetector2DGammaGammaNR(hist2D, fGriffinCrystal, "griffin_crystal_unsup_sceptar_coin_edep_cry_matrix_nr","0RES_Griffin2D");
                }

            } else {
                FillHistDetector1DGamma(hist1D, fGriffinDetector, "griffin_crystal_unsup_sceptar_anticoin_edep", "Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinDetector, "griffin_crystal_unsup_sceptar_anticoin_edep_nr", "0RES_Griffin1D");
                FillHistDetector1DGamma(hist1D, fGriffinCrystal, "griffin_crystal_unsup_sceptar_anticoin_edep_cry", "Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinCrystal, "griffin_crystal_unsup_sceptar_anticoin_edep_cry_nr", "0RES_Griffin1D");
                if(fSettings->Write2DSGGHist()) {
                    FillHistDetector2DGammaGamma(hist2D, fGriffinDetector, "griffin_crystal_unsup_sceptar_anticoin_edep_matrix","Griffin2D");
                    FillHistDetector2DGammaGammaNR(hist2D, fGriffinDetector, "griffin_crystal_unsup_sceptar_anticoin_edep_matrix_nr","0RES_Griffin2D");
                    FillHistDetector2DGammaGamma(hist2D, fGriffinCrystal, "griffin_crystal_unsup_sceptar_anticoin_edep_cry_matrix","Griffin2D");
                    FillHistDetector2DGammaGammaNR(hist2D, fGriffinCrystal, "griffin_crystal_unsup_sceptar_anticoin_edep_cry_matrix_nr","0RES_Griffin2D");
                }
            }

            FillHist2DGriffinHitPattern(hist2D, fGriffinDetector, "griffin_crystal_hit_pattern","Griffin2D");

            FillHistDetector2DGammaGamma(hist2D, fGriffinDetector, "griffin_crystal_unsup_edep_matrix","Griffin2D");
            FillHistDetector2DGammaGammaNR(hist2D, fGriffinDetector, "griffin_crystal_unsup_edep_matrix_nr","0RES_Griffin2D");


            // 3D gamma-gamma corr - Crystal Method
            for(size_t firstDet = 0; firstDet < fGriffinCrystal->size(); ++firstDet) {
                if(fSettings->WriteNDHist()) {
                    // add-back 0 deg hits
                    if(fGriffinCrystal->size()==1) {
                        Double_t fillval[3] = {fGriffinCrystal->at(0).Energy(), fGriffinCrystal->at(0).Energy(),0.0};
                        histND = GetNDHistogram("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_sparse","GriffinND");
                        histND->Fill(fillval); //1.0/64);
                    }
                    for(size_t secondDet = firstDet+1; secondDet < fGriffinCrystal->size(); ++secondDet) {
                        cry1energy  = fGriffinCrystal->at(firstDet).Energy();
                        cry1        = fGriffinCrystal->at(firstDet).CrystalId();
                        cry2energy  = fGriffinCrystal->at(secondDet).Energy();
                        cry2        = fGriffinCrystal->at(secondDet).CrystalId();
                        angle = GriffinCryMap[(int)((4*fGriffinCrystal->at(firstDet).DetectorId())+fGriffinCrystal->at(firstDet).CrystalId())][(int)((4*fGriffinCrystal->at(secondDet).DetectorId())+fGriffinCrystal->at(secondDet).CrystalId())];
                        for(int i = 0; i < 52; i++) {
                            if(GriffinCryMapCombos[i][0] == angle) {
                                norm = (double)GriffinCryMapCombos[i][1];
                                index = i;
                                break;
                            }
                        }
                        if(cry1energy == 0 || cry2energy == 0 || norm == 0) {
                            std::cout << "error, didn't find something" << std::endl;
                            std::cout << "cry1energy = " << cry1energy << std::endl;
                            std::cout << "cry2energy = " << cry2energy << std::endl;
                            std::cout << "norm = " << norm << std::endl;
                            std::cout << "angle = " << angle << std::endl;
                        }
                        Double_t fillval2[3] = {fGriffinCrystal->at(firstDet).Energy(), fGriffinCrystal->at(secondDet).Energy(),(double)index};
                        Double_t fillval3[3] = {fGriffinCrystal->at(secondDet).Energy(), fGriffinCrystal->at(firstDet).Energy(),(double)index};
                        histND = GetNDHistogram("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_sparse","GriffinND");
                        histND->Fill(fillval2); //1.0/64);
                        histND->Fill(fillval3); //1.0/64);
                        cry1 = 0;
                        cry2 = 0;
                        cry1energy = 0;
                        cry2energy = 0;
                        angle = 0;
                        norm = 0;
                    }
                }
            }


            // 3D gamma-gamma corr - Detector Method
            for(size_t firstDet = 0; firstDet < fGriffinDetector->size(); ++firstDet) {
                if(fSettings->WriteNDHist()) {
                    // add-back 0 deg hits
                    if(fGriffinDetector->size()==1) {
                        Double_t fillvalab[3] = {fGriffinDetector->at(0).Energy(), fGriffinDetector->at(0).Energy(),0.0};
                        histND = GetNDHistogram("griffin_crystal_unsup_gamma_gamma_corr_edep_det_sparse","GriffinND");
                        histND->Fill(fillvalab); //1.0/64);
                    }
                    for(size_t secondDet = firstDet+1; secondDet < fGriffinDetector->size(); ++secondDet) {
                        det1        = fGriffinDetector->at(firstDet).DetectorId();
                        det2        = fGriffinDetector->at(secondDet).DetectorId();
                        det1energy  = fGriffinDetector->at(firstDet).Energy();
                        det2energy  = fGriffinDetector->at(secondDet).Energy();
                        angle = GriffinDetMap[(int)((fGriffinDetector->at(firstDet).DetectorId()))][(int)((fGriffinDetector->at(secondDet).DetectorId()))];
                        for(int i = 0; i < 7; i++) {
                            if(GriffinDetMapCombos[i][0] == angle) {
                                norm = (double)GriffinDetMapCombos[i][1];
                                index = i;
                                break;
                            }
                        }
                        if(det1energy == 0 || det2energy == 0 || norm == 0) {
                            std::cout << "error, didn't find something" << std::endl;
                            std::cout << "det1energy = " << det1energy << std::endl;
                            std::cout << "det2energy = " << det2energy << std::endl;
                            std::cout << "det1 = " << det1 << std::endl;
                            std::cout << "det2 = " << det2 << std::endl;
                            std::cout << "norm = " << norm << std::endl;
                            std::cout << "angle = " << angle << std::endl;
                        }
                        Double_t fillval2ab[3] = {fGriffinDetector->at(firstDet).Energy(), fGriffinDetector->at(secondDet).Energy(),(double)index};
                        Double_t fillval3ab[3] = {fGriffinDetector->at(secondDet).Energy(), fGriffinDetector->at(firstDet).Energy(),(double)index};
                        histND = GetNDHistogram("griffin_crystal_unsup_gamma_gamma_corr_edep_det_sparse","GriffinND");
                        histND->Fill(fillval2ab); //1.0/64);
                        histND->Fill(fillval3ab); //1.0/64);
                        det1 = 0;
                        det2 = 0;
                        det1energy = 0;
                        det2energy = 0;
                        angle = 0;
                        norm = 0;
                    }
                }
            }



            // 3D gamma-gamma corr - Add-back Method
            for(size_t firstDet = 0; firstDet < fGriffinDetector->size(); ++firstDet) {
                if(fSettings->WriteNDHist()) {
                    cry1 = 0;
                    cry2 = 0;
                    cry1energy = 0;
                    cry2energy = 0;
                    angle = 0;
                    norm = 0;
                    // add-back 0 deg hits - if there's only one detector, then all the interactions are added back to a zero-degree summed hit
                    if(fGriffinDetector->size()==1) {
                       Double_t fillvalabn[3] = {fGriffinDetector->at(0).Energy(), fGriffinDetector->at(0).Energy(),0.0};
                       histND = GetNDHistogram("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_addback_sparse","GriffinND");
                       histND->Fill(fillvalabn); //1.0/64);
                    } // done 0 deg hits
                    else { // we have interactions in multiple detectors!
                        // iterate over summed detector energies
                        for(size_t secondDet = firstDet+1; secondDet < fGriffinDetector->size(); ++secondDet) {
                            for(size_t thiscry = 0; thiscry < fGriffinCrystal->size(); ++thiscry) { // iterate over all interactions
                                 // if this interaction occurred in the first detector...
                                if(fGriffinCrystal->at(thiscry).DetectorId() == fGriffinDetector->at(firstDet).DetectorId() ) {
                                    //...then compare with cry1energy...
                                    if(fGriffinCrystal->at(thiscry).Energy() > cry1energy){
                                       //...and if the new energy is larger, set the crystal 1 ID and the energy.
                                        cry1energy  = fGriffinCrystal->at(thiscry).Energy();
                                        cry1        = fGriffinCrystal->at(thiscry).CrystalId();
                                    }
                                }
                                 // if this interaction occurred in the second detector...
                                if(fGriffinCrystal->at(thiscry).DetectorId() == fGriffinDetector->at(secondDet).DetectorId() ) {
                                    //...then compare with cry2energy...
                                    if(fGriffinCrystal->at(thiscry).Energy() > cry2energy){
                                       //...and if the new energy is larger, set the crystal 2 ID and the energy.
                                        cry2energy  = fGriffinCrystal->at(thiscry).Energy();
                                        cry2        = fGriffinCrystal->at(thiscry).CrystalId();
                                    }
                                }
                            }
                            angle = GriffinCryMap[(int)((4*fGriffinDetector->at(firstDet).DetectorId())+cry1)][(int)((4*fGriffinDetector->at(secondDet).DetectorId())+cry2)];
                            for(int i = 0; i < 52; i++) {
                                if(GriffinCryMapCombos[i][0] == angle) {
                                    norm = (double)GriffinCryMapCombos[i][1];
                                    index = i;
                                    break;
                                }
                            }
                            if(cry1energy == 0 || cry2energy == 0 || norm == 0) {
										 std::cout << "error, didn't find something" << std::endl;
										 std::cout << "cry1energy = " << cry1energy << std::endl;
										 std::cout << "cry2energy = " << cry2energy << std::endl;
										 std::cout << "norm = " << norm << std::endl;
										 std::cout << "angle = " << angle << std::endl;
                            }
                            Double_t fillval2abn[3] = {fGriffinDetector->at(firstDet).Energy(), fGriffinDetector->at(secondDet).Energy(),(double)index};
                            Double_t fillval3abn[3] = {fGriffinDetector->at(secondDet).Energy(), fGriffinDetector->at(firstDet).Energy(),(double)index};
                            histND = GetNDHistogram("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_addback_sparse","GriffinND");
                            histND->Fill(fillval2abn); //1.0/64);
                            histND->Fill(fillval3abn); //1.0/64);
                            cry1 = 0;
                            cry2 = 0;
                            cry1energy = 0;
                            cry2energy = 0;
                            angle = 0;
                            norm = 0;
                        }
                    }
                }
            }

            // Neighbours
            FillHistDetector1DGamma(hist1D, fGriffinNeighbour, "griffin_crystal_unsup_edep_neigh", "Griffin1D");
            FillHistDetector1DGammaNR(hist1D, fGriffinNeighbour, "griffin_crystal_unsup_edep_neigh_nr", "0RES_Griffin1D");

            if(fSceptarHit) {
                FillHistDetector1DGamma(hist1D, fGriffinNeighbour, "griffin_crystal_unsup_sceptar_coin_edep_neigh", "Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinNeighbour, "griffin_crystal_unsup_sceptar_coin_edep_neigh_nr", "0RES_Griffin1D");
            } else {
                FillHistDetector1DGamma(hist1D, fGriffinNeighbour, "griffin_crystal_unsup_sceptar_anticoin_edep_neigh", "Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinNeighbour, "griffin_crystal_unsup_sceptar_anticoin_edep_neigh_nr", "0RES_Griffin1D");
            }

            // Neighbours Vectors
            FillHistDetector1DGamma(hist1D, fGriffinNeighbourVector, "griffin_crystal_unsup_edep_neighvec", "Griffin1D");
            FillHistDetector1DGammaNR(hist1D, fGriffinNeighbourVector, "griffin_crystal_unsup_edep_neighvec_nr", "0RES_Griffin1D");

            if(fSceptarHit) {
                FillHistDetector1DGamma(hist1D, fGriffinNeighbourVector, "griffin_crystal_unsup_sceptar_coin_edep_neighvec", "Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinNeighbourVector, "griffin_crystal_unsup_sceptar_coin_edep_neighvec_nr", "0RES_Griffin1D");
            } else {
                FillHistDetector1DGamma(hist1D, fGriffinNeighbourVector, "griffin_crystal_unsup_sceptar_anticoin_edep_neighvec", "Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinNeighbourVector, "griffin_crystal_unsup_sceptar_anticoin_edep_neighvec_nr", "0RES_Griffin1D");
            }

            FillHistDetector1DGamma(hist1D, fGriffinArray, "griffin_crystal_unsup_edep_sum", "Griffin1D");
            FillHistDetector1DGammaNR(hist1D, fGriffinArray, "griffin_crystal_unsup_edep_sum_nr", "0RES_Griffin1D");

            // descant coin hits
            if(descantArrayHits == 0) {
                FillHistDetector1DGamma(hist1D, fGriffinCrystal, "griffin_crystal_unsup_descanthit0_coin_edep_cry", "Griffin1D");
                FillHistDetector1DGamma(hist1D, fGriffinDetector, "griffin_crystal_unsup_descanthit0_coin_edep", "Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinCrystal, "griffin_crystal_unsup_descanthit0_coin_edep_cry_nr", "0RES_Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinDetector, "griffin_crystal_unsup_descanthit0_coin_edep_nr", "0RES_Griffin1D");

                FillHistDetector2DGammaGamma(hist2D, fGriffinCrystal, "griffin_crystal_unsup_descanthit0_edep_cry_matrix", "Griffin2D");
                FillHistDetector2DGammaGamma(hist2D, fGriffinDetector, "griffin_crystal_unsup_descanthit0_edep_matrix", "Griffin2D");
                FillHistDetector2DGammaGammaNR(hist2D, fGriffinCrystal, "griffin_crystal_unsup_descanthit0_edep_cry_matrix_nr", "0RES_Griffin2D");
                FillHistDetector2DGammaGammaNR(hist2D, fGriffinDetector, "griffin_crystal_unsup_descanthit0_edep_matrix_nr", "0RES_Griffin2D");

            }
            else if(descantArrayHits == 1) {
                FillHistDetector1DGamma(hist1D, fGriffinCrystal, "griffin_crystal_unsup_descanthit1_coin_edep_cry", "Griffin1D");
                FillHistDetector1DGamma(hist1D, fGriffinDetector, "griffin_crystal_unsup_descanthit1_coin_edep", "Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinCrystal, "griffin_crystal_unsup_descanthit1_coin_edep_cry_nr", "0RES_Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinDetector, "griffin_crystal_unsup_descanthit1_coin_edep_nr", "0RES_Griffin1D");

                FillHistDetector2DGammaGamma(hist2D, fGriffinCrystal, "griffin_crystal_unsup_descanthit1_edep_cry_matrix", "Griffin2D");
                FillHistDetector2DGammaGamma(hist2D, fGriffinDetector, "griffin_crystal_unsup_descanthit1_edep_matrix", "Griffin2D");
                FillHistDetector2DGammaGammaNR(hist2D, fGriffinCrystal, "griffin_crystal_unsup_descanthit1_edep_cry_matrix_nr", "0RES_Griffin2D");
                FillHistDetector2DGammaGammaNR(hist2D, fGriffinDetector, "griffin_crystal_unsup_descanthit1_edep_matrix_nr", "0RES_Griffin2D");
            }
            else if(descantArrayHits == 2) {
                FillHistDetector1DGamma(hist1D, fGriffinCrystal, "griffin_crystal_unsup_descanthit2_coin_edep_cry", "Griffin1D");
                FillHistDetector1DGamma(hist1D, fGriffinDetector, "griffin_crystal_unsup_descanthit2_coin_edep", "Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinCrystal, "griffin_crystal_unsup_descanthit2_coin_edep_cry_nr", "0RES_Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinDetector, "griffin_crystal_unsup_descanthit2_coin_edep_nr", "0RES_Griffin1D");

                FillHistDetector2DGammaGamma(hist2D, fGriffinCrystal, "griffin_crystal_unsup_descanthit2_edep_cry_matrix", "Griffin2D");
                FillHistDetector2DGammaGamma(hist2D, fGriffinDetector, "griffin_crystal_unsup_descanthit2_edep_matrix", "Griffin2D");
                FillHistDetector2DGammaGammaNR(hist2D, fGriffinCrystal, "griffin_crystal_unsup_descanthit2_edep_cry_matrix_nr", "0RES_Griffin2D");
                FillHistDetector2DGammaGammaNR(hist2D, fGriffinDetector, "griffin_crystal_unsup_descanthit2_edep_matrix_nr", "0RES_Griffin2D");
            }
            else {
                FillHistDetector1DGamma(hist1D, fGriffinCrystal, "griffin_crystal_unsup_descanthitn_coin_edep_cry", "Griffin1D");
                FillHistDetector1DGamma(hist1D, fGriffinDetector, "griffin_crystal_unsup_descanthitn_coin_edep", "Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinCrystal, "griffin_crystal_unsup_descanthitn_coin_edep_cry_nr", "0RES_Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinDetector, "griffin_crystal_unsup_descanthitn_coin_edep_nr", "0RES_Griffin1D");

                FillHistDetector2DGammaGamma(hist2D, fGriffinCrystal, "griffin_crystal_unsup_descanthitn_edep_cry_matrix", "Griffin2D");
                FillHistDetector2DGammaGamma(hist2D, fGriffinDetector, "griffin_crystal_unsup_descanthitn_edep_matrix", "Griffin2D");
                FillHistDetector2DGammaGammaNR(hist2D, fGriffinCrystal, "griffin_crystal_unsup_descanthitn_edep_cry_matrix_nr", "0RES_Griffin2D");
                FillHistDetector2DGammaGammaNR(hist2D, fGriffinDetector, "griffin_crystal_unsup_descanthitn_edep_matrix_nr", "0RES_Griffin2D");
            }


            // CLEAR GRIFFIN //
            fGriffinDetector->clear();
            fGriffinNeighbour->clear();
            fGriffinNeighbourVector->clear();
            fGriffinArray->clear();

            //---------------------------------------------------------------------
            // Suppressed GRIFFIN
            //---------------------------------------------------------------------
            SupressGriffin();
            AddbackGriffin();
            if(fSettings->WriteGriffinAddbackVector())
                AddbackGriffinNeighbourVector();

            if(fSettings->WriteTree())
                fTree.Fill(); // Tree contains suppressed data

            //-------------------- crystal histograms
            //multiplicity histogram
            hist1D = Get1DHistogram("GriffinCrystalMultiplicitySup","Statistics");
            hist1D->Fill(fGriffinCrystal->size());
            hist1D = Get1DHistogram("GriffinBgoMultiplicitySup","Statistics");
            hist1D->Fill(fGriffinBgo->size());
            hist1D = Get1DHistogram("GriffinDetectorMultiplicitySup","Statistics");
            hist1D->Fill(fGriffinDetector->size());

            // GRIFFIN Crystal
            FillHistDetector1DGamma(hist1D, fGriffinCrystal, "griffin_crystal_sup_edep_cry", "Griffin1D");
            FillHistDetector1DGammaNR(hist1D, fGriffinCrystal, "griffin_crystal_sup_edep_cry_nr", "0RES_Griffin1D");

            FillHistDetector2DGammaGamma(hist2D, fGriffinCrystal, "griffin_crystal_sup_edep_cry_matrix", "Griffin2D");
            FillHistDetector2DGammaGammaNR(hist2D, fGriffinCrystal, "griffin_crystal_sup_edep_cry_matrix_nr", "0RES_Griffin2D");

            if(fGriffinBgo->size() == 0 && fGriffinBgoBack->size() == 0) {
                FillHistDetector1DGamma(hist1D, fGriffinCrystal, "griffin_crystal_arraysup_edep_cry", "Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinCrystal, "griffin_crystal_arraysup_edep_cry_nr", "0RES_Griffin1D");

                FillHistDetector2DGammaGamma(hist2D, fGriffinCrystal, "griffin_crystal_arraysup_edep_cry_matrix","Griffin2D");
                FillHistDetector2DGammaGammaNR(hist2D, fGriffinCrystal, "griffin_crystal_arraysup_edep_cry_matrix_nr","0RES_Griffin2D");
            }

            // GRIFFIN Detector / Clover
            FillHistDetector1DGamma(hist1D, fGriffinDetector, "griffin_crystal_sup_edep", "Griffin1D");
            FillHistDetector1DGammaNR(hist1D, fGriffinDetector, "griffin_crystal_sup_edep_nr", "0RES_Griffin1D");

            if(fSceptarHit) {
                FillHistDetector1DGamma(hist1D, fGriffinDetector, "griffin_crystal_sup_sceptar_coin_edep", "Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinDetector, "griffin_crystal_sup_sceptar_coin_edep_nr", "0RES_Griffin1D");
                FillHistDetector1DGamma(hist1D, fGriffinCrystal, "griffin_crystal_sup_sceptar_coin_edep_cry", "Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinCrystal, "griffin_crystal_sup_sceptar_coin_edep_cry_nr", "0RES_Griffin1D");
                if(fSettings->Write2DSGGHist()) {
                    FillHistDetector2DGammaGamma(hist2D, fGriffinDetector, "griffin_crystal_sup_sceptar_coin_edep_matrix","Griffin2D");
                    FillHistDetector2DGammaGammaNR(hist2D, fGriffinDetector, "griffin_crystal_sup_sceptar_coin_edep_matrix_nr","0RES_Griffin2D");
                    FillHistDetector2DGammaGamma(hist2D, fGriffinCrystal, "griffin_crystal_sup_sceptar_coin_edep_cry_matrix","Griffin2D");
                    FillHistDetector2DGammaGammaNR(hist2D, fGriffinCrystal, "griffin_crystal_sup_sceptar_coin_edep_cry_matrix_nr","0RES_Griffin2D");
                }
            } else {
                FillHistDetector1DGamma(hist1D, fGriffinDetector, "griffin_crystal_sup_sceptar_anticoin_edep", "Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinDetector, "griffin_crystal_sup_sceptar_anticoin_edep_nr", "0RES_Griffin1D");
                FillHistDetector1DGamma(hist1D, fGriffinCrystal, "griffin_crystal_sup_sceptar_anticoin_edep_cry", "Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinCrystal, "griffin_crystal_sup_sceptar_anticoin_edep_cry_nr", "0RES_Griffin1D");
                if(fSettings->Write2DSGGHist()) {
                    FillHistDetector2DGammaGamma(hist2D, fGriffinDetector, "griffin_crystal_sup_sceptar_anticoin_edep_matrix","Griffin2D");
                    FillHistDetector2DGammaGammaNR(hist2D, fGriffinDetector, "griffin_crystal_sup_sceptar_anticoin_edep_matrix_nr","0RES_Griffin2D");
                    FillHistDetector2DGammaGamma(hist2D, fGriffinCrystal, "griffin_crystal_sup_sceptar_anticoin_edep_cry_matrix","Griffin2D");
                    FillHistDetector2DGammaGammaNR(hist2D, fGriffinCrystal, "griffin_crystal_sup_sceptar_anticoin_edep_cry_matrix_nr","0RES_Griffin2D");

                }
            }

            FillHistDetector2DGammaGamma(hist2D, fGriffinDetector, "griffin_crystal_sup_edep_matrix","Griffin2D");
            FillHistDetector2DGammaGammaNR(hist2D, fGriffinDetector, "griffin_crystal_sup_edep_matrix_nr","0RES_Griffin2D");


            // Neighbours
            FillHistDetector1DGamma(hist1D, fGriffinNeighbour, "griffin_crystal_sup_edep_neigh", "Griffin1D");
            FillHistDetector1DGammaNR(hist1D, fGriffinNeighbour, "griffin_crystal_sup_edep_neigh_nr", "0RES_Griffin1D");

            if(fSceptarHit) {
                FillHistDetector1DGamma(hist1D, fGriffinNeighbour, "griffin_crystal_sup_sceptar_coin_edep_neigh", "Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinNeighbour, "griffin_crystal_sup_sceptar_coin_edep_neigh_nr", "0RES_Griffin1D");
            } else {
                FillHistDetector1DGamma(hist1D, fGriffinNeighbour, "griffin_crystal_sup_sceptar_anticoin_edep_neigh", "Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinNeighbour, "griffin_crystal_sup_sceptar_anticoin_edep_neigh_nr", "0RES_Griffin1D");
            }

            // Neighbours Vectors
            FillHistDetector1DGamma(hist1D, fGriffinNeighbourVector, "griffin_crystal_sup_edep_neighvec", "Griffin1D");
            FillHistDetector1DGammaNR(hist1D, fGriffinNeighbourVector, "griffin_crystal_sup_edep_neighvec_nr", "0RES_Griffin1D");

            if(fSceptarHit) {
                FillHistDetector1DGamma(hist1D, fGriffinNeighbourVector, "griffin_crystal_sup_sceptar_coin_edep_neighvec", "Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinNeighbourVector, "griffin_crystal_sup_sceptar_coin_edep_neighvec_nr", "0RES_Griffin1D");
            } else {
                FillHistDetector1DGamma(hist1D, fGriffinNeighbourVector, "griffin_crystal_sup_sceptar_anticoin_edep_neighvec", "Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinNeighbourVector, "griffin_crystal_sup_sceptar_anticoin_edep_neighvec_nr", "0RES_Griffin1D");
            }

            // GRIFFIN Detector / Clover
            if(fGriffinBgo->size() == 0 && fGriffinBgoBack->size() == 0) {
                FillHistDetector1DGamma(hist1D, fGriffinDetector, "griffin_crystal_arraysup_edep", "Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinDetector, "griffin_crystal_arraysup_edep_nr", "0RES_Griffin1D");

                FillHistDetector2DGammaGamma(hist2D, fGriffinDetector, "griffin_crystal_arraysup_edep_matrix","Griffin2D");
                FillHistDetector2DGammaGammaNR(hist2D, fGriffinDetector, "griffin_crystal_arraysup_edep_matrix_nr","0RES_Griffin2D");
            }

            FillHistDetector1DGamma(hist1D, fGriffinArray, "griffin_crystal_sup_edep_sum", "Griffin1D");
            FillHistDetector1DGammaNR(hist1D, fGriffinArray, "griffin_crystal_sup_edep_sum_nr", "0RES_Griffin1D");

            if(fGriffinBgo->size() == 0 && fGriffinBgoBack->size() == 0 ) {
                FillHistDetector1DGamma(hist1D, fGriffinArray, "griffin_crystal_arraysup_edep_sum", "Griffin1D");
                FillHistDetector1DGammaNR(hist1D, fGriffinArray, "griffin_crystal_arraysup_edep_sum_nr", "0RES_Griffin1D");
            }


            // CLEAR GRIFFIN //
            fGriffinDetector->clear();
            fGriffinNeighbour->clear();
            fGriffinNeighbourVector->clear();
            fGriffinArray->clear();

            // SUPPRESSED GRIFFIN with Ancillary Detectors too
            // now suppress and AddbackGriffin again
            SupressGriffin();
            SupressGriffinByNeighbouringAncillaryBgos();
            AddbackGriffin();


            FillHistDetector1DGamma(hist1D, fGriffinCrystal, "griffin_crystal_ancillaryneighsup_edep_cry", "Griffin1D");
            FillHistDetector1DGammaNR(hist1D, fGriffinCrystal, "griffin_crystal_ancillaryneighsup_edep_cry_nr", "0RES_Griffin1D");

            FillHistDetector1DGamma(hist1D, fGriffinDetector, "griffin_crystal_ancillaryneighsup_edep", "Griffin1D");
            FillHistDetector1DGammaNR(hist1D, fGriffinDetector, "griffin_crystal_ancillaryneighsup_edep_nr", "0RES_Griffin1D");

            // CLEAR GRIFFIN //
            fGriffinDetector->clear();
            fGriffinNeighbour->clear();
            fGriffinNeighbourVector->clear();
            fGriffinArray->clear();

            // SUPPRESSED GRIFFIN with SCEPTAR too
            // now suppress and AddbackGriffin again
            SupressGriffin();
            SupressGriffinBySceptar();
            AddbackGriffin();

            FillHistDetector1DGamma(hist1D, fGriffinCrystal, "griffin_crystal_sceptarsup_edep_cry", "Griffin1D");
            FillHistDetector1DGammaNR(hist1D, fGriffinCrystal, "griffin_crystal_sceptarsup_edep_cry_nr", "0RES_Griffin1D");

            FillHistDetector1DGamma(hist1D, fGriffinDetector, "griffin_crystal_sceptarsup_edep", "Griffin1D");
            FillHistDetector1DGammaNR(hist1D, fGriffinDetector, "griffin_crystal_sceptarsup_edep_nr", "0RES_Griffin1D");

            // LaBr3
            // Unsuppressed
            FillHistDetector1DGamma(hist1D, fLaBrDetector, "labr_crystal_unsup_edep", "LaBr1D");
            FillHistDetector1DGammaNR(hist1D, fLaBrDetector, "labr_crystal_unsup_edep_nr", "0RES_LaBr1D");

            AddbackLaBr();

            FillHistDetector1DGamma(hist1D, fLaBrArray, "labr_crystal_unsup_edep_sum", "LaBr1D");
            FillHistDetector1DGammaNR(hist1D, fLaBrArray, "labr_crystal_unsup_edep_sum_nr", "0RES_LaBr1D");

            AddbackAncillaryBgo();
            SupressLaBr();

            FillHistDetector1DGamma(hist1D, fLaBrDetector, "labr_crystal_sup_edep", "LaBr1D");
            FillHistDetector1DGammaNR(hist1D, fLaBrDetector, "labr_crystal_sup_edep_nr", "0RES_LaBr1D");

            if(fAncillaryBgoCrystal->size() == 0) {
                FillHistDetector1DGamma(hist1D, fLaBrArray, "labr_crystal_sup_edep_sum", "LaBr1D");
                FillHistDetector1DGammaNR(hist1D, fLaBrArray, "labr_crystal_sup_edep_sum_nr", "0RES_LaBr1D");
            }

            SupressLaBrByNeighbouringGriffinShields();

            FillHistDetector1DGamma(hist1D, fLaBrDetector, "labr_crystal_griffinneighsup_edep", "LaBr1D");
            FillHistDetector1DGammaNR(hist1D, fLaBrDetector, "labr_crystal_griffinneighsup_edep_nr", "0RES_LaBr1D");

            if(fAncillaryBgoCrystal->size() == 0) {
                FillHistDetector1DGamma(hist1D, fLaBrArray, "labr_crystal_griffinneighsup_edep_sum", "LaBr1D");
                FillHistDetector1DGammaNR(hist1D, fLaBrArray, "labr_crystal_griffinneighsup_edep_sum_nr", "0RES_LaBr1D");
            }

            if(fGriffinBgo->size() == 0 ) {
                FillHistDetector1DGamma(hist1D, fLaBrDetector, "labr_crystal_griffinanysup_edep", "LaBr1D");
                FillHistDetector1DGammaNR(hist1D, fLaBrDetector, "labr_crystal_anygriffinsup_edep_nr", "0RES_LaBr1D");
                if(fAncillaryBgoCrystal->size() == 0 ) {
                    FillHistDetector1DGamma(hist1D, fLaBrArray, "labr_crystal_griffinanysup_edep_sum", "LaBr1D");
                    FillHistDetector1DGammaNR(hist1D, fLaBrArray, "labr_crystal_anygriffinsup_edep_sum_nr", "0RES_LaBr1D");
                }
            }


            // EightPi3
            // Unsuppressed
            FillHistDetector1DGamma(hist1D, fEightPiDetector, "EightPi_crystal_unsup_edep", "EightPi1D");
            FillHistDetector1DGammaNR(hist1D, fEightPiDetector, "EightPi_crystal_unsup_edep_nr", "0RES_EightPi1D");

            AddbackEightPi();

            FillHistDetector1DGamma(hist1D, fEightPiArray, "EightPi_crystal_unsup_edep_sum", "EightPi1D");
            FillHistDetector1DGammaNR(hist1D, fEightPiArray, "EightPi_crystal_unsup_edep_sum_nr", "0RES_EightPi1D");

            SupressEightPi();

            FillHistDetector1DGamma(hist1D, fEightPiDetector, "EightPi_crystal_sup_edep", "EightPi1D");
            FillHistDetector1DGammaNR(hist1D, fEightPiDetector, "EightPi_crystal_sup_edep_nr", "0RES_EightPi1D");


            // SCEPTAR
            FillHistDetector1DGamma(hist1D, fSceptarDetector, "sceptar_crystal_unsup_edep", "Sceptar1D");
            FillHistDetector1DGammaNR(hist1D, fSceptarDetector, "sceptar_crystal_unsup_edep_nr", "0RES_Sceptar1D");

            AddbackSceptar();

            FillHistDetector1DGamma(hist1D, fSceptarArray, "sceptar_crystal_unsup_edep_sum", "Sceptar1D");
            FillHistDetector1DGammaNR(hist1D, fSceptarArray, "sceptar_crystal_unsup_edep_sum_nr", "0RES_Sceptar1D");


            // DESCANT
            FillHistDetector1DGamma(hist1D, fDescantBlueDetector, "descant_blue_scin_unsup_edep", "Descant1D");
            FillHistDetector1DGammaNR(hist1D, fDescantBlueDetector, "descant_blue_scin_unsup_edep_nr", "0RES_Descant1D");
            FillHistDetector1DGamma(hist1D, fDescantGreenDetector, "descant_green_scin_unsup_edep", "Descant1D");
            FillHistDetector1DGammaNR(hist1D, fDescantGreenDetector, "descant_green_scin_unsup_edep_nr", "0RES_Descant1D");
            FillHistDetector1DGamma(hist1D, fDescantRedDetector, "descant_red_scin_unsup_edep", "Descant1D");
            FillHistDetector1DGammaNR(hist1D, fDescantRedDetector, "descant_red_scin_unsup_edep_nr", "0RES_Descant1D");
            FillHistDetector1DGamma(hist1D, fDescantWhiteDetector, "descant_white_scin_unsup_edep", "Descant1D");
            FillHistDetector1DGammaNR(hist1D, fDescantWhiteDetector, "descant_white_scin_unsup_edep_nr", "0RES_Descant1D");
            FillHistDetector1DGamma(hist1D, fDescantYellowDetector, "descant_yellow_scin_unsup_edep", "Descant1D");
            FillHistDetector1DGammaNR(hist1D, fDescantYellowDetector, "descant_yellow_scin_unsup_edep_nr", "0RES_Descant1D");

            AddbackDescant();

            FillHistDetector1DGamma(hist1D, fDescantArray, "descant_array_scin_unsup_edep_sum", "Descant1D");
            FillHistDetector1DGammaNR(hist1D, fDescantArray, "descant_array_scin_unsup_edep_sum_nr", "0RES_Descant1D");

            // Testcan
            FillHistDetector1DGamma(hist1D, fTestcanDetector, "testcan_scin_unsup_edep", "Testcan1D");
            FillHistDetector1DGammaNR(hist1D, fTestcanDetector, "testcan_scin_unsup_edep_nr", "0RES_Testcan1D");

            // Paces
            FillHistDetector1DGamma(hist1D, fPacesDetector, "paces_crystal_unsup_edep", "Paces1D");
            FillHistDetector1DGammaNR(hist1D, fPacesDetector, "paces_crystal_unsup_edep_nr", "0RES_Paces1D");

            AddbackPaces();

            FillHistDetector1DGamma(hist1D, fPacesArray, "paces_crystal_unsup_edep_sum", "Paces1D");
            FillHistDetector1DGammaNR(hist1D, fPacesArray, "paces_crystal_unsup_edep_sum_nr", "0RES_Paces1D");

            // TI-STAR
            FillHistDetector1DGamma(hist1D, fTISTARArray, "TISTAR_array_unsup_edep", "TISTAR1D");
            FillHistDetector1DGammaNR(hist1D, fTISTARArray, "TISTAR_array_unsup_edep_nr", "0RES_TISTAR1D");
            FillHistDetector1DGamma(hist1D, fTISTARLayer1, "TISTAR_layer1_unsup_edep", "TISTAR1D");
            FillHistDetector1DGammaNR(hist1D, fTISTARLayer1, "TISTAR_layer1_unsup_edep_nr", "0RES_TISTAR1D");
            FillHistDetector1DGamma(hist1D, fTISTARLayer2, "TISTAR_layer2_unsup_edep", "TISTAR1D");
            FillHistDetector1DGammaNR(hist1D, fTISTARLayer2, "TISTAR_layer2_unsup_edep_nr", "0RES_TISTAR1D");
            FillHistDetector1DGamma(hist1D, fTISTARLayer3, "TISTAR_layer3_unsup_edep", "TISTAR1D");
            FillHistDetector1DGammaNR(hist1D, fTISTARLayer3, "TISTAR_layer3_unsup_edep_nr", "0RES_TISTAR1D");


            fGriffinCrystal->clear();
            fGriffinDetector->clear();
            fGriffinNeighbour->clear();
            fGriffinNeighbourVector->clear();
            fGriffinArray->clear();
            fGriffinBgo->clear();
            fGriffinBgoBack->clear();

            fLaBrDetector->clear();
            fLaBrArray->clear();

            fAncillaryBgoCrystal->clear();
            fAncillaryBgoDetector->clear();
            fAncillaryBgoArray->clear();

            fEightPiDetector->clear();
            fEightPiBgoDetector->clear();
            fEightPiArray->clear();

            fSceptarDetector->clear();
            fSceptarArray->clear();

            fDescantArray->clear();
            fDescantBlueDetector->clear();
            fDescantGreenDetector->clear();
            fDescantRedDetector->clear();
            fDescantWhiteDetector->clear();
            fDescantYellowDetector->clear();

            fTestcanDetector->clear();

            fPacesDetector->clear();
            fPacesArray->clear();

            fTISTARArray->clear();
            fTISTARLayer1->clear();            
            fTISTARLayer2->clear();            
            fTISTARLayer3->clear();            

            belowThreshold.clear();
            outsideTimeWindow.clear();

            fSceptarHit = false;
         
            // process TISTAR hits
            // get the entry corresponding to the event we just finished looping over
            status = fTISTARGenChain.GetEntry(eventNumber); 
            if(status == -1) {
                std::cerr<<"Error occured, couldn't read entry "<<fEventNumber<<" from tree "<<fTISTARGenChain.GetName()<<" in file "<<
                           fTISTARGenChain.GetFile()->GetName()<<std::endl;
                continue;
            } else if(status == 0) {
                std::cerr<<"Error occured, entry "<<fEventNumber<<" in tree "<<fTISTARGenChain.GetName()<<" in file "<<
                           fTISTARGenChain.GetFile()->GetName()<<" doesn't exist"<<std::endl;
                return false;
            }

            if(fSettings->VerbosityLevel() > 1) {
                std::cout<<" # of gamma rays in this event: "<<fTISTARGenGammaEnergy->size()<<std::endl;
                for(size_t i=0; i<fTISTARGenGammaEnergy->size(); ++i) {
                    std::cout<<" gamma #: "<<i<<", energy: "<<fTISTARGenGammaEnergy->at(i)<<", theta: "<<fTISTARGenGammaTheta->at(i)<<", phi: "<<fTISTARGenGammaPhi->at(i)<<std::endl;
                }
            }        

            // push the data collected in the tistar vectors that we fill from looping over the normal "hits" saved 
            // in the ntuple (see FillTistarVectors method) into the ParticleMC vectors we use for the analysis
            FillTistarParticleMCs();

            // fill strip/ring number histograms
            for(int panel=0; panel<4; panel++) {
                for(int hit=0; hit<fTISTARFirstDeltaE[panel]->size(); hit++) {
                    hist1D = Get1DHistogram(Form("Layer1Panel%i_nStripZ",panel+1),"TISTAR1D");        
                    hist1D->Fill( fTISTARFirstDeltaE[panel]->at(hit).GetStripNr().at(0) );
                    hist1D = Get1DHistogram(Form("Layer1Panel%i_nStripY",panel+1),"TISTAR1D");        
                    hist1D->Fill( fTISTARFirstDeltaE[panel]->at(hit).GetRingNr().at(0) );
                    hist2D = Get2DHistogram(Form("Layer1Panel%i_nStripZ_vs_Z",panel+1),"TISTAR2D");        
                    hist2D->Fill( fTISTARFirstDeltaE[panel]->at(hit).GetStripNr().at(0),fTISTARFirstDeltaE[panel]->at(hit).GetPosGlobalZ().at(0) );
                    hist2D = Get2DHistogram(Form("Layer1Panel%i_nStripY_vs_Y",panel+1),"TISTAR2D");        
                    hist2D->Fill( fTISTARFirstDeltaE[panel]->at(hit).GetRingNr().at(0),fTISTARFirstDeltaE[panel]->at(hit).GetPosGlobalY().at(0) );
                }
            }     
            for(int panel=0; panel<2; panel++) {
                for(int hit=0; hit<fTISTARSecondDeltaE[panel]->size(); hit++) {
                    hist1D = Get1DHistogram(Form("Layer2Panel%i_nStripZ",panel+1),"TISTAR1D");        
                    hist1D->Fill( fTISTARSecondDeltaE[panel]->at(hit).GetStripNr().at(0) );
                    hist1D = Get1DHistogram(Form("Layer2Panel%i_nStripY",panel+1),"TISTAR1D");        
                    hist1D->Fill( fTISTARSecondDeltaE[panel]->at(hit).GetRingNr().at(0) );
                    hist2D = Get2DHistogram(Form("Layer2Panel%i_nStripZ_vs_Z",panel+1),"TISTAR2D");        
                    hist2D->Fill( fTISTARSecondDeltaE[panel]->at(hit).GetStripNr().at(0),fTISTARSecondDeltaE[panel]->at(hit).GetPosGlobalZ().at(0) );
                    hist2D = Get2DHistogram(Form("Layer2Panel%i_nStripY_vs_Y",panel+1),"TISTAR2D");        
                    hist2D->Fill( fTISTARSecondDeltaE[panel]->at(hit).GetRingNr().at(0),fTISTARSecondDeltaE[panel]->at(hit).GetPosGlobalY().at(0) );
                }
            }     

            hit->Clear();
            ParticleBranch->clear();
            Int_t silicon_mult_first = fTISTARFirstDeltaE[0]->size() + fTISTARFirstDeltaE[1]->size() + fTISTARFirstDeltaE[2]->size() + fTISTARFirstDeltaE[3]->size();
            Int_t silicon_mult_second = fTISTARSecondDeltaE[0]->size()+ fTISTARSecondDeltaE[1]->size();
            TVector3 firstposition;
            TVector3 secondposition;
        
            Int_t index_first = 0;
            Int_t index_second = 0;

            if(fSettings->VerbosityLevel() > 1 && (silicon_mult_first > 1 || silicon_mult_second > 1)) {
                std::cout<<"Warning: Multiple hits in Silicon Tracker! "<<std::endl
                         <<"First layer:  "<<silicon_mult_first<<" ( ";
                for(auto dir : fTISTARFirstDeltaE) {
                    std::cout<<dir->size()<<" ";
                }
                std::cout<<")  Second layer:  "<<silicon_mult_second<<" ( ";
                for(auto dir : fTISTARSecondDeltaE) {
                    std::cout<<dir->size()<<" ";
                }
                std::cout<<")"<<std::endl;
            }

            if(silicon_mult_first == 1 && (silicon_mult_second == 1 || isSolid)) { // begin mult = 1 events
                if(fTISTARFirstDeltaE[0]->size() == 1 ) {
                    fTISTARFirstDeltaE[0]->at(0).ID(0);
                    if(fSettings->VerbosityLevel() > 1) fTISTARFirstDeltaE[0]->at(0).Print();
                    hit->SetFirstDeltaE(fTISTARFirstDeltaE[0]->at(0), kForward);
                    index_first = 0;
                } else if(fTISTARFirstDeltaE[1]->size() == 1 ) {
                    fTISTARFirstDeltaE[1]->at(0).ID(1);
                    if(fSettings->VerbosityLevel() > 1) fTISTARFirstDeltaE[1]->at(0).Print();
                    hit->SetFirstDeltaE(fTISTARFirstDeltaE[1]->at(0), kForward);
                    index_first = 1;
                } else if(fTISTARFirstDeltaE[2]->size() == 1 ) {
                    fTISTARFirstDeltaE[2]->at(0).ID(2);
                    if(fSettings->VerbosityLevel() > 1) fTISTARFirstDeltaE[2]->at(0).Print();
                    hit->SetFirstDeltaE(fTISTARFirstDeltaE[2]->at(0), kBackward);
                    index_first = 2;
                } else {
                    fTISTARFirstDeltaE[3]->at(0).ID(3);
                    if(fSettings->VerbosityLevel() > 1) fTISTARFirstDeltaE[3]->at(0).Print();
                    hit->SetFirstDeltaE(fTISTARFirstDeltaE[3]->at(0), kBackward);
                    index_first = 3;
                }

                if(fTISTARSecondDeltaE[0]->size() == 1 ) {
                    fTISTARSecondDeltaE[0]->at(0).ID(0);
                    if(fSettings->VerbosityLevel() > 1) fTISTARSecondDeltaE[0]->at(0).Print();
                    hit->SetSecondDeltaE(fTISTARSecondDeltaE[0]->at(0), kForward);
                    index_second = 0;
                } else if(fTISTARSecondDeltaE[1]->size() == 1) {
                    fTISTARSecondDeltaE[1]->at(0).ID(1);
                    if(fSettings->VerbosityLevel() > 1) fTISTARSecondDeltaE[1]->at(0).Print();
                    hit->SetSecondDeltaE(fTISTARSecondDeltaE[1]->at(0), kForward);
                    index_second = 1;
                }

                if(silicon_mult_second == 1) {
                    if(fTISTARPad[index_second]->size() == 1 ) {
                        hit->SetPad(fTISTARPad[index_second]->at(0));
                    }

                    if(fSettings->VerbosityLevel() > 1) {
                        std::cout<<"Using pad "<<index_second<<" with "<<fTISTARPad[index_second]->size()<<" detectors"<<std::endl;
                        for(int p = 0; p < 2; ++p) {
                            for(size_t d = 0; d < fTISTARPad[p]->size(); ++d) {
                                std::cout<<p<<": pad "<<fTISTARPad[p]->at(d).GetID()<<" = "<<fTISTARPad[p]->at(d).GetEdet()<<" keV / "<<fTISTARPad[p]->at(d).GetRear()<<" keV"<<std::endl;
                            }
                        }
                    }
                }

                //get position of hit in first layer
                firstposition = hit->FirstPosition(doubleSidedFirstLayer, !dontSmear);

                // get position of hit in second layer
                if(!isSolid) secondposition = hit->SecondPosition(!dontSmear);
                else         secondposition.SetXYZ(0., 0., 0.);

                part.Clear();
            
                // vector between two hits in Siliocn Tracker
                if(isSolid) part.SetPosition(firstposition);
                else        part.SetPosition(secondposition - firstposition);
                if(fSettings->VerbosityLevel() > 1) {
                    std::cout<<"Position to first hit: "<< firstposition.X()<<" "<<firstposition.Y()<<" "<<firstposition.Z()<<std::endl;
                    std::cout<<"Position to second hit: "<< secondposition.X()<<" "<<secondposition.Y()<<" "<<secondposition.Z()<<std::endl;
                    std::cout<<"Position of relative vector: "<< part.GetPosition().X()<<" "<<part.GetPosition().Y()<<" "<<part.GetPosition().Z()<<std::endl;
                }
        
                // error w/ inf's come up when both first/second position are 
                // calculated as (0,0,0) so we just skip them for now...
                if(firstposition == secondposition)  continue;
            

                // reaction angles
                double recoilThetaSim = fTISTARGenRecoilTheta*180./TMath::Pi();
                recoilThetaRec = part.GetPosition().Theta()*180./TMath::Pi();
                double recoilPhiSim = fTISTARGenRecoilPhi*180./TMath::Pi();
                recoilPhiRec = part.GetPosition().Phi()*180./TMath::Pi();
                if(fSettings->VerbosityLevel() > 1) std::cout<<"reaction phi from position: "<<recoilPhiRec<<" - "<<recoilPhiSim<<" = "<<(recoilPhiRec - recoilPhiSim)<<std::endl;
                if(fSettings->VerbosityLevel() > 1) std::cout<<"reaction theta from position: "<<recoilThetaRec<<" - "<<recoilThetaSim<<" = "<<(recoilThetaRec - recoilThetaSim)<<std::endl;
    
                TVector3 vertex;                   //reconstructed vertex
                if(!isSolid) {
                    //find the closest point between beam axis and vector of the two hits in the silicon tracker
                    TVector3 r = part.GetPosition();  //relative vector from first hit to second hit
                    TVector3 r2 = secondposition;     // vector to second hit
                    double t = 0;                          //line parameter to calculate vertex; temp use only
                    if((r*r - r.Z()*r.Z()) != 0 ) t = (r2*r - (r2.Z()*r.Z()))/(r*r - r.Z()*r.Z());
                    vertex = r2 -( t*r);
                } else {
                    vertex.SetXYZ(0., 0., (targetForwardZ + targetBackwardZ)/2.); // middle of target
                }
                // forcing x/y of vertex to be zero
                //vertex.SetX(0.);
                //vertex.SetY(0.);
                
                if(fSettings->VerbosityLevel() > 1) {
                    std::cout<<"Calculated Vertex:\t"<<vertex.X()<<"\t"<<vertex.Y()<<"\t"<<vertex.Z()<<std::endl;
                    std::cout<<"Simulated Vertex: \t"<<fTISTARGenReactionX<<"\t"<<fTISTARGenReactionY<<"\t"<<fTISTARGenReactionZ<<std::endl;
                }

                //update particle information
                if(vertex.Z() > targetForwardZ) {
                    if(fSettings->VerbosityLevel() > 1) std::cout<<"Correcting vertex z from "<<vertex.Z();
                    vertex.SetZ(targetForwardZ);
                    if(fSettings->VerbosityLevel() > 1) std::cout<<" to "<<vertex.Z()<<std::endl;
                }
                if(vertex.Z() < targetBackwardZ) {
                    if(fSettings->VerbosityLevel() > 1) std::cout<<"Correcting vertex z from "<<vertex.Z();
                    vertex.SetZ(targetBackwardZ);
                    if(fSettings->VerbosityLevel() > 1) std::cout<<" to "<<vertex.Z()<<std::endl;
                }
                hist1D = Get1DHistogram("DeltaZ_VertexCorrection","TISTAR1D",2000,-100,100);        
                hist1D->Fill( (fTISTARGenReactionZ-vertex.Z()) );

                // target length at reaction
                double targetThickEvent;
                if(isSolid) targetThickEvent = targetThickness/2.;
                else        targetThickEvent = targetThickness * ( vertex.Z() - targetBackwardZ ) / targetLength;
                if(fSettings->VerbosityLevel()>1) 
                    std::cout<<"Target Thickness at reaction: "<<targetThickEvent<<" = "<<targetThickness<<" * ( "<<vertex.Z()<<" - "<<targetBackwardZ<<" ) / "<<targetLength<<std::endl;
    
                //calculate target thickness for reconstruction of          
                if(targetThickEvent > 0) beamEnergyRec = energyInTarget->Eval(targetThickEvent)/1000.;
                else                     beamEnergyRec = beamEnergy;
                if(fSettings->VerbosityLevel()>1) std::cout<<"Beam Energy at ???????????????????????????????? Reaction: "<<beamEnergyRec<<" MeV"<<std::endl;

                // reconstruct energy of recoil
                recoilEnergyRecdE    =  hit->GetFirstDeltaEEnergy(fSettings->VerbosityLevel()-1) + hit->GetSecondDeltaEEnergy(fSettings->VerbosityLevel()-1);
                recoilEnergyRecErest =  hit->GetPadEnergy();
                recoilEnergyRec = recoilEnergyRecdE + recoilEnergyRecErest;
                //recoilEnergyRec =  recoilEnergyRecErest;
                if(fSettings->VerbosityLevel()>1) { 
                    std::cout<<" 1. layer energy: "<<hit->GetFirstDeltaEEnergy()<<" 2. layer energy: "<<hit->GetSecondDeltaEEnergy()<<" pad energy: "<<hit->GetPadEnergy()
                             <<" => dE: "<<recoilEnergyRecdE<<", Erest: "<<recoilEnergyRecErest<<" => Erec: "<<recoilEnergyRec<<std::endl;
                }
                // reconstruct energy loss in gas and foil
                double sinTheta = TMath::Sin(recoilThetaRec/180.*TMath::Pi());
                double cosTheta = TMath::Cos(recoilThetaRec/180.*TMath::Pi());
                double tmpPhi = recoilPhiRec; while(tmpPhi > 45.) tmpPhi -= 90.; while(tmpPhi < -45.) tmpPhi += 90.;
                double cosPhi   = TMath::Cos(tmpPhi/180.*TMath::Pi());
                double recoilEnergyRecEloss;
 
                if(isSolid) { // FOR SOLID TARGET
                    //for the solid target we only need to reconstruct the energy loss in the foil and the target (no chamber gas)
                    double range = recoilFoilRange->Eval(recoilEnergyRec);
                    //recoilEnergyRecEloss = recoilFoilEnergy->Eval(range + foilThicknessMgCm2/(sinTheta*cosPhi));// original
                    recoilEnergyRecEloss = recoilFoilEnergy->Eval(range + foilThicknessMgCm2/(sinTheta)); // changed bei Dennis
                    range = recoilTargetRange->Eval(recoilEnergyRecEloss);
                    recoilEnergyRecEloss = recoilTargetEnergy->Eval(range + targetThickness/2./TMath::Abs(cosTheta));
                    // due to changing of the target foil from box to cylinder the ernergy loss is corrected by ommitting cosphi. Leila & Dennis
                } 
                else { // FOR GAS TARGET
                    // need to reconstruct energy loss in chamber gas, foil, and target
                    double range;
                    if(fSettings->VerbosityLevel()>1) std::cout<<"theta "<<recoilThetaRec<<", phi "<<recoilPhiRec<<" (sinTheta "<<sinTheta<<", cosTheta "<<cosTheta<<", cosPhi "<<cosPhi<<"): ";
    
                    //*** energy loss through the pad ***
                    if(recoilEnergyRecErest > 0. ) {
                        if(fSettings->VerbosityLevel()>1) std::cout<<"\n\n *** energy loss through the pad *** "<<std::endl;

                        if(fSettings->VerbosityLevel()>1) {
                            std::cout<<"gas thickness: (padDistance - secondLayerDistance)/(sinTheta * cosPhi)"<<std::endl;
                            std::cout<<"gas thickness: ("<<padDistance<<" - "<<secondLayerDistance<<")/("<<sinTheta<<" * "<<cosPhi<<")"<<std::endl;
                            std::cout<<"from pad "<<recoilEnergyRecErest<<" through "<<(padDistance - secondLayerDistance)/(sinTheta*cosPhi)<<" mm gas \n";
                        }
                        range = recoilChamberGasRange->Eval(recoilEnergyRecErest);

                        if(fSettings->VerbosityLevel()>1) std::cout<<"1. range from the pad Erest "<<range<<std::endl;

                        dE2ElossRange = recoilLayerRange->Eval(recoilChamberGasEnergy->Eval(range + thirdGasLayerThicknessMgCm2/(sinTheta*cosPhi)));
                        if(fSettings->VerbosityLevel()>1) std::cout<<"5.a second layer thickness: "<<sett->GetLayerDimensionVector()[1][0].x()<<" range of the second layer: "<<dE2ElossRange<<std::endl;
                        dE2Eloss = recoilLayerEnergy->Eval(dE2ElossRange + secondLayerThicknessMgCm2/(sinTheta*cosPhi)) - recoilChamberGasEnergy->Eval(range + thirdGasLayerThicknessMgCm2/(sinTheta*cosPhi));
                        dE2MeasMinRec = TMath::Abs(hit->GetSecondDeltaEEnergy(fSettings->VerbosityLevel()>1) - dE2Eloss);
                        
                        Get1DHistogram("hdE2MeasMinRec","TistarAnalysis")->Fill(hit->GetSecondDeltaEEnergy(fSettings->VerbosityLevel()-1) - dE2Eloss);                    
                        Get1DHistogram("hdE2ElossRangeWoEpad0","TistarAnalysis")->Fill(dE2ElossRange);                    
                        Get1DHistogram("hdE2Eloss","TistarAnalysis")->Fill(dE2Eloss);
                        Get1DHistogram("hdE2Measured","TistarAnalysis")->Fill(hit->GetSecondDeltaEEnergy(fSettings->VerbosityLevel()-1));
                        Get2DHistogram("hdE2ElossVsMeasuredWoEpad0","TistarAnalysis")->Fill(dE2Eloss,hit->GetSecondDeltaEEnergy(fSettings->VerbosityLevel()-1));

                        recoilEnergyRecEloss = recoilChamberGasEnergy->Eval(range + thirdGasLayerThicknessMgCm2/(sinTheta*cosPhi)) + hit->GetSecondDeltaEEnergy(fSettings->VerbosityLevel()>1);

                        if(fSettings->VerbosityLevel()>1) {
                            std::cout<<"2. energy loss from the pad + de2 "<<recoilEnergyRecEloss<<std::endl;
                        }
                        range = recoilChamberGasRange->Eval(recoilEnergyRecEloss);
                        if(fSettings->VerbosityLevel()>1) {
                            std::cout<<"3. range from the pad Eloss "<<range<<std::endl;
                        }

                        Get1DHistogram("hErestMeasured","TistarAnalysis")->Fill(hit->GetPadEnergy());
                    }
                    else {
                        range = recoilChamberGasRange->Eval(hit->GetSecondDeltaEEnergy(fSettings->VerbosityLevel()>1));
                        if(fSettings->VerbosityLevel()>1) {
                            std::cout<<"\n range from dE2 for Erest=0 "<<range<<std::endl;
                        }
                        dE2ElossRange = recoilLayerRange->Eval(range);
                        dE2Eloss = recoilLayerEnergy->Eval(dE2ElossRange + secondLayerThicknessMgCm2/(sinTheta*cosPhi)) - recoilChamberGasEnergy->Eval(range);
                        Get2DHistogram("hdE2ElossVsMeasured","TistarAnalysis")->Fill(dE2Eloss,hit->GetSecondDeltaEEnergy(fSettings->VerbosityLevel()>1));
                    }

                    //*** energy loss through the second layer ***  
                    dE1ElossRange = recoilLayerRange->Eval(recoilChamberGasEnergy->Eval(range + secondGasLayerThicknessMgCm2/(sinTheta*cosPhi)));
                    if(fSettings->VerbosityLevel()>1) std::cout<<"7.a first layer thickness: "<<sett->GetLayerDimensionVector()[0][0].x()<<" range of the first layer: "<<dE1ElossRange<<std::endl;
                    dE1Eloss = recoilLayerEnergy->Eval(dE1ElossRange + firstLayerThicknessMgCm2/(sinTheta*cosPhi)) - recoilChamberGasEnergy->Eval(range + secondGasLayerThicknessMgCm2/(sinTheta*cosPhi));
                    dE1MeasMinRec = TMath::Abs(hit->GetFirstDeltaEEnergy(fSettings->VerbosityLevel()>1) - dE1Eloss);
                    Get1DHistogram("hdE1ElossRange","TistarAnalysis")->Fill(dE1ElossRange);
                    Get1DHistogram("hdE1Eloss","TistarAnalysis")->Fill(dE1Eloss);
                    Get1DHistogram("hdE1Measured","TistarAnalysis")->Fill(hit->GetFirstDeltaEEnergy(fSettings->VerbosityLevel()>1));
                    Get1DHistogram("hdE1MeasMinRec","TistarAnalysis")->Fill(hit->GetFirstDeltaEEnergy(fSettings->VerbosityLevel()>1) - dE1Eloss);
                    Get2DHistogram("hdE1ElossVsMeasured","TistarAnalysis")->Fill(dE1Eloss,hit->GetFirstDeltaEEnergy(fSettings->VerbosityLevel()>1));
                    if(recoilEnergyRecErest == 0.) Get2DHistogram("hdE1ElossVsMeasuredEpad0","TistarAnalysis")->Fill(dE1Eloss,hit->GetFirstDeltaEEnergy(fSettings->VerbosityLevel()>1));
                    if(recoilEnergyRecErest > 0.)  Get2DHistogram("hdE1ElossVsMeasuredEpadWo0","TistarAnalysis")->Fill(dE1Eloss,hit->GetFirstDeltaEEnergy(fSettings->VerbosityLevel()>1));

                    if(fSettings->VerbosityLevel()>1) {
                        std::cout<<"7.b first layer thickness: "<<sett->GetLayerDimensionVector()[0][0].x()<<" range of the first layer: "<<dE1ElossRange
                                 <<" energi loss "<<dE1Eloss<<" the measured energy: "<<hit->GetFirstDeltaEEnergy(fSettings->VerbosityLevel()>1)
                                 <<" total: "<<dE1MeasuredCorr<<std::endl;
                        std::cout<<"\n\n *** energy loss through the second layer *** "<<std::endl;
                        std::cout<<" with 2. layer "<<hit->GetSecondDeltaEEnergy(fSettings->VerbosityLevel()>1)<<" ("<<recoilEnergyRecEloss
                                 <<") through "<<(secondLayerDistance - firstLayerDistance)/(sinTheta*cosPhi)<<" mm gas \n";
                        std::cout<<"4. range from the pad or second layer "<<range<<std::endl;
                    }


                    recoilEnergyRecEloss = recoilChamberGasEnergy->Eval(range + secondGasLayerThicknessMgCm2/(sinTheta*cosPhi)) + hit->GetFirstDeltaEEnergy(fSettings->VerbosityLevel()>1);
                    if(fSettings->VerbosityLevel()>1) {
                        std::cout<<"5. energy loss from the second layer + de1 "<<recoilEnergyRecEloss<<std::endl;
                    }

                    // for now the foil is box-shaped as well, so we can just continue the same way

                    //*** energy loss through the first layer ***
                    if(fSettings->VerbosityLevel()>1) {
                        std::cout<<"\n\n *** energy loss through the first layer *** "<<std::endl;
                        std::cout<<" with 1. layer "<<hit->GetFirstDeltaEEnergy(fSettings->VerbosityLevel()>1)<<" ("<<recoilEnergyRecEloss<<") through "<<(firstLayerDistance - foilDistance)/(sinTheta*cosPhi)<<" mm gas \n";
                    }
                    range = recoilChamberGasRange->Eval(recoilEnergyRecEloss);

                    if(fSettings->VerbosityLevel()>1) {
                        std::cout<<"6. range from the pad or second layer "<<range<<std::endl;
                    }
                    //recoilEnergyRecEloss = recoilChamberGasEnergy->Eval(range + firstGasLayerThicknessMgCm2/(sinTheta*cosPhi)); // original
                    recoilEnergyRecEloss = recoilChamberGasEnergy->Eval(range + (targetWidthMgCm2*(1 - cosPhi) + 100.*chamberGasMat->GetDensity()*(firstLayerDistance - foilDistance)/(sinTheta*cosPhi))); 
                    // changed by Leila: effectiveLength(firstLayer - vertex)/sinTheta*cosPhi - targetRadii/sinTheta = a0/sinTheta*cosPhi - Rt/sinTheta 

                    if(fSettings->VerbosityLevel()>1) {
                        std::cout<<"7. energy loss from the first layer "<<recoilEnergyRecEloss<<std::endl;
                    }

                    //*** energy loss through the foil ***
                    if(fSettings->VerbosityLevel()>1) {
                        std::cout<<"\n\n *** energy loss through the foil *** "<<std::endl;
                        std::cout<<" at "<<recoilEnergyRecEloss<<" through "<<foilThicknessMgCm2/(sinTheta*cosPhi)<<" mm foil \n";
                    }
                    range = recoilFoilRange->Eval(recoilEnergyRecEloss);

                    if(fSettings->VerbosityLevel()>1) {
                        std::cout<<"8. range from the foil "<<range<<std::endl;
                    }
                    //recoilEnergyRecEloss = recoilFoilEnergy->Eval(range + foilThicknessMgCm2/(sinTheta*cosPhi));// original 
                    recoilEnergyRecEloss = recoilFoilEnergy->Eval(range + foilThicknessMgCm2/(sinTheta));// changed by Leila

                    if(fSettings->VerbosityLevel()>1) {
                        std::cout<<"9. energy loss from the foil "<<recoilEnergyRecEloss<<std::endl;
                        std::cout<<"\n recoilEnergyRec: "<<recoilEnergyRec<<", range: "<<range<<" + foilThicknessMgCm2/(sinTheta*cosPhi) "
                                 <<foilThicknessMgCm2/(sinTheta)<<" => recoilEnergyRecEloss: "<<recoilEnergyRecEloss<<std::endl;
                    }
            
                    //*** energy loss through the target ***
                    // for now assume that the "box" inside the foil is filled with target gas. Not box any more. It is a cylinder --> phi is ommitted!
                    if(fSettings->VerbosityLevel()>1) {
                        std::cout<<"\n\n *** energy loss through the target *** "<<std::endl;
                        std::cout<<" at "<<recoilEnergyRecEloss<<" through "<<targetWidthMgCm2/(sinTheta*cosPhi)<<" mm gas \n";
                    }
                    range = recoilTargetRange->Eval(recoilEnergyRecEloss);

                    if(fSettings->VerbosityLevel()>1) {
                        std::cout<<"10. range from the target "<<range<<std::endl;
                    }
                    //recoilEnergyRecEloss = recoilTargetEnergy->Eval(range + targetWidthMgCm2/(sinTheta*cosPhi));// original
                    recoilEnergyRecEloss = recoilTargetEnergy->Eval(range + targetWidthMgCm2/(sinTheta));// changed by Leila

                    if(fSettings->VerbosityLevel()>1) {
                        std::cout<<"11. energy loss from the target "<<recoilEnergyRecEloss<<std::endl;
                        std::cout<<" => "<<recoilEnergyRecEloss<<std::endl;
                    }


                } // end gas target
                if(fSettings->VerbosityLevel()>1) {
                    std::cout<<"\n theta "<<recoilThetaRec<<", phi "<<recoilPhiRec<<" (sinTheta "<<sinTheta<<", tmpPhi "<<tmpPhi<<", cosPhi "<<cosPhi<<"): "
                             <<foilThicknessMgCm2/(sinTheta)<<" mg/cm^2, "<<recoilEnergyRec<<" => "<<recoilEnergyRecEloss<<", diff "<<recoilEnergyRecEloss-recoilEnergyRec<<std::endl;
                }

                // position has already been set above
                part.SetRecEnergy(recoilEnergyRec);
                part.SetType(2); //proton; this is for one-neutron transfer, only; this sets the mass of the particle
                part.SetReconstructed(); // set TLorentzVector using mass, rec. energy, and position 

                //////////////////////////
                // Q-value reconstruction
                //////////////////////////

                transferP->SetEBeam(beamEnergyRec);
                transferP->Final(recoilThetaRec/180.*TMath::Pi(), 2, true);
                transferP->SetAngles(recoilThetaRec/180.*TMath::Pi(), 2, true);
                double excEnergy = transferP->GetExcEnergy(part.GetReconstructed(), fSettings->VerbosityLevel()-1);

                double recoilThetaCmRec = transferP->GetThetacm(3)/TMath::Pi()*180.;
                if(fSettings->VerbosityLevel()>1) {
                    std::cout<<"\n beamEnergyRec "<<beamEnergyRec<<" => eex = "<<excEnergy<<" (middle spline at "<<recoilThetaRec<<" = "<<middle->Eval(recoilThetaRec)<<", recoilEnergyRec = "<<recoilEnergyRec<<")"<<std::endl;
                    std::cout<<"recoilThetaCmRec = "<<recoilThetaCmRec<<", "<<transferP->GetThetacm(3)<<", "<<transferP->GetThetacm(2)<<", "<<transferP->GetThetacm(1)<<", "<<transferP->GetThetacm(0)<<std::endl;

                    if(excEnergy>5000) std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1********************************1 excEnergy: "<<excEnergy<<" 1. layer E "<< hit->GetFirstDeltaEEnergy(fSettings->VerbosityLevel()>1)<<" 2. layer E "<< hit->GetSecondDeltaEEnergy(fSettings->VerbosityLevel()>1)<<" Epad: "<<recoilEnergyRecErest<<" Epad from hit: "<<hit->GetPadEnergy()<<std::endl;
                }

                ///////////////////////
                // Fill some histograms
                ///////////////////////

                Get1DHistogram("reaction","TistarAnalysis")->Fill(fTISTARGenReaction);
                Get2DHistogram("hitpattern","TistarAnalysis")->Fill(index_first, index_second);
                Get2DHistogram("originXY","TistarAnalysis")->Fill(vertex.X(), vertex.Y());
                Get2DHistogram("originXYErr","TistarAnalysis")->Fill(vertex.X() - fTISTARGenReactionX, vertex.Y() - fTISTARGenReactionY);
                Get2DHistogram("errorOrigin","TistarAnalysis")->Fill(vertex.Z(),  vertex.Z() - fTISTARGenReactionZ );
                Get2DHistogram("errorThetaPhi","TistarAnalysis")->Fill(recoilThetaRec - recoilThetaSim, recoilPhiRec - recoilPhiSim);
                Get2DHistogram("dE12VsPad","TistarAnalysis")->Fill(recoilEnergyRecErest, recoilEnergyRecdE );
                Get2DHistogram("dE12VsE","TistarAnalysis")->Fill(recoilEnergyRec, recoilEnergyRecdE );
                Get2DHistogram("dE1VsE","TistarAnalysis")->Fill(recoilEnergyRec, hit->GetFirstDeltaEEnergy(fSettings->VerbosityLevel()-1));//(firstDeltaE[index_first]->at(0)).GetRear() );
                Get2DHistogram("dE2VsE","TistarAnalysis")->Fill(recoilEnergyRec, hit->GetSecondDeltaEEnergy(fSettings->VerbosityLevel()-1));//(secondDeltaE[index_second]->at(0)).GetRear() );
                Get2DHistogram("dE1VsdE2","TistarAnalysis")->Fill(hit->GetFirstDeltaEEnergy(fSettings->VerbosityLevel()-1), hit->GetSecondDeltaEEnergy(fSettings->VerbosityLevel()-1));//(firstDeltaE[index_first]->at(0)).GetRear() );
                Get2DHistogram("eVsTheta","TistarAnalysis")->Fill(recoilThetaRec, recoilEnergyRec);
                Get2DHistogram("eVsZ","TistarAnalysis")->Fill(vertex.Z(), recoilEnergyRec);

                Get2DHistogram("eRecErrVsESim","TistarAnalysis")->Fill(fTISTARGenRecoilEnergy, recoilEnergyRec - fTISTARGenRecoilEnergy);
                Get2DHistogram("thetaErrorVsZ","TistarAnalysis")->Fill(vertex.Z(), recoilThetaRec - recoilThetaSim);
                //if(hit->GetPadEnergy()>1.00) thetaErrorVsTheta->Fill(recoilThetaSim , recoilThetaRec - recoilThetaSim);
                Get2DHistogram("thetaErrorVsTheta","TistarAnalysis")->Fill(recoilThetaSim , recoilThetaRec - recoilThetaSim); // why twice?
                Get2DHistogram("thetaErrorVsTheta","TistarAnalysis")->Fill(recoilThetaSim , recoilThetaRec - recoilThetaSim);
                if(recoilEnergyRecErest > 0.)  Get2DHistogram("thetaErrorVsThetaEpadCut","TistarAnalysis")->Fill(recoilThetaSim , recoilThetaRec - recoilThetaSim);
                if(fTISTARGenReactionBeamEnergyCM > 0.0) Get2DHistogram("zReactionEnergy","TistarAnalysis")->Fill(vertex.Z(), beamEnergyRec);
                //if(reactionEnergyBeamCM == -1.0) std::cout<<"leila!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<reactionEnergyBeamCM<<std::endl;
                Get1DHistogram("excEnProton","TistarAnalysis")->Fill(excEnergy);
                //if(recoilEnergyRecErest>1.00) excEnProtonVsTheta->Fill(recoilThetaRec, excEnergy); ???????????????????????????
                Get2DHistogram("excEnProtonVsTheta","TistarAnalysis")->Fill(recoilThetaRec, excEnergy);
                Get2DHistogram("excEnProtonVsPhi","TistarAnalysis")->Fill(recoilPhiRec, excEnergy);
                Get2DHistogram("excEnProtonVsThetaCm","TistarAnalysis")->Fill(recoilThetaCmRec, excEnergy);
                Get2DHistogram("excEnProtonVsZ","TistarAnalysis")->Fill(vertex.Z(), excEnergy);
                if(fTISTARGenReaction == 0) {
                    Get2DHistogram("excEnProtonVsThetaGS","TistarAnalysis")->Fill(recoilThetaRec, excEnergy);
                    Get2DHistogram("excEnProtonVsZGS","TistarAnalysis")->Fill(vertex.Z(), excEnergy);
                }
                Get2DHistogram("thetaVsZ","TistarAnalysis")->Fill(vertex.Z(), recoilThetaRec);
                if(index_first == index_second) {
                    Get2DHistogram("eVsZSame","TistarAnalysis")->Fill(vertex.Z(), recoilEnergyRec);
                    Get2DHistogram("thetaVsZSame","TistarAnalysis")->Fill(vertex.Z(), recoilThetaRec);
                } else {
                    Get2DHistogram("eVsZCross","TistarAnalysis")->Fill(vertex.Z(), recoilEnergyRec);
                    Get2DHistogram("thetaVsZCross","TistarAnalysis")->Fill(vertex.Z(), recoilThetaRec);
                }
                Get2DHistogram("phiVsZ","TistarAnalysis")->Fill(vertex.Z(), recoilPhiRec);
                Get2DHistogram("phiErrorVsPhi","TistarAnalysis")->Fill(recoilPhiRec, recoilPhiRec - recoilPhiSim);
                if(fTISTARFirstDeltaE[index_first]->at(0).GetID() == 0) {

                }
                Get2DHistogram("betaCmVsZ","TistarAnalysis")->Fill(vertex.Z(), transferP->GetBetacm());
                Get2DHistogram("eCmVsZ","TistarAnalysis")->Fill(vertex.Z(), transferP->GetCmEnergy()/1000.);
                if(silicon_mult_second == 1) Get2DHistogram("stripPattern","TistarAnalysis")->Fill(index_second*fSettings->GetTISTARnStripsY(0) + fTISTARSecondDeltaE[index_second]->at(0).GetStripNr()[0], fTISTARSecondDeltaE[index_second]->at(0).GetID()*fSettings->GetTISTARnStripsZ(0) + fTISTARSecondDeltaE[index_second]->at(0).GetRingNr()[0]);
                Get2DHistogram("recBeamEnergyErrVsZ","TistarAnalysis")->Fill(vertex.Z(), beamEnergyRec - fTISTARGenReactionBeamEnergy);
                Get2DHistogram("thetaCmVsThetaLab","TistarAnalysis")->Fill(recoilThetaRec, recoilThetaCmRec);
                Get2DHistogram("zErrorVsthetaError","TistarAnalysis")->Fill(recoilThetaRec - recoilThetaSim, vertex.Z() - fTISTARGenReactionZ);
                Get2DHistogram("elossVsTheta","TistarAnalysis")->Fill(recoilThetaRec, recoilEnergyRecEloss - recoilEnergyRec);
                Get2DHistogram("elossVsPhi","TistarAnalysis")->Fill(recoilPhiRec, recoilEnergyRecEloss - recoilEnergyRec);

                Get2DHistogram("dE2VsdE2Pad","TistarAnalysis")->Fill(hit->GetSecondDeltaEEnergy(fSettings->VerbosityLevel())+hit->GetPadEnergy(),hit->GetSecondDeltaEEnergy(fSettings->VerbosityLevel()));
                Get2DHistogram("EPadVsThetaLab","TistarAnalysis")->Fill(recoilThetaRec,hit->GetPadEnergy());
                Get2DHistogram("EPadVsZ","TistarAnalysis")->Fill(vertex.Z(),hit->GetPadEnergy());
                if(recoilEnergyRecErest == 0.) {Get2DHistogram("dE2VsThetaLabEpadCut","TistarAnalysis")->Fill(recoilThetaRec,hit->GetSecondDeltaEEnergy(fSettings->VerbosityLevel()-1));}
                if(recoilThetaRec>0.0 && recoilThetaRec<180.0) {Get2DHistogram("dE2VsEPadThetaCut","TistarAnalysis")->Fill(hit->GetPadEnergy(),hit->GetSecondDeltaEEnergy(fSettings->VerbosityLevel()-1));}
                Get2DHistogram("dE1VsThetaLab","TistarAnalysis")->Fill(recoilThetaRec,hit->GetFirstDeltaEEnergy(fSettings->VerbosityLevel()-1));
                Get2DHistogram("dE2VsThetaLab","TistarAnalysis")->Fill(recoilThetaRec,hit->GetSecondDeltaEEnergy(fSettings->VerbosityLevel()-1));
                Get2DHistogram("dE12VsThetaLab","TistarAnalysis")->Fill(recoilThetaRec,hit->GetFirstDeltaEEnergy(fSettings->VerbosityLevel()-1)+hit->GetSecondDeltaEEnergy(fSettings->VerbosityLevel()-1));
                Get2DHistogram("dE1EpadVsThetaLab","TistarAnalysis")->Fill(recoilThetaRec,hit->GetFirstDeltaEEnergy(fSettings->VerbosityLevel()-1)+hit->GetPadEnergy());
                Get2DHistogram("dE2EpadVsThetaLab","TistarAnalysis")->Fill(recoilThetaRec,hit->GetSecondDeltaEEnergy(fSettings->VerbosityLevel()-1)+hit->GetPadEnergy());

                // ************************* Q-value using the reconstructed energy loss **************************8

                // now we reconstruct the q-value using the reconstructed energy loss
                // position has already been set above
                part.SetRecEnergy(recoilEnergyRecEloss);
                part.SetType(2); //proton; this is for one-neutron transfer, only; this sets the mass of the particle
                part.SetReconstructed(); // set TLorentzVector using mass, rec. energy, and position 
                transferP->Final(recoilThetaRec/180.*TMath::Pi(), 2, true);
                transferP->SetAngles(recoilThetaRec/180.*TMath::Pi(), 2, true);
                excEnergy = transferP->GetExcEnergy(part.GetReconstructed(), fSettings->VerbosityLevel()-1);

                recoilThetaCmRec = transferP->GetThetacm(3)/TMath::Pi()*180.;


        
                Get2DHistogram("hEbeamRecVsSim","TistarAnalysis")->Fill(beamEnergyRec,beamEnergy);

                Get2DHistogram("excEnElossVsTheta","TistarAnalysis")->Fill(recoilThetaRec, excEnergy);
                if(recoilEnergyRecErest > 0. ) {Get2DHistogram("excEnElossVsThetaEpadCut","TistarAnalysis")->Fill(recoilThetaRec, excEnergy);}
                Get1DHistogram("excEnProtonCorr","TistarAnalysis")->Fill(excEnergy);
                if(recoilEnergyRecErest > 0. ) {Get1DHistogram("excEnProtonCorrEpadCut","TistarAnalysis")->Fill(excEnergy);}
                Get2DHistogram("excEnProtonCorrVsX","TistarAnalysis")->Fill(vertex.X(),excEnergy);
                Get2DHistogram("excEnProtonCorrVsY","TistarAnalysis")->Fill(vertex.Y(),excEnergy);
                Get2DHistogram("excEnProtonCorrVsZ","TistarAnalysis")->Fill(vertex.Z(),excEnergy);
                Get2DHistogram("excEnProtonCorrVsT","TistarAnalysis")->Fill(sqrt(vertex.X()*vertex.X()+vertex.Y()*vertex.Y()),excEnergy);
                Get2DHistogram("excEnProtonCorrVsR","TistarAnalysis")->Fill(sqrt(vertex.X()*vertex.X()+vertex.Y()*vertex.Y()+vertex.Z()*vertex.Z()),excEnergy);
                if(dE1MeasMinRec<20.0 && dE2MeasMinRec<40.0) Get1DHistogram("excEnProtonCorrdE1Sigma1","TistarAnalysis")->Fill(excEnergy);
                if(dE1MeasMinRec<40.0 && dE2MeasMinRec<80.0) Get1DHistogram("excEnProtonCorrdE1Sigma2","TistarAnalysis")->Fill(excEnergy);

                //if(0 <= fTISTARGenReaction && fTISTARGenReaction < nofLevels) Get2DHistogram(Form("excEnElossVsThetaLevel_%d",fTISTARGenReaction),"TistarAnalysis")->Fill(recoilThetaRec, excEnergy);
                //if(0 <= reactionSim && reactionSim < nofLevels-1) excEnProtonVsTheta->Fill(recoilThetaRec, excEnergy);//leila 

                // gamma stuff 
                size_t gammaSize = fTISTARGenGammaEnergy->size();
                // doppler correction
                double gamma = (ejectile->GetMass()+fTISTARGenEjectileEnergy/1000.)/ejectile->GetMass();
                double beta = TMath::Sqrt(1.-TMath::Power(1./gamma, 2.));
                double eGammaDoppCorr, resolvedEnergy;
                for(size_t i=0; i<gammaSize; i++) {
                    Get1DHistogram("gammaSpec","TistarAnalysis")->Fill(fTISTARGenGammaEnergy->at(i));
                    Get2DHistogram("excEnProtonVsGamma","TistarAnalysis")->Fill(fTISTARGenGammaEnergy->at(i),excEnergy);
                    
                    eGammaDoppCorr = (1.-beta*TMath::Cos(fTISTARGenGammaTheta->at(i)))/TMath::Sqrt(1.-beta*beta)*fTISTARGenGammaEnergy->at(i); 
                    Get1DHistogram("gammaSpecDoppCorr","TistarAnalysis")->Fill(eGammaDoppCorr);
                    Get2DHistogram("excEnProtonVsGammaDoppCorr","TistarAnalysis")->Fill(eGammaDoppCorr, excEnergy);

                    resolvedEnergy = rndm.Gaus(eGammaDoppCorr,eGammaDoppCorr*0.01/(2.*TMath::Sqrt(2.*TMath::Log(2.))));
                    Get1DHistogram("gammaSpecDoppCorrRes","TistarAnalysis")->Fill(resolvedEnergy);
                    Get2DHistogram("excEnProtonVsGammaDoppCorrRes","TistarAnalysis")->Fill(resolvedEnergy, excEnergy);
                }

            } // end mult = 1 events

            ClearTistarVectors();
            ClearTistarParticleMCs();
            
            // re-set this at the end, as we need the previous hit event number to get
            // the correct entry in the generator tree
            eventNumber = fEventNumber;
        }

        if(fSettings->VerbosityLevel() > 1) {
            std::cout<<"Entry: "<<i<<", Event: "<<fEventNumber<<", Track: "<<fTrackID<<", Det: "<<fDetNumber<<", Cry: "
                     <<fCryNumber<<", Edep: "<<fDepEnergy<<"keV, ParticleID: "<<fParticleType<<", (x,y,z) = ("
                     <<fPosx<<", "<<fPosy<<", "<<fPosz<<" )"<<std::endl;
        }

        // if fSystemID is NOT GRIFFIN, then set fCryNumber to zero
        // This is a quick fix to solve resolution and threshold values from Settings.cc
        // Modified to include TI-STAR, which does use the systemID
        if(fSystemID >= 2000 && fSystemID != 9500) {
            fCryNumber = 0;
        }

        // if fSystemID is TISTAR, decrement the detector and crystal number 
        // this is to solve the numbering problem of the layers/strips, which both
        // start at 1, vs. the vector index which starts at 0
        if(fSystemID == 9500) {
            fDetNumber -= 1;
            fCryNumber -= 1;
        }

        // adding up ti-star hits (checking to see if we have more than 1 in a given strip/ring)
        FillTistarVectors();

        //create energy-resolution smeared energy
        smearedEnergy = fRandom.Gaus(fDepEnergy,fSettings->Resolution(fSystemID,fDetNumber,fCryNumber,fDepEnergy));
        //std::cout << "fDepEnergy= "<<fDepEnergy<<
        //             ", smearedEnergy= "<<smearedEnergy<<
        //             ", fSystemID= "<<fSystemID<<
        //             ", fDetNumber= "<<fDetNumber<<
        //             ", fCryNumber= "<<fCryNumber<<std::endl;

        if((fSettings->SortNumberOfEvents()==0)||(fSettings->SortNumberOfEvents()>=fEventNumber) ) {
            //if the hit is above the threshold, we add it to the vector
            if(AboveThreshold(smearedEnergy, fSystemID)) {
                if(InsideTimeWindow() ) {
                    switch(fSystemID) {
                    case 1000:
                        fGriffinCrystal->push_back(Detector(fEventNumber, fDetNumber, fCryNumber, fDepEnergy, smearedEnergy, TVector3(fPosx,fPosy,fPosz), fTime));
                        break;
                    case 1010:
                    case 1020:
                    case 1030:
                    case 1040:
                        fGriffinBgo->push_back(Detector(fEventNumber, fDetNumber, fCryNumber, fDepEnergy, smearedEnergy, TVector3(fPosx,fPosy,fPosz), fTime));
                        break;
                    case 1050:
                        fGriffinBgoBack->push_back(Detector(fEventNumber, fDetNumber, fCryNumber, fDepEnergy, smearedEnergy, TVector3(fPosx,fPosy,fPosz), fTime));
                        break;
                    case 2000:
                        fLaBrDetector->push_back(Detector(fEventNumber, fDetNumber, fCryNumber, fDepEnergy, smearedEnergy, TVector3(fPosx,fPosy,fPosz), fTime));
                        break;
                    case 3000:
                        fAncillaryBgoCrystal->push_back(Detector(fEventNumber, fDetNumber, fCryNumber, fDepEnergy, smearedEnergy, TVector3(fPosx,fPosy,fPosz), fTime));
                        break;
                    case 5000:
                        fSceptarDetector->push_back(Detector(fEventNumber, fDetNumber, fCryNumber, fDepEnergy, smearedEnergy, TVector3(fPosx,fPosy,fPosz), fTime));
                        fSceptarHit = true;
                        break;
                    case 6000:
                        fEightPiDetector->push_back(Detector(fEventNumber, fDetNumber, fCryNumber, fDepEnergy, smearedEnergy, TVector3(fPosx,fPosy,fPosz), fTime));
                        break;
                    case 6010:
                    case 6020:
                    case 6030:
                        fEightPiBgoDetector->push_back(Detector(fEventNumber, fDetNumber, fCryNumber, fDepEnergy, smearedEnergy, TVector3(fPosx,fPosy,fPosz), fTime));
                        break;
                    case 8010:
                        fDescantBlueDetector->push_back(Detector(fEventNumber, fDetNumber, fCryNumber, fDepEnergy, smearedEnergy, TVector3(fPosx,fPosy,fPosz), fTime));
                         break;
                    case 8020:
                        fDescantGreenDetector->push_back(Detector(fEventNumber, fDetNumber, fCryNumber, fDepEnergy, smearedEnergy, TVector3(fPosx,fPosy,fPosz), fTime));
                        break;
                    case 8030:
                        fDescantRedDetector->push_back(Detector(fEventNumber, fDetNumber, fCryNumber, fDepEnergy, smearedEnergy, TVector3(fPosx,fPosy,fPosz), fTime));
                        break;
                    case 8040:
                        //fDescantWhiteDetector->push_back(Detector(fEventNumber, fDetNumber, fCryNumber, fDepEnergy, smearedEnergy, TVector3(fPosx,fPosy,fPosz), fTime));
                        //fDescantWhiteDetector->push_back(Detector(fEventNumber, fDetNumber, fCryNumber, fEDepD, deuteronSmearedEnergy, TVector3(fPosx,fPosy,fPosz), fTime));
                        //fDescantWhiteDetector->push_back(Detector(fEventNumber, fDetNumber, fCryNumber, fEDepC+fEDepD , carbonSmearedEnergy+deuteronSmearedEnergy, TVector3(fPosx,fPosy,fPosz), fTime));
                        fDescantWhiteDetector->push_back(Detector(fEventNumber, fDetNumber, fCryNumber, fDepEnergy , smearedEnergy, TVector3(fPosx,fPosy,fPosz), fTime));
                        break;
                    case 8050:
                        fDescantYellowDetector->push_back(Detector(fEventNumber, fDetNumber, fCryNumber, fDepEnergy, smearedEnergy, TVector3(fPosx,fPosy,fPosz), fTime));
                        break;

                    case 8500:
                        fTestcanDetector->push_back(Detector(fEventNumber, fDetNumber, fCryNumber, fDepEnergy, smearedEnergy, TVector3(fPosx,fPosy,fPosz), fTime));
                    break;

                    case 9000:
                        fPacesDetector->push_back(Detector(fEventNumber, fDetNumber, fCryNumber, fDepEnergy, smearedEnergy, TVector3(fPosx,fPosy,fPosz), fTime));
                        break;
                    
                    case 9500:
                        fTISTARArray->push_back(Detector(fEventNumber, fDetNumber, fCryNumber, fDepEnergy, smearedEnergy, TVector3(fPosx,fPosy,fPosz), fTime));
                        switch(fDetNumber) {
                            case 0: fTISTARLayer1->push_back(Detector(fEventNumber, fDetNumber, fCryNumber, fDepEnergy, smearedEnergy, TVector3(fPosx,fPosy,fPosz), fTime));
                                    break;
                            case 1: fTISTARLayer2->push_back(Detector(fEventNumber, fDetNumber, fCryNumber, fDepEnergy, smearedEnergy, TVector3(fPosx,fPosy,fPosz), fTime));
                                    break;
                            case 2: fTISTARLayer3->push_back(Detector(fEventNumber, fDetNumber, fCryNumber, fDepEnergy, smearedEnergy, TVector3(fPosx,fPosy,fPosz), fTime));
                                    break;
                        }
                        break;

                    default:
                        std::cerr<<"Unknown detector system ID "<<fSystemID<<std::endl;
                        break;
                    }
                } else {
                    ++outsideTimeWindow[fSystemID];
                }
            } else {
                ++belowThreshold[fSystemID];
            }
        }

        if(i%1000 == 0 && fSettings->VerbosityLevel() > 0) {
            std::cout<<std::setw(3)<<100*i/nEntries<<"% done\r"<<std::flush;
        }
    }

    if(fSettings->VerbosityLevel() > 0) {
        std::cout<<"100% done"<<std::endl;

        if(fSettings->VerbosityLevel() > 1) {
            PrintStatistics();
        }
    }

    return true;
}

bool Converter::AboveThreshold(double energy, int systemID) {
    if(systemID == 5000) {
        // apply hard threshold of 50 keV on Sceptar
        // SCEPTAR in reality saturates at an efficiency of about 80%. In simulation we get an efficiency of 90%
        // 0.9 * 1.11111111 = 100%, 0.8*1.1111111 = 0.888888888
        if(energy > 50.0 && (fRandom.Uniform(0.,1.) < 0.88888888 )) {
            return true;
        }
        else {
            return false;
        }
    }
    else if(energy > fSettings->Threshold(fSystemID,fDetNumber,fCryNumber)+10*fSettings->ThresholdWidth(fSystemID,fDetNumber,fCryNumber)) {
        return true;
    }

    if(fRandom.Uniform(0.,1.) < 0.5*(TMath::Erf((energy-fSettings->Threshold(fSystemID,fDetNumber,fCryNumber))/fSettings->ThresholdWidth(fSystemID,fDetNumber,fCryNumber))+1)) {
        return true;
    }

    return false;
}

bool Converter::InsideTimeWindow() {
    if(fSettings->TimeWindow(fSystemID,fDetNumber,fCryNumber) == 0) {
        return true;
    }
    if(fTime < fSettings->TimeWindow(fSystemID,fDetNumber,fCryNumber)) {
        return true;
    }
    return false;
}

bool Converter::DescantNeutronDiscrimination() { // Assuming perfect gamma-neutron discrimination
    if(fParticleType == 5) { // neutron
        return true;
    }
    return false;
}



void Converter::CheckGriffinCrystalAddback() {
// This method checks that all the "crystal" hits are unique, that is, they have different crystal and detector IDs.
// If they have the same crystal and detector IDs, then we sum the energies together.
// Normally the Geant4 simulation would sum energy deposits on the same volume, but if we ran the code in "step mode",
// or if we merged two ntuples together, this would not be true. This method checks and will do what "hit mode" in Geant4 normally does for us!
    for(auto firstDet = fGriffinCrystal->begin(); firstDet != fGriffinCrystal->end(); ++firstDet) {
        for(auto secondDet = firstDet+1; secondDet != fGriffinCrystal->end();) {
            if((firstDet->DetectorId() == secondDet->DetectorId()) && (firstDet->CrystalId() == secondDet->CrystalId())) {
                firstDet->AddEnergy(secondDet->SimulationEnergy(),secondDet->Energy());
                secondDet = fGriffinCrystal->erase(secondDet);
            }
            else {
                ++secondDet;
            }
        }
    }
}

void Converter::SupressGriffin() {
    //loop over all bgo's and remove all matching germaniums
    for(auto bgo = fGriffinBgo->begin(); bgo != fGriffinBgo->end(); ++bgo) {
        for(auto ge = fGriffinCrystal->begin(); ge != fGriffinCrystal->end();) {
            if(bgo->DetectorId() == ge->DetectorId()) {
                ge = fGriffinCrystal->erase(ge);
            } else {
                ++ge;
            }
        }
    }
    // Now Back Suppressors
    for(auto bgo = fGriffinBgoBack->begin(); bgo != fGriffinBgoBack->end(); ++bgo) {
        for(auto ge = fGriffinCrystal->begin(); ge != fGriffinCrystal->end();) {
            if(bgo->DetectorId() == ge->DetectorId()) {
                ge = fGriffinCrystal->erase(ge);
            } else {
                ++ge;
            }
        }
    }
}

void Converter::SupressGriffinByNeighbouringAncillaryBgos() {
    //loop over all bgo's and remove all matching germaniums
    for(auto bgo = fAncillaryBgoCrystal->begin(); bgo != fAncillaryBgoCrystal->end(); ++bgo) {
        for(auto ge = fGriffinCrystal->begin(); ge != fGriffinCrystal->end();) {
            if(GriffinAncillaryBgoNeighbours_det[ge->DetectorId()][0] == bgo->DetectorId() && GriffinAncillaryBgoNeighbours_cry[ge->DetectorId()][0] == bgo->CrystalId() ) {
                ge = fGriffinCrystal->erase(ge);
            }
            else if(GriffinAncillaryBgoNeighbours_det[ge->DetectorId()][1] == bgo->DetectorId() && GriffinAncillaryBgoNeighbours_cry[ge->DetectorId()][1] == bgo->CrystalId() ) {
                ge = fGriffinCrystal->erase(ge);
            }
            else {
                ++ge;
            }
        }
    }
}

void Converter::SupressGriffinBySceptar() {
    //loop over all bgo's and remove all matching germaniums
    for(auto bgo = fSceptarDetector->begin(); bgo != fSceptarDetector->end(); ++bgo) {

        for(auto ge = fGriffinCrystal->begin(); ge != fGriffinCrystal->end();) {

            if(GriffinSceptarSuppressors_det[ge->DetectorId()][0] == bgo->DetectorId() ) {
                ge = fGriffinCrystal->erase(ge);
            }
            else if(GriffinSceptarSuppressors_det[ge->DetectorId()][1] == bgo->DetectorId() ) {
                ge = fGriffinCrystal->erase(ge);
            }
            else if(GriffinSceptarSuppressors_det[ge->DetectorId()][2] == bgo->DetectorId() ) {
                ge = fGriffinCrystal->erase(ge);
            }
            else if(GriffinSceptarSuppressors_det[ge->DetectorId()][3] == bgo->DetectorId() ) {
                ge = fGriffinCrystal->erase(ge);
            }
            else {
                ++ge;
            }
        }
    }
}

void Converter::AddbackGriffin() {
    std::vector<Detector>::iterator detector;
    for(auto crystal = fGriffinCrystal->begin(); crystal != fGriffinCrystal->end(); ++crystal) {
        //try and find a matching detector to add this crystals energy to
        for(detector = fGriffinDetector->begin(); detector != fGriffinDetector->end(); ++detector) {
            if(crystal->DetectorId() == detector->DetectorId()) {
                detector->AddEnergy(crystal->SimulationEnergy(),crystal->Energy());
                break;
            }
        }
        //if the above loop ended w/o finding a matching detector, we create a new one
        if(detector == fGriffinDetector->end()) {
            fGriffinDetector->push_back(*crystal);
        }
        if(fGriffinArray->size() == 0) {
            fGriffinArray->push_back(*crystal);
        }
        else {
            fGriffinArray->at(0).AddEnergy(crystal->SimulationEnergy(),crystal->Energy());
        }
    }

    // Do the neighbour add-back;
    AddbackGriffinNeighbour();
}

void Converter::AddbackGriffinNeighbour() {
    // This is the idealized neighbour addback method, as we know the order of the Geant4 output.
    // Any scattering should be in order in the output. We could do other things here.
    // For example, we could order the energies from high to low and group them that way.
    std::vector<Detector>::iterator neighbour;
    for(auto crystal = fGriffinCrystal->begin(); crystal != fGriffinCrystal->end(); ++crystal) {
        for(neighbour= fGriffinNeighbour->begin(); neighbour != fGriffinNeighbour->end(); ++neighbour) {
            if(crystal->DetectorId() == neighbour->DetectorId() || GriffinNeighbours_det[neighbour->DetectorId()][0] == crystal->DetectorId() || GriffinNeighbours_det[neighbour->DetectorId()][1] == crystal->DetectorId() || GriffinNeighbours_det[neighbour->DetectorId()][2] == crystal->DetectorId() || GriffinNeighbours_det[neighbour->DetectorId()][3] == crystal->DetectorId() || GriffinNeighbours_det[crystal->DetectorId()][0] == neighbour->DetectorId() || GriffinNeighbours_det[crystal->DetectorId()][1] == neighbour->DetectorId() || GriffinNeighbours_det[crystal->DetectorId()][2] == neighbour->DetectorId() || GriffinNeighbours_det[crystal->DetectorId()][3] == neighbour->DetectorId() ) {
                neighbour->AddEnergy(crystal->SimulationEnergy(),crystal->Energy());
                break;
            }
        }
        //if the above loop ended w/o finding a matching detector, we create a new one
        if(neighbour== fGriffinNeighbour->end()) {
            fGriffinNeighbour->push_back(*crystal);
        }
    }
}

void Converter::AddbackGriffinNeighbourVector() {
    std::vector<Detector>::iterator neighbour;

    double energy_ascend[64] = {0};
    double energy_descend[64] = {0};
    int index;
    bool goodenergy;

    for(auto crystal = fGriffinCrystal->begin(); crystal != fGriffinCrystal->end(); ++crystal) {
        index = crystal->CrystalId() + (crystal->DetectorId() * 4);
        energy_ascend[index] = crystal->SimulationEnergy();
    }

    // Sort energies descending
    int elements = sizeof(energy_ascend) / sizeof(energy_ascend[0]);
    std::sort(energy_ascend, energy_ascend + elements);
    for (int i = 0; i < elements; ++i)
        energy_descend[i] = energy_ascend[elements-1-i];

    for (int i = 0; i < 64; ++i) {
        if(energy_descend[i] == 0)
            break;
        // Now we fill the neighbour vector array according to this list of descending energies
        for(auto crystal1 = fGriffinCrystal->begin(); crystal1 != fGriffinCrystal->end(); ++crystal1) {
            if(crystal1->SimulationEnergy() == energy_descend[i]) { // found crystal with energy energy_descend[i]
                //cout << "new energy = " << crystal1->SimulationEnergy() << " cry = " << crystal1->CrystalId() << " det = " << crystal1->DetectorId() << endl;
                fGriffinNeighbourVector->push_back(*crystal1);
                for(auto crystal2 = fGriffinCrystal->begin(); crystal2 != fGriffinCrystal->end(); ++crystal2) { // loop over crystals again

                    if(crystal1 != crystal2 ) { // make sure we don't add the same crystals together!
                        goodenergy = false;
                        // does energy exist in energy_descend?
                        for (int j = i + 1; j < 64; ++j) {
                            if(energy_descend[j] == crystal2->SimulationEnergy()) {
                                goodenergy = true;
                                index = j;
                                break;
                            }
                            if(energy_descend[j] == 0)
                                break;
                        }

                        if(goodenergy && AreGriffinCrystalCenterPositionsWithinVectorLength(crystal1->CrystalId(), crystal1->DetectorId(), crystal2->CrystalId(), crystal2->DetectorId())) {
                            for(neighbour= fGriffinNeighbourVector->begin(); neighbour != fGriffinNeighbourVector->end(); ++neighbour) { // loop over existing neighbours
                                if(crystal1->DetectorId() == neighbour->DetectorId() && crystal1->CrystalId() == neighbour->CrystalId()) {
                                    neighbour->AddEnergy(crystal2->SimulationEnergy(),crystal2->Energy());
                                    //cout << "neighbour energy = " << neighbour->SimulationEnergy() << " cry = " << neighbour->CrystalId() << " det = " << neighbour->DetectorId() << endl;

                                    // now that we have used this crystal2 energy in the neighbour sum, remove it from the list
                                    if(energy_descend[index] == crystal2->SimulationEnergy()) {
                                        for (int k = index; k < 64-1; ++k) {
                                            energy_descend[k] = energy_descend[k+1];
                                            if(energy_descend[k+1] == 0)
                                                break;
                                        }
                                    }
                                    else {
                                        std::cout << " Something is wrong here! index was not found correctly?" << std::endl;
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }

                // now that we have used this crystal1 energy in the neighbour sum, remove it from the list
                if(energy_descend[i] == crystal1->SimulationEnergy()) {
                    for (int k = i; k < 64-1; ++k) {
                        energy_descend[k] = energy_descend[k+1];
                        if(energy_descend[k+1] == 0)
                            break;
                    }
                }
                else {
                    std::cout << " Something is wrong here! index was not found correctly?" << std::endl;
                }
            }
        }
    }

}

void Converter::CheckLaBrDetectorAddback() {
// This method checks that all the "detector" hits are unique, that is, they have different detector IDs.
// If they have the same detector IDs, then we sum the energies together.
// Normally the Geant4 simulation would sum energy deposits on the same volume, but if we ran the code in "step mode",
// or if we merged two ntuples together, this would not be true. This method checks and will do what "hit mode" in Geant4 normally does for us!
    for(auto firstDet = fLaBrDetector->begin(); firstDet != fLaBrDetector->end(); ++firstDet) {
        for(auto secondDet = firstDet+1; secondDet != fLaBrDetector->end();) {
            if((firstDet->DetectorId() == secondDet->DetectorId())) {
                firstDet->AddEnergy(secondDet->SimulationEnergy(),secondDet->Energy());
                secondDet = fLaBrDetector->erase(secondDet);
            }
            else {
                ++secondDet;
            }
        }
    }
}

void Converter::SupressLaBr() {
    //loop over all bgo's and remove all matching germaniums
    for(auto bgo = fAncillaryBgoDetector->begin(); bgo != fAncillaryBgoDetector->end(); ++bgo) {
        for(auto ge = fLaBrDetector->begin(); ge != fLaBrDetector->end();) {
            if(bgo->DetectorId() == ge->DetectorId()) {
                ge = fLaBrDetector->erase(ge);
            } else {
                ++ge;
            }
        }
    }
}

void Converter::CheckEightPiDetectorAddback() {
// This method checks that all the "detector" hits are unique, that is, they have different detector IDs.
// If they have the same detector IDs, then we sum the energies together.
// Normally the Geant4 simulation would sum energy deposits on the same volume, but if we ran the code in "step mode",
// or if we merged two ntuples together, this would not be true. This method checks and will do what "hit mode" in Geant4 normally does for us!
    for(auto firstDet = fEightPiDetector->begin(); firstDet != fEightPiDetector->end(); ++firstDet) {
        for(auto secondDet = firstDet+1; secondDet != fEightPiDetector->end();) {
            if((firstDet->DetectorId() == secondDet->DetectorId())) {
                firstDet->AddEnergy(secondDet->SimulationEnergy(),secondDet->Energy());
                secondDet = fEightPiDetector->erase(secondDet);
            }
            else {
                ++secondDet;
            }
        }
    }
}

void Converter::SupressEightPi() {
    //loop over all bgo's and remove all matching germaniums
    for(auto bgo = fEightPiBgoDetector->begin(); bgo != fEightPiBgoDetector->end(); ++bgo) {
        for(auto ge = fEightPiDetector->begin(); ge != fEightPiDetector->end();) {
            if(bgo->DetectorId() == ge->DetectorId()) {
                ge = fEightPiDetector->erase(ge);
            } else {
                ++ge;
            }
        }
    }
}

void Converter::SupressLaBrByNeighbouringGriffinShields() {
    //loop over all bgo's and remove all matching germaniums
    for(auto bgo = fGriffinBgo->begin(); bgo != fGriffinBgo->end(); ++bgo) {
        for(auto ge = fLaBrDetector->begin(); ge != fLaBrDetector->end();) {
            if(LaBrGriffinNeighbours_det[ge->DetectorId()][0] == bgo->DetectorId() && LaBrGriffinNeighbours_cry[ge->DetectorId()][0] == bgo->CrystalId() ) {
                ge = fLaBrDetector->erase(ge);
            }
            else if(LaBrGriffinNeighbours_det[ge->DetectorId()][1] == bgo->DetectorId() && LaBrGriffinNeighbours_cry[ge->DetectorId()][1] == bgo->CrystalId() ) {
                ge = fLaBrDetector->erase(ge);
            }
            else if(LaBrGriffinNeighbours_det[ge->DetectorId()][2] == bgo->DetectorId() && LaBrGriffinNeighbours_cry[ge->DetectorId()][2] == bgo->CrystalId() ) {
                ge = fLaBrDetector->erase(ge);
            }
            else {
                ++ge;
            }
        }
    }
}

void Converter::AddbackLaBr() {
    for(auto detector = fLaBrDetector->begin(); detector != fLaBrDetector->end(); ++detector) {
        if(fLaBrArray->size() == 0) {
            fLaBrArray->push_back(*detector);
        } else {
            fLaBrArray->at(0).AddEnergy(detector->SimulationEnergy(),detector->Energy());
        }
    }
}

void Converter::AddbackEightPi() {
    for(auto detector = fEightPiDetector->begin(); detector != fEightPiDetector->end(); ++detector) {
        if(fEightPiArray->size() == 0) {
            fEightPiArray->push_back(*detector);
        } else {
            fEightPiArray->at(0).AddEnergy(detector->SimulationEnergy(),detector->Energy());
        }
    }
}

void Converter::CheckAncillaryBgoCrystalAddback() {
// This method checks that all the "crystal" hits are unique, that is, they have different crystal and detector IDs.
// If they have the same crystal and detector IDs, then we sum the energies together.
// Normally the Geant4 simulation would sum energy deposits on the same volume, but if we ran the code in "step mode",
// or if we merged two ntuples together, this would not be true. This method checks and will do what "hit mode" in Geant4 normally does for us!
    for(auto firstDet = fAncillaryBgoCrystal->begin(); firstDet != fAncillaryBgoCrystal->end(); ++firstDet) {
        for(auto secondDet = firstDet+1; secondDet != fAncillaryBgoCrystal->end();) {
            if((firstDet->DetectorId() == secondDet->DetectorId()) && (firstDet->CrystalId() == secondDet->CrystalId())) {
                firstDet->AddEnergy(secondDet->SimulationEnergy(),secondDet->Energy());
                secondDet = fAncillaryBgoCrystal->erase(secondDet);
            }
            else {
                ++secondDet;
            }
        }
    }
}

void Converter::AddbackAncillaryBgo() {
    std::vector<Detector>::iterator detector;
    for(auto crystal = fAncillaryBgoCrystal->begin(); crystal != fAncillaryBgoCrystal->end(); ++crystal) {
        //try and find a matching detector to add this crystals energy to
        for(detector = fAncillaryBgoDetector->begin(); detector != fAncillaryBgoDetector->end(); ++detector) {
            if(crystal->DetectorId() == detector->DetectorId()) {
                detector->AddEnergy(crystal->SimulationEnergy(),crystal->Energy());
                break;
            }
        }
        //if the above loop ended w/o finding a matching detector, we create a new one
        if(detector == fAncillaryBgoDetector->end()) {
            fAncillaryBgoDetector->push_back(*crystal);
        }
        if(fAncillaryBgoArray->size() == 0 ) {
            fAncillaryBgoArray->push_back(*crystal);
        } else {
            fAncillaryBgoArray->at(0).AddEnergy(crystal->SimulationEnergy(),crystal->Energy());
        }
    }
}

void Converter::CheckSceptarDetectorAddback() {
// This method checks that all the "detector" hits are unique, that is, they have different detector IDs.
// If they have the same detector IDs, then we sum the energies together.
// Normally the Geant4 simulation would sum energy deposits on the same volume, but if we ran the code in "step mode",
// or if we merged two ntuples together, this would not be true. This method checks and will do what "hit mode" in Geant4 normally does for us!
    for(auto firstDet = fSceptarDetector->begin(); firstDet != fSceptarDetector->end(); ++firstDet) {
        for(auto secondDet = firstDet+1; secondDet != fSceptarDetector->end();) {
            if((firstDet->DetectorId() == secondDet->DetectorId())) {
                firstDet->AddEnergy(secondDet->SimulationEnergy(),secondDet->Energy());
                secondDet = fSceptarDetector->erase(secondDet);
            }
            else {
                ++secondDet;
            }
        }
    }
}

void Converter::AddbackSceptar() {
    for(auto detector = fSceptarDetector->begin(); detector != fSceptarDetector->end(); ++detector) {
        if(fSceptarArray->size() == 0) {
            fSceptarArray->push_back(*detector);
        } else {
            fSceptarArray->at(0).AddEnergy(detector->SimulationEnergy(),detector->Energy());
        }
    }
}

void Converter::CheckDescantDetectorAddback() {
// This method checks that all the "detector" hits are unique, that is, they have different detector IDs.
// If they have the same detector IDs, then we sum the energies together.
// Normally the Geant4 simulation would sum energy deposits on the same volume, but if we ran the code in "step mode",
// or if we merged two ntuples together, this would not be true. This method checks and will do what "hit mode" in Geant4 normally does for us!
    for(auto firstDet = fDescantBlueDetector->begin(); firstDet != fDescantBlueDetector->end(); ++firstDet) {
        for(auto secondDet = firstDet+1; secondDet != fDescantBlueDetector->end();) {
            if((firstDet->DetectorId() == secondDet->DetectorId())) {
                firstDet->AddEnergy(secondDet->SimulationEnergy(),secondDet->Energy());
                secondDet = fDescantBlueDetector->erase(secondDet);
            }
            else {
                ++secondDet;
            }
        }
    }
    for(auto firstDet = fDescantGreenDetector->begin(); firstDet != fDescantGreenDetector->end(); ++firstDet) {
        for(auto secondDet = firstDet+1; secondDet != fDescantGreenDetector->end();) {
            if((firstDet->DetectorId() == secondDet->DetectorId())) {
                firstDet->AddEnergy(secondDet->SimulationEnergy(),secondDet->Energy());
                secondDet = fDescantGreenDetector->erase(secondDet);
            }
            else {
                ++secondDet;
            }
        }
    }
    for(auto firstDet = fDescantRedDetector->begin(); firstDet != fDescantRedDetector->end(); ++firstDet) {
        for(auto secondDet = firstDet+1; secondDet != fDescantRedDetector->end();) {
            if((firstDet->DetectorId() == secondDet->DetectorId())) {
                firstDet->AddEnergy(secondDet->SimulationEnergy(),secondDet->Energy());
                secondDet = fDescantRedDetector->erase(secondDet);
            }
            else {
                ++secondDet;
            }
        }
    }
    for(auto firstDet = fDescantWhiteDetector->begin(); firstDet != fDescantWhiteDetector->end(); ++firstDet) {
        for(auto secondDet = firstDet+1; secondDet != fDescantWhiteDetector->end();) {
            if((firstDet->DetectorId() == secondDet->DetectorId())) {
                firstDet->AddEnergy(secondDet->SimulationEnergy(),secondDet->Energy());
                secondDet = fDescantWhiteDetector->erase(secondDet);
            }
            else {
                ++secondDet;
            }
        }
    }
    for(auto firstDet = fDescantYellowDetector->begin(); firstDet != fDescantYellowDetector->end(); ++firstDet) {
        for(auto secondDet = firstDet+1; secondDet != fDescantYellowDetector->end();) {
            if((firstDet->DetectorId() == secondDet->DetectorId())) {
                firstDet->AddEnergy(secondDet->SimulationEnergy(),secondDet->Energy());
                secondDet = fDescantYellowDetector->erase(secondDet);
            }
            else {
                ++secondDet;
            }
        }
    }
}

void Converter::AddbackDescant() {
    for(auto detector = fDescantBlueDetector->begin(); detector != fDescantBlueDetector->end(); ++detector) {
        if(fDescantArray->size() == 0) {
            fDescantArray->push_back(*detector);
        } else {
            fDescantArray->at(0).AddEnergy(detector->SimulationEnergy(),detector->Energy());
            if(fDescantArray->at(0).Time() < detector->Time()) { // If added energy has a later time, set array time to this
                fDescantArray->at(0).SetTime(detector->Time());
            }
        }
    }
    for(auto detector = fDescantGreenDetector->begin(); detector != fDescantGreenDetector->end(); ++detector) {
        if(fDescantArray->size() == 0) {
            fDescantArray->push_back(*detector);
        } else {
            fDescantArray->at(0).AddEnergy(detector->SimulationEnergy(),detector->Energy());
            if(fDescantArray->at(0).Time() < detector->Time()) { // If added energy has a later time, set array time to this
                fDescantArray->at(0).SetTime(detector->Time());
            }
        }
    }
    for(auto detector = fDescantRedDetector->begin(); detector != fDescantRedDetector->end(); ++detector) {
        if(fDescantArray->size() == 0) {
            fDescantArray->push_back(*detector);
        } else {
            fDescantArray->at(0).AddEnergy(detector->SimulationEnergy(),detector->Energy());
            if(fDescantArray->at(0).Time() < detector->Time()) { // If added energy has a later time, set array time to this
                fDescantArray->at(0).SetTime(detector->Time());
            }
        }
    }
    for(auto detector = fDescantWhiteDetector->begin(); detector != fDescantWhiteDetector->end(); ++detector) {
        if(fDescantArray->size() == 0) {
            fDescantArray->push_back(*detector);
        } else {
            fDescantArray->at(0).AddEnergy(detector->SimulationEnergy(),detector->Energy());
            if(fDescantArray->at(0).Time() < detector->Time()) { // If added energy has a later time, set array time to this
                fDescantArray->at(0).SetTime(detector->Time());
            }
        }
    }
    for(auto detector = fDescantYellowDetector->begin(); detector != fDescantYellowDetector->end(); ++detector) {
        if(fDescantArray->size() == 0) {
            fDescantArray->push_back(*detector);
        } else {
            fDescantArray->at(0).AddEnergy(detector->SimulationEnergy(),detector->Energy());
            if(fDescantArray->at(0).Time() < detector->Time()) { // If added energy has a later time, set array time to this
                fDescantArray->at(0).SetTime(detector->Time());
            }
        }
    }
}


void Converter::CheckPacesDetectorAddback() {
// This method checks that all the "detector" hits are unique, that is, they have different detector IDs.
// If they have the same detector IDs, then we sum the energies together.
// Normally the Geant4 simulation would sum energy deposits on the same volume, but if we ran the code in "step mode",
// or if we merged two ntuples together, this would not be true. This method checks and will do what "hit mode" in Geant4 normally does for us!
    for(auto firstDet = fPacesDetector->begin(); firstDet != fPacesDetector->end(); ++firstDet) {
        for(auto secondDet = firstDet+1; secondDet != fPacesDetector->end();) {
            if((firstDet->DetectorId() == secondDet->DetectorId())) {
                firstDet->AddEnergy(secondDet->SimulationEnergy(),secondDet->Energy());
                secondDet = fPacesDetector->erase(secondDet);
            }
            else {
                ++secondDet;
            }
        }
    }
}

void Converter::AddbackPaces() {
    for(auto detector = fPacesDetector->begin(); detector != fPacesDetector->end(); ++detector) {
        if(fPacesArray->size() == 0) {
            fPacesArray->push_back(*detector);
        } else {
            fPacesArray->at(0).AddEnergy(detector->SimulationEnergy(),detector->Energy());
        }
    }
}

void Converter::PrintStatistics() {
}

TH1F* Converter::Get1DHistogram(std::string histogramName, std::string directoryName, int nbins, double lowbin, double highbin) {
    //try and find this histogram
    TH1F* hist = (TH1F*) gDirectory->FindObjectAny(histogramName.c_str());

    if(hist == nullptr){
        //if the histogram doesn't exist, we create it and add it to the histogram list
        hist = new TH1F(histogramName.c_str(),histogramName.c_str(),nbins,lowbin,highbin);
        if(fHistograms.find(directoryName) == fHistograms.end()) {
            fHistograms[directoryName] = new TList;
        }
        fHistograms[directoryName]->Add((TObject*) hist);
    }
    return hist;
}

TH1F* Converter::Get1DHistogram(std::string histogramName, std::string directoryName) {
    //try and find this histogram
    TH1F* hist = (TH1F*) gDirectory->FindObjectAny(histogramName.c_str());

    if(hist == nullptr){
        //if the histogram doesn't exist, we create it and add it to the histogram list
        hist = new TH1F(histogramName.c_str(),histogramName.c_str(),fSettings->NofBins(directoryName),fSettings->RangeLow(directoryName),fSettings->RangeHigh(directoryName));
        if(fHistograms.find(directoryName) == fHistograms.end()) {
            fHistograms[directoryName] = new TList;
        }
        fHistograms[directoryName]->Add((TObject*) hist);
    }
    return hist;
}

TH2F* Converter::Get2DHistogram(std::string histogramName, std::string directoryName, int nbinsX, double lowbinX, double highbinX, int nbinsY, double lowbinY, double highbinY) {
    //try and find this histogram
    TH2F* hist = (TH2F*) gDirectory->FindObjectAny(histogramName.c_str());
    if(hist == nullptr){
        //if the histogram doesn't exist, we create it and add it to the histogram list
        hist = new TH2F(histogramName.c_str(),histogramName.c_str(),
                        nbinsX, lowbinX, highbinX,
                        nbinsY, lowbinY, highbinY);
        if(fHistograms.find(directoryName) == fHistograms.end()) {
            fHistograms[directoryName] = new TList;
        }
        fHistograms[directoryName]->Add((TObject*) hist);
    }
    return hist;
}

TH2F* Converter::Get2DHistogram(std::string histogramName, std::string directoryName) {
    //try and find this histogram
    TH2F* hist = (TH2F*) gDirectory->FindObjectAny(histogramName.c_str());
    if(hist == nullptr){
        //if the histogram doesn't exist, we create it and add it to the histogram list
        hist = new TH2F(histogramName.c_str(),histogramName.c_str(),
                        fSettings->NofBins(directoryName),fSettings->RangeLow(directoryName),fSettings->RangeHigh(directoryName),
                        fSettings->NofBins(directoryName),fSettings->RangeLow(directoryName),fSettings->RangeHigh(directoryName));
        if(fHistograms.find(directoryName) == fHistograms.end()) {
            fHistograms[directoryName] = new TList;
        }
        fHistograms[directoryName]->Add((TObject*) hist);
    }
    return hist;
}

/*TH3I* Converter::Get3DHistogram(std::string histogramName, std::string directoryName) {
    //try and find this histogram
    TH3I* hist = (TH3I*) gDirectory->FindObjectAny(histogramName.c_str());
    if(hist == nullptr){
        //if the histogram doesn't exist, we create it and add it to the histogram list
        hist = new TH3I(histogramName.c_str(),histogramName.c_str(),
                        fSettings->NofBins(directoryName),fSettings->RangeLow(directoryName),fSettings->RangeHigh(directoryName),
                        fSettings->NofBins(directoryName),fSettings->RangeLow(directoryName),fSettings->RangeHigh(directoryName),
                        52,0,52);
        if(fHistograms.find(directoryName) == fHistograms.end()) {
            fHistograms[directoryName] = new TList;
        }
        fHistograms[directoryName]->Add((TObject*) hist);
    }
    return hist;
}*/

THnSparseF* Converter::GetNDHistogram(std::string histogramName, std::string directoryName) {
    //try and find this histogram
    //This method is different for THnSparse if implemented normally with the other method
    //then the histogram would not set the address of fHistograms[directoryName]->FindObject(histogramName.c_str()) and a new histogram would be created for every event
    THnSparseF* hist = nullptr;
    if(fHistograms.find(directoryName) == fHistograms.end() || fHistograms[directoryName]->FindObject(histogramName.c_str()) == nullptr) {
        //if the histogram doesn't exist, we create it and add it to the histogram list
        Double_t min[3] = {fSettings->RangeLow(directoryName), fSettings->RangeLow(directoryName), 0};
        Double_t max[3] = {fSettings->RangeHigh(directoryName), fSettings->RangeHigh(directoryName), 52};
        Int_t Bins[3] = {fSettings->NofBins(directoryName), fSettings->NofBins(directoryName), 52};
        hist = new THnSparseF(histogramName.c_str(),histogramName.c_str(),3, Bins,min,max);
        if(fHistograms.find(directoryName) == fHistograms.end()) {
            fHistograms[directoryName] = new TList;
        }
        fHistograms[directoryName]->Add((TObject*) hist);
    } else {
        hist = (THnSparseF*) fHistograms[directoryName]->FindObject(histogramName.c_str());
    }
    
    return hist;
}

void Converter::FillHistDetector1DGamma(TH1F* hist1D, std::vector<Detector>* detector, std::string hist_name, std::string hist_dir) {
    for(size_t firstDet = 0; firstDet < detector->size(); ++firstDet) {
        hist1D = Get1DHistogram(hist_name,hist_dir);
        hist1D->Fill(detector->at(firstDet).Energy());
    }
}

void Converter::FillHistDetector2DGammaGamma(TH2F* hist2D, std::vector<Detector>* detector, std::string hist_name, std::string hist_dir) {
    if(fSettings->Write2DHist()) {
        for(size_t firstDet = 0; firstDet < detector->size(); ++firstDet) {
            for(size_t secondDet = firstDet+1; secondDet < detector->size(); ++secondDet) {
                hist2D = Get2DHistogram(hist_name,hist_dir);
                // symmetrize!
                hist2D->Fill(detector->at(firstDet).Energy(),detector->at(secondDet).Energy());
                hist2D->Fill(detector->at(secondDet).Energy(),detector->at(firstDet).Energy());
            }
        }
    }
}

void Converter::FillHistDetector1DGammaNR(TH1F* hist1D, std::vector<Detector>* detector, std::string hist_name, std::string hist_dir) {
    for(size_t firstDet = 0; firstDet < detector->size(); ++firstDet) {
        hist1D = Get1DHistogram(hist_name,hist_dir); //
        hist1D->Fill(detector->at(firstDet).SimulationEnergy());
    }
}

void Converter::FillHistDetector2DGammaGammaNR(TH2F* hist2D, std::vector<Detector>* detector, std::string hist_name, std::string hist_dir) {
    if(fSettings->Write2DHist()) {
        for(size_t firstDet = 0; firstDet < detector->size(); ++firstDet) {
            for(size_t secondDet = firstDet+1; secondDet < detector->size(); ++secondDet) {
                hist2D = Get2DHistogram(hist_name,hist_dir);
                // symmetrize!
                hist2D->Fill(detector->at(firstDet).SimulationEnergy(),detector->at(secondDet).SimulationEnergy());
                hist2D->Fill(detector->at(secondDet).SimulationEnergy(),detector->at(firstDet).SimulationEnergy());
            }
        }
    }
}

void Converter::FillHist2DGriffinSceptarHitPattern(TH2F* hist2D, std::vector<Detector>* detector1, std::vector<Detector>* detector2, std::string hist_name, std::string hist_dir) {
    if(fSettings->Write2DHist()) {
        for(size_t firstDet = 0; firstDet < detector1->size(); ++firstDet) {
            for(size_t secondDet = 0; secondDet < detector2->size(); ++secondDet) {
                hist2D = Get2DHistogram(hist_name,hist_dir);
                hist2D->Fill(detector2->at(secondDet).DetectorId(),(4*detector1->at(firstDet).DetectorId()+detector1->at(firstDet).CrystalId()));
            }
        }
    }
}

void Converter::FillHist2DGriffinHitPattern(TH2F* hist2D, std::vector<Detector>* detector, std::string hist_name, std::string hist_dir) {
    if(fSettings->Write2DHist()) {
        for(size_t firstDet = 0; firstDet < detector->size(); ++firstDet) {
            hist2D = Get2DHistogram(hist_name,hist_dir);
            hist2D->Fill((4*detector->at(firstDet).DetectorId()+detector->at(firstDet).CrystalId()),(4*detector->at(firstDet).DetectorId()+detector->at(firstDet).CrystalId()));
            for(size_t secondDet = firstDet+1; secondDet < detector->size(); ++secondDet) {
                hist2D->Fill((4*detector->at(firstDet).DetectorId()+detector->at(firstDet).CrystalId()),(4*detector->at(secondDet).DetectorId()+detector->at(secondDet).CrystalId()));
            }
        }
    }
}


TVector3 Converter::GriffinCrystalCenterPosition(int cry, int det) {

    double theta   = GriffinDetCoords[det][0]*(M_PI/180);
    double phi     = GriffinDetCoords[det][1]*(M_PI/180);

//    double germanium_width                 = 56.5; // mm
//    double germanium_separation            = 0.6; // mm
//    double germanium_length                = 90.0; // mm
    double germanium_dist_from_can_face 	 = 5.5; // mm
    double can_face_thickness              = 1.5; //mm
    double distance_to_can_face            = fSettings->GriffinAddbackVectorCrystalFaceDistancemm(); // mm
//    double germanium_shift                 = 1.05; // mm
//    double germanium_bent_length           = 36.2; // mm
//    double bent_end_angle                  = 22.5*M_PI/180.0; // rad
    // germanium_shift comment from GRIFFIN Geant4 code:
    // this can't be more than 2.75mm. It is the amount by which
    // one side is cut closer to the center than the other
    // the ending length of the cones

    double  depth   = fSettings->GriffinAddbackVectorDepthmm();
//    double  offset  = germanium_bent_length*tan(bent_end_angle)/2.0;
    // this offset is to push the center of the crystal face towards the
    // center of the clover. We do this because the outter edges of the crystal
    // are tappered. We'll take half this value, which is about 7.5 mm.

    // double x0 = (germanium_width + germanium_separation)/2.0 - germanium_shift - offset;
    // double y0 = (germanium_width + germanium_separation)/2.0 - germanium_shift - offset;
    // We push the x and y directions towards the center of the clover by
    // germanium_shift. This is keeping things symmetric along the diagonal,
    // but in reality this doesn`t need to be the case if the core contact is
    // anti-symmetric.

    // From Andrews work, he says the center is x = y = 26 mm.
    double x0 = 26.0; // mm
    double y0 = 26.0; // mm

    double z0 = distance_to_can_face + can_face_thickness + germanium_dist_from_can_face + depth;

    double i = (double)(cry);
    double x = -1*((x0*(pow((-1),(floor((i+1.0)/2.0))))));
    double y = -1*((y0*(pow((-1),(floor((i+2.0)/2.0))))));
    double z = z0;

    TVector3 vec(transX(x,y,z,theta,phi),transY(x,y,z,theta,phi),transZ(x,y,z,theta,phi));
    return vec;
}

bool Converter::AreGriffinCrystalCenterPositionsWithinVectorLength(int cry1, int det1, int cry2, int det2){
    TVector3 vec1 = GriffinCrystalCenterPosition(cry1,det1);
    TVector3 vec2 = GriffinCrystalCenterPosition(cry2,det2);

    double x1 = vec1.X();
    double y1 = vec1.Y();
    double z1 = vec1.Z();
    double x2 = vec2.X();
    double y2 = vec2.Y();
    double z2 = vec2.Z();

    bool result = false;

    if( ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2)) <= fSettings->GriffinAddbackVectorLengthmm()*fSettings->GriffinAddbackVectorLengthmm() ) {
        result = true;
    }

    return result;
}

double Converter::transX(double x, double y, double z, double theta, double phi){
    return (x*cos(theta)+z*sin(theta))*cos(phi)-y*sin(phi);
}

double Converter::transY(double x, double y, double z, double theta, double phi){
    return (x*cos(theta)+z*sin(theta))*sin(phi)+y*cos(phi);
}

double Converter::transZ(double x, double y, double z, double theta, double phi){
    return -x*sin(theta)+z*cos(theta);
}

// ------------- Tistar methods ------------- //

int Converter::CalculateTistarStripNumber(int layerNb, TVector3 particlePos, TVector3 stripPos, TVector3 stripDim) {
    int stripNb = -1;
    TVector3 localPos = particlePos - stripPos;
    //double z = localPos.z() + stripDim.z()/2.;

    double z;
    if(stripPos.z() > 0.) { // forwards
        z = localPos.z() + stripDim.z()/2.;
    } else {                // backwards
        z = localPos.z() - stripDim.z()/2.;
    }
    
    z = fabs(z);
    if(z == stripDim.z()) {
        z -= 1e-5;
    }
    
    stripNb = static_cast<int>(z/fSettings->GetTISTARStripWidthZ(layerNb));

    if(stripNb > static_cast<int>(fSettings->GetTistarSettings()->GetLayerDimensionVector()[layerNb][0].z()/fSettings->GetTISTARStripWidthZ(layerNb)) || stripNb < 0) {
        std::cout<<"Problem: localZ = "<<z<<" , stripNb = "<<stripNb<<std::endl;
        return -1;
    }

    //std::cout << "calculate strip number: " << std::endl;
    //std::cout << "particlePos = "; particlePos.Print();
    //std::cout << "localPos =    "; localPos.Print();
    //std::cout << "stripPos =    "; stripPos.Print()   ;
    //std::cout << "stripDim =    "; stripDim.Print()   ;
    //std::cout << "layerNb =     " << layerNb << std::endl;
    //std::cout << "z =           " << z << std::endl;
    //std::cout << "stripWidthZ = " << fSettings->GetTISTARStripWidthZ(layerNb) << std::endl;
    //std::cout << "stripNb =     " << stripNb << std::endl << std::endl;

    return stripNb;
}

int Converter::CalculateTistarRingNumber(int layerNb, TVector3 particlePos, TVector3 stripPos, TVector3 stripDim) {
    int ringNb = -1;
    TVector3 localPos = particlePos - stripPos;
    double y = localPos.y() + stripDim.y()/2.;

    if(y == stripDim.y()) {
        y -= 1e-5;
    }

    ringNb = static_cast<int>(y/fSettings->GetTISTARStripWidthY(layerNb));
    
    if(ringNb > static_cast<int>(fSettings->GetTistarSettings()->GetLayerDimensionVector()[layerNb][0].y()/fSettings->GetTISTARStripWidthY(layerNb)) || ringNb < 0) {
        std::cout<<"Problem: localY = "<<y<<" , ringNb = "<<ringNb<<std::endl;
        return -1;
    }

    //std::cout << "calculate ring number: " << std::endl;
    //std::cout << "particlePos = "; particlePos.Print();
    //std::cout << "localPos =    "; localPos.Print();
    //std::cout << "stripPos =    "; stripPos.Print()   ;
    //std::cout << "stripDim =    "; stripDim.Print()   ;
    //std::cout << "layerNb =     " << layerNb << std::endl;
    //std::cout << "y =           " << y << std::endl;
    //std::cout << "stripWidthY = " << fSettings->GetTISTARStripWidthY(layerNb) << std::endl;
    //std::cout << "ringNb =     " << ringNb << std::endl << std::endl;

    return ringNb;
}

void Converter::FillTistarVectors() {
    if(fSystemID == 9500) {
        if(fDepEnergy < 10.) return; // try putting a low energy cut on which registers as a hit
        TVector3 particlePos = TVector3(fPosx, fPosy, fPosz);
        TVector3 stripPos = TVector3(fSettings->GetTistarSettings()->GetLayerPositionVector()[fDetNumber][fCryNumber]);
        TVector3 stripDim = TVector3(fSettings->GetTistarSettings()->GetLayerDimensionVector()[fDetNumber][fCryNumber]);
        ParticleMC part;
        int stripNumber, ringNumber;
        switch(fDetNumber) {
            size_t j;
            case 0: // first layer - pixelated w/ 4 panels
                // strip calculations
                stripNumber = CalculateTistarStripNumber(fDetNumber, particlePos, stripPos, stripDim);
                for(j = 0; j < fTISTARFirstLayerStripNb[fCryNumber].size(); j++) { // checking to see if this strip is already activated
                    if(stripNumber == fTISTARFirstLayerStripNb[fCryNumber][j]) { // we found a strip that was already activated!
                        fTISTARFirstLayerStripEnergy[fCryNumber][j] += fDepEnergy; 
                        break;
                    }
                }
                if(j == fTISTARFirstLayerStripNb[fCryNumber].size()) { // we didn't find an activated strip, so we add one
                    fTISTARFirstLayerStripNb[fCryNumber].push_back(stripNumber);
                    fTISTARFirstLayerStripEnergy[fCryNumber].push_back(fDepEnergy);
                    fTISTARFirstLayerStripA[fCryNumber].push_back(fTargetA);
                    fTISTARFirstLayerStripZ[fCryNumber].push_back(fTargetZ);
                    fTISTARFirstLayerStripTrackID[fCryNumber].push_back(fTrackID);
                    fTISTARFirstLayerStripTime[fCryNumber].push_back(fTime);
                    fTISTARFirstLayerStripStopped[fCryNumber].push_back(-1); // might need to actually calculate this using the IsStopped method from TRexBarrelDeltaESingleSensitiveDetector
                    fTISTARFirstLayerStripPos[fCryNumber].push_back(particlePos);
                }
                // ring calculations
                ringNumber = CalculateTistarRingNumber(fDetNumber, particlePos, stripPos, stripDim);
                for(j = 0; j < fTISTARFirstLayerRingNb[fCryNumber].size(); j++) { // checking to see if this Ring is already activated
                    if(ringNumber == fTISTARFirstLayerRingNb[fCryNumber][j]) { // we found a Ring that was already activated!
                        fTISTARFirstLayerRingEnergy[fCryNumber][j] += fDepEnergy; 
                        break;
                    }
                }
                if(j == fTISTARFirstLayerRingNb[fCryNumber].size()) { // we didn't find an activated Ring, so we add one
                    fTISTARFirstLayerRingNb[fCryNumber].push_back(ringNumber);
                    fTISTARFirstLayerRingEnergy[fCryNumber].push_back(fDepEnergy);
                    fTISTARFirstLayerRingA[fCryNumber].push_back(fTargetA);
                    fTISTARFirstLayerRingZ[fCryNumber].push_back(fTargetZ);
                    fTISTARFirstLayerRingTrackID[fCryNumber].push_back(fTrackID);
                    fTISTARFirstLayerRingTime[fCryNumber].push_back(fTime);
                    fTISTARFirstLayerRingStopped[fCryNumber].push_back(-1); // might need to actually calculate this using the IsStopped method from TRexBarrelDeltaESingleSensitiveDetector
                }
                break;
    
            case 1: // second layer - pixelated w/ 2 panels
                // strip calculations
                stripNumber = CalculateTistarStripNumber(fDetNumber, particlePos, stripPos, stripDim);
                for(j = 0; j < fTISTARSecondLayerStripNb[fCryNumber].size(); j++) { // checking to see if this strip is already activated
                    if(stripNumber == fTISTARSecondLayerStripNb[fCryNumber][j]) { // we found a strip that was already activated!
                        fTISTARSecondLayerStripEnergy[fCryNumber][j] += fDepEnergy; 
                        break;
                    }
                }
                if(j == fTISTARSecondLayerStripNb[fCryNumber].size()) { // we didn't find an activated strip, so we add one
                    fTISTARSecondLayerStripNb[fCryNumber].push_back(stripNumber);
                    fTISTARSecondLayerStripEnergy[fCryNumber].push_back(fDepEnergy);
                    fTISTARSecondLayerStripA[fCryNumber].push_back(fTargetA);
                    fTISTARSecondLayerStripZ[fCryNumber].push_back(fTargetZ);
                    fTISTARSecondLayerStripTrackID[fCryNumber].push_back(fTrackID);
                    fTISTARSecondLayerStripTime[fCryNumber].push_back(fTime);
                    fTISTARSecondLayerStripStopped[fCryNumber].push_back(-1); // might need to actually calculate this using the IsStopped method from TRexBarrelDeltaESingleSensitiveDetector
                    fTISTARSecondLayerStripPos[fCryNumber].push_back(particlePos);
                }
                // ring calculations
                ringNumber = CalculateTistarRingNumber(fDetNumber, particlePos, stripPos, stripDim);
                for(j = 0; j < fTISTARSecondLayerRingNb[fCryNumber].size(); j++) { // checking to see if this Ring is already activated
                    if(ringNumber == fTISTARSecondLayerRingNb[fCryNumber][j]) { // we found a Ring that was already activated!
                        fTISTARSecondLayerRingEnergy[fCryNumber][j] += fDepEnergy; 
                        break;
                    }
                }
                if(j == fTISTARSecondLayerRingNb[fCryNumber].size()) { // we didn't find an activated Ring, so we add one
                    fTISTARSecondLayerRingNb[fCryNumber].push_back(ringNumber);
                    fTISTARSecondLayerRingEnergy[fCryNumber].push_back(fDepEnergy);
                    fTISTARSecondLayerRingA[fCryNumber].push_back(fTargetA);
                    fTISTARSecondLayerRingZ[fCryNumber].push_back(fTargetZ);
                    fTISTARSecondLayerRingTrackID[fCryNumber].push_back(fTrackID);
                    fTISTARSecondLayerRingTime[fCryNumber].push_back(fTime);
                    fTISTARSecondLayerRingStopped[fCryNumber].push_back(-1); // might need to actually calculate this using the IsStopped method from TRexBarrelDeltaESingleSensitiveDetector
                }
                break;

            case 2: // third (pad) layer - not pixelated, 2 panels
                // here we can directly create and add the new ParticleMC
                ParticleMC part;
                part.SetEdet(fDepEnergy);
                part.SetTime(fTime);
                part.SetA(fTargetA);
                part.SetZ(fTargetZ);
                part.SetTrackID(fTrackID);
                fTISTARPad[fCryNumber]->push_back(part);
                break;
        }    
    }
}

void Converter::ClearTistarVectors() {
    for(int i=0; i<4; i++) {
        fTISTARFirstLayerStripNb[i].clear();
        fTISTARFirstLayerStripEnergy[i].clear();
        fTISTARFirstLayerStripA[i].clear();
        fTISTARFirstLayerStripZ[i].clear();
        fTISTARFirstLayerStripTrackID[i].clear();
        fTISTARFirstLayerStripTime[i].clear();
        fTISTARFirstLayerStripPos[i].clear(); 
        fTISTARFirstLayerStripStopped[i].clear();
        fTISTARFirstLayerRingNb[i].clear();
        fTISTARFirstLayerRingEnergy[i].clear();
        fTISTARFirstLayerRingA[i].clear();
        fTISTARFirstLayerRingZ[i].clear();
        fTISTARFirstLayerRingTrackID[i].clear();
        fTISTARFirstLayerRingTime[i].clear();
        fTISTARFirstLayerRingStopped[i].clear();
    }
    for(int i=0; i<2; i++) {
        fTISTARSecondLayerStripNb[i].clear();
        fTISTARSecondLayerStripEnergy[i].clear();
        fTISTARSecondLayerStripA[i].clear();
        fTISTARSecondLayerStripZ[i].clear();
        fTISTARSecondLayerStripTrackID[i].clear();
        fTISTARSecondLayerStripTime[i].clear();
        fTISTARSecondLayerStripPos[i].clear(); 
        fTISTARSecondLayerStripStopped[i].clear();
        fTISTARSecondLayerRingNb[i].clear();
        fTISTARSecondLayerRingEnergy[i].clear();
        fTISTARSecondLayerRingA[i].clear();
        fTISTARSecondLayerRingZ[i].clear();
        fTISTARSecondLayerRingTrackID[i].clear();
        fTISTARSecondLayerRingTime[i].clear();
        fTISTARSecondLayerRingStopped[i].clear();
    }
}

void Converter::FillTistarParticleMCs() {
    // loop over all first-layer panels
    for(int panelNb = 0; panelNb < 4; panelNb++) {
        ParticleMC part;
        int nStripsAndRings = 0;
        // loop over all strips that have been hit
        for(size_t i = 0; i < fTISTARFirstLayerStripNb[panelNb].size(); i++) {
            part.AddStrip(fTISTARFirstLayerStripNb[panelNb][i],         // strip number
                          fTISTARFirstLayerStripEnergy[panelNb][i],     // energy
                          fTISTARFirstLayerStripA[panelNb][i],          // particle A
                          fTISTARFirstLayerStripZ[panelNb][i],          // particle Z
                          fTISTARFirstLayerStripTrackID[panelNb][i],    // trackID
                          fTISTARFirstLayerStripTime[panelNb][i],       // time
                          fTISTARFirstLayerStripPos[panelNb][i].x(),    // x position
                          fTISTARFirstLayerStripPos[panelNb][i].y(),    // y position
                          fTISTARFirstLayerStripPos[panelNb][i].z(),    // z position
                          fTISTARFirstLayerStripStopped[panelNb][i]);    // is stopped?
            nStripsAndRings++;
        }
        // loop over all rings that have been hit
        for(size_t i = 0; i < fTISTARFirstLayerRingNb[panelNb].size(); i++) {
            part.AddRing(fTISTARFirstLayerRingNb[panelNb][i],           // ring number
                          fTISTARFirstLayerRingEnergy[panelNb][i],      // energy
                          fTISTARFirstLayerRingA[panelNb][i],           // particle A
                          fTISTARFirstLayerRingZ[panelNb][i],           // particle Z
                          fTISTARFirstLayerRingTrackID[panelNb][i],     // trackID
                          fTISTARFirstLayerRingTime[panelNb][i],        // time
                          fTISTARFirstLayerRingStopped[panelNb][i]);    // is stopped?
            nStripsAndRings++;
        }
        // only fill if we have at least one strip/ring activated
        if(nStripsAndRings > 0) {
            fTISTARFirstDeltaE[panelNb]->push_back(part);
        }
    }
    // loop over all second-layer panels
    for(int panelNb = 0; panelNb < 2; panelNb++) {
        ParticleMC part;
        int nStripsAndRings = 0;
        // loop over all strips that have been hit
        for(size_t i = 0; i < fTISTARSecondLayerStripNb[panelNb].size(); i++) {
            part.AddStrip(fTISTARSecondLayerStripNb[panelNb][i],         // strip number
                          fTISTARSecondLayerStripEnergy[panelNb][i],     // energy
                          fTISTARSecondLayerStripA[panelNb][i],          // particle A
                          fTISTARSecondLayerStripZ[panelNb][i],          // particle Z
                          fTISTARSecondLayerStripTrackID[panelNb][i],    // trackID
                          fTISTARSecondLayerStripTime[panelNb][i],       // time
                          fTISTARSecondLayerStripPos[panelNb][i].x(),    // x position
                          fTISTARSecondLayerStripPos[panelNb][i].y(),    // y position
                          fTISTARSecondLayerStripPos[panelNb][i].z(),    // z position
                          fTISTARSecondLayerStripStopped[panelNb][i]);    // is stopped?
            nStripsAndRings++;
        }
        // loop over all rings that have been hit
        for(size_t i = 0; i < fTISTARSecondLayerRingNb[panelNb].size(); i++) {
            part.AddRing(fTISTARSecondLayerRingNb[panelNb][i],           // ring number
                          fTISTARSecondLayerRingEnergy[panelNb][i],      // energy
                          fTISTARSecondLayerRingA[panelNb][i],           // particle A
                          fTISTARSecondLayerRingZ[panelNb][i],           // particle Z
                          fTISTARSecondLayerRingTrackID[panelNb][i],     // trackID
                          fTISTARSecondLayerRingTime[panelNb][i],        // time
                          fTISTARSecondLayerRingStopped[panelNb][i]);    // is stopped?
            nStripsAndRings++;
        }
        if(nStripsAndRings > 0) {
            fTISTARSecondDeltaE[panelNb]->push_back(part);
        }
    }
}

void Converter::ClearTistarParticleMCs() {
    for(int strip=0; strip<4; strip++) { 
        for(auto particle : *(fTISTARFirstDeltaE[strip])) particle.ClearParticleMC();
        fTISTARFirstDeltaE[strip]->clear();
    }
    for(int strip=0; strip<2; strip++) {
        for(auto particle : *(fTISTARSecondDeltaE[strip])) particle.ClearParticleMC();
        fTISTARSecondDeltaE[strip]->clear();
    }
    for(int strip=0; strip<2; strip++) {
        for(auto particle : *(fTISTARPad[strip])) particle.ClearParticleMC();
        fTISTARPad[strip]->clear();
    }
}

void Converter::CreateTistarHistograms(Kinematics * transferP) {
    // we'll create all the histograms from ReconstructSim here
    TistarSettings * sett = fSettings->GetTistarSettings();
    
    std::string directoryName = "TistarAnalysis";
    fHistograms[directoryName.c_str()] = new TList;
    
    TH2F* originXY = new TH2F("originXY", "y vs. x of reconstructed origin", 200, -10., 10., 200, -10., 10.); fHistograms[directoryName.c_str()]->Add(originXY);
    TH2F* originXYErr = new TH2F("originXYErr", "Error y vs. error x of reconstructed origin - simulated origin", 200, -10., 10., 200, -10., 10.); fHistograms[directoryName.c_str()]->Add(originXYErr);
    TH2F* errorOrigin = new TH2F("errorOrigin", "Error between reconstructed and true origin vs. true origin", 200, -100., 100., 1000, -5, 5); fHistograms[directoryName.c_str()]->Add(errorOrigin);
    TH2F* errorThetaPhi = new TH2F("errorThetaPhi", "Error between reconstructed and true phi vs. error in theta", 600, -30, 30, 720, -360., 360.); fHistograms[directoryName.c_str()]->Add(errorThetaPhi);
    TH1F* excEnProton = new TH1F("excEnProton", "Excitaiton Energy Spectrum from reconstructed Protons", 5000, -20000, 20000); fHistograms[directoryName.c_str()]->Add(excEnProton);
    TH1F* reaction = new TH1F("reaction", "Simulated reaction/level", 10, -0.5, 9.5); fHistograms[directoryName.c_str()]->Add(reaction);
    TH2F* phiErrorVsPhi = new TH2F("phiErrorVsPhi","Error in reconstructed #varphi vs. simulated #varphi", 360, -180., 180., 720, -360., 360.); fHistograms[directoryName.c_str()]->Add(phiErrorVsPhi);

    TH2F* dE12VsPad = new TH2F("dE12VsPad", "energy loss first+second layer vs. pad energy", 2000, 0, 50000, 1000, 0, 10000); fHistograms[directoryName.c_str()]->Add(dE12VsPad);
    TH2F* dE12VsE = new TH2F("dE12VsE", "energy loss first+second layer vs. total energy", 2000, 0, 50000, 1000, 0, 10000); fHistograms[directoryName.c_str()]->Add(dE12VsE);
    TH2F* dE1VsE = new TH2F("dE1VsE", "energy loss first layer vs. total energy", 5000, 0, 50000, 1000, 0, 10000); fHistograms[directoryName.c_str()]->Add(dE1VsE);
    TH2F* dE2VsE = new TH2F("dE2VsE", "energy loss second layer vs. total energy", 5000, 0, 50000, 1000, 0, 10000); fHistograms[directoryName.c_str()]->Add(dE2VsE);
    TH2F* dE1VsdE2 = new TH2F("dE1VsdE2", "energy loss second layer vs. energy loss first layer", 1000, 0, 10000, 1000, 0, 10000); fHistograms[directoryName.c_str()]->Add(dE1VsdE2);
    TH2F* eVsTheta = new TH2F("eVsTheta", "recoil energy vs. theta (lab)", 360, 0, 180, 1000, 0, 25000); fHistograms[directoryName.c_str()]->Add(eVsTheta);
    TH2F* eVsZ = new TH2F("eVsZ", "recoil energy vs. z", 1000, -100., 100., 1000, 0, 25000); fHistograms[directoryName.c_str()]->Add(eVsZ);

    TH2F* eVsZSame = new TH2F("eVsZSame", "recoil energy vs. z, first and second layer both forward or both backward", 1000, -100., 100., 1000, 0, 25000); fHistograms[directoryName.c_str()]->Add(eVsZSame);
    TH2F* eVsZCross = new TH2F("eVsZCross", "recoil energy vs. z, first and second layer over cross", 1000, -100., 100., 1000, 0, 25000); fHistograms[directoryName.c_str()]->Add(eVsZCross);


    TH2F* hEbeamRecVsSim = new TH2F("hEbeamRecVsSim", "beam energy (MeV) in lab reconstructed (x) vs simulated (y)",200,-40,40,200,-40,40); fHistograms[directoryName.c_str()]->Add(hEbeamRecVsSim);

    TH1F* excEnProtonCorr = new TH1F("excEnProtonCorr", "Excitaiton Energy after E loss correction for reconstructed protons", 5000, -5000, 5000); fHistograms[directoryName.c_str()]->Add(excEnProtonCorr);
    TH1F* excEnProtonCorrEpadCut = new TH1F("excEnProtonCorrEpadCut", "Excitaiton Energy after E loss correction for reconstructed protons with Epad>0", 5000, -20000, 20000); fHistograms[directoryName.c_str()]->Add(excEnProtonCorrEpadCut);
    TH2F* excEnProtonCorrVsX = new TH2F("excEnProtonCorrVsX", "Excitation Energy vs Vertex X after E loss correction for reconstructed protons", 5000,-10,10,5000, -20000, 20000); fHistograms[directoryName.c_str()]->Add(excEnProtonCorrVsX);
    TH2F* excEnProtonCorrVsY = new TH2F("excEnProtonCorrVsY", "Excitation Energy vs Vertex Y after E loss correction for reconstructed protons", 5000,-10,10,5000, -20000, 20000); fHistograms[directoryName.c_str()]->Add(excEnProtonCorrVsY);
    TH2F* excEnProtonCorrVsZ = new TH2F("excEnProtonCorrVsZ", "Excitation Energy vs Vertex Z after E loss correction for reconstructed protons", 5000,-100,100,5000, -20000, 20000); fHistograms[directoryName.c_str()]->Add(excEnProtonCorrVsZ);
    TH2F* excEnProtonCorrVsT = new TH2F("excEnProtonCorrVsT", "Excitation Energy vs Vertex T after E loss correction for reconstructed protons", 5000,-10,10,5000, -20000, 20000); fHistograms[directoryName.c_str()]->Add(excEnProtonCorrVsT);
    TH2F* excEnProtonCorrVsR = new TH2F("excEnProtonCorrVsR", "Excitation Energy vs Vertex R after E loss correction for reconstructed protons", 4400,-10,100,5000, -20000, 20000); fHistograms[directoryName.c_str()]->Add(excEnProtonCorrVsR);
    TH1F* excEnProtonCorrdE1Sigma1 = new TH1F("excEnProtonCorrdE1Sigma1", "Excitation Energy for dE1 1#sigma range after Eloss correction for protons", 5000, -5000, 5000); fHistograms[directoryName.c_str()]->Add(excEnProtonCorrdE1Sigma1);
    TH1F* excEnProtonCorrdE1Sigma2 = new TH1F("excEnProtonCorrdE1Sigma2", "Excitation Energy for dE1 2#sigma range after Eloss correction for protons", 5000, -5000, 5000); fHistograms[directoryName.c_str()]->Add(excEnProtonCorrdE1Sigma2);

    TH1F* hdE1ElossRange = new TH1F("hdE1ElossRange", "corrected energy range in the first layer", 1000, -100., 3900.); fHistograms[directoryName.c_str()]->Add(hdE1ElossRange);
    TH1F* hdE2ElossRange = new TH1F("hdE2ElossRange", "corrected energy range in the second layer", 2000, -1000., 9000.); fHistograms[directoryName.c_str()]->Add(hdE2ElossRange);
    TH1F* hdE2ElossRangeWoEpad0 = new TH1F("hdE2ElossRangeWoEpad0", "corrected energy range in the second layer with Epad>0", 2000, -1000., 9000.); fHistograms[directoryName.c_str()]->Add(hdE2ElossRangeWoEpad0);
    TH1F* hdE1Eloss = new TH1F("hdE1Eloss", "corrected energy loss in the first layer", 1000, -100., 9900.); fHistograms[directoryName.c_str()]->Add(hdE1Eloss);
    TH1F* hdE2Eloss = new TH1F("hdE2Eloss", "corrected energy loss in the second layer", 2000, -100., 18900.); fHistograms[directoryName.c_str()]->Add(hdE2Eloss);
    TH1F* hdE1Measured = new TH1F("hdE1Measured", "measured energy loss in the first layer", 1000, -100., 9900.); fHistograms[directoryName.c_str()]->Add(hdE1Measured);
    TH1F* hdE2Measured = new TH1F("hdE2Measured", "measured energy loss in the second layer", 2000, -100., 18900.); fHistograms[directoryName.c_str()]->Add(hdE2Measured);
    TH1F* hErestMeasured = new TH1F("hErestMeasured", "measured energy loss in the pad", 2000, -1000., 39000.); fHistograms[directoryName.c_str()]->Add(hErestMeasured);
    TH2F* hdE1ElossVsMeasured = new TH2F("hdE1ElossVsMeasured", "corrected energy loss in the first layer vs measured", 1000, -100., 2400., 1000,-100.,2400.); fHistograms[directoryName.c_str()]->Add(hdE1ElossVsMeasured);
    TH2F* hdE2ElossVsMeasuredWoEpad0 = new TH2F("hdE2ElossVsMeasuredWoEpad0", "corrected energy loss in the second layer vs measured for Epad>0", 2000, -100., 5400., 2000,-100.,5400.); fHistograms[directoryName.c_str()]->Add(hdE2ElossVsMeasuredWoEpad0);
    TH2F* hdE2ElossVsMeasured = new TH2F("hdE2ElossVsMeasured", "corrected energy loss in the second layer vs measured", 2000, -100., 9400., 2000,-100.,9400.); fHistograms[directoryName.c_str()]->Add(hdE2ElossVsMeasured);
    TH2F* dE2VsdE2Pad = new TH2F("dE2VsdE2Pad", "energy loss in second layer vs. dE2+pad energy", 2000, 0, 50000, 1000, 0, 10000); fHistograms[directoryName.c_str()]->Add(dE2VsdE2Pad);
    TH2F* EPadVsThetaLab = new TH2F("EPadVsThetaLab", "pad energy vs theta lab", 720,0,180.,2000, -1000, 49000); fHistograms[directoryName.c_str()]->Add(EPadVsThetaLab);
    TH2F* EPadVsZ = new TH2F("EPadVsZ", "pad energy vs z", 1000,-100.,100.,2000, -1000, 49000); fHistograms[directoryName.c_str()]->Add(EPadVsZ);
    TH2F* dE2VsThetaLabEpadCut = new TH2F("dE2VsThetaLabEpadCut", "second layer energy vs theta lab for Epad=0", 720,0,180.,1000, -100, 4900); fHistograms[directoryName.c_str()]->Add(dE2VsThetaLabEpadCut);
    TH2F* dE2VsEPadThetaCut = new TH2F("dE2VsEPadThetaCut", "energy loss in second layer vs. pad energy for Theta=40", 2000, -1000, 49000, 1000, -100, 9900); fHistograms[directoryName.c_str()]->Add(dE2VsEPadThetaCut);
    TH2F* dE1VsThetaLab = new TH2F("dE1VsThetaLab", "first layer energy vs theta lab", 720,0,180.,2000, -100, 4900); fHistograms[directoryName.c_str()]->Add(dE1VsThetaLab);
    TH2F* dE2VsThetaLab = new TH2F("dE2VsThetaLab", "second layer energy vs theta lab", 720,0,180.,2000, -100, 4900); fHistograms[directoryName.c_str()]->Add(dE2VsThetaLab);
    TH2F* dE12VsThetaLab = new TH2F("dE12VsThetaLab", "first+second layer energy vs theta lab", 720,0,180.,2000, -100, 7900); fHistograms[directoryName.c_str()]->Add(dE12VsThetaLab);
    TH2F* dE1EpadVsThetaLab = new TH2F("dE1EpadVsThetaLab", "first layer + pad energy vs theta lab", 720,0,180.,2000, -1000, 49900); fHistograms[directoryName.c_str()]->Add(dE1EpadVsThetaLab);
    TH2F* dE2EpadVsThetaLab = new TH2F("dE2EpadVsThetaLab", "second layer + pad energy vs theta lab", 720,0,180.,2000, -1000, 49900); fHistograms[directoryName.c_str()]->Add(dE2EpadVsThetaLab);
    TH2F* hdE1ElossVsMeasuredEpad0 = new TH2F("hdE1ElossVsMeasuredEpad0", "corrected energy loss in the first layer vs measured for Epad=0", 1000, -100., 2400., 1000,-100.,2400.); fHistograms[directoryName.c_str()]->Add(hdE1ElossVsMeasuredEpad0);
    TH2F* hdE1ElossVsMeasuredEpadWo0 = new TH2F("hdE1ElossVsMeasuredEpadWo0", "corrected energy loss in the first layer vs measured for Epad>0", 1000, -100., 2400., 1000,-100.,2400.); fHistograms[directoryName.c_str()]->Add(hdE1ElossVsMeasuredEpadWo0);
    TH1F* hdE1MeasMinRec = new TH1F("hdE1MeasMinRec", "(measured - reconstructed) energy loss in the first layer", 1000, -500., 500.); fHistograms[directoryName.c_str()]->Add(hdE1MeasMinRec);
    TH1F* hdE2MeasMinRec = new TH1F("hdE2MeasMinRec", "(measured - reconstructed) energy loss in the second layer", 1000, -500., 500.); fHistograms[directoryName.c_str()]->Add(hdE2MeasMinRec);

    TH2F* eRecErrVsESim = new TH2F("eRecErrVsESim", "error of reconstructed energy vs. simulated energy of recoil", 1000, 0., 25000., 1000, -5000., 5000.); fHistograms[directoryName.c_str()]->Add(eRecErrVsESim);
    TH2F* thetaErrorVsZ = new TH2F("thetaErrorVsZ", "Error in #vartheta_{lab} reconstruction vs. simulated z-position;z [mm];#Delta#vartheta_{lab} [^{o}]", 200, -100, 100, 100, -15, 15); fHistograms[directoryName.c_str()]->Add(thetaErrorVsZ);
    TH2F* thetaErrorVsTheta = new TH2F("thetaErrorVsTheta", "Error in #vartheta_{lab} reconstruction vs. simulated #vartheta_{lab};#vartheta_{lab} [^{o}];#Delta#vartheta_{lab} [^{o}]", 180, 0, 180, 100, -15, 15); fHistograms[directoryName.c_str()]->Add(thetaErrorVsTheta);
    TH2F* thetaErrorVsThetaEpadCut = new TH2F("thetaErrorVsThetaEpadCut", "Error in #vartheta_{lab} reconstruction vs. simulated #vartheta_{lab};#vartheta_{lab} [^{o}];#Delta#vartheta_{lab} [^{o}] with Epad>0", 180, 0, 180, 100, -15, 15); fHistograms[directoryName.c_str()]->Add(thetaErrorVsThetaEpadCut);
    TH2F* zReactionEnergy = new TH2F("zReactionEnergy", "z position of reaction vs. Beam energy (rec.)", 200, -100, 100, 1000, 0, 1.1*sett->GetBeamEnergy()); fHistograms[directoryName.c_str()]->Add(zReactionEnergy);
    TH2F* excEnProtonVsTheta = new TH2F("excEnProtonVsTheta", "Excitation Energy Spectrum from reconstructed Protons;#vartheta_{lab}[^{o}];E_{exc} [keV]", 180, 0., 180., 5000, -20000, 20000); fHistograms[directoryName.c_str()]->Add(excEnProtonVsTheta);
    TH2F* excEnProtonVsPhi = new TH2F("excEnProtonVsPhi", "Excitation Energy Spectrum from reconstructed Protons;#varphi_{lab}[^{o}];E_{exc} [keV]", 360, -180., 180., 5000, -20000, 20000); fHistograms[directoryName.c_str()]->Add(excEnProtonVsPhi);
    TH2F* excEnProtonVsZ = new TH2F("excEnProtonVsZ", "Excitation Energy Spectrum from reconstructed Protons;z [mm];E_{exc} [keV]", 200, -100., 100., 5000, -20000, 20000); fHistograms[directoryName.c_str()]->Add(excEnProtonVsZ);
    TH2F* excEnProtonVsThetaGS = new TH2F("excEnProtonVsThetaGS", "Excitation Energy Spectrum from reconstructed Protons, ground state only;#vartheta_{lab}[^{o}];E_{exc} [keV]", 180, 0., 180., 5000, -20000, 20000); fHistograms[directoryName.c_str()]->Add(excEnProtonVsThetaGS);
    TH2F* excEnProtonVsZGS = new TH2F("excEnProtonVsZGS", "Excitation Energy Spectrum from reconstructed Protons, ground state only;z [mm];E_{exc} [keV]", 200, -100., 100., 5000, -20000, 20000); fHistograms[directoryName.c_str()]->Add(excEnProtonVsZGS);
    TH2F* excEnProtonVsThetaCm = new TH2F("excEnProtonVsThetaCm", "Excitation Energy Spectrum from reconstructed Protons;#vartheta_{cm}[^{o}];E_{exc} [keV]", 180, 0., 180., 5000, -20000, 20000); fHistograms[directoryName.c_str()]->Add(excEnProtonVsThetaCm);
    TH2F* thetaVsZ = new TH2F("thetaVsZ","#vartheta_{lab} vs. z", 200, -100., 100., 180, 0., 180.); fHistograms[directoryName.c_str()]->Add(thetaVsZ);
    TH2F* thetaVsZSame = new TH2F("thetaVsZSame","#vartheta_{lab} vs. z, first and second layer both forward or both backward", 200, -100., 100., 180, 0., 180.); fHistograms[directoryName.c_str()]->Add(thetaVsZSame);

    TH2F* thetaVsZCross = new TH2F("thetaVsZCross","#vartheta_{lab} vs. z, first and second layer over cross", 200, -100., 100., 180, 0., 180.); fHistograms[directoryName.c_str()]->Add(thetaVsZCross);
    TH2F* phiVsZ = new TH2F("phiVsZ","#varphi_{lab} vs. z", 200, -100., 100., 360, -180., 180.); fHistograms[directoryName.c_str()]->Add(phiVsZ);
    TH2F* hitpattern = new TH2F("hitpattern","detector # of second layer vs. detector # of first layer", 2, -0.5, 1.5, 2, -0.5, 1.5); fHistograms[directoryName.c_str()]->Add(hitpattern);
    TH2F* eCmVsZ = new TH2F("eCmVsZ","energy of cm-system vs. z;z [mm];e-cm [GeV]", 200, -100., 100., 2000, transferP->GetCmEnergy(0.)/1000., transferP->GetCmEnergy(sett->GetBeamEnergy())/1000.); fHistograms[directoryName.c_str()]->Add(eCmVsZ);
    TH2F* betaCmVsZ = new TH2F("betaCmVsZ","#beta of cm-system vs. z", 200, -100., 100., 2000, 0., 0.2); fHistograms[directoryName.c_str()]->Add(betaCmVsZ);
    TH2F* stripPattern = new TH2F("stripPattern","Parallel strip # (#varphi) vs. perpendicular strip # (#vartheta)", 2*fSettings->GetTISTARnStripsY(0), 0., 2.*fSettings->GetTISTARnStripsY(0), 4*fSettings->GetTISTARnStripsZ(0), 0., 4.*fSettings->GetTISTARnStripsZ(0)); fHistograms[directoryName.c_str()]->Add(stripPattern);
    TH2F* recBeamEnergyErrVsZ = new TH2F("recBeamEnergyErrVsZ","Error in reconstructed beam energy vs. z", 200, -100., 100., 1000, -50., 50.); fHistograms[directoryName.c_str()]->Add(recBeamEnergyErrVsZ);
    TH2F* thetaCmVsThetaLab = new TH2F("thetaCmVsThetaLab", "#vartheta_{cm} vs. #vartheta_{lab};#vartheta_{cm} [^{o}];#vartheta_{lab} [^{o}]", 180,0.,180., 180,0.,180.); fHistograms[directoryName.c_str()]->Add(thetaCmVsThetaLab);
    TH2F* zErrorVsThetaSim = new TH2F("zErrorVsThetaSim", "Error between reconstructed and true z-position vs. simulated #vartheta_{lab}", 180, 0., 180., 1000, -5, 5); fHistograms[directoryName.c_str()]->Add(zErrorVsThetaSim);
    TH2F* zErrorVsThetaRec = new TH2F("zErrorVsThetaRec", "Error between reconstructed and true z-position vs. recon structed #vartheta_{lab}", 180, 0., 180., 1000, -5, 5); fHistograms[directoryName.c_str()]->Add(zErrorVsThetaRec);
    TH2F* zErrorVsthetaError = new TH2F("zErrorVstheaError", "z Erorr vs error in theta", 200, -10., 10., 1000, -5., 5.); fHistograms[directoryName.c_str()]->Add(zErrorVsthetaError);
    TH2F* elossVsTheta = new TH2F("elossVsTheta", "reconstructed energy loss vs. #vartheta;#vartheta_lab [^{o}];energy loss [keV]", 360, -180., 180., 10000, -100000., 100000.); fHistograms[directoryName.c_str()]->Add(elossVsTheta);
    TH2F* elossVsPhi = new TH2F("elossVsPhi", "reconstructed energy loss vs. #varphi;#varphi_lab [^{o}];energy loss [keV]", 360, -180., 180., 1000, -10000., 10000.); fHistograms[directoryName.c_str()]->Add(elossVsPhi);
    TH2F* excEnElossVsTheta = new TH2F("excEnElossVsTheta", "Excitation Energy Spectrum from reconstructed energy loss;#vartheta_{lab}[^{o}];E_{exc} [keV]", 360, -180., 180., 5000,-200000, 200000); fHistograms[directoryName.c_str()]->Add(excEnElossVsTheta);
    TH2F* excEnElossVsThetaEpadCut = new TH2F("excEnElossVsThetaEpadCut", "Excitation Energy Spectrum from reconstructed energy loss;#vartheta_{lab}[^{o}];E_{exc} [keV] with Epad>0", 180, 0., 180., 5000, -20000, 20000); fHistograms[directoryName.c_str()]->Add(excEnElossVsThetaEpadCut);

    UInt_t nofLevels = fTISTARGenChain.GetMaximum("reaction")+1;
    if(nofLevels < 1) nofLevels = 1;
    if(nofLevels > 10) nofLevels = 10;
    std::vector<TH2F*> excEnElossVsThetaLevel(nofLevels);
    for(size_t r = 0; r < nofLevels; ++r) {
        excEnElossVsThetaLevel[r] = new TH2F(Form("excEnElossVsThetaLevel_%d",static_cast<int>(r)), Form("Excitation Energy Spectrum from reconstructed energy loss for level %d;#vartheta_{lab}[^{o}];E_{exc} [keV]",static_cast<int>(r)), 180, 0., 180., 5000, -20000, 20000); fHistograms[directoryName.c_str()]->Add(excEnElossVsThetaLevel[r]);
    }

    // particle-gamma matrices
    TH1F* gammaSpec = new TH1F("gammaSpec", "generated gamma-ray spectrum", 5000, 0, 5000); fHistograms[directoryName.c_str()]->Add(gammaSpec);
    TH1F* gammaSpecDoppCorr = new TH1F("gammaSpecDoppCorr", "generated gamma-ray spectrum with doppler correction", 5000, 0, 5000); fHistograms[directoryName.c_str()]->Add(gammaSpecDoppCorr);
    TH1F* gammaSpecDoppCorrRes = new TH1F("gammaSpecDoppCorrRes", "generated gamma-ray spectrum with doppler correction w/ 1% resolution applied", 5000, 0, 5000); fHistograms[directoryName.c_str()]->Add(gammaSpecDoppCorrRes);

    TH2F* excEnProtonVsGamma = new TH2F("excEnProtonVsGamma", "Excitation Energy Spectrum from reconstructed Protons vs gamma ray energy", 2500, 0, 5000, 5000, -5000, 5000); fHistograms[directoryName.c_str()]->Add(excEnProtonVsGamma);
    TH2F* excEnProtonVsGammaDoppCorr = new TH2F("excEnProtonVsGammaDoppCorr", "Excitation Energy Spectrum from reconstructed Protons vs gamma ray energy w/ doppler corrections", 2500, 0, 5000, 5000, -5000, 5000); fHistograms[directoryName.c_str()]->Add(excEnProtonVsGammaDoppCorr);
    TH2F* excEnProtonVsGammaDoppCorrRes = new TH2F("excEnProtonVsGammaDoppCorrRes", "Excitation Energy Spectrum from reconstructed Protons vs gamma ray energy w/ doppler corrections and 1% resolution applied", 2500, 0, 5000, 5000, -5000, 5000); fHistograms[directoryName.c_str()]->Add(excEnProtonVsGammaDoppCorrRes);


}


