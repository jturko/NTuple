#include "Settings.hh"

#include "TEnv.h"
#include "TString.h"

Settings::Settings(std::string fileName, int verbosityLevel) : 
    fVerbosityLevel(verbosityLevel)
    {
    TEnv env;
    env.ReadFile(fileName.c_str(),kEnvLocal);

    //  env.PrintEnv();

    fNtupleName = env.GetValue("NtupleName","/ntuple");

    fBufferSize = env.GetValue("BufferSize",1024000);

    fSortNumberOfEvents = env.GetValue("SortNumberOfEvents",0);

    fWriteTree = env.GetValue("WriteTree",true);

    fWrite2DHist = env.GetValue("Write2DHist",false);

    //fWrite3DHist = env.GetValue("Write3DHist",false);
        
    fWriteNDHist = env.GetValue("WriteNDHist",false);

    fWrite2DSGGHist = env.GetValue("Write2DSGGHist",false);

    fWriteGriffinAddbackVector = env.GetValue("WriteGriffinAddbackVector",false);

    fGriffinAddbackVectorLengthmm = env.GetValue("GriffinAddbackVectorLengthmm",105.0);

    fGriffinAddbackVectorDepthmm = env.GetValue("GriffinAddbackVectorDepthmm",45.0);

    fGriffinAddbackVectorCrystalFaceDistancemm = env.GetValue("GriffinAddbackVectorCrystalFaceDistancemm",110.0);

    // TI-STAR detector/run variables
    // assuming 2 pixelated strips
    fTISTARGenNtupleName =      env.GetValue("TISTAR.NtupleName","/treeGen");
    fTISTARStripWidthX.resize(2); 
    fTISTARStripWidthY.resize(2); 
    fTISTARnStripsX.resize(2);
    fTISTARnStripsY.resize(2);
    fTISTARStripWidthX[0] =  env.GetValue("TISTAR.Layer0.StripWidthX", 0.1); // in mm
    fTISTARStripWidthY[0] =  env.GetValue("TISTAR.Layer0.StripWidthY", 0.1); // in mm
    fTISTARStripWidthX[1] =  env.GetValue("TISTAR.Layer1.StripWidthX", 0.1); // in mm
    fTISTARStripWidthY[1] =  env.GetValue("TISTAR.Layer1.StripWidthY", 0.1); // in mm

    // Griffin
    fResolution[1000].resize(16);
    fThreshold[1000].resize(16,std::vector<double>(4));
    fThresholdWidth[1000].resize(16,std::vector<double>(4));
    fTimeWindow[1000].resize(16,std::vector<double>(4));

    fResolution[1010].resize(16);
    fThreshold[1010].resize(16,std::vector<double>(4));
    fThresholdWidth[1010].resize(16,std::vector<double>(4));
    fTimeWindow[1010].resize(16,std::vector<double>(4));

    fResolution[1020].resize(16);
    fThreshold[1020].resize(16,std::vector<double>(4));
    fThresholdWidth[1020].resize(16,std::vector<double>(4));
    fTimeWindow[1020].resize(16,std::vector<double>(4));

    fResolution[1030].resize(16);
    fThreshold[1030].resize(16,std::vector<double>(4));
    fThresholdWidth[1030].resize(16,std::vector<double>(4));
    fTimeWindow[1030].resize(16,std::vector<double>(4));

    fResolution[1040].resize(16);
    fThreshold[1040].resize(16,std::vector<double>(4));
    fThresholdWidth[1040].resize(16,std::vector<double>(4));
    fTimeWindow[1040].resize(16,std::vector<double>(4));

    fResolution[1050].resize(16);
    fThreshold[1050].resize(16,std::vector<double>(4));
    fThresholdWidth[1050].resize(16,std::vector<double>(4));
    fTimeWindow[1050].resize(16,std::vector<double>(4));

    // LaBr3
    fResolution[2000].resize(16);
    fThreshold[2000].resize(16,std::vector<double>(1));
    fThresholdWidth[2000].resize(16,std::vector<double>(1));
    fTimeWindow[2000].resize(16,std::vector<double>(1));

    // Sceptar
    fResolution[5000].resize(20);
    fThreshold[5000].resize(20,std::vector<double>(1));
    fThresholdWidth[5000].resize(20,std::vector<double>(1));
    fTimeWindow[5000].resize(20,std::vector<double>(1));

    // EightPi
    fResolution[6000].resize(20);
    fThreshold[6000].resize(20,std::vector<double>(4));
    fThresholdWidth[6000].resize(20,std::vector<double>(4));
    fTimeWindow[6000].resize(20,std::vector<double>(4));

    fResolution[6010].resize(20);
    fThreshold[6010].resize(20,std::vector<double>(4));
    fThresholdWidth[6010].resize(20,std::vector<double>(4));
    fTimeWindow[6010].resize(20,std::vector<double>(4));

    fResolution[6020].resize(20);
    fThreshold[6020].resize(20,std::vector<double>(4));
    fThresholdWidth[6020].resize(20,std::vector<double>(4));
    fTimeWindow[6020].resize(20,std::vector<double>(4));

    fResolution[6030].resize(20);
    fThreshold[6030].resize(20,std::vector<double>(4));
    fThresholdWidth[6030].resize(20,std::vector<double>(4));
    fTimeWindow[6030].resize(20,std::vector<double>(4));

    // Descant
    fResolution[8010].resize(15);
    fThreshold[8010].resize(15,std::vector<double>(1));
    fThresholdWidth[8010].resize(15,std::vector<double>(1));
    fTimeWindow[8010].resize(15,std::vector<double>(1));

    fResolution[8020].resize(10);
    fThreshold[8020].resize(10,std::vector<double>(1));
    fThresholdWidth[8020].resize(10,std::vector<double>(1));
    fTimeWindow[8020].resize(10,std::vector<double>(1));

    fResolution[8030].resize(15);
    fThreshold[8030].resize(15,std::vector<double>(1));
    fThresholdWidth[8030].resize(15,std::vector<double>(1));
    fTimeWindow[8030].resize(15,std::vector<double>(1));

    fResolution[8040].resize(20);
    fThreshold[8040].resize(20,std::vector<double>(1));
    fThresholdWidth[8040].resize(20,std::vector<double>(1));
    fTimeWindow[8040].resize(20,std::vector<double>(1));

    fResolution[8050].resize(10);
    fThreshold[8050].resize(10,std::vector<double>(1));
    fThresholdWidth[8050].resize(10,std::vector<double>(1));
    fTimeWindow[8050].resize(10,std::vector<double>(1));

    // Testcan
    fResolution[8500].resize(1);
    fThreshold[8500].resize(1,std::vector<double>(1));
    fThresholdWidth[8500].resize(1,std::vector<double>(1));
    fTimeWindow[8500].resize(1,std::vector<double>(1));
    fProtonCoeff.resize(4);
    fDeuteronCoeff.resize(4);
    fCarbonCoeff.resize(4);
    fBeCoeff.resize(4);
    fBCoeff.resize(4);
    fAlphaCoeff.resize(4);

    // Paces
    fResolution[9000].resize(5);
    fThreshold[9000].resize(5,std::vector<double>(1));
    fThresholdWidth[9000].resize(5,std::vector<double>(1));
    fTimeWindow[9000].resize(5,std::vector<double>(1));

    // TI-STAR
    // resizing to 3 for the each of the layers, the each subvector to 4 for a maximum of 4 strips per layer
    fResolution[9500].resize(3);
    fThreshold[9500].resize(3,std::vector<double>(4));
    fThresholdWidth[9500].resize(3,std::vector<double>(4));
    fTimeWindow[9500].resize(3,std::vector<double>(4));

    double offset, linear, quadratic, cubic;
    double A, B, C;
    //A=0.0826;B=0.2673;C=0.0548; 
    //A=0.15;B=0.1;C=0.02; 
    // Griffin
    for(int detector = 0; detector < 16; ++detector) {
        for(int crystal = 0; crystal < 4; ++crystal) {
            offset = env.GetValue(Form("Griffin.%d.%d.Resolution.Offset",detector,crystal),1.100);
            linear = env.GetValue(Form("Griffin.%d.%d.Resolution.Linear",detector,crystal),0.00183744);
            quadratic = env.GetValue(Form("Griffin.%d.%d.Resolution.Quadratic",detector,crystal),0.0000007);
            cubic = env.GetValue(Form("Griffin.%d.%d.Resolution.Cubic",detector,crystal),0.);
            fResolution[1000][detector].push_back(TF1(Form("Griffin.%d.%d.Resolution",detector,crystal),
                                                      Form("(TMath::Sqrt((%f+%f*x+%f*x*x+%f*x*x*x)))/(2.*TMath::Sqrt(2.*TMath::Log(2.)))",offset, linear, quadratic, cubic),0.,100000.));
            fThreshold[1000][detector][crystal] = env.GetValue(Form("Griffin.%d.%d.Threshold.keV",detector,crystal),10.);
            fThresholdWidth[1000][detector][crystal] = env.GetValue(Form("Griffin.%d.%d.ThresholdWidth.keV",detector,crystal),2.);
            fTimeWindow[1000][detector][crystal] = env.GetValue(Form("Griffin.%d.%d.TimeWindow.sec",detector,crystal),0.);

            offset = env.GetValue(Form("Griffin.BGO.Front.Left.%d.%d.Resolution.Offset",detector,crystal),1.100);
            linear = env.GetValue(Form("Griffin.BGO.Front.Left.%d.%d.Resolution.Linear",detector,crystal),0.00183744);
            quadratic = env.GetValue(Form("Griffin.BGO.Front.Left.%d.%d.Resolution.Quadratic",detector,crystal),0.0000007);
            cubic = env.GetValue(Form("Griffin.BGO.Front.Left.%d.%d.Resolution.Cubic",detector,crystal),0.);
            fResolution[1010][detector].push_back(TF1(Form("Griffin.BGO.Front.Left.%d.%d.Resolution",detector,crystal),
                                                      Form("(TMath::Sqrt((%f+%f*x+%f*x*x+%f*x*x*x)))/(2.*TMath::Sqrt(2.*TMath::Log(2.)))",offset, linear, quadratic, cubic),0.,100000.));
            fThreshold[1010][detector][crystal] = env.GetValue(Form("Griffin.BGO.Front.Left.%d.%d.Threshold.keV",detector,crystal),10.);
            fThresholdWidth[1010][detector][crystal] = env.GetValue(Form("Griffin.BGO.Front.Left.%d.%d.ThresholdWidth.keV",detector,crystal),2.);
            fTimeWindow[1010][detector][crystal] = env.GetValue(Form("Griffin.BGO.Front.Left.%d.%d.TimeWindow.sec",detector,crystal),0.);

            offset = env.GetValue(Form("Griffin.BGO.Front.Right.%d.%d.Resolution.Offset",detector,crystal),1.100);
            linear = env.GetValue(Form("Griffin.BGO.Front.Right.%d.%d.Resolution.Linear",detector,crystal),0.00183744);
            quadratic = env.GetValue(Form("Griffin.BGO.Front.Right.%d.%d.Resolution.Quadratic",detector,crystal),0.0000007);
            cubic = env.GetValue(Form("Griffin.BGO.Front.Right.%d.%d.Resolution.Cubic",detector,crystal),0.);
            fResolution[1020][detector].push_back(TF1(Form("Griffin.BGO.Front.Right.%d.%d.Resolution",detector,crystal),
                                                      Form("(TMath::Sqrt((%f+%f*x+%f*x*x+%f*x*x*x)))/(2.*TMath::Sqrt(2.*TMath::Log(2.)))",offset, linear, quadratic, cubic),0.,100000.));
            fThreshold[1020][detector][crystal] = env.GetValue(Form("Griffin.BGO.Front.Right.%d.%d.Threshold.keV",detector,crystal),10.);
            fThresholdWidth[1020][detector][crystal] = env.GetValue(Form("Griffin.BGO.Front.Right.%d.%d.ThresholdWidth.keV",detector,crystal),2.);
            fTimeWindow[1020][detector][crystal] = env.GetValue(Form("Griffin.BGO.Front.Right.%d.%d.TimeWindow.sec",detector,crystal),0.);

            offset = env.GetValue(Form("Griffin.BGO.Side.Left.%d.%d.Resolution.Offset",detector,crystal),1.100);
            linear = env.GetValue(Form("Griffin.BGO.Side.Left.%d.%d.Resolution.Linear",detector,crystal),0.00183744);
            quadratic = env.GetValue(Form("Griffin.BGO.Side.Left.%d.%d.Resolution.Quadratic",detector,crystal),0.0000007);
            cubic = env.GetValue(Form("Griffin.BGO.Side.Left.%d.%d.Resolution.Cubic",detector,crystal),0.);
            fResolution[1030][detector].push_back(TF1(Form("Griffin.BGO.Side.Left.%d.%d.Resolution",detector,crystal),
                                                      Form("(TMath::Sqrt((%f+%f*x+%f*x*x+%f*x*x*x)))/(2.*TMath::Sqrt(2.*TMath::Log(2.)))",offset, linear, quadratic, cubic),0.,100000.));
            fThreshold[1030][detector][crystal] = env.GetValue(Form("Griffin.BGO.Side.Left.%d.%d.Threshold.keV",detector,crystal),10.);
            fThresholdWidth[1030][detector][crystal] = env.GetValue(Form("Griffin.BGO.Side.Left.%d.%d.ThresholdWidth.keV",detector,crystal),2.);
            fTimeWindow[1030][detector][crystal] = env.GetValue(Form("Griffin.BGO.Side.Left.%d.%d.TimeWindow.sec",detector,crystal),0.);

            offset = env.GetValue(Form("Griffin.BGO.Side.Right.%d.%d.Resolution.Offset",detector,crystal),1.100);
            linear = env.GetValue(Form("Griffin.BGO.Side.Right.%d.%d.Resolution.Linear",detector,crystal),0.00183744);
            quadratic = env.GetValue(Form("Griffin.BGO.Side.Right.%d.%d.Resolution.Quadratic",detector,crystal),0.0000007);
            cubic = env.GetValue(Form("Griffin.BGO.Side.Right.%d.%d.Resolution.Cubic",detector,crystal),0.);
            fResolution[1040][detector].push_back(TF1(Form("Griffin.BGO.Side.Right.%d.%d.Resolution",detector,crystal),
                                                      Form("(TMath::Sqrt((%f+%f*x+%f*x*x+%f*x*x*x)))/(2.*TMath::Sqrt(2.*TMath::Log(2.)))",offset, linear, quadratic, cubic),0.,100000.));
            fThreshold[1040][detector][crystal] = env.GetValue(Form("Griffin.BGO.Side.Right.%d.%d.Threshold.keV",detector,crystal),10.);
            fThresholdWidth[1040][detector][crystal] = env.GetValue(Form("Griffin.BGO.Side.Right.%d.%d.ThresholdWidth.keV",detector,crystal),2.);
            fTimeWindow[1040][detector][crystal] = env.GetValue(Form("Griffin.BGO.Side.Right.%d.%d.TimeWindow.sec",detector,crystal),0.);

            offset = env.GetValue(Form("Griffin.BGO.Back.%d.%d.Resolution.Offset",detector,crystal),1.100);
            linear = env.GetValue(Form("Griffin.BGO.Back.%d.%d.Resolution.Linear",detector,crystal),0.00183744);
            quadratic = env.GetValue(Form("Griffin.BGO.Back.%d.%d.Resolution.Quadratic",detector,crystal),0.0000007);
            cubic = env.GetValue(Form("Griffin.BGO.Back.%d.%d.Resolution.Cubic",detector,crystal),0.);
            fResolution[1050][detector].push_back(TF1(Form("Griffin.BGO.Back.%d.%d.Resolution",detector,crystal),
                                                      Form("(TMath::Sqrt((%f+%f*x+%f*x*x+%f*x*x*x)))/(2.*TMath::Sqrt(2.*TMath::Log(2.)))",offset, linear, quadratic, cubic),0.,100000.));
            fThreshold[1050][detector][crystal] = env.GetValue(Form("Griffin.BGO.Back.%d.%d.Threshold.keV",detector,crystal),10.);
            fThresholdWidth[1050][detector][crystal] = env.GetValue(Form("Griffin.BGO.Back.%d.%d.ThresholdWidth.keV",detector,crystal),2.);
            fTimeWindow[1050][detector][crystal] = env.GetValue(Form("Griffin.BGO.Back.%d.%d.TimeWindow.sec",detector,crystal),0.);
        }
    }

    // LaBr3
    for(int detector = 0; detector < 16; ++detector) {
        offset = env.GetValue(Form("LaBr3.%d.Resolution.Offset",detector),1.7006116);
        linear = env.GetValue(Form("LaBr3.%d.Resolution.Linear",detector),0.5009382);
        quadratic = env.GetValue(Form("LaBr3.%d.Resolution.Quadratic",detector),0.000065451219);
        cubic = env.GetValue(Form("LaBr3.%d.Resolution.Cubic",detector),0.);
        fResolution[2000][detector].push_back(TF1(Form("LaBr3.%d.Resolution",detector),
                                                  Form("(TMath::Sqrt((%f+%f*x+%f*x*x+%f*x*x*x)))/(2.*TMath::Sqrt(2.*TMath::Log(2.)))",offset, linear, quadratic, cubic),0.,100000.));
        fThreshold[2000][detector][0] = env.GetValue(Form("LaBr3.%d.Threshold.keV",detector),10.);
        fThresholdWidth[2000][detector][0] = env.GetValue(Form("LaBr3.%d.ThresholdWidth.keV",detector),2.);
        fTimeWindow[2000][detector][0] = env.GetValue(Form("LaBr3.%d.TimeWindow.sec",detector),0.);
    }

    // Sceptar
    for(int detector = 0; detector < 20; ++detector) {
        offset = env.GetValue(Form("Sceptar.%d.Resolution.Offset",detector),0.0);
        linear = env.GetValue(Form("Sceptar.%d.Resolution.Linear",detector),0.0);
        quadratic = env.GetValue(Form("Sceptar.%d.Resolution.Quadratic",detector),0.0);
        cubic = env.GetValue(Form("Sceptar.%d.Resolution.Cubic",detector),0.0);
        fResolution[5000][detector].push_back(TF1(Form("Sceptar.%d.Resolution",detector),
                                                  Form("(TMath::Sqrt((%f+%f*x+%f*x*x+%f*x*x*x)))/(2.*TMath::Sqrt(2.*TMath::Log(2.)))",offset, linear, quadratic, cubic),0.,100000.));
        fThreshold[5000][detector][0] = env.GetValue(Form("Sceptar.%d.Threshold.keV",detector),0.0);
        fThresholdWidth[5000][detector][0] = env.GetValue(Form("Sceptar.%d.ThresholdWidth.keV",detector),0.0);
        fTimeWindow[5000][detector][0] = env.GetValue(Form("Sceptar.%d.TimeWindow.sec",detector),0.0);
    }

    // EightPi
    for(int detector = 0; detector < 20; ++detector) {
        offset = env.GetValue(Form("EightPi.%d.Resolution.Offset",detector),1.100);
        linear = env.GetValue(Form("EightPi.%d.Resolution.Linear",detector),0.00183744);
        quadratic = env.GetValue(Form("EightPi.%d.Resolution.Quadratic",detector),0.0000007);
        cubic = env.GetValue(Form("EightPi.%d.Resolution.Cubic",detector),0.);
        fResolution[6000][detector].push_back(TF1(Form("EightPi.%d.Resolution",detector),
                                                  Form("(TMath::Sqrt((%f+%f*x+%f*x*x+%f*x*x*x)))/(2.*TMath::Sqrt(2.*TMath::Log(2.)))",offset, linear, quadratic, cubic),0.,100000.));
        fThreshold[6000][detector][0] = env.GetValue(Form("EightPi.%d.Threshold.keV",detector),10.);
        fThresholdWidth[6000][detector][0] = env.GetValue(Form("EightPi.%d.ThresholdWidth.keV",detector),2.);
        fTimeWindow[6000][detector][0] = env.GetValue(Form("EightPi.%d.TimeWindow.sec",detector),0.);

        offset = env.GetValue(Form("EightPi.BGO.%d.Resolution.Offset",detector),1.100);
        linear = env.GetValue(Form("EightPi.BGO.%d.Resolution.Linear",detector),0.00183744);
        quadratic = env.GetValue(Form("EightPi.BGO.%d.Resolution.Quadratic",detector),0.0000007);
        cubic = env.GetValue(Form("EightPi.BGO.%d.Resolution.Cubic",detector),0.);
        fResolution[6010][detector].push_back(TF1(Form("EightPi.BGO.%d.Resolution",detector),
                                                  Form("(TMath::Sqrt((%f+%f*x+%f*x*x+%f*x*x*x)))/(2.*TMath::Sqrt(2.*TMath::Log(2.)))",offset, linear, quadratic, cubic),0.,100000.));
        fThreshold[6010][detector][0] = env.GetValue(Form("EightPi.BGO.%d.Threshold.keV",detector),10.);
        fThresholdWidth[6010][detector][0] = env.GetValue(Form("EightPi.BGO.%d.ThresholdWidth.keV",detector),2.);
        fTimeWindow[6010][detector][0] = env.GetValue(Form("EightPi.BGO.%d.TimeWindow.sec",detector),0.);

        offset = env.GetValue(Form("EightPi.BGO.%d.Resolution.Offset",detector),1.100);
        linear = env.GetValue(Form("EightPi.BGO.%d.Resolution.Linear",detector),0.00183744);
        quadratic = env.GetValue(Form("EightPi.BGO.%d.Resolution.Quadratic",detector),0.0000007);
        cubic = env.GetValue(Form("EightPi.BGO.%d.Resolution.Cubic",detector),0.);
        fResolution[6020][detector].push_back(TF1(Form("EightPi.BGO.%d.Resolution",detector),
                                                  Form("(TMath::Sqrt((%f+%f*x+%f*x*x+%f*x*x*x)))/(2.*TMath::Sqrt(2.*TMath::Log(2.)))",offset, linear, quadratic, cubic),0.,100000.));
        fThreshold[6020][detector][0] = env.GetValue(Form("EightPi.BGO.%d.Threshold.keV",detector),10.);
        fThresholdWidth[6020][detector][0] = env.GetValue(Form("EightPi.BGO.%d.ThresholdWidth.keV",detector),2.);
        fTimeWindow[6020][detector][0] = env.GetValue(Form("EightPi.BGO.%d.TimeWindow.sec",detector),0.);

        offset = env.GetValue(Form("EightPi.BGO.%d.Resolution.Offset",detector),1.100);
        linear = env.GetValue(Form("EightPi.BGO.%d.Resolution.Linear",detector),0.00183744);
        quadratic = env.GetValue(Form("EightPi.BGO.%d.Resolution.Quadratic",detector),0.0000007);
        cubic = env.GetValue(Form("EightPi.BGO.%d.Resolution.Cubic",detector),0.);
        fResolution[6030][detector].push_back(TF1(Form("EightPi.BGO.%d.Resolution",detector),
                                                  Form("(TMath::Sqrt((%f+%f*x+%f*x*x+%f*x*x*x)))/(2.*TMath::Sqrt(2.*TMath::Log(2.)))",offset, linear, quadratic, cubic),0.,100000.));
        fThreshold[6030][detector][0] = env.GetValue(Form("EightPi.BGO.%d.Threshold.keV",detector),10.);
        fThresholdWidth[6030][detector][0] = env.GetValue(Form("EightPi.BGO.%d.ThresholdWidth.keV",detector),2.);
        fTimeWindow[6030][detector][0] = env.GetValue(Form("EightPi.BGO.%d.TimeWindow.sec",detector),0.);
    }

    // DESCANT
    for(int detector = 0; detector < 15; ++detector) {
        offset = env.GetValue(Form("Descant.Blue.%d.Resolution.Offset",detector),0.0);
        linear = env.GetValue(Form("Descant.Blue.%d.Resolution.Linear",detector),0.0);
        quadratic = env.GetValue(Form("Descant.Blue.%d.Resolution.Quadratic",detector),0.009);
        cubic = env.GetValue(Form("Descant.Blue.%d.Resolution.Cubic",detector),0.0);
        fResolution[8010][detector].push_back(TF1(Form("Descant.Blue.%d.Resolution",detector),
                                                  Form("(TMath::Sqrt((%f+%f*x+%f*x*x+%f*x*x*x)))/(2.*TMath::Sqrt(2.*TMath::Log(2.)))",offset, linear, quadratic, cubic),0.,100000.));
        fThreshold[8010][detector][0] = env.GetValue(Form("Descant.Blue.%d.Threshold.keV",detector),0.);
        fThresholdWidth[8010][detector][0] = env.GetValue(Form("Descant.Blue.%d.ThresholdWidth.keV",detector),0.);
        fTimeWindow[8010][detector][0] = env.GetValue(Form("Descant.Blue.%d.TimeWindow.sec",detector),0.);
    }
    for(int detector = 0; detector < 10; ++detector) {
        offset = env.GetValue(Form("Descant.Green.%d.Resolution.Offset",detector),0.0);
        linear = env.GetValue(Form("Descant.Green.%d.Resolution.Linear",detector),0.0);
        quadratic = env.GetValue(Form("Descant.Green.%d.Resolution.Quadratic",detector),0.009);
        cubic = env.GetValue(Form("Descant.Green.%d.Resolution.Cubic",detector),0.0);
        fResolution[8020][detector].push_back(TF1(Form("Descant.Green.%d.Resolution",detector),
                                                  Form("(TMath::Sqrt((%f+%f*x+%f*x*x+%f*x*x*x)))/(2.*TMath::Sqrt(2.*TMath::Log(2.)))",offset, linear, quadratic, cubic),0.,100000.));
        fThreshold[8020][detector][0] = env.GetValue(Form("Descant.Green.%d.Threshold.keV",detector),0.);
        fThresholdWidth[8020][detector][0] = env.GetValue(Form("Descant.Green.%d.ThresholdWidth.keV",detector),0.);
        fTimeWindow[8020][detector][0] = env.GetValue(Form("Descant.Green.%d.TimeWindow.sec",detector),0.);
    }
    for(int detector = 0; detector < 15; ++detector) {
        offset = env.GetValue(Form("Descant.Red.%d.Resolution.Offset",detector),0.0);
        linear = env.GetValue(Form("Descant.Red.%d.Resolution.Linear",detector),0.0);
        quadratic = env.GetValue(Form("Descant.Red.%d.Resolution.Quadratic",detector),0.009);
        cubic = env.GetValue(Form("Descant.Red.%d.Resolution.Cubic",detector),0.0);
        fResolution[8030][detector].push_back(TF1(Form("Descant.Red.%d.Resolution",detector),
                                                  Form("(TMath::Sqrt((%f+%f*x+%f*x*x+%f*x*x*x)))/(2.*TMath::Sqrt(2.*TMath::Log(2.)))",offset, linear, quadratic, cubic),0.,100000.));
        fThreshold[8030][detector][0] = env.GetValue(Form("Descant.Red.%d.Threshold.keV",detector),0.);
        fThresholdWidth[8030][detector][0] = env.GetValue(Form("Descant.Red.%d.ThresholdWidth.keV",detector),0.);
        fTimeWindow[8030][detector][0] = env.GetValue(Form("Descant.Red.%d.TimeWindow.sec",detector),0.);
    }
    for(int detector = 0; detector < 20; ++detector) {
        linear = env.GetValue(Form("Descant.White.%d.Resolution",detector),20.0);
        fResolution[8040][detector].push_back(TF1(Form("Descant.White.%d.Resolution",detector),
                                                  Form("%f*TMath::Sqrt(x)/(2.*TMath::Sqrt(2.*TMath::Log(2.)))", linear),0.,100000.));
        //offset = env.GetValue(Form("Descant.White.%d.Resolution.Offset",detector),0.0);
        //linear = env.GetValue(Form("Descant.White.%d.Resolution.Linear",detector),0.0);
        //quadratic = env.GetValue(Form("Descant.White.%d.Resolution.Quadratic",detector),0.009);
        //cubic = env.GetValue(Form("Descant.White.%d.Resolution.Cubic",detector),0.0);
        //fResolution[8040][detector].push_back(TF1(Form("Descant.White.%d.Resolution",detector),
        //                                          Form("(TMath::Sqrt((%f+%f*x+%f*x*x+%f*x*x*x)))/(2.*TMath::Sqrt(2.*TMath::Log(2.)))",offset, linear, quadratic, cubic),0.,100000.));
        fThreshold[8040][detector][0] = env.GetValue(Form("Descant.White.%d.Threshold.keV",detector),0.);
        fThresholdWidth[8040][detector][0] = env.GetValue(Form("Descant.White.%d.ThresholdWidth.keV",detector),0.);
        fTimeWindow[8040][detector][0] = env.GetValue(Form("Descant.White.%d.TimeWindow.sec",detector),0.);
    }
    for(int detector = 0; detector < 10; ++detector) {
        offset = env.GetValue(Form("Descant.Yellow.%d.Resolution.Offset",detector),0.0);
        linear = env.GetValue(Form("Descant.Yellow.%d.Resolution.Linear",detector),0.0);
        quadratic = env.GetValue(Form("Descant.Yellow.%d.Resolution.Quadratic",detector),0.009);
        cubic = env.GetValue(Form("Descant.Yellow.%d.Resolution.Cubic",detector),0.0);
        fResolution[8050][detector].push_back(TF1(Form("Descant.Yellow.%d.Resolution",detector),
                                                  Form("(TMath::Sqrt((%f+%f*x+%f*x*x+%f*x*x*x)))/(2.*TMath::Sqrt(2.*TMath::Log(2.)))",offset, linear, quadratic, cubic),0.,100000.));
        fThreshold[8050][detector][0] = env.GetValue(Form("Descant.Yellow.%d.Threshold.keV",detector),0.);
        fThresholdWidth[8050][detector][0] = env.GetValue(Form("Descant.Yellow.%d.ThresholdWidth.keV",detector),0.);
        fTimeWindow[8050][detector][0] = env.GetValue(Form("Descant.Yellow.%d.TimeWindow.sec",detector),0.);
    }
        
    // Testcan light
    linear = env.GetValue("Testcan.Resolution",20.0); //std::cout <<"Testcan resolution = " << linear << std::endl;
    A = env.GetValue("Testcan.Resolution.A",0.15);     
    B = env.GetValue("Testcan.Resolution.B",0.10);     
    C = env.GetValue("Testcan.Resolution.C",0.02); 
    for(int i = 1; i <= 4; i++) {
        fProtonCoeff[i-1] = env.GetValue(Form("Testcan.Proton.%d",i), 0.0); //std::cout << "a" << i << " = " << fProtonCoeff[i-1] << " | ";
        fDeuteronCoeff[i-1] = env.GetValue(Form("Testcan.Deuteron.%d",i), 0.0);
        fCarbonCoeff[i-1] = env.GetValue(Form("Testcan.Carbon.%d",i), 0.0);
        fBeCoeff[i-1] = env.GetValue(Form("Testcan.Be.%d",i), 0.0);
        fBCoeff[i-1] = env.GetValue(Form("Testcan.B.%d",i), 0.0);
        fAlphaCoeff[i-1] = env.GetValue(Form("Testcan.Alpha.%d",i), 0.0);
    }    
    //std::cout << "A = " << A << " | B = " << B << " | C = " << C << " | carbon = " << fCarbonCoeff[0] << std::endl;
    //fResolution[8500][0].push_back(TF1("Testcan.Resolution",Form("%f*TMath::Sqrt(x)/(2.*TMath::Sqrt(2.*TMath::Log(2.)))", linear),0.,100000.));
    //fResolution[8500][0].push_back(TF1("Testcan.Resolution",Form("x*TMath::Sqrt(TMath::Power(%f,2)+TMath::Power(%f,2)/x+TMath::Power(%f/x,2))/2.*TMath::Sqrt(2.*TMath::Log(2.)",A,B,C,0.)0.,100000.));
    fResolution[8500][0].push_back(TF1("Testcan.Resolution",Form("x*TMath::Sqrt(TMath::Power(%f,2)+TMath::Power(%f,2)/x+TMath::Power(%f/x,2))",A,B,C),0.,20.)); // From N Desplan Thesis
    //fResolution[8500][0].push_back(TF1("Testcan.Resolution",Form("x*TMath::Sqrt(TMath::Power(%f,2)+TMath::Power(%f,2)/x+TMath::Power(%f/x,2))/(2.*TMath::Sqrt(2.*TMath::Log(2.)))",A,B,C),0.,20.)); // From N Desplan Thesis
    fThreshold[8500][0][0]      = env.GetValue("Testcan.Threshold.keV",0.);
    fThresholdWidth[8500][0][0] = env.GetValue("Testcan.ThresholdWidth.keV",0.);
    fTimeWindow[8500][0][0]     = env.GetValue("Testcan.TimeWindow.sec",0.);
    
    // DESCANT quenching
    fQuenching.resize(7);
    for(int isotope = 0; isotope < Settings::kMax; ++isotope) {
        fQuenching[isotope] = env.GetValue(Form("Descant.Quenching.%d",isotope),1.0);
        //std::cout << "fQuenching[" << isotope << "] = " << fQuenching[isotope] << std::endl;
    }


	 //testcan
	 double fanoFactor = env.GetValue("Testcan.Resolution.FanoFactor",20.);
	 fResolution[8500][0].push_back(TF1("Testcan.Resolution",Form("%f*TMath::Sqrt(x)",fanoFactor), 0., 100000.));
	 fThreshold[8500][0][0] = env.GetValue("Testcan.Threshold.keV",0.);
	 fThresholdWidth[8500][0][0] = env.GetValue("Testcan.ThresholdWidth.keV",0.);
	 fTimeWindow[8500][0][0] = env.GetValue("Testcan.TimeWindow.sec",0.);


    // Paces
    for(int detector = 0; detector < 5; ++detector) {
        offset = env.GetValue(Form("Paces.%d.Resolution.Offset",detector),0.0);
        linear = env.GetValue(Form("Paces.%d.Resolution.Linear",detector),0.0);
        quadratic = env.GetValue(Form("Paces.%d.Resolution.Quadratic",detector),0.0);
        cubic = env.GetValue(Form("Paces.%d.Resolution.Cubic",detector),0.0);
        fResolution[9000][detector].push_back(TF1(Form("Paces.%d.Resolution",detector),
                                                  Form("(TMath::Sqrt((%f+%f*x+%f*x*x+%f*x*x*x)))/(2.*TMath::Sqrt(2.*TMath::Log(2.)))",offset, linear, quadratic, cubic),0.,100000.));
        fThreshold[9000][detector][0] = env.GetValue(Form("Paces.%d.Threshold.keV",detector),0.0);
        fThresholdWidth[9000][detector][0] = env.GetValue(Form("Paces.%d.ThresholdWidth.keV",detector),0.0);
        fTimeWindow[9000][detector][0] = env.GetValue(Form("Paces.%d.TimeWindow.sec",detector),0.0);
    }

    // TI-STAR
    for(int detector = 0; detector < 3; detector++) {
        for(int crystal = 0; crystal < 4; crystal++) {
            offset = env.GetValue(Form("TISTAR.Layer%d.Strip%d.Resolution.Offset",detector,crystal),0.0);
            linear = env.GetValue(Form("TISTAR.Layer%d.Strip%d.Resolution.Linear",detector,crystal),0.0);
            quadratic = env.GetValue(Form("TISTAR.Layer%d.Strip%d.Resolution.Quadratic",detector,crystal),0.0);
            cubic = env.GetValue(Form("TISTAR.Layer%d.Strip%d.Resolution.Cubic",detector,crystal),0.0);
            fResolution[9500][detector].push_back(TF1(Form("TISTAR.Layer%d.Strip%d.Resolution",detector,crystal),
                                                      Form("0.0"),0.,100000.));
            fThreshold[9500][detector][crystal] = env.GetValue(Form("TISTAR.Layer%d.Strip%d.Threshold.keV",detector,crystal),10.);
            fThresholdWidth[9500][detector][crystal] = env.GetValue(Form("TISTAR.Layer%d.Strip%d.ThresholdWidth.keV",detector,crystal),2.);
            fTimeWindow[9500][detector][crystal] = env.GetValue(Form("TISTAR.Layer%d.Strip%d.TimeWindow.sec",detector,crystal),0.);
        }
    }

    fNofBins["Statistics"] = env.GetValue("Histogram.Statistics.NofBins",64);
    fRangeLow["Statistics"] = env.GetValue("Histogram.Statistics.RangeLow.keV",0.);
    fRangeHigh["Statistics"] = env.GetValue("Histogram.Statistics.RangeHigh.keV",64.);

    fNofBins["Griffin1D"] = env.GetValue("Histogram.1D.Griffin.NofBins",4096);
    fRangeLow["Griffin1D"] = env.GetValue("Histogram.1D.Griffin.RangeLow.keV",0.5);
    fRangeHigh["Griffin1D"] = env.GetValue("Histogram.1D.Griffin.RangeHigh.keV",4096.5);

    fNofBins["0RES_Griffin1D"] = env.GetValue("Histogram.1D.Griffin.NofBins",4096);
    fRangeLow["0RES_Griffin1D"] = env.GetValue("Histogram.1D.Griffin.RangeLow.keV",0.5);
    fRangeHigh["0RES_Griffin1D"] = env.GetValue("Histogram.1D.Griffin.RangeHigh.keV",4096.5);

    fNofBins["Griffin2D"] = env.GetValue("Histogram.2D.Griffin.NofBins",4096);
    fRangeLow["Griffin2D"] = env.GetValue("Histogram.2D.Griffin.RangeLow.keV",0.5);
    fRangeHigh["Griffin2D"] = env.GetValue("Histogram.2D.Griffin.RangeHigh.keV",4096.5);

    fNofBins["0RES_Griffin2D"] = env.GetValue("Histogram.2D.Griffin.NofBins",4096);
    fRangeLow["0RES_Griffin2D"] = env.GetValue("Histogram.2D.Griffin.RangeLow.keV",0.5);
    fRangeHigh["0RES_Griffin2D"] = env.GetValue("Histogram.2D.Griffin.RangeHigh.keV",4096.5);

    //fNofBins["Griffin3D"] = env.GetValue("Histogram.3D.Griffin.NofBins",500);
    //fRangeLow["Griffin3D"] = env.GetValue("Histogram.3D.Griffin.RangeLow.keV",0.5);
   // fRangeHigh["Griffin3D"] = env.GetValue("Histogram.3D.Griffin.RangeHigh.keV",500.5);
        
    fNofBins["GriffinND"] = env.GetValue("Histogram.ND.Griffin.NofBins",500);
    fRangeLow["GriffinND"] = env.GetValue("Histogram.ND.Griffin.RangeLow.keV",0.5);
    fRangeHigh["GriffinND"] = env.GetValue("Histogram.ND.Griffin.RangeHigh.keV",500.5);

    fNofBins["Descant1D"] = env.GetValue("Histogram.1D.Descant.NofBins",4096);
    fRangeLow["Descant1D"] = env.GetValue("Histogram.1D.Descant.RangeLow.keV",0.5);
    fRangeHigh["Descant1D"] = env.GetValue("Histogram.1D.Descant.RangeHigh.keV",4096.5);

    fNofBins["0RES_Descant1D"] = env.GetValue("Histogram.1D.Descant.NofBins",4096);
    fRangeLow["0RES_Descant1D"] = env.GetValue("Histogram.1D.Descant.RangeLow.keV",0.5);
    fRangeHigh["0RES_Descant1D"] = env.GetValue("Histogram.1D.Descant.RangeHigh.keV",4096.5);
    
    fNofBins["Testcan1D"] = env.GetValue("Histogram.1D.Testcan.NofBins",4096);
    fRangeLow["Testcan1D"] = env.GetValue("Histogram.1D.Testcan.RangeLow.keV",0.5);
    fRangeHigh["Testcan1D"] = env.GetValue("Histogram.1D.Testcan.RangeHigh.keV",4096.5);

    fNofBins["0RES_Testcan1D"] = env.GetValue("Histogram.1D.Testcan.NofBins",4096);
    fRangeLow["0RES_Testcan1D"] = env.GetValue("Histogram.1D.Testcan.RangeLow.keV",0.5);
    fRangeHigh["0RES_Testcan1D"] = env.GetValue("Histogram.1D.Testcan.RangeHigh.keV",4096.5);

    fNofBins["TISTAR1D"] = env.GetValue("Histogram.1D.TISTAR.NofBins",4096);
    fRangeLow["TISTAR1D"] = env.GetValue("Histogram.1D.TISTAR.RangeLow.keV",0.5);
    fRangeHigh["TISTAR1D"] = env.GetValue("Histogram.1D.TISTAR.RangeHigh.keV",8192.0);
    
    fNofBins["0RES_TISTAR1D"] = env.GetValue("Histogram.1D.TISTAR.NofBins",4096);
    fRangeLow["0RES_TISTAR1D"] = env.GetValue("Histogram.1D.TISTAR.RangeLow.keV",0.5);
    fRangeHigh["0RES_TISTAR1D"] = env.GetValue("Histogram.1D.TISTAR.RangeHigh.keV",8192.0);

}
