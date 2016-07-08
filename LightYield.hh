

double deuteronLight (double E) { return 0.83*E-2.82*(1-TMath::Exp(-0.25*TMath::Power(E,0.93))); }
double protonLight (double E) { return 0.83*E-2.82*(1-TMath::Exp(-0.25*TMath::Power(E,0.93))); }
double carbonLight(double E)    { return 0.017*E; }

double resolutionSigma(double L) { return L*TMath::Sqrt(TMath::Power(0.15,2)+TMath::Power(0.1,2)/L+TMath::Power(0.02/L,2)); }
double resolutionFWHM(double L) { return L*TMath::Sqrt(TMath::Power(0.15,2)+TMath::Power(0.1,2)/L+TMath::Power(0.02/L,2))/(2.*TMath::Sqrt(2.*TMath::Log(2.))); }
