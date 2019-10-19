 #include "Kinematics.hh"      //what is "part"

Kinematics::Kinematics(){
}

Kinematics::Kinematics(Nucleus* projectile, Nucleus* target, double ebeam){
	fParticle[0] = projectile;
	fParticle[1] = target;
	fParticle[2] = NULL;
	fParticle[3] = NULL;
	fM[0]=fParticle[0]->GetMass();
	fM[1]=fParticle[1]->GetMass();
	fEBeam = ebeam;
	fQValue =0;
	Initial(); // kinematic of incoming particles in lab and CM ***LA***
	FinalCm(); // kinematic of outgoing particles in lab and CM ***LA***
}
Kinematics::Kinematics(Nucleus* projectile, Nucleus* target, Nucleus* recoil, Nucleus* ejectile, double ebeam, double ex3){
	fParticle[0] = projectile;
	// std::cout << "projectile " << projectile << std::endl;
	fParticle[1] = target;
	fParticle[2] = recoil;
	fParticle[3] = ejectile;
	for(int i=0;i<4;i++)
		fM[i]=fParticle[i]->GetMass();

	fEBeam = ebeam;
	fQValue = (fM[0]+fM[1])-(fM[2]+fM[3])-ex3;
	//std::cout << "Qvalue :" << fQValue << std::endl;  
	Initial();
	FinalCm();
}

TSpline3* Kinematics::Evslab(double thmin, double thmax, double size, int part, bool upper){
	//std::cout << "Evslab 1 maximum scattering angle: " << GetMaxAngle(fVcm[part])*180./PI  << std::endl;

	if(thmax > GetMaxAngle(fVcm[part])*180./PI){
		thmax = GetMaxAngle(fVcm[part])*180./PI;
	}
	int nBins = (int)((thmax-thmin)/size)+1;
	//std::cout << "Evslab  2 max " << thmax << " min " << thmin << " steps " << (int)((thmax-thmin)/size)+1 << "size " << size << "part " << part << std::endl;
	double* energy = new double[nBins];
	double* angle = new double[nBins];
	int number =0;
	for(int i=0;i<nBins;i++){
		// changed by S.Klupp
		//Final((thmin+i*size)*PI/180.,2);
		Final((thmin+i*size)*PI/180.,part, upper);
		angle[i]=GetThetalab(part)*180./PI;
		energy[i]=GetTlab(part)*1000;
		//std::cout << "Evslab 3 angle " << angle[i] << " energy " << energy[i] << " Final_bla " <<  (thmin+i*size)*PI/180. <<  std::endl;
		//std::cout << "Evslab 4 GetThetalab " << GetThetalab(part)*180./PI << "GetTlab "<< GetTlab(part)*1000 << std::endl;
		if(energy[i]<1e15 && energy[i]>0.0001) number++;
		else break;

		 //std::cout << "Evslab 5 theta " <<(thmin+i*size)<<" res "<< angle[i] << " energy " << energy[i] << std::endl;
	}

	//std::cout << "Evslab 6 nBins " << nBins << " number " << number << std::endl;

	TGraph* graph = new TGraph(number, angle, energy);
	TSpline3* spline = new TSpline3("ETh_lab",graph);
	delete graph;
	delete[] angle;
	delete[] energy;
	return spline;
}

TSpline3* Kinematics::Evscm(double thmin, double thmax, double size, int part){
	double* energy = new double[(int)((thmax-thmin)/size)+1];
	double* angle = new double[(int)((thmax-thmin)/size)+1];
	int number =0;
	for(int i=0;i<((thmax-thmin)/size);i++){
		Final((thmin+i*size)*PI/180.,2);
		angle[i]=GetThetacm(part)*180./PI;
		energy[i]=GetTlab(part)*1000;
		number++;
	}
	TGraph* graph = new TGraph(number, angle, energy);
	TSpline3* spline = new TSpline3("ETh_cm",graph);
	delete graph;
	delete[] angle;
	delete[] energy;
	return spline;
}


double Kinematics::GetExcEnergy(TLorentzVector recoil, bool verbose){    //doing what? DURCHRECHNEN!
	if(verbose) {
		std::cout<<"recoil mass "<<fParticle[2]->GetMass()*1000.<<" keV/c^2, ejectile mass "<<fParticle[3]->GetMass()*1000.<<" keV/c^2"<<std::endl;
		std::cout<<"from recoil tot. energy "<<recoil.E()<<" keV, kin. energy "<<recoil.E()-fParticle[2]->GetMass()*1000.<<" keV"<<std::endl;
	}
	TLorentzVector ejectile;
	recoil.Boost(0,0,-GetBetacm()); //boost to cm system

	if(verbose) std::cout<<"boosted by "<<GetBetacm()<<": recoil e_cm "<<recoil.E()<<" keV"<<std::endl;
	ejectile.SetVect( -recoil.Vect() ); //pr = -pe
	ejectile.SetE(GetCmEnergy()*1000. - recoil.E()); //Ee=Ecm-Er
	if(verbose) std::cout<<GetCmEnergy()*1000.<<": ejectile e_cm "<<ejectile.E()<<std::endl;
	ejectile.Boost(0,0,GetBetacm()); //boost to lab system

	double eex = ejectile.M() - fParticle[3]->GetMass()*1000.;

	if(verbose) std::cout<<"from ejectile "<<ejectile.E()<<", "<<ejectile.M()<<" => eex "<<eex<<std::endl;
	return eex;
}

double Kinematics::GetBeamEnergy(double LabAngle, double LabEnergy){
	double ProjectileMass=fM[0]*1000;
	double TargetMass=fM[1]*1000;

	double ts = pow(TargetMass,2);
	double ps = pow(ProjectileMass,2);
	double cs = pow(cos(LabAngle),2);
	double es = pow(LabEnergy,2);
	double te = TargetMass*LabEnergy;

	return (-8*ProjectileMass - 4*TargetMass + LabEnergy/cs - 2*LabEnergy*tan(LabAngle) 
			+ sqrt(16*ts*cs + LabEnergy*(LabEnergy/cs*pow(cos(LabAngle) - sin(LabAngle),4) + 8*TargetMass*(3 + sin(2*LabAngle))))/cos(LabAngle) 
			+ sqrt((24*ts*LabEnergy*cs + 
					32*ps*TargetMass*pow(cos(LabAngle),4) + 
					2*TargetMass*es*pow(cos(LabAngle) - sin(LabAngle),4) + 
					16*ts*LabEnergy*pow(cos(LabAngle),3)*sin(LabAngle) - 
					8*ps*LabEnergy*cs*(sin(2*LabAngle) - 1) + 
					2*ps*cos(3*LabAngle)*sqrt(2*(4*ts + 12*te + es) + 
						(8*ts - 2*es)*cos(2*LabAngle) + 
						LabEnergy*(LabEnergy/cs + 8*TargetMass*sin(2*LabAngle) - 4*LabEnergy*tan(LabAngle))) + 
					2*cos(LabAngle)*(3*ps + te - te*sin(2*LabAngle))*
					sqrt(2*(4*ts + 12*te + es) + (8*ts - 2*es)*cos(2*LabAngle) + LabEnergy*(LabEnergy/cs + 8*TargetMass*sin(2*LabAngle) - 4*LabEnergy*tan(LabAngle))))/
				(pow(cos(LabAngle),4)*TargetMass)))/8.;
}


void Kinematics::Initial(){
	fT[0]=fEBeam; // projectile kinetic energy in lab ***LA***
	fT[1]=0; // target kinetic energy in lab ***LA***
	fE[0]=fT[0] + fM[0]; // projrctile total energy in lab ***LA***
	fE[1]=fT[1] + fM[1]; // target total energy in lab ***LA***
	fP[0]=sqrt(fT[0]*fT[0]+2*fT[0]*fM[0]); // momentum in lab ***LA***
	fP[1]=0;
	fV[0]=fP[0]/fE[0]; // velocity in lab ***LA***
	fV[1]=fP[1]/fE[1];

	fEcm[0]=GetCmEnergy(fEBeam)/2+(fM[0]*fM[0]-fM[1]*fM[1])/(2*GetCmEnergy(fEBeam));
	fEcm[1]=GetCmEnergy(fEBeam)/2-(fM[0]*fM[0]-fM[1]*fM[1])/(2*GetCmEnergy(fEBeam));
	fTcm[0]=fEcm[0]-fM[0];
	fTcm[1]=fEcm[1]-fM[1];
	fPcm[0]=Pcm_em(fEcm[0],fM[0]);
	fPcm[1]=Pcm_em(fEcm[1],fM[1]);
	fVcm[0]=fPcm[0]/fEcm[0];
	fVcm[1]=fPcm[1]/fEcm[1];
	fBeta_cm = (fP[0]-fP[1])/(fE[0]+fE[1]);
	fBetacm[0] = betacm_tm(fTcm[0],fM[0]);
	fBetacm[1] = -betacm_tm(fTcm[1],fM[1]);
	fGamma_cm = 1/sqrt(1-fBeta_cm*fBeta_cm);
	fTCm_i = GetCmEnergy(fEBeam)-fM[0]-fM[1];
	fTCm_f = fTCm_i + fQValue;

}
void Kinematics::FinalCm(){
	if(fParticle[2]==NULL && fParticle[3]==NULL){
		std::cout << "warning: outgoing particles not defined! Assuming elastic scattering" << std::endl;
		std::cout << "recoil = target" << std::endl;
		std::cout << "ejectile = projectile" << std::endl;
		fM[2]=fParticle[1]->GetMass();
		fM[3]=fParticle[0]->GetMass();     
	}
	fTcm[2]=fTCm_f/2*(fTCm_f+2*fM[3])/GetCmEnergy(fEBeam);
	fTcm[3]=fTCm_f/2*(fTCm_f+2*fM[2])/GetCmEnergy(fEBeam);
	fEcm[2]=fTcm[2] + fM[2];
	fEcm[3]=fTcm[3] + fM[3];
	fPcm[2]=Pcm_em(fEcm[2],fM[2]);
	fPcm[3]=Pcm_em(fEcm[3],fM[3]);
	fVcm[2]=fPcm[2]/fEcm[2];
	fVcm[3]=fPcm[3]/fEcm[3];
	fBetacm[2]=-betacm_tm(fTcm[2],fM[2]);
	fBetacm[3]=betacm_tm(fTcm[3],fM[3]);

	
		for(int i=0;i<4;i++){
		std::cout << "fBetacm["<<i<<"] = " << fBetacm[i]<< "\t";
		}
		std::cout <<std::endl;
		for(int i=0;i<4;i++){
		std::cout << "fPcm["<<i<<"] = " << fPcm[i]<< "\t";
		}
		std::cout <<std::endl;
		std::cout << "fBeta_cm = " << fBeta_cm<<std::endl;
		
}
void Kinematics::Final(double angle, int part, bool upper){             //angle of proton in lab system
	if(angle>GetMaxAngle(fVcm[part])){
		SetAngles(0, part);
		//std::cout<<"Final, angle = "<<angle<<" > MaxAngle = "<<GetMaxAngle(fVcm[part])<<std::endl;
	} else {
		// changed by S.Klupp
		//SetAngles(angle, part); 
		SetAngles(angle, part, upper);
	}
	fE[2] = E_final(2);
	fE[3] = E_final(3);
	fT[2] = T_final(2);
	fT[3] = T_final(3);
	//fP[2]=fGamma_cm*(fPcm[2]+fBeta_cm*fEcm[2]);
	//std::cout << "Vr = " << V_pe(fP[2],fE[2]) << std::endl; 
	//fP[3]=fGamma_cm*(fPcm[3]+fBeta_cm*fEcm[3]);
	//std::cout << "Ve = " << V_pe(fP[3],fE[3]) << std::endl; 
	fP[2] = P_tm(fT[2],fM[2]);
	//fP[2]=fGamma_cm*fPcm[2]*(cos(fThetacm[2])+fBeta_cm/fVcm[2]);
	//  std::cout << fTheta[2]*180./PI << "\t" << GetTlab(2)  << "\tVr = " << V_pe(fP[2],fE[2]) << "\tPr = " << fP[2] << "\tfT[2] "<< fT[2] <<"\tfM[2] " << fM[2] << std::endl; 
	fP[3] = P_tm(fT[3],fM[3]);
	//fP[3]=fGamma_cm*fPcm[3]*(cos(fThetacm[3])+fBeta_cm/fVcm[3]);
	//  std::cout << fTheta[3]*180./PI << "\t" << GetTlab(3)  << "\tVe = " << V_pe(fP[3],fE[3]) << "\tPe = " << fP[3] <<"\tfT[3] "<< fT[3] <<"\tfM[3] " << fM[3] << std::endl; 
	fV[2] = fP[2]/fE[2];
	fV[3] = fP[3]/fE[3];
}
double Kinematics::ELab(double angle_lab, int part){
	Final(angle_lab, part);
	return GetTlab(part);
}
void Kinematics::SetAngles(double angle, int part, bool upper){    
	int given;
	int other;
	if(part==2){
		given =2;
		other =3;
	}
	else if(part==3){
		given =3;
		other =2;
	}
	else{
		std::cout << " error in Kinematics::SetAngles("<<angle<<", "<<part<<") "<<std::endl;
		std::cout << " part must be 2 or 3 " << std::endl;
		exit(4);
	} 
	fTheta[given]=angle;
	//std::cout << "fTheta " << angle *180./TMath::Pi() << std::endl;
	fThetacm[given]=Angle_lab2cm(fVcm[given],fTheta[given]);
	//std::cout << "fThetacm " << fThetacm[given] << std::endl;
	if(given==3&&(fParticle[0]->GetMass()>fParticle[1]->GetMass())){
		//std::cout << "SetAngles: inverse kinematics" << std::endl;
		fThetacm[given]=Angle_lab2cminverse(fVcm[given],fTheta[given],upper);
	}
	//std::cout << "fThetacm " << fThetacm[given] *180./TMath::Pi() << std::endl;
	//std::cout << "------------------" << std::endl;
	fThetacm[other]=PI-fThetacm[given];
	if(fTheta[given]==0)
		fTheta[other]=PI/given;
	else{
		fTheta[other]=Angle_cm2lab(fVcm[other],fThetacm[other]);
	}
	//  fTheta[3]=angle;
	//  fThetacm[3]=Angle_lab2cm(fVcm[3],fTheta[3]);
	//  fThetacm[2]=PI-fThetacm[3];
	//  if(fTheta[3]==0)
	//    fTheta[2]=PI/2;
	//  else
	//    fTheta[2]=Angle_cm2lab(fVcm[2],fThetacm[2]);
}


double Kinematics::GetCmEnergy(double ebeam){
	double ecm;
	ecm = sqrt(fM[0]*fM[0]+fM[1]*fM[1]+2.*fM[1]*(fM[0]+ebeam));
	return ecm;
}
double Kinematics::GetCmEnergy(){
	return GetCmEnergy(fEBeam);
}
double Kinematics::NormalkinEnergy(){                   // doing what? 
	double ENorm = (GetCmEnergy(fEBeam)*GetCmEnergy(fEBeam)-fM[0]*fM[0]-fM[1]*fM[1])/(2*fM[0])-fM[1];
	return ENorm;
}

double Kinematics::GetMaxAngle(double vcm){
	double x;
	x = fBeta_cm/vcm;
	if(x*x<1)
		return PI;
	else
		return atan2(sqrt(1/(x*x-1)),fGamma_cm);
}
double Kinematics::GetMaxAngle(int part){
	return GetMaxAngle(fVcm[part]);
}
bool Kinematics::CheckMaxAngle(double angle, int part){
	return angle <= GetMaxAngle(fVcm[part]);
}
double Kinematics::Angle_lab2cm(double vcm, double angle_lab){
	double tan_lab, gtan,x;
	tan_lab = tan(angle_lab);
	gtan = tan_lab*tan_lab*fGamma_cm*fGamma_cm;
	x = fBeta_cm/vcm;

	if(tan_lab>=0){
		//std::cout << "tan_lab>=0" << std::endl;
		return acos( (-x*gtan+sqrt( 1+gtan*(1-x*x) ))/(1+gtan) );
	}
	else{
		//std::cout << "tan_lab<0" << std::endl;
		return acos( (-x*gtan-sqrt( 1+gtan*(1-x*x) ))/(1+gtan) );
	}
}
double Kinematics::Angle_lab2cminverse(double vcm, double angle_lab, bool upper){
	double tan_lab, gtan,x;
	tan_lab = tan(angle_lab);
	gtan = tan_lab*tan_lab*fGamma_cm*fGamma_cm;
	x = fBeta_cm/vcm;
	//std::cout << "angle_lab " << angle_lab*180./TMath::Pi() << std::endl;

	if(upper){
		return acos( (-x*gtan+sqrt( 1+gtan*(1-x*x) ))/(1+gtan) );
	}
	else{
		return acos( (-x*gtan-sqrt( 1+gtan*(1-x*x) ))/(1+gtan) );
	}
}
void Kinematics::AngleErr_lab2cm(double angle, double &err){
	double angle_lab = angle;
	angle = Angle_lab2cm(fVcm[2], angle_lab);
	err =  fabs(Angle_lab2cm(fVcm[2], angle_lab+err)-Angle_lab2cm(fVcm[2], angle_lab-err))/2.;
	/*
		double tang, tang2, gtang, g2,x;
		tang = tan(angle_lab);
		gtang = tang*tang*fGamma_cm*fGamma_cm;
		tang2 = tang*tang;
		g2= fGamma_cm*fGamma_cm;
		x = fBeta_cm/fVcm[2];

		if(tang>=0){
		err *= fabs((-2*x*tang*g2*(1*tang2)+tang*g2*(1-x*x)*(1+tang2)/sqrt(1+gtang*(1-x*x)))/(1-gtang) + 2*(-x*gtang+sqrt(1+gtang*(1-x*x))*(1+tang2)*tang*g2)/(1-gtang)/(1-gtang) )/sqrt(1-((-x*gtang+sqrt( 1+gtang*(1-x*x) ))/(1+gtang))*((-x*gtang+sqrt( 1+gtang*(1-x*x) ))/(1+gtang)));
		}
		else{
		err *= fabs((-2*x*tang*g2*(1*tang2)-tang*g2*(1-x*x)*(1+tang2)/sqrt(1+gtang*(1-x*x)))/(1-gtang) + 2*(-x*gtang-sqrt(1+gtang*(1-x*x))*(1+tang2)*tang*g2)/(1-gtang)/(1-gtang) )/sqrt(1-((-x*gtang-sqrt( 1+gtang*(1-x*x) ))/(1+gtang))*((-x*gtang-sqrt( 1+gtang*(1-x*x) ))/(1+gtang)));
		if(err>1){
		std::cout << "part1 " << (-2*x*tang*g2*(1*tang2)-tang*g2*(1-x*x)*(1+tang2)/sqrt(1+gtang*(1-x*x)))/(1-gtang) << std::endl;
		std::cout << "part2 " << 2*(-x*gtang-sqrt(1+gtang*(1-x*x))*(1+tang2)*tang*g2)/(1-gtang)/(1-gtang) << std::endl;
		}
		}
		*/
}
double Kinematics::Angle_cm2lab(double vcm, double angle_cm){
	double x;
	x = fBeta_cm/vcm;
	return atan2(sin(angle_cm),fGamma_cm*(cos(angle_cm)+x));
	/*
		std::cout << " old " << atan2(sin(angle_cm),fGamma_cm*(cos(angle_cm)+x)) << " x " << x << std::endl;
		double gam2 = fM[0]*fM[2]/fM[1]/fM[3]*fTCm_i/fTCm_f;
		gam2 = sqrt(gam2);
		double gam3 = fM[0]*fM[3]/fM[1]/fM[2]*fTCm_i/fTCm_f;
		gam3 = sqrt(gam3);
		double y2=1.+gam2*gam2+2.*gam2*cos(angle_cm);
		double y3=1.+gam3*gam3+2.*gam3*cos(PI-angle_cm);
		std::cout << "y2 " << y2 << " y3 " << y3 << "\n";
		if(asin(sin(PI-angle_cm)/sqrt(y2)) < PI-acos(-gam2))
		std::cout << " new " << PI - asin(sin(PI-angle_cm)/sqrt(y2)) << " = asin("<<sin(PI-angle_cm)<<"/"<<sqrt(y2)<<")" << std::endl;
		else 
		std::cout << " new " << asin(sin(PI-angle_cm)/sqrt(y2)) << " = asin("<<sin(PI-angle_cm)<<"/"<<sqrt(y2)<<")" << std::endl;
		*/

}
//x = sqrt(fM[0]*fM[3]/fM[1]/fM[2]*fTCm_i/fTCm_f);
//std::cout << "thorsten\t" << asin(sin(angle_cm)/sqrt(1+x*x+2*x*cos(angle_cm)))*180./PI << std::endl;
//std::cout << "ich\t" << atan2(sin(angle_cm),fGamma_cm*(cos(angle_cm)+x))*180./PI << std::endl;  
TSpline3* Kinematics::labvscm(double thmin, double thmax, double size, int part){
	double* cm = new double[(int)((thmax-thmin)/size)+1];
	double* lab = new double[(int)((thmax-thmin)/size)+1];
	int nr =0;
	for(int i=0;i<((thmax-thmin)/size);i++){
		cm[nr] = i;
		lab[nr] = Angle_cm2lab(fVcm[part],cm[nr]*PI/180.)*180./PI;
		if(lab[nr]>0.01 && lab[nr]<179.99)
			nr++;
	}
	TGraph* graph = new TGraph(nr, cm, lab);
	TSpline3* spline = new TSpline3("Th_cmvslab",graph);
	delete graph;
	delete[] lab;
	delete[] cm;
	return spline;
}
TSpline3* Kinematics::cmvslab(double thmin, double thmax, double size, int part){
	double* cm = new double[(int)((thmax-thmin)/size)+1];
	double* lab = new double[(int)((thmax-thmin)/size)+1];
	int nr =0;
	for(int i=0;i<((thmax-thmin)/size);i++){
		cm[nr] = i;
		lab[nr] = Angle_cm2lab(fVcm[part],cm[nr]*PI/180.)*180./PI;
		if(lab[nr]>0.01 && lab[nr]<179.99)
			nr++;
	}
	TGraph* graph = new TGraph(nr, lab, cm);
	TSpline3* spline = new TSpline3("Th_cmvslab",graph);
	delete graph;
	delete[] lab;
	delete[] cm;
	return spline;
}

double Kinematics::Sigma_cm2lab(double angle_cm, double sigma_cm){
	double gam2 = fM[0]*fM[2]/fM[1]/fM[3]*fTCm_i/fTCm_f;
	gam2 = sqrt(gam2);
	double wurzel=1.+gam2*gam2+2.*gam2*cos(PI-angle_cm);
	wurzel = sqrt(wurzel);
	return sigma_cm*(1/fGamma_cm*wurzel*wurzel*wurzel/(1+gam2*cos(PI-angle_cm)));
}
/*
	double Kinematics::Sigma_cm2lab(double angle_cm, double sigma_cm){
//old
double x;
x = fBeta_cm/fVcm[2];
double wurzel;
wurzel = sqrt(sin(PI-angle_cm)*sin(PI-angle_cm) + fGamma_cm*fGamma_cm*(cos(PI-angle_cm)+x)*(cos(PI-angle_cm)+x));
return sigma_cm*(1/fGamma_cm*wurzel*wurzel*wurzel/(1+x*cos(PI-angle_cm)));  
}
*/
double Kinematics::Sigma_lab2cm(double angle_cm, double sigma_lab){
	double gam2 = fM[0]*fM[2]/fM[1]/fM[3]*fTCm_i/fTCm_f;
	gam2 = sqrt(gam2);
	//test
	//double x;
	//x = fBeta_cm/fVcm[2];
	//std::cout << "x " << x << " gam2 "<< gam2 << " cm "<< angle_cm << " pi - cm "<< PI-angle_cm << std::endl;
	//test

	//double angle_cm = Angle_lab2cm(fVcm[2], angle_lab);
	//orig
	double wurzel=1.+gam2*gam2+2.*gam2*cos(PI-angle_cm);
	wurzel = sqrt(wurzel);
	return sigma_lab/(1/fGamma_cm*wurzel*wurzel*wurzel/(1+gam2*cos(PI-angle_cm)));  
	//test no pi-
	//double wurzel=1.+gam2*gam2+2.*gam2*cos(angle_cm);
	//wurzel = sqrt(wurzel);
	//return sigma_lab/(1/fGamma_cm*wurzel*wurzel*wurzel/(1+gam2*cos(angle_cm)));  
}
void Kinematics::SigmaErr_lab2cm(double angle, double err, double &sigma, double &errsigma){
	double g = fM[0]*fM[2]/fM[1]/fM[3]*fTCm_i/fTCm_f;
	g = sqrt(g);
	double w=1.+g*g+2.*g*cos(PI-angle);
	w = sqrt(w);
	errsigma = fGamma_cm/pow(w,1.5) * sqrt( pow(sigma*g*sin(PI-angle)*(-2+g*g-g*cos(PI-angle))/w * err,2) + pow((1+g*cos(PI-angle))*errsigma,2  ) ); 
	//sigma/=(1/fGamma_cm*w*w*w/(1+g*cos(PI-angle)));
}
void Kinematics::Transform2cm(double &angle, double &sigma){
	//double angle_lab = angle;
	angle = PI-Angle_lab2cm(fVcm[2], angle);
	sigma = Sigma_lab2cm(angle, sigma);
	return;
}
void Kinematics::Transform2cm(double &angle, double &errangle, double &sigma, double &errsigma){
	AngleErr_lab2cm(angle, errangle);
	Transform2cm(angle,sigma);
	SigmaErr_lab2cm(angle, errangle, sigma, errsigma);
	return;
}
double Kinematics::Rutherford(double angle_cm){
	double a = 0.5*1.43997649*fParticle[0]->GetZ()*fParticle[1]->GetZ()/fTCm_i;
	double b = sin(angle_cm/2.)*sin(angle_cm/2.);
	b=b*b;
	return a*a/b*0.0025;//1b=0.01fm
}
TSpline3* Kinematics::Ruthvscm(double thmin, double thmax, double size){
	double* cross = new double[(int)((thmax-thmin)/size)+1];
	double* angle = new double[(int)((thmax-thmin)/size)+1];
	int number =0;
	for(int i=0;i<((thmax-thmin)/size);i++){
		angle[i]=thmin+i*size;
		if(angle[i]>179.99||angle[i]<0.01)
			break;
		cross[i]=Rutherford(angle[i]*PI/180.);
		number++;

		//std::cout << angle[i] << "   " << Rutherford(GetThetacm(2))<<std::endl;
		//std::cout << (thmin+i*size)*PI/180. << "  max angle: " << GetMaxAngle(fVcm[2]) << std::endl;
	}
	TGraph* graph = new TGraph(number, angle, cross);
	TSpline3* spline = new TSpline3("sigmaTh_cm",graph);
	delete graph;
	delete[] angle;
	delete[] cross;
	return spline;
}
TSpline3* Kinematics::Ruthvslab(double thmin, double thmax, double size, int part){
	double* cross = new double[(int)((thmax-thmin)/size)+1];
	double* angle = new double[(int)((thmax-thmin)/size)+1];
	int number =0;
	for(int i=0;i<((thmax-thmin)/size);i++){
		if(part==3||part==2)
			angle[i]=thmin+i*size; //angle[i] is in cm system
		else{
			std::cout << "error " << std::endl;
			exit(1);
		}
		if(angle[i]>179.99||angle[i]<0.01)
			break;
		cross[i]=Rutherford(angle[i]*PI/180.);
		number++;
		//std::cout << " angle cm " << angle[i]; 
		//std::cout << " cs cm " << cross[i];

		cross[i]=Sigma_cm2lab(angle[i]*PI/180., Rutherford(angle[i]*PI/180.));
		if(part==2){
			angle[i]=180-angle[i];
		}
		angle[i]=Angle_cm2lab(fVcm[part],angle[i]*PI/180.)*180./PI;
		//std::cout << fVcm[part] << "fVcm[part]" << std::endl;
		//std::cout << "\t\tangle lab " << angle[i]; 
		//std::cout << " cs lab " << cross[i] << std::endl;; 
	}
	TGraph* graph = new TGraph(number, angle, cross);
	TSpline3* spline = new TSpline3("sigmaTh_lab",graph);
	delete graph;
	delete[] angle;
	delete[] cross;
	return spline;
}

double Kinematics::Pcm_em(double e, double m){     // Impuls cm 
	return sqrt(e*e-m*m);
}
double Kinematics::P_tm(double t, double m){
	return sqrt(t*t+2.*t*m);
}
double Kinematics::E_tm(double t, double m){
	return t+m;
}
double Kinematics::T_em(double e, double m){
	return e-m;
}
double Kinematics::betacm_tm(double t, double m){
	return sqrt(t*t+2*t*m)/(t+m);
}
double Kinematics::V_pe(double p, double e){         // 
	return p/e;
}
double Kinematics::E_final(int i){                  //  calc E_final NACHRECHNEN!  
	return fGamma_cm*(fEcm[i]+fBeta_cm*fPcm[i]);
}
double Kinematics::T_final(int i){                 // NACHRECHNEN!
	//std::cout <<(fGamma_cm-1)*fM[i]<<" + "<<fGamma_cm*fTcm[i]<<" + "<<fGamma_cm*fBeta_cm*fPcm[i]*cos(fThetacm[i])<<std::endl;
	//std::cout <<fGamma_cm<<"*"<<fBeta_cm<<"*"<<fPcm[i]<<"*"<<cos(fThetacm[i])
	return (fGamma_cm-1)*fM[i]+fGamma_cm*fTcm[i]+fGamma_cm*fBeta_cm*fPcm[i]*cos(fThetacm[i]); 
}
