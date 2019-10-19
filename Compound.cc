#include "Compound.hh"
#include "TistarSettings.hh"

#define debug

Compound::Compound() {
}

Compound::Compound(const char* symbol) {
	std::string massfile = TistarSettings::Get()->GetMassFile();
    //std::string massfile = SIMULATION_PATH;
	//massfile.append("/../MassFile.dat");
	
    int length = strlen(symbol);
	if(length == 0) {
		std::cerr<<"Error, length of Material string is zero!"<<std::endl;
		exit(1);
	}
	fSymbol = symbol;
	fDensity = 1.;
	if(strstr(symbol,"PE")) {
		std::cout << "Polyethylene!" << std::endl;
		//H4C2
		fNuclei.resize(2);
		fFrac.resize(2);

		fNuclei[0] = new Nucleus(1,0,massfile.c_str());
		fNuclei[1] = new Nucleus(6,6,massfile.c_str());

		fFrac[0] = 4.*fNuclei[0]->GetMass()/(4.*fNuclei[0]->GetMass() + 2.*fNuclei[1]->GetMass());
		fFrac[1] = 2.*fNuclei[1]->GetMass()/(4.*fNuclei[0]->GetMass() + 2.*fNuclei[1]->GetMass());

		fMass = fNuclei[0]->GetMass()*4. + fNuclei[1]->GetMass()*2.;
		fDensity = 0.94; // from Geant4
	} else if(strstr(symbol,"DPE")) {
		std::cout << "deuterated Polyethylene!" << std::endl;
		//D4C2
		fNuclei.resize(2);
		fFrac.resize(2);

		fNuclei[0] = new Nucleus(1,1,massfile.c_str());
		fNuclei[1] = new Nucleus(6,6,massfile.c_str());

		fFrac[0] = 4.*fNuclei[0]->GetMass()/(4.*fNuclei[0]->GetMass() + 2.*fNuclei[1]->GetMass());
		fFrac[1] = 2.*fNuclei[1]->GetMass()/(4.*fNuclei[0]->GetMass() + 2.*fNuclei[1]->GetMass());

		fMass = fNuclei[0]->GetMass()*4. + fNuclei[1]->GetMass()*2.;
		fDensity = 1.06; //found in a paper. Assuming same density of molecules, the density of 0.94 for PE would correspond to 1.074 for DPE (0.94*32/28).
	} else if(strstr(symbol,"MY")) {
		std::cout << "Mylar!" << std::endl;
		//H8C10O4
		fNuclei.resize(3);
		fFrac.resize(3);

		fNuclei[0] = new Nucleus(1,0,massfile.c_str());
		fNuclei[1] = new Nucleus(6,6,massfile.c_str());
		fNuclei[2] = new Nucleus(8,8,massfile.c_str());

		fFrac[0] = 8.*fNuclei[0]->GetMass()/(8.*fNuclei[0]->GetMass() + 10.*fNuclei[1]->GetMass() + 4.*fNuclei[2]->GetMass());
		fFrac[1] = 10.*fNuclei[1]->GetMass()/(8.*fNuclei[0]->GetMass() + 10.*fNuclei[1]->GetMass() + 4.*fNuclei[2]->GetMass());
		fFrac[2] = 4.*fNuclei[2]->GetMass()/(8.*fNuclei[0]->GetMass() + 10.*fNuclei[1]->GetMass() + 4.*fNuclei[2]->GetMass());

		fMass = fNuclei[0]->GetMass()*8. + fNuclei[1]->GetMass()*10. + fNuclei[2]->GetMass()*4.;
		fDensity = 1.4; // from Geant4 or google 1.39 g/cm3
	} else if(strstr(symbol,"TTI")) {
		std::cout << "Tritiated Titanium Target!" << std::endl;
		// ratioTTI = atomic ratio Tritium/Titanium
		if(isalpha(symbol[0])) {
			std::cerr<<"give atomic ratio of Tritium to Titanium!"<< std::endl;
			exit(1);
		} else{
			double ratio = atof(symbol);
			//std::cout << ratio << std::endl;
			fNuclei.resize(2);
			fFrac.resize(2);
			//     std::cout << "fFrac " << fFrac << std::endl;
			fNuclei[0] = new Nucleus(1,2,massfile.c_str());
			fNuclei[1] = new Nucleus(18,26,massfile.c_str());
			fFrac[0] = ratio*fNuclei[0]->GetMass()/(ratio*fNuclei[0]->GetMass()+fNuclei[1]->GetMass());
			fFrac[1] = fNuclei[1]->GetMass()/(ratio*fNuclei[0]->GetMass()+fNuclei[1]->GetMass());
			fMass = fNuclei[0]->GetMass()*ratio + fNuclei[1]->GetMass();
			fDensity = 4.506; // titanium density
		}
	} else if(strstr(symbol,"DTI")) {
		//std::cout << "Deuterated Titanium Target!" << std::endl;
		if(isalpha(symbol[0])) {
			std::cerr<<"give atomic ratio of Deuterium to Titanium!"<< std::endl;
			exit(1);
		} else{
			double ratio = atof(symbol);
			std::cout << ratio << std::endl;
			fNuclei.resize(2);
			fFrac.resize(2);

			fNuclei[0] = new Nucleus(1,1,massfile.c_str());
			fFrac[0] = ratio/(1+ratio);
			fNuclei[1] = new Nucleus(22,26,massfile.c_str());
			fFrac[1] = 1/(1+ratio);

			fMass = fNuclei[0]->GetMass()*fFrac[0] + fNuclei[1]->GetMass()*fFrac[1];
			fDensity = 4.506; // titanium density
		}
	} else if(strstr(symbol,"2H")) {
		std::cout << "DeuteriumGas!" << std::endl;
		fNuclei.resize(1);
		fFrac.resize(1);

		fNuclei[0] = new Nucleus(1,1,massfile.c_str());

		fFrac[0] = 1;

		fMass = fNuclei[0]->GetMass() * 2.0;              //*2.0 since it's a molecule
		fDensity = 0.000164;
		
	} else if(strstr(symbol,"SolidDeuterium")) {
		std::cout << "Solid Deuterium!" << std::endl;
		fNuclei.resize(1);
		fFrac.resize(1);

		fNuclei[0] = new Nucleus(1,1,massfile.c_str());

		fFrac[0] = 1;

		fMass = fNuclei[0]->GetMass() * 2.0;              //*2.0 since it's a molecule
		fDensity = 0.1967;
	}else if(strstr(symbol,"helium")) {
		std::cout << "Helium Gas!" << std::endl;
		fNuclei.resize(1);
		fFrac.resize(1);

		fNuclei[0] = new Nucleus(2,2,massfile.c_str());

		fFrac[0] = 1;

		fMass = fNuclei[0]->GetMass() * 1.0;
		std::cerr<<"symbol " << symbol<< std::endl;
		//std::cout<<"Helium Mass Leila "<<fMass<<" zero: "<<fNuclei[0]->GetMass()<<" one: "<<fNuclei[1]->GetMass()<<std::endl;
		fDensity = 0.179e-3;// 0.179e-3 g/cm3 @ STP (0 deg/1 bar)
	} 
	
	else if(strstr(symbol,"silicon")) {
		std::cout << "silicon!" << std::endl;
		fNuclei.resize(1);
		fFrac.resize(1);

        std::cout << "massfile = " << massfile << std::endl;
		fNuclei[0] = new Nucleus(14,14,massfile.c_str());
        //std::cout << "just created new nucleus w/ Z = " << fNuclei[0]->GetZ() << ", A = " << fNuclei[0]->GetZ() << std::endl;
        //std::cout << "mass is " << std::flush << fNuclei[0]->GetMass() << std::endl;

		fFrac[0] = 1;

		fMass = fNuclei[0]->GetMass() * 1.0;
		//std::cout<<"silicon Mass Leila "<<fMass<<" zero: "<<fNuclei[0]->GetMass()<<" one: "<<fNuclei[1]->GetMass()<<std::endl;
		fDensity = 2.329;// 2.329 g/cm3 from google
	} 

	else {
		std::cerr<<"Compound \""<<symbol<<"\" not implemented yet!"<< std::endl;
		std::cerr<<"symbol " << symbol<< std::endl;
		exit(1);
	}
}

Compound::Compound(Nucleus* target) {
	fNuclei.resize(1);
	fFrac.resize(1);
	fNuclei[0] = target;
	fFrac[0] = 1;
	fDensity = 1.;
}

Compound::~Compound() {
	// if(fFrac!=NULL)
	//   delete[] fFrac;
	// if(fNuclei!=NULL) {
	//   for(int i=0;i<GetNofElements();i++) {
	//     if(fNuclei[i]!=NULL)
	// 	delete fNuclei[i];
	//   }
	//   delete[] fNuclei;
	// }
	// if(fSymbol!=NULL)
	//   delete[] fSymbol;
}

Nucleus* Compound::GetNucleus(size_t i) {
	if(i < fNuclei.size()) return fNuclei[i];
	else		              return NULL;
}

double Compound::GetFrac(size_t i) {
	// std::cout << "i" << i << "GetNofElements()" << GetNofElements() << std::endl;
	if(i < fFrac.size())	return fFrac[i];
	else           		return 0;
}
