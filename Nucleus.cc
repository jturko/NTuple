#include "Nucleus.hh"

#include <iostream>
#include <sstream>
#define debug

//static double amu = 931.494043; // old value
static double amu = 931.49406121;

Nucleus::Nucleus() {
}

Nucleus::Nucleus(const char* elementsymbol) {
	int length = strlen(elementsymbol);
	if(length == 0) {
		std::cerr<<"error, type numbersymbol, or symbolnumber, i.e. 30Mg oder Mg30"<<std::endl;
		exit(1);
	}
	int first_digit = length;
	int first_letter = length;

	if(isdigit(elementsymbol[0])) {
		first_digit = 0;
	} else if(isalpha(elementsymbol[0])) {
		first_letter = 0;
	} else {
		std::cerr<< "error, type numbersymbol, or symbolnumber, i.e. 30Mg oder Mg30"<<std::endl;
		exit(1);
	}
	int symbol_length = 0;
	int i;
	for(i=0;i<length;i++) {
		if(isdigit(elementsymbol[i]) && first_letter == 0) {
			first_digit = i;
			symbol_length = i;
			break;
		}
		if(isalpha(elementsymbol[i]) && first_digit == 0) {
			first_letter = i;
			symbol_length = length - i;
			break;
		}
	}
	symbol_length +=1; //char one digit more
	if(i==length) {
		if(first_digit == 0 || first_letter == 0)
			std::cerr<< "error, type numbersymbol, or symbolnumber, i.e. 30Mg oder Mg30"<<std::endl;
		exit(1);  
	}
	SetSymbol(elementsymbol+first_letter);
	GetZfromSymbol(elementsymbol+first_letter);
	SetN(atoi(elementsymbol+first_digit)-GetZ());
	SetMass(atof(elementsymbol+first_digit)*amu);
}

Nucleus::Nucleus(int charge, int neutrons, double mass, const char* symbol) {
	fZ = charge;  
	fN = neutrons;
	fSymbol = symbol;
	fMass = mass;
}

Nucleus::Nucleus(int charge, int neutrons, const char* MassFile) {
	fZ = charge;  
	fN = neutrons;
	int n;
	int z;
	double emass;
	double tmp;
	std::ifstream mass_file;
	mass_file.open(MassFile,std::ios::in);
	std::string line;
	std::istringstream str;
	while(getline(mass_file, line).good()) {
		str.clear();
		str.str(line);
		str>>z>>n;
		if(n==fN && z==fZ) {
			str>>emass>>tmp>>fSymbol;
			fMassExcess = emass/1000.;
#ifdef debug
			std::cout << "Symbol " << fSymbol << " tmp " << tmp <<std::endl;
#endif
			SetMass();
			mass_file.close();
			return;
		}
	}
	//max_elements=i;
	std::cout<<"Failed to find z = "<<charge<<", n = "<<neutrons<<" in \""<<MassFile<<"\"!"<<std::endl;
	mass_file.close();
}

void Nucleus::SetZ(int charge) {
	fZ = charge;
}

void Nucleus::SetN(int neutrons) {
	fN = neutrons;
}

void Nucleus::SetMassExcess(double mass_ex) {
	fMassExcess = mass_ex;
}

void Nucleus::SetMass(double mass) {
	fMass = mass;
}

void Nucleus::SetMass() {
	fMass = amu*GetA()+GetMassExcess();
}

void Nucleus::SetSymbol(const char* symbol) {
	fSymbol = symbol;
}

int Nucleus::GetZfromSymbol(const char* symbol) {
	char symbols[105][3] = {"H","HE","LI","BE","B","C","N","O","F","NE","NA","MG","AL","SI","P","S","CL","AR","K","CA","SC","TI","V","CR","MN","FE","CO","NI","CU","ZN","GA","GE","AS","SE","BR","KR","RB","SR","Y","ZR","NB","MO","TC","RU","RH","PD","AG","CD","IN","SN","SB","TE","F","XE","CS","BA","LA","CE","PR","ND","PM","SM","EU","GD","TB","DY","HO","ER","TM","YB","LU","HF","TA","W","RE","OS","IR","PT","AU","HG","TI","PB","BI","PO","AT","RN","FR","RA","AC","TH","PA","U","NP","PU","AM","CM","BK","CF","ES","FM","MD","NO","LR","RF","HA"};
	int length = strlen(symbol);
	//std::cout << symbol << "   " << length << std::endl;
	char* search = new char[length+1];
	for(int i=0;i<length;i++) {
		search[i] = toupper(symbol[i]); // make sure symbol ist in uppercase
	}
	search[length] = '\0';
	for(int i=0;i<105;i++) {
		if(strcmp(search,symbols[i]) == 0) {
			delete[] search;
			SetZ(i+1);
			return i+1;
		}
	}

	delete[] search;
	SetZ(0);
	return 0;

}

int Nucleus::GetZ() {
	return fZ;
}

int Nucleus::GetN() {
	return fN;
}

int Nucleus::GetA() {
	return fN+fZ;
}

double Nucleus::GetMassExcess() {
	return fMassExcess;
}

double Nucleus::GetMass() {
	return fMass;
}

double Nucleus::GetRadius() {
	return 1.12*pow(this->GetA(),1./3.) - 0.94*pow(this->GetA(),-1./3.);
}

const char* Nucleus::GetSymbol() {
	return fSymbol.c_str();
}
