#ifndef __COMPOUND_HH
#define __COMPOUND_HH

#include <iostream>
#include <math.h>
#include <string.h>
#include <iomanip>
#include <stdlib.h>
#include <fstream>
#include <vector>

#include "Nucleus.hh"

class Compound {
	public:
		Compound();
		Compound(const char*);
		Compound(Nucleus*);
		~Compound();
		//Compound(int nofelements, Nucleus* nuclei, double* fracs);
		void SetDensity(double density) { fDensity = density; }
		size_t GetNofElements() { return fNuclei.size(); }
		double GetMass()        { return fMass; }
		double GetDensity()     { return fDensity; }
		const char* GetSymbol() { return fSymbol; }
		Nucleus* GetNucleus(size_t);
		double GetFrac(size_t);
	private:
		std::vector<Nucleus*> fNuclei;
		std::vector<double> fFrac;
		double fMass;
		double fDensity;
		const char* fSymbol;
};
#endif
