/***********************************************************************************
 * S.Klupp: Energy loss calculations e.g. of the beam in the target
 * adapted from irma.f
 **********************************************************************************/

#ifndef __RECONSTRUCTION_HH
#define __RECONSTRUCTION_HH

#include <iostream>
#include <fstream>
#include <string>

#include "TMath.h"
#include "TSpline.h"
#include "TGraph.h"

#include "Nucleus.hh"
#include "Compound.hh"


#ifndef PI
#define PI                       (TMath::Pi())
#endif

class Reconstruction {
	public:
		Reconstruction();
		Reconstruction(Nucleus* projectile, Compound* target);
		Reconstruction(Nucleus* projectile, Compound* target, double thickness);
		void SetTargetThickness(double thickness){
			fTargetThickness = thickness;
		}
		void SetProj(Nucleus* projectile){
			fProj = projectile;
		}
		void SetTarget(Compound* target){
			fTarget = target;
		}
		double StoppingPower(double energy, bool gaseous = true); // total stopping power of the target
		double StoppingPower(Nucleus* target, double energy, bool gaseous = true); // stopping power of a single target element
		double CompoundRange(double energy, int limit, bool gaseous = true); // range of the projectile in the target
		double EnergyAfter(double energy, int limit, bool gaseous = true); // energy after target        //ueberall auf true fuer Gastarget
		double EnergyLoss(double energy, int limit, bool gaseous = true){ // energy loss in target
			return energy - EnergyAfter(energy, limit, gaseous);
		};
		TSpline3* Energy2Range(double emax, double size, bool gaseous = true);
		TSpline3* Range2Energy(double emax, double size, bool gaseous = true);
		TSpline3* Energy2EnergyLoss(double emax, double size, bool gaseous = true);
		TSpline3* Energy2EnergyAfter(double emax, double size, bool gaseous = true);
		TSpline3* Thickness2EnergyAfter(double energy, double maxThickness, double stepSize, bool gaseous);
		TGraph* EnergyAfter2Energy(double emax, double size, bool gaseous = true);
	private:
		double fTargetThickness;
		Nucleus* fProj;
		Compound* fTarget;
		// the following arrays are needed for the calulation of the stopping power
		double a_h(int index, int z); 
		double a_he(int index, int z);
		double b_he(int index, int z);
		double shell_correction(int z); // shell correction for Bethe formula for hydrogen as projectile
		//ClassDef(Reconstruction, 1);
};

#endif
