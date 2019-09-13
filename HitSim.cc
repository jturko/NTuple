#include "HitSim.hh"

HitSim::HitSim(Settings* setting) {
	fSett = setting;
    fTRexSett = fSett->GetTRexSettings();
	fRand = new TRandom();
	Clear();
} // end constructor


void HitSim::Clear() {
	fFirstDeltaE = nullptr;
	fSecondDeltaE = nullptr;
	fPad = nullptr;
	fFirstPosition.SetXYZ(0., 0., 0.);
	fSecondPosition.SetXYZ(0., 0., 0.);
	fFirstEnergy = -99.;
	fSecondEnergy = -99.;
}


void HitSim::SetFirstDeltaE(ParticleMC& firstDeltaE, Direction direction) {
	fFirstDeltaE = &firstDeltaE;
	fFirstDirection = direction;
	fFirstPosition.SetXYZ(0., 0., 0.);
}

void HitSim::SetSecondDeltaE(ParticleMC& secondDeltaE, Direction direction) {                     
	fSecondDeltaE = &secondDeltaE;
	fSecondDirection = direction;
	fSecondPosition.SetXYZ(0., 0., 0.);
}

void HitSim::SetPad(ParticleMC& pad) {
	fPad = &pad;
}

/******************************************************************
 * Barrel hit position 
 * if gaslength in Settingsfile > 0: assume single sided strip detector
 * otherwise assume double-sided strip
 * 
 * returns TVector3(0,0,0) if:
 *                        -> two not neighboring strips are hit
 *                        -> broken strip is hit (not possible in current simulation)
 *
 *****************************************************************/
TVector3 HitSim::FirstPosition(bool doubleSidedFirstLayer, bool smear) {
	if(fFirstDeltaE == nullptr) return TVector3(0., 0., 0.);
	if(fFirstPosition != TVector3(0., 0., 0.)) {
		return fFirstPosition;
	}

	// variables for final hit position
	double x,y,z;

	// quadrant
	int quadr = fFirstDeltaE->GetID();
	//std::cout<<" \n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! first layer no: "<< quadr << std::endl;

	// strip number = perpendicular to beam direction
	double strip = 0;

	// two neighboring strips hit: calculate mean strip number
	if(fFirstDeltaE->GetNeighborStrip()) { 
		for(unsigned int i = 0; i < fFirstDeltaE->GetStripNr().size(); i++) {
			strip += fFirstDeltaE->GetStripNr()[i];
		}
		strip /= fFirstDeltaE->GetStripNr().size();
	} else if(fFirstDeltaE->GetStripNr().size() > 1) { 
		// two not neighbooring strips: ignore them
		//std::cerr<<"found "<<fFirstDeltaE->GetStripNr().size()<<" strips "<<fFirstDeltaE->GetStripNr()[0]<<" and "<<fFirstDeltaE->GetStripNr()[1]<<" but not neighboring!"<<std::endl; 
		return TVector3(0,0,0);
	} else if(fFirstDeltaE->GetStripNr().size() == 1) { 
		// one hit only
		strip = fFirstDeltaE->GetStripNr()[0];
	} else { 
		// no hit
		std::cerr<<"can not find any hit "<<std::endl;
		return TVector3(0,0,0);
	}

	// smear strip position
	if(smear) {
	  strip += (1-fRand->Uniform());  //uniform (0,1] - > 1-uniform [0,1)
	} else { // use mean strip position
	  strip += 0.5;
	}
	
	// running on a gas target? --> use single-sided single strip detector
	if(!doubleSidedFirstLayer) { 
	    // we have no information about phi, so we use the phi of the second layer
	    SecondPosition(smear); //call to make sure the second position has been calculated

	    // strip 0 closest to target
	    if(fFirstDirection == kForward) {//std::cout<<" \n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 1. layer quadr forward: "<< quadr << std::endl;
	      z = fTRexSett->GetLayerPositionZ(0)[quadr] - fTRexSett->GetLayerDimensionsY(0)/2. + strip*fSett->GetTISTARStripWidthY(0);
	    } else { // backward
			//std::cout<<" \n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 1. layer quadr backward: "<< quadr << std::endl;
	      z = fTRexSett->GetLayerPositionZ(0)[quadr+2] + fTRexSett->GetLayerDimensionsY(0)/2. - strip*fSett->GetTISTARStripWidthY(0);
	    }
	    //quadr   0 top   1 lef   2 bot   3 rig
	    //x       +pos    +dtb    -pos    -dtb
	    //y       +dtb    -pos    -dtb    +pos
	    
	    switch(quadr) {
	    case 0:
	      x = fSecondPosition.X()*fTRexSett->GetLayerDistToBeam(0)[quadr]/fTRexSett->GetLayerDistToBeam(1)[quadr];
	      y = fTRexSett->GetLayerDistToBeam(0)[quadr];
	      break; 
	      
	      // exchange case1 and case 2 because it is not any more box shape, just top (0) and bottom (0) Leila 
	      
	    case 1:
	      /*x = fTRexSett->GetLayerDistToBeam(0)[quadr];
	      y = fSecondPosition.Y()*fTRexSett->GetLayerDistToBeam(0)[quadr]/fTRexSett->GetLayerDistToBeam(1)[quadr]; original*/
	      
	      x = fSecondPosition.X()*fTRexSett->GetLayerDistToBeam(0)[quadr]/fTRexSett->GetLayerDistToBeam(1)[quadr];
	      y = -fTRexSett->GetLayerDistToBeam(0)[quadr]; // changed by Leila	      	      
	      break;
	      
	    case 2:
	      /*x = fSecondPosition.X()*fTRexSett->GetLayerDistToBeam(0)[quadr]/fTRexSett->GetLayerDistToBeam(1)[quadr];
	      y = -fTRexSett->GetLayerDistToBeam(0)[quadr]; original */
	      
	      x = fTRexSett->GetLayerDistToBeam(0)[quadr];
	      y = fSecondPosition.Y()*fTRexSett->GetLayerDistToBeam(0)[quadr]/fTRexSett->GetLayerDistToBeam(1)[quadr];// changed by Leila
	      break;
	      
	    case 3:
	      x = -fTRexSett->GetLayerDistToBeam(0)[quadr];
	      y = fSecondPosition.Y()*fTRexSett->GetLayerDistToBeam(0)[quadr]/fTRexSett->GetLayerDistToBeam(1)[quadr];
	      break;
	      
	    default:
	      break;
	    }
	  } else {
		  // running on a solid target --> double-sided strip detector
	    // "ring" number = strips parallel to beam direction
	    double ring = 0;
	    // two neighboring rings hit: calculate mean strip number
	    if(fFirstDeltaE->GetNeighborRing()) {
			 for(unsigned int i = 0; i < fFirstDeltaE->GetRingNr().size(); i++) { 
				 ring += fFirstDeltaE->GetRingNr()[i];
			 }
			 ring /= fFirstDeltaE->GetRingNr().size();
		 } else if(fFirstDeltaE->GetRingNr().size() > 1) { 
			 // two not neighbooring rings: ignore them
	      // std::cerr<<"first : found "<<fFirstDeltaE->GetRingNr().size()<<" \"rings\" "<<fFirstDeltaE->GetRingNr()[0]<<" and "<<fFirstDeltaE->GetRingNr()[1]<<" but not neighboring!"<<std::endl; 
			 return TVector3(0,0,0);
		 } else if(fFirstDeltaE->GetRingNr().size() == 1) { 
			 // one hit only
			 ring = fFirstDeltaE->GetRingNr()[0];
		 } else { 
			 // no hit
			 std::cerr<<"can not find any hit "<<std::endl;
			 return TVector3(0,0,0);
		 }

		 // smear ring position
		 if(smear) {
			 ring+=(1-fRand->Uniform());  //uniform (0,1] - > 1-uniform [0,1)
		 } else { // use mean strip position
			 ring+=0.5;
		 }


		 // strip 0 closest to target
		 if(fFirstDirection == kForward) {
			 z = fTRexSett->GetLayerPositionZ(0)[quadr] - fTRexSett->GetLayerDimensionsY(0)/2. + strip*fSett->GetTISTARStripWidthY(0);
			 
            //change strip # so that the range isn't 0 - (n-1), but -n/2 - n/2
			 ring -= fTRexSett->GetLayerDimensionsX(0)/fSett->GetTISTARStripWidthX(0)/2.;

			 //quadr   0 top   1 lef   2 bot   3 rig
			 //x       +pos    +dtb    -pos    -dtb
			 //y       +dtb    -pos    -dtb    +pos
			 
			 // exchange case1 and case 2 because it is not any more box shape, just top (0) and bottom (0) Leila

			 switch(quadr) {
				 case 0:
					 x = -ring*fSett->GetTISTARStripWidthX(0);
					 y = fTRexSett->GetLayerDistToBeam(0)[quadr];
					 break;
					 
				 case 1:
					 /*x = fTRexSett->GetLayerDistToBeam(0)[quadr];
					 y = ring*fTRexSett->GetFBarrelDeltaESingleStripWidth(); original */
					 
					 x = ring*fSett->GetTISTARStripWidthX(0);
					 y = -fTRexSett->GetLayerDistToBeam(0)[quadr];
					 break;
					 
				 case 2:
					 /*x = ring*fTRexSett->GetFBarrelDeltaESingleStripWidth();
					 y = -fTRexSett->GetLayerDistToBeam(0)[quadr]; original */
					 
					 x = fTRexSett->GetLayerDistToBeam(0)[quadr];
					 y = ring*fSett->GetTISTARStripWidthX(0);
					 break;
					 
				 case 3:
					 x = -fTRexSett->GetLayerDistToBeam(0)[quadr];
					 y = -ring*fSett->GetTISTARStripWidthX(0);
					 break;
				 default:
					 break;
			 }
		 } else { // backward
			 z = fTRexSett->GetLayerPositionZ(0)[quadr+2] + fTRexSett->GetLayerDimensionsY(0)/2. - strip*fSett->GetTISTARStripWidthY(0);
			 //change ring # so that the range isn't 0 - (n-1), but -n/2 - n/2
			 ring -= fTRexSett->GetLayerDimensionsX(0)/fSett->GetTISTARStripWidthX(0)/2.;


			 switch(quadr) {
				 case 0:
					 x = -ring*fSett->GetTISTARStripWidthX(0);
					 y = fTRexSett->GetLayerDistToBeam(0)[quadr+2];
					 break;
				 case 1:
					 /*x = fTRexSett->GetLayerDistToBeam(0)[quadr+2];
					 y = ring*fTRexSett->GetBBarrelDeltaESingleStripWidth(); original */
					 
					 x = ring*fSett->GetTISTARStripWidthX(0);
					 y = -fTRexSett->GetLayerDistToBeam(0)[quadr+2]; // changed by Leila
					 break;
					 
				 case 2:
					 /*x = ring*fTRexSett->GetBBarrelDeltaESingleStripWidth();
					 y = -fTRexSett->GetLayerDistToBeam(0)[quadr+2]; original */
					 
					 x = fTRexSett->GetLayerDistToBeam(0)[quadr+2];
					 y = ring*fSett->GetTISTARStripWidthX(0); // changed by Leila
					 break;
				 case 3:
					 x = -fTRexSett->GetLayerDistToBeam(0)[quadr+2];
					 y = -ring*fSett->GetTISTARStripWidthX(0);
					 break;
				 default:
					 break;
			 }
		 }
	  }

	// final position in lab frame
	fFirstPosition.SetXYZ(x, y, z);

	return fFirstPosition;
}

//#B.Wach

/******************************************************************
 * Second Barrel hit position (assuming double sided strip detector)
 * 
 * returns TVector3(0,0,0) if:
 *                        -> two not neighboring strips are hit
 *                        -> broken strip is hit (not possible in current simulation)
 *
 *****************************************************************/ 
TVector3 HitSim::SecondPosition(bool smear) {                             
	if(fSecondDeltaE == nullptr) return TVector3(0., 0., 0.);
	if(fSecondPosition != TVector3(0., 0., 0.)) {
		return fSecondPosition;
	}

	// variables for final hit position
	double x,y,z;

	// quadrant
	int quadr = fSecondDeltaE->GetID();  
	//std::cout<<" \n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! second layer no: "<< quadr << std::endl;                                                 

	// strip number = strips perpendicular to beam direction
	double strip = 0;

	// two neighboring strips hit: calculate mean strip number
	if(fSecondDeltaE->GetNeighborStrip()) {
		for(unsigned int i = 0; i < fSecondDeltaE->GetStripNr().size(); i++) { 
			strip += fSecondDeltaE->GetStripNr()[i];
		}
		strip /= fSecondDeltaE->GetStripNr().size();
	} else if(fSecondDeltaE->GetStripNr().size() > 1) { 
		// two not neighbooring strips: ignore them
	  //	std::cerr<<"second : found "<<fSecondDeltaE->GetStripNr().size()<<" strips "<<fSecondDeltaE->GetStripNr()[0]<<" and "<<fSecondDeltaE->GetStripNr()[1]<<" but not neighboring!"<<std::endl; 
		return TVector3(0,0,0);
	} else if(fSecondDeltaE->GetStripNr().size() == 1) { 
		// one hit only
		strip = fSecondDeltaE->GetStripNr()[0];
	} else { 
		// no hit
		std::cerr<<"can not find any hit "<<std::endl;
		return TVector3(0,0,0);
	}

	// "ring" number = strips parallel to beam direction
	double ring = 0;

	// two neighboring rings hit: calculate mean strip number
	if(fSecondDeltaE->GetNeighborRing()) {
		for(unsigned int i = 0; i < fSecondDeltaE->GetRingNr().size(); i++) { 
			ring += fSecondDeltaE->GetRingNr()[i];
		}
		ring /= fSecondDeltaE->GetRingNr().size();
	} else if(fSecondDeltaE->GetRingNr().size() > 1) { 
		// two not neighbooring rings: ignore them
		//std::cerr<<"second : found "<<fSecondDeltaE->GetRingNr().size()<<" \"rings\" "<<fSecondDeltaE->GetRingNr()[0]<<" and "<<fSecondDeltaE->GetRingNr()[1]<<" but not neighboring!"<<std::endl; 
		return TVector3(0,0,0);
	} else if(fSecondDeltaE->GetRingNr().size() == 1) { 
		// one hit only
		ring = fSecondDeltaE->GetRingNr()[0];
	} else { 
		// no hit
		std::cerr<<"can not find any hit "<<std::endl;
		return TVector3(0,0,0);
	}

	// smear ring position
	if(smear) {
		ring+=(1-fRand->Uniform());  //uniform (0,1] - > 1-uniform [0,1)
	} else { // use mean strip position
		ring+=0.5;
	}

	// strip 0 closest to target
	if(fSecondDirection == kForward) {//std::cout<<" \n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 2. layer strips/rings quadr forward: "<< quadr << std::endl;
		z = fTRexSett->GetLayerPositionZ(1)[quadr] - fTRexSett->GetLayerDimensionsY(1)/2. + strip*fSett->GetTISTARStripWidthY(1);
		//change ring # so that the range isn't 0 - (n-1), but -n/2 - n/2
		ring -= fTRexSett->GetLayerDimensionsX(1)/fSett->GetTISTARStripWidthX(1)/2.;

		//quadr   0 top   1 lef   2 bot   3 rig
		//x       +pos    +dtb    -pos    -dtb
		//y       +dtb    -pos    -dtb    +pos
		
		 // exchange case1 and case 2 because it is not any more box shape, just top (0) and bottom (0) Leila

		switch(quadr) {
			case 0:
				x = -ring*fSett->GetTISTARStripWidthX(1); // why negative????
				y = fTRexSett->GetLayerDistToBeam(1)[quadr];
				break;
				
			case 1:
				/*x = fTRexSett->GetLayerDistToBeam(1)[quadr];
				y = ring*fTRexSett->GetSecondFBarrelDeltaESingleStripWidth(); original */
				
				x = ring*fSett->GetTISTARStripWidthX(1);
				y = -fTRexSett->GetLayerDistToBeam(1)[quadr]; // changed by Leila
				break;
				
			case 2:
				/*x = ring*fTRexSett->GetSecondFBarrelDeltaESingleStripWidth();
				y = -fTRexSett->GetLayerDistToBeam(1)[quadr]; original*/
				
				x = fTRexSett->GetLayerDistToBeam(1)[quadr];
				y = ring*fSett->GetTISTARStripWidthX(1); // changed by Leila
				break;
				
			case 3:
				x = -fTRexSett->GetLayerDistToBeam(1)[quadr];
				y = -ring*fSett->GetTISTARStripWidthX(1);
				break;
			default:
				break;
		}
	} else { // backward
		//std::cout<<" \n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 2. layer strips/rings quadr backward: "<< quadr << std::endl;
		z = fTRexSett->GetLayerPositionZ(0)[quadr+2] + fTRexSett->GetLayerDimensionsY(1)/2. - strip*fSett->GetTISTARStripWidthY(1);
		//change ring # so that the range isn't 0 - (n-1), but -n/2 - n/2
		ring -= fTRexSett->GetLayerDimensionsX(1)/fSett->GetTISTARStripWidthX(1)/2.;


		switch(quadr) {
			case 0:
				x = -ring*fSett->GetTISTARStripWidthX(1);
				y = fTRexSett->GetLayerDistToBeam(1)[quadr];
				break;
				
			case 1:
				/*x = fTRexSett->GetLayerDistToBeam(1)[quadr];
				y = ring*fTRexSett->GetSecondBBarrelDeltaESingleStripWidth(); original */
				
				x = ring*fSett->GetTISTARStripWidthX(1);
				y = -fTRexSett->GetLayerDistToBeam(1)[quadr]; // changed by Leila
				break;
			case 2:
				/*x = ring*fTRexSett->GetSecondBBarrelDeltaESingleStripWidth();
				y = -fTRexSett->GetLayerDistToBeam(1)[quadr]; original */ 
				
				x = fTRexSett->GetLayerDistToBeam(1)[quadr];
				y = ring*fSett->GetTISTARStripWidthX(1); // changed by Leila
				break;
				
			case 3:
				x = -fTRexSett->GetLayerDistToBeam(1)[quadr];
				y = -ring*fSett->GetTISTARStripWidthX(1);
				break;
			default:
				break;
		}
	}

	// final position in lab frame
	fSecondPosition.SetXYZ(x, y, z);

	return fSecondPosition;
}

double HitSim::GetFirstDeltaEEnergy(bool verbose) {
	// loop over all strips and return their added energies
	if(fFirstEnergy == -99.) {
		if(fFirstDeltaE == nullptr) return 0.;
		fFirstEnergy = 0.;
		if(verbose) std::cout<<"first: "<<fFirstEnergy<<" -> ";
		for(auto en : fFirstDeltaE->GetStripEnergy()) {
			fFirstEnergy += en;
			if(verbose) std::cout<<fFirstEnergy<<" -> ";
		}
		if(verbose) std::cout<<fFirstEnergy<<std::endl;
	}

	return fFirstEnergy;
}

double HitSim::GetSecondDeltaEEnergy(bool verbose) {
	// loop over all strips and return their added energies
	if(fSecondEnergy == -99.) {
		if(fSecondDeltaE == nullptr) return 0.;
		fSecondEnergy = 0.;
		if(verbose) std::cout<<"second: "<<fSecondEnergy<<" -> ";
		for(auto en : fSecondDeltaE->GetStripEnergy()) {
			fSecondEnergy += en;
			if(verbose) std::cout<<fSecondEnergy<<" -> ";
		}
		if(verbose) std::cout<<fSecondEnergy<<std::endl;
	}

	return fSecondEnergy;
}

double HitSim::GetPadEnergy() {
	if(fPad != nullptr) return fPad->GetEdet();
	return 0.;
}
