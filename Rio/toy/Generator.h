#include <stdio.h>
#include <TMath.h>
#include <TROOT.h>
#include <TRandom1.h>
#include <math.h>
#include <stdio.h>
#include "../src/GenericFunctions.h"

using namespace std;

void Generator(int &final_state, bool &is_gaussian, double &Mass_min, double &Mass_max, bool &is_bkg, double &Bkg_par1, double &Bkg_par2, double &M, double &s12, double &s13, double &s23, double &cos13_12, double &cos12_13, double &cos31_23, TRandom1 * Random, int &npoints, double * points, double * weights){
	
	double s12min, s12max, s13min, s13max, s23min, s23max, s, m1sq, m2sq, m3sq, deltas12, deltas13, random_number_0, random_number_1, m1, m2, m3, mother_mass, mother_width;

	bool boundary = 0;	

	// Checks which is the final state and defines masses acording to it
	
	if (final_state == 0) {
		m1 = mK;
		m2 = mK;
		m3 = mpi;
	}
	else if (final_state == 1) {
		m1 = mK;
		m2 = mpi;
		m3 = mpi;
	}
	else if (final_state == 2) {
		m1 = mpi;
		m2 = mpi;
		m3 = mpi;
	}
	else if (final_state == 3) {
		m1 = mK;
		m2 = mK;
		m3 = mK;
	}
	else{
		cout << "Generator - ERROR: Incorrect final state inserted";
		exit(-1);
	}

	// Generates mass within the mass window with signal or background mass distribution shape

	if (is_gaussian == 1) {
		if (is_bkg) {
			M = Get_Bkg_Mass(Random, Bkg_par1, Bkg_par2, Mass_min, Mass_max);
            mother_mass = D_Mass;
            mother_width = D_width;
		}else {
			M = Get_Gaussian_Mass(Random, mother_mass, mother_width, Mass_min, Mass_max, npoints, points, weights);
		}
	}else {
		M = D_Mass;
	}

	
	
	// Calculation of the limits of the Dalitz plot
	s = M*M;
	s12min = (m1 + m2)*(m1 + m2);
	s12max = (M - m3)*(M - m3);
	s13min = (m1 + m3)*(m1 + m3);
	s13max = (M - m2)*(M - m2);
	s23min = (m2 + m3)*(m2 + m3);
	s23max = (M - m1)*(M - m1);
	
	// Calculation of the square of the masses   
	m1sq = m1*m1;
	m2sq = m2*m2;
	m3sq = m3*m3;
	
	// Calculation of s12 and s13 range
	deltas12 = s12max - s12min;
	deltas13 = s13max - s13min;
	
	while(!boundary){

   	   // Setting values to 0
	   s12 = 0;
	   s13 = 0;
	   s23 = 0;
	   cos13_12 = 0;
	   cos12_13 = 0;
	   cos31_23 = 0;
	
       while(s23 < s23min || s23 > s23max){
			
	      // Generation of s12 an s13
			random_number_0 = Random->Rndm();
			random_number_1 = Random->Rndm();

			s12 = s12min + random_number_0*deltas12;
			s13 = s13min + random_number_1*deltas13;
			
			// Calculation of s23 through the kinematical contraint, checking if s23 is inside the limits
			s23 = s + m1sq + m2sq + m3sq - s12 - s13;
                        
       }
		
	   // Calculation of cos13_12, the cosine of the angle between particles 1 and 3 in (1,2) rest frame
	   cos13_12 = ((s - s12 - m3sq)*(s12 + m1sq - m2sq) + 2*s12*(m1sq + m3sq - s13))/sqrt(lambda(s, s12, m3sq)*lambda(s12, m1sq, m2sq));
		
	   // Calculation of cos12_13, the cosine of the angle between particles 1 and 2 in (1,3) rest frame
	   cos12_13 = ((s - s13 - m2sq)*(s13 + m1sq - m3sq) + 2*s13*(m1sq + m2sq - s12))/sqrt(lambda (s, s13, m2sq)*lambda (s13, m3sq, m1sq));
		
	   // Calculation of cos31_23, the cosine of the angle between particles 3 and 1 in (2,3) rest frame
	   cos31_23 = ((s - s23 - m1sq)*(s23 + m2sq - m3sq) + 2*s23*(m1sq + m3sq - s13))/sqrt(lambda(s, s23, m1sq)*lambda (s23, m2sq, m3sq));
		
	   // Check if the generated point is inside de Dalitz plot boundary. In this case, |cos13| <=1 
	   if(fabs(cos13_12)<1) boundary = 1;
	}
	return;
}

