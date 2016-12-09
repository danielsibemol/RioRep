#include <stdio.h>
#include <math.h>
#include <string>
#include <TComplex.h>
#include <TEnv.h>

using namespace std;

void Read_Number_of_Samples(string input_txt_file_name, int &number_of_samples){

	TEnv *input_txt_file = new TEnv(input_txt_file_name.c_str());

	// This function reads the number of samples to be used in the fit

	// Checks if the input .txt file was loaded correctly

	if (input_txt_file->IsZombie()) {
		std::cout << "ReadParameters - ERROR: File not found - " << input_txt_file_name << std::endl;
		exit(-1);
	}

    // Reads the number of samples
    
    number_of_samples = input_txt_file->GetValue("number_of_samples", 1);
    
    cout << "number of samples: " << number_of_samples << endl;
    cout << endl;
    cout << "------------------------------------------------------------------------" << endl;
    cout << endl;
}
