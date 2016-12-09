#include <stdio.h>
#include <iostream>
#include <fstream>
#define MAIN_FILE
#include "shared.h"
#include "ToyMCGenerator.h"


using namespace std;

int main(int argc,char *argv[]){
	cout << endl;

	cout << "###################################################" << endl;
	cout << "#                                                 #" << endl;
	cout << "#             Welcome to ToyMCGenerator           #" << endl;
	cout << "#                                                 #" << endl;
	cout << "#                                                 #" << endl;
	cout << "#       Authors: Sandra Amato, Rafael Aoude       #" << endl;
	cout << "#   Carla Gobel, Josue Molina, Érica Polycarpo,   #" << endl;
	cout << "#         Alberto Reis, Danielle Tostes,          #" << endl;   
	cout << "#                and Daniel Vieira                #" << endl;
	cout << "#                                                 #" << endl;
	cout << "#            contact: dvieira@if.ufrj.br          #" << endl;
	cout << "#                                                 #" << endl;
	cout << "#                    Have fun!                    #" << endl;
	cout << "#                                                 #" << endl;
	cout << "###################################################" << endl;

	cout << endl;
	cout << "GLOBAL_TOY " << GLOBAL_TOY << endl;

	//////////////////////////////
	if ( argc !=2 ) {
		cout << "main - ERROR: Cannot open input text file, no argument found" << endl;
		return EXIT_FAILURE;
	}
	//////////////////////////////

	string string=argv[1];
	cout << "Reading input file name:" << string<< endl;
	ToyMCGenerator(string);

	return 0;
}

