//
// This file contains the main function of the project
//
#include <iostream>
#include "Process.h"

#define DM_VERSION "2.0.0"

using namespace std;

void ShowHelp(){
	cout << "USAGE: dm Type DiffIndex source [target]\n"
			"WHERE:\n"
			"Type\t\tMAW|RAW\n"
			"DiffIndex\tMAW_LWI_SDIFF|MAW_LWI_INTERSECT|MAW_GCC_SDIFF|MAW_GCC_INTERSECT|MAW_JD|MAW_TVD|RAW_LWI|RAW_GCC";
	cout << endl;
}

int main(int argc, char *argv[]){
    cout << "DM version " << DM_VERSION << endl;
	//
	// First check whether necessary arguments are supplied
	//
	if(argc < 4){
		cout << "Insufficient argument(s) supplied!" << endl;
		ShowHelp();
		return EXIT_FAILURE; 	// Not enough arguments to work with, so exit FAILURE
	}

	//
	// Change the argv to readable values
	//
	string word_type  = argv[1];  // Absent word type
	string diff_index = argv[2];  // Diff Index
	string input_dir  = argv[3];
	string output_dir = (argc > 5) ? argv[4] : input_dir; // Set input dir as output dir if 4th argument isn't set

    //
    // Change the following 2 variables to run with different
    // absent word type (RAW) and diff indices.
    //
    int absWordType = (word_type == "RAW") ? Process::RAW : Process::MAW;
    int diffIndex = Process::MAW_TVD;    // Default for RAW: RAW_GCC
	
	if(absWordType == Process::RAW){
		if(diff_index == "RAW_LWI") diffIndex = Process::RAW_LWI;
		else diffIndex = Process::RAW_GCC;	// Default for RAW
	}else{	// Process::AW::MAW
		if(diff_index == "MAW_LWI_SDIFF") diffIndex = Process::MAW_LWI_SDIFF;
		else if(diff_index == "MAW_LWI_INTERSECT") diffIndex = Process::MAW_LWI_INTERSECT;
		else if(diff_index == "MAW_GCC_SDIFF") diffIndex = Process::MAW_GCC_SDIFF;
		else if(diff_index == "MAW_GCC_INTERSECT") diffIndex = Process::MAW_GCC_INTERSECT;
		else if(diff_index == "MAW_JD") diffIndex = Process::MAW_JD;
		// else diffIndex = MAW_TVD; // No need, already set
	}

    cout << "Input Dir  : " << input_dir << endl;
    cout << "Output Dir : " << output_dir << endl;

    try {
        Process p(absWordType, diffIndex, input_dir, output_dir);
        p.printDiffMatrix();
        p.printSortedSpeciesRelations();
    }catch(exception){ // NOLINT
        return EXIT_FAILURE;
    };
    return EXIT_SUCCESS;
}
