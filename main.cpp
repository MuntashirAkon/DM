//
// This file contains the main function of the project
//
#include <iostream>		// <Mun>
#include <string>		// <Mun>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cassert>

#include "Set.h"
#include "Tests_1.cpp"
#include "Constants.cpp"
#include "Globals.cpp"
#include "DllMain.cpp"

using namespace std;

//
// Debug mode
//
#define DEBUG true

//
//	<Mun> CLI Return Values
//
#define SUCCESS 0
#define FAILURE 1

//
// Input files
//
#define SPECIES_FULL "SpeciesFull.txt"


//
// Output files
//
#define SPECIES_RELATION_TXT  "SpeciesRelation.txt"
#define SPECIES_RELATION_JSON "SpeciesRelation.json"
#define DISTANCE_MATRIX       "DistanceMatrix.txt"

int _gLargestSpeciesLength = 0;

extern int  g_numGenes;
extern char g_strSpeciesShortName[NUM_GENE][MAX_SPECIES_NAME_LEN];
extern char g_strSpeciesFullName[NUM_GENE][MAX_SPECIES_NAME_LEN];
extern char g_strDataDir[MAX_PATH];

extern int Initialize(
    char **geneFullNames,
    char **geneShortNames,
    int nGenes,
    char *dataDir
    );

extern int getDiffMatrix(
    double diffMatrix[][NUM_GENE],
    int absWordType,
    int diffIndex
    );

extern int getRanks(
    int rank[][NUM_GENE],
    int absWordType,
    int diffIndex
    );

void PrintSpeciesRelations(int ranks[][NUM_GENE], ofstream *text, ofstream *json){
    int i, j;
    
    if(!text->is_open() && !json->is_open()) return;
    
    *json << "{" << endl;
    for (i = 0; i < g_numGenes; i++){
        *text << setw(_gLargestSpeciesLength) << g_strSpeciesFullName[i] << ":";
        *json << "    \"" << g_strSpeciesFullName[i] << "\": [" << endl;
        for (j = 0; j < g_numGenes; j++){
            *text << " -> " << setw(_gLargestSpeciesLength) << g_strSpeciesFullName[ranks[i][j]] << ((j != g_numGenes - 1) ? " ": "\n");
            *json << "        \"" << g_strSpeciesFullName[ranks[i][j]] << "\"" << ((j != g_numGenes - 1) ?  ",": "") << endl;
        }
        *json << "    ]" << ((i != g_numGenes - 1) ? ", " : "") << endl;
    }
    *json << "}" << endl;
    text->close();
    json->close();
}

//
// Java array format, so that it can be tested directly using Match7
//
void PrintDiffMatrix(double diffMatrix[][NUM_GENE], FILE* f1, FILE* f2){
    int i, j;

    for (i = 0; i < g_numGenes; i++)
    {
        cout << "{ ";

        for (j = 0; j < i; j++)
        {
            printf("%5.2lf", diffMatrix[i][j]);
            fprintf(f2,"%5.2lf", diffMatrix[i][j]);
            fprintf(f1, "%5.4lf\n", diffMatrix[i][j]);
            printf("%s", (j == i - 1) ? "" : ", ");
            fprintf(f2,"%s", (j == i - 1) ? "" : ", ");
        }

        cout << " }" << ((i == g_numGenes - 1) ? "" : ",") << endl;
        fprintf(f2," %s\n", (i == g_numGenes - 1) ? "" : ",");
    }

}

void ShowHelp(){
	cout << "Usage: dm Type DiffIndex source [target]" << endl;
	cout << "Where:" << endl;
	cout << "Type\t\tMAW|RAW" << endl;
	cout << "DiffIndex\tMAW_LWI_SDIFF|MAW_LWI_INTERSECT|"
		"MAW_GCC_SDIFF|MAW_GCC_INTERSECT|MAW_JD|MAW_TVD|"
			"RAW_LWI|RAW_GCC" << endl;
}

int main(int argc, char *argv[]){
	//
	// <Mun> First check if necessary arguments are supplied
	//
	if(argc < 4){
		cout << "Insufficient argument(s) supplied!" << endl;
		ShowHelp();
		return FAILURE; 	// Not enough arguments to work with, so exit FAILURE
	}

	//
	// <Mun> Change the argv to readable values
	//
	string word_type  = argv[1];  // Absent word type
	string diff_index = argv[2];  // Diff Index
	string input_dir  = argv[3];
	string output_dir = (argc > 5) ? argv[4] : input_dir; // Set input dir as output dir if 4th argument isn't set

    //
    // Change the following 2 variables to run with different
    // absent word type (RAW) and diff indices.
    //
    int absWordType = (word_type == "RAW") ? RAW : MAW;
    int diffIndex = MAW_TVD;    // Default for RAW: RAW_GCC
	
	// <Mun>
	if(absWordType == RAW){
		if(diff_index == "RAW_LWI") diffIndex = RAW_LWI;
		else diffIndex = RAW_GCC;	// Default for RAW
	}else{	// MAW
		if(diff_index == "MAW_LWI_SDIFF") diffIndex = MAW_LWI_SDIFF;
		else if(diff_index == "MAW_LWI_INTERSECT") diffIndex = MAW_LWI_INTERSECT;
		else if(diff_index == "MAW_GCC_SDIFF") diffIndex = MAW_GCC_SDIFF;
		else if(diff_index == "MAW_GCC_INTERSECT") diffIndex = MAW_GCC_INTERSECT;
		else if(diff_index == "MAW_JD") diffIndex = MAW_JD;
		// else diffIndex = MAW_TVD; // No need, already set
	}
	
    if(DEBUG) cout << "Data Dir   : " << input_dir << endl;

    //
    // Change the following variables for running a different set of species/gene sequences
    // FIXME: Transform the bellow into pure C++ code
    //

    
    //
    // Short name code used for each species.
    char* strSpeciesShortName[MAX_SPECIES_NAME_LEN];

    //
    // Full name of each species.
    char* strSpeciesFullName[MAX_SPECIES_NAME_LEN];

    if(DEBUG) cout << "SpeciesFull: " << input_dir << '/' << SPECIES_FULL << endl;

    int nGenes = 0;
    ifstream speciesFull(input_dir + '/' + SPECIES_FULL);
    if(!speciesFull.is_open()){
        cerr << "The source directory doesn't contain \"SpeciesFull.txt\"!" << endl;
        return FAILURE;
    }
    // Get line has problems with carriage return
    for (string line; getline(speciesFull, line); ++nGenes){
        if(_gLargestSpeciesLength < line.length()) _gLargestSpeciesLength = line.length();
        strcpy(g_strSpeciesFullName[nGenes], line.c_str());
    }
    g_numGenes = nGenes;
	
    double diffMatrix[NUM_GENE][NUM_GENE] = {0};
    int rank[NUM_GENE][NUM_GENE] = {0};

    if(DEBUG) cout << "Output file: " << output_dir << '/' << DISTANCE_MATRIX << endl;

    FILE *f3 = fopen((output_dir + '/' + DISTANCE_MATRIX).c_str(), "w+");

    // Formatted for Match7
    FILE *f4 = fopen((output_dir + '/' + "Output.txt").c_str(), "w+");
    strcpy(g_strDataDir, input_dir.c_str());

    //Initialize(strSpeciesFullName, strSpeciesShortName, nGenes, g_strDataDir);

    getDiffMatrix(diffMatrix, absWordType, diffIndex);
    PrintDiffMatrix(diffMatrix, f3, f4);

    getRanks(rank, absWordType, diffIndex);
    PrintSpeciesRelations(rank, new ofstream(output_dir + '/' + SPECIES_RELATION_TXT), new ofstream(output_dir + '/' + SPECIES_RELATION_JSON));

    return SUCCESS;
}
