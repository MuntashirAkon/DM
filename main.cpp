//
// This file contains the main function of the project
//
#include <iostream>		// <Mun>
#include <string>		// <Mun>
#include <cstdio>
#include <cassert>

#include "Set_1.cpp"
#include "Tests_1.cpp"
#include "Constants.cpp"
#include "Globals.cpp"
#include "DllMain.cpp"

using namespace std;

//
//	<Mun> CLI Return Values
//
#define SUCCESS 0
#define FAILURE 1

extern  int g_numGenes;
extern  char g_strSpeciesShortName[NUM_GENE][MAX_SPECIES_NAME_LEN];
extern  char g_strSpeciesFullName[NUM_GENE][MAX_SPECIES_NAME_LEN];
extern  char g_strDataDir[MAX_PATH];

extern  int Initialize(
    char **geneFullNames,
    char **geneShortNames,
    int nGenes,
    char *dataDir
    );

extern  int getDiffMatrix(
    double diffMatrix[][NUM_GENE],
    int absWordType,
    int diffIndex
    );

extern  int getRanks(
    int rank[][NUM_GENE],
    int absWordType,
    int diffIndex
    );

void PrintSpeciesRelations(int ranks[][NUM_GENE], FILE* text, FILE* json)
{
    int i, j;

    fprintf(json, "{\n");
    for (i = 0; i < g_numGenes; i++)
    {
        fprintf(text, "%7s:", g_strSpeciesFullName[i]);
        fprintf(json, "    \"%s\": [\n", g_strSpeciesFullName[i]);
        for (j = 0; j < g_numGenes; j++)
        {
            fprintf(text, " -> %s%c", g_strSpeciesFullName[ranks[i][j]], (j == g_numGenes-1 ? '\n' : ' '));
            fprintf(json, "        \"%s\"", g_strSpeciesFullName[ranks[i][j]]);
            if(j != g_numGenes - 1) fprintf(json, ",");
            fprintf(json, "\n");
        }
        fprintf(json, "    ]");
        if(i != g_numGenes - 1) fprintf(json, ", ");
        fprintf(json, "\n");
    }
    fprintf(json, "}\n");
}

void PrintDiffMatrix(double diffMatrix[][NUM_GENE], FILE* f1, FILE* f2)
{
    int i, j;

    for (i = 0; i < g_numGenes; i++)
    {
        printf("{ ");

        for (j = 0; j < i; j++)
        {
            printf("%5.2lf", diffMatrix[i][j]);
            fprintf(f2,"%5.2lf", diffMatrix[i][j]);
            fprintf(f1, "%5.4lf\n", diffMatrix[i][j]);
            printf("%s", (j == i - 1) ? "" : ", ");
            fprintf(f2,"%s", (j == i - 1) ? "" : ", ");
        }

        printf(" }%s\n", (i == g_numGenes - 1) ? "" : ",");
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

int main(int argc, char *argv[])
{
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
	string word_type 	= argv[1];
	string diff_index 	= argv[2];
	string input_dir	= argv[3];
	string output_dir	= (argc > 5) ? argv[4] : input_dir; // Set input dir as output dir if 4th argument isn't set
	
	//
	// <Mun> Add trailing slash to the directories
	//
        //if(input_dir.back() != '/')
	input_dir   += "/";
	//if(output_dir.back() != '/')
	output_dir += "/";
		
    //
    // Change the following 2 variables to run with different
    // absent word type (RAW) and diff indices.
    //
    int absWordType = MAW;		// Default: RAW
    int diffIndex = MAW_TVD;	// Default: RAW_LWI
    // char url[300];				// <Mun> To store enough data

	
	if(word_type == "RAW") absWordType = RAW;	// <Mun> No need for else, since MAW is already set

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
	
    //
    // Change the following variable to change the data folder location
    //
    char *strDataDir;
	
    strDataDir = (char *) malloc( input_dir.length() + 1 );
    strcpy(strDataDir, input_dir.c_str());
    cout << "Data Dir   : " << strDataDir << endl;

    //
    // Change the following variables for running a different set of species/gene sequences
    //

    int nGenes = 0;

    //
    // Short name code used for each species.
    char* strSpeciesShortName[MAX_SPECIES_NAME_LEN];

    //
    // Full name of each species.
    char* strSpeciesFullName[MAX_SPECIES_NAME_LEN];

    char line[MAX_SPECIES_NAME_LEN + 3];
    cout << "SpeciesFull: " << input_dir << "SpeciesFull.txt\n";

    FILE *f2 = fopen((input_dir + "SpeciesFull.txt").c_str(), "r");
	if(!f2){
		cout << "The source directory doesn't contain \"SpeciesFull.txt\"!" << endl;
		return FAILURE;
	}

    while( fgets(line, sizeof line, f2) != NULL ){
        line[strlen(line)-1] = '\0';
        strcpy(g_strSpeciesFullName[nGenes], line);
        nGenes++;
    }
    g_numGenes = nGenes;
	
    double diffMatrix[NUM_GENE][NUM_GENE] = {0};
    int rank[NUM_GENE][NUM_GENE] = {0};

    printf("%sDistanceMatrix.txt\n", output_dir.c_str());


    FILE *f3 = fopen((input_dir + "DistanceMatrix.txt").c_str(), "w+");

    FILE *f4 = fopen((output_dir + "Output.txt").c_str(), "w+");
    strcpy(g_strDataDir, strDataDir);

    FILE *SR_text = fopen((output_dir + "SpeciesRelation.txt").c_str(), "w+");
    FILE *SR_json = fopen((output_dir + "SpeciesRelation.json").c_str(), "w+");


    Initialize(strSpeciesFullName, strSpeciesShortName, nGenes, strDataDir);

    getDiffMatrix(diffMatrix, absWordType, diffIndex);
    PrintDiffMatrix(diffMatrix, f3, f4);



    getRanks(rank, absWordType, diffIndex);
    PrintSpeciesRelations(rank, SR_text, SR_json);

    return 0;
}
