//
// Created by Muntashir Al-Islam on 3/13/18.
//

#ifndef DM_PROCESS_H
#define DM_PROCESS_H

#include <vector>
#include <string>
#include <cassert>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Set.h"

using namespace std;

//
// Input file
//
#define SPECIES_FULL    "SpeciesFull.txt"

//
// Output files
//
#define SPECIES_RELATION_TXT    "SpeciesRelation.txt"
#define SPECIES_RELATION_JSON   "SpeciesRelation.json"
#define DISTANCE_MATRIX         "DistanceMatrix.txt"
#define DISTANCE_MATRIX_PHYLIP  "DistanceMatrix.dist"

class Process{
    int total_genes;
    vector<string> species;
    int largest_species_len;
    int aw_type, di_type;
    string input_dir, output_dir;
    vector< vector<double> > diffMatrix;
    vector< vector<int> > rank;

public:
    Process(int absent_word_type, int dissimilarity_index_type, const string& input_dir, const string& output_dir);
    void printSortedSpeciesRelations();
    void printDiffMatrix();
    /**
     * Absent Word types
     */
    enum AW{RAW, MAW};
    /**
     * Defined constants for different MAW based similarity matrices
     */
    enum MAW{MAW_LWI_SDIFF, MAW_LWI_INTERSECT, MAW_GCC_SDIFF, MAW_GCC_INTERSECT, MAW_JD, MAW_TVD};
    /**
     * Defined constants for different RAW based similarity matrices
     */
    enum RAW{RAW_LWI, RAW_GCC};

private:
    /**
     * Gets a list of species from SpeciesFull.txt
     */
    void getSpecies();
    /**
     * Gets the different matrix
     * @return
     */
    int getDiffMatrix();
    /**
     * Ranks the different species relative to one species based on the difference matrix.
     * It does this analysis for each species and generates a 2D rank table.
     */
    void getRanks();
    int runMaw();
    double calculateTVD(vector<pair<int, double> > rgA, vector<pair<int, double> > rgB);
    int runRaw();
    double GetRawBasedDiff(int i, int j, int diffIndex);
};

#endif //DM_PROCESS_H
