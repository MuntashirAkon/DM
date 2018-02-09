//
// This file contains different global variable definitions.
//

#pragma once

#include "Constants.cpp"

//
// Total number of genes/species
//
int  g_numGenes = 0;

//
// Short name code used for each species.
// (Used in generating the filename for Raw/Maw sets for respective species)
//
char g_strSpeciesShortName[NUM_GENE][MAX_SPECIES_NAME_LEN];

//
// Full name of each species.
//
char g_strSpeciesFullName[NUM_GENE][MAX_SPECIES_NAME_LEN];

//
// Location of the folder that contains necessary data files
//
char g_strDataDir[MAX_PATH] = "";

