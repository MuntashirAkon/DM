//
// This file contains different global variable definitions.
//

#pragma once

#include "Constants.h"

//
// Total number of genes/species
//
extern __declspec(selectany) int g_numGenes;

//
// Short name code used for each species.
// (Used in generating the filename for Raw/Maw sets for respective species)
//
extern __declspec(selectany) char g_strSpeciesShortName[NUM_GENE][MAX_SPECIES_NAME_LEN];

//
// Full name of each species.
//
extern __declspec(selectany) char g_strSpeciesFullName[NUM_GENE][MAX_SPECIES_NAME_LEN];

//
// Location of the folder that contains necessary data files
//
extern __declspec(selectany) char g_strDataDir[MAX_PATH];
