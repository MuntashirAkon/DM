//
// This file contains different constant definitions.
//

#pragma once

//
// Disable warning generated due to using depricated (unsafe) functions such as
// freopen, sprintf etc.
//
#pragma warning(disable:4996)

//
// Maximum number of genes that can be used in experiment
//
#define NUM_GENE 20

//
// Maximum length of a species name Number of genes used in experiment
//

#define MAX_SPECIES_NAME_LEN 41

//
// Maximum possible length of an absent word. To be safe this should be set to
// a length that is greater than the longest gene used in experiment.
// Since we use a gene with maximum length of 105, and concatenate reverse
// complement, so the size can be 210. So, choosing limit accordingly
//
#define MAX_WORD_LEN 256

//
// Maximum length of a file location
//
#define MAX_PATH 256

#define MAW     1   // Minimum abset word
#define RAW     2   // Relative absent word

//
// Defined constants for different MAW based similarity matrices
//
#define MAW_LWI_SDIFF       1  // Length weighted index of symmetric difference
#define MAW_LWI_INTERSECT   2  // Length weighted index of intersection
#define MAW_GCC_SDIFF       3  // GC content of symmetric difference
#define MAW_GCC_INTERSECT   4  // GC content of intersection
#define MAW_JD              5  // Jaccard Distance
#define MAW_TVD             6  // Jaccard Distance

//
// Defined constants for different RAW based similarity matrices
//
#define RAW_LWI            1   // Length weighted index
#define RAW_GCC            2   // GC content
