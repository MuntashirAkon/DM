//
// This file contains necessary code for Relative Absent Word (RAW)
// based difference matrix calculation
//
// The lower triangle of the diff matrix is filled. This is because
// the code used to generate phylogenetic tree needs to lower triangle
//

#include <iostream>
#include <string>
#include <cassert>
#include "Constants.cpp"
#include "Globals.cpp"
#include "Set.h"

using namespace std;

double GetRawBasedDiff(int i, int j, int diffIndex){
    double res = 0.0;
    char strFileName1[MAX_PATH], strFileName2[MAX_PATH];
    Set raw1, raw2;

    assert(i >= 0 && i < sizeof(g_strSpeciesFullName));
    assert(j >= 0 && j < sizeof(g_strSpeciesFullName));

    //
    // Make the file name from the gene indices
    //
    string filename1 = string(g_strDataDir) + '/'
        + string(g_strSpeciesFullName[i]) + '/'
        + string(g_strSpeciesFullName[i]) + '_' + string(g_strSpeciesFullName[j]) + ".raw.txt";
    string filename2 = string(g_strDataDir) + '/'
        + string(g_strSpeciesFullName[j]) + '/'
        + string(g_strSpeciesFullName[j]) + '_' + string(g_strSpeciesFullName[i]) + ".raw.txt";

    raw1 = Set::CreateFromFile(filename1);
    raw2 = Set::CreateFromFile(filename2);

    if (diffIndex == RAW_LWI){
        res = 0.5 * (raw1.LengthWeightedIndex() + raw2.LengthWeightedIndex());
    }else if (diffIndex == RAW_GCC){
        res = 0.5 * (raw1.GCContent() + raw2.GCContent());
    }

    return res;
}

int runRaw(double diffMatrix[][NUM_GENE], int diffIndex){
    int i, j;

    if (diffIndex < RAW_LWI || diffIndex > RAW_GCC)
        return -1;

    for (i = 0; i < g_numGenes; i++)
    {
        diffMatrix[i][i] = 0;
        for (j = 0; j < i; j++)
        {
            diffMatrix[i][j] = diffMatrix[j][i] = GetRawBasedDiff(i, j, diffIndex);
        }
    }

    return 0; // Is there any way to filter it? I think it can be by checking if all the DM are 0.0
}
