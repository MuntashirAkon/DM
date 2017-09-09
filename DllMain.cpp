//
// This file defines entry point functions of the MAW DLL.
//
// Initialize() must be called before any other function is called.
// It can be called as many times as needed, with different data sets
//

#include <string.h>

#include "Constants.cpp"
#include "Globals.cpp"
#include "Raw.cpp"
#include "Maw.cpp"
extern int runMaw(double diffMatrix[][NUM_GENE], int diffIndex);
extern int runRaw(double diffMatrix[][NUM_GENE], int diffIndex);
extern void getRanks(double diffMatrix[][NUM_GENE], int rank[][NUM_GENE]);

extern int Initialize(
    char **geneFullNames,
    char **geneShortNames,
    int nGenes,
    char *dataDir
    )
{
    if (!geneFullNames || !geneShortNames || !dataDir)
    {
        return -1;
    }

    if (nGenes <= 0)
    {
        return -1;
    }

    for (int i = 0; i < nGenes; i++)
    {
        if (!geneShortNames[i] || !geneFullNames[i])
            return -1;
    }

    int i;

    g_numGenes = nGenes;

    for (i = 0; i < nGenes; i++)
    {
        strcpy(g_strSpeciesShortName[i], geneShortNames[i]);
        strcpy(g_strSpeciesFullName[i], geneFullNames[i]);
    }

    i = strlen(dataDir);
    if (i == 0)
    {
        strcpy(g_strDataDir, ".");
    }
    else
    {
        strcpy(g_strDataDir, dataDir);
        if (g_strDataDir[i - 1] == '\\')
        {
            g_strDataDir[i - 1] = '\0';
        }
    }

    return 0;
}

extern  int getDiffMatrix(
    double diffMatrix[][NUM_GENE],
    int absWordType,
    int diffIndex
    )
{
    int ret = 0;

    if (absWordType != MAW && absWordType != RAW)
        return -1;

    if (absWordType == RAW)
    {
        ret = runRaw(diffMatrix, diffIndex);
    }
    else if (absWordType == MAW)
    {
        ret = runMaw(diffMatrix, diffIndex);
    }

    return 0;
}

extern  int getRanks(int rank[][NUM_GENE], int absWordType, int diffIndex)
{
    double diffMatrix[NUM_GENE][NUM_GENE] = {0};
    int ret = 0;

    ret = getDiffMatrix(diffMatrix, absWordType, diffIndex);


    int i, j;
    for (i = 0; i < g_numGenes; i++)
    {
        for (j = i; j < g_numGenes; j++)
        {
            diffMatrix[i][j] = diffMatrix[j][i];
        }


    }

    double mn, temp;
    int index, k;
    for ( k = 0 ; k < g_numGenes ; k++ )
    {
        for ( i = 0 ; i < g_numGenes ; i++ )
        {
            mn = diffMatrix[k][i];
            index = i;
            for ( j = 0 ; j < g_numGenes ; j++ )
            {
                if ( mn > diffMatrix[k][j] && diffMatrix[k][j] != 999999999 )
                {
                    mn = diffMatrix[k][j];
                    index = j;
                }
            }
            diffMatrix[k][index] = 999999999;
            rank[k][i] = index;
        }
    }

    return ret;
}
