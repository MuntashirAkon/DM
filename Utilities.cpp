//
// This file contains some utility code.
// 
// The getRanks method ranks the different species relative to one species
// based on the difference matrix. It does this analysis for each species and
// generates a 2D rank table.
//
// In addition, this file contains a method PrintDiffMatrix_LatexFormat
// that was used to suitably format the difference matrix into latex style
// table (with higher triangle)
//

#include <cstdio>
#include <algorithm>
#include <numeric>
#include "Constants.h"
#include "Globals.h"

using namespace std;

class IndexComparer
{
protected:
    double *m_pVal;
    
public:
    IndexComparer(double *pVal)
        : m_pVal(pVal)
    {
    }

    bool operator() (size_t lindex, size_t rindex)
    {
        return m_pVal[lindex] < m_pVal[rindex];
    }
};

void getRanks(double diffMatrix[][NUM_GENE], int rank[][NUM_GENE])
{
    int i;

    for (i = 0; i < g_numGenes; i++)
    {
        iota(rank[i], rank[i] + g_numGenes, 0);
        sort(rank[i], rank[i] + g_numGenes, IndexComparer(diffMatrix[i]));
    }
}

//void PrintDiffMatrix_LatexFormat(double diffMatrix[NUM_GENE][NUM_GENE])
//{
//    int i, j;
//
//    printf("\\hline\n");
//    printf("Species");
//    for (i = 0; i < NUM_GENE; i++)
//    {
//        printf(" & %7s", g_strSpeciesFullName[i]);
//    }
//    printf("\\\\ \\hline\n");
//    
//    for (i = 0; i < NUM_GENE; i++)
//    {
//        printf("%7s", g_strSpeciesFullName[i]);
//        for (j = 0; j <= i; j++)
//        {
//            printf(" &        ");
//        }
//
//        for (; j < NUM_GENE; j++)
//        {
//            printf(" & %7.2lf", diffMatrix[i][j]);
//        }
//        printf("\\\\\n");
//    }
//
//    printf("\\hline\n");
//}

