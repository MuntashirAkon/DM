//
// This file contains some utility code.
//
// In addition, this file contains a method PrintDiffMatrix_LatexFormat
// that was used to suitably format the difference matrix into latex style
// table (with higher triangle)
//

#include <cstdio>
#include <algorithm>
#include <numeric>

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

