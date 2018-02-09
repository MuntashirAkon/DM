//
// This file contains necessary code for Minimal Absent Word (MAW)
// based difference matrix calculation
//
// The lower triangle of the diff matrix is filled. This is because
// the code used to generate phylogenetic tree needs the lower triangle
//

#include <string>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include "Set.h"
#include "Constants.cpp"
#include "Globals.cpp"

double calculateTVD(
    vector<pair<int, double> > rgA,
    vector<pair<int, double> > rgB
    );

int runMaw(double diffMatrix[][NUM_GENE], int diffIndex){
    int i, j;
    Set maw[NUM_GENE];
    Set diff, a, b;
    string filename;

    //
    // 1. Read Input
    // 2. Compute Symmetric Difference
    // 3. Calculate index
    // 4. Produce the difference matrix for the 11 genes
    //
    for (i = 0; i < g_numGenes; i++){
        filename = string(g_strDataDir) + '/' + string(g_strSpeciesFullName[i]) + ".maw.txt";
        maw[i] = Set::CreateFromFile(filename);
    }

    for (i = 0; i < g_numGenes; i++){
        for (j = 0; j < i; j++){
            switch(diffIndex){
                case MAW_LWI_SDIFF:
                    diff = maw[i].SymmetricDifference(maw[j]);
                    diffMatrix[i][j] = diffMatrix[j][i] = diff.LengthWeightedIndex();
                    break;

                case MAW_LWI_INTERSECT:
                    diff = maw[i].Intersection(maw[j]);
                    diffMatrix[i][j] = diffMatrix[j][i] = -diff.LengthWeightedIndex();
                    break;

                case MAW_GCC_SDIFF:
                    diff = maw[i].SymmetricDifference(maw[j]);
                    diffMatrix[i][j] = diffMatrix[j][i] = diff.GCContent();
                    break;

                case MAW_GCC_INTERSECT:
                    diff = maw[i].Intersection(maw[j]);
                    diffMatrix[i][j] = diffMatrix[j][i] = 1.0 - diff.GCContent();
                    break;

                case MAW_JD:

                    a = maw[i].Union(maw[j]);

                    b = maw[i].Intersection(maw[j]);

                    diffMatrix[i][j] = diffMatrix[j][i] = 1.0 - 1.0 * b.Cardinality() / a.Cardinality();

                    break;

                case MAW_TVD:
                    diffMatrix[i][j] = diffMatrix[j][i]
                                    = calculateTVD(
                                            maw[i].LengthDistribution(),
                                            maw[j].LengthDistribution()
                                            );
                    break;
            }

        }
    }

    return 0;
}

//
// Calculate the Total Variation Distance (TVD) between
// 2 sets of (minimal absent) words.
//
// TVD(P, Q) = 0.5 * sum(abs(P(i) - Q(i))), over all lengths i
//
// Assumption: the vectors are sorted in length order
//
double calculateTVD(
    vector<pair<int, double> > rgA,
    vector<pair<int, double> > rgB
    )
{
    double res = 0.0;

    int indexA, indexB;
    double distA, distB;
    int lenA, lenB;

    for (indexA = 0, indexB = 0; indexA < (int)rgA.size() || indexB < (int)rgB.size(); )
    {
        distA = distB = 0.0;
        lenA = lenB = 999999;

        if (indexA < (int)rgA.size() && indexB < (int)rgB.size())
        {
            lenA = rgA[indexA].first;
            lenB = rgB[indexB].first;

            if (lenA == lenB)
            {
                distA = rgA[indexA++].second;
                distB = rgB[indexB++].second;
            }
            else if (lenA < lenB)
            {
                distA = rgA[indexA++].second;
            }
            else
            {
                distB = rgB[indexB++].second;
            }
        }
        else if (indexA < (int)rgA.size())
        {
            distA = rgA[indexA++].second;
        }
        else
        {
            distB = rgB[indexB++].second;
        }
        if (distA > distB) { res += (distA - distB); }
        else { res += (distB - distA); }
    }

    res /= 2.0;

    return res;
}

