//
// This file contains declaration of the Set class. This
// class is used extensively in different set operations
// used in difference matrix calculations.
//

#pragma once

#include <vector>
#include <string>
using namespace std;

class Set
{
    vector<string> rgElements;
    bool m_fBulkAddElements;

public:
    Set();
    ~Set();
    Set operator=(const Set& B);
    Set(const Set& B);

    static Set CreateFromFile(char *pstrFileName);

    void StartBulkAddElements();
    void EndBulkAddElements();

    void AddElement(string strElement);
    Set Union(Set &B);
    Set Intersection(Set &B);
    Set Difference(Set B);
    Set Difference(Set &B);
    Set SymmetricDifference(Set &B);

    void Print();

    double LengthWeightedIndex();
    double GCContent();
    double GCContent_2();
    int Cardinality();

    vector<pair<int, double> > LengthDistribution();

    //
    // Compare 2 elements of the set for ordering
    // Used for sorting purposes.
    //
    // Returns true if a is smaller than b
    //
    static bool cmp (string a, string b)
    {
        return a.compare(b) < 0;
    }

private:
    void CopyInternal(const Set& B);
};
