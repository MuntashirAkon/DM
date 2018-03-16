//
// This file contains definition/implementation of the Set class.
// This class is used extensively in different set operations used
// in difference matrix calculations.
//

#include "Set.h"

using namespace std;

Set Set::CreateFromFile(string filename){
    // Create new Set
    Set set;
    cout << filename << endl;
    ifstream species(filename);
    if (!species.is_open()){
        cerr << "The file " << filename << " is not found!" << endl;
        throw exception();
    }
    // Add elements
    set.StartBulkAddElements();
    for(string line; getline(species, line); ) set.AddElement(line);
    set.EndBulkAddElements();
    species.close();
    return set;
}

Set::Set(){ m_fBulkAddElements = false; }

Set::Set(const Set& B) { CopyInternal(B); } // NOLINT

Set& Set::operator=(const Set& B){
    CopyInternal(B);
    return *this;
}

Set::~Set(){
    rgElements.clear();

    //
    // Release the memory
    //
    //rgElements.shrink_to_fit();
}

void Set::CopyInternal(const Set& B){
    m_fBulkAddElements = false;
    rgElements.assign(
        B.rgElements.begin(),
        B.rgElements.end()
        );
}


void Set::StartBulkAddElements(){
    m_fBulkAddElements = true;
}

void Set::EndBulkAddElements(){
    m_fBulkAddElements = false;

    //
    // bulk addition of elements in done. Let's sort the list.
    // By sorting the list, later set operations can be performed
    // efficiently
    //
    sort(rgElements.begin(), rgElements.end(), cmp);
}

void Set::AddElement(string strElement){
    unsigned long i, j;

    if (m_fBulkAddElements){
        rgElements.push_back(strElement);
    }else{
        //
        // We ensure that the vector remains sorted after adding an element
        // O(n) operations
        //
        for (i = 0; i < rgElements.size() && cmp(rgElements[i], strElement); i++);

        if (i == rgElements.size()){
            // We simply add the new element at the end
            rgElements.push_back(strElement);
        }else{
            // We make room for the new element at the right spot
            // by displacing the larger items 1 position to the right

            //
            // Save the last element, we will add it at the end of the
            // vector later on.
            //
            string strTemp = rgElements.back();

            //
            // Shift the elements to the right, one spot
            //
            for (j = rgElements.size() - 1 ; j > i ; j--){
                rgElements[j] = rgElements[j - 1];
            }

            rgElements[i] = strElement;

            rgElements.push_back(strTemp);
        }
    }
}

Set Set::Union(Set &B){
    assert(!m_fBulkAddElements);

    Set C;
    unsigned long i, j;

    C.rgElements.assign(rgElements.begin(), rgElements.end());

    for (i = 0; i < (int)B.rgElements.size(); i++){
        C.rgElements.push_back(B.rgElements[i]);
    }

    sort(C.rgElements.begin(), C.rgElements.end(), cmp);

    //
    // Duplicate Removal
    //
    for (i = 0, j = 0; i < (int)C.rgElements.size(); j++){
        string str = C.rgElements[i++];
        while (i < (int)C.rgElements.size() && 0 == str.compare(C.rgElements[i])){
            i++;
        }

        C.rgElements[j] = str;
    }

    //
    // Now remove the redundant elements from the back
    // (which have already been shifted to earlier positions)
    // so that the vector size remains correct
    //
    i = C.rgElements.size();
    while (j < i){
        C.rgElements.pop_back();
        j++;
    }

    return C;
}

Set Set::Intersection(Set &B){
    assert(!m_fBulkAddElements);

    Set C;
    int i;

    for (i = 0; i < (int)rgElements.size(); i++)
    {
        string str = rgElements[i];
        if (binary_search(B.rgElements.begin(), B.rgElements.end(), str, cmp))
        {
            C.rgElements.push_back(str);
        }
    }

    return C;
}

Set Set::Difference(Set &B){
    assert(!m_fBulkAddElements);

    Set C;
    int i;

    for (i = 0; i < (int)rgElements.size(); i++){
        string str = rgElements[i];
        if (!binary_search(B.rgElements.begin(), B.rgElements.end(), str, cmp)){
            C.rgElements.push_back(str);
        }
    }

    return C;
}

Set Set::Difference(Set B){
    assert(!m_fBulkAddElements);

    Set C;
    int i;

    for (i = 0; i < (int)rgElements.size(); i++){
        string str = rgElements[i];
        if (!binary_search(B.rgElements.begin(), B.rgElements.end(), str, cmp)){
            C.rgElements.push_back(str);
        }
    }

    return C;
}

Set Set::SymmetricDifference(Set &B){
    assert(!m_fBulkAddElements);

    return Union(B).Difference(Intersection(B));
}

void Set::Print(){
    assert(!m_fBulkAddElements);

    int i;

    for (i = 0; i < (int)rgElements.size(); i++){
        printf("%s\n", rgElements[i].c_str());
    }
    printf("\n");
}

double Set::LengthWeightedIndex(){
    double res = 0.0;
    int i;

    for (i = 0; i < (int)rgElements.size(); i++){
        unsigned long len = rgElements[i].length();

        res += 1.0 / (len * len);
    }

    return res;
}

//
// Computes the ratio of 'G' and 'C' (combined), compared to
// the total number of base pairs in the entire set
//
double Set::GCContent(){
    int i, j;
    double res = 0.0;
    double numerator = 0.0;
    double denominator = 0.0;

    for (i = 0; i < (int)rgElements.size(); i++){
        denominator += rgElements[i].length();

        for (j = 0; j < (int)rgElements[i].length(); j++){
            char ch = rgElements[i][j];
            if ('G' == ch || 'C' == ch){
                numerator++;
            }
        }
    }

    //
    // 0 divide should not occur here, since we expect
    // the set to contain at least 1 non-empty strand
    //
    assert(denominator > 0);

    res = numerator / denominator;

    return res;
}

int Set::Cardinality(){
    return (int)rgElements.size();
}

//
// Generates the distribution of words for different lengths.
// Each element in the vector specifies a length, followed by
// how many words were there of that length.
//
vector<pair<int, double> > Set::LengthDistribution(){
    map<int, int> mapDist;
    vector<pair<int, double> > rgDist;
    int i, len;
    i = 0;
    while (i < rgElements.size()){
        len = static_cast<int>(rgElements[i++].length());
        if (mapDist.find(len) == mapDist.end()){
            mapDist[len] = 1;
        }else{
            mapDist[len]++;
        }
    }

    map<int, int>::iterator it;
    for (it = mapDist.begin(); it != mapDist.end(); it++){
        rgDist.push_back(pair<int, double>(it->first, 1.0 * it->second / rgElements.size()));
    }
    sort(rgDist.begin(), rgDist.end());

    //
    // Double check that the probability distribution
    // sums to 1.0
    //
#if 1
    double p = 0.0;
    for (i = 0; i < (int)rgDist.size(); i++){
        p += rgDist[i].second;
    }
#endif

    return rgDist;
}
