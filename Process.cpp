//
// Created by Muntashir Al-Islam on 3/13/18.
//
#import "Process.h"

Process::Process(int absent_word_type, int dissimilarity_index_type, const string& input_dir, const string& output_dir) {
    this->aw_type = absent_word_type;
    this->di_type = dissimilarity_index_type;
    this->input_dir = input_dir;
    this->output_dir = output_dir;
    this->total_genes = 0;
    this->largest_species_len = 0;
    this->getSpecies();
    for (int i = 0; i < this->total_genes; i++) {
        this->diffMatrix.push_back(vector<double>(static_cast<unsigned int>(this->total_genes)));
        this->rank.push_back(vector<int>(static_cast<unsigned int>(this->total_genes)));
    }
    this->getDiffMatrix();
}

void Process::getSpecies() {
    ifstream speciesFull(this->input_dir + '/' + SPECIES_FULL);
    if(!speciesFull.is_open()){
        cerr << "The source directory doesn't contain \"SpeciesFull.txt\"!" << endl;
        throw exception();
    }
    // Get line has problems with carriage return
    for (string line; getline(speciesFull, line); ++(this->total_genes)){
        if(largest_species_len < line.length()) largest_species_len = static_cast<int>(line.length());
        this->species.push_back(line);
        line.clear();
    }
}

int Process::getDiffMatrix(){
    return (this->aw_type == RAW ? this->runRaw() : this->runMaw());
}

void Process::getRanks(){
    int i, j;
    for (i = 0; i < this->total_genes; i++) {
        for (j = i; j < this->total_genes; j++) {
            diffMatrix[i][j] = diffMatrix[j][i];
        }
    }

    double mn;
    int index, k;
    for ( k = 0 ; k < this->total_genes ; k++ ) {
        for ( i = 0 ; i < this->total_genes ; i++ ) {
            mn = diffMatrix[k][i];
            index = i;
            for ( j = 0 ; j < this->total_genes ; j++ ) {
                if ( mn > diffMatrix[k][j] && diffMatrix[k][j] != INFINITY ) {
                    mn = diffMatrix[k][j];
                    index = j;
                }
            }
            diffMatrix[k][index] = INFINITY;
            rank[k][i] = index;
        }
    }
}

void Process::printDiffMatrix(){
    int i, j;
    FILE *f1 = fopen((this->output_dir + '/' + DISTANCE_MATRIX).c_str(), "w+");
    // Formatted for Match7
    FILE *f2 = fopen((this->output_dir + '/' + "Output.txt").c_str(), "w+");

    for (i = 0; i < this->total_genes; i++) {
        cout << "{ ";
        for (j = 0; j < i; j++) {
            printf("%5.2lf", diffMatrix[i][j]);
            fprintf(f2,"%5.2lf", diffMatrix[i][j]);
            fprintf(f1, "%5.4lf\n", diffMatrix[i][j]);
            printf("%s", (j == i - 1) ? "" : ", ");
            fprintf(f2,"%s", (j == i - 1) ? "" : ", ");
        }
        cout << " }" << ((i == this->total_genes - 1) ? "" : ",") << '\n';
        fprintf(f2," %s\n", (i == this->total_genes - 1) ? "" : ",");
    }
}

void Process::printSortedSpeciesRelations(){
    this->getRanks();
    int i, j;
    ofstream text(this->output_dir + '/' + SPECIES_RELATION_TXT);
    ofstream json(this->output_dir + '/' + SPECIES_RELATION_JSON);

    if(!text.is_open() && !json.is_open()) return;

    json << "{" << '\n';
    for (i = 0; i < this->total_genes; i++){
        text << setw(this->largest_species_len) << this->species[i] << ":";
        json << "    \"" << this->species[i] << "\": [" << '\n';
        for (j = 0; j < this->total_genes; j++){
            text << " -> " << setw(this->largest_species_len) << this->species[rank[i][j]] << ((j != this->total_genes - 1) ? " ": "\n");
            json << "        \"" << this->species[rank[i][j]] << "\"" << ((j != this->total_genes - 1) ?  ",": "") << '\n';
        }
        json << "    ]" << ((i != this->total_genes - 1) ? ", " : "") << '\n';
    }
    json << "}" << '\n';
    text.close();
    json.close();
}

int Process::runRaw(){
    int i, j;

    if (this->di_type < RAW_LWI || this->di_type > RAW_GCC)
        return -1;

    for (i = 0; i < this->total_genes; i++) {
        diffMatrix[i][i] = 0;
        for (j = 0; j < i; j++) {
            diffMatrix[i][j] = diffMatrix[j][i] = GetRawBasedDiff(i, j, this->di_type);
        }
    }

    return 0; // Is there any way to filter it? I think it can be by checking if all the DM are 0.0
}

int Process::runMaw(){
    int i, j;
    vector<Set> maw(static_cast<unsigned int>(this->total_genes));
    Set diff, a, b;
    string filename;

    //
    // 1. Read Input
    // 2. Compute Symmetric Difference
    // 3. Calculate index
    // 4. Produce the difference matrix for the 11 genes
    //
    for (i = 0; i < this->total_genes; i++){
        filename = this->input_dir + '/' + this->species[i] + ".maw.txt";
        maw[i] = Set::CreateFromFile(filename);
    }

    for (i = 0; i < this->total_genes; i++){
        for (j = 0; j < i; j++){
            switch(this->di_type){
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
                            = calculateTVD(maw[i].LengthDistribution(), maw[j].LengthDistribution());
                    break;
                default:
                    cerr << "No matching dissimilarity index!" << endl;
                    throw exception();
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
double Process::calculateTVD(vector<pair<int, double> > rgA, vector<pair<int, double> > rgB) {
    double res = 0.0;

    int indexA, indexB;
    double distA, distB;
    int lenA, lenB;

    for (indexA = 0, indexB = 0; indexA < (int)rgA.size() || indexB < (int)rgB.size(); ) {
        distA = distB = 0.0;
//        lenA = lenB = INFINITY;

        if (indexA < (int)rgA.size() && indexB < (int)rgB.size()) {
            lenA = rgA[indexA].first;
            lenB = rgB[indexB].first;

            if (lenA == lenB) {
                distA = rgA[indexA++].second;
                distB = rgB[indexB++].second;
            } else if (lenA < lenB) {
                distA = rgA[indexA++].second;
            } else {
                distB = rgB[indexB++].second;
            }
        } else if (indexA < (int)rgA.size()) {
            distA = rgA[indexA++].second;
        } else {
            distB = rgB[indexB++].second;
        }

        if (distA > distB) {
            res += (distA - distB);
        } else {
            res += (distB - distA);
        }
    }

    res /= 2.0;

    return res;
}

double Process::GetRawBasedDiff(int i, int j, int diffIndex){
    int ref, main;
    double res = 0.0;
    Set raw1, raw2;

    //
    // Make the file name from the gene indices
    //
    ref = i, main = j;
    string filename1 = this->input_dir + '/' + this->species[ref] + '/'
                       + this->species[ref] + '_' + this->species[main] + ".raw.txt";
    ref = j, main = i;
    string filename2 = this->input_dir + '/' + this->species[ref] + '/'
                       + this->species[ref] + '_' + this->species[main] + ".raw.txt";

    raw1 = Set::CreateFromFile(filename1);
    raw2 = Set::CreateFromFile(filename2);

    if (diffIndex == RAW_LWI){
        res = 0.5 * (raw1.LengthWeightedIndex() + raw2.LengthWeightedIndex());
    }else if (diffIndex == RAW_GCC){
        res = 0.5 * (raw1.GCContent() + raw2.GCContent());
    }

    return res;
}
