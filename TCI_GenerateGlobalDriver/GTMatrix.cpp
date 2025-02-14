/* 
 * File:   GTMatrix.cpp
 * Author: kevin
 * 
 * Created on May 21, 2015, 9:36 AM
 */

#include <math.h>
#include <algorithm>    // std::find
#include <vector>
#include <fstream> 
#include <iostream>
#include <string>
#include <string.h>
#include <sstream>
#include <float.h>

//#include "TDIMatrix.h"
#include "GTMatrix.h"

using namespace std;

GTMatrix::GTMatrix() : TDIMatrix() {
    //    allAlphaPriors = NULL;

}

GTMatrix::GTMatrix(string gtfileName, string priorfileName) {
    SGApriors = NULL;
    this->load(gtfileName, priorfileName);
    //    calcGlobalAlphaPriors();
}

GTMatrix::~GTMatrix() {
    //    if(allAlphaPriors) delete allAlphaPriors;
    if (SGApriors) delete SGApriors;
}

/**
 * Create array representing prior probability for all GTs.
 */
//void GTMatrix::calcGlobalAlphaPriors()

void GTMatrix::calcGlobalAlphaPriors(const float v0, vector<float>& allAlphaPriors) {

    /*
     * Iterate through all tumors of the GTMatrix.
     * For each tumor, find GTs that are mutated for that tumor.
     * Then calculate a normalized weight for the probability of one of these
     * GTs being the driving mutation. This weight is 1 / (number of mutated GTs).
     * Finally, update indeces in "allAlphaPriors" corresponding to these mutated
     * GTs by adding this normalized weight to the existing value in allAlphaPriors[GT].
     *
     */
    vector<int> geneIndices;
    for (int i = 0; i < nTumors; i++) {
        findGeneWithOnesInTumor(i, geneIndices);
        int numMutatedGTsForTumor = geneIndices.size();
        float gtWeight = 1.0 / (float) numMutatedGTsForTumor;
        for (int j = 0; j < numMutatedGTsForTumor; j++) {
            allAlphaPriors[geneIndices[j]] += gtWeight;
        }
        geneIndices.clear();
    }

    // by mxj start *********

    float sumAlphaPrior = 0.0;
    for (int i = 1; i < nGenes; i++) //we only use the alphas of real SGAs to calculate
        sumAlphaPrior += allAlphaPriors[i];

    for (int i = 0; i < nGenes; i++) {
        if (i == 0)
            allAlphaPriors[i] = log(v0);
        else
            if (allAlphaPriors[i] == 0)
            allAlphaPriors[i] = -FLT_MAX;
        else
            allAlphaPriors[i] = log((1 - v0) * allAlphaPriors[i] / sumAlphaPrior);
    }
    // by mxj end ***********
}
/**
 * 
 * @param gtIndx List of GT indices for which you want to calculate the log prior.
 * @param v0 
 * @param lnTumorPriorAlphas List of normalized, log prior probabilities for each of the GTs in gtIndx 
 *         (with original sequence intact).
 */


/**
 * Load GT or GE data from raw data (CSV file).
 * @param fileName Filepath for a CSV file representing a GT or GE matrix.
 */
void GTMatrix::load(string gtfileName, string priorfileName) {
    std::stringstream ss;
    std::string line;
    ifstream inFileStream;
    vector<int*> matrixAsVec;
    int nCol = 0, nRow = 0;

    try {
        inFileStream.open(gtfileName.c_str());
        if ((inFileStream.rdstate() & std::ifstream::failbit) != 0) {
            std::cerr << "Error opening file when loading GT matrix, quit.\n";
            inFileStream.close();
            exit(EXIT_FAILURE);
        }
    } catch (...) { //std::ifstream::failure e
        cerr << "Fail to open file " << gtfileName;

    }

    // read in first line to figure out the number of columns and column Names;
    getline(inFileStream, line);
    line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

    string tmp;
    stringstream firstSS(line);

    bool firstColFlag = true;
    geneNames.push_back("A0");
    while (getline(firstSS, tmp, ',')) {
        if (firstColFlag) {
            firstColFlag = false;
            continue;
        }
        geneNames.push_back(tmp);


        //geneIndxMap.insert(std::pair<string, int>(tmp, nCol));
        nCol++;
    }

    //User may accidently input GTMatrix as PanCanGTMatrix, so check the input file is actual GTMatrix
    //Check the last column to see if the column name is something like 'Cancel Type', if it is , then exit program
    // Cancel type column name can be 'Cancer Type' or 'Cancel Types' or 'Can Type' or 'Can Types' or 'Cancer Type Code' or 'CanType code'. 
    //There can be no or more space between each word. Each character in the word is not case sensitive, can be lower or upper case.

    string colName;
    colName = geneNames.back();
    //Convert to lower case
    transform(colName.begin(), colName.end(), colName.begin(), ::tolower);
    //remove all spaces
    colName.erase(remove(colName.begin(), colName.end(), ' '), colName.end());
    bool hasCanTypeCol = true;
    if (colName == "cancertype" || colName == "cancertypes" || colName == "cantype" || colName == "cantypes"
            || colName == "cantypecode" || colName == "cancertypecode" || colName == "cantypescode" || colName == "cancertypescode") {
        cout << "Input PanCanGTMatrix has cancer type \n ";
        geneNames.pop_back();
        nCol--;
    } else
        hasCanTypeCol = false;

    while (getline(inFileStream, line)) {
        firstColFlag = true;

        stringstream anotherss(line);
        string tmp;
        int curCol = 0;
        while (getline(anotherss, tmp, ',')) {
            if (firstColFlag) {
                firstColFlag = false;
                tumorNames.push_back(tmp);
                matrixAsVec.push_back(new int[nCol]);
                continue;
            }
            if (hasCanTypeCol) {
                if (curCol != nCol) {
                    matrixAsVec[nRow][curCol] = atoi(tmp.c_str());
                    curCol++;
                }
            } else {
                matrixAsVec[nRow][curCol] = atoi(tmp.c_str());
                curCol++;
            }
        }
        nRow++;
    }
    inFileStream.close();

    // transform the vector of inter arrays to a consecutive array so that 
    // computation can be done efficiently
    int indx = 0;
    mat = new bool [(nCol + 1) * nRow ]();
    for (int i = 0; i < nRow; i++) {
        mat[indx++] = 1;
    }

    for (int g = 0; g < nCol; g++) {
        for (int t = 0; t < nRow; t++) {
            mat[indx++] = matrixAsVec[t][g];
        }
    }

    nTumors = tumorNames.size();
    nGenes = geneNames.size();

    // free temporary matrix
    for (int i = 0; i < matrixAsVec.size(); i++)
        delete[] matrixAsVec[i];

    //////////////////////////////////////////////////////////////////////////////////
    // This section read in prior probability matrix.
    ifstream inFileStream_p;
    try {
        inFileStream_p.open(priorfileName.c_str());
        if ((inFileStream_p.rdstate() & std::ifstream::failbit) != 0) {
            std::cerr << "Error opening file when loading global prior file, quit.\n";
            inFileStream_p.close();
            exit(EXIT_FAILURE);
        }
    } catch (...) { //std::ifstream::failure e
        cerr << "Fail to open file " << priorfileName;

    }

    // read in first line which is the SGA names, and do nothing;
    // There are two lines in prior file. we read in the second line to get the prior data
    SGApriors = new float[nCol+1]();//includes A0
    std::string line_p;
        
    string tmp_p;
    int curCol_p = 0;
    double f_prior;

    getline(inFileStream_p, line_p);
    //read in the second line
    getline(inFileStream_p, line_p);
    std::stringstream ss_p(line_p);
 
    while (getline(ss_p, tmp_p, ',')) {
        stringstream s_tmp_p;
        s_tmp_p << tmp_p;
        s_tmp_p>>f_prior;
        SGApriors[curCol_p++] = f_prior;
    }
    if (curCol_p != nCol+1){
        cerr << "Global prior has different number of genes than GtMatrix. #prior:" << curCol_p <<", #GtMatrix:"<< nCol + 1;
        exit(1);
    }
        
    inFileStream_p.close();
}

