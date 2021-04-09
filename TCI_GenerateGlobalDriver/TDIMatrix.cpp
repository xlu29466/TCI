/* 
 * File:   TDIMatrix.cpp
 * Author: Kevin Lu
 * 
 * Created on April 25, 2015, 6:13 PM
 */
using namespace std;

#include "TDIMatrix.h"
#include <algorithm>    // std::find
#include <vector>
#include <fstream> 
#include <iostream>
#include <string>
#include <string.h>
#include <sstream>

//using namespace std;

TDIMatrix::TDIMatrix(){
    mat = NULL;
    // other initialization stuff
}

TDIMatrix::TDIMatrix(string fileName){
    this->load(fileName);    
    
}

TDIMatrix::TDIMatrix(const TDIMatrix& orig) {
    // copy constructor not implemented.
}


TDIMatrix::~TDIMatrix() {
    if (mat) delete [] mat; 
}


/**
 * 
 * @param geneIndx index of gene of interest
 * @param tumorIndx index of tumor of interest
 * @return mutation or expresssion value of the given gene for the given tumor
 */
int TDIMatrix::valueAt(int geneIndx, int tumorIndx){
    return mat[geneIndx * nTumors + tumorIndx];
}



/**
 * 
 * @param inGeneNames list of tumors
 * @param outGeneIndices corresponding indeces for the tumors in "inGeneNames"
 */


/**
 * 
 * @param fileName filepath of CSV file to be converted into a TDIMatrix object
 */
void TDIMatrix::load(string fileName)
{
    std::stringstream ss;
    std::string line;
    ifstream inFileStream;   
    vector<bool*> matrixAsVec;
    int nCol = 0, nRow = 0;

    try
    {
        inFileStream.open(fileName.c_str());
        if ( (inFileStream.rdstate() & std::ifstream::failbit ) != 0 )
        {
            std::cerr << "Error opening file when loading TDI matrix, quit.\n";
            inFileStream.close();
            exit(EXIT_FAILURE);
        }
    }
    catch (...) 
    { //std::ifstream::failure e
        cerr << "Fail to open file " << fileName;
        
    }

    //cout << "Opened GE file ok.\n";
            
    // read in first line to figure out the number of columns and column Names;
    getline(inFileStream, line);
    line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

    string tmp;
    stringstream firstSS(line);

    bool firstColFlag = true;
    while (getline(firstSS, tmp, ',' ))
    {
        if(firstColFlag)
        {
            firstColFlag = false;
            continue;
        }
        geneNames.push_back(tmp);

        //geneIndxMap.insert(std::pair<string, int>(tmp, nCol));
        nCol++;
    }

    while (getline(inFileStream, line))
    {
        firstColFlag = true; 

        stringstream anotherss(line);
        string tmp;
        int curCol = 0;
        while (getline(anotherss, tmp, ','))
        {
            if(firstColFlag)
            {
                firstColFlag = false;
                tumorNames.push_back(tmp);
                
 
                matrixAsVec.push_back(new bool[nCol]());
                continue;
            }

            matrixAsVec[nRow][curCol] = atoi(tmp.c_str());
            curCol++;
        }        
        nRow++;		
    }

    inFileStream.close();          
    
    // transform the vector of inter arrays to a consecutive array so that 
    // computation can be done efficiently
    mat = new bool [nCol * nRow]();
    int indx = 0;    
    for (int g = 0; g < nCol; g++){
        for (int t = 0; t < nRow; t++){
            mat[indx++] = matrixAsVec[t][g];
        }
    }
    
    nTumors = tumorNames.size();
    nGenes = geneNames.size();
    
    // free temporary matrix
    for (int i = 0; i < matrixAsVec.size(); i++) 
        delete [] matrixAsVec[i];
}


/**
 * 
 * @param geneID index of gene
 * @param tumorIndices vector of indices representing tumors that have a mutation in the given gene
 */
//void TDIMatrix::findTumorsWithOnesPerGene(int geneID, vector<int>& tumorIndices)
//{
//    for (int t = 0; t < nTumors; t++){
//        if (mat[nTumors * geneID + t] == 1) tumorIndices.push_back(t);
//    }    
//}
//
///**
// * 
// * @param tumorID index of tumor
// * @param geneIndices vector of indices representing genes that are mutated in the given tumor
// */
// //should this return A0 if it is a 1?
void TDIMatrix::findGeneWithOnesInTumor(int tumorID, vector<int>& geneIndices){
    for (int g = 0; g < nGenes; g ++ ){
        if (mat[g * nTumors + tumorID] == 1) geneIndices.push_back(g);
    }
}

bool TDIMatrix::writeToCSV(string outFilePath)
{
    ofstream outFile;
    try
    {
        outFile.open(outFilePath.c_str());
    }
    catch(ofstream::failure e)
    {
        cout << "Exception opening output file. Please ensure you have an existing directory for file.\n";
        return false;
    }
    
    //start writing CSV representation of TDIMatrix
    
    //write column headers
    for(int i = 0; i < geneNames.size(); i++)
    {
        outFile << "," << geneNames[i];
    }
    outFile << "\n";
    
    for(int i = 0; i < tumorNames.size(); i++)
    {
        outFile << tumorNames[i];
        for(int j = 0; j < geneNames.size(); j++)
        {
            outFile << "," << mat[j * nTumors + i];
        }
        outFile << "\n";
    }
    
    outFile.close();
    return true;
}
