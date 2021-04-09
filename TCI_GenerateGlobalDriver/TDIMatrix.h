/* 
 * File:   TDIMatrix.h
 * 
 * Author: Kevin Lu
 *
 * Created on April 25, 2015, 6:13 PM
 * 
 * This matrix class is developed to contain the SGA and DEG data for TDI analysis.
 * 
 * This class consists of a matrix hold cancer genomic  data and a set of functions 
 * to extract information from the matrix.   The matrix is organized in tumor-by-gene
 * format, which simplify the pre-processing of the data.   
 */



#include <string>
#include <vector>
#include <map>


using namespace std;

#ifndef TDIMATRIX_H
#define	TDIMATRIX_H

class TDIMatrix {
public:
    // constructors
    TDIMatrix();
    TDIMatrix(const TDIMatrix& orig);    
    TDIMatrix(string fileName);
    virtual ~TDIMatrix();
    
    // Basic functions
    void load(string );
    int valueAt(int, int);
    
    bool* getMatPtr(void) const {return mat;};
    
    void findGeneWithOnesInTumor(int, vector<int>& );
    
    // other look ups
    vector<string>& getGeneNames(void) {return geneNames;}; 
    int getNTumors(void){return tumorNames.size();};
    int getNGenes(void){return geneNames.size();};


    //write functions
    bool writeToCSV(string outFilePath);


protected:
    unsigned int nTumors;
    unsigned int nGenes;
    vector<string> tumorNames;
    vector<string> geneNames;
    bool* mat;   // we represent the matrix as a consecutive liner array 
                // with  gene-by-tumor.  This way, it is more efficient to 
                // collect statistics from tumors because we loop through tumors
    //map<string, int> geneIndxMap;
    
};

#endif	/* MATRIX_H */

