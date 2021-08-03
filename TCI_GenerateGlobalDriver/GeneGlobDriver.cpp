/* 
 * File:   GENEGLOBDRIVER.cpp
 */

using namespace std;
#include "GeneGlobDriver.h"
#include <fstream>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <string>
#include <vector>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <algorithm>


/**
 * This function first calculate the global driver score, and then generate global driver  
 *  
 * @param gtMatrix     An TDIMatrix reference, in which SGA data of a collection of tumors are stored
 *                              in a tumor-by-gene matrix
 * @param geMatrix        An TDIMatrix reference, in which DEG data of a collection of tumors are stored
 *                              in a tumor-by-gene matrix
 * @param mapGlobDrivers        A reference to a dictionary of global drivers of DEGs in the "geMatrix",
                                Each entry is keyed by the DEG gene name, and the value is a reference of a 
                                vector<string> containing the two 2 global driver.  
 * @param outFileName       A string that contains the output file name include the path 
 * @param v0   A constant float
 */
//void TDIC(GTMatrix& gtMatrix, TDIMatrix& geMatrix, map<string, 
//        string>& mapGlobDrivers, const int tumorID, const string outPath, const float v0){
void GeneGlobDriver(GTMatrix& gtMatrix, TDIMatrix& geMatrix, const string outFileName, const float v0){
    //Calculate global driver score
    
    bool * gtDataMatrix = gtMatrix.getMatPtr();
    bool * geDataMatrix = geMatrix.getMatPtr();
    
     // Allocate memory for nGT x nGE matrix

    unsigned int nGT = gtMatrix.getNGenes();
    unsigned int nGE = geMatrix.getNGenes();
    int nTumors = gtMatrix.getNTumors();
    if (nTumors != geMatrix.getNTumors()) // number of tumors should agree between gtMatrix and geMatrix
    {
        cerr << "Bug: gtMatrix and geMatrix contain different number of tumors ";
    }

    /*Get names of mutated GTs and GEs for this tumor
     */
    vector<string> gtNames = gtMatrix.getGeneNames();
    vector<string> geNames = geMatrix.getGeneNames();
 
    float* tumorPosteriorMatrix = new float[nGT * nGE]();

    //ofstream file;
    ofstream outMatrixFile, outPosteriorMatrixFile;
    string marginalMatrixFileName = outFileName + ".marginal.matrix.csv";
    string posteriorMatrixFileName = outFileName + ".posterior.matrix.csv";
    try
    {
        outMatrixFile.open(marginalMatrixFileName.c_str());
        outPosteriorMatrixFile.open(posteriorMatrixFileName.c_str());
    }
    catch(ofstream::failure e)
    {
        cout << "Exception opening output file. Please ensure you have an existing directory for file.\n";
    }
    
    // Output the marginal likelihood table
    //output header
    for (int gt = 0; gt < nGT; gt ++ ){
        outMatrixFile << "," << gtNames[gt];
        outPosteriorMatrixFile << "," << gtNames[gt];
    }
    outMatrixFile << "\n";
    outPosteriorMatrixFile << "\n";

    #pragma omp parallel for
    for(unsigned int ge = 0; ge < nGE; ge++)
    {
        if (ge % 1000 == 0)
                printf("TCI_GD processing %dth phenotype.\n", ge);
        
        float normalizer = 0;
         unsigned int rowStartForGE = ge * nTumors; 
        
        // loop through each GT in the tumor
        for (unsigned int gt = 0; gt < nGT; gt++)
        {
            // variables to collect state of gt and ge|gt.  
            float T[2] = {0.0};
            float TE[4] = {0.0};

            int gtRowStart = gt * nTumors;
            
            for(int t = 0; t < nTumors; t++)
            {
                int tVal = gtDataMatrix[gtRowStart + t];
                int eVal = geDataMatrix[rowStartForGE + t];
                
                TE[tVal*2 + eVal] ++;  //conditional statistics is stored as a binary tree
            }

            //T[0xT] T is the gt value
            //e.g. T[1] save the count when gt value T = 1 
            T[0] = TE[0] + TE[1]; //T0 = T0E0 + T0E1
            T[1] = TE[2] + TE[3]; //T1 = T1E0 + T1E1  
                              
            float TFscore;
            if (gt == 0)
            {
                TFscore = calcA0Fscore(T[1],  TE[3], TE[2], T[0],  TE[1], TE[0]);
            }
            else 
            {
                TFscore = calcFscore( T[1],  TE[3], TE[2], T[0],  TE[1], TE[0] );
            }

            // 
            float gtPriorRaw = gtMatrix.getPriorPtr()[gt];
            float gtPrior = 0.0;
            if (gtPriorRaw == 0.0)
                gtPrior = -FLT_MAX;
            else
                gtPrior = log(gtPriorRaw);
            
            float lnData = TFscore + gtPrior;

            tumorPosteriorMatrix[gt * nGE + ge] = lnData;

            float pGT1GE1, pGT0GE1;
            if(gt == 0)
            {
                pGT1GE1 = (ALPHANULL + TE[3]) / (ALPHANULL + ALPHANULL + T[1]);
                pGT0GE1 = (ALPHANULL + TE[1]) / (ALPHANULL + ALPHANULL + T[0]);
                
            }
            else
            {
                pGT1GE1 = (ALPHAIJK11 + TE[3]) / (ALPHAIJK11 + ALPHAIJK10 + T[1]);
                pGT0GE1 = (ALPHAIJK01 + TE[1]) / (ALPHAIJK01 + ALPHAIJK00 + T[0]);                      
            }

            if(pGT1GE1 <= pGT0GE1)
            {
                tumorPosteriorMatrix[gt* nGE + ge] = -FLT_MAX;
            }
        }

        for(unsigned int gt = 0; gt < nGT; gt++)
        {
            if(gt == 0)
            {
                normalizer = tumorPosteriorMatrix[gt * nGE + ge];
            }
            else
            {
                normalizer = logSum(normalizer, tumorPosteriorMatrix[gt * nGE + ge]);
            }
        }
        
        // finished populating a marginal likelihood of GTs with respect to a given GE, normalize so that the column sum to 1
        // output results to files
        outMatrixFile << geNames[ge]  ;
        outPosteriorMatrixFile << geNames[ge];
        for (unsigned int gt = 0; gt < nGT; gt++)
        {
            outMatrixFile << "," << tumorPosteriorMatrix[gt * nGE + ge];
            tumorPosteriorMatrix[gt * nGE + ge] = exp(tumorPosteriorMatrix[gt * nGE + ge] - normalizer); 
            outPosteriorMatrixFile << "," << tumorPosteriorMatrix[gt * nGE + ge];
        }
        outMatrixFile << "\n";
        outPosteriorMatrixFile << "\n";
    }
    outPosteriorMatrixFile.close();
    outMatrixFile.close();

    // Output global driver list
    ofstream outFile;
    try
    {
        outFile.open(outFileName.c_str());
    }
    catch(ofstream::failure e)
    {
        cout << "Exception opening output file. Please ensure you have an existing directory for file.\n";
    }
 
       
    // Generate global driver from tumorPosteriorMatrix
    int* maxGtIndex = new int[nGE]() ;
    float* maxGtProb = new float [nGE]();
    for (int ge = 0; ge < nGE; ge++)
    {   
        float maxGtValue = 0;
        for (int gt = 0; gt < nGT; gt++ )
        {
            if(tumorPosteriorMatrix[gt * nGE + ge] > maxGtValue)
            {
                maxGtValue = tumorPosteriorMatrix[gt * nGE + ge];
                maxGtIndex[ge] = gt;   
                maxGtProb[ge] = tumorPosteriorMatrix[gt * nGE + ge];
            }                
        }
    }
    
    //write column headers
    for(int ge = 0; ge < nGE; ge++)
    {
        outFile << geNames[ge] << "," ;
        outFile << gtNames[maxGtIndex[ge]] << "," << maxGtProb[ge];
        outFile << "\n";
    }

    outFile.close();

    delete [] maxGtIndex;
    delete [] tumorPosteriorMatrix;
}

/********** logSum *********************************************************/ 
/**
 * Evaluate Ln(x + y)
 * @param lnx ln(x)
 * @param lny ln(y)
 * @return ln(x + y)
 */
float logSum(float lnx, float lny){
    float maxExp = -4950.0;

    if(lny > lnx){                
        float tmp = lnx;
        lnx = lny;
        lny = tmp;
    }

    float lnyMinusLnX = lny - lnx;
    float lnXplusLnY;

    if(lnyMinusLnX < maxExp)
        lnXplusLnY = lnx;
    else
        lnXplusLnY = log(1 + exp(lnyMinusLnX)) + lnx;

    return (lnXplusLnY); 
}


/***************** calcSingleGtFscore  **************************************/
/**
 * 
 * @param gt1 
 * @param gt1ge1
 * @param gt1ge0
 * @param gt0
 * @param gt0ge1
 * @param gt0ge0
 * @return 
 */
float calcFscore(float gt1,  float gt1ge1, float gt1ge0, 
    float gt0, float gt0ge1, float gt0ge0 )
{
    // Calculation of Equation 7    
    float glnNi0 = lgamma(ALPHAIJK00 + ALPHAIJK01) - lgamma(gt0 + ALPHAIJK00 + ALPHAIJK01);
    float glnNi1 = lgamma(ALPHAIJK10 + ALPHAIJK11) - lgamma(gt1 + ALPHAIJK10 + ALPHAIJK11);

    float fscore = glnNi0 + glnNi1;   
    fscore += lgamma(gt0ge0 + ALPHAIJK00) - lgamma(ALPHAIJK00);
    fscore += lgamma(gt0ge1 + ALPHAIJK01) - lgamma(ALPHAIJK01);
    fscore += lgamma(gt1ge0 + ALPHAIJK10) - lgamma(ALPHAIJK10);
    fscore += lgamma(gt1ge1 + ALPHAIJK11) - lgamma(ALPHAIJK11);

    return (fscore);
}


/***************** calcSingleGtFscore  **************************************/
/**
 * 
 * @param gt1
 * @param gt1ge1
 * @param gt1ge0
 * @param gt0
 * @param gt0ge1
 * @param gt0ge0
 * @return 
 */
float calcA0Fscore(float gt1,  float gt1ge1, float gt1ge0, 
    float gt0, float gt0ge1, float gt0ge0 )
{

    // Calculation of Equation 7    
    float glnNi0 = lgamma( ALPHANULL + ALPHANULL) - lgamma(gt0 + ALPHANULL + ALPHANULL);
    float glnNi1 = lgamma(ALPHAIJK10 + ALPHANULL) - lgamma(gt1 + ALPHANULL + ALPHANULL);

    float fscore = glnNi0 + glnNi1;   
    fscore += lgamma(gt0ge0 + ALPHANULL) - lgamma(ALPHANULL);
    fscore += lgamma(gt0ge1 + ALPHANULL) - lgamma(ALPHANULL);
    fscore += lgamma(gt1ge0 + ALPHANULL) - lgamma(ALPHANULL);
    fscore += lgamma(gt1ge1 + ALPHANULL) - lgamma(ALPHANULL);

    return (fscore);
}
