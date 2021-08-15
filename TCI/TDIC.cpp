/* 
 * File:   TDIC.cpp
 * Author: Kevin Lu
 * 
 * Created on April 25, 2015, 6:13 PM
 */

using namespace std;
#include "TDIC.h"
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
 * This function performs tumor-specific driver identification.  It calculate the causal score for all 
 * GT-vs-GE pairs observed in a given tumor, populate a GT-by-GE score matrix and output to file.  
 * @param gtMatrix     An TDIMatrix reference, in which SGA data of a collection of tumors are stored
 *                              in a tumor-by-gene matrix
 * @param geMatrix        An TDIMatrix reference, in which DEG data of a collection of tumors are stored
 *                              in a tumor-by-gene matrix
 * @param mapGlobDrivers        A reference to a dictionary of global drivers of DEGs in the "geMatrix",
                                Each entry is keyed by the DEG gene name, and the value is a reference of a 
                                vector<string> containing the two 2 global driver.  
 * @param tumorID               The tumor to be process
 */
void TDIC(GTMatrix& gtMatrix, TDIMatrix& geMatrix, map<string, 
        vector<string> > & mapGlobDrivers, const int tumorID, const string outPath, const float v0){
    // initializations 
 
    string curTumorName = gtMatrix.getTumorNameById(tumorID);
   
    // Find the GTs and GEs of the current tumor
    vector<int> tumorGtIndices, tumorGeIndices;
    gtMatrix.findGeneWithOnesInTumor(tumorID, tumorGtIndices);
    geMatrix.findGeneWithOnesInTumor(tumorID, tumorGeIndices);
    
    bool* gtDataMatrix = gtMatrix.getMatPtr();
    bool* geDataMatrix = geMatrix.getMatPtr();

    float* priorMatrix = gtMatrix.getPriorPtr();
    // Find the global drivers corresponding to the
    vector<int> tumorGlobDriverIndices1, tumorGlobDriverIndices2;
    //map <string, string> globalDriversMap;
    if (!getDEGGlobDriverIndices(gtMatrix, geMatrix, mapGlobDrivers, tumorGeIndices, tumorGlobDriverIndices1, tumorGlobDriverIndices2))
    {
        cout << "Error occurred when retrieving global drivers";
    }

    // Allocate memory for nGT x nGE matrix
    unsigned int nGT = tumorGtIndices.size();
    unsigned int nGE = tumorGeIndices.size();
    int nTumors = gtMatrix.getNTumors();
    if (nTumors != geMatrix.getNTumors()) // number of tumors should agree between gtMatrix and geMatrix
    {
        cerr << "Bug: gtMatrix and geMatrix contain different number of tumors ";
    }

    cout << "Processing tumor " << curTumorName << " with " << nGT << " GAs, and " << nGE << " GEs" << "\n";

    /*Get names of mutated GTs and GEs for this tumor
     */
    vector<string> gtNames;
    vector<string> geNames;
    
    gtMatrix.getGeneNamesByIndices(tumorGtIndices, gtNames);
    geMatrix.getGeneNamesByIndices(tumorGeIndices, geNames);
  
    float* tumorPosteriorMatrix = new float[nGT * nGE]();    

    // loop through each GE
    #pragma omp parallel for
    for(unsigned int ge = 0; ge < nGE; ge++)
    {
        float normalizer = 0;
        unsigned int curGeIndx = tumorGeIndices[ge];
        unsigned int rowStartForGE = curGeIndx * nTumors; 
        
        // find the globDriver for this give ge   
        unsigned int curGDriverIndx = tumorGlobDriverIndices1[ge]; //curGDriverIndx is found by ge indx    
        unsigned int rowStartForGlobDriver = curGDriverIndx * nTumors;

        // loop through each GT in the tumor
        for (unsigned int gt = 0; gt < nGT; gt++)
        {   
            // we use a binary tree to keep the statistics
            float T[2] = {0.0};
            float TE[4] = {0.0};
            float TD[4] = {0.0};
            float TDE[8] = {0.0};            
            
            int curGTIndx = tumorGtIndices[gt];
            int gtRowStart = curGTIndx * nTumors;

            //Check if current gt is the same as GD, if yes, switch GD 
            if (curGTIndx == curGDriverIndx){
                unsigned int GD2Indx = tumorGlobDriverIndices2[ge];
                rowStartForGlobDriver = GD2Indx * nTumors;
            }

            float gtPrior;
            if (priorMatrix[gtRowStart + tumorID] == 0)
            {
                gtPrior = -FLT_MAX;
            }
            else
            {
                gtPrior = log(priorMatrix[gtRowStart + tumorID]);
            }

            for(int t = 0; t < nTumors; t++)
            {
                int tVal = gtDataMatrix[gtRowStart + t];
                int eVal = geDataMatrix[rowStartForGE + t];
                int dVal = gtDataMatrix[rowStartForGlobDriver + t];

                //TDE[0xTDE] T is the gt value, D is the global driver value, E is the ge value
                //e.g. TDE[7] means TDE[0x110] save the count when T=1 and D=1 and E=1
                TDE[tVal*4+dVal*2+eVal]++;
            }

            //TD[0xTD] T is the gt value, D is the global driver value
            //e.g. TD[2] means TD[ox10] save the count when T=1 D=0
            TD[0] = TDE[0] + TDE[1]; //T0D0 = T0D0E0 + T0D0E1 
            TD[1] = TDE[2] + TDE[3]; //T0D1 = T0D1E0 + T0D1E1
            TD[2] = TDE[4] + TDE[5]; //T1D0 = T1D0E0 + T1D0E1 
            TD[3] = TDE[6] + TDE[7]; //T0D1 = T1D1E0 + T1D1E1 
            //TE[0xTE]] T is the gt value, E is the ge value
            //e.g. TE[3] means TE[0x11] save the count when T=1 and E=1
            TE[0] = TDE[0] + TDE[2]; //T0E0 = T0D0E0 + T0D1E0
            TE[1] = TDE[1] + TDE[3]; //T0E1 = T0D0E1 + T0D1E1 
            TE[2] = TDE[4] + TDE[6]; //T1E0 = T1D0E0 + T1D1E0
            TE[3] = TDE[5] + TDE[7]; //T1E1 = T1D0E1 + T1D1E1
            //T[0xT] T is the gt value
            //e.g. T[1] save the count when gt value T = 1 
            T[0] = TE[0] + TE[1]; //T0 = T0E0 + T0E1
            T[1] = TE[2] + TE[3]; //T1 = T1E0 + T1E1  
            
            //Therr is no count for T0ge0, T0ge1 and T0
            TE[0]=TE[1] = 0.0;
            T[0] = 0.0;                    
                    
            float TFscore;
            if(curGTIndx == 0)
            {
                //TFscore = calcA0Fscore(T1,  T1ge1, T1ge0, T0,  T0ge1, T0ge0);
                TFscore = calcA0Fscore(T[1],  TE[3], TE[2], T[0],  TE[1], TE[0]);
            }
            else 
            {
                //TFscore = calcFscore( T1,  T1ge1, T1ge0, T0,  T0ge1, T0ge0 );
                TFscore = calcFscore( T[1],  TE[3], TE[2], T[0],  TE[1], TE[0] );
            }

            //float DFscore = calcFscore( D1, D1ge1, D1ge0, D0, D0ge1, D0ge0 );
            float DFscore = calcFscore( TD[1], TDE[3], TDE[2], TD[0], TDE[1], TDE[0] );

            float lnData = TFscore + DFscore + gtPrior;

            tumorPosteriorMatrix[gt * nGE + ge] = lnData;

            float pGT1GE1, pGT0GE1;
            if(gt == 0)
            {
                pGT1GE1 = (ALPHANULL + TE[3]) / (ALPHANULL + ALPHANULL + T[1]);
                pGT0GE1 = (ALPHANULL + TDE[1] + TDE[3]) / (ALPHANULL + ALPHANULL + nTumors - T[1]);
            }
            else
            {
                pGT1GE1 = (ALPHAIJK11 + TE[3]) / (ALPHAIJK11 + ALPHAIJK10 + T[1]);
                pGT0GE1 = (ALPHAIJK01 + TDE[1] + TDE[3]) / (ALPHAIJK01 + ALPHAIJK00 + nTumors - T[1]);                      
            }

            if(pGT1GE1 <= pGT0GE1)
            {
                tumorPosteriorMatrix[gt* nGE + ge] = -FLT_MAX;
            }

            // restore GD after processing current gt == GD
            if (curGTIndx == curGDriverIndx)
                rowStartForGlobDriver = curGDriverIndx * nTumors;
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
        
        // finished populating a column of GTs with respect to a given GE, normalize so that the column sum to 1
        for (unsigned int gt = 0; gt < nGT; gt++)
            tumorPosteriorMatrix[gt * nGE + ge] = exp(tumorPosteriorMatrix[gt * nGE + ge] - normalizer);  
        
    }
  
    // save results to file
     string outFileName;
    if (*outPath.end() != '/')
    {
        outFileName = outPath + "/" +  curTumorName + ".csv";
    }
    else
    {
        outFileName = outPath + curTumorName + ".csv";
    }

    //ofstream file;
    ofstream outFile;
    try
    {
        outFile.open(outFileName.c_str());
    }
    catch(ofstream::failure e)
    {
        cout << "Exception opening output file. Please ensure you have an existing directory for file.\n";
    }
    
    //start writing CSV representation of TDIMatrix    
    //write column headers
    for(int i = 0; i < nGE; i++)
    {
        outFile << "," << geNames[i];
    }
    outFile << "\n";
    
    for(int i = 0; i < nGT; i++)
    {
        outFile << gtNames[i];
        for(int j = 0; j < nGE; j++)
        {
            outFile << "," << tumorPosteriorMatrix[i * nGE + j];
        }
        outFile << "\n";
    }    
    outFile.close();

    delete [] tumorPosteriorMatrix;
}

/**
 * This function calculate marginal likelihood using TDI algorithm.  It calculate the causal score for all 
 * GT-vs-GE pairs observed in a given tumor, populate a GT-by-GE score matrix and output to file.  
 * @param gtMatrix     An TDIMatrix reference, in which SGA data of a collection of tumors are stored
 *                              in a tumor-by-gene matrix
 * @param geMatrix        An TDIMatrix reference, in which DEG data of a collection of tumors are stored
 *                              in a tumor-by-gene matrix
 * @param mapGlobDrivers        A reference to a dictionary of global drivers of DEGs in the "geMatrix",
                                Each entry is keyed by the DEG gene name, and the value is a reference of a 
                                vector<string> containing the two 2 global driver.  

 */
void TDIC_marginal(GTMatrix& gtMatrix, TDIMatrix& geMatrix, map<string, 
        vector<string> >& mapGlobDrivers, const string outPath, const float v0)
{
    // initializations 
   
    bool* gtDataMatrix = gtMatrix.getMatPtr();
    bool* geDataMatrix = geMatrix.getMatPtr();    
    
     // Allocate memory for nGT x nGE matrix
    unsigned int nGT = gtMatrix.getNGenes();
    unsigned int nGE = geMatrix.getNGenes();
    int nTumors = gtMatrix.getNTumors();
    if (nTumors != geMatrix.getNTumors()) // number of tumors should agree between gtMatrix and geMatrix
    {
        cerr << "Bug: gtMatrix and geMatrix contain different number of tumors ";
    }
    float* priorMatrix = gtMatrix.getPriorPtr();

    /*Get names of all GE and GTs     */    
    vector<string> geNames = geMatrix.getGeneNames();    
    vector<string> gtNames = gtMatrix.getGeneNames();

    // Populate the global drivers corresponding to the GEs
    vector<int> tumorGeIndices;
    for (int i = 0; i < nGE; i++)
    {
        tumorGeIndices.push_back(i);
    }

    // Declare 2 vectors to hold the 1st and 2nd global drivers for DEGs
    vector<int> tumorGlobDriverIndices1, tumorGlobDriverIndices2;
    //map <string, string> globalDriversMap;
    if (!getDEGGlobDriverIndices(gtMatrix, geMatrix, mapGlobDrivers, tumorGeIndices, tumorGlobDriverIndices1, tumorGlobDriverIndices2))
    {
        cout << "Error occurred when retrieving global drivers";
    }
   
    // allocate memory to hold the matrix
    float* tumorPosteriorMatrix = new float[nGT * nGE]();    

    // loop through each GE
    #pragma omp parallel for
    for(unsigned int ge = 0; ge < nGE; ge++)
    {
        float normalizer = 0;
        unsigned int curGeIndx = ge;
        unsigned int rowStartForGE = curGeIndx * nTumors; 
        
        // find the globDriver for this give ge   
        unsigned int curGDriverIndx = tumorGlobDriverIndices1[ge]; //curGDriverIndx is found by ge indx     
        unsigned int rowStartForGlobDriver = curGDriverIndx * nTumors;
        
        // loop through each GT in the tumor
        for (unsigned int gt = 0; gt < nGT; gt++)
        { 
            //Check if current gt is the same as GD, if yes, switch GD
            unsigned int oldRowStartForGlobDriver = rowStartForGlobDriver; 
            if (gt == curGDriverIndx){
                unsigned int GD2Indx = tumorGlobDriverIndices2[ge];
                rowStartForGlobDriver = GD2Indx * nTumors;
            }

            // we use a binary tree to keep the statistics
            float T[2] = {0.0};
            float TE[4] = {0.0};
            float TD[4] = {0.0};
            float TDE[8] = {0.0};            
            
            int curGTIndx = gt;  //This is special case because we are going through all GTs

            int gtRowStart = curGTIndx * nTumors;

            for(int t = 0; t < nTumors; t++)
            {
                int tVal = gtDataMatrix[gtRowStart + t];
                int eVal = geDataMatrix[rowStartForGE + t];
                int dVal = gtDataMatrix[rowStartForGlobDriver + t];

                //TDE[0xTDE] T is the gt value, D is the global driver value, E is the ge value
                //e.g. TDE[7] means TDE[0x110] save the count when T=1 and D=1 and E=1
                TDE[tVal*4+dVal*2+eVal]++;
            }

            //TD[0xTD] T is the gt value, D is the global driver value
            //e.g. TD[2] means TD[ox10] save the count when T=1 D=0
            TD[0] = TDE[0] + TDE[1]; //T0D0 = T0D0E0 + T0D0E1 
            TD[1] = TDE[2] + TDE[3]; //T0D1 = T0D1E0 + T0D1E1
            TD[2] = TDE[4] + TDE[5]; //T1D0 = T1D0E0 + T1D0E1 
            TD[3] = TDE[6] + TDE[7]; //T0D1 = T1D1E0 + T1D1E1 
            //TE[0xTE]] T is the gt value, E is the ge value
            //e.g. TE[3] means TE[0x11] save the count when T=1 and E=1
            TE[0] = TDE[0] + TDE[2]; //T0E0 = T0D0E0 + T0D1E0
            TE[1] = TDE[1] + TDE[3]; //T0E1 = T0D0E1 + T0D1E1 
            TE[2] = TDE[4] + TDE[6]; //T1E0 = T1D0E0 + T1D1E0
            TE[3] = TDE[5] + TDE[7]; //T1E1 = T1D0E1 + T1D1E1
            //T[0xT] T is the gt value
            //e.g. T[1] save the count when gt value T = 1 
            T[0] = TE[0] + TE[1]; //T0 = T0E0 + T0E1
            T[1] = TE[2] + TE[3]; //T1 = T1E0 + T1E1  
            
            //Therr is no count for T0ge0, T0ge1 and T0
            TE[0]=TE[1] = 0.0;
            T[0] = 0.0;                    
                    
            float TFscore;
            if(curGTIndx == 0)
            {
                //TFscore = calcA0Fscore(T1,  T1ge1, T1ge0, T0,  T0ge1, T0ge0);
                TFscore = calcA0Fscore(T[1],  TE[3], TE[2], T[0],  TE[1], TE[0]);
            }
            else 
            {
                //TFscore = calcFscore( T1,  T1ge1, T1ge0, T0,  T0ge1, T0ge0 );
                TFscore = calcFscore( T[1],  TE[3], TE[2], T[0],  TE[1], TE[0] );
            }

            //float DFscore = calcFscore( D1, D1ge1, D1ge0, D0, D0ge1, D0ge0 );
            float DFscore = calcFscore( TD[1], TDE[3], TDE[2], TD[0], TDE[1], TDE[0] );

            float lnData = TFscore + DFscore;

            tumorPosteriorMatrix[gt * nGE + ge] = lnData;

            float pGT1GE1, pGT0GE1;
            if(gt == 0)
            {
                pGT1GE1 = (ALPHANULL + TE[3]) / (ALPHANULL + ALPHANULL + T[1]);
                pGT0GE1 = (ALPHANULL + TDE[1] + TDE[3]) / (ALPHANULL + ALPHANULL + nTumors - T[1]);
            }
            else
            {
                pGT1GE1 = (ALPHAIJK11 + TE[3]) / (ALPHAIJK11 + ALPHAIJK10 + T[1]);
                pGT0GE1 = (ALPHAIJK01 + TDE[1] + TDE[3]) / (ALPHAIJK01 + ALPHAIJK00 + nTumors - T[1]);                      
            }

            if(pGT1GE1 <= pGT0GE1)
            {
                tumorPosteriorMatrix[gt* nGE + ge] = -FLT_MAX;
            }

            // restore GD after processing current gt == GD
            if (gt == curGDriverIndx)
                rowStartForGlobDriver = oldRowStartForGlobDriver;
        }        
    }
  
    // save results to file
     string outFileName;
    if (*outPath.end() != '/')
    {
        outFileName = outPath + "/" +  "ICI_discrete_all_marginal" + ".csv";
    }
    else
    {
        outFileName = outPath + "ICI_discrete_all_marginal" + ".csv";
    }

    //ofstream file;
    ofstream outFile;
    try
    {
        outFile.open(outFileName.c_str());
    }
    catch(ofstream::failure e)
    {
        cout << "Exception opening output file. Please ensure you have an existing directory for file.\n";
    }
    
    //start writing CSV representation of TDIMatrix    
    //write column headers
    for(int i = 0; i < nGE; i++)
    {
        outFile << "," << geNames[i];
    }
    outFile << "\n";
    
    for(int i = 0; i < nGT; i++)
    {
        outFile << gtNames[i];
        for(int j = 0; j < nGE; j++)
        {
            outFile << "," << tumorPosteriorMatrix[i * nGE + j];
        }
        outFile << "\n";
    }    
    outFile.close();

    delete [] tumorPosteriorMatrix;
}

/**
 * This function parse the text file that list top 2 global drivers for each of 
 * DEGs observed in a DEG matrix. 
 * @param A string fileName
 * @return A boolean value indicating the success  
 */
bool parseGlobDriverDict(string fileName, map<string, vector<string> > & globDriverMap){
    ifstream inFileStream;
    string line;
    vector<string> fields;  
    vector<string> drivers;  
   
    try {
        inFileStream.open(fileName.c_str()); 
 
        while(getline(inFileStream, line))
        {   
            line.erase(std::remove(line.begin(), line.end(), '\n'), line.end()); 
            line.erase(std::remove(line.begin(), line.end(), '\r'), line.end()); 
            fields = split(line, ',');
            
            drivers.push_back(fields.at(1));
            drivers.push_back(fields.at(2));
            
            globDriverMap.insert(std::pair<string, vector<string> > (fields.at(0), drivers));                
        }
        inFileStream.close();
    }
    catch (ifstream::failure e) {
        cout << "Fail to open file " << fileName;
        return false;
    }   
    return true; 
}

/**
 * Split a string by a given delimiter.
 * @param s String to be split.
 * @param delim Single character delimiter by which to split the string.
 * @param elems List containing each split substring of 's' split by 'delim'.
 * @return 
 */
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

/**
 * This split function calls '&split'. User calls this function.
 * @param s String to be split by 'delim'.
 * @param delim Character delimiter to split the string 's'.
 * @return List of substrings resulting from the split.
 */
std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}


/**
 * 
 * @param gtMat GT matr
 * @param geMat
 * @param mapGlobDrivers
 * @param inDEGIndices
 * @param OutGlobDriverVec
 * @return 
 */

bool getDEGGlobDriverIndices(GTMatrix& gtMat, TDIMatrix& geMat, map<string, vector<string> >& mapGlobDrivers, 
                            vector<int>& inDEGIndices, vector<int>& OutGlobDriverVec1, vector<int>& OutGlobDriverVec2)                            
{
    /*
     * First we must get the names of the DEGs corresponding to the indices in "inDEGIndices".
     * Then, using these DEG names, we can access their global driver through our map "mapGlobDrivers" 
     * and push them onto 'OutGlobDriverVec'.
     */
    //cout << "Inside getDEGGlobDriver.\n";
    vector<string> inDEGNames;
    geMat.getGeneNamesByIndices(inDEGIndices, inDEGNames);

    vector<string> globalDriverNames1, globalDriverNames2;
    for(int i = 0; i < inDEGNames.size(); i++)
    {
        string geneName = inDEGNames[i];
        vector<string> topDrivers = mapGlobDrivers[geneName]; 
        globalDriverNames1.push_back(topDrivers[0]);
        globalDriverNames2.push_back(topDrivers[1]);
    }
    
    gtMat.getGeneIndicesByNames(globalDriverNames1, OutGlobDriverVec1);
    gtMat.getGeneIndicesByNames(globalDriverNames2, OutGlobDriverVec2);
    return true;
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


