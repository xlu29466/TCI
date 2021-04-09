/* 
 * File:   PanCanTDIC.cpp
 *
 * 
 * 
 */

//#include "TDIC.h"
#include "PanCanTDIC.h"
//#include "TDIMatrix.h"
//#include "GTMatrix.h"
//#include "PanCanGTMatrix.h"
#include <fstream>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <float.h>

using namespace std;
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
 * @param tumorID               The tumor to be processed
 */
void PanCanTDIC(PanCanGTMatrix& gtMatrix, TDIMatrix& geMatrix, map<string, 
        string>& mapGlobDrivers, const int tumorID, const string outPath, const float v0){
    // initializations 

    int tumorCanType = gtMatrix.getCanTypeByTumorId(tumorID);
    string curTumorName = gtMatrix.getTumorNameById(tumorID);
   
    // Find the GTs and GEs of the current tumor
    vector<int> tumorGtIndices, tumorGeIndices;
    gtMatrix.findGeneWithOnesInTumor(tumorID, tumorGtIndices);
    geMatrix.findGeneWithOnesInTumor(tumorID, tumorGeIndices);
    
    bool* gtDataMatrix = gtMatrix.getMatPtr();
    bool* geDataMatrix = geMatrix.getMatPtr();

    // Find the global drivers corresponding to the
    vector<int> tumorGlobDriverIndices;
    //map <string, string> globalDriversMap;
    if (!getDEGGlobDriverIndices(gtMatrix, geMatrix, mapGlobDrivers, tumorGeIndices, tumorGlobDriverIndices))
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
  
    vector<float> lntumorMutPriors;
    gtMatrix.calcLnTumorPriors(tumorGtIndices, v0, lntumorMutPriors);

    float* tumorPosteriorMatrix = new float[nGT * nGE]();
    float* tumorPosteriorMatrixC0 = new float[nGT * nGE]();
    float* tumorPosteriorMatrixC1 = new float[nGT * nGE]();
     
    
    // loop through each GE
    #pragma omp parallel for
    for(unsigned int ge = 0; ge < nGE; ge++)
    {

        //nC1 is the number of tumors of the same cancel type
        float nC1 = 0.0;
        unsigned int curGeIndx = tumorGeIndices[ge];
        unsigned int rowStartForGE = curGeIndx * nTumors; 
        
        // find the globDriver for this give ge   
        unsigned int curGDriverIndx = tumorGlobDriverIndices[ge]; //curGDriverIndx is found by ge indx        
        unsigned int rowStartForGlobDriver = curGDriverIndx * nTumors;
        
        // loop through each GT in the tumor
        for (unsigned int gt = 0; gt < nGT; gt++)
        {

            //Variables to save the count considering the cancel type
            float CT[4] = {0.0};
            float CTE[8] = {0.0};
            float CTD[8] = {0.0};
            float CTDE[16] = {0.0};
            
            
            int curGTIndx = tumorGtIndices[gt];

            int gtRowStart = curGTIndx * nTumors;
            
            for(int t = 0; t < nTumors; t++)
            {
                
                int compCanType = gtMatrix.getCanTypeByTumorId(t);
                int theSameCanType = (tumorCanType == compCanType);//if cancer type is the same, the value = 1
                
                int tVal = gtDataMatrix[gtRowStart + t];
                int eVal = geDataMatrix[rowStartForGE + t];
                int dVal = gtDataMatrix[rowStartForGlobDriver + t];
                
//                //Only need to count CTDE, all other variables can be calculated from CTDE
//                
//                //count while NOT considering cancel type
//                T[tVal]++;
//                TE[tVal*2+eVal]++;
//                TD[tVal*2+dVal]++;
//                TDE[tVal*4+dVal*2+eVal]++;
//                
//                //count considering cancel type                
//                CT[theSameCanType*2+tVal]++;
//                CTE[theSameCanType*4+tVal*2+eVal]++;
//                CTD[theSameCanType*4+tVal*2+dVal]++;
                
                
                CTDE[theSameCanType*8+tVal*4+dVal*2+eVal]++;
            }
            
            //Calculate CTE,CTD,TDE from CTDE
            for(int a=0; a<=1; a++)
                for(int b=0; b<=1; b++)
                    for(int c=0; c<=1; c++)
                    {
                        CTE[ a*4 + b*2 + c ] = CTDE[ a*8 + b*4 + 0*2 + c ] + CTDE[ a*8 + b*4 + 1*2 + c ];
                        CTD[ a*4 + b*2 + c ] = CTDE[ a*8 + b*4 + c*2 + 0 ] + CTDE[ a*8 + b*4 + c*2 + 1 ];
                    }
                
            //Calculate CT from CTE
            for(int a=0; a<=1; a++)
                for(int b=0; b<=1; b++)
                    CT[ a*2 + b ] = CTE[ a*4 + b*2 + 0 ] + CTE[ a*4 + b*2 + 1 ];
                            
            //nC1 is the number of tumors of the same cancel type
            nC1 = CT[2]+CT[3];  // nC1= C1T0 + C1T1=CT[2]+CT[3]
            CT[0] = CT[2] = 0.0;
            //There is no count for C0T0E0, C0T0E1 C1TOE0 C1T0E1
            CTE[0] = CTE[1] = CTE[4] = CTE[5] = 0.0;
            
            float TFscore ;
            float DFscore ;
            float pGT1GE1 ;
            float pGT0GE1 ; 
           /*
            *Calculate lnData while cancel types are the same, save the lnData into lnDataC1
            */
 
            float lnDataC1 = 0.0;
            if(curGTIndx == 0)
            {
                TFscore = calcA0Fscore(CT[3],  CTE[7], CTE[6], CT[2],  CTE[5], CTE[4]);
            }
            else 
            {
                //TFscore = calcFscore( T1,  T1ge1, T1ge0, T0,  T0ge1, T0ge0 );
                TFscore = calcFscore( CT[3],  CTE[7], CTE[6], CT[2],  CTE[5], CTE[4] );
            }
            
            DFscore = calcFscore( CTD[5], CTDE[11], CTDE[10], CTD[4], CTDE[9], CTDE[8] );
            
            lnDataC1 = TFscore + DFscore + lntumorMutPriors[gt];

            //tumorPosteriorMatrix[gt * nGE + ge] = lnData;    Comment out for cancel type calculation; Not ready to save yet. Save lnData at the last step
            
            if(gt == 0)
            {
                 pGT1GE1 = (ALPHANULL + CTE[7]) / (ALPHANULL + ALPHANULL + CT[3]);
                pGT0GE1 = (ALPHANULL + CTDE[9] + CTDE[11]) / (ALPHANULL + ALPHANULL + nTumors - CT[3]);
            }
            else
            {
                 pGT1GE1 = (ALPHAIJK11 + CTE[7]) / (ALPHAIJK11 + ALPHAIJK10 + CT[3]);
                pGT0GE1 = (ALPHAIJK01 + CTDE[9] + CTDE[11]) / (ALPHAIJK01 + ALPHAIJK00 + nTumors - CT[3]);         
            }
         
            if(pGT1GE1 <= pGT0GE1)
            {
                //tumorPosteriorMatrix[gt* nGE + ge] = -FLT_MAX;      Changed for cancel type calculation. Save lnData at the last step
                lnDataC1 = -FLT_MAX;
            }
            
            
            /*
            *Calculate lnData while cancel types are NOT the same, save the lnData into lnDataC0
            */
 
            float lnDataC0 = 0.0;
            if(curGTIndx == 0)
            {
                TFscore = calcA0Fscore(CT[1],  CTE[3], CTE[2], CT[0],  CTE[1], CTE[0]);
            }
            else 
            {
                //TFscore = calcFscore( T1,  T1ge1, T1ge0, T0,  T0ge1, T0ge0 );
                TFscore = calcFscore( CT[1],  CTE[3], CTE[2], CT[0],  CTE[1], CTE[0] );
            }

            DFscore = calcFscore( CTD[1], CTDE[3], CTDE[2], CTD[0], CTDE[1], CTDE[0] );
            
                       
            
            lnDataC0 = TFscore + DFscore + lntumorMutPriors[gt];

              //tumorPosteriorMatrix[gt * nGE + ge] = lnData;    Comment out for cancel type calculation; Not ready to save yet. Save lnData at the last step
            if(gt == 0)
            {
                 pGT1GE1 = (ALPHANULL + CTE[3]) / (ALPHANULL + ALPHANULL + CT[1]);
                pGT0GE1 = (ALPHANULL + CTDE[1] + CTDE[3]) / (ALPHANULL + ALPHANULL + nTumors - CT[1]);
            }
            else
            {
                pGT1GE1 = (ALPHAIJK11 + CTE[3]) / (ALPHAIJK11 + ALPHAIJK10 + CT[1]);
                pGT0GE1 = (ALPHAIJK01 + CTDE[1] + CTDE[3]) / (ALPHAIJK01 + ALPHAIJK00 + nTumors - CT[1]);         
            }
         
            if(pGT1GE1 <= pGT0GE1)
            {
                //tumorPosteriorMatrix[gt* nGE + ge] = -FLT_MAX;      Changed for cancel type calculation. Save lnData at the last step
                lnDataC0 = -FLT_MAX;
            }
            
             
            //save lnData lnDataC0 lnDataC1
//            tumorPosteriorMatrix[gt* nGE + ge] = lnData; 
            tumorPosteriorMatrixC0[gt* nGE + ge] = lnDataC0;
            tumorPosteriorMatrixC1[gt* nGE + ge] = lnDataC1;

    
        }
        
        float normalizerC0 = 0;
        float normalizerC1 = 0;
        
        for(unsigned int gt = 0; gt < nGT; gt++)
        {
            if(gt == 0)
            {
                normalizerC0 = tumorPosteriorMatrixC0[gt * nGE + ge];
                normalizerC1 = tumorPosteriorMatrixC1[gt * nGE + ge];
            }
            else
            {
                normalizerC0 = logSum(normalizerC0, tumorPosteriorMatrixC0[gt * nGE + ge]);
                normalizerC1 = logSum(normalizerC1, tumorPosteriorMatrixC1[gt * nGE + ge]);
            }
        }
        
        // finished populating a column of GTs with respect to a given GE, normalize so that the column sum to 1
        for (unsigned int gt = 0; gt < nGT; gt++)
        {
            tumorPosteriorMatrixC0[gt * nGE + ge] = exp(tumorPosteriorMatrixC0[gt * nGE + ge] - normalizerC0);  
            tumorPosteriorMatrixC1[gt * nGE + ge] = exp(tumorPosteriorMatrixC1[gt * nGE + ge] - normalizerC1); 
        }
        
        float normC0;  //cancer types are not the same
        float normC1; //cancer type is the same
        
        float sumC0C1 = 0.0; //for test purpose

        
        for (unsigned int gt = 0; gt < nGT; gt++)
        {
            normC0 = tumorPosteriorMatrixC0[gt * nGE + ge] ;
            normC1 = tumorPosteriorMatrixC1[gt * nGE + ge];
            
            tumorPosteriorMatrix[gt * nGE + ge] = normC0 * (nTumors - nC1)/nTumors + normC1 * nC1/nTumors;
        }    
    }
  

    // save results to file
    //vector<string> gtNames = gtMatrix.getGeneNames();

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
    delete [] tumorPosteriorMatrixC0;
    delete [] tumorPosteriorMatrixC1;
}