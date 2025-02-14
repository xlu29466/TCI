/* 
 * File:   TDIC.h
 * Author: Kevin Lu
 *
 * Created on April 25, 2015, 10:21 PM
 */
    
    
using namespace std;
#include <map>
#include <vector>
#include <string>

//#include "TDIMatrix.h"
#include "GTMatrix.h"

#ifndef TDIC_H
#define	TDIC_H


// define constant hyper parameters
#define ALPHANULL 1.0

#define ALPHAIJK00 2.0
#define ALPHAIJK01 1.0
#define ALPHAIJK10 1.0
#define ALPHAIJK11 2.0


//Function declarations
void TDIC(GTMatrix&, TDIMatrix&, map<string, vector<string> > & , const int, const string , const float );
void TDIC_marginal(GTMatrix& , TDIMatrix& , map<string, vector<string> >& , const string , const float );

//bool parseGlobDriverDict(string fileName, map<string, string> globDriverMap);
bool parseGlobDriverDict(string fileName, map<string, vector<string> > & globDriverMap);
bool getDEGGlobDriverIndices(GTMatrix& gtMat, TDIMatrix& geMat, map<string, vector<string> >& mapGlobDrivers, 
                            vector<int>& inDEGIndices, vector<int>& OutGlobDriverVec1, vector<int>& OutGlobDriverVec2);

//vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
vector<std::string> split(const std::string &s, char delim);

float logSum(float lnx, float lny);
float calcFscore(float gt1,  float gt1ge1, float gt1ge0, float gt0, float gt0ge1, float gt0ge0 );
float calcA0Fscore(float gt1,  float gt1ge1, float gt1ge0, float gt0, float gt0ge1, float gt0ge0 );

#endif	/* TDIC_H */

