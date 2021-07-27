/* 
 * File:   GeneGlobDriver.h
 * Author: Kevin Lu
 *
 * Created on April 25, 2015, 10:21 PM
 */
    
    
using namespace std;
#include <map>
#include <vector>
#include <string>

#include "GTMatrix.h"

#ifndef GENEGLOBDRIVER_H
#define	GENEGLOBDRIVER_H

// define constant hyper parameters
#define ALPHANULL 1.0

#define ALPHAIJK00 2.0
#define ALPHAIJK01 1.0
#define ALPHAIJK10 1.0
#define ALPHAIJK11 2.0

//Function declarations
void GeneGlobDriver(GTMatrix&, TDIMatrix&, const string outPath, const float v0, bool );

//vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
vector<std::string> split(const std::string &s, char delim);

float logSum(float lnx, float lny);
float calcFscore(float gt1,  float gt1ge1, float gt1ge0, float gt0, float gt0ge1, float gt0ge0 );
float calcA0Fscore(float gt1,  float gt1ge1, float gt1ge0, float gt0, float gt0ge1, float gt0ge0 );

#endif	/* GENEGLOBDRIVER_H */

