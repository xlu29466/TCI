// This version of TDIC reads in tumor-SGA prior matrix. Prior is not calculated from GtMatrix anymore.
// When loading GtMatrix, at the same time load the prior matrix as well

#include <cstdlib>
#include <unistd.h>
#include <getopt.h>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <ctype.h>
#include <fstream>
#include <time.h>
#include "TDIC.h"

//#include "TDIMatrix.h"
//#include "GTMatrix.h"
//#include "PanCanGTMatrix.h"


using namespace std;

int main(int argc, char** argv) {
    
    // parse arguments 
    //extern char *optarg;
    //extern int optind, opterr, optopt;
    
    
    time_t t_start,t_end;
    time (&t_start);
    float v0 = 0.05;    
    int rowStart = -1;
    int rowEnd = -1;
    int hasOpt;
    int nTumors;
    GTMatrix* gtMatrix;
    
    // values of argment 
    string gtFilePath,  globalDriverPath, degFilePath, outPath, priorFilePath, strv0;
    bool outPutMarginal = false;

    while((hasOpt = getopt(argc, argv, "mhs:e:f:d:g:o:p:v:")) != -1)
    {
        switch(hasOpt)
        {
            case 'm' :
                outPutMarginal = true;
                break;

            case 'p':
                priorFilePath = optarg;
                break;
                
            case 'f':
                gtFilePath = optarg;
                break;

            case 'v':
                strv0 = optarg ;
                v0 = atof(strv0.c_str());
                break;
            
            case 'd':
                degFilePath = optarg;
                break;

            case 'g':
                globalDriverPath = optarg ;
                break;
            
            case 's':
                if(atoi(optarg) >= 0)       
                    rowStart = atoi(optarg);
                else
                {
                    cout << "rowStart given is less than zero. Exiting out.\n";
                    exit(1);
                }
                break;
            
            case 'e':
                rowEnd = atoi(optarg);
                if(rowEnd < 0)
                {
                    cout << "rowEnd given is less than zero. Exiting out.\n";
                    exit(1);                    
                }
                break;

            case 'o':
                outPath = optarg;
                break;

            case 'h':
                cerr << "Usage: TDIC -p tumorPriorFile -f inputGaMatrix -d inputGeMatrix -g inputGlobDriverDictionary -o pathForOutputResults [-s rowStart index -e rowEnd index]\n";
                exit(1);
                break;
                
            case '?':
                if(optopt == 'p')
                {
                  cout << "Option -" << optopt << " requires an argument.\n";
                  return 0;
                }
                if(optopt == 'f')
                {
                  cout << "Option -" << optopt << " requires an argument.\n";
                  return 0;
                }
                else if(optopt == 'g')
                {
                  cout << "Option -" << optopt << " requires an argument.\n";
                  return 0;          
                }
                else if(optopt == 'd')
                {
                  cout << "Option -" << optopt << " requires an argument.\n";
                  return 0;          
                }
                else if(isprint(optopt))
                {
                  cout << "Unknown option -" << optopt << ".\n";
                  return 0;
                }
                else
                {
                  cout << "Unknown option character.\n";
                  return 0;
                }

            default:
                cerr << "Usage: TDIC -p tumorPriorFile -f inputGaMatrix -d inputGeMatrix -g inputGlobDriverDictionary -o pathForOutputResults\n";
                abort();
        }
    }

    if (priorFilePath.empty())
    {
        cerr <<"Must input tumor-SGA prior matrix";
        exit(1);
    }
    
    if (gtFilePath.empty() )
     {
        //both gtFileePath and gtcFilePath not exist
        cerr << "Must input GtMatrix \n";
        exit(1);            
    }    
    else 
    {           //read in GT matrices
        cout << "Reading SGA matrix: " << gtFilePath << " and prior file: " << priorFilePath << "\n";
        gtMatrix = new GTMatrix(gtFilePath, priorFilePath);
        nTumors = gtMatrix->getNTumors();
    }
 
    //read in GE matrices       
    cout << "Reading DEG matrix. " << degFilePath << "\n";
    TDIMatrix* geMatrix = new TDIMatrix(degFilePath);
   
    cout << "Reading global driver file.\n";
    map<string, vector<string> > globalDriverMap;
    parseGlobDriverDict(globalDriverPath, globalDriverMap);

    if(rowStart == -1)
        rowStart = 0;
    if(rowEnd == -1 || rowEnd > nTumors)
        rowEnd = nTumors;

	if(rowStart > rowEnd)
    {
        cout << "Given rowEnd index is smaller than given rowStart. Exiting out.\n";
        exit(1);
    }

    cout << "v0 = " << v0 << "\n";
    if (!gtFilePath.empty() && !gtFilePath.empty())
    {
        for(int i = rowStart; i < rowEnd; i++)
        {
            if (i % 50 == 0)
                printf("TDIC processed %d tumors.\n", i);
            TDIC(*gtMatrix, *geMatrix, globalDriverMap, i, outPath, v0);
                    
        }
    }

    //check if the flag for output overall marginal table
    if (outPutMarginal && !gtFilePath.empty() && !gtFilePath.empty())
    {
        TDIC_marginal( *gtMatrix, *geMatrix, globalDriverMap, outPath, v0);
    }

    delete gtMatrix;
    delete geMatrix;  

    time (&t_end);
    long seconds = difftime (t_end,t_start);
    
    int hours, minutes;
 
    minutes = seconds / 60;
    hours = minutes / 60;
    
    cout <<  " Elasped time is  " << hours << " hours " << minutes%60 << " minutes " << seconds%60 << " seconds." << "\n";

    return 0;
}

