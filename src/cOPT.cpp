#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "CenData.h"
#include "Region.h"
#include "WKM.h"
#include "OPT.h"
#include "PHI.h"
#include "MyMath.h"
#include <Rcpp.h>

double pINF = 1.0;

using namespace std;
using namespace Rcpp;

bool parse_roi(char *_str,Region *roi)
{
  char str[10240];
  strcpy(str,_str);
  char *p[32];
  int np = 1, n = strlen(str);
  p[0] = str;
  for( int i=0 ; i<n ; i++ ) {
    if( str[i] == '-' || str[i] == ',' ) { str[i] = '\0'; p[np] = &str[i+1]; np++; }
  }
  if( np%2 > 0 ) return false;
  int nd = np/2;

  region_t r;
  for( int i=0 ; i<nd ; i++ ) {
    r.begin = atof(p[2*i]);
    r.end = atof(p[2*i+1]);
    roi->push_back(r);
  }

  return true;
}

// [[Rcpp::export]]
NumericVector cOPT(NumericMatrix SurvD)
{
  // default parameters
  //char *in_file = NULL;
  char *out_file = NULL;
  char *phi_file = NULL;
  char *prior_file = NULL;
  char *roi_str = NULL;
  int max_depth = 5;
  double min_num_points = 1.0;
  int max_iter = 5;
  double min_error = 0.01;
  double rho = 0.5;
  double alpha = 0.5;
  bool int_part = false;
  bool do_once = false;
  bool mixed_mode = false;
  bool verb_mode = false;

  string out_file_s = "output.surv.txt";
  out_file = &out_file_s[0];
  opterr = 0;

  // read data
  CenData Data;
  if( !Data.readData(SurvD) ) { exit(0); }

  // region of interests
  Region ROI, *pROI = NULL;
  if( roi_str != NULL ) {
    if( !parse_roi(roi_str,&ROI) ) {
      cerr << "Invalid ROI format" << endl;
    }
    if( (int)(ROI.size()) != Data.getDimension() ) {
      cerr << "The dimension of the given ROI does not agree with the data dimension" << endl;
      exit(0);
    }
    pROI = &ROI;
  }

  // OPT
  OPT OPT(&Data,pROI);
  OPT.setMaxDepth(max_depth);
  OPT.setMinNumPoints(min_num_points);
  OPT.setMaxIter(max_iter);
  OPT.setMinErr(min_error);
  OPT.setRho(rho);
  OPT.setAlpha(alpha);
  OPT.setIntPart(int_part);
  OPT.setMixedMode(mixed_mode);
  if( verb_mode ) OPT.printOptOptions();

  // prior density
  Density Sp, *pSp = NULL;
  if( prior_file && do_once ) {
    Sp.readDensityFromFile(prior_file);
    pSp = &Sp;
  }

  // run opt
  if( do_once ) OPT.runOptOne(pSp);
  else OPT.runOpt();

  // output
  Density *Den = OPT.getDensity();
  if( out_file ) Den->writeDensityToFile(out_file);
  else Den->printDensity();

  // optional PHI table output
  PHI *Phi = OPT.getPHI();
  if( phi_file ) Phi->writePhiToFile(phi_file);


  NumericVector res = 0;
  return res;
}



