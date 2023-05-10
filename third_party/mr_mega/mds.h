//
//  mds.hpp
//  MR-MEGA
//
//  Created by reedik on 03/08/2016.
//  Copyright © 2016 Estonian Genome Center. All rights reserved.
//

#ifndef mds_hpp
#define mds_hpp

#include <stdio.h>
#include <vector>
#include "structures.h"

template<class T>
inline const T SQR(const T a) {return a*a;}

template<class T>
inline const T MAX(const T &a, const T &b)
{return b > a ? (b) : (a);}

template<class T>
inline const T MIN(const T &a, const T &b)
{return b < a ? (b) : (a);}

template<class T>
inline const T SIGN(const T &a, const T &b)
{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}



matrixD
calculateMDS(matrixD D, int pcCount);


double
pythag(const double a, const double b);
bool svd(matrixD & u, arrayD &w, matrixD &v);
void tqli(vector<double> &d, vector<double>&e,
          vector<vector<double> > &z);
void tred2(vector<vector<double> > & a,
           vector<double> & d,
           vector<double> &e);
vector<double> eigenvalues(vector<vector<double> > & a);

void svdvar(vector<vector<double> > & v,
            vector<double> & w,
            vector<vector<double> > & cvm);

//bool svdcmp(matrixD & a,
//            arrayD & w,
//            matrixD &v);

bool svdcmp(vector<vector<double> > & a,
            vector<double> & w,
            vector<vector<double> > & v);

#endif /* mds_hpp */
