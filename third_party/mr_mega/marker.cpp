//
//  marker.cpp
//  MR-MEGA
//
//  Created by reedik on 28/07/2016.
//  Copyright © 2016 Estonian Genome Center. All rights reserved.
//

#include "marker.h"
#include <iostream>
#include <vector>
#include <stdlib.h>
#include "global.h"
#include "structures.h"

marker::marker(std::string name, global & G)
{
    _name = name;
    _chr = 0;
    _pos = 0;
    _ea = "N";
    _nea = "N";
    _cohortCount=(unsigned int)G.fileList.size();
    _eaf.resize(_cohortCount, -1);
    _n.resize(_cohortCount, 0);
    _beta.resize(_cohortCount, -666);
    _se.resize(_cohortCount, -666);
}


int
marker::addBetaSE(int i, float beta, float se, std::ofstream & LOG)
{
    if (se<0)return 1;
    _beta[i]=beta;
    _se[i]=se;
    return 0;
}



double
marker::getAverageEAF()
{
    double _sum=0;
    int _count=0;
    for (int i=0; i<_cohortCount;i++)
    {
        if (_eaf[i]!=-1) {_sum+=_eaf[i]*_n[i];_count += _n[i];}
    }
    if (_count)return _sum/_count;
    return 0;
}

double
marker::getN()
{
    double _count=0;
    for (int i=0; i<_cohortCount;i++)
    {
        if (_n[i]>0) {_count += _n[i];}
    }
    return _count;
}


int
marker::getCohortCount()
{
    double _sum=0;
    for (int i=0; i<_cohortCount;i++)
    {
        if (_eaf[i]!=-1) {_sum++;}
    }
    return _sum;
}

std::string
marker::getCohortFlags()
{
    std::string x = "";
    for (int i=0; i<_cohortCount;i++)
    {
        if (_beta[i]==-666) {x.push_back('?');}
        else if (_beta[i]==0) {x.push_back('0');}
        else if (_beta[i]>0) {x.push_back('+');}
        else {x.push_back('-');}
        
    }
    return x;
}

bool
marker::getSE(arrayD * W, bool gc, vector <cohort> & cohorts)
{
    int j=0;
    for (int i=0; i<_cohortCount;i++)
    {
        if (_beta[i]!=-666)
        {
            if (gc && cohorts[i].lambda>1){_se[i]=_se[i]*sqrt(cohorts[i].lambda);}
            W->put(j,pow(1/_se[i],2));j++;
        }
    }
    return 0;
}

bool marker::getBeta(arrayD * Y)
{
    int j=0;
    for (int i=0; i<_cohortCount;i++)
    {
        if (_beta[i]!=-666) {Y->put(j,_beta[i]);j++;}
    }
    return 0;
}

bool marker::getPC(matrixD * X, matrixD & PCs, int countPC)
{
    int j=0;
    for (int i=0; i<_cohortCount;i++)
    {
        if (_beta[i]!=-666)
        {
            X->put(j, 0,1.0);
            for (int k = 0; k < countPC;k++)
            {
                X->put(j, k+1,PCs.get(i,k));
                
            }
            j++;
        }
    }
    return 0;
}

bool marker::getnullPC(matrixD * X, matrixD & PCs, int countPC)
{
    int j=0;
    for (int i=0; i<_cohortCount;i++)
    {
        if (_beta[i]!=-666)
        {
            X->put(j, 0,0.0);
            for (int k = 0; k < countPC;k++)
            {
                X->put(j, k+1,PCs.get(i,k));
                
            }
            j++;
        }
    }
    return 0;
}

bool marker::getonlyPC(matrixD * X, matrixD & PCs, int countPC)
{
    int j=0;
    for (int i=0; i<_cohortCount;i++)
    {
        if (_beta[i]!=-666)
        {
            for (int k = 0; k < countPC;k++)
            {
                X->put(j, k,PCs.get(i,k));
                
            }
            j++;
        }
    }
    return 0;
}


bool
marker::allAboveMAF()
{
    for (int i=0; i<_cohortCount;i++)
    {
        // also takes care of -1's
        if (_eaf[i]<0.01)return false;
        if(_eaf[i]>0.99)return false;
    }
    return true;
}



int
marker::addPos(std::string pos, std::ofstream & LOG)
{
    int mypos = atoi(pos.c_str());
    if (!(mypos>0)){return 2;}
    if (mypos != _pos && _pos !=0){return 1;}
    _pos=mypos;
    return 0;
}


int
marker::addChr(std::string chr, std::ofstream & LOG)
{
    int mychr;
    if (chr=="X")mychr=23;
    else if (chr=="Y")mychr=24;
    else if (chr=="XY")mychr=25;
    else if (chr=="MT")mychr=26;
    else mychr = atoi(chr.c_str());
    if (!(mychr>0 && mychr<26)){return 2;}
    if (mychr != _chr && _chr !=0){return 1;}
    _chr=mychr;
    return 0;
}

int
marker::addAlleles(std::string ea, std::string nea, std::ofstream & LOG)
{
    if (ea==nea)return 1;
    if (_ea=="N"){_ea = ea; _nea=nea;return 0;}
    if (!((_ea == ea && _nea ==nea)||(_ea==nea && _nea ==ea)))return 2;
    return 0;
}


int
marker::pushEAF(std::string ea, std::string nea, double eaf,int n, int i, std::ofstream & LOG)
{
    double myeaf = eaf;
    _n[i]=n;
    if (myeaf<0 || myeaf>1)return 1;
    if (_ea == ea && _nea ==nea){_eaf[i]=myeaf;}
    else if (_ea == nea && _nea == ea){_eaf[i]=(1-myeaf);}
    else return 2;
    return 0;
}