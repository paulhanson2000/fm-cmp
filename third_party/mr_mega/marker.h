//
//  marker.hpp
//  MR-MEGA
//
//  Created by reedik on 28/07/2016.
//  Copyright © 2016 Estonian Genome Center. All rights reserved.
//

#ifndef marker_h
#define marker_h

#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "global.h"
#include "structures.h"

class global;

class marker{
private:
    std::string _name;
    int _pos;
    int _chr;
    unsigned int _cohortCount;
    std::string _ea;
    std::string _nea;
    std::vector <float> _eaf;
    std::vector <int> _n;
    
    std::vector <float> _beta;
    std::vector <float> _se;
    
    
    
    
public:
    marker(std::string, global & G);

    int addPos(std::string, std::ofstream &);
    int addChr(std::string, std::ofstream &);
    int addAlleles(std::string, std::string, std::ofstream &);
    int addBetaSE(int i, float beta,float se, std::ofstream &);
    
    int pushEAF(std::string, std::string,double,int, int, std::ofstream &);

    std::string getName(){return _name;}
    int getPos(){return _pos;}
    int getChr(){return _chr;}
    std::string getEA(){return _ea;}
    std::string getNEA(){return _nea;}
    std::vector <float> getEAF(){return _eaf;}
    double getAverageEAF();
    double getN();
    bool allAboveMAF();
    int getCohortCount();
    std::string getCohortFlags();
    bool getSE(arrayD * W, bool gc, vector <cohort> & cohorts);
    bool getBeta(arrayD * Y);
    bool getPC(matrixD * X, matrixD & PCs, int countPC);
    bool getonlyPC(matrixD * X, matrixD & PCs, int countPC);
    bool getnullPC(matrixD * X, matrixD & PCs, int countPC);

};




#endif /* marker_h */
