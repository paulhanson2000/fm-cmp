//
//  cohort.h
//  MR-MEGA
//
//  Created by reedik on 27/07/2016.
//  Copyright © 2016 Estonian Genome Center. All rights reserved.
//

#ifndef cohort_h
#define cohort_h
#include <map>
#include <vector>
#include <iostream>
#include <stdio.h>
#include "filter.h"
#include "tools.h"


class cohort
{
private:
    
public:
        cohort();
        double lambda;
//        std::map <std::string, double> mafList;
        std::vector <std::string> header;
        std::vector <int> headerPos;
        std::vector <int> filterPos;
        bool createHeader(std::string);
};

#endif /* cohort_h */
