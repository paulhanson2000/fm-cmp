//
//  cohort.cpp
//  MR-MEGA
//
//  Created by reedik on 27/07/2016.
//  Copyright © 2016 Estonian Genome Center. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <vector>
#include "cohort.h"
#include "tools.h"
#include "filter.h"

cohort::cohort()
{
    for (unsigned int i = 0; i<=11;i++)headerPos.push_back(-1);
}


bool
cohort::createHeader(std::string x)
{
    x = uc(x);
    if (!Tokenize(x, header, " "))return false;
    return true;
}


