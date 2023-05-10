//
//  filter.h
//  MR-MEGA
//
//  Created by reedik on 08/07/2016.
//  Copyright © 2016 Estonian Genome Center. All rights reserved.
//

#ifndef filter_h
#define filter_h

#include <iostream>
#include <vector>
#include "cohort.h"


class filterD
{
private:
    
    
public:
    filterD(std::string name,std::string equation,double value);
    
    std::string _name;
    std::string _equation;
    double _value;
};




#endif /* filter_h */
