//
//  filter.cpp
//  MR-MEGA
//
//  Created by reedik on 27/07/2016.
//  Copyright © 2016 Estonian Genome Center. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include "filter.h"

filterD::filterD(std::string name,std::string equation,double value)
{
    _name = name;
    _equation = equation;
    _value = value;
}