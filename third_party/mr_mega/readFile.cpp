//
//  readFile.cpp
//  MR-MEGA
//
//  Created by reedik on 27/07/2016.
//  Copyright © 2016 Estonian Genome Center. All rights reserved.
//

#include <stdio.h>
#include "readFile.h"
#include "tools.h"

readFile::readFile(std::string fileName)
{
    buffer=new char[1000000];
    eof=false;
    isOK=true;
    if (fileName.substr(fileName.length()-2)=="gz")
    {
        _F1 = gzopen(fileName.c_str(),"r");
        if (_F1==NULL)isOK=false;
    }
    else
    {
        _F2.open (fileName.c_str());
        if (!_F2.is_open())isOK=false;
    }
}

readFile::~readFile()
{
    delete[] buffer;
}

std::string
readFile::getLine(std::string fileName)
{
    
    if (eof==true)return NULL;
    
    if (fileName.substr(fileName.length()-2)=="gz")
    {
        while(0!=gzgets(_F1,buffer,1000000))
        {
            std::string s = buffer;
            s.erase(s.find_last_not_of(" \n\r\t")+1);
            return s;
        }
        eof=true;
    }
    else
    {
        while (! _F2.eof() )
        {
            std::string s;
            getline (_F2,s);
            s.erase(s.find_last_not_of(" \n\r\t")+1);
            return s;
        }
        eof=true;
    }
    return "EOF";
}
