//
//  readFile.h
//  MR-MEGA
//
//  Created by reedik on 27/07/2016.
//  Copyright © 2016 Estonian Genome Center. All rights reserved.
//

#ifndef readFile_h
#define readFile_h


#include <zlib.h>
#include <iostream>
#include <fstream>
using namespace std;


class readFile
{
private:
    gzFile _F1;
    ifstream _F2;
    string fileName;
    char *buffer;
    
public:
    readFile(string);
    ~readFile();
    string getLine(string);
    bool eof;
    bool isOK;
};

#endif /* readFile_h */
