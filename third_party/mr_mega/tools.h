/*************************************************************************
 SCOPA software:  March, 2016
 
 Contributors:
 * Andrew P Morris A.P.Morris@liverpool.ac.uk
 * Reedik Magi reedik.magi@ut.ee
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *************************************************************************/

#ifndef _TOOLS_H_
#define _TOOLS_H_

#include <iostream>
#include <map>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <cctype> // std::toupper
#include <string>
#include <algorithm>
class cohort;
class filterD;

using namespace std;

double
chisqr(int Dof, double Cv);

vector <bool>
passFilter(std::vector <std::string> _tokens, cohort & _cohort, vector <filterD> _filters);

bool
setFilters(cohort & _cohort, std::vector <filterD> & _filters);

bool beta2or(double _beta, double _se, double & _or, double & _or95l, double & _or95u);
bool or2beta(double & _beta, double & _se, double _or, double _or95l, double _or95u);


void sortVec(vector <double>& x, int size);
string stoi(int x);
int Tokenize(const  string& str1,
                      vector<string>& tokens,
			const string& delimiters);
string uc(string s);	//uppercase
bool checkAlleles(string & s1, string & s2);	//check if alleles are ok and change numbers to letters if necessary
string flip(string s);	//flip the alleles if 
vector <bool> phenoMasker(int, int);

string  HWE(double aa, double aA, double AA);



#endif
