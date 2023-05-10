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


#include <iostream>
#include <map>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <cctype> // std::toupper
#include <string>
#include <algorithm>
#include "tools.h"
#include "studenttdistr.h"
#include "chisquaredistr.h"
#include "cohort.h"
#include "filter.h"
#include "gammafunc.h"
#include <iomanip> 

using namespace std;




static long double log_igf(long double S, long double Z);
static long double KM(long double S, long double Z);
double gamma(double N);

#define A 15 // 15

double approx_gamma(double Z)
{
    const double RECIP_E = 0.36787944117144232159552377016147;  // RECIP_E = (E^-1) = (1.0 / E)
    const double TWOPI = 6.283185307179586476925286766559;  // TWOPI = 2.0 * PI
    
    double D = 1.0 / (10.0 * Z);
    D = 1.0 / ((12 * Z) - D);
    D = (D + Z) * RECIP_E;
    D = pow(D, Z);
    D *= sqrt(TWOPI / Z);
    
    return D;
}

long double log_gamma(double N)
{
    /*
     The constant SQRT2PI is defined as sqrt(2.0 * PI);
     For speed the constant is already defined in decimal
     form.  However, if you wish to ensure that you achieve
     maximum precision on your own machine, you can calculate
     it yourself using (sqrt(atan(1.0) * 8.0))
     */
    
    //const long double SQRT2PI = sqrtl(atanl(1.0) * 8.0);
    const long double SQRT2PI = (sqrt(atan(1.0) * 8.0));
    
    long double Z = (long double)N;
    long double Sc;
    
    Sc = (logl(Z + A) * (Z + 0.5)) - (Z + A) - logl(Z);
    
    long double F = 1.0;
    long double Ck;
    long double Sum = SQRT2PI;
    
    
    for(int K = 1; K < A; K++)
    {
        Z++;
        Ck = powl(A - K, K - 0.5);
        Ck *= expl(A - K);
        Ck /= F;
        
        Sum += (Ck / Z);
        
        F *= (-1.0 * K);
    }
    
    return logl(Sum) + Sc;
}


double chisqr(int Dof, double Cv)
{
    //printf("Dof:  %i\n", Dof);
    //printf("Cv:  %f\n", Cv);
    if(Cv < 0 || Dof < 1)
    {
        return 0.0;
    }
    double K = ((double)Dof) * 0.5;
    double X = Cv * 0.5;
    if(Dof == 2)
    {
        return exp(-1.0 * X);
    }
    long double PValue, Gam;
    long double ln_PV;
    ln_PV = log_igf(K, X);
    Gam = log_gamma(K);
    ln_PV -= Gam;
    PValue = 1.0 - expl(ln_PV);
    
    return (double)PValue;
    
}


static long double log_igf(long double S, long double Z)
{
    if(Z < 0.0)
    {
        return 0.0;
    }
    long double Sc, K;
    Sc = (logl(Z) * S) - Z - logl(S);
    
    K = KM(S, Z);
    
    return logl(K) + Sc;
}


static long double KM(long double S, long double Z)
{
    long double Sum = 1.0;
    long double Nom = 1.0;
    long double Denom = 1.0;
    
    for(int I = 0; I < 1000; I++) // Loops for 1000 iterations
    {
        Nom *= Z;
        S++;
        Denom *= S;
        Sum += (Nom / Denom);
    }
    return Sum;
    
}




bool beta2or(double _beta, double _se, double & _or, double & _or95l, double & _or95u)
{
    if (_se<=0)return false;
    _or = exp(_beta);
    _or95l = exp(_beta - 1.96 * _se);
    _or95u = exp(_beta + 1.96 * _se);
    return true;
}
bool or2beta(double & _beta, double & _se, double _or, double _or95l, double _or95u)
{
    if (_or<=0)return false;
    _beta = log(_or);
    _se = ((log(_or)-log(_or95l))/1.96);
    return true;
}



string stoi(int x)
{
    ostringstream convert;
    convert << x;
    return convert.str();
}

vector <bool>
passFilter(std::vector <std::string> _tokens, cohort & _cohort, std::vector <filterD> _filters)
{
    vector <bool> result;
    for (unsigned int i=0; i<_filters.size();i++)result.push_back(true);
    for (unsigned int i=0; i<_filters.size();i++)
    {
        if(_cohort.filterPos[i]!=-1)
        {
            if (_tokens.size()<_cohort.filterPos[i])
            {
                cerr << "Not enough columns in file - must be a broken line. Exit program." << endl; exit(1);
            }
            double value = atof(_tokens[_cohort.filterPos[i]].c_str());
            if (_filters[i]._equation=="eq"){if (value!=_filters[i]._value)result[i]=false;}
            else if (_filters[i]._equation=="ge"){if (value<_filters[i]._value)result[i]=false;}
            else if (_filters[i]._equation=="le"){if (value>_filters[i]._value)result[i]=false;}
            else if (_filters[i]._equation=="gt"){if (value<=_filters[i]._value)result[i]=false;}
            else if (_filters[i]._equation=="lt"){if (value>=_filters[i]._value)result[i]=false;}
            else if (_filters[i]._equation=="ne"){if (value==_filters[i]._value)result[i]=false;}
        }
    }
    return result;
}


bool
setFilters(cohort & _cohort, std::vector <filterD> & _filters)
{
    for (unsigned int i =0; i< _filters.size();i++)
    {
        _cohort.filterPos.push_back(-1);
    }
    for (unsigned int i =0; i < _filters.size(); i++)
    {
        for (unsigned int j=0; j < _cohort.header.size(); j++)
        {
            if (uc(_cohort.header[j])==uc(_filters[i]._name)){_cohort.filterPos[i]=j;}
        }
    }
    
    return true;
}

int Tokenize(const  string& str1,
                      vector<string>& tokens,
			const string& delimiters = " ")
{
    int cnt = 0;
    string emptyStr = "";
    string str = str1;
    std::replace( str.begin(), str.end(), '\r', ' ' );
    std::replace( str.begin(), str.end(), '\t', ' ' );



    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos )
    {
	if (str.substr(lastPos, pos - lastPos) != emptyStr)
	{
        	// Found a token, add it to the vector.
        	tokens.push_back(str.substr(lastPos, pos - lastPos));
        	// Skip delimiters.  Note the "not_of"
        	lastPos = str.find_first_not_of(delimiters, pos);
        	// Find next "non-delimiter"
        	pos = str.find_first_of(delimiters, lastPos);
		//
		cnt++;
	}
    }
    return cnt;
}
void
sortVec(vector <double>& x, int size)
{
	 std::sort(x.begin(), x.end());
}



string uc(string s)
{
  const int length = (int)s.length();
  for(int i=0; i!=length ; ++i)
  {
    s[i] = std::toupper(s[i]);
  }
  return s;
}

bool checkAlleles(string & s1, string & s2)
{
	if (s1.compare("1")==0) s1 = "A";
	if (s1.compare("2")==0) s1 = "C";
	if (s1.compare("3")==0) s1 = "G";
	if (s1.compare("4")==0) s1 = "T";
	if (s2.compare("1")==0) s2 = "A";
	if (s2.compare("2")==0) s2 = "C";
	if (s2.compare("3")==0) s2 = "G";
	if (s2.compare("4")==0) s2 = "T";
	if (!(s1.compare("A")==0 || s1.compare("C")==0 || s1.compare("G")==0 || s1.compare("T")==0))return false;
	if (!(s2.compare("A")==0 || s2.compare("C")==0 || s2.compare("G")==0 || s2.compare("T")==0))return false;
	if (s1==s2)return false;
	return true;
}
string flip(string s)
{
	if (s.compare("A")==0) return "T";
	if (s.compare("C")==0) return "G";
	if (s.compare("G")==0) return "C";
	if (s.compare("T")==0) return "A";
	return "N";
}

vector <bool> phenoMasker(int variant, int size)
{
    vector<bool> X;
    for (int i = 0; i<size; i++){X.push_back(0);}
    int fap = pow(2,size-1);
    int place = 0;
    while(fap>=1)
    {
        if (variant>=fap)
        {
            X[place]=1;
            variant = variant - fap;
        }
        fap = fap/2;
        place ++;        
    }
    return X;
}

string
HWE(double aa, double aA, double AA)
{
    if (aa+aA+AA<=0 || aa<0 || aA<0 || AA<0){return "NA";}
    double sum = aa+aA+AA;
    double a = (2*aa + aA)/(2*sum);
    double eaa = a * a * sum;
    double eaA = 2 * a * (1-a) * sum;
    double eAA = (1-a) * (1-a) * sum;
    double chi = (((aa-eaa)*(aa-eaa))/(eaa))+(((aA-eaA)*(aA-eaA))/(eaA))+(((AA-eAA)*(AA-eAA))/(eAA));
    double p = 1-chisquaredistribution(1,chi);
    std::stringstream strs;
    strs << p;
    //if (p>0 & p<=1)
    return strs.str();
    return "NA";
 }


