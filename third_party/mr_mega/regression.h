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


#pragma once

#include "structures.h"

class lr
{

private:
	matrixD * V;			// Least squares and var/covar matrix
	double RYSQ;			// Multiple correlation coefficient
	double SDV;				// Standard deviation of errors
	double FReg;			// Fisher F statistic for regression
    double VarG;            // Variance of Y for linear regression null likelihood
    double YBAR;            // Defining YBAR
	arrayD * Ycalc;			// Calculated values of Y
	arrayD * DY;			// Residual values of Y
//    lr(const lr& v){};
    void operator=(const lr& v){};
    
	bool isC, isSEC, isV, isYcalc, isDY, isLogistic, isCOV;

		
public:
	arrayD * Cstat;				// Coefficients	
	arrayD * SECstat;			// Std Error of coefficients	
    matrixD * covariance;			//covar matrix  
	bool lr_w(arrayD * Y, matrixD * X, arrayD * W);	//linear regression
	bool lg_w(arrayD * Y, matrixD * X, arrayD * W);  //logistic regression
    double getLnLk(arrayD * Y, matrixD * X, arrayD * W); // Return -2 * sample log-likelihood
    const double nullLikelihood (const vector<double>& COPYpheno, arrayD * W);
    double getFReg(){return FReg;}
    double TSS, TSS0, RSS;
	~lr();
};

