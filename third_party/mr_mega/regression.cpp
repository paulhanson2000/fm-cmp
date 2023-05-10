/*************************************************************************
 SCOPA software:  March, 2016
 
 Contributors:
 * Andrew P Morris A.P.Morris@liverpool.ac.uk
 * Reedik Magi reedik.magi@ut.eek
    * Reedik Magi reedik@well.ox.ac.uk

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
// Current file contains modified methods implemented in plink (logistic regression) and
// by Walt Fair, Jr. (weighted linear regression)


#include "regression.h"
#include "math.h"
#include <iostream>

double mabs(double x){if (x>0) return x; return x*-1.0;} // absolut value for double


using namespace std;
bool
multip(matrixD & a,
		matrixD * b, matrixD & c)                        // matrix multiplication
{
  int ar = a.getRows();
  int br = b->getRows();
  if (ar == 0 || br == 0)
    return false;
  int ac = a.getCols();
  int bc = b->getCols();
  if ( ac != br )
    return false;
  
//  int cr = ar;
//  int cc = bc;
  
  for (int i = 0; i < ar; i++)
    for (int j = 0; j < bc; j++)
      for (int k = 0; k < ac; k++)
	c.put(i,j, c.get(i,j) + a.get(i,k) * b->get(k,j));

  return true;
}

bool
lr::lr_w(arrayD * Y, matrixD * X, arrayD * W)
{
//////////////////////////////////////////////////////////
// Following code is weighted linear regression model
// This algorithm is implemented by Walt Fair, Jr. 
// (http://www.codeproject.com/KB/recipes/LinReg.aspx)
// and been rewritten in C++ by Reedik Mägi
////////////////////////////////////////////////////////
   // Y->print();
   // X->print();
   // W->print();

    // Y[j]   = j-th observed data point
    // X[i,j] = j-th value of the i-th independent varialble
    // W[j]   = j-th weight value
	isC = false;
	isSEC = false;
	isV = false;
	isYcalc = false;
	isDY = false;
    isLogistic = false;

    int M = Y->size();             // M = Number of data points
    int N = X->size() / M;         // N = Number of linear terms
    int NDF = M - N;			  // Degrees of freedom
    
  //  cout << M  << " " << N << " " << NDF << endl;
    
	Ycalc = new arrayD(M); isYcalc = true;
	DY = new arrayD(M); isDY = true;
    V= new matrixD(N,N); isV = true;
    Cstat= new arrayD(N); isC = true;
    SECstat= new arrayD(N); isSEC = true;
    arrayD B(N);   // Vector for LSQ
    covariance = new matrixD(N,N);
    
    
    // If not enough data, don't attempt regression
    if (NDF < 1)
    {
        return false;
    }


    // Form Least Squares Matrix
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
			V->put(i,j,0);
            for (int k = 0; k < M; k++)
			{
                V->put(i, j, (V->get(i,j) + W->get(k) * X->get(k, i) * X->get(k, j)));

			}
        }
        B.put(i,0);
        for (int k = 0; k < M; k++)
			B.put(i, B.get(i) + W->get(k) * X->get(k, i) * Y->get(k));
    }

   // B.print();
   // V->print();
    
    // V now contains the raw least squares matrix
    if (!V->InvertS())
    {
        /*
        delete V;
        delete Cstat;
        delete SECstat;
        return false;
         */
    }
    // V now contains the inverted least square matrix
    // Matrix multpily to get coefficients C = VB
    for (int i = 0; i < N; i++)
    {
        Cstat->put(i, 0);
        for (int j = 0; j < N; j++)
            Cstat->put(i, (Cstat->get(i) + V->get(i, j) * B.get(j)));
    }
    //V->print();
    //Cstat->print();
    
    // Calculate statistics
    TSS = 0;
    TSS0 = 0;
    RSS = 0;
    YBAR = 0;
    double WSUM = 0;
    for (int k = 0; k < M; k++)
    {
        YBAR = YBAR + W->get(k) * Y->get(k);
        WSUM = WSUM + W->get(k);
    }
    YBAR = YBAR / WSUM;
    
  //  cout << YBAR << endl;
    
    for (int k = 0; k < M; k++)
    {
        Ycalc->put(k, 0);
        for (int i = 0; i < N; i++)
            Ycalc->put(k, Ycalc->get(k) + Cstat->get(i) * X->get(k, i));
        DY->put(k, Ycalc->get(k) - Y->get(k));
        TSS = TSS + W->get(k) * (Y->get(k) - YBAR) * (Y->get(k) - YBAR);
        TSS0 = TSS0 + W->get(k) * (Y->get(k)) * (Y->get(k));
        
        RSS = RSS + W->get(k) * DY->get(k) * DY->get(k);
    }
  //  Ycalc->print();
  //  DY->print();
  //  cout << TSS << " " << RSS << endl;
    
    double SSQ = RSS / NDF;
    RYSQ = 1 - RSS / TSS;
    FReg = 9999999;
    if (RYSQ < 0.9999999)
        FReg = RYSQ / (1 - RYSQ) * NDF / (N - 1);
    SDV = sqrt(SSQ);

    VarG = TSS/(M-1);
    

    
    // Calculate var-covar matrix and std error of coefficients
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            V->put(i, j, V->get(i, j) * SSQ);
            covariance->put(i, j, V->get(i, j));
        }
        // SECstat->put(i, sqrt(V->get(i, i)));
    }
    
    //lets try this
    for (int i = 0; i < N; i++)SECstat->put(i, sqrt(covariance->get(i, i)));

//covariance->print();
    
   // SECstat->print();
    

	return true;
}

bool
lr::lg_w(arrayD * Y, matrixD * X, arrayD * W)
{
    // Y[j]   = j-th observed data point
    // X[i,j] = j-th value of the i-th independent varialble
    // W[j]   = j-th weight value
	isC = false;
	isSEC = false;
	isV = false;
	isYcalc = false;
	isDY = false;
    isLogistic = true;
    
    int M = Y->size();             // M = Number of data points
    int N = X->size() / M;         // N = Number of linear terms
    int NDF = M - N;			  // Degrees of freedom


    // If not enough data, don't attempt regression
    if (NDF < 1)
    {
        return false;
    }

	V= new matrixD(N,N); isV = true;
	Cstat= new arrayD(N); isC = true;
	SECstat= new arrayD(N); isSEC = true;

    
	//START PLINK
	
	arrayD p(M);
	arrayD V1(M);
///////////////////////////////////////
// Newton-Raphson to fit logistic model.
// Following code is edited version of PLINK
// (http://pngu.mgh.harvard.edu/~purcell/plink)
// source code.
///////////////////////////////////////    
    bool converge = false;
    int it = 0;

    while ( ! converge && it < 20) 
      {
	
	// Determine p and V
	for (int i = 0; i < M; i++)
	  {
	    double t = 0;
	    for (int j = 0; j < N; j++)
	      t += Cstat->get(j) *  X->get(i,j);	    
	    p.put(i, 1/(1+exp(-t)));
		V1.put(i,  p.get(i) * (1 - p.get(i)));
	  }
	
	// Update coefficients
	// b <- b +  solve( t(X) %*% V %*% X ) %*% t(X) %*% ( y - p ) 
	
	matrixD T(N,N);

	for (int j = 0; j < N; j++)
	  for (int k = j; k < N; k++) 
	    {
	      double sum = 0;
	      for (int i = 0; i < M; i++)
				sum += X->get(i, j) * V1.get(i) *  X->get(i, k);
	      T.put(j, k, sum);
		  T.put(k, j, sum);	      
	    }
	if (!T.InvertS())
	{
 //       delete V;
 //       delete Cstat;
 //       delete SECstat;
		return false;
	}
	matrixD T2(N,M);
	
	// note implicit transpose of X
	for (int i = 0; i < N; i++)
	  for (int j = 0; j < M; j++)
	    for (int k = 0; k < N; k++)
	      T2.put(i, j, T2.get(i,j) + T.get(i, k) *  X->get(j, k));  
		
	arrayD t3(M);
	for (int i = 0; i < M; i++) 
	  t3.put(i,  (Y->get(i) - p.get(i))*W->get(i));
	
	arrayD ncoef(N);
	for (int j = 0; j < N; j++) 
	  for (int i = 0; i < M; i++) 
	    ncoef.put(j, ncoef.get(j) + T2.get(j, i) * t3.get(i));

	// Update coefficients, and check for 
	// convergence
	double delta = 0;
	for (int j = 0; j < N; j++) 	
	  {
	    delta += mabs(ncoef.get(j));
	    Cstat->put(j, Cstat->get(j) + ncoef.get(j));
	  }

	if ( delta < 1e-6 )
	  converge = true;

	// Next iteration
	it++;
      }
    /////////////////////////////////////////
    // Obtain covariance matrix of estimates
    // S <- solve( t(X) %*% V %*% X )        
    // Transpose X and multiple by diagonal V
    matrixD Xt(N,M);
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++) 
		Xt.put(j, i,  X->get(i,j)* V1.get(i)*W->get(i));
	matrixD V2(N,N);
	if (!multip(Xt, X, V2))
	{
 //       delete V;
 //       delete Cstat;
 //       delete SECstat;
		return false;
	}
	//end PLINK code

	if (!V2.InvertS())
	{
//        delete V;
//        delete Cstat;
//        delete SECstat;
		return false;
	}


    // Calculate var-covar matrix and std error of coefficients
    for (int i = 0; i < N; i++)
    {
        SECstat->put(i, sqrt(V2.get(i, i)));
    }
    covariance = new matrixD(N,N);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            covariance->put(i,j,V2.get(i,j));
        }    
    }
    return true;

}

//getLnLk is implemented by Kanishka Bhattacharya for INTRAPID software

double 
lr::getLnLk(arrayD * Y, matrixD * X, arrayD * W)
{

        // Return -2 * sample log-likelihood
        // We assume the model is fit, and all Y's are either 0 or 1
        double lnlk = 0;
        int nind = Y->size();             // M = Number of data points
        int np = X->size() / nind;         // N = Number of linear terms
      //  cout << "INDS:" << nind << endl;
      //  cout << "LINEAR TERMS:" << np << endl;
    if (isLogistic)
    {
        
        for (int i=0; i<nind; i++)                 
        {
            double t = 0;
            for (int j=0; j<np; j++)       
                t += Cstat->get(j) * X->get(i, j);
            lnlk += Y->get(i) == 1 ? log( 1/(1+exp(-t))) : log(1 - (1/(1+exp(-t))) ); //logistic regression
            
        }                                
        
        return -2 * lnlk;              
    }
    else
    {

        double sigma = SDV;
        
        cout << sigma << endl;
        
        for (int i=0; i<nind; i++)
        {
            double t = 0;
            for (int j=0; j<np; j++)
                t += Cstat->get(j) * X->get(i, j);
            lnlk +=
            -(
               ((Y->get(i)-t)*(Y->get(i)-t))/(2*sigma*sigma)+log(sigma)
//              (W->get(i)*pow(Y->get(i)-t,2))/(2*sigma*sigma)+(W->get(i)*log(sigma))
             );
        }
        return -2 * lnlk;
    }
    return -9999;
}

const double
lr::nullLikelihood (const vector<double>& COPYpheno, arrayD * W)
{
    if (isLogistic)
    {
        double output = 0.0;
        double sumPheno = 0;
        double probSuccess = 0;
        const double vectorSize = COPYpheno.size();
        for (int i=0; i<vectorSize; ++i)
        {
            sumPheno = sumPheno + COPYpheno[i];
        }
        probSuccess = sumPheno/vectorSize ;
        output = sumPheno * log (probSuccess) + (vectorSize - sumPheno) * log (1 - probSuccess); //logistic regression
        return -2 * output;
    }
    else
    {
        
    
        int nind = (unsigned int) COPYpheno.size();
        double lnlk =0;
        double t = YBAR;
        double sigma = sqrt(VarG);

        cout << t << " " << sigma << endl;
        
        for (int i=0; i<nind; i++)
        {
            lnlk +=
            -(
              ((COPYpheno[i]-t)*(COPYpheno[i]-t))/(2*sigma*sigma)+log(sigma)
 //             (W->get(i)*pow(COPYpheno[i]-t,2))/(2*sigma*sigma)+(W->get(i)*log(sigma))
              );
        }
        return -2 * lnlk;
    }
    return -9999;
}


lr::~lr()
{
    if (isV == true)
    {
        delete V;
        V= NULL;
        delete Cstat;
        Cstat = NULL;
        delete SECstat;
        SECstat = NULL;
        delete covariance;      //not yet implemented in linear regression
        covariance = NULL;
        delete Ycalc;  // added by Jani Heikkinen 20/05/2015
        Ycalc = NULL;  // see
        delete DY;     // above
        DY = NULL;     // this as well.
    }
}



