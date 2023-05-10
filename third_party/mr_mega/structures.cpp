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



#include <math.h>
#include <iostream>
#include <vector>
#include "structures.h"
using namespace std;


double myabs(double x){if (x>0) return x; return x*-1.0;}

matrixD::matrixD(int cRow, int cColumn) : _rows(0), _cols(0)
{
//	cout << "Printing array " << cRow*cColumn << endl;
    if (cRow>0 && cColumn>0)
	{
		
		_rows=cRow;
		_cols=cColumn;
		_M.resize(_rows*_cols,0);
	}
}



matrixD
matrixD::correl()
{
    matrixD m(_cols, _cols);
    for (int i = 0 ; i < _cols-1; i++)
    {
        m.put(i,i,1);
        for (int j = i+1; j < _cols; j++)
        {
            double x=0, y=0, x2=0, y2=0, xy=0;
            for (int k = 0; k < _rows; k++)
            {
                x+=get(i,k);
                y+=get(j,k);
                x2+=pow(get(i,k),2);
                y2+=pow(get(j,k),2);
                xy+=get(i,k)*get(j,k);
            }
            double correl = ((_rows * xy) - (x * y))/(sqrt(_rows*pow(x,2)*x2)*sqrt(_rows*pow(y,2)*y2));
            m.put(i,j,correl);
            m.put(j,i,correl);
        }
    }
    return m;
}
matrixD
matrixD::distance()
{
    matrixD m(_cols, _cols);
    for (int i = 0 ; i < _cols-1; i++)
    {
//        m.put(i,i,0);
        for (int j = i+1; j < _cols; j++)
        {
            double x=0;
            for (int k = 0; k < _rows; k++)
            {
                x+=pow(get(k,i)-get(k,j),2);

            }
            m.put(i,j,x/_rows);
            m.put(j,i,x/_rows);
        }
    }

    return m;
}



bool
matrixD::print()
{
//    cout << "Printing matrix " << _rows << " x " << _cols << endl;
	for (int i = 0; i < _rows; i++)
	{
        cout << i;
		for (int j = 0; j < _cols; j++)
            cout << "\t" << _M[j+_cols*i];
			cout << endl;
	}
	return true;
}

bool
arrayD::print()
{
	for (int i = 0; i < _rows; i++)
	{
			cout << i << "\t" <<  _M[i]<< endl;
	}
	return true;
}
//matrixD::~matrixD(void)
//{
//	if (_rows !=0 && _cols!=0)
//	{
//		delete _M;
//	}
//}

void
matrixD::put(int row, int column, double value)
{
	if (row<=_rows && row>=0 && column<=_cols && column>=0 && _rows>0 && _cols>0)
		_M[row*_cols+column] = value;
}

double 
matrixD::get(int row, int column)
{
	if (row<=_rows && row>=0 && column<=_cols && column>=0 && _rows>0 && _cols>0)
		return _M[row*_cols+column];
	else 
		return 0;
}

vector <vector <double> >
matrixD::toVector()
{
    vector <vector <double> > x(_cols, vector<double> (_rows));
    for (int i = 0; i < _cols; i++)
        for (int j = 0; j < _rows; j++)
            x[i][j]=_M[j*_cols+i];

    return x;
}

void
matrixD::putFromVector(vector <vector <double> > x)
{
    for (int i = 0; i < x.size(); i++)
        for (int j = 0; j < x[0].size(); j++)
            _M[j*x.size()+i]= x[i][j];
    _cols=(int) x.size();
    _rows=(int) x[0].size();
}

void
arrayD::putFromVector(vector <double> x)
{
    _M = x;
    _rows = (int) x.size();
}


//Method for inverting symmetric matrix (www.codeproject.com/KB/recipes/LinReg.aspx) by Walt Fair, Jr.

bool 
matrixD::InvertS()
{
	if (_rows < 1 || _cols < 1 || _rows != _cols) return false;
	int n = _rows;
	double * t = new double [n];
	double * Q = new double [n];
	double * R = new double [n];
	double ab;
	int k, l, m;

	for (m = 0; m < n; m++) 
		R[m]=1;
	k = 0;
	for (m = 0; m < n; m++)
	{
		double big = 0;
		for (l=0; l < n; l++)
		{
			ab = myabs(_M[l*_rows+l]);
			if (ab > big && R[l] !=0)
			{
				big = ab;
				k = l;
			}
		}
		if (big == 0)
        {
            delete [] t;
            delete [] Q;
            delete [] R;
            return false;
        }
		R[k] = 0;
		Q[k] = 1 / _M[k*_rows+k];
		t[k] = 1;
		_M[k*_rows+k] = 0;
		if (k !=0)
		{
			for (l = 0; l < k; l++)
			{
				t[l] = _M[l*_rows+k];
				if (R[l] == 0)
					Q[l] = _M[l*_rows+k] * Q[k];
				else
					Q[l] = -_M[l*_rows+k] * Q[k];
				_M[l*_rows+k]=0;
			}
		}
		if ((k + 1) < n)
		{
			for (l = k + 1; l < n; l++)
			{
				if (R[l] !=0)
					t[l] = _M[k*_rows+l];
				else
					t[l] = -_M[k*_rows+l];
				Q[l] = -_M[k*_rows+l] * Q[k];
				_M[k*_rows+l]=0;
			}
		}
		for (l = 0; l < n; l++)
			for (k = l; k < n ; k++)
				_M[l*_rows+k] = _M[l*_rows+k] + t[l] * Q[k];
	}
	m = n;
	l = n - 1;
	for (k = 1; k < n; k++) 
	{
		m = m - 1;
		l = l - 1;
		for (int j = 0; j <= l; j++)
			_M[m*_rows+j] = _M[j*_rows+m];
	}
    delete [] t;
    delete [] Q;
    delete [] R;
	return true;
}



arrayD::arrayD(int cRow) : _rows(0)
{
	if (cRow>0)
	{
		_rows=cRow;
		for (int i = 0; i < cRow; i++) _M.push_back(0);
	}
}
//arrayD::~arrayD(void)
//{
//	if (_rows>0)
//	{
//		delete _M;
//	}
//}

void
arrayD::put(int row, double value)
{
	if (row<=_rows && row>=0 && _rows>0)
		_M[row] = value;
}

double 
arrayD::get(int row)
{
	if (row<=_rows && row>=0 && _rows>0)
		return _M[row];
	else 
		return 0;
}

