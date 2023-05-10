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

#ifndef structures_h
#define structures_h

#include <vector>
using namespace std;

class matrixD
{
private:
	vector <double> _M;
	int _rows;
	int _cols;

public:
	matrixD(int cRow, int cColumn);

//	~matrixD(void);
	void put(int row, int column, double value);
	double get(int row, int column);
	int size(){return _rows*_cols;}
	int getRows(){return _rows;}
	int getCols(){return _cols;}
    vector <vector <double> > toVector();
    void putFromVector(vector <vector <double> > x);

	bool InvertS(); //invert symmetric matrix
	bool print();
    matrixD correl();
    matrixD distance();
};

class arrayD
{
private:
	vector <double> _M;
	int _rows;

public:
	arrayD(int cRow);
//	~arrayD(void);
	void put(int row,  double value);
	double get(int row);
	int size(){return _rows;}
    bool print();
    vector <double> toVector(){return _M;}
    void putFromVector(vector <double> x);
};

#endif
