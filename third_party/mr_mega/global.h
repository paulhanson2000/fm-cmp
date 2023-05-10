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

#ifndef __GLOBAL_H__
#define __GLOBAL_H__

#include <string>
#include <vector>
#include <map>
#include <stdlib.h>

#include "cohort.h"
#include "filter.h"
#include "marker.h"
class marker;

class global
{
public:
	global(void);
//	~global(void);
	void createOutput(void);
    bool pushStandards(void);
    bool parseFilters(void);

	std::string version ;

    std::string inputFile;
    std::string outputRoot;
    std::string mapFile;
    
    std::string outputResult;
    std::string outputLog;
    

    bool gc;
    bool gco;
    bool qt;
    bool debugMode;
    bool stdnames;
    bool noalleles;
    bool pcsOnly;
    double threshold;
    bool precalculatedPC;
    
    std::vector <std::string> filterList;
    std::vector <std::string> name_markerList;
    std::vector <std::string> name_eaList;
    std::vector <std::string> name_neaList;
    std::vector <std::string> name_eafList;
    std::vector <std::string> name_betaList;
    std::vector <std::string> name_seList;
    std::vector <std::string> name_orList;
    std::vector <std::string> name_or95lList;
    std::vector <std::string> name_or95uList;
    std::vector <std::string> name_strandList;
    std::vector <std::string> name_nList;
    std::vector <std::string> name_chrList;
    std::vector <std::string> name_posList;
    
    
    std::vector <filterD> filters;
    
    std::vector <std::string> fileList;
    
    std::vector <marker> markerList;
    std::map <std::string, unsigned int> markerNum;
    int pcCount;
    std::vector <cohort> cohorts;
    
    static double CHI_3_P[1447];
};
#endif

