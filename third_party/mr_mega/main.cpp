/*************************************************************************
 MRMEGA software:  July, 2016
 
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
#include <algorithm>
#include <string>
#include <vector>
#include <fstream>
#include <cctype> // std::toupper
#include <map>
#include <zlib.h>
#include "CmdLine.h"
#include "global.h"
#include "tools.h"
#include "structures.h"
#include "regression.h"
#include "studenttdistr.h"
#include "chisquaredistr.h"
#include "filter.h"
#include "readFile.h"
#include "mds.h"



#define LENS 1000000
using namespace TCLAP;
using namespace std;
//double ddabs(double d){if(d<0)return d*-1;return d;}

//bool readSampleFile(global & G, ofstream & LOG);
//bool readGenoFile(global & G, ofstream & LOG);
//bool readExclFile(global & G, ofstream & LOG);
bool readInputFile(global & GLOBAL, ofstream & LOG);
bool calcLambdaGetMAF(int i, global & GLOBAL, ofstream & LOG);
bool readBetaSE(int i, global & GLOBAL, ofstream & LOG);
bool filterOK(vector<bool>_filter){for (unsigned int i=0; i<_filter.size();i++)if(!_filter[i])return false; return true;}
bool headerOK(std::vector <std::string> & _header, int cohN, global & G, ofstream & LOG);
double ddabs(double d){if(d<0)return d*-1;return d;}
void getPCfromInputFiles(global & G, ofstream & L, matrixD & M);

int main (int argc,  char * argv[])
{
    global GLOBAL;


        CmdLine cmd("For more info: http://www.geenivaramu.ee/en/tools/mrmega", ' ', GLOBAL.version );
//FILE NAMES
        ValueArg<string> inputfArg("i","filelist","Specify studies' result files. Default = mrmega.in",false,"mrmega.in","string", cmd);
        ValueArg<string> outfArg("o","out","This specifies output root. By default mrmega",false,"","string", cmd);
        ValueArg<string> mapfArg("m","map","This specifies map file",false,"","string", cmd);
//BOOLEANS
        SwitchArg noallelesArg("", "no_alleles","No allele information has been given. Expecting always the same EA", cmd);
        SwitchArg gcArg("", "gc","Use genomic control correction on input files", cmd);
        SwitchArg gcoArg("", "gco","Use second genomic control correction on output file", cmd);
        SwitchArg qtArg("", "qt","Use this option, if trait is quantitative (columns BETA & SE). Default is binary trait (columns OR, OR95_U, OR_95_L)", cmd);
        SwitchArg debugArg("", "debug","Debug mode on (default OFF)", cmd);
        SwitchArg pcsonlyArg("", "print_pcs_and_die","Print only PCs (default OFF)", cmd);
        SwitchArg stdnamesArg("", "no_std_names","Default column names are not used. All columns must be be defined by user", cmd);
        SwitchArg precalcPCArg("", "precalculated","Using precalculated PCs. PCs must be added into input filelist file starting from second column. The number of PCs must be defined using --pc option (default OFF)", cmd);
        
//p-value print threshold
        ValueArg<double> thresholdArg("t", "threshold", "The p-value threshold for showing direction. Default = 1", false, 1,"double" , cmd);
        ValueArg<int> pcfArg("","pc","This specifies the number of PC to use in regression. Default = 3",false,3,"int", cmd);
        
//LISTS
        MultiArg<string> filterArg("f","filter", "Set a filtering based on column name. It needs 3 arguments: column name, equation [>,<,>=,<=,==,!=], numeric filter value. Multiple filters can be set. Please note that UNIX may require using '\\' before '<' and '>' signs. Column names are not case sensitive. (Example: INFO\\>0.4)", false, "string", cmd);
        MultiArg<string> name_markerArg("","name_marker", "Alternative header to marker name column. Default MARKERNAME", false, "string", cmd);
        MultiArg<string> name_eaArg("","name_ea", "Alternative header to effect allele column. Default EA", false, "string", cmd);
        MultiArg<string> name_neaArg("","name_nea", "Alternative header to other allele column. Default NEA", false, "string", cmd);
        MultiArg<string> name_eafArg("","name_eaf", "Alternative header to effect allele frequency column. Default EAF", false, "string", cmd);
        MultiArg<string> name_betaArg("","name_beta", "Alternative header to effect column. Default BETA", false, "string", cmd);
        MultiArg<string> name_seArg("","name_se", "Alternative header to standard error column. Default SE", false, "string", cmd);
        MultiArg<string> name_orArg("","name_or", "Alternative header to odds ratio column. Default OR", false, "string", cmd);
        MultiArg<string> name_or95lArg("","name_or_95l", "Alternative header to lower 95 CI of odds ratio column. Default OR_95L", false, "string", cmd);
        MultiArg<string> name_or95uArg("","name_or_95u", "Alternative header to upper 95 CI of odds ratio column. Default OR_95U", false, "string", cmd);
        MultiArg<string> name_strandArg("","name_strand", "Alternative header to strand column. Default STRAND", false, "string", cmd);
        MultiArg<string> name_nArg("","name_n", "Alternative header to sample size column. Default N", false, "string", cmd);
        MultiArg<string> name_chrArg("","name_chr", "Alternative header to chromosome column. Default CHROMOSOME", false, "string", cmd);
        MultiArg<string> name_posArg("","name_pos", "Alternative header to position column. Default POSITION", false, "string", cmd);
        
        cmd.parse(argc,argv);
        
        GLOBAL.inputFile = inputfArg.getValue();
        GLOBAL.outputRoot = outfArg.getValue();
        GLOBAL.mapFile = mapfArg.getValue();
        GLOBAL.pcsOnly = pcsonlyArg.getValue();
        GLOBAL.debugMode = debugArg.getValue();
        GLOBAL.gc = gcArg.getValue();
        GLOBAL.gco = gcoArg.getValue();
        GLOBAL.qt = qtArg.getValue();
        GLOBAL.stdnames = stdnamesArg.getValue();
        GLOBAL.noalleles = noallelesArg.getValue();
        GLOBAL.pcCount = pcfArg.getValue();
        GLOBAL.precalculatedPC=precalcPCArg.getValue();

        
        GLOBAL.threshold = thresholdArg.getValue();
        if (GLOBAL.threshold<0 || GLOBAL.threshold>1)
        {
            cout << "P-value threshold out of range. Must be between 0  and 1. Exit program.";
            exit(1);
        }
        GLOBAL.filterList=filterArg.getValue();
        GLOBAL.name_markerList=name_markerArg.getValue();
        GLOBAL.name_eaList=name_eaArg.getValue();
        GLOBAL.name_neaList=name_neaArg.getValue();
        GLOBAL.name_eafList=name_eafArg.getValue();
        GLOBAL.name_betaList=name_betaArg.getValue();
        GLOBAL.name_seList=name_seArg.getValue();
        GLOBAL.name_orList=name_orArg.getValue();
        GLOBAL.name_or95lList=name_or95lArg.getValue();
        GLOBAL.name_or95uList=name_or95uArg.getValue();
        GLOBAL.name_strandList=name_strandArg.getValue();
        GLOBAL.name_nList=name_nArg.getValue();
        GLOBAL.name_chrList=name_chrArg.getValue();
        GLOBAL.name_posList=name_posArg.getValue();
        
        if (!GLOBAL.stdnames){GLOBAL.pushStandards();}
        
        GLOBAL.parseFilters();
        
        GLOBAL.createOutput();
        ofstream LOG (GLOBAL.outputLog.c_str());
        
        cout << "###################\n# MR-MEGA v." << GLOBAL.version << "\n###################\n" << endl;
        LOG << "###################\n# MR-MEGA v." << GLOBAL.version << "\n###################\n" << endl;
        
        cout << "Using following command line options:" << endl;
        LOG << "Using following command line options:" << endl;

        cout << "Input file: " << GLOBAL.inputFile << endl;
        LOG << "Input file: " << GLOBAL.inputFile << endl;
        
        cout << "Output result file: " << GLOBAL.outputResult << endl;
        LOG << "Output result file: " << GLOBAL.outputResult << endl;

        cout << "Output log file: " << GLOBAL.outputLog << endl;
        LOG << "Output log file: " << GLOBAL.outputLog << endl;
        
        cout << "Number of PC-s in regression: " << GLOBAL.pcCount << endl;
        LOG << "Number of PC-s in regression: " << GLOBAL.pcCount << endl;

        if (GLOBAL.precalculatedPC)
        {
            cout << "Program will use precalculated PCs" << endl;
            LOG << "Program will use precalculated PCs" << endl;
        }
        if (GLOBAL.pcsOnly)
        {
            cout << "Program will only print the PCs and die without running meta-analysis!" << endl;
            LOG << "Program will only print the PCs and die without running meta-analysis!" << endl;
        }
        
        if (GLOBAL.gc)
        {
            cout << "Using genomic correction on input files" << endl;
            LOG << "Using genomic correction on input files" << endl;
        }

        if (GLOBAL.gco)
        {
            cout << "Using genomic correction on output file" << endl;
            LOG << "Using genomic correction on output file" << endl;
        }
        if (!GLOBAL.gc && GLOBAL.gco)
        {
            cout << "WARNING: you are using second genomic control only onoutput files and not on the input files, which is not commonly done." << endl;
        }
        if (GLOBAL.qt)
        {
            cout << "Quantitative trait analysis (expecting columns BETA and SE)" << endl;
            LOG << "Quantitative trait analysis (expecting columns BETA and SE)" << endl;
        }
        else
        {
            cout << "Binary trait analysis (expecting columns OR, OR_95L, OR_95U)" << endl;
            LOG << "Binary trait analysis (expecting columns OR, OR_95L, OR_95U)" << endl;
        }
        if (GLOBAL.stdnames)
        {
            cout << "Using default column names has been switched off" << endl;
            LOG << "Using default column names has been switched off" << endl;
        }
        if (GLOBAL.noalleles)
        {
            cout << "Not using allele info and expecting all effects to be measured using the same reference allele" << endl;
            LOG << "Not using allele info and expecting all effects to be measured using the same reference allele" << endl;
        }
        cout << "P-value threshold for showing cohort effect direction: " << GLOBAL.threshold << endl;
        LOG << "P-value threshold for showing cohort effect direction: " << GLOBAL.threshold << endl;
        
        if (GLOBAL.filters.size()==0)
        {
            cout << "No column filters set" << endl;
            LOG << "No column filters set" << endl;
        }
        else
        {
            cout << "Filters set:" << endl;
            LOG << "Filters set:" << endl;
            for (int i = 0;i<GLOBAL.filters.size(); i++)
            {
                cout << "\t" << GLOBAL.filters[i]._name << " " << GLOBAL.filters[i]._equation << " " << GLOBAL.filters[i]._value << endl;
                LOG << "\t" << GLOBAL.filters[i]._name << " " << GLOBAL.filters[i]._equation << " " << GLOBAL.filters[i]._value << endl;
            }
        }
        cout << "Column names:" << endl;
        LOG << "Column names:" << endl;

        cout << "\tMarker name:"; LOG << "\tMarker name:";
        for (int i = 0;i<GLOBAL.name_markerList.size(); i++){cout << " " << uc(GLOBAL.name_markerList[i]);LOG << " " << uc(GLOBAL.name_markerList[i]);}
        cout << endl; LOG << endl;
        
        if (!GLOBAL.noalleles)
        {
            cout << "\tEffect allele:"; LOG << "\tEffect allele:";
            for (int i = 0;i<GLOBAL.name_eaList.size(); i++){cout << " " << uc(GLOBAL.name_eaList[i]);LOG << " " << uc(GLOBAL.name_eaList[i]);}
            cout << endl; LOG << endl;
            
            cout << "\tOther allele:"; LOG << "\tOther allele:";
            for (int i = 0;i<GLOBAL.name_neaList.size(); i++){cout << " " << uc(GLOBAL.name_neaList[i]);LOG << " " << uc(GLOBAL.name_neaList[i]);}
            cout << endl; LOG << endl;
        }
        
        if (GLOBAL.name_eafList.size()>0)
        {
            cout << "\tEffect allele frequency:"; LOG << "\tEffect allele frequency:";
            for (int i = 0;i<GLOBAL.name_eafList.size(); i++){cout << " " << uc(GLOBAL.name_eafList[i]);LOG << " " << uc(GLOBAL.name_eafList[i]);}
            cout << endl; LOG << endl;
        }
        if (GLOBAL.qt)
        {
            cout << "\tEffect (beta):"; LOG << "\tEffect (beta):";
            for (int i = 0;i<GLOBAL.name_betaList.size(); i++){cout << " " << uc(GLOBAL.name_betaList[i]);LOG << " " << uc(GLOBAL.name_betaList[i]);}
            cout << endl; LOG << endl;
            cout << "\tStd. error of effect:"; LOG << "\tStd. error of effect:";
            for (int i = 0;i<GLOBAL.name_seList.size(); i++){cout << " " << uc(GLOBAL.name_seList[i]);LOG << " " << uc(GLOBAL.name_seList[i]);}
            cout << endl; LOG << endl;
            
        }
        else
        {
            cout << "\tEffect (OR):"; LOG << "\tEffect (OR):";
            for (int i = 0;i<GLOBAL.name_orList.size(); i++){cout << " " << uc(GLOBAL.name_orList[i]);LOG << " " << uc(GLOBAL.name_orList[i]);}
            cout << endl; LOG << endl;

            cout << "\tUpper CI of effect:"; LOG << "\tUpper CI of effect:";
            for (int i = 0;i<GLOBAL.name_or95uList.size(); i++){cout << " " << uc(GLOBAL.name_or95uList[i]);LOG << " " << uc(GLOBAL.name_or95uList[i]);}
            cout << endl; LOG << endl;
            
            cout << "\tLower CI of effect:"; LOG << "\tLower CI of effect:";
            for (int i = 0;i<GLOBAL.name_or95lList.size(); i++){cout << " " << uc(GLOBAL.name_or95lList[i]);LOG << " " << uc(GLOBAL.name_or95lList[i]);}
            cout << endl; LOG << endl;
        }
        
        if (GLOBAL.name_strandList.size()>0)
        {
            cout << "\tStrand:"; LOG << "\tStrand:";
            for (int i = 0;i<GLOBAL.name_strandList.size(); i++){cout << " " << uc(GLOBAL.name_strandList[i]);LOG << " " << uc(GLOBAL.name_strandList[i]);}
            cout << endl; LOG << endl;
        }
        
        if (GLOBAL.name_nList.size()>0)
        {
            cout << "\tSample size:"; LOG << "\tSample size:";
            for (int i = 0;i<GLOBAL.name_nList.size(); i++){cout << " " << uc(GLOBAL.name_nList[i]);LOG << " " << uc(GLOBAL.name_nList[i]);}
            cout << endl; LOG << endl;
        }
        cout << "\tChromosome:"; LOG << "\tChromosome:";
        for (int i = 0;i<GLOBAL.name_chrList.size(); i++){cout << " " << uc(GLOBAL.name_chrList[i]);LOG << " " << uc(GLOBAL.name_chrList[i]);}
        cout << endl; LOG << endl;
        cout << "\tPosition:"; LOG << "\tPosition:";
        for (int i = 0;i<GLOBAL.name_posList.size(); i++){cout << " " << uc(GLOBAL.name_posList[i]);LOG << " " << uc(GLOBAL.name_posList[i]);}
        cout << endl; LOG << endl;
        
        if (GLOBAL.debugMode) {LOG << "DEBUG MODE ON"<<endl;}

        cout << "Reading cohort list file..." << endl; 
        readInputFile(GLOBAL, LOG);
        
        cout << "Processing files..."<<endl;
        for (unsigned int i = 0; i<GLOBAL.fileList.size();i++)
        {
            cohort k;
            GLOBAL.cohorts.push_back(k);
            calcLambdaGetMAF(i, GLOBAL, LOG);
        }
    
        matrixD _pcs(int(GLOBAL.fileList.size()), GLOBAL.pcCount);
        if (!GLOBAL.precalculatedPC)
        {
            cout << "Calculating the number of markers present in all cohorts and with MAF>1%..." << endl;
            int _goodMarkers=0;
            matrixD _matrix(30,300);
            for (int i = 0; i< GLOBAL.markerList.size();i++)
            {
                if (GLOBAL.markerList[i].allAboveMAF() && GLOBAL.markerList[i].getChr()>0 && GLOBAL.markerList[i].getPos()>0)
                {
                    _matrix.put(GLOBAL.markerList[i].getChr(), (int)GLOBAL.markerList[i].getPos()/1000000, GLOBAL.markerNum[GLOBAL.markerList[i].getName()]-1);
                    _goodMarkers++;
                }
            }
            cout << "Altogether " << _goodMarkers << " good markers." << endl;
            vector <int> _selectedMarkers;
            for (int i=1; i<=23;i++)
            {
                for (int j=0; j<299; j++)
                {
                    if (_matrix.get(i,j)>0)
                    {
                        _selectedMarkers.push_back(_matrix.get(i,j));
                    }
                }
            }
            cout << "Selected " << _selectedMarkers.size() << " independent variants for EAF correlation calculation" << endl;
        
            cout << "Calculating distance matrix..." << endl;
            matrixD _forDistance(int(_selectedMarkers.size()), int(GLOBAL.fileList.size()));
            for (int i = 0; i < GLOBAL.fileList.size(); i++)
            {
                for (int j = 0; j < _selectedMarkers.size(); j++) _forDistance.put(j,i,GLOBAL.markerList[_selectedMarkers[j]].getEAF()[i]);
            }
            if (GLOBAL.debugMode)
            {
                LOG << "Selected markers: " << endl;
                for (int j = 0; j < _selectedMarkers.size(); j++) LOG << "\t" << GLOBAL.markerList[_selectedMarkers[j]].getName() << endl;

                LOG << "eaf values:" << endl;
                _forDistance.print();
                LOG << endl << endl;
            }
            matrixD _dist = _forDistance.distance();
            if (GLOBAL.debugMode)
            {
                LOG << "Distance matrix between cohorts:" << endl;
                _dist.print();
                LOG << endl << endl;
            }
            
            
            
            cout << "Calculating MDS between cohorts..." << endl;
            _pcs = calculateMDS(_dist, GLOBAL.pcCount);
        }
    else
    {
        getPCfromInputFiles(GLOBAL, LOG, _pcs);
    }
    
        LOG << "Principal components:" << endl;
    LOG << "PCs";
    for (int i = 0; i < _pcs.getCols();i++)
    {
        LOG << " PC" << i;
    }
    LOG << endl;
    for (int j = 0; j < _pcs.getRows();j++)
    {
        LOG << GLOBAL.fileList[j];
        for (int i = 0; i < _pcs.getCols();i++)
        {
            LOG << " " << _pcs.get(j, i);
        }
        LOG << endl;
        
    }
        LOG << endl;
        if (GLOBAL.pcsOnly)
        {
            cout << "PCs";
            for (int i = 0; i < _pcs.getCols();i++)
            {
                cout << " PC" << i;
            }
            cout << endl;
            
            for (int j = 0; j < _pcs.getRows();j++)
            {
                cout << GLOBAL.fileList[j];
                for (int i = 0; i < _pcs.getCols();i++)
                {
                    cout << " " << _pcs.get(j, i);
                }
                cout << endl;
                
            }
            cout << endl <<  "Time to die..." << endl;
            exit(0);
        }
        
  /*
        //we have everything now for the calculation
        GLOBAL.markerList.clear();
        GLOBAL.markerNum.clear();
        cout << "Reading effect sizes from files..." << endl;

        for (int i = 0; i < GLOBAL.fileList.size(); i++)
        {
            readBetaSE(i, GLOBAL, LOG);
        }
   */
        //run regressions
        cout << "Preparing output..." << endl;
        ofstream OUT (GLOBAL.outputResult.c_str());
        OUT << "MarkerName\tChromosome\tPosition\tEA\tNEA\tEAF\tNsample\tNcohort\tEffects";
        for (int j = 0; j < GLOBAL.pcCount+1; j++) OUT << "\tbeta_" << j << "\tse_" << j;
        OUT << "\tchisq_association\tndf_association\tP-value_association\tchisq_ancestry_het\tndf_ancestry_het\tP-value_ancestry_het\tchisq_residual_het\tndf_residual_het\tP-value_residual_het\tlnBF\tComments" <<endl;
        
        string _tmp = "";
        for (int j = 0; j < GLOBAL.pcCount+1; j++) _tmp = _tmp + "\tNA\tNA";
        vector <string> _output1;       //BETAS & SEs
        _output1.resize(GLOBAL.markerList.size(), _tmp);
        
        _tmp = "\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
        vector <string> _output2;     //ancestry & het results
        _output2.resize(GLOBAL.markerList.size(), _tmp);
        
        vector <string> _remarks;     //remarks
        _remarks.resize(GLOBAL.markerList.size(), "\tNA");
        
        vector<double> _chisq;        //assoc chi
        _chisq.resize(GLOBAL.markerList.size());
        
        vector<double> _good_chisq;        //assoc chi for gco

        
        for (int i=0; i<GLOBAL.markerList.size();i++)
        {
            
 /*           OUT << GLOBAL.markerList[i].getName() << "\t" << GLOBAL.markerList[i].getChr() << "\t" << GLOBAL.markerList[i].getPos() << "\t" << GLOBAL.markerList[i].getEA();
            OUT << "\t" << GLOBAL.markerList[i].getNEA() << "\t" << GLOBAL.markerList[i].getAverageEAF() << "\t" << GLOBAL.markerList[i].getN() << "\t"<< GLOBAL.markerList[i].getCohortCount() << "\t" << GLOBAL.markerList[i].getCohortFlags();
   */
            if (GLOBAL.markerList[i].getCohortCount()-2>GLOBAL.pcCount)
            {
            
                arrayD *Y = new arrayD(GLOBAL.markerList[i].getCohortCount());
                matrixD *X2 = new matrixD(GLOBAL.markerList[i].getCohortCount(), GLOBAL.pcCount+1);
                    arrayD *W = new arrayD(GLOBAL.markerList[i].getCohortCount());
                GLOBAL.markerList[i].getSE(W, GLOBAL.gc, GLOBAL.cohorts);
                GLOBAL.markerList[i].getBeta(Y);
                GLOBAL.markerList[i].getPC(X2, _pcs, GLOBAL.pcCount);
  //              GLOBAL.markerList[i].getnullPC(X3, _pcs, GLOBAL.pcCount);
                vector<double> copyPheno = Y->toVector();
                lr LR2;

                    double TSS_RSS=0;
                    double TSS0_RSS=0;
                    double RSS=0;

                    
//                    if (LR3.lr_w(Y, X3, W))
//                    {
//                        nullinterceptLogLikelihood = -LR3.getLnLk(Y, X3, W)/2;
                        
                        
                    
                    if (LR2.lr_w(Y, X2, W))
                    {
                        TSS_RSS = LR2.TSS-LR2.RSS;
                        TSS0_RSS = LR2.TSS0-LR2.RSS;
                        RSS = LR2.RSS;

                        


  //                      double _pModelTest;
  //                      if (ddabs(TSS0_RSS)>0){_pModelTest = chisqr(GLOBAL.pcCount+1,ddabs(TSS0_RSS));}
  //                      else {_pModelTest = NAN;}
                        
                    
                        double _pModelHet;
                        if (ddabs(TSS_RSS)>0){_pModelHet = chisqr(GLOBAL.pcCount,ddabs(TSS_RSS));}
                        else {_pModelHet = NAN;}
                        
                        
      

                        double _pResidHet;
                        int df = GLOBAL.markerList[i].getCohortCount() - GLOBAL.pcCount - 1;
                        if (ddabs(RSS)>0 && df>0){_pResidHet = chisqr(df,ddabs(RSS));}
                        else {_pResidHet = NAN;}
                        
                       // double _BayesFactor = ((LR2.TSS0 - RSS)-(GLOBAL.pcCount+1)*log(GLOBAL.markerList[i].getCohortCount()))/2;
                       // not sure if LR2.TSS0 - RSS is wrong - must be chsq of association through
                        double _BayesFactor = ((ddabs(TSS0_RSS))-((GLOBAL.pcCount+1)*log(GLOBAL.markerList[i].getCohortCount())))/2;
                        
                        
                        std::ostringstream _tmp2;
                        for (int j = 0; j < GLOBAL.pcCount+1; j++){ _tmp2 <<  "\t" << LR2.Cstat->get(j) << "\t" << LR2.SECstat->get(j);
                        }
                        _output1[i]= _tmp2.str();
                        
                        _chisq[i] = ddabs(TSS0_RSS);
                        _good_chisq.push_back(ddabs(TSS0_RSS));
                        std::ostringstream _tmp3;
                        _tmp3 << "\t" << ddabs(TSS_RSS)
                        << "\t" << GLOBAL.pcCount << "\t" <<  _pModelHet << "\t" << ddabs(RSS) << "\t" << df
                        << "\t" << _pResidHet << "\t" << _BayesFactor;
                        _output2[i]= _tmp3.str();

                    }
//                 }
                
                    else
                    {
                        _remarks[i] = "\tcollinear";
                    }

                X2=0;
                delete X2;
                Y=0;
                delete Y;
                W=0;
                delete W;
            }
            else
            {
                _remarks[i] = "\tSmallCohortCount";
            }
        }
        
        
        double _gco_lambda = 1;
        
        if (GLOBAL.gco)
        {
            cout << "Calculating gco lambda: ";
            sortVec(_good_chisq,(int) _good_chisq.size());
            double _median = 0;
            if ((_good_chisq.size())%2!=0) _median = _good_chisq[((_good_chisq.size())/2)-1];
            else _median = (_good_chisq[(_good_chisq.size())/2 -1] + _good_chisq[(_good_chisq.size())/2])/2;
            _gco_lambda = _median/invchisquaredistribution(GLOBAL.pcCount+1,0.5);
            cout << _gco_lambda << endl;
            LOG << "GCO Lambda:" <<  _gco_lambda << endl;
        }
        

        
        for (int i=0; i<GLOBAL.markerList.size();i++)
        {
            
            OUT << GLOBAL.markerList[i].getName() << "\t" << GLOBAL.markerList[i].getChr() << "\t" << GLOBAL.markerList[i].getPos() << "\t" << GLOBAL.markerList[i].getEA();
            OUT << "\t" << GLOBAL.markerList[i].getNEA() << "\t" << GLOBAL.markerList[i].getAverageEAF() << "\t" << GLOBAL.markerList[i].getN() << "\t"<< GLOBAL.markerList[i].getCohortCount() << "\t" << GLOBAL.markerList[i].getCohortFlags();
            
            
            
            double _pModelTest;
            if (ddabs(_chisq[i])>0){_pModelTest = chisqr(GLOBAL.pcCount+1,_chisq[i]/_gco_lambda);}
            else {_pModelTest = NAN;}

            
            
            
            OUT << _output1[i];
            if (_remarks[i]=="\tNA") OUT << "\t" << _chisq[i]/_gco_lambda << "\t" << GLOBAL.pcCount+1 << "\t" << _pModelTest;
            else OUT << "\tNA\tNA\tNA";
            OUT << _output2[i];
            OUT << _remarks[i] << endl;

        }
        
        
        
        LOG << "Analysis finished." <<endl;
        cout << "Analysis finished." <<endl;
        
        
        
        

    return 0;
}
bool readInputFile(global & G, ofstream & L)
{
    char *buffer = new char[LENS];
    if (G.inputFile.substr(G.inputFile.length()-2)=="gz")
    {
        //      	ifstream F (G.inputGenFile.c_str());
        gzFile F =gzopen(G.inputFile.c_str(),"r");
        L << "Cohorts list:" << endl;
        while(0!=gzgets(F,buffer,LENS))
        {
            string line;
            vector<string> tokens;
            string currentmarker = "";
            if (Tokenize(buffer, tokens, " ")>0)
            {
                L << tokens[0] << endl;
                G.fileList.push_back(tokens[0]);
            }
        }
    }
    else
    {
        ifstream F (G.inputFile.c_str());

        if (F.is_open())
        {
            L << "Cohorts list:" << endl;
            while (! F.eof() )
            {
                string line;
                vector<string> tokens;
                getline (F,line);
                string currentmarker = "";
                if (Tokenize(line, tokens, " ")>0)
                {
                    L << tokens[0] << endl;
                    G.fileList.push_back(tokens[0]);
                }

            }
        }
    }
    return true;
}

void getPCfromInputFiles(global & G, ofstream & L, matrixD & _output)
{
    char *buffer = new char[LENS];
    if (G.inputFile.substr(G.inputFile.length()-2)=="gz")
    {
        //          ifstream F (G.inputGenFile.c_str());
        gzFile F =gzopen(G.inputFile.c_str(),"r");
        L << "Reading PCs from input file." << endl;
        int j = 0;
        while(0!=gzgets(F,buffer,LENS))
        {
            string line;
            vector<string> tokens;
            string currentmarker = "";
            if (Tokenize(buffer, tokens, " ")>=G.pcCount+1)
            {
                for (int i=1;i<G.pcCount+1;i++)_output.put(j,i-1,atof(tokens[i].c_str()));
            }
            j++;
        }
    }
    else
    {
        ifstream F (G.inputFile.c_str());

        if (F.is_open())
        {
            L << "Reading PCs from input file." << endl;
            int j = 0;
            while (! F.eof() )
            {
                string line;
                vector<string> tokens;
                getline (F,line);
                string currentmarker = "";
                if (Tokenize(line, tokens, " ")>=G.pcCount+1)
                {
                    for (int i=1;i<G.pcCount+1;i++)
                    {
                        _output.put(j,i-1,atof(tokens[i].c_str()));
                    }
                }
                j++;
 
            }
        }
    }

}


bool calcLambdaGetMAF(int i, global & G, ofstream & L)
{
    unsigned int _lineNr=0;
    unsigned int _goodMarker=0;
    
    std::vector <int> _errorM;
    std::vector <std::string> _errorMlog;
    for (unsigned j = 0; j<9;j++){_errorM.push_back(0);_errorMlog.push_back("");}
    std::vector <int> _filterM;
    for (unsigned j = 0; j<G.filterList.size();j++){_filterM.push_back(0);}
    std::vector <double> _chi;
    cout << "Reading file " << G.fileList[i]<< endl;
    readFile FILE (G.fileList[i]);
    if (FILE.isOK)
    {
        if (G.debugMode) L << "File is ok" << endl;
        map <std::string, unsigned int> _duplicates;
         while(!FILE.eof)
        {
            if (_lineNr==0)
            {
                if (G.debugMode) L << "Reading first line" << endl;

                G.cohorts[i].createHeader(FILE.getLine(G.fileList[i]));
                if (!headerOK( G.cohorts[i].header, i, G,  L))
                {cerr << "Header of " << G.fileList[i] << " is missing mandatory columns. Exit program."<<endl; exit(1);}
                setFilters(G.cohorts[i],G.filters);
            }
            else
            {
                vector <std::string> tokens;
                Tokenize(FILE.getLine(G.fileList[i]),tokens," ");
                if (tokens.size()>1)
                {
                vector <bool> _filter = passFilter(tokens, G.cohorts[i], G.filters);
                if (filterOK(_filter))
                {
                    //get marker from line
                    std::string _name, _chr, _pos, _ea, _nea, _eaf, _beta, _se, _or, _or95u, _or95l, _n;
                    if (G.cohorts[i].headerPos[0]>=0 && G.cohorts[i].headerPos[0]<=tokens.size()) {_name = tokens[G.cohorts[i].headerPos[0]];}
                    else {cerr << "Error reading marker name from the file " << G.fileList[i] << " at line " << _lineNr+1 << ". Exit program." << endl; exit(1);}
                    
                    if (G.cohorts[i].headerPos[3]>=0 && G.cohorts[i].headerPos[3]<=tokens.size()) {_eaf = tokens[G.cohorts[i].headerPos[3]];}
                    else {cerr << "Error reading effect allele frequency from the file " << G.fileList[i] << " at line " << _lineNr+1 << ". Exit program." << endl; exit(1);}
                    
                    if (G.cohorts[i].headerPos[9]>=0 && G.cohorts[i].headerPos[9]<=tokens.size()) {_chr = tokens[G.cohorts[i].headerPos[9]];}
                    else {cerr << "Error reading chromosome from the file " << G.fileList[i] << " at line " << _lineNr+1 << ". Exit program." << endl; exit(1);}
                    
                    if (G.cohorts[i].headerPos[10]>=0 && G.cohorts[i].headerPos[10]<=tokens.size()) {_pos = tokens[G.cohorts[i].headerPos[10]];}
                    else {cerr << "Error reading position from the file " << G.fileList[i] << " at line " << _lineNr+1 << ". Exit program." << endl; exit(1);}

                    if (G.cohorts[i].headerPos[11]>=0 && G.cohorts[i].headerPos[11]<=tokens.size()) {_n = tokens[G.cohorts[i].headerPos[11]];}
                    else {cerr << "Error reading sample size from the file " << G.fileList[i] << " at line " << _lineNr+1 << ". Exit program." << endl; exit(1);}
                    
                    if (!G.noalleles)
                    {
                        if (G.cohorts[i].headerPos[1]>=0 && G.cohorts[i].headerPos[1]<=tokens.size()) {_ea = tokens[G.cohorts[i].headerPos[1]];}
                        else {cerr << "Error reading effect allele from the file " << G.fileList[i] << " at line " << _lineNr+1 << ". Exit program." << endl; exit(1);}
                        
                        if (G.cohorts[i].headerPos[2]>=0 && G.cohorts[i].headerPos[2]<=tokens.size()) {_nea = tokens[G.cohorts[i].headerPos[2]];}
                        else {cerr << "Error reading other allele from the file " << G.fileList[i] << " at line " << _lineNr+1 << ". Exit program." << endl; exit(1);}
                    }
                    
                    if (G.qt)
                    {
                        if (G.cohorts[i].headerPos[4]>=0 && G.cohorts[i].headerPos[4]<=tokens.size()) {_beta = tokens[G.cohorts[i].headerPos[4]];}
                        else {cerr << "Error reading beta from the file " << G.fileList[i] << " at line " << _lineNr+1 << ". Exit program." << endl; exit(1);}
                        
                        if (G.cohorts[i].headerPos[5]>=0 && G.cohorts[i].headerPos[5]<=tokens.size()) {_se = tokens[G.cohorts[i].headerPos[5]];}
                        else {cerr << "Error reading std. err from the file " << G.fileList[i] << " at line " << _lineNr+1 << ". Exit program." << endl; exit(1);}
                    }
                    else
                    {
                        if (G.cohorts[i].headerPos[6]>=0 && G.cohorts[i].headerPos[6]<=tokens.size()) {_or = tokens[G.cohorts[i].headerPos[6]];}
                        else {cerr << "Error reading odds ratio from the file " << G.fileList[i] << " at line " << _lineNr+1 << ". Exit program." << endl; exit(1);}
                        
                        if (G.cohorts[i].headerPos[7]>=0 && G.cohorts[i].headerPos[7]<=tokens.size()) {_or95u = tokens[G.cohorts[i].headerPos[7]];}
                        else {cerr << "Error reading upper confidence interval for odds ratio from the file " << G.fileList[i] << " at line " << _lineNr+1 << ". Exit program." << endl; exit(1);}
                        
                        if (G.cohorts[i].headerPos[8]>=0 && G.cohorts[i].headerPos[8]<=tokens.size()) {_or95l = tokens[G.cohorts[i].headerPos[8]];}
                        else {cerr << "Error reading lower confidence interval for odds ratio for odds ratio from the file " << G.fileList[i] << " at line " << _lineNr+1 << ". Exit program." << endl; exit(1);}
                    }
                    
                    if (G.markerNum[_name]==0) // new marker and first cohort
                    {
                        marker x(_name, G);
                        int e=0,f=0,g=0,h=0, d=0;
                        double beta,se;
                        if (!G.noalleles)
                        {
                            e = x.addAlleles(_ea,_nea,L);
                            if (e==1)
                            {
                                _errorM[4]++;
                                if (_errorM[4]<=5)_errorMlog[4] += "\tEXAMPLE: Effect allele is the same as other allele: " + _name + " (" + _ea + "," + _nea + ")\n";
                            }
                        }
                        f = x.addChr(_chr, L);
                        if (f==2)
                        {
                                _errorM[3]++;
                                if (_errorM[3]<=5)_errorMlog[3] += "\tEXAMPLE: Unvalid chromosome number: " + _name + " (" + _chr + ")\n";
                        }
                        g = x.addPos(_pos, L);
                        if (g==2)
                        {
                            _errorM[1]++;
                            if (_errorM[1]<=5)_errorMlog[1] += "\tEXAMPLE: Unvalid position: " + _name + " (" + _pos + ")\n";
                        }
                        if (!G.noalleles) h = x.pushEAF(_ea, _nea, atof(_eaf.c_str()), atoi(_n.c_str()),i, L);
                        else h = x.pushEAF("N", "N", atof(_eaf.c_str()), atoi(_n.c_str()),i, L);
                        if (h==1)
                        {
                            _errorM[6]++;
                            if (_errorM[6]<=5)_errorMlog[6] += "\tEXAMPLE: Unvalid effect allele frequency: " + _name + " (" + _eaf + ")\n";
                        }
                        if (!G.qt)
                        {
                            if (!(or2beta(beta,se,atof(_or.c_str()),atof(_or95l.c_str()),atof(_or95u.c_str())))){d=1;}
                            if (x.addBetaSE(i,(float) beta, (float) se, L)){d=1;}
                        }
                        else
                        {
                            beta=atof(_beta.c_str());
                            se=atof(_se.c_str());
                            if (x.addBetaSE(i,(float) beta, (float) se, L)){d=1;}
                        }

                        if (e!=1 && h!=1)
                        {
                            G.markerList.push_back(x);
                            G.markerNum[_name]=(unsigned int) G.markerList.size();
                            _goodMarker++;
                        }
   //                     else x.markerDestroy();
                        
                        
                        if (G.qt)
                        {
                            _chi.push_back(pow(beta/se,2));
                        }
                        else
                        {
                            double beta,se;
                            if (or2beta(beta,se,atof(_or.c_str()),atof(_or95l.c_str()),atof(_or95u.c_str())))
                            {
                            }
                            else
                            {
                                _errorM[8]++;
                                if (_errorM[8]<=5) _errorMlog[8]+= "\tEXAMPLE: Problem with effect size for marker " + _name + "(OR: " + _or + ", CI_95L: " + _or95l + ", CI_95U: " + _or95u + "\n";
                            }
                            _chi.push_back(pow(beta/se,2));
                        }
                        
                    }
                    else if (_duplicates[_name]>0)
                    {
                        //duplicated marker - ignore and count
                        _errorM[7]++;
                        if (_errorM[7]<=5){_errorMlog[7] += "\tEXAMPLE: Duplicated marker: " + _name + "\n";}
                    }
                    else
                    {
                        int e = 0,f = 0,g = 0,h = 0;
                        if (!G.noalleles)
                        {
 //                           string testn = _name;
  //                          int testi  =G.markerNum[_name];
                            e = G.markerList[G.markerNum[_name]-1].addAlleles(_ea,_nea, L);
                        }
                        if (e==1)
                        {
                            _errorM[4]++;
                            if (_errorM[4]<=5)_errorMlog[4] += "\tEXAMPLE: Effect allele is the same as other allele: " + _name + " (" + _ea + "," + _nea + ")\n";
                        }
                        if (e==2)
                        {
                            _errorM[5]++;
                            if (_errorM[5]<=5)_errorMlog[5] += "\tEXAMPLE: Alleles dont match with previous cohorts: " + _name + " (" + _ea + "," + _nea + " vs. (previous) " + G.markerList[G.markerNum[_name]-1].getEA() + "," + G.markerList[G.markerNum[_name]-1].getNEA()+")\n";
                        }
                        f = G.markerList[G.markerNum[_name]-1].addChr(_chr,L);
                        if (f==1)
                        {
                            _errorM[2]++;
                            if (_errorM[2]<=5)_errorMlog[2] += "\tEXAMPLE: Chromosome mismatch with previous cohorts: " + _name + " (" + _chr + " vs. previous " + stoi(G.markerList[G.markerNum[_name]-1].getChr()) + ")\n";
                        }
                        if (f==2)
                        {
                            _errorM[3]++;
                            if (_errorM[3]<=5)_errorMlog[3] += "\tEXAMPLE: Unvalid chromosome number: " + _name + " (" + _chr + ")\n";
                        }
                        g = G.markerList[G.markerNum[_name]-1].addPos(_pos,L);
                        if (g==1)
                        {
                            _errorM[0]++;
                            if (_errorM[0]<=5)_errorMlog[0] += "\tEXAMPLE: Position mismatch with previous cohorts: " + _name + " (" + _pos + " vs. previous " + stoi(G.markerList[G.markerNum[_name]-1].getPos()) + ")\n";
                        }
                        if (g==2)
                        {
                            _errorM[1]++;
                            if (_errorM[1]<=5)_errorMlog[1] += "\tEXAMPLE: Unvalid position: " + _name + " (" + _pos + ")\n";
                        }
                        if (!G.noalleles)
                        {
                            h = G.markerList[G.markerNum[_name]-1].pushEAF(_ea, _nea, atof(_eaf.c_str()), atoi(_n.c_str()),i, L);
                            
                        }
                        else h = G.markerList[G.markerNum[_name]-1].pushEAF("N", "N", atof(_eaf.c_str()), atoi(_n.c_str()),i, L);
                        if (h==1)
                        {
                            _errorM[6]++;
                            if (_errorM[6]<=5)_errorMlog[6] += "\tEXAMPLE: Unvalid effect allele frequency: " + _name + " (" + _eaf + ")\n";
                        }
                        double beta, se;
                        if (!G.qt)
                        {
                            or2beta(beta,se,atof(_or.c_str()),atof(_or95l.c_str()),atof(_or95u.c_str()));
//                            if (G.gc && G.cohorts[i].lambda>1){se=se*sqrt(G.cohorts[i].lambda);}
                        }
                        else
                        {
                            beta=atof(_beta.c_str());
                            se=atof(_se.c_str());
//                            if (G.gc && G.cohorts[i].lambda>1){se=se*sqrt(G.cohorts[i].lambda);}
                        }
                        
                        if (e==0 && h==0)
                        {

                                if (_ea == G.markerList[G.markerNum[_name]-1].getEA()) G.markerList[G.markerNum[_name]-1].addBetaSE(i,(float) beta, (float) se, L);
                                else G.markerList[G.markerNum[_name]-1].addBetaSE(i,(float) -1 * beta, (float) se, L);

                                _chi.push_back(pow(beta/se,2));
                                _goodMarker++;
                        }
                        
                        
                        
                    }
                    

                    
                }
                else
                {
                    //create filter problem counts
                    for (unsigned j = 0; j<G.filterList.size();j++){if (!_filter[j])_filterM[j]++;}
                    
                }
                }
            }
            _lineNr++;
        }

    }
    else
    {
        if (G.debugMode) L << "Error reading file" << endl;
        cerr << "Error reading file: " << G.fileList[i] << ". Exit program!"; exit(1);
        return false;
    }
    
    // should i remove markers not present in following cohorts? it will screw up the positions
    
    
    cout << "\tMarker problems:" << endl;
    for (unsigned j = 0; j<9;j++)
    {
        if (_errorM[j]>5)
        {
            cout << _errorMlog[j] << "\tAltogether " << _errorM[j] << " similar errors" << endl;
        }
    }
    
    
    cout << "\tMarker filtering:"<< endl;
    for (unsigned j = 0; j<G.filterList.size();j++)
    {
        cout << "\t" << G.filterList[j] << ": "  <<_filterM[j]<< endl;
    }
    
    if (G.debugMode)
    {
        L << "Chisquare value count " << _chi.size() << endl;
    
    }
    
    cout << "Lambda: ";
    sortVec(_chi,(int) _chi.size());
    double _median = 0;
    if ((_chi.size())%2!=0) _median = _chi[((_chi.size())/2)-1];
    else _median = (_chi[(_chi.size())/2 -1] + _chi[(_chi.size())/2])/2;
    double _lambda = _median/0.4549364;		//median of chi-sq from R ... qchisq(0.5, df= 1)
    cout << _lambda << endl;
    L << "Lambda:" <<  _lambda << endl;
    if (G.gc && _lambda>1){cout << "Using GC correction for this cohort" << endl;L << "Using GC correction for this cohort" << endl;  }
    G.cohorts[i].lambda=_lambda;
    return true;
}

bool headerOK(std::vector <std::string> & _header, int cohN, global & G, ofstream & LOG)
{
    unsigned int _marker=0;
    unsigned int _ea = 0;
    unsigned int _nea = 0;
    unsigned int _eaf = 0;
    unsigned int _beta = 0;
    unsigned int _se = 0;
    unsigned int _or = 0;
    unsigned int _or95u = 0;
    unsigned int _or95l = 0;
    unsigned int _chr = 0;
    unsigned int _pos = 0;
    unsigned int _n = 0;
    for (unsigned int i = 0; i<_header.size();i++)
    {
        for (unsigned int j = 0; j < G.name_markerList.size();j++){if (uc(_header[i])==uc(G.name_markerList[j])){_marker++;G.cohorts[cohN].headerPos[0]=i;}}
        for (unsigned int j = 0; j < G.name_eaList.size();j++){if (uc(_header[i])==uc(G.name_eaList[j])){_ea++;G.cohorts[cohN].headerPos[1]=i;}}
        for (unsigned int j = 0; j < G.name_neaList.size();j++){if (uc(_header[i])==uc(G.name_neaList[j])){_nea++;G.cohorts[cohN].headerPos[2]=i;}}
        for (unsigned int j = 0; j < G.name_eafList.size();j++){if (uc(_header[i])==uc(G.name_eafList[j])){_eaf++;G.cohorts[cohN].headerPos[3]=i;}}
        for (unsigned int j = 0; j < G.name_betaList.size();j++){if (uc(_header[i])==uc(G.name_betaList[j])){_beta++;G.cohorts[cohN].headerPos[4]=i;}}
        for (unsigned int j = 0; j < G.name_seList.size();j++){if (uc(_header[i])==uc(G.name_seList[j])){_se++;G.cohorts[cohN].headerPos[5]=i;}}
        for (unsigned int j = 0; j < G.name_orList.size();j++){if (uc(_header[i])==uc(G.name_orList[j])){_or++;G.cohorts[cohN].headerPos[6]=i;}}
        for (unsigned int j = 0; j < G.name_or95uList.size();j++){if (uc(_header[i])==uc(G.name_or95uList[j])){_or95u++;G.cohorts[cohN].headerPos[7]=i;}}
        for (unsigned int j = 0; j < G.name_or95lList.size();j++){if (uc(_header[i])==uc(G.name_or95lList[j])){_or95l++;G.cohorts[cohN].headerPos[8]=i;}}
        for (unsigned int j = 0; j < G.name_chrList.size();j++){if (uc(_header[i])==uc(G.name_chrList[j])){_chr++;G.cohorts[cohN].headerPos[9]=i;}}
        for (unsigned int j = 0; j < G.name_posList.size();j++){if (uc(_header[i])==uc(G.name_posList[j])){_pos++;G.cohorts[cohN].headerPos[10]=i;}}
        for (unsigned int j = 0; j < G.name_nList.size();j++){if (uc(_header[i])==uc(G.name_nList[j])){_n++;G.cohorts[cohN].headerPos[11]=i;}}
        }
    if (_marker==0){cerr << "MarkerName column is missing. Exit program." << endl;exit(1);}else if (_marker>1){cout << "WARNING: Several columns are labeled as marker name column. Selecting the last one." << endl;LOG << "WARNING: Several columns are labeled as marker name column. Selecting the last one." << endl;}
    
    
    if (!G.noalleles)
    {
        if (_ea==0){cerr << "Effect allele column is missing. If all files have the effects measured from the same allele, please use --no_alleles option. Exit program." << endl;exit(1);}else if (_ea>1){cout << "WARNING: Several columns are labeled as effect allele column. Selecting the last one." << endl;LOG << "WARNING: Several columns are labeled as effect allele column. Selecting the last one." << endl;}
        if (_nea==0){cerr << "Other allele column is missing. If all files have the effects measured from the same allele, please use --no_alleles option. Exit program." << endl;exit(1);}else if (_nea>1){cout << "WARNING: Several columns are labeled as other allele column. Selecting the last one." << endl;LOG << "WARNING: Several columns are labeled as other allele column. Selecting the last one." << endl;}
    }
    
    
    if (_eaf==0){cerr << "Effect allele frequency column is missing. Exit program." << endl;exit(1);}else if (_eaf>1){cout << "WARNING: Several columns are labeled as effect allele frequency column. Selecting the last one." << endl;LOG << "WARNING: Several columns are labeled as effect allele frequency column. Selecting the last one." << endl;}
    
    if (G.qt)
    {
        if (_beta==0){cerr << "Effect size column is missing. Exit program." << endl;exit(1);}else if (_beta>1){cout << "WARNING: Several columns are labeled as effect size column. Selecting the last one." << endl;LOG << "WARNING: Several columns are labeled as effect size column. Selecting the last one." << endl;}
        if (_se==0){cerr << "Std. error column is missing. Exit program." << endl;exit(1);}else if (_se>1){cout << "WARNING: Several columns are labeled as std. error column. Selecting the last one." << endl;LOG << "WARNING: Several columns are labeled as std. error column. Selecting the last one." << endl;}
    }
    else
    {
        if (_or==0){cerr << "Odds ratio column is missing. Exit program." << endl;exit(1);}else if (_or>1){cout << "WARNING: Several columns are labeled as odds ratio column. Selecting the last one." << endl;LOG << "WARNING: Several columns are labeled as odds ratio column. Selecting the last one." << endl;}
        if (_or95u==0){cerr << "Upper confidence interval column is missing. Exit program." << endl;exit(1);}else if (_or95u>1){cout << "WARNING: Several columns are labeled as upper confidence interval column. Selecting the last one." << endl;LOG << "WARNING: Several columns are labeled as upper confidence interval column. Selecting the last one." << endl;}
        if (_or95l==0){cerr << "Lower confidence interval column is missing. Exit program." << endl;exit(1);}else if (_or95l>1){cout << "WARNING: Several columns are labeled as lower confidence interval column. Selecting the last one." << endl;LOG << "WARNING: Several columns are labeled as lower confidence interval column. Selecting the last one." << endl;}
    }
    if (_chr==0){cerr << "Chromosome column is missing. Exit program." << endl;exit(1);}else if (_chr>1){cout << "WARNING: Several columns are labeled as chromosome column. Selecting the last one." << endl;LOG << "WARNING: Several columns are labeled as chromosome column. Selecting the last one." << endl;}
    if (_pos==0){cerr << "Position column is missing. Exit program." << endl;exit(1);}else if (_pos>1){cout << "WARNING: Several columns are labeled as position column. Selecting the last one." << endl;LOG << "WARNING: Several columns are labeled as position column. Selecting the last one." << endl;}
    if (_n==0){cerr << "Sample size column is missing. Exit program." << endl;exit(1);}else if (_pos>1){cout << "WARNING: Several columns are labeled as position column. Selecting the last one." << endl;LOG << "WARNING: Several columns are labeled as position column. Selecting the last one." << endl;}
    return true;
    }





/* MDS TEST DATA
 matrixD testM(13,13);
 testM.put(0,0,0.000000);
 testM.put(0,1,4.473483);
 testM.put(0,2,2.624622);
 testM.put(0,3,3.440069);
 testM.put(0,4,4.316154);
 testM.put(0,5,2.992642);
 testM.put(0,6,101.601725);
 testM.put(0,7,6.145833);
 testM.put(0,8,102.573851);
 testM.put(0,9,1.541608);
 testM.put(0,10,2.318835);
 testM.put(0,11,1.5248490);
 testM.put(0,12,2.0597933);
 testM.put(1,0,4.473483);
 testM.put(1,1,0.000000);
 testM.put(1,2,5.943080);
 testM.put(1,3,4.685793);
 testM.put(1,4,8.515545);
 testM.put(1,5,7.236852);
 testM.put(1,6,105.395296);
 testM.put(1,7,4.058136);
 testM.put(1,8,105.884743);
 testM.put(1,9,4.768057);
 testM.put(1,10,5.458447);
 testM.put(1,11,5.4188154);
 testM.put(1,12,6.2983955);
 testM.put(2,0,2.624622);
 testM.put(2,1,5.943080);
 testM.put(2,2,0.000000);
 testM.put(2,3,5.128179);
 testM.put(2,4,5.101069);
 testM.put(2,5,3.835862);
 testM.put(2,6,110.672061);
 testM.put(2,7,8.110299);
 testM.put(2,8,111.407110);
 testM.put(2,9,1.974524);
 testM.put(2,10,2.719314);
 testM.put(2,11,1.9456327);
 testM.put(2,12,3.0355746);
 testM.put(3,0,3.440069);
 testM.put(3,1,4.685793);
 testM.put(3,2,5.128179);
 testM.put(3,3,0.000000);
 testM.put(3,4,7.221607);
 testM.put(3,5,5.818194);
 testM.put(3,6,106.549569);
 testM.put(3,7,4.288650);
 testM.put(3,8,106.969368);
 testM.put(3,9,3.833129);
 testM.put(3,10,4.192902);
 testM.put(3,11,4.1995164);
 testM.put(3,12,4.9301374);
 testM.put(4,0,4.316154);
 testM.put(4,1,8.515545);
 testM.put(4,2,5.101069);
 testM.put(4,3,7.221607);
 testM.put(4,4,0.000000);
 testM.put(4,5,4.458763);
 testM.put(4,6,112.359114);
 testM.put(4,7,10.427909);
 testM.put(4,8,112.663735);
 testM.put(4,9,4.222986);
 testM.put(4,10,4.925352);
 testM.put(4,11,3.2692857);
 testM.put(4,12,4.0800916);
 testM.put(5,0,2.992642);
 testM.put(5,1,7.236852);
 testM.put(5,2,3.835862);
 testM.put(5,3,5.818194);
 testM.put(5,4,4.458763);
 testM.put(5,5,0.000000);
 testM.put(5,6,110.332137);
 testM.put(5,7,9.491245);
 testM.put(5,8,111.081128);
 testM.put(5,9,3.104981);
 testM.put(5,10,3.521602);
 testM.put(5,11,1.3905075);
 testM.put(5,12,1.5490503);
 testM.put(6,0,101.601725);
 testM.put(6,1,105.395296);
 testM.put(6,2,110.672061);
 testM.put(6,3,106.549569);
 testM.put(6,4,112.359114);
 testM.put(6,5,110.332137);
 testM.put(6,6,0.000000);
 testM.put(6,7,103.291470);
 testM.put(6,8,7.487637);
 testM.put(6,9,107.963120);
 testM.put(6,10,109.310930);
 testM.put(6,11,109.1271259);
 testM.put(6,12,106.5551132);
 testM.put(7,0,6.145833);
 testM.put(7,1,4.058136);
 testM.put(7,2,8.110299);
 testM.put(7,3,4.288650);
 testM.put(7,4,10.427909);
 testM.put(7,5,9.491245);
 testM.put(7,6,103.291470);
 testM.put(7,7,0.000000);
 testM.put(7,8,103.449329);
 testM.put(7,9,6.504795);
 testM.put(7,10,6.828160);
 testM.put(7,11,7.5897106);
 testM.put(7,12,8.4393112);
 testM.put(8,0,102.573851);
 testM.put(8,1,105.884743);
 testM.put(8,2,111.407110);
 testM.put(8,3,106.969368);
 testM.put(8,4,112.663735);
 testM.put(8,5,111.081128);
 testM.put(8,6,7.487637);
 testM.put(8,7,103.449329);
 testM.put(8,8,0.000000);
 testM.put(8,9,108.821772);
 testM.put(8,10,109.834397);
 testM.put(8,11,109.8344609);
 testM.put(8,12,107.3102939);
 testM.put(9,0,1.541608);
 testM.put(9,1,4.768057);
 testM.put(9,2,1.974524);
 testM.put(9,3,3.833129);
 testM.put(9,4,4.222986);
 testM.put(9,5,3.104981);
 testM.put(9,6,107.963120);
 testM.put(9,7,6.504795);
 testM.put(9,8,108.821772);
 testM.put(9,9,0.000000);
 testM.put(9,10,1.756661);
 testM.put(9,11,1.3540491);
 testM.put(9,12,2.3649700);
 testM.put(10,0,2.318835);
 testM.put(10,1,5.458447);
 testM.put(10,2,2.719314);
 testM.put(10,3,4.192902);
 testM.put(10,4,4.925352);
 testM.put(10,5,3.521602);
 testM.put(10,6,109.310930);
 testM.put(10,7,6.828160);
 testM.put(10,8,109.834397);
 testM.put(10,9,1.756661);
 testM.put(10,10,0.000000);
 testM.put(10,11,1.9406277);
 testM.put(10,12,2.8489843);
 testM.put(11,0,1.524849);
 testM.put(11,1,5.418815);
 testM.put(11,2,1.945633);
 testM.put(11,3,4.199516);
 testM.put(11,4,3.269286);
 testM.put(11,5,1.390508);
 testM.put(11,6,109.127126);
 testM.put(11,7,7.589711);
 testM.put(11,8,109.834461);
 testM.put(11,9,1.354049);
 testM.put(11,10,1.940628);
 testM.put(11,11,0.0000000);
 testM.put(11,12,0.8387904);
 testM.put(12,0,2.059793);
 testM.put(12,1,6.298395);
 testM.put(12,2,3.035575);
 testM.put(12,3,4.930137);
 testM.put(12,4,4.080092);
 testM.put(12,5,1.549050);
 testM.put(12,6,106.555113);
 testM.put(12,7,8.439311);
 testM.put(12,8,107.310294);
 testM.put(12,9,2.364970);
 testM.put(12,10,2.848984);
 testM.put(12,11,0.8387904);
 testM.put(12,12,0.0000000);
 //testM.print();
 
 //matrixD pcs = calculateMDS(testM, 3);
 //pcs.print();
 
 //exit(0);
 
 
 */
