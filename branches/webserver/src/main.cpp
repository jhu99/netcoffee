/* main.cpp
Author: Jialu Hu
Data: 18.06.2012*/

#include <iostream>
//#include <fstream>
//#include <vector>
#include <omp.h>
//#include <lemon/list_graph.h>
//#include <lemon/smart_graph.h>
//#include <lemon/time_measure.h>
//#include "verbose.h"
//#include "input/networkpool.h"
//#include "input/recordstore.h"
//#include "input/kp_graph.h"
//#include "input/processprofile.h"
//#include "algorithms/solution_x.h"
//#include "algorithms/simulatedannealing.h"
//#include "analyze/semantic_go.h"
//#include "analyze/score_model.h"
//#include "format/format.h"
#include "format/output_html.h"

using namespace std;
/*typedef lemon::ListGraph Graph;
typedef lemon::SmartGraph BpGraph;
typedef NetworkPool<Graph,BpGraph> InputGraph;
typedef KpGraph<InputGraph> InputKpGraph;*/

typedef struct _Option
{
  // Main opitions
  bool out;
  bool create_records;
  bool analyze;
  bool model;
  bool bscore;
  bool tcoffee;
  int task;
  double edgefactor;
  double alpha;
  double beta;
  double eta;
  int numspecies;
  int nmax;
  int numthreads;
  std::string alignmentfile;
  std::string avefunsimfile;
  std::string recordsfile;
  std::string scorefile;/// log-ratio scoring model
  std::string profile;
  std::string logfile;
  std::string formatfile;
  std::string orthologyfile;
  std::string randomfile;
  std::string distributionfile;
  std::string resultfolder;
  std::vector<std::string> blastfiles;
  std::vector<std::string> nullfiles;
  std::vector<std::string> networkfiles;
  std::vector<std::string> associationfiles;/// gene ontology association
  _Option()
  : out(false), create_records(false),
  analyze(false), model(false),bscore(false),task(1),edgefactor(0.1),alpha(0.5),beta(1.0),eta(1.0),numspecies(4),nmax(2000)
  {
    profile="./test_profile.input";
	numthreads=omp_get_max_threads();
	resultfolder="./result/four_species/netcoffee/";
  }
}Option;

typedef Output_html<Option> OutputHtml;
/*
typedef Solution_X<Graph, InputKpGraph, BpGraph, Option> Solution;
typedef RecordStore<InputKpGraph, Option> InputRecord;
typedef SimulatedAnnealing<InputRecord,Solution,InputGraph,Option> SimulatedAnnealingType;
typedef Format<InputGraph, Option> FormatType;
typedef GoList<InputGraph, Option> GoAnalyzer;

bool runParser(Option& myoption)
{
  std::string filename;
  ProcessProfile<Option> myprofile(myoption.profile);
  //myprofile.getOption(myoption);
  filename.append(myoption.resultfolder);
  filename.append(myoption.avefunsimfile);
  myoption.avefunsimfile=filename;
  return true;
}*/

int main()
{
  OutputHtml webpage;

  //string recordsfile = "./dataset/bldata/alignment_records_blast.data";
  //string hmodelfile= "./dataset/homology.model";
  //string nmodelfile= "./dataset/null.model";
  //vector<string> filenames;/// strings for input networks filename
  //string outputfile("./result/alignment.network.data");
  Option myoption;
  //g_verbosity = VERBOSE_DEBUG;
  
  webpage.set_header();
  webpage.get_data(myoption);
  //runParser(myoption);
  std::cout <<"test"<<std::endl;
  webpage.set_footer();

  //InputKpGraph kpgraph(myoption.blastfiles,myoption.numspecies);
  //InputRecord records(hmodelfile,nmodelfile,myoption);
  //std::ofstream mylog(myoption.logfile,std::ios_base::out|std::ios_base::app);
  //InputGraph networks;
  //Timer t;
     
    //t.restart();
    //SimulatedAnnealingType simAnnealing(myoption);
    //networks.initNetworkPool(myoption.networkfiles,myoption.numthreads);

     ////run simulated annealing on multiple networks with tcoffee technique.
    //records.createBpGraphAll(kpgraph,networks);
    //simAnnealing.run_t(kpgraph,records,networks);
    //mylog <<"------------------------NETCOFFEE------------------------------"<<std::endl;
    //mylog <<"It takes "<<t <<" seconds to obtain the alignment with alpha "<<myoption.alpha<<"."<<std::endl;

    //t.stop();
    //if(myoption.out)
    //{
        //simAnnealing.printAlignment_t(myoption.alignmentfile, networks);
    //}
    //mylog.close();
    return 0;
}
