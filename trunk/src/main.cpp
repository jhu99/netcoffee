/* main.cpp
Author: Jialu Hu
Data: 18.06.2012*/

#include <iostream>
#include <fstream>
#include <vector>
#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include <lemon/arg_parser.h>
#include <lemon/time_measure.h>
#include "verbose.h"
#include "input/shrinknetwork.h"
#include "input/networkpool.h"
#include "input/recordstore.h"
#include "input/kp_graph.h"
#include "input/processprofile.h"
#include "analyze/mask_whitespace.h"
#include "algorithms/solution_x.h"
#include "algorithms/simulatedannealing.h"
#include "analyze/semantic_go.h"
#include "analyze/score_model.h"
#include "format/format.h"

using namespace std;
typedef lemon::ListGraph Graph;
typedef lemon::SmartGraph BpGraph;
typedef NetworkPool<Graph,BpGraph> InputGraph;
typedef KpGraph<InputGraph> InputKpGraph;

typedef struct _Option
{
  // Main opitions
  bool out;
  bool create_records;
  bool analyze;
  bool model;
  bool blastformat;
  bool tcoffee;
  int task;
  double edgefactor;
  double alpha;
  double beta;
  double threshold;
  int numspecies;
  int nmax;
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
  analyze(false), model(false),blastformat(true),task(0),edgefactor(0.1),alpha(0.5),beta(1.0),threshold(0.4),numspecies(5),nmax(2000)
  {
    profile="./profile.input";
  }
}Option;

typedef Solution_X<Graph, InputKpGraph, BpGraph, Option> Solution;
typedef RecordStore<InputKpGraph, Option> InputRecord;
typedef SimulatedAnnealing<InputRecord,Solution,InputGraph,Option> SimulatedAnnealingType;
typedef Format<InputGraph, Option> FormatType;
typedef GoList<InputGraph, Option> GoAnalyzer;


bool setParser(ArgParser& parser, Option& myoption)
{
  parser
  /// Option for package information
  .boolOption("version","Show the version number.")
  /// Option for four tasks needed to perform on
  .boolOption("records","Generate alignment records file.")
  .boolOption("alignment","Execute the alignment algorithm.")
  .boolOption("analyse","Make analysis on alignmenr result.")
  .boolOption("format","Process input or output file into proper format.")
  .optionGroup("method","version")
  .optionGroup("method","records")
  .optionGroup("method","alignment")
  .optionGroup("method","analyse")
  .optionGroup("method","format")
  .onlyOneGroup("method")
  .mandatoryGroup("method")
  /// Option for variables
  .refOption("task","Complete different tasks. It must be used with -alignment option together.",myoption.task)
  .refOption("alignmentfile","The filename of alignment which is required to either create or analyse.",myoption.alignmentfile)
  .refOption("avefunsimfile", "The filename for functional similarity of alignment.", myoption.avefunsimfile)
  .refOption("recordsfile", "Records file for writing and reading. It is used to store the triplet edges in NetCoffee program.", myoption.recordsfile)
  .refOption("alpha", "Prameter controlling how much sequence similarity contributes. Default is 0.5.", myoption.alpha)
  .refOption("edgefactor", "The factor of the power law normalization. Default is 0.1.", myoption.alpha)
  .refOption("numspecies","Number of the species compared. Default is 4.", myoption.numspecies)
  .refOption("profile","Profile of input parameters.", myoption.profile)
  .refOption("formatfile","Profile of input parameters.", myoption.formatfile)
  .refOption("orthologyfile","Training data for orthology model.", myoption.orthologyfile)
  .refOption("randomfile","Training data for random model.", myoption.randomfile)
  .refOption("distributionfile","Output file for log ratio distribution.", myoption.distributionfile)
  .refOption("resultfolder","The folder which was used as active folder in the data process.", myoption.resultfolder)
  .refOption("beta","Probability for randomly picking out the alignment records. Default is 1.0",myoption.beta)
  .refOption("threshold","Threshold gamma which is used to exclude match-sets with score below a certain value. Default is 0.0.",myoption.threshold)
  .refOption("nmax","The parameter for SA algorithm N.",myoption.nmax)
  .refOption("out","Print the alignment result into file.",myoption.out)
  ;
  return true;
}

bool runParser(ArgParser& myparser, Option& myoption)
{
  std::string filename;
  ProcessProfile<Option> myprofile(myoption.profile);
  myprofile.getOption(myoption);
  myparser.run();
  filename.append(myoption.resultfolder);
  filename.append(myoption.avefunsimfile);
  myoption.avefunsimfile=filename;
  return true;
}

int main(int argc, char** argv)
{
  ArgParser myparser(argc,argv);

  string recordsfile = "./dataset/bldata/alignment_records_blast.data";
  string hmodelfile= "./dataset/homology.model";
  string nmodelfile= "./dataset/null.model";
  //string scorefile= "./dataset/bldata/score_composit.model";/// normalised log-ration likelyhood
  vector<string> filenames;/// strings for input networks filename
  string outputfile("./result/alignment.network.data");
  Option myoption;
  g_verbosity = VERBOSE_DEBUG;

  /// Parse command line.
  setParser(myparser, myoption);
  runParser(myparser, myoption);

  InputKpGraph kpgraph(myoption.blastfiles,myoption.numspecies);
  InputRecord records(hmodelfile,nmodelfile,myoption);
  //InputRecord records(hmodelfile,nmodelfile,myoption.scorefile,myoption.recordsfile,myoption.numspecies,myoption.beta,myoption.alpha);
  std::ofstream mylog(myoption.logfile,std::ios_base::out|std::ios_base::app);
  InputGraph networks;
  Timer t;
  
  /// Print the version of M-NetAligner. 
 if(myparser.given("version"))
  {
    std::cout <<"NetCoffee Version 1.0\
    Author: Jialu Hu\
    email: jialu.hu@fu-berlin.de"<<std::endl;
  }

  if(myparser.given("records"))
  {
    if(myoption.task==0)
    {
      // Read training data for null model and homology model.
        Score_Model mymodel(myoption.orthologyfile,myoption.randomfile,myoption.distributionfile);
        mymodel.run(myoption.blastformat);
    }
    else if(myoption.task==1)
    {
      t.restart();
      networks.initNetworkPool(myoption.networkfiles);
      records.createRecords_MNetAligner(myoption.blastfiles,networks);
      t.stop();
      mylog <<"------------------------M-NETALIGNER-------------------------"<<std::endl;
      mylog <<"It takes "<<t <<" seconds to create alignment records!";
    }
    else if(myoption.task==2)
    {
      records.readRecords2(myoption.recordsfile.c_str());/// Read records for multiple networks.
    }
    else if(myoption.task==4)
    {
      /// Construct triple nodes using tcoffee.
      networks.initNetworkPool(myoption.networkfiles);
      kpgraph.constructGraph();
      //records.createRecords_t(kpgraph,networks);
    }
  }

  if(myparser.given("alignment"))
  {
    t.restart();
    SimulatedAnnealingType simAnnealing(myoption);
    networks.initNetworkPool(myoption.networkfiles);
    if(myoption.task==0)
    {
      // run simulated annealing on multiple networks.
      records.readRecords2(myoption.recordsfile.c_str());
      simAnnealing.run2(records,networks);
      mylog <<"------------------------M-NETALIGNER-------------------------"<<std::endl;
      mylog <<"It takes "<<t <<" seconds to obtain the alignment!"<<std::endl;
    }
    else
    {
      // run simulated annealing on multiple networks with tcoffee technique.
      records.createBpGraphAll(kpgraph,networks);
      simAnnealing.run_t(kpgraph,records,networks);
      mylog <<"------------------------NETCOFFEE------------------------------"<<std::endl;
      mylog <<"It takes "<<t <<" seconds to obtain the alignment with alpha "<<myoption.alpha<<"."<<std::endl;
    }
    t.stop();
    if(myoption.out)
    {
      if(myoption.task==0)
      {
        simAnnealing.printAlignment2(myoption.alignmentfile,records);
      }
      else
      {
        simAnnealing.printAlignment_t(myoption.alignmentfile, networks);
      }
    }
  }

  if(myparser.given("analyse"))
  {
    GoAnalyzer analyzer(myoption);
    analyzer.goInitial();
    if(myoption.task==0)
    {
      // Figure out annotation information of proteins.
      // Figure out coverage performance of the alignment.
      networks.initNetworkPool(myoption.networkfiles);
      std::cout << "Annotation information of proteins in " << myoption.alignmentfile <<std::endl;
      analyzer.readAlignment(myoption.alignmentfile.c_str());/// alignment file must strictly on the format
      analyzer.getAlignmentCoverage(networks);
      //analyzer.getNetworkAnnotation(networks);
    }
    else if(myoption.task==1)
    {
      // Analyze the alignment results in term of semantic similarity with GO.
      // Get average score.
      analyzer.readAlignment(myoption.alignmentfile.c_str());/// alignment file must strictly on the format
      analyzer.convert_fsst();
      analyzer.deleteRedundancy();
      analyzer.getMulFunSim(myoption.avefunsimfile);// avefunsimfile is fsst.result file
    }
    else if(myoption.task==2)
    {
      // Get average score.
      analyzer.readAlignment(myoption.alignmentfile.c_str());/// alignment file must strictly on the format 1
      analyzer.getMulFunSim(myoption.avefunsimfile);// avefunsimfile is fsst.result file
    }
    else if(myoption.task==3)
    {
      // Extract p-value for alignment
      analyzer.extractPValue();
    }
    else{}
  }

  if(myparser.given("format"))
  {
    if(myoption.task==0)
    {
      // discard interactions whose involving proteins donnot appear in SwissUniprot.
      ShrinkNetwork sn(myoption.formatfile);
      sn.run();
    }
    else if(myoption.task==1)
    {
      MaskWhiteSpace maskspace(myoption.resultfolder);
      maskspace.run();
    }
    else if(myoption.task==2)
    {
      // discard bipartite edges whose involving proteins donnot appear in the ppi networks.
      FormatType myformat(myoption);
      networks.initNetworkPool(myoption.networkfiles);
      myformat.removeBiEdges(networks);
    }else if(myoption.task==3)
    {	
      FormatType myformat(myoption);
      std::string outfile("./result/alignment_netcoffee.data");
      myformat.formatAlignment(myoption.alignmentfile,outfile);
    }else
    {
	  // discard these replicate interaction and self-loop in PPI networks.
	  FormatType myformat(myoption);
	}
  }
  mylog.close();
  return 0;
}