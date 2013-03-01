/* semantic_go.h
Author: Jialu Hu
Data: 23.07.2012*/
#ifndef SEMANTIC_GO_H_
#define SEMANTIC_GO_H_
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <unordered_map>
#include <cassert>
#include "macro.h"
#include "stdlib.h"

template<typename NetworkPoolType,typename Option>
class GoList
{
private:
  typedef std::vector<std::string> GoTerms;
  typedef std::unordered_map<std::string,GoTerms*> GoMap;
  typedef typename NetworkPoolType::Graph Graph;
  typedef typename NetworkPoolType::InvOrigLabelNodeMap InvOrigLabelNodeMap;
  typedef typename NetworkPoolType::WeightEdgeMap WeightEdgeMap;
  typedef typename Graph::template NodeMap<int> OrigLabelNodeMap;
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  typedef struct _GOntology
  {
    GoTerms MF;
    GoTerms BP;
    GoTerms CC;
  }GOntology;
  typedef std::unordered_map<std::string,GOntology> GOntologyMap;
  typedef struct
  {
    std::string protein1;
    std::string protein2;
    std::string protein3;
  }AlignmentNode;
  typedef std::vector<std::string> AlignmentNodeVector;
  typedef struct _FunScore
  {
    float rfunsim;
    float funsim;
    float rfunsimall;
    float funsimall;
    float mf;
    float bp;
    float cc;
    _FunScore()
    :rfunsim(0)
    ,funsim(0)
    ,rfunsimall()
    ,funsimall()
    ,mf(0)
    ,bp(0)
    ,cc(0){}
    _FunScore& operator+= (const _FunScore& x){
      rfunsim+=x.rfunsim;
      funsim+=x.funsim;
      rfunsimall+=x.rfunsimall;
      funsimall+=x.funsimall;
      mf+=x.mf;
      bp+=x.bp;
      cc+=x.cc;
      return *this;
    }
    _FunScore& operator/= (const unsigned x)
    {
      rfunsim/=x;
      funsim/=x;
      rfunsimall/=x;
      funsimall/=x;
      mf/=x;
      bp/=x;
      cc/=x;
      return *this;
    }
  }FunScore;
  typedef std::array<AlignmentNodeVector*,NUM_ALIGNMENT_NODE> NewAlignmentList;
  typedef std::array<AlignmentNode,NUM_ALIGNMENT_NODE> AlignmentList;
  unsigned length;
  typedef std::unordered_map<std::string,FunScore> FunsimMap;
  int _numSpecies;
  int _numQualified;
  FunsimMap fun_map;
public:
  std::array<GoMap,NUM_SPECIES> go_list;
  GOntologyMap go_map;
  NewAlignmentList *newalignment;
  AlignmentList *alignment;
  std::string alignmentfile;
  std::string avefunsimfile;
  std::string resultfolder;
  std::vector<std::string> associationfiles;
  Graph _g;
  OrigLabelNodeMap node2label;
  InvOrigLabelNodeMap label2node;
  WeightEdgeMap weight;
  int numMF,numBP,numCC;/// Keep the number of annotated proteins in alignment graph
  GoList(Option& myoption)
  :length(0)
  ,_numSpecies(myoption.numspecies)
  ,_numQualified(0)
  ,fun_map()
  ,go_list()
  ,go_map()
  ,_g()
  ,node2label(_g)
  ,weight(_g)
  {
    newalignment=new NewAlignmentList();
    alignment = new AlignmentList();
    associationfiles = myoption.associationfiles;
    alignmentfile=myoption.alignmentfile;
    avefunsimfile=myoption.avefunsimfile;
    resultfolder=myoption.resultfolder;
    numMF=0;numBP=0;numCC=0;
  }
  ~GoList()
  {
    for(unsigned i=0;i<length;i++)
      delete (*newalignment)[i];
    delete newalignment;
    delete alignment;
  }
  bool readGeneOntology(const char*);/// Read gene ontology file for one species.
  bool convert_fsst();/// Output fsst format for alignment output.
  bool outputGOterm(std::ofstream&,std::string);
  //bool readList(const char*,int);
  bool readAlignment(const char*);
  bool goInitial();
  bool analyse();
  std::string getProtein(unsigned,short);
  GoTerms& getGoTerm(std::ofstream&,GOntology&,short);
  bool maskComment();
  bool deleteRedundancy();
  bool getNetworkAnnotation(NetworkPoolType&);
  bool getAlignmentCoverage(NetworkPoolType&);
  bool getMulFunSim(std::string&);
  bool readFsstResult(std::string&);
  bool isAnnotated(AlignmentNodeVector* pRecords);
  bool getAveFunSim(AlignmentNodeVector* pRecords,std::ofstream&,std::ofstream&,std::ofstream&);
  bool getAlignmentEdge(NetworkPoolType&);
  bool extractPValue();
};

template<typename NetworkPoolType,typename Option>
bool GoList<NetworkPoolType,Option>::readFsstResult(std::string& filename)
{
  if(!fun_map.empty())
    fun_map.clear();
  std::ifstream input(filename.c_str());
  std::string line;
  std::string proteinpair;
  std::getline(input,line);/// Ignore the header of input file.
  while(std::getline(input,line))
  {
    std::stringstream streamline(line);
    FunScore funscore;
    streamline >> proteinpair >> funscore.rfunsim >> funscore.funsim >> funscore.rfunsimall >> funscore.funsimall >> funscore.mf >> funscore.bp >> funscore.cc;
    fun_map[proteinpair] = funscore;  
  }
  return true;
}

template<typename NetworkPoolType,typename Option>
bool GoList<NetworkPoolType,Option>::getMulFunSim(std::string& filename)
{
  std::string filename1(resultfolder);
  std::string filename2(resultfolder);
  std::string filename3(resultfolder);
  filename1.append("aveFunSim.result");
  filename2.append("aveQualified.result");
  filename3.append("aveUnknown.result");
  std::ofstream output(filename1.c_str());
  std::ofstream outputQualified(filename2.c_str());
  std::ofstream outputUnknown(filename3.c_str());
  readFsstResult(filename);
  int numNAN=0;
  for(unsigned i=0;i<length;i++)
  {
    if(!isAnnotated((*newalignment)[i]))/// Count the number of "unknown records", any of proteins's annotation is unknown.
    {
      numNAN++;///
      int j=0;
      int recordSize = (*newalignment)[i]->size();
      outputUnknown <<numNAN<<"\t";
      for(; j<recordSize-1;++j)
      {
        if(go_map.find((*(*newalignment)[i])[j])==go_map.end())
          outputUnknown <<"*" << (*(*newalignment)[i])[j] <<"\t";
        else
          outputUnknown <<(*(*newalignment)[i])[j]<<"\t";
      }
      if(go_map.find((*(*newalignment)[i])[j])==go_map.end())
        outputUnknown << "*" <<(*(*newalignment)[i])[j] <<std::endl;
      else
        outputUnknown << (*(*newalignment)[i])[j] <<std::endl;
      //while(it!=pRecords->end())
      //{
        //if(go_map.find(*it)==go_map.end())
        //return false;
        //++it;
      //}
      std::string termfinder("./bin/gotermfinder.sh ");
      std::stringstream ss;
      ss << termfinder << numNAN;
      termfinder = ss.str();
      //system(termfinder.c_str());
      continue;
    }else
    {
	  std::string filename4("./result/genes.txt");
      std::ofstream outputGenes(filename4.c_str());
      getAveFunSim((*newalignment)[i],output,outputQualified,outputGenes);
      outputGenes.close();
    }
  }
  output << "#The number of unknown alignment records: "<< numNAN <<std::endl;
  output << "#The number of proteins annotated with MF,BP and CC in alginment graph respectively: "<< numMF <<"\t" <<numBP <<"\t" <<numCC <<std::endl;
  output.close();
  outputQualified.close();
  outputUnknown.close();
  return true;
}

template<typename NetworkPoolType,typename Option>
bool GoList<NetworkPoolType,Option>::getAveFunSim(AlignmentNodeVector* pRecords,std::ofstream& output,std::ofstream& outputQualified,std::ofstream& outputGenes)
{
  AlignmentNodeVector::iterator it,it2;
  std::string protein,protein1,protein2;
  FunScore ave;
  for(it=pRecords->begin();it!=pRecords->end()-1;++it)
  {
    for(it2=it+1;it2!=pRecords->end();++it2)
    {
      protein1=*it;protein2=*it2;
      if(g_verbosity>VERBOSE_ESSENTIAL)
        assert(protein1.compare(protein2)!=0);
      if(protein1.compare(protein2)>0)
      {
        protein.append(protein1);
        protein.append("_");
        protein.append(protein2);
      }else
      {
        protein.append(protein2);
        protein.append("_");
        protein.append(protein1);
      }
      if(g_verbosity>VERBOSE_ESSENTIAL)
        assert(fun_map.find(protein)!=fun_map.end());
      ave+=fun_map[protein];
      protein.clear();
    }
  }
  unsigned s=pRecords->size();
  ave/=s*(s-1)/2;
  for(it=pRecords->begin();it!=pRecords->end();++it)
  {
    output << "#" <<*it <<"\t";
    outputGenes << *it << std::endl;
    std::string protein=*it;
    if(go_map[protein].MF.size()!=0)
      numMF++;
    if(go_map[protein].BP.size()!=0)
      numBP++;
    if(go_map[protein].CC.size()!=0)
      numCC++;
  }
  output  << std::endl <<ave.rfunsim <<"\t"<<ave.funsim<<"\t"<< ave.rfunsimall<<"\t"
    <<ave.funsimall<<"\t"<< ave.mf<<"\t"<<ave.bp<<"\t"<<ave.cc<<std::endl;
    if(ave.mf >= 0.8 || ave.bp >= 0.6)
    {
      _numQualified++;
      outputQualified << _numQualified <<"\t";
      for(it=pRecords->begin();it!=pRecords->end();++it)
        outputQualified <<*it<<"\t";
      outputQualified << ave.mf <<"\t" << ave.bp <<std::endl;
      std::string termfinder("./bin/gotermfinder.sh ");
      std::stringstream ss;
      ss << termfinder << _numQualified;
      termfinder = ss.str();
      //system(termfinder.c_str());
    }
  //system("./bin/termfinderClient.sh ");
  return true;
}

template<typename NetworkPoolType,typename Option>
bool GoList<NetworkPoolType,Option>::isAnnotated(AlignmentNodeVector* pRecords)
{
  AlignmentNodeVector::iterator it=pRecords->begin();
  for(;it!=pRecords->end();++it)
  {
    if(go_map.find(*it)==go_map.end())
      return false;
  }
  return true;
}

template<typename NetworkPoolType,typename Option>
bool GoList<NetworkPoolType,Option>::readGeneOntology(const char* filename)
{
  /// Format of input file should be like this: Acc GO_term Inferrence Ontologies 
  std::ifstream input(filename);
  std::string line;
  while(std::getline(input,line))
  {
    if(line[0]=='!')continue;
    std::array<std::string,4> gorecord;
    std::stringstream streamline(line);
    for(int i=0;i<4;i++)
    {
      if(streamline.good())
      {
        streamline >> gorecord[i];
        //std::cout << gorecord[i] <<"\t";
      }
      else break;
    }

    /// If GO_term was inferred from IEA or ISS, ignore it.
    if(gorecord[2].compare("IEA")==0 ||gorecord[2].compare("ISS")==0)
      continue;
    GOntology& record=go_map[gorecord[0]];
    if(gorecord[3].compare("C")==0)
    {
       record.CC.push_back(gorecord[1]);
      /// Cellular Component annotation is not needed in our analysis
    }
    else if(gorecord[3].compare("F")==0)
    {
      record.MF.push_back(gorecord[1]);
    }
    else if(gorecord[3].compare("P")==0)
    {
      record.BP.push_back(gorecord[1]);
    }
    else{}
  }
  input.close();
  return true;
}

//bool GoList::readList(const char* filename, int i)
//{
  //std::ifstream input(filename);
  //std::string line;
  //while(std::getline(input,line))
  //{
    //std::stringstream streamline(line);
    //std::string protein;
    //GoTerms *goTerms=new GoTerms();
    //streamline >> protein;
    //while(streamline.good())
    //{
      //std::string goterm;
      //streamline >> goterm;
      //(*goTerms).insert(goterm);
    //}
    //go_list[i][protein] = goTerms;
  //}
  //return true;
//}

template<typename NetworkPoolType,typename Option>
bool GoList<NetworkPoolType,Option>::goInitial()
{
  for(unsigned i=0; i<associationfiles.size(); i++)
  {
    readGeneOntology(associationfiles[i].c_str());
  }
  ///readAlignment(alignmentfile.c_str());
  return true;
}

template<typename NetworkPoolType,typename Option>
bool GoList<NetworkPoolType,Option>::convert_fsst()
{
  std::string protein;
  std::string filename;
  filename.append(resultfolder);
  filename.append("fsst_revised.out");
  std::ofstream output2(filename.c_str());
  std::string line,commandline;
  filename.clear();
  filename.append(resultfolder);filename.append("alignment_proteins.fsst");
  commandline.append("./bin/fsst.sh ");commandline.append(resultfolder);
  for(unsigned i=0;i<length;i++)
  {
    std::ofstream output(filename.c_str());
    AlignmentNodeVector::iterator it,itEnd;
    it=(*newalignment)[i]->begin();
    itEnd=(*newalignment)[i]->end();
    while(it!=itEnd)
    {
      protein=*it;
      outputGOterm(output,protein);
      ++it;
    }
    output.close();
    system(commandline.c_str());
  }

  filename.clear();
  filename.append(resultfolder);filename.append("fsst.out");
  std::ifstream input;
  input.open(filename.c_str());
  while(std::getline(input,line))//acc	rank	acc	rfunsim (max)	funsim (max)	rfunsimAll (max)	funsimAll (max)	mfscore (max)	bpscore (max)	ccscore (max)	matching pfams	non-matching pfams	rank	acc	rfunsim (avg)	funsim (avg)	rfunsimAll (avg)	funsimAll (avg)	mfscore (avg)	bpscore (avg)	ccscore (avg)	matching pfams	non-matching pfams	
  {
    std::stringstream streamline(line);
    while(streamline.good())
    {
      std::string mystring;
      streamline >> mystring;
      if(mystring.compare("NaN")==0)
        output2 << "0" << "\t";
      else
        output2 << mystring << "\t";
    }
    output2 << std::endl;
  }
  input.close();
  output2.close();
  return true;
}

template<typename NetworkPoolType,typename Option>
bool GoList<NetworkPoolType,Option>::maskComment()
{
  std::ifstream input("./result/fsst.out");
  std::ofstream output("./result/fsst.temp.out");
  std::string line;
  while(std::getline(input,line))
  {
    std::string firststring;
    std::stringstream streamline(line);
    streamline >> firststring;
    if(firststring.compare("acc")!=0 || line[0]!='-')
    {
      output << line;
    }
  }
  input.close();
  output.close();
  return true;
}

template<typename NetworkPoolType,typename Option>
bool GoList<NetworkPoolType,Option>::deleteRedundancy()
{
  std::string filename1(resultfolder);filename1.append("fsst_revised.out");
  std::string filename2(avefunsimfile);//filename2.append("fsst.result");
  std::ifstream input(filename1.c_str());
  std::ofstream output(filename2.c_str());
  std::string line;
  
  while(std::getline(input,line))//acc	rank	acc	rfunsim (max)	funsim (max)	rfunsimAll (max)	funsimAll (max)	mfscore (max)	bpscore (max)	ccscore (max)	matching pfams	non-matching pfams	rank	acc	rfunsim (avg)	funsim (avg)	rfunsimAll (avg)	funsimAll (avg)	mfscore (avg)	bpscore (avg)	ccscore (avg)	matching pfams	non-matching pfams	
  {
    std::string protein1,protein2;
    int rank;
    FunScore funscore;
    std::stringstream streamline(line);
    streamline >> protein1 >> rank >> protein2 >> funscore.rfunsim >> funscore.funsim >> funscore.rfunsimall >> funscore.funsimall >> funscore.mf >> funscore.bp >> funscore.cc;
    if(protein1.compare(protein2)>0)
    {
      protein1.append("_");
      protein1.append(protein2);
      fun_map[protein1]=funscore;
    }
    else if(protein1.compare(protein2)<0)
    {
      protein2.append("_");
      protein2.append(protein1);
      fun_map[protein2]=funscore;
    }
    else
    {
    }
  }
  output <<"acc_acc\trfunsim\tfunsim\trfunsimall\tfunsimall\tmf\tbp\tcc"<<std::endl;
  for(typename FunsimMap::iterator it=fun_map.begin();it!=fun_map.end();++it)
  {
    output << it->first <<"\t"
    << it->second.rfunsim <<"\t"
    << it->second.funsim << "\t"
    << it->second.rfunsimall <<"\t"
    << it->second.funsimall << "\t"
    << it->second.mf <<"\t"
    << it->second.bp <<"\t"
    << it->second.cc <<std::endl;
  }
  input.close();
  output.close();
  return true;
}

template<typename NetworkPoolType,typename Option>
bool GoList<NetworkPoolType,Option>::outputGOterm(std::ofstream& output,std::string protein)
{
  if(go_map.find(protein)==go_map.end())
    return false;
  output << ">"<<protein << std::endl;
  GOntology& go=go_map[protein];
  for(short i=0;i<2;i++)
  {
    GoTerms& goterm=getGoTerm(output,go,i);
    for(GoTerms::iterator it=goterm.begin(); it!=goterm.end(); it++)
    {
      output << *it <<" ";
    }
    output << std::endl;
  }
  return true;
}

template<typename NetworkPoolType,typename Option>
typename GoList<NetworkPoolType,Option>::GoTerms& GoList<NetworkPoolType,Option>::getGoTerm(std::ofstream& output,GOntology& go, short i)
{
  GoTerms* goterm;
  switch(i)
  {
    case 0:goterm=&go.MF;output <<"MF:";break;
    case 1:goterm=&go.BP;output <<"BP:";break;
    case 2:goterm=&go.CC;output <<"CC:";break;
  }
  return *goterm;
}

template<typename NetworkPoolType,typename Option>
bool GoList<NetworkPoolType,Option>::readAlignment(const char* filename)
{
  std::ifstream input(filename);
  std::string line;
  unsigned len=0;
  std::array<unsigned, 6> numRecords;
  numRecords.fill(0);
  while(std::getline(input,line))
  {
    std::string str;
    std::stringstream streamline(line);
    if(line[0]=='#') continue;
    (*newalignment)[len]= new AlignmentNodeVector();
    Node node=_g.addNode();
    while(streamline.good())// reference of natelia
    {
      streamline >> str;
      if(!str.compare("NA"))continue;
      (*newalignment)[len]->push_back(str);
      label2node[str] = node;
      node2label[node] = len;
    }
    
    unsigned temp = (*newalignment)[len]->size();
    if(temp>=5)
      numRecords[5]++;
    else
      numRecords[temp]++;
    len++;
  }
  //std::cout << "The number of matching sets with 3 proteins: "<<numRecords[3] <<"\t"<<numRecords[4]<<"\t"<<numRecords[5]<< std::endl;
  length=len;
  return true;
}

template<typename NetworkPoolType,typename Option>
std::string GoList<NetworkPoolType,Option>::getProtein(unsigned len, short i)
{
  std::string pro;
  switch(i)
  {
    case 0: pro = (*alignment)[len].protein1; break;
    case 1: pro = (*alignment)[len].protein2; break;
    case 2: pro = (*alignment)[len].protein3; break;
    default : break;
  }
  return pro;
}

template<typename NetworkPoolType,typename Option>
bool GoList<NetworkPoolType,Option>::analyse()
{
  typedef std::unordered_map<std::string,unsigned> _GoMap;
  std::array<int, NUM_GO_SHARED> sharedNum;
  readAlignment(alignmentfile.c_str());
  for(int i=0;i<NUM_GO_SHARED;i++)
  {
    sharedNum[i]=0;
  }
  for(unsigned i=0;i<length;++i)
  {
    short numShared=0;
    _GoMap gomap;
    for(unsigned j=0; j<NUM_SPECIES; ++j)
    {
      std::string protein = getProtein(i,j);
      GoTerms* goterm = go_list[j][protein];
      for(GoTerms::iterator it=goterm->begin();it!=goterm->end();++it)
      {
        if(gomap.find(*it)!=gomap.end())
          ++gomap[*it];
        else
          gomap[*it]=1;
      }
    }
    for(_GoMap::iterator it=gomap.begin();it!=gomap.end();++it)
    {
      if(it->second==NUM_SPECIES) ++numShared;
      // std::cout << it->first <<"\t";
    }
    sharedNum[numShared]++;
  }
  return true;
}

template<typename NetworkPoolType,typename Option>
bool GoList<NetworkPoolType,Option>::getAlignmentEdge(NetworkPoolType& networks)
{
  for(NodeIt it(_g);it!=lemon::INVALID;++it)
  {
    for(NodeIt pt=it+1;pt!=lemon::INVALID;++pt)
    {
      /// conserved edges might be more than the number of species. So it makes no sense to count the number of conserved edges.
    }
  }
  return true;
}

template<typename NetworkPoolType,typename Option>
bool GoList<NetworkPoolType,Option>::getAlignmentCoverage(NetworkPoolType& networks)
{
  unsigned *numProtein = new unsigned[_numSpecies];
  unsigned *mySpecies = new unsigned[_numSpecies+1];
  unsigned *protein_k = new unsigned[_numSpecies+1];
  unsigned alignmentProtein=0;
  std::unordered_map<std::string,int> coveredProtein;
  for(unsigned i=0;i<length;i++)
  {
    unsigned mysize = (*newalignment)[i]->size();
    unsigned numSpecies=0;
    for(unsigned j=0;j<mysize;j++)
    {
      std::string myprotein = (*(*newalignment)[i])[j];
      if(coveredProtein.find(myprotein)!=coveredProtein.end())
      {
		coveredProtein[myprotein]++;
	  }
	  else
	  {
		  coveredProtein[myprotein]=1;
		  alignmentProtein++;
	  }
      //std::cout << i <<" " <<j <<" " << myprotein << std::endl;
      numProtein[networks.getHost(myprotein)]++;
    }
    for(short k=0;k<_numSpecies;k++)
    {
      if(numProtein[k]>0)
      {
        numSpecies++;
        numProtein[k]=0;
      }
    }
    mySpecies[numSpecies]++;
    protein_k[numSpecies]+=mysize;
  }
  std::cout <<"The number of proteins coverved by the alignment: "<<alignmentProtein<<std::endl;
  for(short k=1;k<=_numSpecies;k++)
  {
    std::cout <<"The number of clusters contains proteins from exactly"<<k<<" species:"<<mySpecies[k] <<"\\" << protein_k[k] <<std::endl;
  }
  delete [] numProtein;
  delete [] mySpecies;
  return true;
}

template<typename NetworkPoolType,typename Option>
bool GoList<NetworkPoolType,Option>::getNetworkAnnotation(NetworkPoolType& networks)
{
  Graph *graph;
  unsigned numMF,numBP,numCC,numGO;// number of proteins annotated by MF BP CC
  std::cout <<"MF\tBP\tCC\tANY\n";
  for(int i=0;i<_numSpecies;i++)
  {
    numMF=numBP=numCC=numGO=0;
    graph = networks.getGraph(i)->g;
    for(typename Graph::NodeIt it(*graph);it!=lemon::INVALID;++it)
    {
      std::string protein=(*networks.getGraph(i)->label)[it];
      if(go_map.find(protein)!=go_map.end())
      {
        GOntology& go=go_map[protein];
        if(go.MF.size()!=0)
          numMF++;
        if(go.BP.size()!=0)
          numBP++;
        if(go.CC.size()!=0)
          numCC++;
        numGO++;
      }
    }
    std::cout << numMF <<"\t"<< numBP<<"\t"<<numCC<<"\t"<<numGO <<std::endl;
  }
  return true;
}
template<typename NetworkPoolType,typename Option>
bool GoList<NetworkPoolType,Option>::extractPValue()
{
  std::ifstream input("./result/alignment_raw_pvalue.txt");
  std::ofstream output("./result/alignment_pvalue.txt");
  std::string line;
  double pvalue=0;
  while(std::getline(input,line))
  {
    std::string str;
    std::stringstream streamline(line);
    while(streamline.good())// reference of natelia
    {
      streamline >> str;
      if(!str.compare("CORRECTED"))
      {
        streamline >> str;
        streamline >> pvalue;
        output << pvalue << std::endl;
      }
    }
  }
  input.close();
  output.close();
  return true;
}
#endif /// SEMANTIC_GO_H_
