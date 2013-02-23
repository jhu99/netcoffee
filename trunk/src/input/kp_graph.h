/* kp_graph.h
Author: Jialu Hu
Date: 28.06.2012*/

#ifndef KP_GRAPH_H_
#define KP_GRAPH_H_

#include <lemon/core.h>
#include <lemon/bits/graph_extender.h>
#include <iostream>
#include <string>
#include <array>
#include <unordered_map>
#include <utility>
#include "macro.h"
#include "verbose.h"
#include <assert.h>

template<typename NetworkPool>
class KpGraph
{
public:
  typedef NetworkPool NetwookPool_Type;
  typedef std::unordered_multimap<std::string,std::string> EdgeMap;
  typedef std::unordered_map<std::string,double> WeightMap;
  typedef std::unordered_map<std::string,unsigned> StWeightMap;
  typedef typename NetworkPool::Graph Graph;
  typedef typename NetworkPool::GraphData GraphData;
  typedef typename NetworkPool::OrigLabelNodeMap OrigLabelNodeMap;
  typedef typename NetworkPool::InvOrigLabelNodeMap InvOrigLabelNodeMap;
  typedef typename Graph::template EdgeMap<float> EdgeWeight;
  typedef struct _BpGraph
  {
    EdgeMap redBlue;
    EdgeMap blueRed;
    WeightMap seWeight;
    StWeightMap stWeight;
  }BpGraph;
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  
  std::vector<std::string> fileList;
  std::vector<BpGraph*> graphs;
  int numSpecies;
  unsigned maxStrWeight;
  
  KpGraph(std::vector<std::string>&,int);//{fileList=files;}
  int  getBpIndex(int,int);
  bool readHomoList(std::string&,BpGraph*,int,int);
  bool constructGraph();
  bool reweighting(NetworkPool&,int,int);
  bool reweightingAll(NetworkPool&);
  bool isEdge(std::string protein,GraphData*);
  bool createBpGraph(Graph&,OrigLabelNodeMap&,InvOrigLabelNodeMap&,EdgeWeight&,NetworkPool&);
};

template<typename NetworkPool>
KpGraph<NetworkPool>::KpGraph(std::vector<std::string>& files,int species)
:numSpecies(species)
,maxStrWeight(0)
{
  fileList=files;
}
template<typename NetworkPool>
bool KpGraph<NetworkPool>::isEdge(std::string protein,GraphData* network)
{
  if( network->invIdNodeMap->find(protein)==network->invIdNodeMap->end())
  {
    return false;
  }
  return true;
}
template<typename NetworkPool>
bool KpGraph<NetworkPool>::createBpGraph(Graph& gr,
                                         OrigLabelNodeMap& node2label,
                                         InvOrigLabelNodeMap& label2node,
                                         EdgeWeight& edgemap,
                                         NetworkPool& networkpool)
{
  BpGraph* bp12=graphs[0];
  EdgeMap::iterator it;
  typename NetworkPool::GraphData* network_1;
  typename NetworkPool::GraphData* network_2;
  network_1 = networkpool.getGraph(0);// The first network
  network_2 = networkpool.getGraph(1);// The second network
  std::string protein1,protein2;
  for(it=bp12->redBlue.begin();it!=bp12->redBlue.end();++it)
  {
    protein1=it->first;
    protein2=it->second;
    if(!isEdge(protein1,network_1)||!isEdge(protein2,network_2))continue;
    Node node1,node2;
    Edge e;
    if(label2node.find(protein1)==label2node.end())
    {
      node1 = gr.addNode();
      node2label.set(node1,protein1);
      label2node[protein1] = node1;
    } 
    else
    {
      node1 = label2node[protein1];
    }
    if(label2node.find(protein2)==label2node.end())
    {
      node2 = gr.addNode();
      node2label.set(node2,protein2);
      label2node[protein2] = node2;
    }
    else
    {
      node2 = label2node[protein2];
    }
    e = gr.addEdge(node1,node2);
    edgemap.set(e,0.0);
  }
  return true;
}

template<typename NetworkPool>
int KpGraph<NetworkPool>::getBpIndex(int i,int j)
{
  if (i>=j || j>=numSpecies)
    std::cerr << "Triplets doesn't exist!" << std::endl;
    return (i*(2*numSpecies-i-1))/2+j;
}

template<typename NetworkPool>
bool KpGraph<NetworkPool>::constructGraph()
{
  graphs.clear();
  int i=0;
  for(int ni=0;ni<numSpecies;ni++)
  {
    for(int nj=ni;nj<numSpecies;nj++,i++)
    {
		graphs.push_back(new BpGraph());
		readHomoList(fileList[i],graphs[i],ni,nj);
	}
  }
  return true;
}

//template<typename NetworkPool>
//float KpGraph<NetworkPool>::scoreCalculate(double e1,double e2,double e3)
//{
  //float score;
  //score=getScore(e1)+getScore(e2)+getScore(e3);
  //return score;
//}

//template<typename NetworkPool>
//float KpGraph<NetworkPool>::getScore(double evalue)
//{
  //int index=0;
  //if(evalue>1.0e-180)
  //{
    //index = static_cast<int>(ceil(log10(evalue))/NUM_PRO_INTERVAL)+NUM_OFFSET_BARS;
    ////std::cout << evalue << index <<std::endl;
  //}
  //if(index>=NUM_PRO_BARS)
  //{
    //std::cerr <<"E-Value "<<evalue<<" is out of range!"<<std::endl;
    //return 0.0;
  //}
  //return (*resultDistr)[index];
//}

template<typename NetworkPool>
bool KpGraph<NetworkPool>::readHomoList(std::string& filename,BpGraph* graph,int ni,int nj)
{
  std::ifstream input(filename.c_str());
  std::string line;
  std::string protein1="";
  std::string protein2="";
  double evalue;
  if(!input.is_open())
  {
    std::cerr << filename << "cannot be opened" <<std::endl;
    return false;
  }
  while(std::getline(input,line))
  {
    std::stringstream lineStream(line);
    lineStream >> protein1 >> protein2 >> evalue;
    std::string keystring;
    if(ni==nj)
    {
		if(protein1.compare(protein2)>0)
		{
			std::string temp = protein1;
			protein1 = protein2;
			protein2 = temp;
		}
	}
	keystring.append(protein1);keystring.append(protein2);
    if(graph->seWeight.find(keystring)!=graph->seWeight.end())
    {
      if(graph->seWeight[keystring]>evalue)
        graph->seWeight[keystring]=evalue;
    }
    else
    {
      graph->seWeight[keystring]=evalue;
      graph->stWeight[keystring]=0;
      graph->redBlue.insert(std::make_pair(protein1,protein2));
      graph->blueRed.insert(std::make_pair(protein2,protein1));
    }
  }
  return true;
}

template<typename NetworkPool>
bool KpGraph<NetworkPool>::reweightingAll(NetworkPool& networkpool)
{
  constructGraph();
  for(int ni=0;ni<numSpecies-1;ni++)
    for(int nj=ni+1;nj<numSpecies;nj++)
      reweighting(networkpool,ni,nj);
  return true;
}

template<typename NetworkPool>
bool KpGraph<NetworkPool>::reweighting(NetworkPool& networkpool,int ni,int nj)
{
  int myindex=0;
  myindex=getBpIndex(ni,nj);
  EdgeMap::iterator it;
  typename NetworkPool::GraphData* network_1;
  typename NetworkPool::GraphData* network_2;
  typename NetworkPool::GraphData* network_k;
  network_1 = networkpool.getGraph(ni);// The first network
  network_2 = networkpool.getGraph(nj);// The second network
  BpGraph* bp12=graphs[myindex];
  Node nodeA1,nodeA2,nodeA3,nodeB1,nodeB2,nodeB3;
  std::string kst1,kst2,kst3,kst4,kst5,kst6;
  int nc=0;
  for(it=bp12->redBlue.begin();it!=bp12->redBlue.end();++it,++nc)
  {
    //std::cerr <<"The number of steps: "<<nc<<std::endl;
    std::string protein1,protein2,protein3,protein4,protein5,protein6;//Protein in left and right circles.
    protein1=it->first;
    protein2=it->second;
    if(!isEdge(protein1,network_1)||!isEdge(protein2,network_2))continue;
    kst1.append(protein1);kst1.append(protein2);
    nodeA1=(*network_1->invIdNodeMap)[protein1];
    nodeA2=(*network_2->invIdNodeMap)[protein2];
    for(int x=nj+1;x<numSpecies;++x)
    {
      network_k = networkpool.getGraph(x);// The third network
      BpGraph* bp1k=graphs[getBpIndex(ni,x)];
      BpGraph* bp2k=graphs[getBpIndex(nj,x)];
      auto range=bp1k->redBlue.equal_range(protein1);
      EdgeMap::iterator pt;
      for(pt=range.first;pt!=range.second;++pt)
      {
        protein3=pt->second;
        if(!isEdge(protein3,network_k))continue;
        nodeA3 =(*network_k->invIdNodeMap)[protein3];
        kst2.append(protein1);
        kst2.append(protein3);
        kst3.append(protein2);
        kst3.append(protein3);
        if(bp2k->stWeight.find(kst3)!=bp2k->stWeight.end())
        {
          /// Find a left circle. Then we need to check whether there exist a right circle.
          for(IncEdgeIt ie((*network_1->g),nodeA1);ie!=lemon::INVALID;++ie)
          {
            nodeB1=network_1->g->runningNode(ie);
            protein4=(*network_1->label)[nodeB1];
            if(!isEdge(protein4,network_1))continue;
            for(IncEdgeIt pe((*network_2->g),nodeA2);pe!=lemon::INVALID;++pe)
            {
              nodeB2=network_2->g->runningNode(pe);
              protein5=(*network_2->label)[nodeB2];
              if(!isEdge(protein5,network_2))continue;
              kst4.append(protein4);
              kst4.append(protein5);
              for(IncEdgeIt me((*network_k->g),nodeA3);me!=lemon::INVALID;++me)
              {
                nodeB3=network_k->g->runningNode(me);
                protein6=(*network_k->label)[nodeB3];
                if(!isEdge(protein6,network_k))continue;
                kst5.append(protein4);
                kst5.append(protein6);
                kst6.append(protein5);
                kst6.append(protein6);
                if(bp1k->stWeight.find(kst5)!=bp1k->stWeight.end()
                  && bp2k->stWeight.find(kst6)!=bp2k->stWeight.end())
                  {
                    /// reweight on match edges.
                    bp1k->stWeight[kst2]++;
                    bp1k->stWeight[kst5]++;
                    bp2k->stWeight[kst3]++;
                    bp2k->stWeight[kst6]++;
                    bp12->stWeight[kst1]++;
                    if(maxStrWeight < bp12->stWeight[kst1])
                      maxStrWeight = bp12->stWeight[kst1];
                    if(bp12->stWeight.find(kst4)!=bp12->stWeight.end())
                    {
                      bp12->stWeight[kst4]++;
                      if(maxStrWeight < bp12->stWeight[kst4])
                        maxStrWeight = bp12->stWeight[kst4];
                    }
                  }
                  kst5.clear();kst6.clear();
              }
              kst4.clear();
            }
          }
        }
        kst2.clear();kst3.clear();
      }
    }
    kst1.clear();
  }
  return true;
}
#endif /// KP_GRAPH_H_
