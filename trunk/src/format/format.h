/* format.h
Author: Jialu Hu
Data: 10.12.2012*/
#ifndef FORMAT_H_
#define FORMAT_H_
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <unordered_map>

template<typename NetworksType,typename MyOption>
class Format
{
public:
  Format(MyOption&);
  ~Format(){};
  std::vector<std::string> blastfile;
  std::unordered_multimap<std::string,std::string>  KOmap;
  std::unordered_map<std::string,unsigned> KOgroup;
  //std::array<
  bool removeBiEdges(NetworksType&);
  bool formatAlignment(std::string,std::string);
  bool removeRodundantInteraction();
  bool retrieveKOnumber(std::string);
  bool extractHomologyProteins(std::string,NetworksType&);
  bool extractGoAssociation(std::string);
  bool extractGraemlinAlignment(std::string,std::string);
  bool readTrainingData(std::string);
  bool retrieveKOgroups(std::string);
};

template<typename NetworksType,typename MyOption>
Format<NetworksType,MyOption>::Format(MyOption& myoption)
{
  blastfile = myoption.blastfiles;
}

template<typename NetworksType,typename MyOption>
bool Format<NetworksType,MyOption>::removeBiEdges(NetworksType& networks)
{
  std::vector<std::string>::iterator it;
  for(it=blastfile.begin();it!=blastfile.end();++it)
  {
    std::string filename = *it;
    std::string outname=filename;
    outname.append(".sw");
    std::ifstream input(filename);
    std::ofstream output(outname);
    if(!input.is_open() && !output.is_open())
      std::cout << "Error in opening the files" << std::endl;
    std::string line="";
    while(std::getline(input,line))
    {
      std::string protein1,protein2;
      std::stringstream streamline(line);
      streamline >> protein1 >> protein2;
      if(networks.existNode(protein1) && networks.existNode(protein2))
      {
        output << line <<std::endl;
      }
    }
  }
  return true;
}

template<typename NetworksType,typename MyOption>
bool Format<NetworksType,MyOption>::formatAlignment(std::string alignmentfile,std::string outfile)
{
  std::ifstream input(alignmentfile.c_str());
  std::ofstream output(outfile.c_str());
  std::string line;
  while(std::getline(input,line))
  {
    unsigned i=0;
    std::stringstream streamline(line);
    std::vector<std::string> proteins;
    while(streamline.good())
    {
      std::string protein;
      streamline >> protein;
      if(protein.compare("NA")!=0)
      proteins.push_back(protein);
    }
    if(proteins.size()>=2)
    {
      for(i=0;i<proteins.size()-1;i++)
      {
        output << proteins[i] << "\t";
      }
      output << proteins[i] << std::endl;
    }
  }
  return true;	
}

template<typename NetworksType,typename MyOption>
bool Format<NetworksType,MyOption>::removeRodundantInteraction()
{
	return true;
}

template<typename NetworksType,typename MyOption>
bool Format<NetworksType,MyOption>::readTrainingData(std::string formatfile)
{
	std::ifstream input(formatfile.c_str());
	std::string line;
	unsigned index=1;
	while(std::getline(input,line))
	{
		std::stringstream streamline(line);
		while(streamline.good())
		{
			std::string term;
			streamline >> term;
			KOgroup[term] = index;
		}
		index++;
	}
	return true;
}

template<typename NetworksType,typename MyOption>
bool Format<NetworksType,MyOption>::retrieveKOgroups(std::string formatfile)
{
	std::ifstream input(formatfile.c_str());
	std::string line;
	std::string trainingfile("./benchmark/graemlin/graemlin-2.0_test_files/test_cases/6way/training.txt");
	unsigned qualifiedNum=0;
	unsigned linenum=0;
	readTrainingData(trainingfile.c_str());
	while(std::getline(input,line))
	{
		std::vector<std::string> matchset;
		std::unordered_map<unsigned,unsigned> indexMap;
		std::stringstream streamline(line);
		unsigned maxKOnum=0;
		unsigned matchsetSize=0;
		while(streamline.good())
		{
			std::string term;
			streamline >> term;
			matchset.push_back(term);
			if(KOgroup.find(term)==KOgroup.end())continue;
			unsigned myindex=KOgroup[term];//discard all nodes in the alignment without a KO group
			matchsetSize++;
			if(indexMap.find(myindex)!=indexMap.end())
			{
				indexMap[myindex]++;
			}else
			{
				indexMap[myindex]=1;
			}
			if(indexMap[myindex]>maxKOnum) maxKOnum=indexMap[myindex];
		}
		float conservedRate=maxKOnum/(1.0*matchsetSize);
		if(conservedRate>=0.6) qualifiedNum+=matchsetSize;
		linenum++;
	}
	std::cout <<"#The number of qualified match-sets, match-sets, conrectness in "<<formatfile <<" is:"<<std::endl;
	std::cout << qualifiedNum <<"\t" << linenum << "\t" <<qualifiedNum/(1.0*linenum)<< std::endl;
	return true;
}
template<typename NetworksType,typename MyOption>
bool Format<NetworksType,MyOption>::retrieveKOnumber(std::string formatfile)
{
	std::ifstream input(formatfile.c_str());
	std::string line;
	std::string module;
	std::string KOnumber;
	while(std::getline(input,line))
	{
	  std::stringstream streamline(line);
	  std::string koindex;
      streamline >> koindex;
      if(koindex.compare("D")==0)
      {
		  streamline >> module;
	  }else if(koindex.compare("E")==0)
	  {
		  streamline >> KOnumber;
		  std::cout << KOnumber << std::endl;
		  KOmap.insert(std::make_pair(module,KOnumber));
	  }
	  else continue;
	}
	return true;
}

template<typename NetworksType,typename MyOption>
bool Format<NetworksType,MyOption>::extractHomologyProteins(std::string filename,NetworksType& networks)
{
	std::ifstream input(filename.c_str());
	std::string inputNet[]={"Celeg20130131","Dmela20130131","Ecoli20130131","Hpylo20130131","Hsapi20130131","Mmusc20130131","Rnorv20130131","Scere20130131"};
	std::string line;
	std::getline(input,line);// Header line
	std::unordered_map<std::string,int> checklist;
	while(std::getline(input,line))
	{
		std::stringstream streamline(line);
		std::string protein1,protein2;
		double evalue;
		streamline >> protein1 >> protein2 >> evalue;
		unsigned i=networks.getHost(protein1);
		unsigned j=networks.getHost(protein2);
		if(i==100 || j==100) continue;
		if(i>j)
		{
			unsigned temp=i;
			std::string protein=protein1;
			i=j;
			j=temp;
			protein1=protein2;
			protein2=protein;
		}
		if(checklist.find(protein1)==checklist.end())
		{
		  checklist[protein1]=1;
		  std::cout << protein1 <<std::endl;
	    }
		if(checklist.find(protein2)==checklist.end())
		{
		  checklist[protein2]=1;
		  std::cout << protein2 <<std::endl;
	    }
		std::string outfilename("./dataset/dip/");
    outfilename.append(inputNet[i]);
		outfilename.append("-");
		outfilename.append(inputNet[j]);
		outfilename.append(".evals");
		std::ofstream output(outfilename.c_str(),std::ios_base::out|std::ios_base::app);
		output << protein1 <<"\t" << protein2 << "\t" << evalue << std::endl;
		output.close();		
	}
	return true;
}

template<typename NetworksType,typename MyOption>
bool Format<NetworksType,MyOption>::extractGoAssociation(std::string formatfile)
{
	std::ifstream input(formatfile.c_str());
	std::ofstream output("./dataset/goa/gene_association.goa_target");
	std::unordered_map<std::string,int> checklist;
	std::string line;
	while(std::getline(input,line))
	{
		std::stringstream streamline(line);
		std::string protein;
		streamline >> protein;
		if(checklist.find(protein)==checklist.end())
		{
		  checklist[protein]=1;
	    }				
	}
	input.close();
	std::ifstream input2("./dataset/goa/gene_association.goa_uniprot");
	if(!input2.is_open())
	{
		std::cerr << " Gene association file can't open." << std::endl;
	}
	while(std::getline(input2,line))
	{
		if(line[0]=='!')continue;
		std::vector<std::string> goterm;
		std::stringstream streamline(line);
		while(streamline.good())
		{
			std::string term;
			streamline >> term;
			goterm.push_back(term);
		}
		if(checklist.find(goterm[1])!=checklist.end())
			output << line << std::endl;
		 //output << goterm[1] << "\t" << goterm[3] << "\t" << goterm[5] << "\t" << goterm[7] << std::endl;
	}
	input2.close();
	output.close();
	return true;
}

template<typename NetworksType,typename MyOption>
bool Format<NetworksType,MyOption>::extractGraemlinAlignment(std::string formatfile,std::string outfile)
{
	std::ifstream input(formatfile.c_str());
	std::ifstream input2("./maptable.txt");
	std::ofstream output(outfile.c_str());
	std::string line;
	std::unordered_map<std::string,std::string> idmap;
	
	while(std::getline(input2,line))
	{
		std::stringstream streamline(line);
		std::string id1,id2;
		streamline >> id1 >> id2;
		idmap[id1]=id2;
	}

    /// The next line goes for header line. If there is no header line, set it as comment.
	std::getline(input, line);
	while(std::getline(input,line))
	{
		std::vector<std::string> matchset;
		std::stringstream streamline(line);
		while(streamline.good())
		{
			std::string term;
			streamline >> term;
			matchset.push_back(term);
		}
		if(matchset.size()<2)continue;
		unsigned i=0;
		for(i=0;i<matchset.size()-1;i++)
		{
		   std::string uniprot_id;
		   if(idmap.find(matchset[i])!=idmap.end())
		   {
			  uniprot_id=idmap[matchset[i]];
			  output << uniprot_id <<"\t";
		   }else
		   {
			   output << matchset[i] <<"\t";
		   }
	   }
		if(idmap.find(matchset[i])!=idmap.end())
		{
		  output << idmap[matchset[i]] <<std::endl;
	   }else
	   {
		   output << matchset[i] <<std::endl;
	   }
	}
	return true;
}
#endif
