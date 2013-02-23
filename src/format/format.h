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

template<typename NetworksType,typename MyOption>
class Format
{
public:
  Format(MyOption&);
  ~Format(){};
  std::vector<std::string> blastfile;
  bool removeBiEdges(NetworksType&);
  bool formatAlignment(std::string,std::string);
  bool removeRodundantInteraction();
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
#endif
