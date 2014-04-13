/* format.h
Author: Jialu Hu
Data: 10.04.2014*/
#ifndef OUTPUT_HTML_H_
#define OUTPUT_HTML_H_

#include <exception>
#include <iostream>
#include <string>
#include <vector>
#include "getpost.h"
#include "function.h"
//#include "mysql_connection.h"
//#include <cppconn/driver.h>
//#include <cppconn/exception.h>
//#include <cppconn/resultset.h>
//#include <cppconn/statement.h>

using namespace std;

template<typename Option>
class Output_html
{
public:
	std::string jobid;
	Output_html()
	{
	}
	void set_header();
	void set_footer();
	void get_data(Option&);
};

template<typename Option>
void Output_html<Option>::set_header()
{
	std::cout <<"Content-Type: text/html;charset=us-ascii\n\n<html>\n<head><title>A Web Server for Network Alignment Problem</title></head>\n<body bgcolor=\"white\">\
	<h1 style=\"text-align:center;color:blue\">NetCoffee & LocalAli</h1>\
	<p style=\"text-align:center;font-size:20px\"><i>A Web Server for Multiple (Global & Local) Network Alignment</i></p><hr>\n\
	<table style=\"width:800px\" bgcolor=#aoffa0>\
		<tr>\
			<td><a href=\"../index.php\">Introduction</a></td>\
			<td><a href=\"../application.php\">Application</a></td>\
			<td><a href=\"../job.php\">Submitted Jobs</a></td>\
			<td><a href=\"../contact.php\">Contact Us</a></td>\
		</tr>\
	</table>";
}

template<typename Option>
void Output_html<Option>::set_footer()
{
	std::cout <<"</body></html>\n";
}

template<typename Option>
void Output_html<Option>::get_data(Option& myoption)
{
	std::string jobid,numspecies;
	if(!getPost(jobid,numspecies))return;
	std::cout <<"Job Id:\t"<<jobid<<std::endl;
	std::cout <<"<br>numspecies:\t"<<numspecies<<std::endl;
	int num=std::stoi(numspecies);
	for(int i=0;i<num;i++)
	{
		std::string ppifilename="./uploadfiles/";
		ppifilename.append(jobid);ppifilename.append("/ppi");
		ppifilename.append(convert_num2str(i));
		ppifilename.append(".tab");
		std::cout <<"<br>"<<ppifilename<<std::endl;
		myoption.networkfiles.push_back(ppifilename);
		for(int j=i;j<num;j++)
		{
			std::string ssfilename="./uploadfiles/";
			ssfilename.append(jobid);ssfilename.append("/ppi");
			ssfilename.append(convert_num2str(i));
			ssfilename.append("-ppi");ssfilename.append(convert_num2str(j));
			ssfilename.append(".cf");
			std::cout <<"<br>"<<ssfilename<<std::endl;
		}
	}
}

#endif //OUTPUT_HTML_H_
