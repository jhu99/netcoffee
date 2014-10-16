/* format.h
Author: Jialu Hu
Data: 10.04.2014*/
#ifndef OUTPUT_HTML_H_
#define OUTPUT_HTML_H_

#include <exception>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include "getpost.h"
#include "function.h"
#include "mysql_connection.h"
#include <cppconn/driver.h>
#include <cppconn/exception.h>
#include <cppconn/resultset.h>
#include <cppconn/statement.h>


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
	void update_table();
};

template<typename Option>
void Output_html<Option>::set_header()
{
	std::cout <<"Content-Type: text/html;charset=us-ascii\n\n";
	/*<html><head><title>A Web Server for Network Alignment Problem</title></head>\n<body bgcolor=\"white\">\
	<h1 style=\"text-align:center;color:blue\">NetCoffee & LocalAli</h1>\
	<p style=\"text-align:center;font-size:20px\"><i>A Web Server for Multiple (Global & Local) Network Alignment</i></p><hr>\n\
	<table style=\"width:800px\" bgcolor=#aoffa0>\
		<tr>\
			<td><a href=\"index.php\">Introduction</a></td>\
			<td><a href=\"application.php\">Application</a></td>\
			<td><a href=\"job.php\">Submitted Jobs</a></td>\
			<td><a href=\"contact.php\">Contact Us</a></td>\
		</tr>\
	</table>";*/
}

template<typename Option>
void Output_html<Option>::set_footer()
{
	std::cout <<"</body></html>\n";
}

template<typename Option>
void Output_html<Option>::get_data(Option& myoption)
{
	std::unordered_map<std::string,std::string> formData;
	std::string numspecies,alpha,eta,model,algorithm;
	if(!getPost(formData))return;
	jobid=formData["jobid"];
	numspecies=formData["numspecies"];
	model=formData["model"];
	alpha=formData["alpha"];
	algorithm=formData["algorithm"];
	eta=formData["eta"];
	std::cout <<"Job id:\t"<<jobid<<"<br>"<<std::endl;
	std::cout <<"numspecies:\t"<<numspecies<<"<br>"<<std::endl;
	std::cout <<"model:\t"<<model<<"<br>"<<std::endl;
	int num=std::stoi(numspecies);
	myoption.numspecies=num;
	myoption.alpha=std::stof(alpha);
	myoption.eta=std::stof(eta);
	std::cout <<"alpha:\t"<<myoption.alpha<<"<br>"<<std::endl;
	std::cout <<"eta:\t"<<myoption.eta<<"<br>"<<std::endl;
	std::cout <<"algorithm:\t"<<algorithm<<"<br>"<<std::endl;
	
	for(int i=0;i<num;i++)
	{
		std::string ppifilename="../../data-rw/uploadfiles/";
		ppifilename.append(jobid);ppifilename.append("/ppi");
		ppifilename.append(convert_num2str(i));
		ppifilename.append(".tab");
		myoption.networkfiles.push_back(ppifilename);
		for(int j=i;j<num;j++)
		{
			std::string ssfilename="../../data-rw/uploadfiles/";
			ssfilename.append(jobid);ssfilename.append("/ppi");
			ssfilename.append(convert_num2str(i));
			ssfilename.append("-ppi");ssfilename.append(convert_num2str(j));
			ssfilename.append(".cf");
			myoption.blastfiles.push_back(ssfilename);
		}
	}
	if(model.compare("evalue")==0)
	{
		myoption.scorefile="../../data-rw/model/score_composit.model";
	}
	else
	{
		myoption.bscore=true;
		myoption.scorefile="../../data-rw/model/score_composit.bmodel";
	}
	myoption.resultfolder.append(jobid);
	myoption.resultfolder.append("/");
	myoption.alignmentfile.append(myoption.resultfolder);
	myoption.alignmentfile.append("alignment_netcoffee.data");
	myoption.logfile.append(myoption.resultfolder);
	myoption.logfile.append("log.txt");
}

template<typename Option>
void Output_html<Option>::update_table()
{
	sql::Driver *driver;
	sql::Connection *con;
	sql::Statement *stmt;
	
	driver = get_driver_instance();
	//con = driver->connect("mysql.imp.fu-berlin.de","jhu","ADcVBRwLfy4XTfFh");
	//con->setSchema("jhu");
	//stmt = con->createStatement();
	//std::string sqlstatement("UPDATE job_description SET status='completed' WHERE id='");
	//sqlstatement.append(jobid);
	//sqlstatement.append("'");
	//std::cout << sqlstatement <<"<br>\n";
	//stmt->execute(sqlstatement.c_str());
	delete stmt;
	delete con;
}

#endif //OUTPUT_HTML_H_
