#ifndef __GETPOST_H__
#define __GETPOST_H__

#include <string.h>
#include <iostream>
#include <vector>

bool getPost(std::string &jobid,std::string &numspecies)
{
	char *buffer=NULL;
	char *strlength = getenv("CONTENT_LENGTH");
	if (strlength == NULL) return false;
	int content_length = atoi(strlength);
	if(content_length==0)return false;
	buffer = new char[(content_length+1)*sizeof(char)];
	if(!fread(buffer, sizeof(char), content_length, stdin))return false;
	buffer[content_length]='\0';
	char *value,*ivalue;
	int i=0;
	while(buffer[i]!='\0')
	{
		if(buffer[i]=='=')
		{
			buffer[i]='\0';
			value=buffer+i+1;
			i++;
		}else if(buffer[i]=='&')
		{
			buffer[i]='\0';
			ivalue=value;
			i++;
			jobid=ivalue;
		}
		else
		{
			i++;
		}
	}	
	ivalue=value;
	numspecies=ivalue;
	return true;
}
#endif /*__GETPOST_H__*/
