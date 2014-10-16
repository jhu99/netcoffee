#ifndef __GETPOST_H__
#define __GETPOST_H__

#include <string.h>
#include <iostream>
#include <unordered_map>
#include <vector>

bool getPost(std::unordered_map<std::string,std::string> &formData)
{
	char *buffer=NULL;
	char *strlength = getenv("CONTENT_LENGTH");
	if (strlength == NULL) return false;
	int content_length = atoi(strlength);
	if(content_length==0)return false;
	buffer = new char[(content_length+1)*sizeof(char)];
	if(!fread(buffer, sizeof(char), content_length, stdin))return false;
	buffer[content_length]='\0';
	//std::cout << buffer<<"<br>"<<std::endl;
	char *key,*value;
	int i=0;
	key=buffer;
	std::string ikey,ivalue;
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
			ikey=key;
			i++;
			key=buffer+i;
			formData[ikey]=ivalue;
		}
		else
		{
			i++;
		}
	}	
	ikey=key;
	ivalue=value;
	formData[ikey]=ivalue;
	return true;
}
#endif /*__GETPOST_H__*/
