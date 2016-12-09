#include <limits> // for std::numeric_limits
#include <string> // for std::string
#include <iostream> // I/O
#include <fstream>// file I/O


ifstream& SkipComment(ifstream &strm)
{
   char c;
   const int c_MAX_LINESIZE=std::numeric_limits<int>::max();

   top:
   
   if(!strm) return strm;
   
   while((isspace(c=strm.get()))||(c==',')||(c==';')) ;

   if(c=='#'||c=='%'||((c=='/') &&( strm.peek()=='/')))
   {
	//skip the rest of the line
        
	strm.ignore(c_MAX_LINESIZE, '\n');
        goto top;		                  
   }
   else if(c=='/' && strm.peek()=='*')
   {
	   //skip everything in the comment block
	   c=strm.get();          //skip the first '*'
	   char last='\0';
	   while(!(strm.eof())&&strm.good())
	   {
		   c=strm.get();
		   if(c=='/'&&last=='*')break;
		   else last=c;
	   }
	   return strm;

   }
   else if(c!=EOF)
   {
	   strm.putback(c);
   }

     if(c=='#'||c=='%'||((c=='/') &&( strm.peek()=='/')))
        {
		//skip the rest of the line
		strm.ignore(c_MAX_LINESIZE, '\n');

	}
		

   return strm;
}
