#include <iostream>
#include <cmath>
#include <cstring>
#include <fstream>
#include <locale.h>

#include <vector>
#include <map>

using namespace std;

int main (int argc, char ** argv)
{
  //setlocale( LC_ALL, "Russian"); // для вывода на экран русского текста
  //setlocale(LC_CTYPE, "Russian_Russia.1251"); 
  ifstream fin;
  fin.open("data.dat", ios::in);
  char c;
  map<char, int> freqmap;
  std::vector<char> vecchar;
  std::vector<int>  vecfreq(0);
  while (!fin.eof())
  {
    fin.get(c);
    cout<<c<<" ";
    
    int I=-1;
    for(int i=0; i<vecchar.size(); ++i)
    {
      if(vecchar.at(i)==c) I=i;
    }

    if(I>-1) ++vecfreq.at(I);
    else      vecchar.push_back(c); vecfreq.push_back(1);
   
    if(freqmap.count(c)>0) freqmap[c]++;
    else freqmap[c]=1;
    //freqmap[c]=1;
   
  }
  
  /*
  c='a';
  freqmap[c]=1;
  c='b';
  freqmap[c]=2;
  */

  map<char,int>::iterator it = freqmap.begin();
  int i=0;
  for(; it != freqmap.end(); it++)
  {
    cout<<i<<"  "<<"f="<<(int)it->first<<" f="<<it->first<<"  s="<<it->second<<endl;
    ++i;  
  }

  cout<<"Size of the map is="<<freqmap.size()<<endl;

  std::cout<<"SIZE="<<vecchar.size()<<std::endl;
  for(int i=0; i<vecchar.size(); ++i)
  {
    std::cout<<vecchar.at(i)<<"   "<<vecfreq.at(i)<<std::endl;
  }

  
  cout<<"After sort:"<<endl;
  for(int i=0; i<vecchar.size(); ++i)
  {
    for(int j=0; j<vecchar.size()-i; ++j)
    {
      if(vecfreq.at(j)<vecfreq.at(j+1))
      {
        int temp=vecfreq.at(j+1);
        vecfreq.at(j+1)=vecfreq.at(j);
        vecfreq.at(j)=temp;
         
        char ctemp=vecchar.at(j+1);
        vecchar.at(j+1)=vecchar.at(j);
        vecchar.at(j)=ctemp;
      }
    }
  }

  std::cout<<"SIZE="<<vecchar.size()<<std::endl;
  for(int i=0; i<vecchar.size(); ++i)
  {
    std::cout<<vecchar.at(i)<<"   "<<vecfreq.at(i)<<std::endl;
  }
 
  cout<<endl;
  return 0;
}
