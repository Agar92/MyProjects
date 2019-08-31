#ifndef FILE1_H
#define FILE1_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <iomanip>

using namespace std;
void print(){cout<<"END OF PRINT\n";}

template <typename T, typename... Types>
void print(T var1, Types&&... var2)
{
    cout<<var1<<endl;
    print(var2...) ;
}

int countWords(string str)
{
    stringstream s(str);
    string word;
    int count = 0;
    while (s >> word)
        count++;
    return count;
}

void printFrequency(string s)
{
    map<string, int> FW;
    stringstream ss(s);
    string Word;
    while (ss >> Word) FW[Word]++;
    map<string, int>::iterator m;
    for (m = FW.begin(); m != FW.end(); m++)
        cout<<m->first<<" -> "<< m->second<<"\n";
}

//Function to remove spaces from a string
string removeSpaces(string str)
{
    stringstream ss;
    string temp;
    //store the whole string into string stream
    ss << str;
    //make the string empty
    str = "";
    //run loop till end of stream
    while (!ss.eof())
    {
        ss >> temp;
        str = str + temp;
    }
    return str;
}

template <typename T>
bool StrToNum(const std::string& s, T &t)
{
    std::istringstream is(s);
    return (is>>t)?true:false;
}

template <typename T>
string NumToStr(T &&t)
{
    std::ostringstream os;
    os.precision(10);
    os<<t;
    return os.str();
}

void PrintUintDecimalToHexaDecimal(unsigned int i)
{
    stringstream ss;
    ss<<hex<<i;
    string res=ss.str();
    cout<<"0x"<<res<<endl;
}

void PrintHexadecimalToUint(string hexStr="0x3ae")
{
    unsigned int x;
    stringstream ss;
    ss<<std::hex<<hexStr;
    ss>>x;
    cout<<x<<endl;
}

void PrintUintDecimalToOctal(unsigned int i)
{
    stringstream ss;
    ss<<oct<<i;
    string res=ss.str();
    cout<<"0"<<res<<endl;
}

void PrintOctalToUint(string octStr="01656")
{
    unsigned int x;
    stringstream ss;
    ss<<std::oct<<octStr;
    ss>>x;
    cout<<x<<endl;
}

//function to find all words in a string that have length > K
void findWords(string str, int K)
{
    string word;
    stringstream ss(str);
    int count = 0;
    while (ss >> word) {
        if (word.size() > K) {
            cout << word << " ";
            count++;
        }
    }
}

//works only from console (not from IDE) for English
void ReadALineToAFile(string WhatToWriteToAFile)
{
    std::streambuf *psbuf, *backup;
    std::ofstream filestr;
    filestr.open ("test.txt", ios::app);
    backup=std::cout.rdbuf();// back up cout's streambuf
    psbuf=filestr.rdbuf();   // get file's streambuf
    std::cout.rdbuf(psbuf);  // assign streambuf to cout
    std::cout<<"WhatToWriteToAFile is written to the file"<<endl;
    std::cout<<WhatToWriteToAFile<<endl;
    std::cout.rdbuf(backup); // restore cout's original streambuf
    filestr.close();
}


#endif // FILE1_H
