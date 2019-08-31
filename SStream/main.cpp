#include <iostream>
#include "File1.h"

using namespace std;

int main()
{
    print(1, 234, 67.5, "MAY");
    string s="dfsgdgdfg sdfsdf dsfsdfs";
    cout<<countWords(s)<<endl;
    printFrequency(s);
    cout<<removeSpaces(s)<<endl;
    float s1;
    StrToNum("9.5", s1);
    cout<<s1<<endl;
    string s2=NumToStr(3456.789);
    cout<<s2<<endl;
    PrintUintDecimalToHexaDecimal(942);
    PrintHexadecimalToUint();
    PrintUintDecimalToOctal(942);
    PrintOctalToUint();
    findWords(" 323 sdfsd sdfsdfsdfsdf 23234 67567567657567fdgdfgd", 10);
    cout<<"\nPlease, type what to write to a file:"<<endl;
    string str;
    //works only from console (not from IDE) for English
    cin>>str;
    cout<<str<<endl;
    //1st write
    ReadALineToAFile(str);
    str.clear();
    //2nd write
    cin>>str;
    ReadALineToAFile(str);
    return 0;
}
