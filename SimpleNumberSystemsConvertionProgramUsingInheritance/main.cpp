#include <iostream>
#include <cstring>
#include <cmath>
#include <vector>

using namespace std;

class NumberSystem
{
public:
  NumberSystem()=delete;
  NumberSystem(int _number, int _base):number(_number), base(base){}
  NumberSystem(int _base):base(_base)
  {
    digit_count = 0;
    std::cin>>number;
    int num = number;
    std::cout<<"Typed number="<<number<<" base="<<base<<std::endl;
    while(num != 0)
    {
      int rem = num%10;
      if(rem >= base) std::cout<<"***ERROR: Incorrect number input"<<std::endl;
      num /= 10;
      ++digit_count;
    }
  }
public:
  int ToDecimal()
  {
    int i=0, rem=0, decimalNum=0;
    int num = number;
    while(num != 0)
    {
      rem = num % 10;
      num /= 10;
      decimalNum += rem * pow(base, i);
      ++i;
    }
    cout << "Equivalent Decimal number: " << decimalNum << endl;
    return decimalNum;
  }
  int ToBase(int base)
  {
    int DecimalNumber = ToDecimal();
    int BaseNum = 0;
    int i=0, rem=0;

    std::vector<std::string> symvec;
    std::string NUMBER="";
    while(DecimalNumber != 0)
    {
      if(base != 16)
      {
        rem = DecimalNumber % base;
        DecimalNumber /= base;
        BaseNum += rem * pow(10, i);
        ++i;
        std::cout<<"rem="<<rem<<" base="<<base<<" DecimalNumber="<<DecimalNumber<<" BaseNum="<<BaseNum<<std::endl;
      }
      else
      {
        rem = DecimalNumber % base;
        std::string symbol;
        switch(rem)
        {
          case 15:
            symbol="F";
            break;
          case 14:
            symbol="E";
            break;
          case 13:
            symbol="D";
            break;
          case 12:
            symbol="C";
            break;
          case 11:
            symbol="B";
            break;
          case 10:
            symbol="A";
            break;
          default:
            symbol=std::to_string(rem);
            break;
        }
        DecimalNumber /= base;
        ++i;
        symvec.push_back(symbol);
        std::cout<<"rem="<<rem<<" base="<<base<<" DecimalNumber="<<DecimalNumber<<" BaseNum="<<BaseNum<<std::endl;
      }
    }
    if(base!=16)
      cout << "Equivalent "<<base<<"-base number: " << BaseNum << endl;
    else
    {
      for(int i=0; i<symvec.size(); ++i) NUMBER += symvec.at(symvec.size()-1-i);
      cout << "Equivalent "<<base<<"-base number: " << NUMBER << endl;
    }
    return BaseNum;
  }
  int ToBinary(){return ToBase(2);}
  int ToOctal(){return ToBase(8);}
  int ToHexaDecimal(){return ToBase(16);}
private:
  int number;
  int base;
  int digit_count;
};

class BinaryNS : public NumberSystem
{
public:
  BinaryNS() : NumberSystem(2){}
private:
};

class OctalNS : public NumberSystem
{
public:
  OctalNS() : NumberSystem(8){}
private:
};

class HexadecimalNS : public NumberSystem
{
public:
  HexadecimalNS() : NumberSystem(16){}
private:
};

int main(int argc, char ** argv)
{
  NumberSystem n = NumberSystem(2);
  std::cout<<"Decimal="<<n.ToDecimal()<<" Octal="<<n.ToBase(8)<<" Hexadecimal="<<n.ToBase(16)<<std::endl;
  std::cout<<"END"<<std::endl;
  return 0;
}

