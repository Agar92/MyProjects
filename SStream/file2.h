#ifndef FILE2_H
#define FILE2_H

#include <cmath>
#include <fstream>
#include <string>

int power(int x, unsigned int y)//y - power
{
    if( y == 0) return 1;
    if (y%2 == 0) return power(x, y/2)*power(x, y/2);
    return x*power(x, y/2)*power(x, y/2);
}

//function to calculate the number of digits in the given number x
int order(int x)
{
    int n = 0;
    while(x)
    {
        n++;
        x = x/10;
    }
    return n;
}

//function to check whether the given number is Armstrong number or not
bool isArmstrong(int x)
{
    int n = order(x);
    int temp = x, sum = 0;
    while(temp)
    {
        int r = temp%10;
        sum += power(r, n);
        temp = temp/10;
    }
    //satisfies or not the Armstrong condition
    return (sum == x);
}

bool IsPrime(int x)
{
    for(int i=2; i<(int)std::sqrt(x)+1; ++i)
        if(x%i==0) return false;
    return true;
}

//greatest common divisor
int GCD(int a,int b){return b?GCD(b,a%b):a;}
//least common multiple
int LCM(int a,int b){return a*b/GCD(a,b);}

void CopyFileToAFile(std::string source_filename, std::string dest_filename)
{
    std::ifstream src(source_filename, std::ios::binary);
    std::ofstream dst(dest_filename,   std::ios::binary);
    dst << src.rdbuf();
}
#endif // FILE2_H
