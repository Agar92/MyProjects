#include <conio.h>
#include <cstdlib>
#include <iostream>
#include <random>
#include <unistd.h>
#include "windows.h"
#include <ctime>
enum Direction {STOP=0, LEFT, RIGHT, UP, DOWN};
Direction dir;
bool gameOver;
int x, y, FruitX, FruitY, width=20, height=20;
using namespace std;
int level=3;//1 - quick, 10 - slow.
int score=0;
int tailX[100], tailY[100];
int nTail=0;
long int Time=0;

void Setup()
{
   gameOver=false;
   dir=STOP;
   x=width/2;
   y=height/2;
   FruitX=rand() % width-1;
   FruitY=rand() % height;
   if(FruitX<1) FruitX+=1;
   if(FruitY<1) FruitY+=1;
   if(FruitX>width-2) FruitX-=1;
   if(FruitY>height-2) FruitY-=1;
}

void Draw()
{
    system("cls");
    for(int i=0; i<width; ++i) std::cout<<"#";
    std::cout<<std::endl;
    //
    //std::cout<<"x="<<x<<" y="<<y<<std::endl;
    for(int i=0; i<height; ++i)
    {
        for(int j=0; j<width; ++j)
        {
            if(x==j && y==i) std::cout<<"0";
            else if(FruitX==j && FruitY==i) std::cout<<"F";
            else if(j==0 || j==width-1) std::cout<<"#";
            else
            {
                bool print=false;
                for(int k=0; k<nTail; ++k)
                {
                    if(tailX[k]==j && tailY[k]==i)
                    {
                        cout<<"o";
                        print=true;
                    }
                }
                if(!print) cout<<" ";
            }
        }
        cout<<endl;
    }
    for(int i=0; i<width; ++i) std::cout<<"#";
    std::cout<<std::endl;
    std::cout<<"Score="<<score<<" Time="<<Time<<" s"<<std::endl;
    //std::cout<<"x="<<x<<std::endl;
    //cout<<"%x="<<x<<" %y="<<y<<endl;
    Sleep(level*100);//delay=100ms*level.
}

void Logic()
{
    int prevX=tailX[0];
    int prevY=tailY[0];
    tailX[0]=x;
    tailY[0]=y;
    //std::cout<<"1@@@ x="<<x<<" y="<<y<<" t0="<<tailX[0]<<" t1="<<tailX[1]<<endl;
    switch(dir)
    {
        case LEFT:
                --x;
               break;
        case RIGHT:
                ++x;
               break;
        case UP:
                --y;
               break;
        case DOWN:
                ++y;
               break;
        case 'x':
            gameOver = true;
               break;
    }
    if(x<1)     x=width-2;
    else if(x>width-1) x=1;
    if(y<1)     y=height-2;
    else if(y>height-1) y=1;
    //
    cout<<"x="<<x<<" y="<<y<<" Fx="<<FruitX<<" Fy="<<FruitY<<endl;
    if(x==FruitX && y==FruitY)
    {
        //cout<<"x==FruitX && y==FruitY"<<endl;
        //if(FruitX=x) FruitX=1, FruitY=1;
        FruitX=tailX[0];
        FruitY=tailY[0];
        bool Flag=true;
        while(Flag)
        {
            cout<<"Flag"<<endl;
            FruitX=rand()%width;
            FruitY=rand()%height;
            bool Flag1=false;
            for(int k=0; k<nTail; ++k)
            {
                if(FruitX==tailX[k] && FruitY==tailY[k]) Flag1=true;
            }
            if(FruitX<1 || FruitX>width-1 || FruitY<1 || FruitY>height-1) Flag1=true;
            Flag=Flag1;
        }
        score+=10;
        ++nTail;
    }
    //
    int prev2X, prev2Y;
    for(int i=1; i<nTail; ++i)
    {
        prev2X=tailX[i];
        prev2Y=tailY[i];
        tailX[i]=prevX;
        tailY[i]=prevY;
        prevX=prev2X;
        prevY=prev2Y;
    }
    //std::cout<<"2@@@ x="<<x<<" y="<<y<<" t0="<<tailX[0]<<" t1="<<tailX[1]<<" Fx="<<FruitX<<" Fy="<<FruitY<<endl;
    for(int k=0; k<nTail; ++k) if(x==tailX[k] && tailY[k]==y) gameOver=true;
    //std::cout<<"x="<<x<<" nTail="<<nTail<<" xt0="<<tailX[0]<<" "<<" xt1="<<tailX[1]<<endl;
    //Sleep(1000);
}

void Input()
{
    bool flag=_kbhit();
    cout<<"flag="<<flag<<endl;
    if(flag)
    {
        //std::cout<<"WE ARE HERE"<<endl;
        //cout<<"_________________________________________________________________"<<endl;
        switch(getch())
        {
            //getchar();
            case 'a':
                dir = LEFT;
                //cout<<"&&&LEFT";
                break;
            case 'd':
                dir = RIGHT;
                //cout<<"&&&RIGHT";
                break;
            case 'w':
                dir = UP;
                //cout<<"&&&UP";
                break;
            case 's':
                dir = DOWN;
                //cout<<"&&&DOWN";
                break;
            case 'x':
                gameOver = true;
                break;
        }
    }
    //cout<<"@x="<<x<<" @y="<<y<<endl;
    //Sleep(3000);

}
int main(int argc, char *argv[])
{
    Setup();
    int count=0;
    clock_t begin = clock();
    while(!gameOver)//count<1000)//!gameOver)
    {
        //cout<<"---------------------"<<endl;
        //Draw();
        Input();
        Logic();
        clock_t end = clock();
        Time = double(end - begin) / CLOCKS_PER_SEC;//time in seconds
        //getchar();
        Draw();
        //cout<<"dir="<<dir<<" LEFT="<<LEFT<<" UP="<<UP<<endl;
        ++count;
        //getchar();
        //sleep(0.1);
    }
    if(gameOver) cout<<"Game Over"<<endl;

    return 0;
}
