#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>

#include <vector>
#include <list>
#include <map>
#include <queue>

//чтобы можно было использовать sleep() на Linux:
//#include "unistd.h"
//чтобы вывести размер файла в байтах:
#include <sys/stat.h>
 
using namespace std;

//узел бинарного дерева Хаффмана:
class Node
{
public:
  int  freq;
  char c;
  Node *left, *right;
  Node(char data, int _freq):c(data),freq(_freq),left(NULL),right(NULL){}
private:
  Node(){}
};

//функция печати дерева Хоффмана в виде графа
//k нужны для отступов.
void print(Node * root, unsigned int k=0)
{
  if(root!=NULL)
  {
    print(root->left, k+3);
    for(unsigned int i=0; i<k; ++i) cout<<" ";
    if(root->c) cout<<root->freq<<" ("<<root->c<<")"<<endl;
    else cout<<root->freq<<endl;
    print(root->right, k+3);
  }
}

//возвращает размер файла в байтах:
int getFileSize(const char * fileName)
{
  struct stat file_stat;
  stat(fileName, &file_stat);
  return file_stat.st_size;
}

//вспомогательные контейнеры:
//code - вспомогательный вектор для записи бинарного кода char символов
//в таблицу (символ char -> его двоичный код) table.
std::vector<bool> code;
std::map<char, vector<bool> > table;
//функция построения таблицы (символ char -> его двоичный код):
void BuildTableCharCode(Node * root)
{
  if(root->left!=NULL)
  {
    code.push_back(0);
    BuildTableCharCode(root->left);
  }
  if(root->right!=NULL)
  {
    code.push_back(1);
    BuildTableCharCode(root->right);
  }
  if(root->c) table[root->c]=code;
  code.pop_back();
}

//эта вспомогательная функция печатает таблицу символ char->его бинарный код,
//если требуетсяя её вывести.
void printTableCharBinaryCode()
{
  std::cout<<"Check char symbols binary codes:"<<std::endl;
  for(std::map<char, vector<bool> >::iterator iter=table.begin(); iter!=table.end(); ++iter)
  {
    std::cout<<iter->first<<" ";
    for(int j=0; j<iter->second.size(); ++j) std:cout<<iter->second.at(j)<<" ";
    std::cout<<std::endl;
  }
}

//эта структура необходима для правильной сортировки priority_queue:
struct Compare
{
  bool operator ()(Node * l, Node * r) const
  {return l->freq > r->freq;}
};

//освобождает все динамически выделенные до этого для построения бинарного
//дерева Хоффмана указатели типа Node*.
void DeleteNodesOfHuffmanTreeUsingPostOrderTraversal(Node*& root)
{
  if(root!=NULL)
  {
    DeleteNodesOfHuffmanTreeUsingPostOrderTraversal(root->left);
    DeleteNodesOfHuffmanTreeUsingPostOrderTraversal(root->right);
    delete root;
    root = NULL;
  }
}

//функия подсчитывает сколько раз каждый символ char в тексте в нём встречался:
map<char, int> CalculateCharFrequencies(ifstream & f)
{
  map<char, int> m;
  while(!f.eof())
  {
    //нужно использовать get(), потому что f>>c не считывает символы пробелов:
    char c=f.get();
    m[c]++;
  }
  //печать получившейся таблицы (символ char -> сколько раз он встретился в тексте):
  for(std::map<char, int>::iterator itr=m.begin(); itr!=m.end(); ++itr) cout<<itr->first<<" "<<itr->second<<endl;
  //печать размера получившейся таблицы (символ char -> сколько раз он встретился в тексте)
  //и сколько байт в памяти она занимает:
  cout<<"TABLE_SIZE="<<m.size()<<" Number of bytes required for the table of frequencies="<<m.size()*(1+4)+4<<" bytes"<<endl;
  
  return m;
}

//функция построения бинарного дерева Хаффмана.
//возвращает root узел дерева Хаффмана:
Node * BuildHuffmanBinaryTree(map<char, int> m)
{
  priority_queue<Node*, vector<Node*>, Compare> t;
  //заполнение priority_queue:
  for(map<char, int>::iterator itr=m.begin(); itr!=m.end(); ++itr) t.push(new Node(itr->first, itr->second));
  //построение дерева Хаффмана:
  while(t.size()!=1)
  {
    Node * sonL=t.top();
    t.pop();
    Node * sonR=t.top();
    t.pop();
    Node * top=new Node('$', sonL->freq+sonR->freq);
    top->left=sonL;
    top->right=sonR;
    t.push(top);
  }
  return t.top();
}

//записывает таблицу (символ char -> сколько раз он встретился в тексте) и
//последовательность бинарных кодов символов в бинарный файл:
void SaveToBinaryFile(map<char, int> m, ifstream & f)
{
  //файл для бинарного вывода по умолчанию:
  ofstream out("output.bin");
  //размер таблицы (символ char->сколько раз он втретился в тексте):
  const int MAP_SIZE=m.size();
  char chararr[MAP_SIZE]{0};
  int  freqarr[MAP_SIZE]{0};
  int index=0;
  for(std::map<char, int>::iterator iter=m.begin(); iter!=m.end(); ++iter)
  {
    chararr[index]=iter->first;
    freqarr[index]=iter->second;
    ++index;
  }
  
  /*
  cout<<"1 chararr:"<<endl;
  cout<<"MAP_SIZE="<<MAP_SIZE<<endl;
  for(int i=0; i<MAP_SIZE; ++i) cout<<chararr[i]<<" ";
  cout<<endl;

  cout<<"1 freqarr:"<<endl;
  for(int i=0; i<MAP_SIZE; ++i) cout<<freqarr[i]<<" ";
  cout<<endl;
  */

  //записываем таблицу (символ char->сколько раз он втретился в тексте) в файл для бинарного вывода:
  out.write((const char *)&MAP_SIZE, sizeof(int));
  out.write((const char *)chararr, MAP_SIZE*sizeof(char));
  out.write((const char *)freqarr, MAP_SIZE*sizeof(int));
  int count=0;
  char buf=0;
  f.clear();
  //сдвигаем курсор в файле ввода в начало файла:
  f.seekg(0);
  while(!f.eof())
  {
    char c=f.get();
    vector<bool> x=table[c];
    //побитово записываем бинарный код:
    for(int n=0; n<x.size(); ++n)
    {
      buf = buf | x[n]<<(7-count);
      count++;
      if(count==8)
      {
        count=0;
        out.put(buf);
        buf=0;
      }
    }
  }
  f.close();
  out.close();
}

//читает таблицу (символ char -> сколько раз он встретился в тексте) и
//последовательность бинарных кодов символов из бинарного файла:
map<char, int> ReadFromBinaryFile(ifstream & F)
{ 
  int NEW_MAP_SIZE;
  F.read((char *)&NEW_MAP_SIZE, sizeof(int));
  const int SIZE=NEW_MAP_SIZE;
  char newchararr[SIZE]{0};
  int  newfreqarr[SIZE]{0};
  F.read((char *)newchararr, SIZE*sizeof(char));
  F.read((char *)newfreqarr, SIZE*sizeof(int));
  map<char, int> m;
  for(int i=0; i<SIZE; ++i) m[newchararr[i]]=newfreqarr[i];

  /*
  cout<<"2 chararr:"<<endl;
  cout<<"NEW_MAP_SIZE="<<NEW_MAP_SIZE<<endl;
  for(int i=0; i<NEW_MAP_SIZE; ++i) cout<<newchararr[i]<<" ";
  cout<<endl;
  cout<<"2 freqarr:"<<endl;
  for(int i=0; i<NEW_MAP_SIZE; ++i) cout<<newfreqarr[i]<<" ";
  cout<<endl;
  */

  return m;
}

//функция, выводящая последовательность закодированных в бинарном файле символов:
void decodeHelpFunction(Node * root, ifstream & F)
{
  Node * p = root;
  int count=0;
  char byte=F.get();
  if(!F.good())
  {
    std::cout<<"***ERROR: F.good()="<<F.good()<<std::endl;
    exit(0);
  }
  while(!F.eof())
  {
    bool b = byte & (1 << (7-count) );
    if(b) p=p->right;
    else  p=p->left;
    if(p->left==NULL && p->right==NULL) 
    {
      cout<<p->c;
      p=root;
    }
    count++;
    if(count==8)
    {
      count=0;
      byte=F.get();
    }
  }
  F.close();
}

//Кодирование-сжатие:
void Encode(const char * filename)
{
  ifstream f(filename);
  map<char, int> m=CalculateCharFrequencies(f);
  Node * root=BuildHuffmanBinaryTree(m);
  //печать получившегося дерева Хаффмана, если необходимо:
  //print(root);
  BuildTableCharCode(root);
  //печать получившейся таблицы (символ char -> его двоичный код), если необходимо:
  printTableCharBinaryCode();
  SaveToBinaryFile(m, f);
  DeleteNodesOfHuffmanTreeUsingPostOrderTraversal(root);
}

//Раскодирование/извлечение:
void Decode()
{
  ifstream F("output.bin");
  map<char, int> m2=ReadFromBinaryFile(F);
  //очищаем code и table, т. к. их необходимо использовать заново:
  code.clear();
  table.clear();
  //построение дерева Хаффмана по таблице (символ char->сколько раз он втретился в тексте),
  //прочитанной в output.bin:
  Node * root2=BuildHuffmanBinaryTree(m2);
  BuildTableCharCode(root2);
  decodeHelpFunction(root2, F);  
  DeleteNodesOfHuffmanTreeUsingPostOrderTraversal(root2);
}


//Т. к. в условии сказано, что входной текстовый файл весит до 10 МБ,
//то число типа int, принимающее значения от -2^32 до 2^32-1, для сохранения
//частот появления символов в тексте подойдёт. 
int main ()
{ 
  //текстовый файл, который нужно сжать:
  const char * filename=/*"BigText.txt";*/"SmallText.txt";
  //по умолчанию, закодированная последовательность
  //сохраняется в бинарный файл output.bin.
  Encode(filename);
  //по умолчанию, закодированная последовательность
  //читается из бинарного файла output.bin.
  Decode(); 
  //вычиcляем коэффициент сжатия:
  int size1=getFileSize(filename);
  cout<<"Size of "<<filename<<"="<<size1<<" bytes"<<endl;
  //по умолчанию, закодированная последовательность
  //сохраняется в бинарный файл output.bin.
  int size2=getFileSize("output.bin");
  cout<<"Size of output.bin="<<size2<<" bytes"<<endl;
  cout<<"Coefficient of compression="<<static_cast<double>(size1)/size2<<endl;
  return 0;
}
