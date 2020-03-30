#include <cmath>
#include <iostream>
#include <chrono>
#include <cassert>

#include "LinkedList.h"
#include "BinaryTree.h"

using namespace std;


int main(int argc, char ** argv)
{

    BinaryTree<double> tree;
    tree.Insert(1);
    tree.Insert(10);
    tree.Insert(100);
    tree.Insert(-1);
    tree.Insert(5000);
    tree.Insert(-5000);
    tree.PrintPreOrder();
    cout<<"\ntree.countNodes()="<<tree.countNodes()<<endl;
    cout<<"tree.search(5000)->data="<<tree.search(5000)->data<<endl;
    tree.Delete(-1);
    tree.PrintPreOrder();
    exit(0);
    List<double> list;
    cout<<list.AddNode(1)<<endl;
    list.AddNode(99);
    list.AddNode(101);
    list.AddNode(150);
    cout<<list.GetSize()<<endl;
    list.AddNode(-1, 4);
    list.AddNode(1000);
    list.AddNode(2000);
    list.AddNode(3000);
    list.AddNode(100000);
    list.AddNode(-100000);
    list.PrintAll();
    cout<<"\n"<<list.GetSize()<<endl;
    list.RemoveNodeP(0);
    list.PrintAll();
    cout<<"\n"<<list.GetSize()<<endl;
    list.Display(list.GetHead());
    list.SortAsc();
    list.SortDesc();
    return 0;
}
