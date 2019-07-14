#include <iostream>
#include "treenode.h"
#include "binstree.h"

//double G(double t, double f(double x)){return t * f(t);}
//double SQR(double x){return x*x;}
int main(int argc, char *argv[])
{
    auto visit=[](int& item){std::cout<<1<<" ";};

    TreeNode<double> * root=new TreeNode<double>(50.0);
    BinSTree<double> * bstr=new BinSTree<double>();
    bstr->SetRoot(root);
    bstr->Print(root);
    TreeNode<double> **p=&root;
    bstr->Insert(root, 30);
    std::cout<<std::endl;
    bstr->Print(root);
    bstr->Insert(root, 2);
    bstr->Insert(root, 300);
    bstr->Insert(root, 70);
    std::cout<<std::endl;
    bstr->Print(root);

    //std::cout<<G(6,SQR)<<std::endl;
    return 0;
}
