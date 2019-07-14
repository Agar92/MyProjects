#ifndef BINSTREE_H
#define BINSTREE_H

#include <iostream>
#include <stdlib.h>
#include "Treenode.h"

template <class T>
class BinSTree
{
  protected:  // требуется для наследования
    TreeNode<T> *root;// указатели на корень и на текущий узел
    TreeNode<T> *current;
    int size;// число элементов дерева
    TreeNode<T> *GetTreeNode(const T& item,
                      TreeNode<T> *lptr=NULL, TreeNode<T> *rptr=NULL);
    void FreeTreeNode(TreeNode<T> *p);
    void DeleteTree(TreeNode<T> *t);
  public:
    // конструкторы и деструктор
    BinSTree(){}
    BinSTree(const BinSTree<T>& tree){}
    ~BinSTree(void){}
     // оператор присваивания
     BinSTree<T>& operator= (const BinSTree<T>& rhs);
     TreeNode<T> * getNodeByValue(TreeNode<T> *root, T value);
     void Insert(TreeNode<T> *head, T value);
     // Возвращает указатель на корень.
     TreeNode<T> *GetRoot(void) const {return root;}
     void SetRoot(TreeNode<T> * _root) {root=_root;}
     TreeNode<T>* getMinNode(TreeNode<T> *root) const;
     TreeNode<T>* getMaxNode(TreeNode<T> *root) const;
     void Print(TreeNode<T> * root);
};

template <class T>
TreeNode<T> * BinSTree<T>::GetTreeNode(const T & item, TreeNode<T> *lptr,
                           TreeNode<T> *rptr)
{
  TreeNode<T> *p;
  p = new TreeNode<T> (item, lptr, rptr);
  if (p == NULL)//если памяти недостаточно, завершить программу сообщением об ошибке
  {
    std::cerr << "Ошибка при выделении памяти!\n";
    exit(1);
  }
  return p;
}

template <class T>
void BinSTree<T>::Print(TreeNode<T> * root)
{
   TreeNode<T> * t=root;
   // рекурсивное прохождение завершается на пустом поддереве
   if (t != NULL)
   {
      Print(t->Left());
      std::cout<<t->data<<" ";
      Print(t->Right());
   }
}

// оператор присваивания
template <class T>
BinSTree<T>& BinSTree<T>::operator = (const BinSTree<T>& rhs)
{
  // нельзя копировать дерево в само себя
  if (this == &rhs)
    return *this;
  root = CopyTree(rhs.root);
  // присвоить текущему указателю значение корня и задать
  // размер дерева
  current = root;
  size = rhs.size;
  // возвратить ссылку на текущий объект
  return *this;
}

template <class T>
TreeNode<T> * BinSTree<T>::getNodeByValue(TreeNode<T> *root, T value)
{
  while (root)
  {
    if (root->data > value)
    {
      root = root->left;
      continue;
    }
    else if (root->data < value)
    {
      root = root->right;
      continue;
    }
    else return root;
  }
  return NULL;
}

template <class T>
void BinSTree<T>::Insert(TreeNode<T> *head, T value)
{
    TreeNode<T> *tmp = NULL;
    TreeNode<T> *ins = NULL;
    if (head == NULL)
    {
      head = GetTreeNode(value);
      return;
    }
    tmp = head;
    while (tmp)
    {
      if (value>tmp->data)
      {
        if (tmp->right)
        {
          tmp = tmp->right;
          continue;
        }
        else
        {
          tmp->right = GetTreeNode(value);
          return;
        }
      }
      else if (value<tmp->data)
      {
        if (tmp->left)
        {
          tmp = tmp->left;
          continue;
        }
        else
        {
          tmp->left = GetTreeNode(value);
          return;
        }
      }
      else exit(2);
    }
}

template <class T>
TreeNode<T>* BinSTree<T>::getMinNode(TreeNode<T> *root) const
{
  while (root->left) root = root->left;
  return root;
}

template <class T>
TreeNode<T>* BinSTree<T>::getMaxNode(TreeNode<T> *root) const
{
  while (root->right) root = root->right;
  return root;
}

#endif // BINSTREE_H
