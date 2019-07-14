#ifndef TREENODE_H
#define TREENODE_H

#include <iostream>

// Следующее объявление понадобится в дальнейшем
template <class T>
class BinSTree;

//объявление объекта узла бинарного дерева
template <class T>
class TreeNode
{
  private:
    //указатели левого и правого дочерних узлов
    TreeNode<T> *left;
    TreeNode<T> *right;
  public:
    T data;
    TreeNode(const T& item, TreeNode<T> *lptr=NULL,
      TreeNode<T> *rptr=NULL):data(item), left(lptr), right(rptr){}
    TreeNode<T>* Left() const{return left;}
    TreeNode<T>* Right() const{return right;}
    //сделать класс BinSTree дружественным, поскольку необходим доступ к полям left и right
    friend class BinSTree<T>;
};

template <class T>
void FreeTreeNode(TreeNode<T> *p){delete p;}

//симметричное рекурсивное прохождение узлов дерева
template <class T>
void Inorder (TreeNode<T> *t, void visit(T& item))
{
   // рекурсивное прохождение завершает-ся на пустом поддереве
   if (t != NULL)
   {
      Inorder(t->Left(), visit);
      visit(t->data);
      Inorder(t->Right(), visit);
   }
}

// при посещении узла проверяется, является ли он листовым
template <class T>
void CountLeaf(TreeNode<T> *t, int& count)
{
  // Использовать обратный метод прохождения
  if (t != NULL)
  {
    CountLeaf(t->Left(), count);
    CountLeaf(t->Right(), count);
    if (t->Left() == NULL && t->Right() == NULL)
      count++;
  }
}

template <class T>
int Depth (TreeNode<T> *t)
{
  int depthLeft, depthRight, depthval;

  if (t == NULL)
    depthval = -1;
  else
  {
    depthLeft = Depth(t->Left());
    depthRight = Depth(t->Right());
    depthval = 1 + (depthLeft > depthRight ? depthLeft : depthRight);
  }
  return depthval;
}

// создать дубликат дерева t и возвратить корень нового дерева
template <class T>
TreeNode<T> *CopyTree(TreeNode<T> *t)
{
   // переменная newnode указывает на новый узел, создаваемый
   // посредством вызова GetTreeNode и присоединяемый в дальнейшем
   // к новому дереву. указатели newlptr и newrptr адресуют сыновей
   // нового узла и передаются в качестве параметров в GetTreeNode
   TreeNode<T> *newlptr, *newrptr, *newnode;
   // остановить рекурсивное прохождение при достижении пустого дерева
   if (t == NULL)
      return NULL;
   // CopyTree строит новое дерево в процессе прохождения
   // узлов дерева t. в каждом узле это-го дерева функция
   // CopyTree проверяет наличие левого сына. если он есть,
   // создается его копия. в противном случае возвращается
   // NULL. CopyTree создает копию узла с помощью GetTreeNode
   // и подвешивает к нему копии сыновей.
   if (t->Left() != NULL)
      newlptr = CopyTree(t->Left());
   else
      newlptr = NULL;
   if (t->Right() != NULL)
      newrptr = CopyTree(t->Right());
   else
      newrptr = NULL;
   // построить новое дерево снизу вверх, сначала создавая
   // двух сыновей, а затем их родителя
   newnode = GetTreeNode(t->data, newlptr, newrptr);
   // вернуть указатель на вновь созданное дерево
   return newnode;
}

// использовать обратный алгоритм для прохождения узлов дерева
// и удалить каждый узел при его посещении
template <class T>
void DeleteTree(TreeNode<T> *t)
{
   if(t != NULL)
   {
      DeleteTree(t->Left());
      DeleteTree(t->Right());
      FreeTreeNode(t);
   }
}

#endif // TREENODE_H
