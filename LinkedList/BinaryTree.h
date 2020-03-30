#ifndef BINARYTREE_H
#define BINARYTREE_H

#include <climits>

template <typename T>
struct TreeNode
{
  T data;
  TreeNode * left;
  TreeNode * right;
  TreeNode(T val):data(val),left(NULL),right(NULL){}
};

template <typename T>
class BinaryTree
{
public:
  BinaryTree():root(NULL){}
  ~BinaryTree(){DestroyRecursive(root);}
  void DestroyRecursive(TreeNode<T> * node)
  {
    if(node)
    {
        DestroyRecursive(node->left);
        DestroyRecursive(node->right);
        delete node;
    }
  }
  int countNodes();
  void PrintPreOrder(){preorderPrint(root);}
  TreeNode<T> * search(T key);
  TreeNode<T> * newNode(T item);
  void Insert(T val);
  TreeNode<T> *& Delete(T data);
private:
  int countNodes(TreeNode<T> * root);
  void preorderPrint(TreeNode<T> * root);
  TreeNode<T> * search(TreeNode<T> * root, T key);
  void Insert(T val, TreeNode<T> *& node);
  TreeNode<T> *& Delete(TreeNode<T> *& root, T data);
  TreeNode<T> *& FindMin(TreeNode<T> *& root);   
  TreeNode<T> *root;
};

template <typename T>
int BinaryTree<T>::countNodes(TreeNode<T> * root)
{
  if(root == NULL) return 0;// The tree is empty.  It contains no nodes.
  else
  {
    int count = 1;
    count += countNodes(root->left);
    count += countNodes(root->right);
    return count;  // Return the total.
  }
}//end countNodes()

template <typename T>
int BinaryTree<T>::countNodes()
{
  return countNodes(root);
}//end countNodes()

template <typename T>
void BinaryTree<T>::preorderPrint(TreeNode<T> *root)
{
  if(root!=NULL)
  {
    cout<<root->data<< " ";
    preorderPrint(root->left);
    preorderPrint(root->right);
  }
  //cout<<endl;
}

template <typename T>
TreeNode<T>* BinaryTree<T>::search(TreeNode<T>* root, T key) 
{ 
  if (root == NULL || root->data == key) 
     return root; 
  if (root->data < key) 
     return search(root->right, key); 
  return search(root->left, key); 
}

template <typename T>
TreeNode<T>* BinaryTree<T>::search(T key) 
{ 
  return search(root, key); 
}

template <typename T>
TreeNode<T> * BinaryTree<T>::newNode(T item) 
{ 
  TreeNode<T> *temp = new TreeNode<T>(item); 
  temp->left = temp->right = NULL; 
  return temp; 
}

/// Insert a new value into the subtree starting at node
template <typename T>
void BinaryTree<T>::Insert(T val, TreeNode<T> *& node)
{
  //Check if node's value equals val
  //If so, warn the user and then exit function
  if(val == node->data)
  {
    std::cout << "Warning: Value already exists, so nothing will be done." << std::endl;
    return;
  }
  if(val < node->data)
  {
    if(node->left == nullptr)  node->left = new TreeNode(val);
    else                       this->Insert(val, node->left);
  }
  else
  {
    if(node->right == nullptr) node->right = new TreeNode(val);
    else                       this->Insert(val, node->right);
    }
}

template <typename T>
void BinaryTree<T>::Insert(T val)
{
  if(root == nullptr) this->root = new TreeNode(val);
  else                this->Insert(val, this->root);
}

template <typename T>
TreeNode<T> *& BinaryTree<T>::Delete(TreeNode<T> *& root, T data)
{
  //if(root==NULL)              return NULL;
  if(data<root->data)         root->left=Delete(root->left, data);
  else if (data > root->data) root->right=Delete(root->right, data);
  else
  {
    if (root->left==NULL && root->right==NULL)
    {
      delete(root);
      root = NULL;
    }
    else if(root->left==NULL)
    {
      TreeNode<T> *temp = root;
      root = root->right;
      delete temp;
    }
    else if(root->right==NULL)
    {
      TreeNode<T> *temp = root; // save current node as a backup
      root = root->left;
      delete temp;
    }
    else
    {
      TreeNode<T> *temp = FindMin(root->right); // find minimal value of right sub tree
      root->data = temp->data; // duplicate the node
      root->right = Delete(root->right, temp->data); // delete the duplicate node
    }
  }
  return root;
}

template <typename T>
TreeNode<T> *& BinaryTree<T>::Delete(T data)
{
  return Delete(root, data);
}

template <typename T>
TreeNode<T> *& BinaryTree<T>::FindMin(TreeNode<T> *& root)
{
  if(root==NULL)         throw "Min value not found";
  if(root->left != NULL) return FindMin(root->left);//left tree is smaller
  else return root;
}
#endif//BINARYTREE_H
