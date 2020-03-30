#ifndef LINKEDLIST_H
#define LINKEDLIST_H

#include <cmath>
#include <iostream>
#include <chrono>
#include <cassert>

using namespace std;

template <typename T>
struct Node
{
    T data;
    Node * next;
    Node(T _data=T(), Node * _next=new Node()):data(_data)
    {
        next=_next;
        _next=NULL;
    }
    //~Node(){delete next;}
};

template <typename T>
class List
{
public:
    List():size(0){head=NULL;}
    ~List()
    {
        Node<T> * prevnode=head;
        Node<T> * curnode=head;
        while(curnode->next!=NULL)
        {
            prevnode=curnode;
            curnode=curnode->next;
            prevnode->next=NULL;
            delete prevnode;
            --size;
        }
    }
    Node<T> * CreateNode(T data);
    int AddNode(T data);
    int AddNode(T data, int position);
    bool RemoveNode(T data);
    bool RemoveNodeP(int position);
    int GetSize() const {return size;}
    void SortAsc();
    void SortDesc();
    void PrintAll()
    {
        Node<T> * curnode=head;
        int i=0;
        if(head!=NULL)
        {
            cout<<"("<<i<<","<<curnode->data<<") ";
            curnode=curnode->next;
            ++i;
        }
        while(curnode!=NULL)
        {
            cout<<"("<<i<<","<<curnode->data<<") ";
            curnode=curnode->next;
            ++i;
        }
    }
    void Display(Node<T> * head)
    {
        if(head==NULL) return;
        else
        {
            cout<<head->data<<" ";
            Display(head->next);
        }
        cout<<endl;
    }
    Node<T> * GetHead() const {return head;}

private:
  int size;
  Node<T> * head;
};

template <typename T>
Node<T> * List<T>::CreateNode(T data)
{
    Node<T> * newnode=new Node<T>(data, NULL);
    return newnode;
}

template <typename T>
int List<T>::AddNode(T data)
{
    Node<T> * newnode=new Node<T>(data, NULL);
    if(head==NULL) head=newnode;
    else
    {
        Node<T> * lastnode=head;
        while(lastnode->next!=NULL) lastnode=lastnode->next;
        lastnode->next=newnode;
    }
    ++size;
    return size-1;
}

template <typename T>
int List<T>::AddNode(T data, int position)
{
    assert(position<=size && "position>size");
    int i=0;
    if(position==0)
    {
        Node<T> * newnode=new Node<T>(data, NULL);
        newnode->next=head;
        head=newnode;
        ++i;
    }
    else
    {
        Node<T> *prevnode=NULL;
        Node<T> *curnode=head;
        int i=0;
        while(i<position)
        {
            prevnode=curnode;
            curnode=curnode->next;
            ++i;
        }
        Node<T> * newnode=new Node<T>(data, NULL);
        prevnode->next=newnode;
        newnode->next=curnode;
        cout<<"i="<<i<<" i-1="<<i-1<<" position="<<position<<endl;
        assert(i==position && "i!=position");
    }
    ++size;
    return i-1;
}

template <typename T>
bool List<T>::RemoveNode(T data)
{
    if(head==NULL) return false;
    int i=0;
    int index=-1;
    Node<T> * prevnode=head;
    Node<T> * curnode=head;
    bool found=false;
    while(curnode!=NULL)
    {
        if(curnode->data==data)
        {
            index=i;
            found=true;
            //prevnode=curnode;
            //curnode=curnode->next;
            cout<<"2 i="<<i<<" found="<<found<<" data="<<data
                <<" prevnode->data="<<prevnode->data
                <<" curnode->data="<<curnode->data<<endl;
            break;
        }
        else
        {
            cout<<"1 i="<<i<<" found="<<found<<" data="<<data
                <<" curnode->data="<<curnode->data<<endl;
            prevnode=curnode;
            curnode=curnode->next;
        }
        ++i;
    }
    if(found)
    {
       if(curnode==head)            head=head->next;
       else if(curnode->next==NULL) prevnode->next=NULL;
       else
       {
            cout<<"A: prevnode->data="<<prevnode->data
                <<" prevnode->next->data="<<prevnode->next->data
                <<" curnode->data="<<curnode->data
                <<" curnode->next->data="<<curnode->next->data
                <<endl;
            prevnode->next=curnode->next;
            cout<<"curnode="<<curnode<<" curnode->data="<<curnode->data<<endl;
       }
       delete curnode;
       --size;
    }
    return true;
}

template <typename T>
bool List<T>::RemoveNodeP(int position)
{
    assert(position<=size-1 && "position>size-1");
    int i=0;
    int index=-1;
    Node<T> * prevnode=head;
    Node<T> * curnode=head;
    bool found=false;
    while(curnode!=NULL)
    {
        if(i==position)
        {
            index=i;
            found=true;
            break;
        }
        prevnode=curnode;
        curnode=curnode->next;
        ++i;
    }
    if(found)
    {
        if(curnode==head)            head=head->next;
        else if(curnode->next==NULL) prevnode->next=NULL;
        else
        {
            cout<<"A: prevnode->data="<<prevnode->data
                <<" prevnode->next->data="<<prevnode->next->data
                <<" curnode->data="<<curnode->data
                <<" curnode->next->data="<<curnode->next->data
                <<endl;
            prevnode->next=curnode->next;
            cout<<"curnode="<<curnode<<" curnode->data="<<curnode->data<<endl;
        }
        delete curnode;
        --size;
    }
    return true;
}

template <typename T>
void List<T>::SortAsc()
{
    Node<T> * curnode=head;
    cout<<"Sort: size="<<size<<endl;
    cout<<"-1@: "<<head->data<<endl;
    Display(head);
    cout<<"0@: "<<head->data<<endl;
    for(int i=0; i<size; ++i)
    {
        curnode=head;
        for(int k=0; k<i; ++k) curnode=curnode->next;
        int ind_min=i;
        T min=curnode->data;
        bool change=false;
        cout<<"ITER#"<<i<<" curnode->data="<<curnode->data<<endl;
        for(int j=i; j<size; ++j)
        {
            cout<<"2@: "<<curnode->data<<endl;
            if(curnode->data<min)
            {
                min=curnode->data;
                ind_min=j;
            }
            curnode=curnode->next;
        }
        cout<<"ind_min="<<ind_min<<std::endl;
        Node<T> * node1=head;
        Node<T> * node2=head;
        for(int i1=0; i1<i; ++i1)       node1=node1->next;
        cout<<"###: node1->data="<<node1->data<<endl;
        for(int i2=0; i2<ind_min; ++i2) node2=node2->next;
        cout<<"###: node1->data="<<node1->data<<" node2->data="<<node2->data<<endl;
        T temp=node1->data;
        node1->data=node2->data;
        node2->data=temp;
        cout<<"i="<<i<<":"<<endl;
        Display(head);
    }
    cout<<endl;
}

template <typename T>
void List<T>::SortDesc()
{
    Node<T> * curnode=head;
    cout<<"Sort: size="<<size<<endl;
    cout<<"-1@: "<<head->data<<endl;
    Display(head);
    cout<<"0@: "<<head->data<<endl;
    for(int i=0; i<size; ++i)
    {
        curnode=head;
        for(int k=0; k<i; ++k) curnode=curnode->next;
        int ind_max=i;
        T max=curnode->data;
        bool change=false;
        cout<<"ITER#"<<i<<" curnode->data="<<curnode->data<<endl;
        for(int j=i; j<size; ++j)
        {
            cout<<"2@: "<<curnode->data<<endl;
            if(curnode->data>max)
            {
                max=curnode->data;
                ind_max=j;
            }
            curnode=curnode->next;
        }
        cout<<"ind_min="<<ind_max<<std::endl;
        Node<T> * node1=head;
        Node<T> * node2=head;
        for(int i1=0; i1<i; ++i1)       node1=node1->next;
        cout<<"###: node1->data="<<node1->data<<endl;
        for(int i2=0; i2<ind_max; ++i2) node2=node2->next;
        cout<<"###: node1->data="<<node1->data<<" node2->data="<<node2->data<<endl;
        T temp=node1->data;
        node1->data=node2->data;
        node2->data=temp;
        cout<<"i="<<i<<":"<<endl;
        Display(head);
    }
    cout<<endl;
}

#endif//LINKEDLIST_H
