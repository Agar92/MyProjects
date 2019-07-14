#include <iostream>
#include <random>
#include <algorithm>
#include <chrono>

void BubbleSort(int a[], int size)
{
  for(int i = size - 1; i >= 0; i--)
  {
    for(int j = 0; j < i; j++)
    {
      if(a[j] > a[j+1])
      {
        int tmp = a[j];
        a[j] = a[j + 1];
        a[j + 1] = tmp;
      }
    }
  }
}

void ParallelBubbleSort(int* ar, int size)
{
#pragma omp parallel
  for(int i = 0; i < size; i++)
  {
    int v = 0;
    if(i % 2 == 0)
    {
    #pragma omp for private(v)
      for(int j = 0; j < size; j += 2)
      {
        if(j < size - 1)
        {
          if(ar[j] > ar[j+1])
          {
            v = ar[j];
            ar[j] = ar[j+1];
            ar[j+1] = v;
          }
        }
      }
    }
    else
    {
    #pragma omp for private(v)
      for(int j = 1; j < size; j += 2)
      {
        if(j < size - 1)
        {
          if(ar[j] > ar[j+1])
          {
            v = ar[j];
            ar[j] = ar[j+1];
            ar[j+1] = v;
          }
        }
      }
    }
  }
}

void InsertionSort(int data[], int size)
{
  int key = 0;
  int i = 0;
  for(int j = 1; j < size; j++){
    key = data[j];
    i = j-1;
    while(i >= 0 && data[i] > key){
      data[i+1] = data[i];
      i = i-1;
      data[i+1]=key;
    }
  }
}

void InsertionSortSimple(int data[], int size)
{
  int key = 0;
  int j = 0;
  for(int i = 1; i < size; i++){
    j = i-1;
    while(j >= 0 && data[j] > data[j+1]){
      std::swap(data[j], data[j+1]);
      --j;
    }
  }
}

//сортировка Шелла
void ShellSort (int data[], int size){
  int h, i, j;//h - distance between elements
  for(h = size/2 ; h > 0 ; h = h/2)
  {
    for(i = 0 ; i < size-h ; i++)
    {
      for(j = i ; j >= 0 ; j = j - h)
      {
        //std::cout<<"h="<<h<<" i="<<i
        //         <<" j="<<j<<std::endl;
        if(data[j] > data[j+h])
          std::swap(data[j], data[j+h]);
        else j = 0;
      }
    }
  }
}
//gnome sort begin:
// 4 7 2 8 9 1 5 7 1 - original array.
// 4 2 7 8 9 1 5 7 1
// 2 4 7 8 9 1 5 7 1
// 2 4 7 8 1 9 5 7 1
// 2 4 7 1 8 9 5 7 1
// 2 4 1 7 8 9 5 7 1
// 2 1 4 7 8 9 5 7 1
// 1 2 4 7 8 9 5 7 1
// 1 2 4 7 8 5 9 7 1
// 1 2 4 7 5 8 9 7 1
// 1 2 4 5 7 8 9 7 1
// 1 2 4 5 7 8 7 9 1
// 1 2 4 5 7 7 8 9 1
// 1 2 4 5 7 7 8 1 9
// 1 2 4 5 7 7 1 8 9
// 1 2 4 5 7 1 7 8 9
// 1 2 4 5 1 7 7 8 9
// 1 2 4 1 5 7 7 8 9
// 1 2 1 4 5 7 7 8 9
// 1 1 2 4 5 7 7 8 9 - sorted array.
//gnome sort end.
//гномья сортировка
void GnomeSort(int A[], int size)
{
  for(int i = 0; i < size-1; ++i)
  {
    if(A[i] > A[i + 1])
    {
      std::swap(A[i], A[i+1]);
      if(i != 0) i -= 2; //вычитается два и потом прибавляется один
    }
  }
}

//сортировка слиянием:
// Merges two subarrays of arr[].
// First subarray is arr[l..m]
// Second subarray is arr[m+1..r]
void Merge(int arr[], int l, int m, int r)
{
  int i, j, k;
  int n1 = m - l + 1;
  int n2 =  r - m;
  /* create temp arrays */
  int L[n1], R[n2];
  /* Copy data to temp arrays L[] and R[] */
  for(i = 0; i < n1; i++) L[i] = arr[l + i];
  for(j = 0; j < n2; j++) R[j] = arr[m + 1+ j];
  /* Merge the temp arrays back into arr[l..r]*/
  i = 0; // Initial index of first subarray
  j = 0; // Initial index of second subarray
  k = l; // Initial index of merged subarray
  while(i < n1 && j < n2)
  {
    if(L[i] <= R[j])
    {
      arr[k] = L[i];
      i++;
    }
    else
    {
      arr[k] = R[j];
      j++;
    }
    k++;
  }
  /* Copy the remaining elements of L[], if there are any */
  while(i < n1)
  {
    arr[k] = L[i];
    i++;
    k++;
  }
  /* Copy the remaining elements of R[], if there are any */
  while(j < n2)
  {
    arr[k] = R[j];
    j++;
    k++;
  }
}

/* l is for left index and r is right index of the
   sub-array of arr to be sorted */
void MergeSort(int arr[], int l, int r)
{
  if(l < r)
  {
    // Same as (l+r)/2, but avoids overflow for
    // large l and h
    int m = l+(r-l)/2;
    // Sort first and second halves
    MergeSort(arr, l, m);
    MergeSort(arr, m+1, r);
    Merge(arr, l, m, r);
  }
}

//быстрая сортировка (опорный элемент - последний):

//A utility function to swap two elements
void swap(int* a, int* b)
{
  int t = *a;
  *a = *b;
  *b = t;
}

/* This function takes last element as pivot, places
the pivot element at its correct position in sorted
array, and places all smaller (smaller than pivot)
to left of pivot and all greater elements to right
of pivot */
int Partition (int arr[], int low, int high)
{
  int pivot = arr[high]; // pivot - опорный элемент
  int i = (low - 1); // Index of smaller element
  for(int j = low; j <= high - 1; j++)
  {
    // If current element is smaller than or equal to pivot
    if(arr[j] <= pivot)
    {
      i++; // increment index of smaller element
      swap(&arr[i], &arr[j]);
    }
  }
  swap(&arr[i + 1], &arr[high]);
  return (i + 1);
}

/* The main function that implements QuickSort
arr[] --> Array to be sorted,
low --> Starting index,
high --> Ending index */
void QuickSort(int arr[], int low, int high)
{
    if(low < high)
    {
        /* pi is partitioning index, arr[p] is now
        at right place */
        int pi = Partition(arr, low, high);

        // Separately sort elements before
        // partition and after partition
        QuickSort(arr, low, pi - 1);
        QuickSort(arr, pi + 1, high);
    }
}

//сортировка кучей:
//To heapify a subtree rooted with node i which is
//an index in arr[]. size is size of heap (size of array).
void Heapify(int arr[], int n, int i)//делаем 2-ичную кучу
{
  int largest = i; //Initialize largest as root
  int l = 2*i + 1; //left = 2*i + 1
  int r = 2*i + 2; //right = 2*i + 2
  //If left child is larger than root
  if(l < n && arr[l] > arr[largest]) largest = l;
  //If right child is larger than largest so far
  if(r < n && arr[r] > arr[largest]) largest = r;
  //If largest is not root
  if(largest != i)
  {
    std::swap(arr[i], arr[largest]);
    //Recursively heapify the affected sub-tree
    Heapify(arr, n, largest);
  }
}

//main function to do heap sort
void HeapSort(int arr[], int n)
{
  //Build heap (rearrange array)
  for(int i = n/2 - 1; i >= 0; i--) Heapify(arr, n, i);
  //One by one extract an element from heap
  for(int i=n-1; i>=0; i--)
  {
    //Move current root to end
    std::swap(arr[0], arr[i]);
    //call max heapify on the reduced heap
    Heapify(arr, i, 0);
  }
}

//**************************************************************//
int main(int argc, char *argv[])
{
    const int N=50000;//не больше 50000.
//N=50000:   T1=6558 ms T2=16765 ms T3=10576 ms T4=5689 ms
//           T5=16 ms T6=10 ms T7=129 ms T8=35 ms.
    int a1[N]{0};
    for(int i=0; i<N; ++i) a1[i]=rand()%100+1;//1-:-100 - possible values.
    //std::cout<<"Before sort 1:"<<std::endl;
    //for(int i=0; i<N; ++i) std::cout<<a1[i]<<" ";
    //std::cout<<std::endl;
//===========================================================================//
//сортировка выбором                                                         //
//===========================================================================//
  /*********** Начало сортировки **************/
    auto start1 = std::chrono::steady_clock::now();
    for (int i = 0; i < N; i++)
    {
      int min = i;
        for (int j = i + 1; j < N; j++) min = ( a1[j] < a1[min] ) ? j : min;
      if (i != min) std::swap(a1[min] ,a1[i]);
    }
    auto end1 = std::chrono::steady_clock::now();
  /*********** Конец сортировки **************/
    //std::cout<<"After sort 1:"<<std::endl;
    //for(int i=0; i<N; ++i) std::cout<<a1[i]<<" ";
    //std::cout<<std::endl;
//===========================================================================//
//пузырьковая сортировка                                                     //
//===========================================================================//
    int a2[N]{0};
    for(int i=0; i<N; ++i) a2[i]=rand()%100+1;//1-:-100 - possible values.
    //std::cout<<"Before sort 2:"<<std::endl;
    //for(int i=0; i<N; ++i) std::cout<<a2[i]<<" ";
    //std::cout<<std::endl;
    auto start2 = std::chrono::steady_clock::now();
    for (int i = N - 1; i >= 0; i--)
    {
      for (int j = 0; j < i; j++)
      {
        if (a2[j] > a2[j+1])
        {
          int tmp = a2[j];
          a2[j] = a2[j + 1];
          a2[j + 1] = tmp;
        }
      }
    }
    auto end2 = std::chrono::steady_clock::now();
    //std::cout<<"After sort 2:"<<std::endl;
    //for(int i=0; i<N; ++i) std::cout<<a2[i]<<" ";
    //std::cout<<std::endl;
//===========================================================================//
//параллельная OpenMP пузырьковая сортировка                                 //
//===========================================================================//
    int a3[N]{0};
    for(int i=0; i<N; ++i) a3[i]=rand()%100+1;//1-:-100 - possible values.
    //std::cout<<"Before sort 3:"<<std::endl;
    //for(int i=0; i<N; ++i) std::cout<<a3[i]<<" ";
    //std::cout<<std::endl;
    auto start3 = std::chrono::steady_clock::now();
    ParallelBubbleSort(a3, N);
    auto end3 = std::chrono::steady_clock::now();
    //std::cout<<"After sort 3:"<<std::endl;
    //for(int i=0; i<N; ++i) std::cout<<a3[i]<<" ";
    //std::cout<<std::endl;
//===========================================================================//
//сортировка вставками                                                       //
//===========================================================================//

    int a4[N]{0};
    for(int i=0; i<N; ++i) a4[i]=rand()%100+1;//1-:-100 - possible values.
    //std::cout<<"Before sort 4:"<<std::endl;
    //for(int i=0; i<N; ++i) std::cout<<a4[i]<<" ";
    //std::cout<<std::endl;
    auto start4 = std::chrono::steady_clock::now();
    InsertionSort(a4, N);
    auto end4 = std::chrono::steady_clock::now();
    //std::cout<<"After sort 4:"<<std::endl;
    //for(int i=0; i<N; ++i) std::cout<<a4[i]<<" ";
    //std::cout<<std::endl;
//===========================================================================//
//сортировка Шелла                                                           //
//===========================================================================//
    int a5[N]{0};
    for(int i=0; i<N; ++i) a5[i]=rand()%100+1;//1-:-100 - possible values.
    //std::cout<<"Before sort 5:"<<std::endl;
    //for(int i=0; i<N; ++i) std::cout<<a5[i]<<" ";
    //std::cout<<std::endl;
    auto start5 = std::chrono::steady_clock::now();
    ShellSort(a5, N);
    auto end5 = std::chrono::steady_clock::now();
    //std::cout<<"After sort 5:"<<std::endl;
    //for(int i=0; i<N; ++i) std::cout<<a5[i]<<" ";
    //std::cout<<std::endl;
//===========================================================================//
//сортировка слиянием                                                        //
//===========================================================================//
    int a6[N]{0};
    for(int i=0; i<N; ++i) a6[i]=rand()%100+1;//1-:-100 - possible values.
    //std::cout<<"Before sort 6:"<<std::endl;
    //for(int i=0; i<N; ++i) std::cout<<a6[i]<<" ";
    //std::cout<<std::endl;
    auto start6 = std::chrono::steady_clock::now();
    MergeSort(a6, 0, N-1);
    auto end6 = std::chrono::steady_clock::now();
    //std::cout<<"After sort 6:"<<std::endl;
    //for(int i=0; i<N; ++i) std::cout<<a6[i]<<" ";
    //std::cout<<std::endl;
//===========================================================================//
//быстрая сортировка                                                         //
//===========================================================================//
    int a7[N]{0};
    for(int i=0; i<N; ++i) a7[i]=rand()%100+1;//1-:-100 - possible values.
    //std::cout<<"Before sort 7:"<<std::endl;
    //for(int i=0; i<N; ++i) std::cout<<a7[i]<<" ";
    //std::cout<<std::endl;
    auto start7 = std::chrono::steady_clock::now();
    QuickSort(a7, 0, N-1);
    auto end7 = std::chrono::steady_clock::now();
    //std::cout<<"After sort 7:"<<std::endl;
    //for(int i=0; i<N; ++i) std::cout<<a7[i]<<" ";
    //std::cout<<std::endl;
//===========================================================================//
//сортировка кучей                                                           //
//===========================================================================//
    int a8[N]{0};
    for(int i=0; i<N; ++i) a8[i]=rand()%100+1;//1-:-100 - possible values.
    //std::cout<<"Before sort 8:"<<std::endl;
    //for(int i=0; i<N; ++i) std::cout<<a8[i]<<" ";
    //std::cout<<std::endl;
    auto start8 = std::chrono::steady_clock::now();
    HeapSort(a8, N);
    auto end8 = std::chrono::steady_clock::now();
    //std::cout<<"After sort 8:"<<std::endl;
    //for(int i=0; i<N; ++i) std::cout<<a8[i]<<" ";
    //std::cout<<std::endl;

    auto const T1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1-start1).count();
    auto const T2 = std::chrono::duration_cast<std::chrono::milliseconds>(end2-start2).count();
    auto const T3 = std::chrono::duration_cast<std::chrono::milliseconds>(end3-start3).count();
    auto const T4 = std::chrono::duration_cast<std::chrono::milliseconds>(end4-start4).count();
    auto const T5 = std::chrono::duration_cast<std::chrono::milliseconds>(end5-start5).count();
    auto const T6 = std::chrono::duration_cast<std::chrono::milliseconds>(end6-start6).count();
    auto const T7 = std::chrono::duration_cast<std::chrono::milliseconds>(end7-start7).count();
    auto const T8 = std::chrono::duration_cast<std::chrono::milliseconds>(end8-start8).count();
    std::cout<<"T1="<<T1<<" ms T2="<<T2<<" ms T3="<<T3<<" ms T4="<<T4<<" ms"
             <<" T5="<<T5<<" ms T6="<<T6<<" ms T7="<<T7<<" ms T8="<<T8<<" ms"<<std::endl;

    return 0;
}
