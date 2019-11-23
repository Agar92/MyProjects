#include <chrono>
#include <iostream>
#include "unistd.h"
#include <typeinfo>
#include <cmath>
#include <fstream>
#include <iomanip>

//#include "T3Globals.h"
#include "T3Defs.h"
#include "T3RNG.h"

//
#include "T2NSGangular_RW.hh"
//

#include <iostream>
#include "unistd.h"
#include <iterator>

//#include "T3NSGangular_RW.h"

using namespace t3;



#ifdef OPENACC
#include <accelmath.h>
#include <openacc.h>
#ifdef CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif
#else
#include <omp.h>
  struct float3{
    float x, y, z;
  };
  struct float4{
    float x, y, z, w;
  };
#endif



//#include "T3DataHolder.h"

//build options on i7:
// cmake3 . -DCMAKE_C_COMPILER=pgcc -DCMAKE_CXX_COMPILER=pgc++
//-DCMAKE_C_FLAGS="-acc -Minfo=all -mcmodel=medium -ta=tesla:cc30
//-Mnollvm -Minline -Mcuda=cuda10.1"
//-DCMAKE_CXX_STANDARD=17 -DACC=ON -DCUDA=ON

//build options on KNL:
//cmake . -DCMAKE_C_COMPILER=pgcc -DCMAKE_CXX_COMPILER=pgc++
//-DCMAKE_C_FLAGS="-acc -Minfo=acc -mcmodel=medium -ta=tesla:cc70
//-tp=haswell -Mnollvm -Minline -Mcuda=cuda10.1"
//-DCMAKE_CXX_FLAGS="-acc -Minfo=acc -mcmodel=medium -ta=tesla:cc70
//-tp=haswell -Mnollvm -Minline -Mcuda=cuda10.1"
//-DCMAKE_CXX_STANDARD=17 -DACC=ON -DCUDA=ON




////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!//
//*************************************************************************//
//Before beginning to work with this program every day, You SHOULD MAKE IN THE TERMINAL
//$ source ~/TPT_BUILD_DIRECTORY_18_09_2019/install/sourcemegcc-xsystem-g4v10.5.1-tptmaster-Debug.sh
//IT SETS THE ENVIRONMENT VARIABLES OF T2_DATA AND OTHERS, NECESSARY FOR READING DATABASES FROM
///home/70-gaa/TPT_BUILD_DIRECTORY_18_09_2019/data/T2_DATA.
//*************************************************************************//




//TOTAL NUMBER OF INJECTED PARTICLES:
const int NPARTICLES=100000000;//1 million

//the indexes in the rw.at(0).at(index) database,
//which i will use for taking partial sums at necessary energies:
//10   MeV index=66
//5    MeV index=56
//2    MeV index=39
//1    MeV index=29
//0.1  MeV index=11
//0.01 MeV index=2


//!!!!!!!!!!!!!!!!!!!!!!!!!!!1//
//const double Energy=0.01;
double Ee=2.0;//0.1;
//!!!!!!!!!!!!!!!!!!!!!!!!!!!1//


//const int PSindex=2;//39;//66;

//the step is 2/128 and there are 129 points in Ox (cos(theta_CM)) array:
//from -1 to 1, there are 129 points:

//Found that PAW does not allow to plot more than 1024 points
//when i try to plot 2048 bins, it plots on ly the first 1024
//=>i can see only the left zoom and the left side of the center xoom.
constexpr int Nbin=2048;//1024;//2048;//1024;
constexpr int NpointsInPS=129;

const double delta=static_cast<double>(2.0)/128;

constexpr int NpointsInCosCM=Nbin+1;
//cos(theta_CM):
const double deltacoscm=2.0/Nbin;
double coscm[129];
//differentiated partial sums (dps):
//there are  NpointsInCosCM-1 elements because can not take the derivative of x=1.0.
double dpsmain[NpointsInPS-1]{0.0};
//dpsexp:
double dpsexp[NpointsInPS-1]{0.0};

//Ox (cos(theta_cm)) in the middle of bins (128 bins=>128 elements in Ox):
double Ox[NpointsInCosCM-1]{0.0};

const int CHANNELS=2;//2 channels: exponent and main distribution.

int HistogramNiCosCM[NpointsInCosCM-1]{0};


int InHistogramRange=0;
int OutOfHistogramRange=0;


//counters:
int COUNTCHANNEL1=0;
int COUNTCHANNEL1INBIN=0;
int COUNTCHANNEL1OUTOFBIN=0;

int COUNTCHANNEL2=0;
int COUNTCHANNEL2INBIN=0;
int COUNTCHANNEL2OUTOFBIN=0;

//end of counters.

//**********************************************************************//
//here the partial sums from ~/TPT_BUILD_DIRECTORY_18_09_2019/data/T2_DATA/
//should be stored.

//number of partial sums, values of E, values of pr and sl:
const int N1=67;
//number of elements in partial sums (in SDI's partial sums there are only 127 elements,
//but i will add 0.0 and 1.0 and have 129 elements):
const int N2=129;
double PS[N1][N2]{0.0};
double En[N1]{0.0};
double PR[N1]{0.0};
double SL[N1]{0.0};
//**********************************************************************//
//this is a function for filling PS, E, PR, SL arrays from SDI's
//database at ~/TPT_BUILD_DIRECTORY_18_09_2019/data/T2_DATA/:

/*
void FillT2_DATA_DD_Database()
{
  //the object for reading the database:
  T2NSGangular_RW rw;
  //target deuteron Z=1
  static constexpr auto tgZ = 1;
  //target deuteron A=2
  static constexpr auto tgA = 2;
  //reaction product=n+3He
  //this is a neutron from a reaction product:
  static constexpr auto sPDG = 2112;
  //1000*Z+A//inc deuteron
  static constexpr auto incZA = 1002;
  //MT index of reaction in ENDF "Mt=50".
  static constexpr auto rid = "50";
  //load the data from the binary file:
  rw.load_binary(tgZ, tgA, rid, sPDG, incZA);
  //print the data records (e, partical sums):
  //std::cout << rw << std::endl;
  
  std::cout<<"rw size="<<rw.size()<<std::endl;
  std::cout<<"SEE:"<<std::endl;

  std::cout<<"SIZE="<<rw.size()<<" SIZE OF THE RECORD="<<rw.at(0).size()<<std::endl;

  for(int i=0; i<N1; ++i)
  {
    //print the node:
    //std::cout<<rw.at(0).at(i)<<std::endl;
    //the energy of the partial sums distribution:
    En[i]=rw.at(0).at(i).Get_E();
    //get partial sums array:
    std::array<double, T2NSGangular_RWnode::_num_point> ps0=rw.at(0).at(i).Get_V();
    PS[i][NpointsInCosCM-1]=0.0;
    for(int j=0; j<T2NSGangular_RWnode::_num_point; ++j) PS[i][j+1]=ps0.at(j);
    PS[i][NpointsInCosCM-1]=1.0;
    //exponent contribution:
    PR[i]=rw.at(0).at(i).Get_pr();
    //exponent slope: A*e^(-B*x), here B is the slope of the exponent.
    SL[i]=rw.at(0).at(i).Get_sl();
  }

  std::cout<<"Check the arrays:"<<std::endl;
  std::cout<<"Check array En:"<<std::endl;
  for(int i=0; i<N1; ++i) std::cout<<En[i]<<" ";
  std::cout<<std::endl;

  std::cout<<"Check array PR:"<<std::endl;
  for(int i=0; i<N1; ++i) std::cout<<PR[i]<<" ";
  std::cout<<std::endl;

  std::cout<<"Check array SL:"<<std::endl;
  for(int i=0; i<N1; ++i) std::cout<<SL[i]<<" ";
  std::cout<<std::endl;
  
}
*/




int SaveDDPartialSumsToAFile()
{
  std::ofstream out_stream("t2ddpartialsums.dat", std::ios::binary);
  if(out_stream.good())
  {
    for(int i=0; i<N1; ++i)
    {
      out_stream.write((const char*) &En[i], sizeof(double) );
      out_stream.write((const char*) &PR[i], sizeof(double) );
      out_stream.write((const char*) &SL[i], sizeof(double) );
      out_stream.write((const char*) &PS[i][0], N2 * sizeof(double));
    }
  }
  else
  {
    std::cout<<"***Warning: Bad stream!"<<std::endl;
    return 0;
  }
  out_stream.close();
  return 1;
}

int ReadDDPartialSumsFromAFile()
{
  std::ifstream in_stream("t2ddpartialsums.dat", std::ios::binary);
  if(in_stream.good())
  {
    for(int i=0; i<N1; ++i)
    {
      in_stream.read((char*) &En[i], sizeof(double) );
      in_stream.read((char*) &PR[i], sizeof(double) );
      in_stream.read((char*) &SL[i], sizeof(double) );
      in_stream.read((char*) &PS[i][0], N2 * sizeof(double));
    }
  }
  else
  {
    std::cout<<"***Warning: Bad stream!"<<std::endl;
    return 0;
  }
  in_stream.close();
  return 1;
}


//**********************************************************************//



//there are totally 128 bins.
//the length from -1 to 1 is 2
//points -1 and 1 are included in Ox (cos(theta_CM)),
void MakeEquidistantBinsCosThetaCM()
{
#ifdef OPENACC
#pragma acc parallel loop gang vector present(coscm) copy(delta)
#endif
  for(int i=0; i<129; ++i) coscm[i]=-1.0+i*delta;
}

void FillDifferentiatedPartialSums(double * ps/*read from T2_DATA/ partial sums - 67 elements in ps*/,
                                   double PS_CONTRIBUTION/*contribution of main partial sums distribution=1-pr, where pr is the contribution of the exponent*/
                                  )
{
#ifdef OPENACC
#pragma acc parallel loop gang vector present(dpsmain,ps,CosCM) 
#endif
  for(int m=0; m<NpointsInPS-1; ++m)
  {
    //differentiated partial sums:
    //std::cout<<"m="<<m<<std::endl;
    dpsmain[m]=PS_CONTRIBUTION*((ps[m+1]*ps[m+1]-ps[m]*ps[m]))/delta;

    std::cout<<"delta="<<delta<<" PS_CONTRIBUTION="<<PS_CONTRIBUTION
             <<std::endl;
    for(int j=0; j<129; ++j) std::cout<<ps[j]<<" ";
    std::cout<<std::endl;
    //sleep(1);
    
    //std::cout<<"ps.at(m+1)="<<ps.at(m+1)<<" ps.at(m)="<<ps.at(m)<<" ps.at(m+1)-ps.at(m)="<<ps.at(m+1)-ps.at(m)<<" CosCM.at(m+1)-CosCM.at(m)="<<CosCM.at(m+1)-CosCM.at(m)
    //         <<" dpsmain.at(m)="<<dpsmain.at(m)<<std::endl;
  }
}

//in the middle of bins:
void FillOxAxisCosCM() {
  //in the middle of the bin:
  //std::cout<<"Filling Ox:"<<std::endl;

#ifdef OPENACC
#pragma acc parallel loop gang present(Ox,coscm) copy(deltacoscm)
#endif
  for(int m=0; m<NpointsInCosCM; ++m)
  {
    Ox[m]=coscm[0]+(m+0.5)*deltacoscm;
  }
}

void Histogram_OxAxisCosCM() {
  std::ofstream foutne_theta;
  foutne_theta.open("ne_thetacoscm.dat");
#ifdef OPENACC
#pragma acc update host(Ox)
#endif
  for(int m=0; m<NpointsInCosCM-1; ++m)
  {
    foutne_theta<<m<<"   ";
    //in the middle of the bin:
    foutne_theta<<std::setw(8)<<Ox[m]<<"   ";
    foutne_theta<<std::endl;
  }
  foutne_theta.close();
}

void Histogram_DifferentiatedPartialSums() {
  std::ofstream foutne_theta;
  foutne_theta.open("ne_thetadps.dat");
#ifdef OPENACC
#pragma acc update host(dpsmain)
#endif
  for(int m=0; m<NpointsInPS-1; ++m)
  {
    const double deltax=2.0/128;
    const double x=-1.0+deltax*(m+0.5);
    foutne_theta<<m<<"   "<<std::setw(8)<<x<<"   "<<dpsmain[m]<<"   ";
    foutne_theta<<std::endl;
  }
  foutne_theta.close();
}

void FillExponent(double sl, double pr/*pr*/) {

#ifdef OPENACC
#pragma acc parallel loop gang vector present(Ox,dpsexp) copy(sl,pr) 
#endif
  for(int m=0; m<NpointsInPS-1; ++m)
  {
    const double deltax=2.0/128;
    const double x=-1.0+deltax*(m+0.5);
    const double norm=sl*(1.0-exp(-2.0/sl));
    /*
      std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
      std::cout<<"slope="<<sl<<std::endl;
      std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
    */
    const double fexp=pr*exp(-(1-x)/sl)/norm;
    dpsexp[m]=fexp;
    //std::cout<<"x="<<Ox.at(m)<<" norm="<<norm<<" fexp="<<fexp
    //         <<" dpsexp.at(m)="<<dpsexp.at(m)<<std::endl;
  }
}

void Histogram_Exponent() {
  std::ofstream foutne_theta;
  foutne_theta.open("ne_thetaexp.dat");
#ifdef OPENACC
#pragma acc update host(dpsexp)
#endif
  for(int m=0; m<NpointsInPS-1; ++m)
  {
    foutne_theta<<m<<"   "<<std::setw(8)<<dpsexp[m]<<"   ";
    foutne_theta<<std::endl;
  }
  foutne_theta.close();
}

void FillHistogramNi()
{
  std::ofstream foutne_theta;
  foutne_theta.open("ne_thetani.dat");
#ifdef OPENACC
#pragma acc update host(HistogramNiCosCM)
#endif
  for(int m=0; m<NpointsInCosCM-1; ++m)
  {
    foutne_theta<<m<<"   ";
    //1)number of events in each bin, normalized by the total number of events.
    //if the surface under the graphic of the function is 1, then this is equal to
    //the surface under the function graphic in each bin:
    //2)Since we take the normalized partial sums from T2_DATA/ and
    //differentiate them (divide by deltacoscm=2.0/Nbin),
    //here we also need to divide normalized Ni/Ntot by deltacoscm=2.0/Nbin,
    //to obtain a derivative. 
    //\\//const double Hm=static_cast<double>(HistogramNiCosCM.at(j).at(m))/NPARTICLES/deltacoscm;
    foutne_theta<<std::setw(8)<<static_cast<double>(HistogramNiCosCM[m])/NPARTICLES/deltacoscm<<"   ";
    foutne_theta<<std::endl;
  }
  foutne_theta.close();
}


int BinarySearch(double array[NpointsInCosCM], double key)
{
  int l=0, r=NpointsInCosCM-1;
  int mid;
  int i=0;
  while(l!=r-1)
  {
    //std::cout<<"Iteration #"<<i<<std::endl;
    mid=(l+r)/2;
    /*
    std::cout<<"l="<<l<<" r="<<r<<" mid="<<mid
             <<" array[mid]="<<array[mid]<<" key="<<key<<std::endl;
    */
    if(array[mid]==key) return mid;
    if(array[mid]>key)  r=mid;
    else                l=mid;

    //std::cout<<"l="<<l<<" r="<<r<<" mid="<<mid<<std::endl;

    ++i;
  }
  //????????????//
  if(key>=array[NpointsInCosCM-1]) return r;
  return l;
}

int BinarySearch(double * array, const int N, double key)
{
  int l=0, r=N-1;
  int mid;
  int i=0;
  while(l!=r-1)
  {
    //std::cout<<"Iteration #"<<i<<std::endl;
    mid=(l+r)/2;
    /*
    std::cout<<"l="<<l<<" r="<<r<<" mid="<<mid
             <<" array[mid]="<<array[mid]<<" key="<<key<<std::endl;
    */
    if(array[mid]==key) return mid;
    if(array[mid]>key)  r=mid;
    else                l=mid;

    //std::cout<<"l="<<l<<" r="<<r<<" mid="<<mid<<std::endl;

    ++i;
  }

  //\\//std::cout<<"key="<<key<<" array[l]="<<array[l]<<" array[r]="<<array[r]<<std::endl;
  
  //????????????//
  if(key>=array[N-1]) return r;
  return l;
}



double GetCosThetaCMFromDDInelasticPartialSums(double E/*Energy*/, double RNDnumber/*random number from 0 to 1*/)
{
  //std::cout<<"Enter GetBinNumberForT2_DATA_DD_FROM_PARTIAL_SUMS(): E="<<E<<" RNDnumber="<<RNDnumber<<std::endl;
  if(0.0001<=E && E<=10.0)
  {
    //std::cout<<"A"<<std::endl;
    if(RNDnumber==1) return 1.0;
    int binE=BinarySearch(En, N1, E);
    double costhetacm=0.0;
    //interpolate partial sums:
    //if E<10 MeV
    if(binE!=N1-1)
    {
      //std::cout<<"AA"<<std::endl;
      const int binPSEb=BinarySearch(PS[binE], N2, RNDnumber);
      const int binPSEb1=BinarySearch(PS[binE+1], N2, RNDnumber);      
      //double

      /*
      std::cout<<"Check 3 binary searches:"<<std::endl;
      std::cout<<"binE="<<binE<<" RNDnumber="<<RNDnumber<<" E="<<E<<" Eb="<<En[binE]<<" Eb1="<<En[binE+1]<<std::endl;
      std::cout<<"binPSEb="<<binPSEb<<" "<<PS[binE][binPSEb]<<" "<<PS[binE][binPSEb+1]<<std::endl;
      std::cout<<"binPSEb1="<<binPSEb1<<" "<<PS[binE+1][binPSEb1]<<" "<<PS[binE+1][binPSEb1+1]<<std::endl;

      std::cout<<"Check PS[binE] last bin #"<<N2-1<<" "<<PS[binE][N2-1]<<" "<<PS[binE][N2]<<std::endl;
      std::cout<<"Check PS[binE+1] last bin #"<<N2-1<<" "<<PS[binE+1][N2-1]<<" "<<PS[binE+1][N2]<<std::endl;
      
      
      std::cout<<"Check PS[binE]:"<<std::endl;
      for(int i=0; i<N2; ++i) std::cout<<PS[binE][i]<<" ";
      std::cout<<std::endl;
      
      std::cout<<"Check PS[binE+1]:"<<std::endl;
      for(int i=0; i<N2; ++i) std::cout<<PS[binE+1][i]<<" ";
      std::cout<<std::endl;
      */

      //const double deltaX1=Ox[binPSEb+1]-Ox[binPSEb];  //=deltacoscm
      //const double deltaX2=Ox[binPSEb1+1]-Ox[binPSEb1];//=deltacoscm
      const double dPS1=PS[binE][binPSEb+1]*PS[binE][binPSEb+1]-PS[binE][binPSEb]*PS[binE][binPSEb];
      const double dPS2=PS[binE+1][binPSEb1+1]*PS[binE+1][binPSEb1+1]-PS[binE+1][binPSEb1]*PS[binE+1][binPSEb1];
      const double DeltaXCosCM=(double)2.0/128;
      const double costhetacm1=coscm[binPSEb]+(RNDnumber*RNDnumber-PS[binE][binPSEb]*PS[binE][binPSEb])*DeltaXCosCM/dPS1;
      const double costhetacm2=coscm[binPSEb1]+(RNDnumber*RNDnumber-PS[binE+1][binPSEb1]*PS[binE+1][binPSEb1])*DeltaXCosCM/dPS2;
      //interpolation of cos(theta_cm) for Ei and Ei+1:
      costhetacm=costhetacm1+(costhetacm2-costhetacm1)*(E-En[binE])/(En[binE+1]-En[binE]);

      /*
      std::cout<<"E="<<E<<" RNDnumber="<<RNDnumber<<" binE="<<binE<<std::endl;
      std::cout<<"Energies:"<<std::endl;
      for(int g=0; g<N1; ++g) std::cout<<En[g]<<" ";
      std::cout<<std::endl;

      std::cout<<"DeltaXCosCM="<<DeltaXCosCM<<std::endl;

      std::cout<<"binPSEb="<<binPSEb<<" binPSEb1="<<binPSEb1<<" R="<<RNDnumber<<std::endl;
      std::cout<<"binPSEb="<<binPSEb<<" "<<PS[binE][binPSEb]<<" "<<PS[binE][binPSEb+1]<<std::endl;
      std::cout<<"binPSEb1="<<binPSEb1<<" "<<PS[binE+1][binPSEb1]<<" "<<PS[binE+1][binPSEb1+1]<<std::endl;
      //sleep(1);
      */

      /*
      std::cout<<"deltacosthetacm="<<deltacoscm<<std::endl;
      std::cout<<"!!!costhetacm="<<costhetacm<<" costhetacm1="<<costhetacm1<<" costhetacm2="<<costhetacm2
               <<" coscm[binPSEb]="<<coscm[binPSEb]<<" coscm[binPSEb1]="<<coscm[binPSEb1]
               <<" 1-coscm[binPSEb]="<<1.0-coscm[binPSEb]<<" 1-coscm[binPSEb1]="<<1.0-coscm[binPSEb1]<<std::endl;
      std::cout<<"RNDnumber="<<RNDnumber<<std::endl;
      */
      //sleep(2);
      
    }
    else//if E==10MeV
    {
      //std::cout<<"AB"<<std::endl;
      const int binPSEb=BinarySearch(PS[binE], N2, RNDnumber);
      //const double deltaX1=Ox[binPSEb+1]-Ox[binPSEb];//=deltacoscm
      const double dPS1=PS[binE][binPSEb+1]-PS[binE][binPSEb];
      //interpolation of cos(theta_cm) for Ei=10 MeV:
      const double DeltaXCosCM=(double)2.0/128;
      costhetacm=coscm[binPSEb]+(RNDnumber-PS[binE][binPSEb])*DeltaXCosCM/dPS1;
    }

    /*//THIS DOES NOT WORK ON GPU!!!
    //THERE ARE 2 WAYS HOW TO INTERPOLATE cos(theta_cm):
    //1) AS SDI ADVISED
    //2) INTERPOLATION BY 4 POINTS, AS KMV TOLD.
    double tempPS[N2];
    for(int j=0; j<N2; ++j)
    {
      tempPS[j]=PS[binE][j]+(E-En[binE])*(PS[binE+1][j]-PS[binE][j])/(En[binE+1]-En[binE]);
    }
    */
    return costhetacm;
  }
  else
  {
    //std::cout<<"B"<<std::endl;
    ///////std::cout<<"***ERROR: E is out of energy range!"<<std::endl;
    return -100.0;
  }  
}


double GetPR(double E)
{
  if(0.0001<=E && E<=10.0)
  {
    int binE=BinarySearch(En, N1, E);
    //interpolate PR:
    double pr=PR[binE]+(PR[binE+1]-PR[binE])*(E-En[binE])/(En[binE+1]-En[binE]);
    return pr;
  }
  else
  {
    //////std::cout<<"***ERROR: E is out of energy range!"<<std::endl;
    return -100.0;
  }  
}

double GetSL(double E)
{
  if(0.0001<=E && E<=10.0)
  {
    int binE=BinarySearch(En, N1, E);
    //interpolate SL:
    double sl=SL[binE]+(SL[binE+1]-SL[binE])*(E-En[binE])/(En[binE+1]-En[binE]);
    return sl;
  }
  else
  {
    //////std::cout<<"***ERROR: E is out of energy range!"<<std::endl;
    return -100.0;
  }  
}

void FillHistograms()
{


//BEGIN OF PROGRAM
//*******************************************************************************************  

  COUNTCHANNEL1=0;
  COUNTCHANNEL1INBIN=0;
  COUNTCHANNEL1OUTOFBIN=0;

  COUNTCHANNEL2=0;
  COUNTCHANNEL2INBIN=0;
  COUNTCHANNEL2OUTOFBIN=0;

  //make differentiated partial sums in dps array:
  const int PSINDEX=BinarySearch(En, N1, Ee);

  std::cout<<"En="<<Ee<<" PSINDEX="<<PSINDEX<<" "<<En[PSINDEX]
           <<" "<<En[PSINDEX+1]<<std::endl;

  //for(int i=0; i<N1; ++i) std::cout<<En[i]

  /////sleep(3);


  std::cout<<"After FillDifferentiatedPartialSums():"<<std::endl;
  
  /*
  std::cout<<"Check PartialSums:"<<std::endl;
  for(int i=0; i<NpointsInCosCM; ++i) std::cout<<PS[i]<<" ";
  std::cout<<std::endl;
    
  //fill exponent dpsexp:
  FillExponent(SL[PSindex], PR[PSindex], ind);
  std::cout<<"Print Ox:"<<std::endl;
  for(int i=0; i<NpointsInCosCM-1; ++i) std::cout<<Ox[i]<<" ";
  std::cout<<std::endl;


  std::cout<<"Print PartialSums: size="<<std::endl;
  for(int i=0; i<NpointsInCosCM; ++i) std::cout<<PS[i]<<" ";
  std::cout<<std::endl;
  */

  
//-----------------------------------------------------------------------------------------------//
//Fill in the histogram:
//-----------------------------------------------------------------------------------------------//
  //prepare the border array:
  double array[CHANNELS+1]{0.0};
  array[0]=0.0;
  array[1]=array[0]+1.0-GetPR(Ee);
  array[2]=array[1]+GetPR(Ee);
  for(int i=0; i<CHANNELS+1; ++i) array[i]=array[i]/array[2];
  std::cout<<"Check border array: CHANNELS="<<CHANNELS<<std::endl;
  for(int i=0; i<CHANNELS+1; ++i) std::cout<<array[i]<<" ";
  std::cout<<std::endl;
  
  //exit(0);
  
  //


  std::cout<<"CHECK BEFORE FOR LOOP IN FillHistograms():"<<std::endl;
  
  int count=0;
  RNDGenerator rand(1);
#ifdef OPENACC
//#pragma acc parallel loop num_gangs(1) vector_length(1) present(PS,HistogramNiCosCM) copy(array,rand)
#pragma acc serial loop present(HistogramNiCosCM) copy(array,rand,COUNTCHANNEL1,COUNTCHANNEL1INBIN,COUNTCHANNEL1OUTOFBIN,COUNTCHANNEL2,COUNTCHANNEL2INBIN,COUNTCHANNEL2OUTOFBIN,OutOfHistogramRange)
#endif
  for(int i=0; i<NPARTICLES; ++i)
  {
    const double P=GenerateSubCanonical<double>(rand);
    /*
    int channel=0;
    for(int j=0; j<CHANNELS; ++j)
    {
      if(array[j]<=R && R<=array[j+1]) channel=j+1;
    }
    */
    
    int bin=-1;
    double r=0.0;

    //std::cout<<"i="<<i<<std::endl;
    
    double costhetacm=-1.0;
    const double pr=GetPR(Ee);
    //const double sl=GetSL(Ee);
    //std::cout<<"pr="<<pr<<" sl="<<sl<<std::endl;
    //sleep(1);
    if(P < pr)
    {//if(channel==2)
      // use exponent
      r=GenerateSubCanonical<double>(rand);
      if(0.0<=r && r<=1.0)
      {
        const double sl=GetSL(Ee);
        costhetacm=1. + sl * log(1. - ( 1. - exp(-2 / sl) ) * r );
        bin=(costhetacm+1.0)/deltacoscm;

        //std::cout<<"1 bin="<<bin<<" costhetacm="<<costhetacm<<" costhetacm+1.0="<<costhetacm+1.0
        //         <<" deltacoscm="<<deltacoscm<<std::endl;
        
      #ifdef OPENACC
      #pragma acc atomic update
      #endif
        ++HistogramNiCosCM[bin];
      #ifdef OPENACC
      #pragma acc atomic update
      #endif
        ++COUNTCHANNEL2INBIN;
      #ifdef OPENACC
      #pragma acc atomic update
      #endif
        ++COUNTCHANNEL2;
      }
    }
    else
    {//if(channel==1)
    #ifdef OPENACC
    #pragma acc atomic update
    #endif
      ++COUNTCHANNEL1;
      r=GenerateSubCanonical<double>(rand);
//\\//!!!!!!!!!!!!!!!!//
      r=sqrt(r);
//\\//!!!!!!!!!!!!!!!!//            
      //thetams=std::sqrt(-D1(ZZ/t3::units::um)*std::log(r));
      if(0.0<=r && r<=1.0)
      {
        ///\\\///bin=BinarySearch(PS[PSindex], N2, r);

        ////////bin=GetBinNumberForT2_DATA_DD_FROM_PARTIAL_SUMS(10.0/*MeV*/, r);
        const double CosThetaCM=GetCosThetaCMFromDDInelasticPartialSums(Ee, r);
        // std::cout<<"Exit GetBinNumberForT2_DATA_DD_FROM_PARTIAL_SUMS():"<<std::endl;
        bin=(CosThetaCM+1.0)/deltacoscm;

        //std::cout<<"2 bin="<<bin<<" CosThetaCM="<<CosThetaCM<<" CosThetaCM+1.0="<<CosThetaCM+1.0
        //         <<" deltacoscm="<<deltacoscm<<std::endl;
        
        //std::cout<<"bin="<<bin<<" CosThetaCM="<<CosThetaCM<<" r="<<r<<std::endl;
        //sleep(1);
        
        
        //std::cout<<"channel "<<1<<std::endl;
      #ifdef OPENACC
      #pragma acc atomic update
      #endif
        ++HistogramNiCosCM[bin];
      #ifdef OPENACC
      #pragma acc atomic update
      #endif
        ++COUNTCHANNEL1INBIN;
        
      }
      else
      {
      #ifdef OPENACC
      #pragma acc atomic update
      #endif
        ++COUNTCHANNEL1OUTOFBIN;
      }
    }

    /*
    else
    {
    #ifdef OPENACC
    #pragma acc atomic update
    #endif
      ++OutOfHistogramRange;

      printf("No channel!");
    }
    */

    //std::cout<<"i="<<i<<" channel="<<channel<<std::endl;
    //std::cout<<"costhetacm="<<costhetacm<<std::endl;
    //sleep(1);

        

      /*
      if(bin==-1) bin=1;
      
      std::cout<<"i="<<i<<std::endl;
      std::cout<<"array: "<<array[0]<<" "<<array[1]<<" "<<array[2]<<std::endl;
      std::cout<<"R="<<R<<" channel="<<channel<<std::endl;
      std::cout<<"Print ps0:"<<std::endl;
      for(int j=0; j<NpointsInCosCM; ++j) std::cout<<ps0.at(j)<<" ";
      std::cout<<std::endl;
      std::cout<<"bin="<<bin<<" r="<<r<<std::endl;
      std::cout<<"ps0.at(0)="<<ps0.at(0)<<" ps0.at(NpointsInCosCM-1)="<<ps0.at(NpointsInCosCM-1)<<std::endl;
      if(bin>0 && bin<NpointsInCosCM-1)
      std::cout<<"ps0.at("<<bin-1<<")="<<ps0.at(bin-1)
               <<"ps0.at("<<bin<<")="<<ps0.at(bin)
               <<"ps0.at("<<bin+1<<")="<<ps0.at(bin+1)
               <<std::endl;
      //sleep(1);
      */
  #ifdef OPENACC
  #pragma acc atomic update
  #endif
    ++count;
  }

  std::cout<<"deltacoscm="<<deltacoscm<<std::endl;
  std::cout<<"COUNT="<<count<<std::endl;

  std::cout<<"COUNTCHANNEL1="<<COUNTCHANNEL1<<" COUNTCHANNEL1INBIN="<<COUNTCHANNEL1INBIN<<" COUNTCHANNEL1OUTOFBIN="<<COUNTCHANNEL1OUTOFBIN
           <<" COUNTCHANNEL2="<<COUNTCHANNEL2<<" COUNTCHANNEL2INBIN="<<COUNTCHANNEL2INBIN<<" COUNTCHANNEL2OUTOFBIN="<<COUNTCHANNEL2OUTOFBIN
           <<std::endl;
  std::cout<<"Check balance:"<<std::endl;
  std::cout<<"COUNTCHANNEL1="<<COUNTCHANNEL1<<" COUNTCHANNEL1INBIN+COUNTCHANNEL2OUTOFBIN="<<COUNTCHANNEL1INBIN+COUNTCHANNEL2OUTOFBIN
           <<" COUNTCHANNEL2="<<COUNTCHANNEL2<<" COUNTCHANNEL2INBIN+COUNTCHANNEL2OUTOFBIN="<<COUNTCHANNEL2INBIN+COUNTCHANNEL2OUTOFBIN
           <<std::endl;
  std::cout<<"COUNTCHANNEL1+COUNTCHANNEL2="<<COUNTCHANNEL1+COUNTCHANNEL2<<std::endl;

  std::cout<<"Print Histogram:"<<std::endl;
  int TotalParticles=0;
#ifdef OPENACC
#pragma acc parallel loop gang vector reduction(+:TotalParticles)
#endif
  for(int h=0; h<NpointsInCosCM-1; ++h)
  {
    //#pragma acc atomic update
    TotalParticles+=HistogramNiCosCM[h];
  }

#ifdef OPENACC
#pragma acc update host(HistogramNiCosCM)
#endif
  for(int h=0; h<NpointsInCosCM-1; ++h)
  {
    std::cout<<HistogramNiCosCM[h]<<" ";
  }
  
  std::cout<<std::endl;
  std::cout<<"There are "<<TotalParticles<<" particles in the histogram"<<std::endl;
  std::cout<<"There are "<<OutOfHistogramRange<<" particles ouf of histogram range"<<std::endl;


  std::cout<<"pr="<<GetPR(Ee)<<std::endl;
  std::cout<<"sl="<<GetSL(Ee)<<std::endl;
  std::cout<<"EXPONENT_CONTRIBUTION="<<GetPR(Ee)<<std::endl;
  std::cout<<"PS_CONTRIBUTION="<<1.0-GetPR(Ee)<<std::endl;
  std::cout<<"EXPONENT_SLOPE="<<GetSL(Ee)<<std::endl;
  
  

  /*
  //Check of coscm:
  std::cout<<"Print coscm:"<<std::endl;
  for(int i=0; i<T2NSGangular_RWnode::_num_point; ++i) std::cout<<coscm.at(i)<<" ";
  std::cout<<std::endl;
  std::cout<<"Check coscm:"<<std::endl;
  const double d=2.0/128;
  for(int i=0; i<T2NSGangular_RWnode::_num_point; ++i) std::cout<<-1.0+d+i*d<<" ";
  std::cout<<std::endl;
  exit(0);
  */

} 







//END OF WHAT I WILL USE IN THE NEW PROGRAM


//////////////////////////////////////////////////////////////////////////////
// Program main
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
#ifdef OPENACC
  #ifdef CUDA
    std::cout<<"OPENACC IS DEFINED, CUDA IS DEFINED"<<std::endl;
  #else
    std::cout<<"OPENACC IS DEFINED, CUDA IS NOT DEFINED"<<std::endl;
  #endif
#else
  #ifdef CUDA
    std::cout<<"OPENACC IS NOT DEFINED, CUDA IS DEFINED"<<std::endl;
  #else
    std::cout<<"OPENACC IS NOT DEFINED, CUDA IS NOT DEFINED"<<std::endl;
  #endif
#endif

  //sleep(1);

  auto begin=std::chrono::steady_clock::now();


  //FillT2_DATA_DD_Database();
  //SaveDDPartialSumsToAFile();

//****************************************************************************************//
//Begin of read the partial sums form T2_DATA:
//****************************************************************************************//  
  T2NSGangular_RW rw;  
  //target deuteron Z=1
  static constexpr auto tgZ = 1;
  //target deuteron A=2
  static constexpr auto tgA = 2;
  //reaction product=n+3He
  //this is a neutron from a reaction product:
  static constexpr auto sPDG = 2112;
  //1000*Z+A//inc deuteron
  static constexpr auto incZA = 1002;
  //MT index of reaction in ENDF
  static constexpr auto rid = "50";
  //load the data from the binary file:
  rw.load_binary(tgZ, tgA, rid, sPDG, incZA);
  //print the data records (e, partical sums):
  //std::cout << rw << std::endl;
  std::cout<<"rw size="<<rw.size()<<std::endl;
  std::cout<<"SEE:"<<std::endl;
  std::cout<<rw<<std::endl;


  std::cout<<"Check 0.1 MeV:"<<std::endl;
  std::cout<<rw.at(0).at(11)<<std::endl;

  T2NSGangular_RWnode node=rw.at(0).interpolate(Ee);

  std::cout<<"node:"<<std::endl;
  std::cout<<node<<std::endl;

//******************************************************//
//Fill the arrays:
  for(int i=0; i<N1; ++i)
  {
    En[i]=rw.at(0).at(i).Get_E();
    PR[i]=rw.at(0).at(i).Get_pr();
    SL[i]=rw.at(0).at(i).Get_sl();
    PS[i][0]=0.0;
    PS[i][N2-1]=1.0;
    std::array<G4double, 127> psarray=rw.at(0).at(i).Get_V();
    for(int j=0; j<N2-2; ++j)
    {
      PS[i][j+1]=psarray.at(j);
    }
  }
//End of fill the arrays:
//******************************************************//
  
//******************************************************//   
//End if read the partial sums from T2_DATA.
//****************************************************************************************//

  
//Old read from the file:
  //\\//ReadDDPartialSumsFromAFile();

  std::cout<<"Check data: 2:"<<std::endl;
  std::cout<<"Energies:"<<std::endl;
  for(int i=0; i<N1; ++i) std::cout<<En[i]<<" ";
  std::cout<<std::endl;
  std::cout<<"PR:"<<std::endl;
  for(int i=0; i<N1; ++i) std::cout<<PR[i]<<" ";
  std::cout<<std::endl;
  std::cout<<"SL:"<<std::endl;
  for(int i=0; i<N1; ++i) std::cout<<SL[i]<<" ";
  std::cout<<std::endl;
  std::cout<<"Partial Sums:"<<std::endl;
  for(int i=0; i<N1; ++i)
  {
    std::cout<<"E="<<En[i]<<" pr="<<PR[i]<<" sl="<<SL[i]<<" ";
    for(int j=0; j<N2; ++j) std::cout<<PS[i][j]<<" ";
    std::cout<<std::endl;
  }
  std::cout<<std::endl;


  T2NSGangular_RWnode curnode=rw.at(0).interpolate(Ee);
  std::array<G4double, 127> psarray=curnode.Get_V();
  double PsArray[129];
  PsArray[0]=0.0;
  PsArray[128]=1.0;
  for(int k=0; k<127; ++k) PsArray[k+1]=psarray.at(k);

  MakeEquidistantBinsCosThetaCM();
  FillOxAxisCosCM();
  //output bin number/Ox axis points to "ne_thetacoscm.dat" file.
  Histogram_OxAxisCosCM();

  std::cout<<"coscm:"<<std::endl;
  for(int i=0; i<129; ++i) std::cout<<coscm[i]<<" ";
  std::cout<<std::endl;
  std::cout<<"Ox:"<<std::endl;
  for(int i=0; i<128; ++i) std::cout<<Ox[i]<<" ";
  std::cout<<std::endl;

  std::cout<<"PsArray:"<<" E="<<Ee<<" pr="<<curnode.Get_pr()<<" sl="<<curnode.Get_sl()<<std::endl;
  for(int i=0; i<128; ++i) std::cout<<PsArray[i]<<" ";
  std::cout<<std::endl;
  
  FillDifferentiatedPartialSums(PsArray, 1.0-curnode.Get_pr());
  Histogram_DifferentiatedPartialSums();
  //Fill dpsexp:
  FillExponent(curnode.Get_sl(), curnode.Get_pr());
  //output dpsexp to "ne_thetexp.dat":
  Histogram_Exponent();
  std::cout<<"E="<<Ee<<" pr="<<curnode.Get_pr()<<" sl="<<curnode.Get_sl()<<std::endl;
  std::cout<<"curnode:"<<std::endl;
  std::cout<<curnode<<std::endl;

  //sleep(1);
  //exit(0);


#ifdef OPENACC
#pragma acc data copy(coscm,Ox) copy(dpsmain,dpsexp,HistogramNiCosCM) copyin(coscm,PS,En,PR,SL)
{
#endif
  
  std::cout<<"MakeEquidistantBinsCosThetaCM():"<<std::endl;
  MakeEquidistantBinsCosThetaCM();
#ifdef OPENACC
#pragma acc update host(coscm)
#endif
  std::cout<<"Check coscm:"<<std::endl;
  for(int i=0; i<NpointsInCosCM; ++i) std::cout<<coscm[i]<<" ";
  std::cout<<std::endl;

  std::cout<<"CHECK 1"<<std::endl;

  FillOxAxisCosCM();

  std::cout<<"Check coscm and Ox:"<<std::endl;
  std::cout<<"coscm:"<<std::endl;
  for(int i=0; i<NpointsInCosCM; ++i) std::cout<<coscm[i]<<" ";

  std::cout<<std::endl;
  std::cout<<"Ox:"<<std::endl;
  for(int i=0; i<NpointsInCosCM-1; ++i) std::cout<<Ox[i]<<" ";
  std::cout<<std::endl;

  //sleep(10);

//*
//Here will check energy extrapolation:

  /*
  const double E=8.2;
  RNDGenerator gen(1);
  const double Rand=GenerateSubCanonical<double>(gen);
  double costhetacm=GetCosThetaCMFromDDInelasticPartialSums(E, Rand);
  //std::cout.precision(8);
  std::cout<<"costhetacm="<<costhetacm<<std::endl;
  */

  //exit(0);
  
//End of where will check energy extrapolation:  
//*/
  


  std::cout<<"CHECK 2"<<std::endl;
  
  
  //output bin number/Ox axis points to "ne_thetacoscm.dat" file.
  //Histogram_OxAxisCosCM();

  std::cout<<"CHECK 3"<<std::endl;

 
  //output bin number/differentiated partial sums to "ne_thetadps.dat".
  //Histogram_DifferentiatedPartialSums();

  std::cout<<"CHECK 4"<<std::endl;
  
  //fill NE histograms for NE energies:
  FillHistograms();

  std::cout<<"CHECK 5"<<std::endl;

  //exit(0);

  FillHistogramNi();

  
  Histogram_DifferentiatedPartialSums();

  //Fill dpsexp:
  //FillExponent(GetSL(Ee), GetPR(Ee));

  //std::cout<<"pr="<<GetPR(Ee)<<" sl="<<GetSL(Ee)<<std::endl;
  
  //output dpsexp to "ne_thetexp.dat":
  //Histogram_Exponent();

#ifdef OPENACC
 }
#endif





  

  


//END OF PROGRAM

  auto end=std::chrono::steady_clock::now();
	auto elapsed_ms=std::chrono::duration_cast<std::chrono::milliseconds>(end-begin);
	std::cout<<"time="<<elapsed_ms.count()<<" ms"<<std::endl;
  std::cout<<"Nbin="<<Nbin<<std::endl;
}
