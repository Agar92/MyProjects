#include "T2NSGangular_RW.hh"
#include <iostream>
#include "unistd.h"
#include <iterator>
#include <fstream>
#include <iomanip>


//the step is 2/128, but there are 127 points in Ox (cos(theta_CM)) array:
//from -1+2/128 to 1-2/128, there are 127 points:
constexpr int Nbin=128;
constexpr int NpointsInCosCM=Nbin-1;
std::array<double, NpointsInCosCM> coscm;

//there are totally 128 bins.
//the length from -1 to 1 is 2
//points -1 and 1 are not included in Ox (cos(theta_CM)),
//So, the 1st point of Ox is -1+2/128,
//the last point of Ox is 1-2/128.
void MakeEquidistantBinsCosThetaCM()
{
  const double delta=static_cast<double>(2.0)/Nbin;
  const double coscm1=-1.0+delta;
  for(int i=0; i<NpointsInCosCM; ++i) coscm.at(i)=coscm1+i*delta;
}

int main()
{
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
    //sleep(3);
    //std::cout<<rw.Get_vector_of_T2NSGangular_RWrecords().at(0)<<std::endl;

    /*
    for(T2NSGangular_RW::const_iterator it = rw.begin(); it != rw.end(); ++it)
    {
      std::cout<<*it<<std::endl;
    }
    */

    //std::cout<<"SIZE="<<rw.size()<<" SIZE OF THE RECORD="<<rw.at(0).at(11).Get_E()<<std::endl;
    //std::cout<<rw.at(0).at(11)<<std::endl;

    std::cout<<rw.at(0)<<std::endl;

    std::array<G4double, T2NSGangular_RWnode::_num_point> arr= rw.at(0).at(11).Get_V();
    
    MakeEquidistantBinsCosThetaCM();


    std::ofstream out1;
    out1.open("output1.dat");
    for(int i=0; i<(int)T2NSGangular_RWnode::_num_point; ++i)
    {
      out1<<std::setw(10)<<coscm.at(i)<<"   "<<arr.at(i)<<std::endl;
    }
    out1<<std::endl;
    out1.close();

    std::ofstream out2;
    out2.open("output2.dat");
    const double delta=static_cast<double>(2.0)/Nbin;
    for(int i=0; i<(int)T2NSGangular_RWnode::_num_point-1; ++i)
    {
      const double v=(arr.at(i+1)-arr.at(i))/(delta);
      out2<<std::setw(10)<<coscm.at(i)<<"   "<<v<<std::endl;
    }
    out2<<std::endl;
    out2.close();
    

    /*
    std::cout<<"Print coscm:"<<std::endl;
    for(int i=0; i<NpointsInCosCM; ++i) std::cout<<coscm.at(i)<<" ";
    std::cout<<std::endl;
    */


    /*
    std::array<double, T2NSGangular_RWnode::_num_point> ps0=rw.at(0).at(66).Get_V();

    for(int i=0; i<T2NSGangular_RWnode::_num_point; ++i) std::cout<<ps0.at(i)<<" ";
    std::cout<<std::endl;
    */
    
    
    
    
    /*
    rw.MakeEquidistantBinsCosThetaCM();
    int Nbin=rw.GetNpointsInCosCM();
    std::array<double, T2NSGangular_RW::NpointsInCosCM> CosCM;
    CosCM=rw.GetEquidistantBinsCosThetaCM();
    std::cout<<"Print cos(theta_CM):"<<" Number Of Elements="
             <<CosCM.size()<<std::endl;
    for(int i=0; i<T2NSGangular_RW::NpointsInCosCM; ++i)
      std::cout<<CosCM.at(i)<<" ";
    std::cout<<std::endl;
    */
    
    
    return 0;
}
