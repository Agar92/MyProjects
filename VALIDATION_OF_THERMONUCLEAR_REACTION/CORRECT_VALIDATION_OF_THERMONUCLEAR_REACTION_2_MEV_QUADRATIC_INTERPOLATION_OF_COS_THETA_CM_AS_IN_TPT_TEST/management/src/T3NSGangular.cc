#include "T3NSGangular.hh"
#include <stdexcept>
//\\//#include "Randomize.hh"
#include "T3NSGangular_record.hh"
#include <map>
#include <mutex>
#include <array>
#include <string>

#include "unistd.h"

//#define debug

namespace T3NSGangular
{
  class redChanID
  {
    // FIXME private
public:
    redChanID(T3int sPDG, T3int rid, T3int level): secPDG(sPDG), rID(rid), levN(level){}
    T3int secPDG;
    T3int rID;
    T3int levN;
  };

  namespace
  {
    class redMidar : public std::map<const redChanID, T3NSGangular_record*>
    {
    public:
      ~redMidar()
      {
        // TODO the same as in energy-angular, unify
#ifdef pdebug
        T3cout << "Angular data destructor is called" << T3endl;
#endif
        for(auto& it : *this)
        {
          if (it.second)
          {
            delete it.second;
            it.second = nullptr;
          }
        }
      }
    };

    static const size_t totN = 256;
    static const size_t totZ = 128;
    static const size_t totZN = totZ*totN;
    T3int ZNindex(T3int Z, T3int N) {return totN*Z + N;}
    static std::array<redMidar, totZN> nodes;
//    static std::array<std::once_flag, totZN> initOnce;
//    static std::array<G4bool, totZN> hasInitialized;
    static std::array<std::mutex, totZN> nodeMutexes;
  }


  //\\//G4double RandomizeCost(G4int tgZ, G4int tgA, G4int sPDG, G4int rID,
  //\\//                       G4int level, G4double Einc, G4int incZA)
  T3double RandomizeCost(RNDGenerator & generator, T3int tgZ, T3int tgA, T3int sPDG, T3int rID,
                         T3int level, T3double Einc, T3int incZA)
  {

    //\\//std::cout<<"Begin T2NSGangular::RandomizeCost():"<<std::endl;
    
    // TODO profile, maybe load all levels for a reaction at once
    // FIXED multiple adding for the same chid: at() throws, insert() returns false
//    static std::mutex totalmutex;
    auto find = [tgZ, tgA, sPDG, rID, level, incZA](){
//      const chanID chid(tgZ, tgA, sPDG, rID, level);
      const redChanID rchid(sPDG, rID, level);
      auto const znindex = ZNindex(tgZ, tgA - tgZ);
      auto& node = nodes.at(znindex);
      auto& nodeMutex = nodeMutexes.at(znindex);
//      auto& hIn = hasInitialized.at(znindex);
      // FIXME don't lock for reads
      //
      // TODO implement readers-writers approach
      // DONE independent load for each element
      std::lock_guard<std::mutex> nodeLock(nodeMutex);
      auto r = node.insert(std::make_pair(rchid, nullptr));

      //\\//std::cout<<"r.second="<<r.second<<std::endl;
      
      if(r.second)
      {
        T3NSGangular_RW rw;
        T3NSGangular_RWrecord rwrec;
        
        // TODO process inside rw
        T3String rid;
        // TODO other reactions
        if(rID == 4) rid = "inel";
        else if(rID == 2) rid = "elastic";
        else rid = std::to_string(rID);

        rw.load_binary(tgZ, tgA, rid, sPDG, incZA);

        size_t ind = 0;
        for(; ind < rw.size(); ++ind) if (level == rw.at(ind).Get_levN()) break;
        if(ind != rw.size()) rwrec = rw.at(ind);

////\\\\////Here square the values of partial sums:
        for(T3NSGangular_RWrecord::iterator it = rwrec.begin(); it != rwrec.end(); ++it)
        {
          for(int i=0; i<T3NSGangular_RWnode::_num_point; ++i)
            it->Set_V(i+1, it->Get_V(i+1)*it->Get_V(i+1));
        }
        ////std::cout<<"Check rwrec:"<<std::endl;
        ////std::cout<<rwrec.at(39)<<std::endl;        
////\\\\////End of here square the values of partial sums.                
        
        r.first->second = rwrec.size() ? new T3NSGangular_record(rwrec) : NULL;

#ifdef debug
        T3cout << "RandomizeCost: adding Z=" << tgZ << " A=" << tgA << " rID=" << rID
               << " level=" << level << " sPDG=" << sPDG << " rwrecsize="
               << rwrec.size() << T3endl;
#endif
      }
#ifdef debug
      else T3cout << "RandomizeCost: found Z=" << tgZ << " A=" << tgA << " rID=" << rID
                  << " level=" << level << " sPDG=" << sPDG << T3endl;
#endif
      return r;
    };
    auto r = find();
    T3NSGangular_record* rec = r.first->second;

    //\\//std::cout<<"End T2NSGangular::RandomizeCost():"<<std::endl;
    
    //\\//return rec ? rec->RandomizeCost(Einc) : 1 - 2*G4UniformRand();;
//\\//***********************************************//
    return rec ? rec->RandomizeCost(generator, Einc) : 1 - 2*t3::GenerateSubCanonical<FloatingType>(generator);;
//\\//***********************************************//        
  }

  inline T3bool operator<(const T3NSGangular::redChanID& lhs,
                          const T3NSGangular::redChanID& rhs)
  {
    if (lhs.secPDG < rhs.secPDG) return true;
    else if (lhs.rID < rhs.rID) return true;
    else if (lhs.levN < rhs.levN) return true;
    else return false;
  }

  inline T3bool operator==(const T3NSGangular::redChanID& lhs,
                           const T3NSGangular::redChanID& rhs)
  {
    return (lhs.secPDG == rhs.secPDG && lhs.rID == rhs.rID &&
            lhs.levN == rhs.levN);
  }
}
