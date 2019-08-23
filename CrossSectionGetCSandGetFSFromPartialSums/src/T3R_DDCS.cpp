// ---------------------------------------------------------------------------
//  In this file there is a d-d elastic scattering approximation
//  from /users/ALPHA_ELAST/FINAL_CORRECT_dd.
// ---------------------------------------------------------------------------

// #define debug
// #define pdebug

#include "T3R_DDCS.h"
#include <cmath>
#include <map>

//********************************************************************************************************************//
//This is a file for writing and reading |t| and partial sums for each energy for d-d elastic scattering approximation//
//********************************************************************************************************************//

//const double md=1.8756*GeV;
//const double md2=md*md;//GeV^2

namespace t3 {

T3R_DDCS::T3R_DDCS():lnE{0.0}, lncos{0.0}, Energy{0.0},
                     CS{0.0}, Sk{0.0}, Fk{0.0}, ecs{}, tps{}, trw{}, TailCS{0.0}, RutherfordTailCS{0.0}
{}

void T3R_DDCS::Fill()
{
//prepare E axis for E in partial sums db (Bin2):
  const double lnEmin=std::log(Emin);
  const double lnEmax=std::log(Emax);
  deltalnE=(lnEmax-lnEmin)/Bin2;
  for(int b=0; b<Bin2+1; ++b) lnE.at(b)=lnEmin+b*deltalnE;
  //for min value of E (Ede - displacement energy) calculate min value of |t| for d-d:
  const double tmin=2*md*Edisplace_deuteron;//|t|min=2*M*Ede, Ede for d-d = 10 eV.
  //for max value of E calculate max value of tmax=4pcm^2 for d-d:
  const double mpls2=Emax*(2*md+Emax);//pls^2 of the inc deuteron.
  //for d-d rat=rat2=1.0
  const double mpcm2=mpls2/(2.0+2.*std::sqrt(1.0+mpls2/md2));
  //INTEGRATE dsigma/d|t| from |t|min to (|t|max-|t|min), because for scattering of equal particles d-d
  //dsigma/d|t| is symmetric at |t|=2*pcm^2 and multiple scattering cross section is int_0_|t|min(dsigma/d|t|)+
  //int_(|t|max-|t|min)_|t|max(dsigma/d|t|)
  const double tmaxEmax=4*mpcm2;
  const double tmEmax=2*mpcm2;
  //make logarithmic |t|min and |t|m:
  const double lnTmin=std::log(tmin);
  const double lnTm=std::log(tmEmax);
//we decided that the minimum of 1-cos(theta_cm) axis (it must be > 0, because the scale
//is logarithmic) must be tmin/2*pcm^2(Emax).
  MinLnCosGlobal=std::log(tmin/tmEmax);
  //prepare 1-cos(theta_cm) axis:
  deltalncos=(std::log(1.0)-MinLnCosGlobal)/Bin3;
  //we take the values of 1-cos(theta_cm) at each border of the bin from right to left (we sum partial sums from right to left)
  //If Nbin=100, number of Ei will be 101, (if there are Nbin bins, the number of walls is Nbin+1).
  for(int b=0; b<Bin3+1; ++b) lncos.at(b)=MinLnCosGlobal+b*deltalncos;
  //here we divide the 1-cos(theta_cm) intevals by SubBin3 subbins:
  std::array<std::array<double, SubBin3+1>, Bin3> SubBinsCS{0.0};//ln(1-cos(theta_cm)) subbins array
  //here we divide the 1-cos(theta_cm) intervals by SubBin3 subbins:
  const double deltalncosSubBins=deltalncos/SubBin3;
  //\\//for(int b=1; b<Bin3+1; ++b)
  for(int b=0; b<Bin3; ++b)
  {
    //we fill sub bins from right to left, because we make partial sums from right to left
    for(int sb=0; sb<SubBin3+1; ++sb)
        //\\//SubBinsCS.at(b-1).at(sb)=lncos.at(b-1)-sb*deltalncosSubBins;
        SubBinsCS.at(b).at(sb)=lncos.at(b)+sb*deltalncosSubBins;
  }
  std::cout<<"MinLnCosGlobal="<<std::exp(MinLnCosGlobal)<<" Check 1-cos(theta_cm):   ";
  for(int b=0; b<Bin3+1; ++b) std::cout<<std::exp(lncos.at(b))<<" ";
  std::cout<<std::endl;
  std::cout<<"Check SubBinsCS:"<<std::endl;
  for(int b=0; b<Bin3; ++b)
  {
    std::cout<<std::exp(lncos.at(b))<<":   ";
    std::array<double, SubBin3+1> sbarray=SubBinsCS.at(b);
    for(int sb=0; sb<SubBin3+1; ++sb) std::cout<<std::exp(sbarray.at(sb))<<" ";
    std::cout<<std::endl;
  }

  //iterate through all bins in E[0:Bin2]:
  for(int i=0; i<Bin2+1; ++i)
  {
    const double Tls=std::exp(lnE.at(i));//kinetic energy of the particle E[i] in inner units (MeV).
    //calculate tmax:
    const double pls2 = Tls*(2*md+Tls);
    const double pcm2 = pls2/(2.+2.*std::sqrt(1.+pls2/md2));//rat=md/md=1. rat2=1.
    //FOR d-d |t|max>Ede, because Tls (LS kin energy of the inc particle) should be >=2*Ede.
    const double tm =   2*pcm2;
    //real 1-cos(theta_cm_min) for this energy = |t|min/(2*pcm^2)
    const double MinLnCos=std::log(tmin/tm);

    std::cout<<"Tls("<<i<<")="<<Tls<<" pls2="<<pls2/GeV/GeV<<" pcm2="<<pcm2/GeV/GeV<<" tm="<<tm/GeV/GeV
             <<" MinLnCos="<<MinLnCos<<std::endl;

//fill partial sums.
    double psumb=0.0;//here we will sum sum_0_Bin(dsigma/d|t|*d|t|).   
    double psumb_Rutherford=0.0;//here we will sum sum_0_Bin(dsigma/d|t|*d|t|) for Rutherford.
    Fk.at(i).at(0)=0.0;
    for(int b=1; b<Bin3+1; ++b)
    {
//WE SUM  PARTIAL SUMS FROM RIGHT TO LEFT, BECAUSE RUTHERFORD FUNCTION ~ 1/|t|^2,
//AND TO SAVE (MORE ACCURATE CALCULATIONS) SMALL VALUES FROM THE BINS NEAR |t|=2*pcm^2,
//WE SUM FROM RIGHT TO LEFT (from smaller values of dsigma/d|t| to bigger).
      std::array<double, SubBin3+1> subbinscs=SubBinsCS.at(Bin3-b);
      std::cout<<"Bin #"<<b<<" "<<std::exp(lncos.at(Bin3-b))<<" "<<std::exp(lncos.at(Bin3-b+1))
               <<" tmin="<<tmin/GeV/GeV<<" pls2="<<pls2/GeV/GeV<<" pcm2="<<pcm2/GeV/GeV<<" tm="<<tm/GeV/GeV<<std::endl;
      std::cout<<"Check subbins:   ";
      for(int sb=0; sb<SubBin3+1; ++sb) std::cout<<std::exp(subbinscs.at(SubBin3-sb))<<" ";
      std::cout<<std::endl;
      std::cout<<"Check |t| subbins:   ";
      for(int sb=0; sb<SubBin3+1; ++sb) std::cout<<tm*std::exp(subbinscs.at(SubBin3-sb))/GeV/GeV<<" ";
      std::cout<<" |t|min="<<tmin/GeV/GeV<<std::endl;
      for(int sb=0; sb<SubBin3; ++sb)
      {
        int path=0;
        int stop_bin=-1;
        const double tl=tm*std::exp(subbinscs.at(SubBin3-sb-1));
        const double tr=tm*std::exp(subbinscs.at(SubBin3-sb));
        double dleft=tl*GetdSigmadt(std::exp(lnE.at(i)), tl);
        double dright=tr*GetdSigmadt(std::exp(lnE.at(i)), tr);

        double dleftR=tl*GetdSigmadt_DD_Rutherford(std::exp(lnE.at(i)), tl);
        double drightR=tr*GetdSigmadt_DD_Rutherford(std::exp(lnE.at(i)), tr);

        double dlncos=deltalncosSubBins;//=dx
        const double DeltaLnCos=MinLnCos-subbinscs.at(SubBin3-sb);
        if(DeltaLnCos>=0.0)
        {
          dright=0.0;
          dleft=0.0;

          drightR=0.0;
          dleftR=0.0;

          path=2;
        }
        else if(MinLnCos>subbinscs.at(SubBin3-sb-1) && MinLnCos<subbinscs.at(SubBin3-sb))
        {
          const double tleft=tmin;
          dleft=tleft*GetdSigmadt(std::exp(lnE.at(i)), tleft);
          dleftR=tleft*GetdSigmadt_DD_Rutherford(std::exp(lnE.at(i)), tleft);
          dlncos=subbinscs.at(SubBin3-sb)-MinLnCos;
          stop_bin=SubBin3-sb-1;
          path=3;
        }
        const double dsigmai=(dleft+dright)*dlncos/2;
        psumb+=dsigmai;
        const double dsigmai_Rutherford=(dleftR+drightR)*dlncos/2;
        psumb_Rutherford+=dsigmai_Rutherford;
        if(5.0e-1*tm<tr) TailCS.at(i)+=dsigmai;
        if(5.0e-1*tm<tr) RutherfordTailCS.at(i)+=dsigmai_Rutherford;
        std::cout<<"tl="<<tl/GeV/GeV<<" tr="<<tr/GeV/GeV<<" 1.0e-2*GeV^2"
                 <<" TailCS="<<TailCS.at(i)/mbarn<<" RuTailCS="<<RutherfordTailCS.at(i)<<std::endl;

        std::cout<<"SubBin #"<<sb<<" "<<subbinscs.at(SubBin3-sb-1)<<" "<<subbinscs.at(SubBin3-sb)
                 <<" "<<std::exp(subbinscs.at(SubBin3-sb-1))<<" "<<std::exp(subbinscs.at(SubBin3-sb))
                 <<" tl="<<tl/GeV/GeV<<" tr="<<tr/GeV/GeV<<" ds/dt_l="<<GetdSigmadt(std::exp(lnE.at(i)),tl)/mbarn*GeV*GeV
                 <<" ds/dt_r="<<GetdSigmadt(std::exp(lnE.at(i)),tr)/mbarn*GeV*GeV
                 <<" dleft="<<dleft/mbarn<<" dright="<<dright/mbarn<<" dlncos="<<dlncos
                 <<" path="<<path<<" dsigmai="<<dsigmai/mbarn<<" ps="<<psumb/mbarn<<std::endl;
      }
      //BECAUSE THE PARTICLES ARE EQUAL (WE CAN NOT TELL THE INCIDENT PARTICLE FROM THE SCATTERED)
      //WE GET THE RANDOM VALUE OF |t| from |t|_{min} to |t|_{m}=2*pcm^2 (NOT FROM |t|_{min} to |t|_{max}).
      //THIS WILL WILL BE |t| of the particle flying forward (from 0 to 90 grad) in CM and u=-2*m_D*Tls-t.
      //THE LS MOMENTUM OF THE OTHER PARTICLE (from 90 to 180 grad) CAN BE FOUND FROM THE ENERGY CONSERVATION LAW.

      Fk.at(i).at(b)=psumb;
      TailFk.at(i).at(b)=psumb;
      RutherfordTailFk.at(i).at(b)=psumb_Rutherford;
    }
  }

  std::cout<<"Check TailsCS:"<<std::endl;
  for(int i=0; i<Bin2+1; ++i) std::cout<<TailCS.at(i)/mbarn<<" ";
  std::cout<<std::endl;
  std::cout<<"Check RutherfordTailCS:"<<std::endl;
  for(int i=0; i<Bin2+1; ++i) std::cout<<RutherfordTailCS.at(i)/mbarn<<" ";
  std::cout<<std::endl;

  std::cout<<"Check Fk:"<<std::endl;
  for(int i=0; i<Bin2+1; ++i)
  {
    std::cout<<"E="<<std::exp(lnE.at(i))/MeV<<" MeV"<<std::endl;
    for(int b=0; b<Bin3+1; ++b) std::cout<<Fk.at(i).at(b)<<" ";
    std::cout<<std::endl;
  }
  std::cout<<std::endl;
//--------------------------------------------------------------------------//
//prepare data for the file (E, partial sums):
//--------------------------------------------------------------------------//
  tps.resize(0);
  for(int i=0; i<Bin2+1; ++i)
  {
    std::cout<<"Bin #"<<i<<" E="<<std::exp(lnE.at(i))/MeV<<std::endl;
    Energy.at(i)=std::exp(lnE.at(i));
    std::cout<<"Fk:   ";
    for(int b=0; b<Bin3+1; ++b) std::cout<<Fk.at(i).at(b)<<" ";
    std::cout<<std::endl;
    T3NSGangular_RWnode newnode =
      T3NSGangular_RWnode(std::exp(lnE.at(i)), Fk.at(i));
    tps.push_back(newnode);
  }
  trw.resize(0);
  trw.push_back(tps);
  std::cout<<"End of Fill():"<<std::endl;
  std::cout<<" trw size="<<trw.size()<<" tps size="<<tps.size()<<" tps: "<<std::endl;
  std::cout<<trw<<std::endl;
}

//returns the differential cross section from d-d approximation in inner units
double T3R_DDCS::GetdSigmadt(double Tls/*LS kinetic energy of the incident deuteron in inner units*/,
                   double t/*|t| - Mandelstahm variable |t| in inner units*/) const
{
  //1. in aproximation Tls is in GeV. calculate tmax:
  const double pls2 = Tls*(2*md+Tls);
  const double pcm2 = pls2/(2.+2.*std::sqrt(1.+pls2/md2));//rat=md/md=1. rat2=1.
  const double tmax = 4*pcm2;//tmax=4*pcm2.
//if |t|>|t|max, exit:///
  if(t>tmax)
  {
    return 0.0;
  }
/////////////////////////
  //if(|t|>2*pcm^2) |t|=4*pcm^2-|t|.
  //d-d scattering of equal particles=>dsigma/d|t| is symmetric at 2*pcm^2=>
  //to get cross section at |t|>2*pcm^2 we should make
  //if(|t|>2*pcm^2) |t|=4*pcm^2-|t| and dsigma/d|t| will be dsigma/d|t|(4*pcm^2-|t|).
  if(t>2*pcm2) t=tmax-t;
  //Ar0 without hc
  const double Ar0=md / 137. * std::sqrt(10 * M_PI / pcm2);
  const double pcm4 = 4 * pcm2;
  const double sin2 = t / pcm4;
  const double cos2 = 1. - sin2;
  const double tg2 = sin2 / cos2;
  //Rutherford amplitude for equal particles d-d:
//---------------------------------------------------------------------------
//when pcm/Ecm ~ alpha:
//---------------------------------------------------------------------------  
  const double alpha=1.0/137.0;
  const double pcm=std::sqrt(pcm2);
// ////////////////////////////////////////////////
  //Here is an error: Ecm - CM energy of the inc particle or the reduced particle:
  //1) if inc particle - Ecm=std::sqrt(md2+pcm2)
  //2) if reduced particle - Ecm=std::sqrt(md2/4+pcm2)
// ////////////////////////////////////////////////
  const double Ecm=std::sqrt(md2+pcm2);
  const double CA=std::cos(alpha*Ecm/pcm*std::log(tg2));
//---------------------------------------------------------------------------  
  const double Ar = Ar0 * std::sqrt(1. + tg2 * tg2 + (2. / 3.) * tg2 * CA);
  const double Tcm = Tls / 2;//for d-d Tcm=Tls/2 in MeV.
  const double s = 4 * md2 * (1. + Tcm / md);//here Tcm in MeV
  //here Tcm in MeV
  const double g=std::exp(-0.5*std::sqrt(0.986 * MeV / Tcm));//Gamow factor//0.986 MeV - Gamow energy for d-d.
  const double As = 49000.0 * (1. + std::pow(Tcm/1.173, 3.06)) / (1. + std::pow(Tcm/0.64, 2.0)) /
                    (1. + std::pow(Tcm/5.45, 2.0)) /  (1.+std::pow(Tcm/23.7, 1.775)) * g;
//  As=0.49E5/(1.+(Tcm/0.64)**2.0)*
// *(1.+(Tcm/1.173)**3.06)/(1.+(Tcm/5.45)**2.)/
// /(1.+(Tcm/23.7)**1.775)*g
  const double Atu = 1041.0 * (1. + std::pow(Tcm/0.97, 3.89)) /
    (1. + std::pow(Tcm/0.056, 1.29)) / (1. + std::pow(Tcm/1.574, 2.6)) /
    (1. + std::pow(Tcm/16.75, 1.67)) * g;
//  Atu=1041.0/(1.+(Tcm/0.056)**1.29)*
// *(1.+(Tcm/0.97)**(1.29+2.6))/
// /(1.+(Tcm/1.574)**2.6)/
// /(1.+(Tcm/16.75)**1.67)*g  
  const double phases=std::pow(Tcm/0.93, 2.76) / (1. + std::pow(Tcm/1.63, 2.8));
//  phases=(Tcm/0.93)**2.76/
// /(1.+(Tcm/1.63)**2.8)
  const double phasetu=std::pow(Tcm/2.65, 3.1) / (1. + std::pow(Tcm/3.22, 3.245));
//  phasetu=(Tcm/2.65)**3.1/
// /(1.+(Tcm/3.22)**3.245)
  //squared mass of pi-0 mezon = 0.01822 GeV^2.
  const double mpi02=0.01822 * GeV * GeV;
  const double dtg = mpi02 / pcm4;
  const double tg2t = (sin2 + dtg) / (cos2 + dtg);
  const double termAtu = Atu * std::sqrt(1. + tg2t * tg2t + (2. / 3.) * tg2t * CA);
  const double term11 = Ar / t + As * std::cos(phases) / s;
  const double term12 = termAtu * std::cos(phasetu) / (t + mpi02);
  const double term1 = term11 + term12;
  const double term21 = As * std::sin(phases) / s;
  const double term22 = termAtu * std::sin(phasetu) / (t + mpi02);
  const double term2=term21+term22;
  const double hc = 0.2 * GeV * fm;//=200 MeV*fm.  
  double dsigmadt = hc*hc*(term1*term1+term2*term2);// =/mbarn
  //In approximation Ar0=amd/137.*SQRT(10*pi/pcm2).
  //10 is because the result is in fm^2
  //and 1fm^2=10mbarn, and to convert the result in mbarn/GeV^2 from fm^2/mbarn,
  //we multiply pi by 10. That is why, to get the result not in mbarn/GeV^2, but
  //in natural units, in fm^2/GeV^2,
  //we should divide the result of approximation by 10.
  dsigmadt/=10;//because std::sqrt(10 * M_PI / pcm2);.

  /*
  std::cout<<"md="<<md/GeV<<" md2="<<md2/GeV/GeV<<" Tls="<<Tls/MeV<<" MeV |t|="<<t/GeV/GeV<<" GeV^2"
           <<" pls2="<<pls2/GeV/GeV<<" GeV^2"<<" pcm2="<<pcm2/GeV/GeV<<" GeV^2"
           <<" tmax="<<tmax/GeV/GeV<<" GeV^2"<<" Ar0="<<Ar0<<" pcm4="<<pcm4/GeV/GeV<<" GeV^2"
           <<" sin2="<<sin2<<" cos2="<<cos2<<" tg2="<<tg2<<" alpha="<<alpha
           <<" pcm="<<pcm/GeV<<" GeV"<<" Ecm="<<Ecm/GeV<<" GeV"<<" CA="<<CA
           <<" Ar="<<Ar<<" Tcm="<<Tcm/MeV<<" MeV"<<" s="<<s/GeV/GeV<<" GeV^2"<<" g="<<g
           <<" As="<<As<<" Atu="<<Atu<<" fs="<<phases<<" ftu="<<phasetu<<" dtg="<<dtg<<" tg2t="<<tg2t
           <<" s="<<s/GeV/GeV<<" termAtu="<<termAtu<<" |t|="<<t/GeV/GeV
           <<" t11="<<term11*GeV*GeV<<" t12="<<term12*GeV*GeV
           <<" t21="<<term21*GeV*GeV<<" t22="<<term22*GeV*GeV
           <<" t1="<<term1*GeV*GeV<<" t2="<<term2*GeV*GeV
           <<" sum="<<(term1*term1+term2*term2)*GeV*GeV*GeV*GeV<<std::endl;
  */           
  ////std::cout<<"CS: dsigma/d|t|="<<dsigmadt/mbarn*GeV*GeV<<" mbarn/GeV^2"<<std::endl;
  return dsigmadt;
}

//returns the Rutherford differential cross section for d-d (Rutherford formula for equla particles)
double T3R_DDCS::GetdSigmadt_DD_Rutherford(double Tls/*LS kinetic energy of the incident deuteron in inner units*/,
                   double t/*|t| - Mandelstahm variable |t| in inner units*/) const
{
  //1. in aproximation Tls is in GeV. calculate tmax:
  const double pls2 = Tls*(2*md+Tls);
  const double pcm2 = pls2/(2.+2.*std::sqrt(1.+pls2/md2));//rat=md/md=1. rat2=1.
  const double tmax = 4*pcm2;//tmax=4*pcm2.
  //Ar0 without hc
  const double Ar0=md / 137. * std::sqrt(M_PI / pcm2);
  //Rutherford amplitude for equal particles d-d:
//---------------------------------------------------------------------------
//when pcm/Ecm ~ alpha:
//---------------------------------------------------------------------------
  const double alpha=1.0/137.0;
  const double pcm=std::sqrt(pcm2);
// ////////////////////////////////////////////////
  //Here is an error: Ecm - CM energy of the inc particle or the reduced particle:
  //1) if inc particle - Ecm=std::sqrt(md2+pcm2)
  //2) if reduced particle - Ecm=std::sqrt(md2/4+pcm2)
// ////////////////////////////////////////////////
  const double mr=md/2;
  const double mr2=mr*mr;
  const double Ecm=std::sqrt(mr2+pcm2);
  const double argln=t/(tmax-t);
  const double CA=std::cos(alpha*Ecm/pcm*std::log(argln));
//---------------------------------------------------------------------------
  const double term1 = 1.0/t/t;
  const double term2=1.0/(tmax-t)/(tmax-t);
  const double term3 = 2.0/3.0/t/(tmax-t)*CA;
  const double hc = 0.2 * GeV * fm;//=200 MeV*fm.
  const double coef=Ar0*Ar0;
  double dsigmadt = hc*hc*coef*(term1+term2+term3);// =/mbarn
  return dsigmadt;
}

void T3R_DDCS::Save_to_CS(int tgZ, int tgA/*deuteron Z=1 and A=2*/, T3String suffix)
{
  //target - deuteron:
  std::cout<<"Check ecs in Save_to_CS():"<<ecs<<std::endl;
  ecs.SetZ(tgZ);
  ecs.SetA(tgA);
  ecs.save_binary(1,2,1,2,suffix);
}

bool T3R_DDCS::Load_from_CS(int tgZ, int tgA/*deuteron Z=1 and A=2*/, T3String suffix)
{
  //target - deuteron:
  ecs.SetZ(tgZ);
  ecs.SetA(tgA);
  ecs.load_binary(1,2,1,2,suffix);
  return true;
}

void T3R_DDCS::Save_to_PS(/*target Z and A*/int tgZ, int tgA, const T3String & rid,
                          T3int pdg, T3int incZA)
{
  //target - deuteron:
  std::cout<<"Check trw in Save_to_PS():"<<trw<<std::endl;
  trw.Set_Z(tgZ);
  trw.Set_A(tgA);
  trw.save_binary(1,2,"elastic",1002,1002);
}

bool T3R_DDCS::Load_from_PS(/*target Z and A*/int tgZ, int tgA, const T3String & rid,
                          T3int pdg, T3int incZA)
{
  //target - deuteron:
  trw.Set_Z(tgZ);
  trw.Set_A(tgA);
  trw.load_binary(1,2,"elastic",1002,1002);
  std::cout<<"records in trw="<<trw.Get_vector_of_T3NSGangular_RWrecords().size()<<std::endl;
  for(int i=0; i<trw.Get_vector_of_T3NSGangular_RWrecords().size(); ++i)
  {
    std::cout<<"nodes in record["<<i<<"] in trw="<<trw.Get_vector_of_T3NSGangular_RWrecords().at(i).size()<<std::endl;
    std::cout<<"record["<<i<<"]:"<<std::endl;
    for(int j=0; j<trw.Get_vector_of_T3NSGangular_RWrecords().at(i).size(); ++j)
    {
      std::array<T3double, NumberOfBinsInDDApproximation3+1> vv=
              trw.Get_vector_of_T3NSGangular_RWrecords().at(i).Get_vector_of_T3NSGangular_RWnodes().at(j).Get_V();
      for(int k=0; k<NumberOfBinsInDDApproximation3+1; ++k) std::cout<<vv.at(k)<<" ";
      std::cout<<std::endl;

    }
  }

  std::cout<<"Check trw in Load_from_PS():"<<trw<<std::endl;
  return true;
}

//returns the integral cross section in inner units
//to get in mbarn - /mbarn
double T3R_DDCS::GetCS(double Tls/*LS kinetic energy of the inc deuteron*/) const
{
  if(Tls<100.0*eV || Tls>20.0*MeV) std::cout<<"***Error: Out of the dd energy range!"<<std::endl;
  const double lnTls=std::log(Tls)-lnE.at(0);
  const int bin=lnTls/deltalnE;
  std::cout<<"bin="<<bin<<std::endl;
  if(bin<0 || bin>Bin2-1) std::cout<<"***Error: bin="<<bin<<" bin+1="<<bin+1
                                 <<" Bin2="<<Bin2<<std::endl;
  std::cout<<"bin="<<bin<<std::endl;
  const double Ei=std::exp(lnE.at(bin));
  const double Ei1=std::exp(lnE.at(bin+1));
//for d-d we filled partial sums from |t|=0 to |t|=2*pcm^2
//because the particles are equal, and the differential
//cross section is symmetric at |t|=2*pcm^2. So the value
//of the partial sum for any energy in the last bin (b==Bin)
//is the half value of the integral cross section
//(from 0 to 4*pcm^2), because the differential cross section
//is symmetric at |t|=2*pcm^2.
  const double si=2*Fk.at(bin).at(Bin3);
  const double si1=2*Fk.at(bin+1).at(Bin3);
  const double cs = si + (Tls - Ei) * (si1 - si) / (Ei1 - Ei);
  std::cout<<"Tls="<<Tls/MeV<<" MeV lnTls="<<lnTls<<std::endl;
  std::cout<<"Check lnE bins: deltalnE="<<deltalnE<<std::endl;
  for(int b=0; b<Bin2+1; ++b) std::cout<<lnE.at(b)<<" ";
  std::cout<<std::endl;

  std::cout<<"Ei="<<Ei/MeV<<" MeV Ei1="<<Ei1/MeV<<" MeV si="<<si/barn<<" barn si1="<<si1/barn
           <<" barn sigma="<<cs/barn<<" barn"<<" Ei1-Ei="<<(Ei1 - Ei)/MeV<<" MeV"
           <<" E-Ei="<<(Tls-Ei)/MeV<<" MeV si1-si="<<(si1-si)/barn<<" barn"<<std::endl;
  return cs;
}

std::ostream& operator<<(std::ostream& os, const T3R_DDCS& inst)
{
  os<<"<< Bin2="<<T3R_DDCS::GetNumberOfBins2()<<" Bin3="<<inst.GetNumberOfBins3()<<std::endl;
  os<<"<< lnE:"<<std::endl;
  for(int b=0; b<inst.Get_size_2(); ++b) os<<inst.lnE.at(b)<<" ";
  os<<std::endl;
  os<<"<< partial sums:"<<std::endl;
  for(int i=0; i<inst.Get_size_2(); ++i)
  {
    for(int b=0; b<inst.Get_size_3(); ++b) os<<inst.Fk.at(i).at(b)/mbarn<<" ";
    os<<std::endl;
  }
  return os;
}

} // namespace t3
