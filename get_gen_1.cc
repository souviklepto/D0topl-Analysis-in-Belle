//___________________________________________________________________
// genT class to get generated information
// This class is useful in making the code very small for large combinations
// Also, you can Print the daughters of the particle.
// 
//   
// AUTHOR: Vishal  Place: Nara W. University 
//
// Version 0.9 : bug fixed for repetition
//               bug fixed in the returned HepLorentzVector 
//               Addition of momentum
//		 Simplify the logic to get daughters
//               Addition of more decays upto 13 daughters now
//               Addition of  PHOTOS, till 3 gammas are supported here
//                         These gammas can be reduced.
//               Printing all daughters option available
//               Rank option added to get the same charged particles
//               Absolute option added to make it easy
// Functions description:
// :o: cx_body( )
//          compares cda1 with da1, cda2 with da2,... and so on
// :o: dau (Gen_hepevt a)
//      Return the number of daughters
// :o: MyChild(Gen_hepevt a, int rank)
//      Return the daugther corresponsing to 1st, 2nd  and so on
// :o: checkDecay(Gen_hepevt a, int da1, int da2,..)
//      Check the lund id (da1,da2,...) are daughter of a
// :o: PcheckDecay(Gen_hepevt a, int da1, int da2,..)
//      Photos added and Check the lund id (da1,da2,...) are daughter of a
//      upto three gammas are tested in addition to da1,da2,...
// :o: AbscheckDecay(Gen_hepevt a, int da1, int da2,..)
//      Check the lund id of daughters but insensitive to their charge
// :o: PAbscheckDecay(Gen_hepevt a, int da1, int da2,..)
//      Photos and check the lund id of daughter but insensitive to their charge
// :o: AmId(Gen_hepevt a);
//      Gives the idhep() of the particle a
// :o: PrintDau(Gen_hepevt a);
//      Print all the daughters of a, Useful for debug (?)
// :o: DaME(Gen_hepevt a, int id);
//      Return daughter of a, after matching with its lund id.
// :o: DaME(Gen_hepevt a, int id, int rank);
//      Return daughter of a after matching its lund id and taking repetiton
//      Useful for pi+ B0-> chic1 K+ pi- pi+ pi- pi+ pi+, last pi+ has rank 3
// :o: AbsDaME(Gen_hepevt a, int id);
//      insensitive to chg
//      Return daughter of a, after matching with its lund id.
// :o: AbsDaME(Gen_hepevt a, int id, int rank);
//      Insensitive to chg
//      Return daughter of a after matching its lund id and taking repetiton
//      Useful for pi+ B0-> chic1 K+ pi- pi+ pi- pi+ pi+, last pi+ has rank 3
// :o: Mom4DAU(Gen_hepevt a, int id)
//      Returns HepLorentzVector of child of a having lund id
// :o: Mom4DAU(Gen_hepevt a, int id, int rank)
//      Returns HepLorentzVector of child of a having lund id
//      rank take care of repetition
// :o: AbsMom4DAU(Gen_hepevt a, int id)
//      Charge insensitive,Returns HepLorentzVector of child of a having lund id
// :o: AbsMom4DAU(Gen_hepevt a, int id, int rank)
//      Charge insenitive,Returns HepLorentzVector of child of a having lund id
//      rank take care of repetition//
// :o: Mom4(Gen_hepevt a);
//     Return HepLorentzVector of a
// :o: MomDAU(Gen_hepevt a, int id)
//      Returns momentum of child of a having lund id
// :o: MomDAU(Gen_hepevt a, int id, int rank)
//      Returns momentum of child of a having lund id
//      rank take care of repetition
// :o: AbsMomDAU(Gen_hepevt a, int id)
//      Charge insensitive,Returns momentum of child of a having lund id
// :o: AbsMomDAU(Gen_hepevt a, int id, int rank)
//      Charge insenitive,Returns momentum of child of a having lund id
//      rank take care of repetition//
// :o: Mom(Gen_hepevt a);
//     Return momentum of a//
// :o: pTDAU(Gen_hepevt a, int id)
//      Returns pT of child of a having lund id
// :o: pTDAU(Gen_hepevt a, int id, int rank)
//      Returns pT of child of a having lund id
//      rank take care of repetition
// :o: AbspTDAU(Gen_hepevt a, int id)
//      Charge insensitive,Returns pT of child of a having lund id
// :o: AbspTDAU(Gen_hepevt a, int id, int rank)
//      Charge insenitive,Returns pT of child of a having lund id
//      rank take care of repetition//
// :o: pT(Gen_hepevt a);
//     Return pT of a//
//
//___________________________________________________________________

#include "belle.h"
#include "get_gen.h"
#include HEPEVT_H
#include MDST_H
#include "particle/Ptype.h"
#include <math.h>

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif


  int genT::AmId(Gen_hepevt a){
    return a.idhep();
  }

  void genT::PrintDAU(Gen_hepevt a){
    std::cout << AmId(a) << " has " << dau(a) << " daughters" << std::endl;
    for(int j=0; j<dau(a); j++){
      std::cout << " Daughter " << j << " : " <<  AmId(MyChild(a,j+1)) << 
	std::endl;
    }
  }
 
  //Gives pT of HepLorentzVector
  double genT::pT(HepLorentzVector a){
     return sqrt(a.px()*a.px() + a.py()*a.py());
    
  }

 //Gives pT of HepLorentzVector
  double genT::pT(Gen_hepevt a){
     return sqrt(a.PX()*a.PX() + a.PY()*a.PY());
    
  }
  
 //Gives mom of HepLorentzVector
  double genT::Mom(HepLorentzVector a){
    return sqrt(a.px()*a.px() + a.py()*a.py() + a.pz()*a.pz());
  }
  

 //Gives mom of HepLorentzVector
  double genT::Mom(Gen_hepevt a){
    return sqrt(a.PX()*a.PX() + a.PY()*a.PY() + a.PZ()*a.PZ());
  }
 

  
  //Gives number of daughters of Generated particle
  int genT::dau(Gen_hepevt a){
    return a.da(1) - a.da(0) + 1;
  }
  

  HepLorentzVector genT::Mom4(Gen_hepevt a){
    HepLorentzVector p(a.PX(), a.PY(), a.PZ(),a.E());
    return p;
  }  

  
  //Gives the child of particle
  // rank is the number of kid
  // rank 1 will give first daughter
  // rank 2 will give second daughter
  // and so on...
  Gen_hepevt genT::MyChild(Gen_hepevt a, int rank){
    Gen_hepevt_Manager & hepevt_mag = Gen_hepevt_Manager::get_manager();
    Gen_hepevt &child = *(hepevt_mag.begin() + a.da(0) + (rank -2));
    return child;
  }
  
  
  // Returns the daughter after matching with its id,
  // WARNING :please use it only after checking that there is a matching 
  // daughter otherwise it will give segmentation violation :(
  Gen_hepevt genT::DaME(Gen_hepevt a, int id, int rank){
    int i=1;
    for(int j=0; j<dau(a); j++){
      if(MyChild(a,j+1).idhep()==id){if (rank==i){return MyChild(a,j+1);}
	else { i=i+1;} }
    }
  }
  
  //DAME without need of rank.
  // This function lives on the previous only the rank is put 1 here.
  Gen_hepevt genT::DaME(Gen_hepevt a, int id){
    return DaME(a,id,1);
  }

  // For Absolute, dont care about charge

  Gen_hepevt genT::AbsDaME(Gen_hepevt a, int id, int rank){
    int i=1;
    for(int j=0; j<dau(a); j++){
      if(abs(MyChild(a,j+1).idhep())==abs(id)){if (rank==i){return MyChild(a,j+1);}
	else { i=i+1;} }
    }
  }
  
  //DAME without need of rank.
  // This function lives on the previous only the rank is put 1 here.
  Gen_hepevt genT::AbsDaME(Gen_hepevt a, int id){
    return AbsDaME(a,id,1);
  }



  // Returns the HepLoretnzVector of daughter after matching with its id,
  // WARNING :please use it only after checking that there is a matching 
  // daughter otherwise it will give segmentation violation :(
  HepLorentzVector genT::Mom4DAU(Gen_hepevt a, int id, int rank){
    int i=1;
    for(int j=0; j<dau(a); j++){
      if(MyChild(a,j+1).idhep()==id){if (rank==i){return Mom4(MyChild(a,j+1));}
	else { i=i+1;} }
    }
  }
  
  //DAME without need of rank.
  // This function lives on the previous only the rank is put 1 here.
  HepLorentzVector genT::Mom4DAU(Gen_hepevt a, int id){
    return Mom4DAU(a,id,1);
  }
  
  // Abs doesn't care about the sign..
  // ...
  HepLorentzVector genT::AbsMom4DAU(Gen_hepevt a, int id, int rank){
    int i=1;
    for(int j=0; j<dau(a); j++){
      if(abs(MyChild(a,j+1).idhep())==abs(id)){if (rank==i){
	  return Mom4(MyChild(a,j+1));}
	else { i=i+1;} }
    }
  }
  
  //DAME without need of rank.
  // This function lives on the previous only the rank is put 1 here.
  HepLorentzVector genT::AbsMom4DAU(Gen_hepevt a, int id){
    return AbsMom4DAU(a,abs(id),1);
  }

  


 


  // Returns the momentum of daughter after matching with its id,
  // WARNING :please use it only after checking that there is a matching 
  // daughter otherwise it will give segmentation violation :(
  double genT::MomDAU(Gen_hepevt a, int id, int rank){
    return Mom(Mom4DAU(a,id,rank));
  }
  
  //DAME without need of rank.
  // This function lives on the previous only the rank is put 1 here.
  double genT::MomDAU(Gen_hepevt a, int id){
    return MomDAU(a,id,1);
  }
  
  //With Absolute
   double genT::AbsMomDAU(Gen_hepevt a, int id, int rank){
    return Mom(AbsMom4DAU(a,id,rank));
  }
  double genT::AbsMomDAU(Gen_hepevt a, int id){
    return AbsMomDAU(a,id,1);
  }

  // Returns the PT of daughter after matching with its id,
  // WARNING :please use it only after checking that there is a matching 
  // daughter otherwise it will give segmentation violation :(
  double genT::pTDAU(Gen_hepevt a, int id, int rank){
    return pT(Mom4DAU(a,id,rank));
  }
  
  // without need of rank.
  // This function lives on the previous only the rank is put 1 here.
  double genT::pTDAU(Gen_hepevt a, int id){
    return pTDAU(a,id,1);
  }

  //with Absolute
  double genT::AbspTDAU(Gen_hepevt a, int id, int rank){
    return pT(AbsMom4DAU(a,id,rank));
  }
  double genT::AbspTDAU(Gen_hepevt a, int id){
    return AbspTDAU(a,id,1);
  }


  int genT::checkDecay(Gen_hepevt a, int da1, int da2){
    if (dau(a)!=2){ return 0;}
    else {
      return c2_body(MyChild(a,1).idhep(), MyChild(a,2).idhep(), da1, da2);
    }
  }
  
  
  
  int genT::checkDecay(Gen_hepevt a, int da1, int da2, int da3){
    if (dau(a)!=3){ return 0;}
    else {
      return c3_body(MyChild(a,1).idhep(), MyChild(a,2).idhep(), 
		     MyChild(a,3).idhep(),da1, da2, da3);
    }
  }
  
  
  int genT::checkDecay(Gen_hepevt a, int da1, int da2, int da3, int da4){
    if (dau(a)!=4){ return 0;}
    else {
      return c4_body(MyChild(a,1).idhep(), MyChild(a,2).idhep(), 
		     MyChild(a,3).idhep(), MyChild(a,4).idhep(),
		     da1, da2, da3, da4);
    }
  }
  
  int genT::checkDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
		       int da5){
    if (dau(a)!=5){ return 0;}
    else {
      return c5_body(MyChild(a,1).idhep(), MyChild(a,2).idhep(), 
		     MyChild(a,3).idhep(), MyChild(a,4).idhep(),
		     MyChild(a,5).idhep(), da1, da2, da3, da4, da5);
    }
  }
  
  
 int genT::checkDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
		      int da5, int da6){
   if (dau(a)!=6){ return 0;}
   else {
     return c6_body(MyChild(a,1).idhep(), MyChild(a,2).idhep(), 
		    MyChild(a,3).idhep(), MyChild(a,4).idhep(),
		    MyChild(a,5).idhep(), MyChild(a,6).idhep(),
		    da1, da2, da3, da4, da5, da6);
   }
 }
  
  int genT::checkDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
		       int da5, int da6, int da7){
    if (dau(a)!=7){ return 0;}
    else {
      return c7_body(MyChild(a,1).idhep(), MyChild(a,2).idhep(), 
		     MyChild(a,3).idhep(), MyChild(a,4).idhep(),
		     MyChild(a,5).idhep(), MyChild(a,6).idhep(),
		     MyChild(a,7).idhep(), da1, da2, da3, da4, da5, 
		     da6, da7);
    }
  }
  
  int genT::checkDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
		       int da5, int da6, int da7, int da8){
    if (dau(a)!=8){ return 0;}
    else {
      return c8_body(MyChild(a,1).idhep(), MyChild(a,2).idhep(), 
		     MyChild(a,3).idhep(), MyChild(a,4).idhep(),
		     MyChild(a,5).idhep(), MyChild(a,6).idhep(),
		     MyChild(a,7).idhep(), MyChild(a,8).idhep(),
		     da1, da2, da3, da4, da5, da6, da7, da8);
    }
  }
  
  int genT::checkDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
		       int da5, int da6, int da7, int da8, int da9){
    if (dau(a)!=9){ return 0;}
    else {
      return c9_body(MyChild(a,1).idhep(), MyChild(a,2).idhep(), 
		     MyChild(a,3).idhep(), MyChild(a,4).idhep(),
		     MyChild(a,5).idhep(), MyChild(a,6).idhep(),
		     MyChild(a,7).idhep(), MyChild(a,8).idhep(),
		     MyChild(a,9).idhep(), da1, da2, da3, da4, da5, 
		     da6, da7, da8, da9);
    }
  }
  
  
  int genT::checkDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
		       int da5, int da6, int da7, int da8, int da9, int da10){
    if (dau(a)!=10){ return 0;}
    else {
      return c10_body(MyChild(a,1).idhep(), MyChild(a,2).idhep(), 
		      MyChild(a,3).idhep(), MyChild(a,4).idhep(),
		      MyChild(a,5).idhep(), MyChild(a,6).idhep(),
		      MyChild(a,7).idhep(), MyChild(a,8).idhep(),
		      MyChild(a,9).idhep(), MyChild(a,10).idhep(),
		      da1, da2, da3, da4, da5, da6, da7, da8, da9, da10);
    }
  }
  
  
  int genT::checkDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
		       int da5, int da6, int da7, int da8, int da9, int da10,
		       int da11 ){
    if (dau(a)!=11){ return 0;}
    else {
      return c11_body(MyChild(a,1).idhep(), MyChild(a,2).idhep(), 
		      MyChild(a,3).idhep(), MyChild(a,4).idhep(),
		      MyChild(a,5).idhep(), MyChild(a,6).idhep(),
		      MyChild(a,7).idhep(), MyChild(a,8).idhep(),
		      MyChild(a,9).idhep(), MyChild(a,10).idhep(),
		      MyChild(a,11).idhep(), da1, da2, da3, da4, da5, 
		      da6, da7, da8, da9, da10, da11);
    }
  }
  
  int genT::checkDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
		       int da5, int da6, int da7, int da8, int da9, int da10,
		       int da11, int da12 ){
    if (dau(a)!=12){ return 0;}
    else {
      return c12_body(MyChild(a,1).idhep(), MyChild(a,2).idhep(), 
		      MyChild(a,3).idhep(), MyChild(a,4).idhep(),
		      MyChild(a,5).idhep(), MyChild(a,6).idhep(),
		      MyChild(a,7).idhep(), MyChild(a,8).idhep(),
		      MyChild(a,9).idhep(), MyChild(a,10).idhep(),
		      MyChild(a,11).idhep(), MyChild(a,12).idhep(),
		      da1, da2, da3, da4, da5, da6, da7, da8, da9, 
		      da10, da11, da12);
    }
  }
  
  int genT::checkDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
		       int da5, int da6, int da7, int da8, int da9, int da10,
		       int da11, int da12, int da13 ){
    if (dau(a)!=13){ return 0;}
    else {
      return c13_body(MyChild(a,1).idhep(), MyChild(a,2).idhep(), 
		      MyChild(a,3).idhep(), MyChild(a,4).idhep(),
		      MyChild(a,5).idhep(), MyChild(a,6).idhep(),
		      MyChild(a,7).idhep(), MyChild(a,8).idhep(),
		      MyChild(a,9).idhep(), MyChild(a,10).idhep(),
		      MyChild(a,11).idhep(), MyChild(a,12).idhep(),
		      MyChild(a,13).idhep(), da1, da2, da3, da4, da5, 
		      da6, da7, da8, da9, da10, da11, da12, da13);
    }
  }
  
  
  
  // FOR ABSOLUTE don't care about any signs
  
  int genT::AbscheckDecay(Gen_hepevt a, int da1, int da2){
    if (dau(a)!=2){ return 0;}
    else {
      return c2_body(abs(MyChild(a,1).idhep()), abs(MyChild(a,2).idhep()), 
		     da1, da2);
    }
  }
  
  int genT::AbscheckDecay(Gen_hepevt a, int da1, int da2, int da3){
    if (dau(a)!=3){ return 0;}
    else {
      return c3_body(abs(MyChild(a,1).idhep()), abs(MyChild(a,2).idhep()), 
		     abs(MyChild(a,3).idhep()), da1, da2, da3);
    }
  }
  
  
  
  int genT::AbscheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4){
    if (dau(a)!=4){ return 0;}
    else {
      return c4_body(abs(MyChild(a,1).idhep()), abs(MyChild(a,2).idhep()), 
		     abs(MyChild(a,3).idhep()), abs(MyChild(a,4).idhep()),
		     da1, da2, da3, da4);
    }
  }
  
  int genT::AbscheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
			  int da5){
    if (dau(a)!=5){ return 0;}
    else {
      return c5_body(abs(MyChild(a,1).idhep()), abs(MyChild(a,2).idhep()), 
		     abs(MyChild(a,3).idhep()), abs(MyChild(a,4).idhep()),
		     abs(MyChild(a,5).idhep()), da1, da2, da3, da4, da5);
    }
  }
  
  
  int genT::AbscheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
		      int da5, int da6){
    if (dau(a)!=6){ return 0;}
    else {
      return c6_body(abs(MyChild(a,1).idhep()), abs(MyChild(a,2).idhep()), 
		     abs(MyChild(a,3).idhep()), abs(MyChild(a,4).idhep()),
		     abs(MyChild(a,5).idhep()), abs(MyChild(a,6).idhep()),
		     da1, da2, da3, da4, da5, da6);
    }
  }
  
  int genT::AbscheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
			  int da5, int da6, int da7){
    if (dau(a)!=7){ return 0;}
    else {
      return c7_body(abs(MyChild(a,1).idhep()), abs(MyChild(a,2).idhep()), 
		     abs(MyChild(a,3).idhep()), abs(MyChild(a,4).idhep()),
		     abs(MyChild(a,5).idhep()), abs(MyChild(a,6).idhep()),
		     abs(MyChild(a,7).idhep()), da1, da2, da3, da4, da5, 
		     da6, da7);
    }
  }
  
  int genT::AbscheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
			  int da5, int da6, int da7, int da8){
    if (dau(a)!=8){ return 0;}
    else {
      return c8_body(abs(MyChild(a,1).idhep()), abs(MyChild(a,2).idhep()), 
		     abs(MyChild(a,3).idhep()), abs(MyChild(a,4).idhep()),
		     abs(MyChild(a,5).idhep()), abs(MyChild(a,6).idhep()),
		     abs(MyChild(a,7).idhep()), abs(MyChild(a,8).idhep()),
		     da1, da2, da3, da4, da5, da6, da7, da8);
    }
  }
  
  int genT::AbscheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
			  int da5, int da6, int da7, int da8, int da9){
    if (dau(a)!=9){ return 0;}
    else {
      return c9_body(abs(MyChild(a,1).idhep()), abs(MyChild(a,2).idhep()), 
		     abs(MyChild(a,3).idhep()), abs(MyChild(a,4).idhep()),
		     abs(MyChild(a,5).idhep()), abs(MyChild(a,6).idhep()),
		     abs(MyChild(a,7).idhep()), abs(MyChild(a,8).idhep()),
		     abs(MyChild(a,9).idhep()), da1, da2, da3, da4, da5, 
		     da6, da7, da8, da9);
    }
  }
  
  
  int genT::AbscheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
			  int da5, int da6, int da7, int da8, int da9, int da10){
    if (dau(a)!=10){ return 0;}
    else {
      return c10_body(abs(MyChild(a,1).idhep()), abs(MyChild(a,2).idhep()), 
		      abs(MyChild(a,3).idhep()), abs(MyChild(a,4).idhep()),
		      abs(MyChild(a,5).idhep()), abs(MyChild(a,6).idhep()),
		      abs(MyChild(a,7).idhep()), abs(MyChild(a,8).idhep()),
		      abs(MyChild(a,9).idhep()), abs(MyChild(a,10).idhep()),
		      da1, da2, da3, da4, da5, da6, da7, da8, da9, da10);
    }
  }
  
  
  int genT::AbscheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
			  int da5, int da6, int da7, int da8, int da9, int da10,
			  int da11 ){
    if (dau(a)!=11){ return 0;}
    else {
      return c11_body(abs(MyChild(a,1).idhep()), abs(MyChild(a,2).idhep()), 
		      abs(MyChild(a,3).idhep()), abs(MyChild(a,4).idhep()),
		      abs(MyChild(a,5).idhep()), abs(MyChild(a,6).idhep()),
		      abs(MyChild(a,7).idhep()), abs(MyChild(a,8).idhep()),
		      abs(MyChild(a,9).idhep()), abs(MyChild(a,10).idhep()),
		      abs(MyChild(a,11).idhep()), da1, da2, da3, da4, da5, 
		      da6, da7, da8, da9, da10, da11);
    }
  }
  
  int genT::AbscheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
			  int da5, int da6, int da7, int da8, int da9, int da10,
			  int da11, int da12 ){
    if (dau(a)!=12){ return 0;}
    else {
      return c12_body(abs(MyChild(a,1).idhep()), abs(MyChild(a,2).idhep()), 
		      abs(MyChild(a,3).idhep()), abs(MyChild(a,4).idhep()),
		      abs(MyChild(a,5).idhep()), abs(MyChild(a,6).idhep()),
		      abs(MyChild(a,7).idhep()), abs(MyChild(a,8).idhep()),
		      abs(MyChild(a,9).idhep()), abs(MyChild(a,10).idhep()),
		      abs(MyChild(a,11).idhep()), abs(MyChild(a,12).idhep()),
		      da1, da2, da3, da4, da5, da6, da7, da8, da9, 
		      da10, da11, da12);
    }
  }
  
int genT::AbscheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
			int da5, int da6, int da7, int da8, int da9, int da10,
			int da11, int da12, int da13 ){
  if (dau(a)!=13){ return 0;}
  else {
    return c13_body(abs(MyChild(a,1).idhep()), abs(MyChild(a,2).idhep()), 
		    abs(MyChild(a,3).idhep()), abs(MyChild(a,4).idhep()),
		    abs(MyChild(a,5).idhep()), abs(MyChild(a,6).idhep()),
		    abs(MyChild(a,7).idhep()), abs(MyChild(a,8).idhep()),
		    abs(MyChild(a,9).idhep()), abs(MyChild(a,10).idhep()),
		    abs(MyChild(a,11).idhep()), abs(MyChild(a,12).idhep()),
		    abs(MyChild(a,13).idhep()), da1, da2, da3, da4, da5, 
		    da6, da7, da8, da9, da10, da11, da12, da13);
  }
}
  
  
  
  
  // PHOTOS
  
  //Check the number of daughters and add upto three in oder to 
  // take Photos into account
  int genT::Pcheck(Gen_hepevt a, int num){
    if (dau(a)!=num && dau(a)!=num+1 && dau(a)!=num+2 &&dau(a)!=num+3){
      return 0;}
    else return 1;
  }

  int genT::PcheckDecay(Gen_hepevt a, int da1, int da2){
    if (Pcheck(a,2)==0){ return 0;}
    else {
      if (dau(a)==2) { return checkDecay(a, da1, da2);}
      else if (dau(a)==3){return checkDecay(a,da1,da2,22);}
      else if (dau(a)==4){return checkDecay(a,da1,da2,22,22);}
      else if (dau(a)==5){return checkDecay(a,da1,da2,22,22,22);}
      else return 0;
    }
  }
  
  
  
  int genT::PcheckDecay(Gen_hepevt a, int da1, int da2, int da3){
   if (Pcheck(a,3)==0){ return 0;}
    else {
      if (dau(a)==3) { return checkDecay(a, da1, da2,da3);}
      else if (dau(a)==4){return checkDecay(a,da1,da2,da3,22);}
      else if (dau(a)==5){return checkDecay(a,da1,da2,da3,22,22);}
      else if (dau(a)==6){return checkDecay(a,da1,da2,da3,22,22,22);}
      else return 0;
    }
  }
  
  
  int genT::PcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4){
    if (Pcheck(a,4)==0){ return 0;}
    else {
      if (dau(a)==4) { return checkDecay(a, da1, da2,da3,da4);}
      else if (dau(a)==5){return checkDecay(a,da1,da2,da3,da4,22);}
      else if (dau(a)==6){return checkDecay(a,da1,da2,da3,da4,22,22);}
      else if (dau(a)==7){return checkDecay(a,da1,da2,da3,da4,22,22,22);}
      else return 0;
    }
  }
  
  int genT::PcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
		       int da5){
    if (Pcheck(a,5)==0){ return 0;}
    else {
      if (dau(a)==5) { return checkDecay(a, da1, da2,da3,da4,da5);}
      else if (dau(a)==6){return checkDecay(a,da1,da2,da3,da4,da5,22);}
      else if (dau(a)==7){return checkDecay(a,da1,da2,da3,da4,da5,22,22);}
      else if (dau(a)==8){return checkDecay(a,da1,da2,da3,da4,da5,22,22,22);}
      else return 0;
    }
  }
  
  
 int genT::PcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
		      int da5, int da6){
    if (Pcheck(a,6)==0){ return 0;}
    else {
      if (dau(a)==6) { return checkDecay(a, da1, da2,da3,da4,da5,da6);}
      else if (dau(a)==7){return checkDecay(a,da1,da2,da3,da4,da5,da6,22);}
      else if (dau(a)==8){return checkDecay(a,da1,da2,da3,da4,da5,da6,22,22);}
      else if (dau(a)==9){return checkDecay(a,da1,da2,da3,da4,da5,da6,
					    22,22,22);}
      else return 0;
   }
 }
  
  int genT::PcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
		       int da5, int da6, int da7){
    if (Pcheck(a,7)==0){ return 0;}
    else {
      if (dau(a)==7) { return checkDecay(a, da1, da2,da3,da4,da5,da6,da7);}
      else if (dau(a)==8){return checkDecay(a,da1,da2,da3,da4,da5,da6,da7,22);}
      else if (dau(a)==9){return checkDecay(a,da1,da2,da3,da4,da5,da6,da7,
					    22,22);}
      else if (dau(a)==10){return checkDecay(a,da1,da2,da3,da4,da5,da6,da7,
					    22,22,22);}
      else return 0;
    }
  }
  
  int genT::PcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
		       int da5, int da6, int da7, int da8){
    if (Pcheck(a,8)==0){ return 0;}
    else {
      if (dau(a)==8) { return checkDecay(a, da1, da2,da3,da4,da5,da6,da7,da8);}
      else if (dau(a)==9){return checkDecay(a,da1,da2,da3,da4,da5,da6,da7,da8,
					    22);}
      else if (dau(a)==10){return checkDecay(a,da1,da2,da3,da4,da5,da6,da7,da8,
					    22,22);}
      else if (dau(a)==11){return checkDecay(a,da1,da2,da3,da4,da5,da6,da7,da8,
					     22,22,22);}
      else return 0;
    }
  }
  
  int genT::PcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
			int da5, int da6, int da7, int da8, int da9,int da10){
    if (Pcheck(a,10)==0){ return 0;}
    else {
      if (dau(a)==10) { return checkDecay(a, da1, da2,da3,da4,da5,da6,da7,da8,
					  da9,da10,22);}
      else if (dau(a)==11){return checkDecay(a,da1,da2,da3,da4,da5,da6,da7,da8,
					     da9,da10,22);}
      else if (dau(a)==12){return checkDecay(a,da1,da2,da3,da4,da5,da6,da7,da8,
					     da9,da10,22,22);}
      else if (dau(a)==13){return checkDecay(a,da1,da2,da3,da4,da5,da6,da7,da8,
					     da9,da10,22,22,22);}
      else return 0;
    }
  }
  
  
  int genT::PcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
		       int da5, int da6, int da7, int da8, int da9){
    if (Pcheck(a,9)==0){ return 0;}
    else {
      if (dau(a)==9) { return checkDecay(a, da1, da2,da3,da4,da5,da6,da7,da8,
					 da9);}
      else if (dau(a)==10){return checkDecay(a,da1,da2,da3,da4,da5,da6,da7,da8,
					     da9,22);}
      else if (dau(a)==11){return checkDecay(a,da1,da2,da3,da4,da5,da6,da7,da8,
					     da9,22,22);}
      else if (dau(a)==12){return checkDecay(a,da1,da2,da3,da4,da5,da6,da7,da8,
					     da9,22,22,22);}
      else return 0;
    }
  }
  
  
  // Absolute Photos

   int genT::AbsPcheckDecay(Gen_hepevt a, int da1, int da2){
    if (Pcheck(a,2)==0){ return 0;}
    else {
      if (dau(a)==2) { return AbscheckDecay(a, da1, da2);}
      else if (dau(a)==3){return AbscheckDecay(a,da1,da2,22);}
      else if (dau(a)==4){return AbscheckDecay(a,da1,da2,22,22);}
      else if (dau(a)==5){return AbscheckDecay(a,da1,da2,22,22,22);}
      else return 0;
    }
  }
  
  
  
  int genT::AbsPcheckDecay(Gen_hepevt a, int da1, int da2, int da3){
   if (Pcheck(a,3)==0){ return 0;}
    else {
      if (dau(a)==3) { return AbscheckDecay(a, da1, da2,da3);}
      else if (dau(a)==4){return AbscheckDecay(a,da1,da2,da3,22);}
      else if (dau(a)==5){return AbscheckDecay(a,da1,da2,da3,22,22);}
      else if (dau(a)==6){return AbscheckDecay(a,da1,da2,da3,22,22,22);}
      else return 0;
    }
  }
  
  
  int genT::AbsPcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4){
    if (Pcheck(a,4)==0){ return 0;}
    else {
      if (dau(a)==4) { return AbscheckDecay(a, da1, da2,da3,da4);}
      else if (dau(a)==5){return AbscheckDecay(a,da1,da2,da3,da4,22);}
      else if (dau(a)==6){return AbscheckDecay(a,da1,da2,da3,da4,22,22);}
      else if (dau(a)==7){return AbscheckDecay(a,da1,da2,da3,da4,22,22,22);}
      else return 0;
    }
  }
  
  int genT::AbsPcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
		       int da5){
    if (Pcheck(a,5)==0){ return 0;}
    else {
      if (dau(a)==5) { return AbscheckDecay(a, da1, da2,da3,da4,da5);}
      else if (dau(a)==6){return AbscheckDecay(a,da1,da2,da3,da4,da5,22);}
      else if (dau(a)==7){return AbscheckDecay(a,da1,da2,da3,da4,da5,22,22);}
      else if (dau(a)==8){return AbscheckDecay(a,da1,da2,da3,da4,da5,22,22,22);}
      else return 0;
    }
  }
  
  
 int genT::AbsPcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
		      int da5, int da6){
    if (Pcheck(a,6)==0){ return 0;}
    else {
      if (dau(a)==6) { return AbscheckDecay(a, da1, da2,da3,da4,da5,da6);}
      else if (dau(a)==7){return AbscheckDecay(a,da1,da2,da3,da4,da5,da6,22);}
      else if (dau(a)==8){return AbscheckDecay(a,da1,da2,da3,da4,da5,da6,22,22);}
      else if (dau(a)==9){return AbscheckDecay(a,da1,da2,da3,da4,da5,da6,
					    22,22,22);}
      else return 0;
   }
 }
  
  int genT::AbsPcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
		       int da5, int da6, int da7){
    if (Pcheck(a,7)==0){ return 0;}
    else {
      if (dau(a)==7) { return AbscheckDecay(a, da1, da2,da3,da4,da5,da6,da7);}
      else if (dau(a)==8){return AbscheckDecay(a,da1,da2,da3,da4,da5,da6,da7,
					       22);}
      else if (dau(a)==9){return AbscheckDecay(a,da1,da2,da3,da4,da5,da6,da7,
					    22,22);}
      else if (dau(a)==10){return AbscheckDecay(a,da1,da2,da3,da4,da5,da6,da7,
					    22,22,22);}
      else return 0;
    }
  }
  
  int genT::AbsPcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
		       int da5, int da6, int da7, int da8){
    if (Pcheck(a,8)==0){ return 0;}
    else {
      if (dau(a)==8) { return AbscheckDecay(a, da1, da2,da3,da4,da5,da6,da7,da8);}
      else if (dau(a)==9){return AbscheckDecay(a,da1,da2,da3,da4,da5,da6,da7,da8,
					    22);}
      else if (dau(a)==10){return AbscheckDecay(a,da1,da2,da3,da4,da5,da6,da7,da8,
					    22,22);}
      else if (dau(a)==11){return AbscheckDecay(a,da1,da2,da3,da4,da5,da6,da7,da8,
					     22,22,22);}
      else return 0;
    }
  }
  
  int genT::AbsPcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
		       int da5, int da6, int da7, int da8, int da9){
    if (Pcheck(a,10)==0){ return 0;}
    else {
      if (dau(a)==10) { return AbscheckDecay(a, da1, da2,da3,da4,da5,da6,da7,
					     da8, da9);}
      else if (dau(a)==11){return AbscheckDecay(a,da1,da2,da3,da4,da5,da6,da7,
						da8, da9,22);}
      else if (dau(a)==12){return AbscheckDecay(a,da1,da2,da3,da4,da5,da6,da7,
						da8, da9,22,22);}
      else if (dau(a)==13){return AbscheckDecay(a,da1,da2,da3,da4,da5,da6,da7,
						da8, da9,22,22,22);}
      else return 0;
    }
  }
  
  
  int genT::AbsPcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, 
		       int da5, int da6, int da7, int da8, int da9, int da10){
    if (Pcheck(a,9)==0){ return 0;}
    else {
      if (dau(a)==9) { return AbscheckDecay(a, da1, da2,da3,da4,da5,da6,da7,da8,
					 da9);}
      else if (dau(a)==10){return AbscheckDecay(a,da1,da2,da3,da4,da5,da6,da7,
						da8, da9,22);}
      else if (dau(a)==11){return AbscheckDecay(a,da1,da2,da3,da4,da5,da6,da7,
						da8,da9,22,22);}
      else if (dau(a)==12){return AbscheckDecay(a,da1,da2,da3,da4,da5,da6,da7,
						da8, da9,22,22,22);}
      else return 0;
    }
  }
  



  int genT::c2_body(int cda1, int cda2, int da1, int da2 ){
    if ((cda1==da1 && cda2==da2) || (cda1==da2 && cda2==da1)) 
      { return 1;}
    else return 0;
  }


  int genT::c3_body(int cda1, int cda2, int cda3,
		      int da1, int da2, int da3 ){
    if (cda1==da1 && cda2==da2 && cda3==da3 ||
	cda1==da1 && cda3==da2 && cda2==da3 ||
	cda2==da1 && cda1==da2 && cda3==da3 ||
	cda2==da1 && cda3==da2 && cda1==da3 ||
	cda3==da1 && cda1==da2 && cda2==da3 || 
	cda3==da1 && cda2==da2 && cda1==da3) 
      { return 1;}
    else return 0;
  }
  

  int genT::c4_body(int cda1, int cda2, int cda3, int cda4,
		      int da1, int da2, int da3, int da4 ){

    if (cda1==da1 && cda2==da2 && cda3==da3 && cda4==da4  ||
	cda1==da1 && cda2==da2 && cda4==da3 && cda3==da4  ||
	cda1==da1 && cda3==da2 && cda2==da3 && cda4==da4  ||
	cda1==da1 && cda3==da2 && cda4==da3 && cda2==da4  ||
	cda1==da1 && cda4==da2 && cda3==da3 && cda2==da4  ||
	cda1==da1 && cda4==da2 && cda2==da3 && cda3==da4  ||
	
	cda2==da1 && cda1==da2 && cda3==da3 && cda4==da4  ||
	cda2==da1 && cda1==da2 && cda4==da3 && cda3==da4  ||
	cda2==da1 && cda3==da2 && cda1==da3 && cda4==da4  ||
	cda2==da1 && cda3==da2 && cda4==da3 && cda1==da4  ||
	cda2==da1 && cda4==da2 && cda3==da3 && cda1==da4  ||
	cda2==da1 && cda4==da2 && cda1==da3 && cda3==da4  ||
	
	cda3==da1 && cda2==da2 && cda1==da3 && cda4==da4  ||
	cda3==da1 && cda2==da2 && cda4==da3 && cda1==da4  ||
	cda3==da1 && cda1==da2 && cda2==da3 && cda4==da4  ||
	cda3==da1 && cda1==da2 && cda4==da3 && cda2==da4  ||
	cda3==da1 && cda4==da2 && cda1==da3 && cda2==da4  ||
	cda3==da1 && cda4==da2 && cda2==da3 && cda1==da4  ||
		
	cda4==da1 && cda2==da2 && cda3==da3 && cda1==da4  ||
	cda4==da1 && cda2==da2 && cda1==da3 && cda3==da4  ||
	cda4==da1 && cda3==da2 && cda2==da3 && cda1==da4  ||
	cda4==da1 && cda3==da2 && cda1==da3 && cda2==da4  ||
	cda4==da1 && cda1==da2 && cda3==da3 && cda2==da4  ||
	cda4==da1 && cda1==da2 && cda2==da3 && cda3==da4  )
      { return 1;}
    else return 0;
  }





  int genT::c5_body(int cda1, int cda2, int cda3, int cda4, int cda5,
		       int da1, int da2, int da3, int da4, int da5 ){
    int c_da[5]={da1,da2,da3,da4,da5};
    int flag=0;
    for (int i1=0; i1<5; i1++){
      for( int i2=0; i2<5; i2++){
	if (i1==i2) continue;
	for( int i3=0; i3<5; i3++){
	  if (i3==i1) continue;
	  if (i3==i2) continue;
	  for(int i4=0; i4<5; i4++){
	    if(i4==i1) continue;
	    if(i4==i2) continue;
	    if(i4==i3) continue;
	    for(int i5=0; i5<5; i5++){
	      if(i5==i1) continue;
	      if(i5==i2) continue;
	      if(i5==i3) continue;
	      if(i5==i4) continue;
	      if (cda1==c_da[i1] && cda2==c_da[i2] && cda3==c_da[i3] &&
		  cda4==c_da[i4] && cda5==c_da[i5] )
		    { flag=1; return 1; }
	    }
	  }
	}
      }
    }
    if(flag !=1){
      return 0;}
  }
  


  int genT::c6_body(int cda1, int cda2, int cda3, int cda4, int cda5, 
		       int cda6,
		       int da1, int da2, int da3, int da4, int da5, int da6 ){
    
    int c_da[6]={da1,da2,da3,da4,da5,da6};
    int flag=0;
    for (int i1=0; i1<6; i1++){
      for( int i2=0; i2<6; i2++){
	if (i1==i2) continue;
	for( int i3=0; i3<6; i3++){
	  if (i3==i1) continue;
	  if (i3==i2) continue;
	  for(int i4=0; i4<6; i4++){
	    if(i4==i1) continue;
	    if(i4==i2) continue;
	    if(i4==i3) continue;
	    for(int i5=0; i5<6; i5++){
	      if(i5==i1) continue;
	      if(i5==i2) continue;
	      if(i5==i3) continue;
	      if(i5==i4) continue;
	      for(int i6=0; i6<6; i6++){
		if(i6==i1) continue;
		if(i6==i2) continue;
		if(i6==i3) continue;
		if(i6==i4) continue;
		if(i6==i5) continue;
		if (cda1==c_da[i1] && cda2==c_da[i2] && cda3==c_da[i3] &&
		    cda4==c_da[i4] && cda5==c_da[i5] && cda6==c_da[i6])
		  { flag=1; return 1; }
	      }
	    }
	  }
	}
      }
    }
    if(flag !=1){
      return 0;}
  }
  




  int genT::c7_body(int cda1, int cda2, int cda3, int cda4, int cda5, int cda6, int cda7,
		    int da1, int da2, int da3, int da4, int da5, int da6, int da7 ){
    int c_da[7]={da1,da2,da3,da4,da5,da6,da7};
    int flag=0;
    for (int i1=0; i1<7; i1++){
      for( int i2=0; i2<7; i2++){
	if (i1 == i2) continue;
	for( int i3=0; i3<7; i3++){
	  if (i3== i1) continue;
	  if (i3== i2) continue;
	  for(int i4=0; i4<7; i4++){
	    if(i4==i1) continue;
	    if(i4==i2) continue;
	    if(i4==i3) continue;
	    for(int i5=0; i5<7; i5++){
	      if(i5==i1) continue;
	      if(i5==i2) continue;
	      if(i5==i3) continue;
	      if(i5==i4) continue;
	      for(int i6=0; i6<7; i6++){
		if(i6==i1) continue;
		if(i6==i2) continue;
		if(i6==i3) continue;
		if(i6==i4) continue;
		if(i6==i5) continue;
		for(int i7=0; i7<7; i7++){
		  if(i7==i1) continue;
		  if(i7==i2) continue;
		  if(i7==i3) continue;
		  if(i7==i4) continue;
		  if(i7==i5) continue;
		  if(i7==i6) continue;
		  if (cda1==c_da[i1] && cda2==c_da[i2] && cda3==c_da[i3] &&
		      cda4==c_da[i4] && cda5==c_da[i5] && cda6==c_da[i6] &&  
		      cda7==c_da[i7])
		    { flag=1; return 1; }
		}
	      }
	    }
	  }
	}
      }
    }
    if(flag !=1){
      return 0;}
  }
  









  int genT::c8_body(int cda1, int cda2, int cda3, int cda4, int cda5, int cda6, int cda7, int cda8,
		    int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8 ){
    int c_da[8]={da1,da2,da3,da4,da5,da6,da7, da8};
    int flag=0;
    for (int i1=0; i1<8; i1++){
      for( int i2=0; i2<8; i2++){
	if (i1 == i2) continue;
	for( int i3=0; i3<8; i3++){
	  if (i3== i1) continue;
	  if (i3== i2) continue;
	  for(int i4=0; i4<8; i4++){
	    if(i4==i1) continue;
	    if(i4==i2) continue;
	    if(i4==i3) continue;
	    for(int i5=0; i5<8; i5++){
	      if(i5==i1) continue;
	      if(i5==i2) continue;
	      if(i5==i3) continue;
	      if(i5==i4) continue;
	      for(int i6=0; i6<8; i6++){
		if(i6==i1) continue;
		if(i6==i2) continue;
		if(i6==i3) continue;
		if(i6==i4) continue;
		if(i6==i5) continue;
		for(int i7=0; i7<8; i7++){
		  if(i7==i1) continue;
		  if(i7==i2) continue;
		  if(i7==i3) continue;
		  if(i7==i4) continue;
		  if(i7==i5) continue;
		  if(i7==i6) continue;
		  for(int i8=0; i8<8; i8++){
		  if(i8==i1) continue;
		  if(i8==i2) continue;
		  if(i8==i3) continue;
		  if(i8==i4) continue;
		  if(i8==i5) continue;
		  if(i8==i6) continue;
		  if(i8==i7) continue;
		  if (cda1==c_da[i1] && cda2==c_da[i2] && cda3==c_da[i3] &&
		      cda4==c_da[i4] && cda5==c_da[i5] && cda6==c_da[i6] &&  
		      cda7==c_da[i7] && cda8==c_da[i8])
		    { flag=1; return 1; }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    if(flag !=1){
      return 0;}
  }
  





  int genT::c9_body(int cda1, int cda2, int cda3, int cda4, int cda5, int cda6, int cda7, int cda8, int cda9,
		    int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8, int da9 ){
    int c_da[9]={da1,da2,da3,da4,da5,da6,da7,da8,da9};
    int flag=0; int maX = 9;
    for (int i1=0; i1<maX; i1++){
      for( int i2=0; i2<maX; i2++){
	if (i1 == i2) continue;
	for( int i3=0; i3<maX; i3++){
	  if (i3== i1) continue;
	  if (i3== i2) continue;
	  for(int i4=0; i4<maX; i4++){
	    if(i4==i1) continue;
	    if(i4==i2) continue;
	    if(i4==i3) continue;
	    for(int i5=0; i5<maX; i5++){
	      if(i5==i1) continue;
	      if(i5==i2) continue;
	      if(i5==i3) continue;
	      if(i5==i4) continue;
	      for(int i6=0; i6<maX; i6++){
		if(i6==i1) continue;
		if(i6==i2) continue;
		if(i6==i3) continue;
		if(i6==i4) continue;
		if(i6==i5) continue;
		for(int i7=0; i7<maX; i7++){
		  if(i7==i1) continue;
		  if(i7==i2) continue;
		  if(i7==i3) continue;
		  if(i7==i4) continue;
		  if(i7==i5) continue;
		  if(i7==i6) continue;
		  for(int i8=0; i8<maX; i8++){
		  if(i8==i1) continue;
		  if(i8==i2) continue;
		  if(i8==i3) continue;
		  if(i8==i4) continue;
		  if(i8==i5) continue;
		  if(i8==i6) continue;
		  if(i8==i7) continue;
		  for(int i9=0; i9<maX; i9++){
		    if(i9==i1) continue;
		    if(i9==i2) continue;
		    if(i9==i3) continue;
		    if(i9==i4) continue;
		    if(i9==i5) continue;
		    if(i9==i6) continue;
		    if(i9==i7) continue;
		    if(i9==i8) continue;
		    if (cda1==c_da[i1] && cda2==c_da[i2] && cda3==c_da[i3] &&
			cda4==c_da[i4] && cda5==c_da[i5] && cda6==c_da[i6] &&  
			cda7==c_da[i7] && cda8==c_da[i8] && cda9==c_da[i9])
		      { flag=1; return 1; }
		  }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    if(flag !=1){
      return 0;}
  }
  



  int genT::c10_body(int cda1, int cda2, int cda3, int cda4, int cda5, int cda6, int cda7, int cda8, int cda9, int cda10,
		     int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8, int da9, int da10 ){
    int c_da[10]={da1,da2,da3,da4,da5,da6,da7,da8,da9,da10};
    int flag=0; int maX = 10;
    for (int i1=0; i1<maX; i1++){
      for( int i2=0; i2<maX; i2++){
	if (i1 == i2) continue;
	for( int i3=0; i3<maX; i3++){
	  if (i3== i1) continue;
	  if (i3== i2) continue;
	  for(int i4=0; i4<maX; i4++){
	    if(i4==i1) continue;
	    if(i4==i2) continue;
	    if(i4==i3) continue;
	    for(int i5=0; i5<maX; i5++){
	      if(i5==i1) continue;
	      if(i5==i2) continue;
	      if(i5==i3) continue;
	      if(i5==i4) continue;
	      for(int i6=0; i6<maX; i6++){
		if(i6==i1) continue;
		if(i6==i2) continue;
		if(i6==i3) continue;
		if(i6==i4) continue;
		if(i6==i5) continue;
		for(int i7=0; i7<maX; i7++){
		  if(i7==i1) continue;
		  if(i7==i2) continue;
		  if(i7==i3) continue;
		  if(i7==i4) continue;
		  if(i7==i5) continue;
		  if(i7==i6) continue;
		  for(int i8=0; i8<maX; i8++){
		  if(i8==i1) continue;
		  if(i8==i2) continue;
		  if(i8==i3) continue;
		  if(i8==i4) continue;
		  if(i8==i5) continue;
		  if(i8==i6) continue;
		  if(i8==i7) continue;
		  for(int i9=0; i9<maX; i9++){
		    if(i9==i1) continue;
		    if(i9==i2) continue;
		    if(i9==i3) continue;
		    if(i9==i4) continue;
		    if(i9==i5) continue;
		    if(i9==i6) continue;
		    if(i9==i7) continue;
		    if(i9==i8) continue;
		    for(int i10=0; i10<maX; i10++){
		    if(i10==i1) continue;
		    if(i10==i2) continue;
		    if(i10==i3) continue;
		    if(i10==i4) continue;
		    if(i10==i5) continue;
		    if(i10==i6) continue;
		    if(i10==i7) continue;
		    if(i10==i8) continue;
		    if(i10==i9) continue;
		    if (cda1==c_da[i1] && cda2==c_da[i2] && cda3==c_da[i3] &&
			cda4==c_da[i4] && cda5==c_da[i5] && cda6==c_da[i6] &&  
			cda7==c_da[i7] && cda8==c_da[i8] && cda9==c_da[i9] &&
			cda10==c_da[i10])
		      { flag=1; return 1; }
		    }
		  }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    if(flag !=1){
      return 0;}
  }
  



  int genT::c11_body(int cda1, int cda2, int cda3, int cda4, int cda5, 
		     int cda6, int cda7, int cda8, int cda9, int cda10, 
		     int cda11, int da1, int da2, int da3, int da4, 
		     int da5, int da6, int da7, int da8, int da9, 
		     int da10, int da11 ){
    int c_da[11]={da1,da2,da3,da4,da5,da6,da7,da8,da9,da10,da11};
    int flag=0; int maX = 11;
    for (int i1=0; i1<maX; i1++){
      for( int i2=0; i2<maX; i2++){
	if (i1 == i2) continue;
	for( int i3=0; i3<maX; i3++){
	  if (i3== i1) continue;
	  if (i3== i2) continue;
	  for(int i4=0; i4<maX; i4++){
	    if(i4==i1) continue;
	    if(i4==i2) continue;
	    if(i4==i3) continue;
	    for(int i5=0; i5<maX; i5++){
	      if(i5==i1) continue;
	      if(i5==i2) continue;
	      if(i5==i3) continue;
	      if(i5==i4) continue;
	      for(int i6=0; i6<maX; i6++){
		if(i6==i1) continue;
		if(i6==i2) continue;
		if(i6==i3) continue;
		if(i6==i4) continue;
		if(i6==i5) continue;
		for(int i7=0; i7<maX; i7++){
		  if(i7==i1) continue;
		  if(i7==i2) continue;
		  if(i7==i3) continue;
		  if(i7==i4) continue;
		  if(i7==i5) continue;
		  if(i7==i6) continue;
		  for(int i8=0; i8<maX; i8++){
		  if(i8==i1) continue;
		  if(i8==i2) continue;
		  if(i8==i3) continue;
		  if(i8==i4) continue;
		  if(i8==i5) continue;
		  if(i8==i6) continue;
		  if(i8==i7) continue;
		  for(int i9=0; i9<maX; i9++){
		    if(i9==i1) continue;
		    if(i9==i2) continue;
		    if(i9==i3) continue;
		    if(i9==i4) continue;
		    if(i9==i5) continue;
		    if(i9==i6) continue;
		    if(i9==i7) continue;
		    if(i9==i8) continue;
		    for(int i10=0; i10<maX; i10++){
		    if(i10==i1) continue;
		    if(i10==i2) continue;
		    if(i10==i3) continue;
		    if(i10==i4) continue;
		    if(i10==i5) continue;
		    if(i10==i6) continue;
		    if(i10==i7) continue;
		    if(i10==i8) continue;
		    if(i10==i9) continue;
		    for(int i11=0; i11<maX; i11++){
		      if(i11==i1) continue;
		      if(i11==i2) continue;
		      if(i11==i3) continue;
		      if(i11==i4) continue;
		      if(i11==i5) continue;
		      if(i11==i6) continue;
		      if(i11==i7) continue;
		      if(i11==i8) continue;
		      if(i11==i9) continue;
		      if(i11==i10) continue;
		      
		    if (cda1==c_da[i1] && cda2==c_da[i2] && cda3==c_da[i3] &&
			cda4==c_da[i4] && cda5==c_da[i5] && cda6==c_da[i6] &&  
			cda7==c_da[i7] && cda8==c_da[i8] && cda9==c_da[i9] &&
			cda10==c_da[i10] && cda11==c_da[i11])
		      { flag=1; return 1; }
		    }
		    }
		  }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    if(flag !=1){
      return 0;}
  }
  


  int genT::c12_body(int cda1, int cda2, int cda3, int cda4, int cda5, 
		     int cda6, int cda7, int cda8, int cda9, int cda10, 
		     int cda11, int cda12, int da1, int da2, int da3, 
		     int da4, int da5, int da6, int da7, int da8, int da9, 
		     int da10, int da11, int da12){
    int c_da[12]={da1,da2,da3,da4,da5,da6,da7,da8,da9,da10,da11,da12};
    int flag=0; int maX = 12;
    for (int i1=0; i1<maX; i1++){
      for( int i2=0; i2<maX; i2++){
	if (i1 == i2) continue;
	for( int i3=0; i3<maX; i3++){
	  if (i3== i1) continue;
	  if (i3== i2) continue;
	  for(int i4=0; i4<maX; i4++){
	    if(i4==i1) continue;
	    if(i4==i2) continue;
	    if(i4==i3) continue;
	    for(int i5=0; i5<maX; i5++){
	      if(i5==i1) continue;
	      if(i5==i2) continue;
	      if(i5==i3) continue;
	      if(i5==i4) continue;
	      for(int i6=0; i6<maX; i6++){
		if(i6==i1) continue;
		if(i6==i2) continue;
		if(i6==i3) continue;
		if(i6==i4) continue;
		if(i6==i5) continue;
		for(int i7=0; i7<maX; i7++){
		  if(i7==i1) continue;
		  if(i7==i2) continue;
		  if(i7==i3) continue;
		  if(i7==i4) continue;
		  if(i7==i5) continue;
		  if(i7==i6) continue;
		  for(int i8=0; i8<maX; i8++){
		  if(i8==i1) continue;
		  if(i8==i2) continue;
		  if(i8==i3) continue;
		  if(i8==i4) continue;
		  if(i8==i5) continue;
		  if(i8==i6) continue;
		  if(i8==i7) continue;
		  for(int i9=0; i9<maX; i9++){
		    if(i9==i1) continue;
		    if(i9==i2) continue;
		    if(i9==i3) continue;
		    if(i9==i4) continue;
		    if(i9==i5) continue;
		    if(i9==i6) continue;
		    if(i9==i7) continue;
		    if(i9==i8) continue;
		    for(int i10=0; i10<maX; i10++){
		    if(i10==i1) continue;
		    if(i10==i2) continue;
		    if(i10==i3) continue;
		    if(i10==i4) continue;
		    if(i10==i5) continue;
		    if(i10==i6) continue;
		    if(i10==i7) continue;
		    if(i10==i8) continue;
		    if(i10==i9) continue;
		    for(int i11=0; i11<maX; i11++){
		      if(i11==i1) continue;
		      if(i11==i2) continue;
		      if(i11==i3) continue;
		      if(i11==i4) continue;
		      if(i11==i5) continue;
		      if(i11==i6) continue;
		      if(i11==i7) continue;
		      if(i11==i8) continue;
		      if(i11==i9) continue;
		      if(i11==i10) continue;
		      for(int i12=0; i12<maX; i12++){
			if(i12==i1) continue;
			if(i12==i2) continue;
			if(i12==i3) continue;
			if(i12==i4) continue;
			if(i12==i5) continue;
			if(i12==i6) continue;
			if(i12==i7) continue;
			if(i12==i8) continue;
			if(i12==i9) continue;
			if(i12==i10) continue;
			if(i12==i11) continue;
			if (cda1==c_da[i1]&&cda2==c_da[i2]&&cda3==c_da[i3] &&
			    cda4==c_da[i4]&&cda5==c_da[i5]&&cda6==c_da[i6] &&  
			    cda7==c_da[i7]&&cda8==c_da[i8]&&cda9==c_da[i9] &&
			    cda10==c_da[i10]&&cda11==c_da[i11]
			    &&cda12==c_da[i12])
			  { flag=1; return 1; }
		      }
		    }
		    }
		  }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    if(flag !=1){
      return 0;}
  }
  
  
  
  
  
  
  int genT::c13_body(int cda1, int cda2, int cda3, int cda4, int cda5, 
		     int cda6, int cda7, int cda8, int cda9, int cda10, 
		     int cda11, int cda12,int cda13, int da1, int da2, 
		     int da3, int da4, int da5, int da6, int da7, 
		     int da8, int da9, int da10, int da11, int da12, int da13){
    int c_da[13]={da1,da2,da3,da4,da5,da6,da7,da8,da9,da10,da11,da12,da13};
    int flag=0; int maX = 13;
    for (int i1=0; i1<maX; i1++){
      for( int i2=0; i2<maX; i2++){
	if (i1 == i2) continue;
	for( int i3=0; i3<maX; i3++){
	  if (i3== i1) continue;
	  if (i3== i2) continue;
	  for(int i4=0; i4<maX; i4++){
	    if(i4==i1) continue;
	    if(i4==i2) continue;
	    if(i4==i3) continue;
	    for(int i5=0; i5<maX; i5++){
	      if(i5==i1) continue;
	      if(i5==i2) continue;
	      if(i5==i3) continue;
	      if(i5==i4) continue;
	      for(int i6=0; i6<maX; i6++){
		if(i6==i1) continue;
		if(i6==i2) continue;
		if(i6==i3) continue;
		if(i6==i4) continue;
		if(i6==i5) continue;
		for(int i7=0; i7<maX; i7++){
		  if(i7==i1) continue;
		  if(i7==i2) continue;
		  if(i7==i3) continue;
		  if(i7==i4) continue;
		  if(i7==i5) continue;
		  if(i7==i6) continue;
		  for(int i8=0; i8<maX; i8++){
		  if(i8==i1) continue;
		  if(i8==i2) continue;
		  if(i8==i3) continue;
		  if(i8==i4) continue;
		  if(i8==i5) continue;
		  if(i8==i6) continue;
		  if(i8==i7) continue;
		  for(int i9=0; i9<maX; i9++){
		    if(i9==i1) continue;
		    if(i9==i2) continue;
		    if(i9==i3) continue;
		    if(i9==i4) continue;
		    if(i9==i5) continue;
		    if(i9==i6) continue;
		    if(i9==i7) continue;
		    if(i9==i8) continue;
		    for(int i10=0; i10<maX; i10++){
		    if(i10==i1) continue;
		    if(i10==i2) continue;
		    if(i10==i3) continue;
		    if(i10==i4) continue;
		    if(i10==i5) continue;
		    if(i10==i6) continue;
		    if(i10==i7) continue;
		    if(i10==i8) continue;
		    if(i10==i9) continue;
		    for(int i11=0; i11<maX; i11++){
		      if(i11==i1) continue;
		      if(i11==i2) continue;
		      if(i11==i3) continue;
		      if(i11==i4) continue;
		      if(i11==i5) continue;
		      if(i11==i6) continue;
		      if(i11==i7) continue;
		      if(i11==i8) continue;
		      if(i11==i9) continue;
		      if(i11==i10) continue;
		      for(int i12=0; i12<maX; i12++){
			if(i12==i1) continue;
			if(i12==i2) continue;
			if(i12==i3) continue;
			if(i12==i4) continue;
			if(i12==i5) continue;
			if(i12==i6) continue;
			if(i12==i7) continue;
			if(i12==i8) continue;
			if(i12==i9) continue;
			if(i12==i10) continue;
			if(i12==i11) continue;
			for(int i13=0; i13<maX; i13++){
			  if(i13==i1) continue;
			  if(i13==i2) continue;
			  if(i13==i3) continue;
			  if(i13==i4) continue;
			  if(i13==i5) continue;
			  if(i13==i6) continue;
			  if(i13==i7) continue;
			  if(i13==i8) continue;
			  if(i13==i9) continue;
			  if(i13==i10) continue;
			  if(i13==i11) continue;
			  if(i13==i12) continue;
			  
			  if (cda1==c_da[i1]&&cda2==c_da[i2]&&cda3==c_da[i3]&&
			      cda4==c_da[i4]&&cda5==c_da[i5]&&cda6==c_da[i6]&&  
			      cda7==c_da[i7]&&cda8==c_da[i8]&&cda9==c_da[i9]&&
			      cda10==c_da[i10]&&cda11==c_da[i11]&&
			      cda12==c_da[i12]&&cda13==c_da[i13])
			    { flag=1; return 1; }
			}
		      }
		    }
		    }
		  }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    if(flag !=1){
      return 0;}
  }
  


  // Vish: Motheid class defined from here

int motherid::oiD(Particle &pa) {
    Gen_hepevt_Manager& hepevt_mag = Gen_hepevt_Manager::get_manager();
    testid =  pa.relation().genHepevt().get_ID();
    if (abs(testid) != 0) {

      return pa.relation().genHepevt().idhep(); }
    else return 9876789; // default if we don't get anything
}


  int motherid::MiD( Particle &pa, int gen) {
    Gen_hepevt_Manager& hepevt_mag = Gen_hepevt_Manager::get_manager();

       testid = pa.relation().genHepevt().get_ID();
      if (abs(testid) !=0) {        ma_id = pa.relation().genHepevt().idhep();
       } else { ma_id = 9876789; }
      for ( int p=1; p<= gen; p++) {
    if ( abs(testid) !=0 && abs(ma_id) !=300553 && abs(ma_id) != 10022  ) {
// test for mother, if we have not reached the Y(4S) or the virtual photon

   for(std::vector<Gen_hepevt>::iterator it = hepevt_mag.begin();
        it != hepevt_mag.end(); ++it){
      Gen_hepevt& particle = *it;
      if (particle.get_ID()== testid) {
        ma_id = particle.mother().idhep();
        testid =particle.mother().get_ID();
      }
    }
    }
      }
      return ma_id;
  }


// CHECK PI0
// return 2; if both gammas are from pi0
// return 1; if one gamma is correct
// return 0; if both gamma are not of pi0

 int motherid::pi0trumaa(Particle &PI0){
    Mdst_sim_ecl_Manager & sime = Mdst_sim_ecl_Manager::get_manager();
    int ma1 =0; int ma2 =0;
    int id1 =0; int id2 =0;
    int eclid1 = PI0.mdstPi0().gamma(0).ecl().get_ID();
    int eclid2 = PI0.mdstPi0().gamma(1).ecl().get_ID();
    for(std::vector<Mdst_sim_ecl>::iterator i=sime.begin();
        i!=sime.end(); i++){
      Mdst_sim_ecl& simec = *i;
      if( simec.ecl().get_ID()== eclid1) {id1= simec.hepevt().get_ID();
        if (id1 !=0 && simec.hepevt().mother().get_ID()!=0) {
          ma1= simec.hepevt().mother().idhep();} }
      if( simec.ecl().get_ID()== eclid2) {id2= simec.hepevt().get_ID();
        if (id2 !=0&& simec.hepevt().mother().get_ID()!=0) {
          ma2= simec.hepevt().mother().idhep();} }
    }

    if( (abs(ma1) ==111  && abs(ma2)!=111 ) || (abs(ma2) ==111 && abs(ma1)
!=111))
      {return 1; }
    else if(abs(ma1)==111 && abs(ma2)==111) {return 2; }
    else {return 0; }
  }


//mother of the gamma daughter of pi0
//Particle Pi0
//  nkid is kid 0, 1
//gen is the generation of mother
// 1 mother
// 2 mother
//3 mother
// you can go till infnite mother as permitted or code allow
// make good logic inorder to avoid Segmentation violation

 int motherid::pi0Gmaa(Particle &PI0, int nkid, int gen){
    Mdst_sim_ecl_Manager & sime = Mdst_sim_ecl_Manager::get_manager();
    Gen_hepevt_Manager& hepevt_mag = Gen_hepevt_Manager::get_manager();
    int eclid = PI0.mdstPi0().gamma(nkid).ecl().get_ID();
    int maaid=0; int maa =0;
    for(std::vector<Mdst_sim_ecl>::iterator i=sime.begin();
        i!=sime.end(); i++){
      Mdst_sim_ecl& simec = *i;
      if( simec.ecl().get_ID()== eclid) {
        maaid = simec.hepevt().mother().get_ID();
        if (maaid !=0) { maa = simec.hepevt().mother().idhep(); }
      }
    }

    if(maaid ==0) return 9876789;

    if (gen >1) {
      for ( int p=0; p< gen-1; p++) {
        if ( abs(maaid) !=0 && abs(maa) !=300553 && abs(maa) != 10022) {
          for(std::vector<Gen_hepevt>::iterator it = hepevt_mag.begin();
              it != hepevt_mag.end(); ++it){
            Gen_hepevt& particle = *it;
            if (particle.get_ID()== maaid) {
              maa = particle.mother().idhep();
              maaid = particle.mother().get_ID();
              //      for(int ss=0; ss<gen;  ss++){  // for fun
              //                std::cout << "::::::::"; }
              //  std::cout << endl;
            }
          }
        } else{ return 9876789;}      }
    }

    return maa;
  }



/// To get mother of Gamma using ECLID
//USAGE : eclid is
//      int eclid = gam.ecl().get_ID();
//generation of the mother

 int motherid::GMiD(int eclid, int gen){
    Mdst_sim_ecl_Manager & sime = Mdst_sim_ecl_Manager::get_manager();
    Gen_hepevt_Manager& hepevt_mag = Gen_hepevt_Manager::get_manager();
    int maaid=0; int maa =0;
    for(std::vector<Mdst_sim_ecl>::iterator i=sime.begin();
        i!=sime.end(); i++){
      Mdst_sim_ecl& simec = *i;
      if( simec.ecl().get_ID()== eclid) {
        maaid = simec.hepevt().mother().get_ID();
        if (maaid !=0) { maa = simec.hepevt().mother().idhep(); }
      }
    }

    if(maaid ==0) return 9876789;

    if (gen >1) {
      for ( int p=0; p< gen-1; p++) {
        if ( abs(maaid) !=0 && abs(maa) !=300553 && abs(maa) != 10022) {
          for(std::vector<Gen_hepevt>::iterator it = hepevt_mag.begin();
              it != hepevt_mag.end(); ++it){
            Gen_hepevt& particle = *it;
            if (particle.get_ID()== maaid) {
              maa = particle.mother().idhep();
              maaid = particle.mother().get_ID();
              //      for(int ss=0; ss<gen;  ss++){  // for fun
              //                std::cout << "::::::::"; }
              //  std::cout << endl;
            }
          }
        } else{ return 9876789;}      }
    }

    return maa;
  }


  int motherid::p0Info(Particle &D0){

    Particle pi0_1 = D0.child(0);
    Particle pi0_2 = D0.child(1);


    Mdst_sim_ecl_Manager & sime = Mdst_sim_ecl_Manager::get_manager();
    int ma1_1 =0; int ma1_2 =0;
    int id1_1 =0; int id1_2 =0;
    int eclid1_1 = pi0_1.mdstPi0().gamma(0).ecl().get_ID();
    int eclid1_2 = pi0_1.mdstPi0().gamma(1).ecl().get_ID();
    int ma2_1 =0; int ma2_2 =0;
    int id2_1 =0; int id2_2 =0;
    int eclid2_1 = pi0_2.mdstPi0().gamma(0).ecl().get_ID();
    int eclid2_2 = pi0_2.mdstPi0().gamma(1).ecl().get_ID();
    int maid1_1, maid1_2, maid2_1, maid2_2;
    maid1_1=maid1_2=maid2_1=maid2_2=0;




    for(std::vector<Mdst_sim_ecl>::iterator i=sime.begin();
        i!=sime.end(); i++){
      Mdst_sim_ecl& simec = *i;
      if( simec.ecl().get_ID()== eclid1_1) {id1_1= simec.hepevt().get_ID();
        if (id1_1 !=0 && simec.hepevt().mother().get_ID()!=0) {
          ma1_1= simec.hepevt().mother().idhep();
          maid1_1= simec.hepevt().mother().get_ID();
	} }
      if( simec.ecl().get_ID()== eclid1_2) {id1_2= simec.hepevt().get_ID();
	if (id1_2 !=0&& simec.hepevt().mother().get_ID()!=0) {
          ma1_2= simec.hepevt().mother().idhep();
          maid1_2= simec.hepevt().mother().get_ID();
	} }
    }

    for(std::vector<Mdst_sim_ecl>::iterator i=sime.begin();
        i!=sime.end(); i++){
      Mdst_sim_ecl& simec = *i;
      if( simec.ecl().get_ID()== eclid2_1) {id2_1= simec.hepevt().get_ID();
        if (id2_1 !=0 && simec.hepevt().mother().get_ID()!=0) {
          ma2_1= simec.hepevt().mother().idhep();
	  maid2_1= simec.hepevt().mother().get_ID();
	} }
      if( simec.ecl().get_ID()== eclid2_2) {id2_2= simec.hepevt().get_ID();
        if (id2_2 !=0&& simec.hepevt().mother().get_ID()!=0) {
          ma2_2= simec.hepevt().mother().idhep();
          maid2_2= simec.hepevt().mother().get_ID();
	} }
    }


    if (ma1_1 == 111 && ma1_2 ==111 && ma2_1 ==111 && ma2_2 ==111){
          
	/*  std::cout << maid1_1 << " : " << maid1_2 
		<< " : " 
		<< maid2_1 << " : " << maid2_2 
		<< " Why " << std::endl; 
      */

      if (maid1_1 == maid1_2 && maid2_1 == maid2_2){
	return 41111;}
      else if (maid1_1 == maid2_2 && maid1_2 == maid2_1) {
	return 43131; }
      else if (maid1_1 == maid2_1 && maid1_2 == maid2_2) {
	return 43113; }
      //else if (maid1_1 == maid2_2 && maid1_2 == maid2_1) {
	//return 43131; }
      else  if (maid1_1 == maid1_2 && maid2_1 != maid2_2){
        return 41133;}
      else  if (maid1_1 != maid1_2 && maid2_1 == maid2_2){
        return 43311;}
      else  if (maid1_1 != maid1_2 && maid2_1 != maid2_2){
        return 43333;}
      else { return -49999;}
    }
    else if (ma1_1 == 111 && ma1_2 !=111 && ma2_1==111 && ma2_2 == 111){
      return 31011;}
    else if (ma1_1 != 111 && ma1_2 ==111 && ma2_1==111 && ma2_2 == 111){
      return 30111;}
    else if (ma1_1 == 111 && ma1_2 ==111 && ma2_1==111 && ma2_2 != 111){
      return 31110;}
    else if (ma1_1 == 111 && ma1_2 ==111 && ma2_1!=111 && ma2_2 == 111){
      return 31101;}
    else if (ma1_1 == 111 && ma1_2 ==111 && ma2_1==111 && ma2_2 != 111){
      return 31110;}
    else if (ma1_1 == 111 && ma1_2 ==111 && ma2_1!=111 && ma2_2 != 111){
      return 21100;}
    else if (ma1_1 != 111 && ma1_2 !=111 && ma2_1==111 && ma2_2 == 111){
      return 20011;}
    else if (ma1_1 == 111 && ma1_2 !=111 && ma2_1==111 && ma2_2 != 111){
      return 21010;}
    else if (ma1_1 == 111 && ma1_2 !=111 && ma2_1!=111 && ma2_2 == 111){
      return 21001;}
    else if (ma1_1 != 111 && ma1_2 ==111 && ma2_1!=111 && ma2_2 == 111){
      return 20101;}
    else if (ma1_1 != 111 && ma1_2 ==111 && ma2_1==111 && ma2_2 != 111){
      return 20110;}
    else if (ma1_1 == 111 && ma1_2 !=111 && ma2_1!=111 && ma2_2 != 111){
      return 11000;}
    else if (ma1_1 != 111 && ma1_2 ==111 && ma2_1!=111 && ma2_2 != 111){
      return 10100;}
    else if (ma1_1 != 111 && ma1_2 !=111 && ma2_1==111 && ma2_2 != 111){
      return 10010;}
    else if (ma1_1 != 111 && ma1_2 !=111 && ma2_1!=111 && ma2_2 == 111){
      return 10001;}
    else if (ma1_1 != 111 && ma1_2 !=111 && ma2_1!=111 && ma2_2 != 111){
      return 90000;}
    return -99999;

  }

    int  motherid::Kstrumaa(Particle &Ks){
    Mdst_sim_trk_Manager & sime = Mdst_sim_trk_Manager::get_manager();
    Gen_hepevt_Manager& hepevt_mag = Gen_hepevt_Manager::get_manager();
    int Ksid1_1 = Ks.mdstVee2().chgd(0).trk().get_ID();
    int Ksid1_2 = Ks.mdstVee2().chgd(1).trk().get_ID();

    int maaid1_1=-98789; int maa1_1=-98789;
    int maaid1_2=-98789; int maa1_2=-98789;

    for(std::vector<Mdst_sim_trk>::iterator i=sime.begin();
        i!=sime.end(); i++){
      Mdst_sim_trk& simec = *i;
      if( simec.trk().get_ID()== Ksid1_1) {
        maaid1_1 = simec.hepevt().mother().get_ID();
        if (maaid1_1 !=0) { maa1_1 = simec.hepevt().mother().idhep(); }
      }
      if( simec.trk().get_ID()== Ksid1_2) {
        maaid1_2 = simec.hepevt().mother().get_ID();
        if (maaid1_2 !=0) { maa1_2 = simec.hepevt().mother().idhep(); }
      }
    }

    if (maa1_1==310 && maa1_2==310){
      if (maaid1_1==maaid1_2) return 2011;
      else return 1011;
      }

    else if (maa1_1==310 && maa1_2!=310){
      return 1010; // First one from Ks,
     }	
    else if (maa1_1!=310 && maa1_2==310){
      return 1001;// Secondone from Ks
    }
    else return -9000; // Rest of case, no one seems from Ks

  }

  //check swappig of pion in D0
  // X00YY  X is dummy and Y tell whcih pion true
  // 90099 nothing is true.
  // 90001 if one pion true
  // 90002 if two pion true of same pion2
  // -90002 swapped
  // 90020  if two pion true of same pion1
  // -90020 swapped
  // 90011  if two pion true but of different
  // -90011 swapped
  // 90021  if three pion true
  //-90021 swapped
  // 90012  if three pion true
  // -90012 swapped
  // 90022  all four pion from same 
  //-92022 pion swapped
  
  
  int motherid::D0Info(Particle &D0){
    
    Particle K_S1 = D0.child(0);
    Particle K_S2 = D0.child(1);
    Mdst_sim_trk_Manager & sime = Mdst_sim_trk_Manager::get_manager();
    
    Gen_hepevt_Manager& hepevt_mag = Gen_hepevt_Manager::get_manager();
    int Ksid1_1 = K_S1.mdstVee2().chgd(0).trk().get_ID();
    int Ksid1_2 = K_S1.mdstVee2().chgd(1).trk().get_ID();
    int Ksid2_1 = K_S2.mdstVee2().chgd(0).trk().get_ID();
    int Ksid2_2 = K_S2.mdstVee2().chgd(1).trk().get_ID();

    
    // 98789 magic number with negative
    int maaid1_1=-98789; int maa1_1=-98789;
    int maaid1_2=-98789; int maa1_2=-98789;
    int maaid2_1=-98789; int maa2_1=-98789;
    int maaid2_2=-98789; int maa2_2=-98789;
    
    
    for(std::vector<Mdst_sim_trk>::iterator i=sime.begin();
	i!=sime.end(); i++){
      Mdst_sim_trk& simec = *i;
      if( simec.trk().get_ID()== Ksid1_1) {
	maaid1_1 = simec.hepevt().mother().get_ID();
	if (maaid1_1 !=0) { maa1_1 = simec.hepevt().mother().idhep(); }
      }
      if( simec.trk().get_ID()== Ksid1_2) {
 	maaid1_2 = simec.hepevt().mother().get_ID();
	if (maaid1_2 !=0) { maa1_2 = simec.hepevt().mother().idhep(); }
      }
      if( simec.trk().get_ID()== Ksid2_1) {
	maaid2_1 = simec.hepevt().mother().get_ID();
	if (maaid2_1 !=0) { maa2_1 = simec.hepevt().mother().idhep(); }
      }
      if( simec.trk().get_ID()== Ksid2_2) {
	maaid2_2 = simec.hepevt().mother().get_ID();
	if (maaid2_2 !=0) { maa2_2 = simec.hepevt().mother().idhep(); }
      }
    }
    /*
      std::cout << "Start of PUPPY \n" << 
      maa1_1 << "\t" << maaid1_1 << "\n" <<
      maa1_2 << "\t" << maaid1_2 << "\n" <<
      maa2_1 << "\t" << maaid2_1 << "\n" <<
      maa2_2 << "\t" << maaid2_2 << "\n" << 
      " Bye Bye PUPPY \n";
    */
    
    if (maa1_1==310 && maa1_2==310 &&
	maa2_1==310 && maa2_2==310){
      if (maaid1_1==maaid1_2 && maaid2_1==maaid2_2) return 90022;
      else if (maaid1_1==maaid2_2 && maaid1_2==maaid2_1) return -92022;
      else if (maaid1_1==maaid2_1 && maaid1_2==maaid2_2) return -92022;
    }

    // CHeck for 90021
    else if (maa1_1==310 && maa1_2==310 &&
	     maa2_1==310 && maa2_2!=310){
      if (maaid1_1==maaid1_2) return 90021;
      else if(maaid1_1==maaid2_1) return -90021;
      else if(maaid1_2==maaid2_1) return -90021;
    }
    else if (maa1_1==310 && maa1_2==310 &&
	     maa2_1!=310 && maa2_2==310){
      if (maaid1_1==maaid2_2) return 90021;
      else if(maaid1_1==maaid2_2) return -90021;
      else if(maaid2_2==maaid2_1) return -90021;
    }
    else if (maa1_1==310 && maa1_2!=310 &&
	     maa2_1==310 && maa2_2==310){
      if (maaid2_1==maaid2_2) return 90012;
      else if(maaid2_1==maaid1_1) return -90012;
      else if(maaid2_2==maaid1_1) return -90012;
    }
    else if (maa1_1!=310 && maa1_2==310 &&
	     maa2_1==310 && maa2_2==310){
      if (maaid2_1==maaid2_2) return 90012;
      else if(maaid2_1==maaid1_2) return -90012;
      else if(maaid2_2==maaid1_2) return -90012;
    }
    
    else if (maa1_1!=310 && maa1_2==310 &&
	     maa2_1==310 && maa2_2!=310){
      if (maaid1_2==maaid2_1) return -92011;
      else return 90011;}
    else if (maa1_1!=310 && maa1_2==310 &&
	     maa2_1!=310 && maa2_2==310){
      if (maaid1_2==maaid2_2) return -92011;
      else return 90011;}
    else if (maa1_1==310 && maa1_2!=310 &&
	     maa2_1!=310 && maa2_2==310){
      if (maaid1_1==maaid2_2) return -92011;
      else return 90011;
    }
    else if(maa1_1==310 && maa1_2!=310 &&
	    maa2_1==310 && maa2_2!=310){
      if (maaid1_1==maaid2_1) return -92011;
      else return 90011;
    }
    
    else if(maa1_1==310&&maa1_2==310&&
	    maa2_1!=310 && maa2_2!=310){
      if (maaid1_1==maaid1_2) return 90020;
      else return -90020;
    }
     
    else if(maa1_1!=310&&maa1_2!=310&&
	    maa2_1==310 && maa2_2==310){
      if (maaid1_1==maaid1_2) return 90002;
      else return -90002;
    }
    else if( (maa1_1==310&&maa1_2!=310&&
	      maa2_1!=310&&maa2_2!=310) ||
	     (maa1_1!=310&&maa1_2==310&&
	      maa2_1!=310&&maa2_2!=310) ||
	     (maa1_1!=310&&maa1_2!=310&&
	      maa2_1==310&&maa2_2!=310) ||
	     (maa1_1!=310&&maa1_2!=310&&
	      maa2_1!=310&&maa2_2==310)) return 90001;
    else 90099;
  }




//mother of the pi+ of KS
//Particle Ks0
//  nkid is kid 0, 1
//gen is the generation of mother
// 1 mother
// 2 mother
// 3 mother
// you can go till infnite mother as permitted or code allow
// make good logic inorder to avoid Segmentation violation
  
  int motherid::KsPimaa(Particle &K_S, int nkid, int gen){
    Mdst_sim_trk_Manager & sime = Mdst_sim_trk_Manager::get_manager();
    Gen_hepevt_Manager& hepevt_mag = Gen_hepevt_Manager::get_manager();
    int eclid = K_S.mdstVee2().chgd(nkid).trk().get_ID();
    int maaid=0; int maa =0;
    for(std::vector<Mdst_sim_trk>::iterator i=sime.begin();
        i!=sime.end(); i++){
      Mdst_sim_trk& simec = *i;
      if( simec.trk().get_ID()== eclid) {
        maaid = simec.hepevt().mother().get_ID();
        if (maaid !=0) { maa = simec.hepevt().mother().idhep(); }
      }
    }
    
    if(maaid ==0) return 9876789;
    
    if (gen >1) {
      for ( int p=0; p< gen-1; p++) {
        if ( abs(maaid) !=0 && abs(maa) !=300553 && abs(maa) != 10022) {
          for(std::vector<Gen_hepevt>::iterator it = hepevt_mag.begin();
              it != hepevt_mag.end(); ++it){
            Gen_hepevt& particle = *it;
            if (particle.get_ID()== maaid) {
              maa = particle.mother().idhep();
              maaid = particle.mother().get_ID();
              //      for(int ss=0; ss<gen;  ss++){  // for fun
              //                std::cout << "::::::::"; }
              //  std::cout << endl;
            }
          }
        } else{ return 9876789;}      }
    }
    
    return maa;
  }
  





#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
