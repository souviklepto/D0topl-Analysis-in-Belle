

//___________________________________________________________________
// genT class to get generated information
// This class is useful in making the code very small for large combinations
// Also, you can Print the daughters of the particle.
// 
//   
// AUTHOR: Vishal  Place: Nara W. University, U. South Carolina
//                  
// Version 0.10 : bug fixed for repetition
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
// :o: Vec4DAU(Gen_hepevt a, int id)
//      Returns HepLorentzVector of  position of child of a having lund id
// :o: Vec4DAU(Gen_hepevt a, int id, int rank)
//      Returns HepLorentzVector position of child of a having lund id
//      rank take care of repetition
// :o: AbsVec4DAU(Gen_hepevt a, int id)
//      Charge insensitive,Returns HepLorentzVectorPosition of child of a having lund id
// :o: AbsVec4DAU(Gen_hepevt a, int id, int rank)
//      Charge insenitive,Returns HepLorentzVectorPosition of child of a having lund id
//      rank take care of repetition//
// :o: Vec(Gen_hepevt a);
//     Return HepLorentzVector ofPosition of a
// known bug for Bmode ==> 209, 387, 465,650,1008,1021,1042,1388
//
//___________________________________________________________________

#include "belle.h"
#include "get_gen_tag.h"
#include HEPEVT_H
#include BELLETDF_H
#include MDST_H

#include "particle/Ptype.h"
#include <math.h>

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif


  //Mothers


  int genT::oiD(Particle &pa) {
    Gen_hepevt_Manager& hepevt_mag = Gen_hepevt_Manager::get_manager();
    testid =  pa.relation().genHepevt().get_ID();
    if (abs(testid) != 0) {
      
      return pa.relation().genHepevt().idhep(); }
    else return 9876789; // default if we don't get anything
  }
  


  int genT::MiD( Particle &pa, int gen) {
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
  

  int genT::MiDNum( Particle &pa, int gen) {
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
    return testid;
  }
  



  // CHECK PI0
  // return 2; if both gammas are from pi0
  // return 1; if one gamma is correct
  // return 0; if both gamma are not of pi0

  int genT::pi0trumaa(Particle &PI0){
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

int genT::pi0Gmaa(Particle &PI0, int nkid, int gen){
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
  
  int genT::GMiD(int eclid, int gen){
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
  
  


  //mother of the pi+ of KS                                                                                                                                      
  //Particle Ks0                                                                                                                                                 
  //  nkid is kid 0, 1                                                                                                                                           
  //gen is the generation of mother                                                                                                                              
  // 1 mother                                                                                                                                                    
  // 2 mother                                                                                                                                                    
  // 3 mother                                                                                                                                                    
  // you can go till infnite mother as permitted or code allow                                                                                                   
  // make good logic inorder to avoid Segmentation violation                                                                                                     

  int genT::KsPimaa(Particle &K_S, int nkid, int gen){
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
  


  int genT::Gisthep( int eclid){
    Mdst_ecl_aux_Manager &EcL_AuX = Mdst_ecl_aux_Manager::get_manager();
    Mdst_sim_ecl_Manager & sime = Mdst_sim_ecl_Manager::get_manager();
    Gen_hepevt_Manager& hepevt_mag = Gen_hepevt_Manager::get_manager();
    int ownid=0; int isthep =0;
    for(std::vector<Mdst_sim_ecl>::iterator i=sime.begin();
        i!=sime.end(); i++){
      Mdst_sim_ecl& simec = *i;
      if( simec.ecl().get_ID()== eclid) {
        ownid = simec.hepevt().get_ID();
        if (ownid !=0) { isthep = simec.hepevt().isthep(); }
      }
    }
    return isthep;
  }


  int genT::GTinfo( int eclid, int noTri) {
    Mdst_ecl_aux_Manager & EcL_AuX = Mdst_ecl_aux_Manager::get_manager();

    double nhit  = EcL_AuX[eclid -1].property(6);

    float propy1 = EcL_AuX[eclid-1].property(1);
    float propy2 = EcL_AuX[eclid-1].property(2);
    float propy3 = EcL_AuX[eclid-1].property(3);

    int pop1 = *(int*)&propy1;
    int pop2 = *(int*)&propy2;
    int pop3 = *(int*)&propy3;

    int tcid1 = (pop1 & 0x3ff);
    int tsid1 = (pop1 >> 10) & 0x1;
    int tcid2 = (pop1 >> 11) & 0x3ff;
    int tsid2 = (pop1 >> 21) & 0x1;
    int tdccount1 = pop2 & 0xffff;
    int elinkcls1 = (pop2 >> 16) & 0x7f;
    int elinktc1 = (pop2 >> 23) & 0x7f;
    int hecls1 = (pop2 >> 30) & 0x1;
    int hetc1 = (pop2 >> 31) & 0x1;

    int tdccount2 = pop3 & 0xffff;
    int elinkcls2 = (pop3 >> 16) & 0x7f;
    int elinktc2 = (pop3 >> 23) & 0x7f;
    int hecls2 = (pop3 >> 30) & 0x1;
    int hetc2 = (pop3 >> 31) & 0x1;

    if(noTri==1){
      return tdccount1;}
    else if(noTri ==2){
      return tdccount2;}
    else{
      std::cout << "FATAL ERROR MAA says\t PLEASE PUT CORRECT WHICH TRIGGER you want\n correct input 1 or 2 " << std::endl;
      return 0; }
    // if( tcid1!=0 && tcid1==tcid2 && tsid1==tsid2){                                                                                                          

    //return nhit;                                                                                                                                             
  }
  
  
  //
  



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

  // get 4 Vector of X,Y,Z of vertex and production time in mm/c
  HepLorentzVector genT::Vec4(Gen_hepevt a){
    HepLorentzVector p(a.VX(), a.VY(), a.VZ(),a.T());
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

  /***************/


  // Returns the HepLoretnzVector of Vertex daughter after matching with its id,
  // WARNING :please use it only after checking that there is a matching 
  // daughter otherwise it will give segmentation violation :(
  HepLorentzVector genT::Vec4DAU(Gen_hepevt a, int id, int rank){
    int i=1;
    for(int j=0; j<dau(a); j++){
      if(MyChild(a,j+1).idhep()==id){if (rank==i){return Vec4(MyChild(a,j+1));}
	else { i=i+1;} }
    }
  }
  
  //DAME without need of rank.
  // This function lives on the previous only the rank is put 1 here.
  HepLorentzVector genT::Vec4DAU(Gen_hepevt a, int id){
    return Vec4DAU(a,id,1);
  }
  
  // Abs doesn't care about the sign..
  // ...
  HepLorentzVector genT::AbsVec4DAU(Gen_hepevt a, int id, int rank){
    int i=1;
    for(int j=0; j<dau(a); j++){
      if(abs(MyChild(a,j+1).idhep())==abs(id)){if (rank==i){
	  return Vec4(MyChild(a,j+1));}
	else { i=i+1;} }
    }
  }
  
  //DAME without need of rank.
  // This function lives on the previous only the rank is put 1 here.
  HepLorentzVector genT::AbsVec4DAU(Gen_hepevt a, int id){
    return AbsVec4DAU(a,abs(id),1);
  }

  /***/


 


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

  int genT::Dmode(Gen_hepevt genpart){
    if(AmId(genpart)==-421)
{
  if ( PcheckDecay(genpart,321, 11, -12)==1){
    return 1;}//anti-D0 decays to K+ e- anti-nu_e 
  
  if ( PcheckDecay(genpart,323, 11, -12)==1){
    return 2;}//anti-D0 decays to K*+ e- anti-nu_e 
  
  if ( PcheckDecay(genpart,325, 11, -12)==1){
    return 3;}//anti-D0 decays to K_2*+ e- anti-nu_e 
  
  if ( PcheckDecay(genpart,211, 11, -12)==1){
    return 4;}//anti-D0 decays to pi+ e- anti-nu_e 
  
  if ( PcheckDecay(genpart,213, 11, -12)==1){
    return 5;}//anti-D0 decays to rho+ e- anti-nu_e 
  
  if ( PcheckDecay(genpart,10213, 11, -12)==1){
    return 6;}//anti-D0 decays to b_1+ e- anti-nu_e 
  
  if ( PcheckDecay(genpart,321, 111, 11, -12)==1){
    return 7;}//anti-D0 decays to K+ pi0 e- anti-nu_e 
  
  if ( PcheckDecay(genpart,311, 211, 11, -12)==1){
    return 8;}//anti-D0 decays to K0 pi+ e- anti-nu_e 
  
  //if ( PcheckDecay(genpart,return 9;}//anti-D0 decays to 
  
  if ( PcheckDecay(genpart,321, 13, -14)==1){
    return 10;}//anti-D0 decays to K+ mu- anti-nu_mu 
  
  if ( PcheckDecay(genpart,323, 13, -14)==1){
    return 11;}//anti-D0 decays to K*+ mu- anti-nu_mu 
  
  if ( PcheckDecay(genpart,325, 13, -14)==1){
    return 12;}//anti-D0 decays to K_2*+ mu- anti-nu_mu 
  
  if ( PcheckDecay(genpart,211, 13, -14)==1){
    return 13;}//anti-D0 decays to pi+ mu- anti-nu_mu 
  
  if ( PcheckDecay(genpart,213, 13, -14)==1){
    return 14;}//anti-D0 decays to rho+ mu- anti-nu_mu 
  
  if ( PcheckDecay(genpart,10213, 13, -14)==1){
    return 15;}//anti-D0 decays to b_1+ mu- anti-nu_mu 
  
  if ( PcheckDecay(genpart,321, 111, 13, -14)==1){
    return 16;}//anti-D0 decays to K+ pi0 mu- anti-nu_mu 

if ( PcheckDecay(genpart,311, 211, 13, -14)==1){
return 17;}//anti-D0 decays to K0 pi+ mu- anti-nu_mu 

if ( PcheckDecay(genpart,321, -211)==1){
return 18;}//anti-D0 decays to K+ pi- 

if ( PcheckDecay(genpart,311, 111)==1){
return 19;}//anti-D0 decays to K0 pi0 

if ( PcheckDecay(genpart,311, 221)==1){
return 20;}//anti-D0 decays to K0 eta 

if ( PcheckDecay(genpart,311, 331)==1){
return 21;}//anti-D0 decays to K0 eta' 

if ( PcheckDecay(genpart,113, 311)==1){
return 22;}//anti-D0 decays to rho0 K0 

if ( PcheckDecay(genpart,-213, 321)==1){
return 23;}//anti-D0 decays to rho- K+ 

if ( PcheckDecay(genpart,223, 311)==1){
return 24;}//anti-D0 decays to omega K0 

if ( PcheckDecay(genpart,313, 221)==1){
return 25;}//anti-D0 decays to K*0 eta 

if ( PcheckDecay(genpart,-20213, 321)==1){
return 26;}//anti-D0 decays to a_1- K+ 

if ( PcheckDecay(genpart,323, -213)==1){
return 27;}//anti-D0 decays to K*+ rho- 

if ( PcheckDecay(genpart,313, 113)==1){
return 28;}//anti-D0 decays to K*0 rho0 

if ( PcheckDecay(genpart,313, 223)==1){
return 29;}//anti-D0 decays to K*0 omega 

if ( PcheckDecay(genpart,313, 111)==1){
return 30;}//anti-D0 decays to K*0 pi0 

if ( PcheckDecay(genpart,323, -211)==1){
return 31;}//anti-D0 decays to K*+ pi- 

if ( PcheckDecay(genpart,10323, -211)==1){
return 32;}//anti-D0 decays to K_1+ pi- 

if ( PcheckDecay(genpart,10221, 311)==1){
return 33;}//anti-D0 decays to f_0 K0 

if ( PcheckDecay(genpart,10331, 311)==1){
return 34;}//anti-D0 decays to f'_0 K0 

if ( PcheckDecay(genpart,10321, -211)==1){
return 35;}//anti-D0 decays to K_0*+ pi- 

if ( PcheckDecay(genpart,325, -211)==1){
return 36;}//anti-D0 decays to K_2*+ pi- 

if ( PcheckDecay(genpart,10311, 111)==1){
return 37;}//anti-D0 decays to K_0*0 pi0 

if ( PcheckDecay(genpart,321, -211, 111)==1){
return 38;}//anti-D0 decays to K+ pi- pi0 

if ( PcheckDecay(genpart,311, 111, 111)==1){
return 39;}//anti-D0 decays to K0 pi0 pi0 

if ( PcheckDecay(genpart,313, 211, -211)==1){
return 40;}//anti-D0 decays to K*0 pi+ pi- 

if ( PcheckDecay(genpart,313, 111, 111)==1){
return 41;}//anti-D0 decays to K*0 pi0 pi0 

if ( PcheckDecay(genpart,321, -213, 111)==1){
return 42;}//anti-D0 decays to K+ rho- pi0 

if ( PcheckDecay(genpart,321, -211, 113)==1){
return 43;}//anti-D0 decays to K+ pi- rho0 

if ( PcheckDecay(genpart,321, -211, 223)==1){
return 44;}//anti-D0 decays to K+ pi- omega 

if ( PcheckDecay(genpart,321, -211, 221)==1){
return 45;}//anti-D0 decays to K+ pi- eta 

if ( PcheckDecay(genpart,321, -211, 331)==1){
return 46;}//anti-D0 decays to K+ pi- eta' 

if ( PcheckDecay(genpart,313, -211, 211, 111)==1){
return 47;}//anti-D0 decays to K*0 pi- pi+ pi0 

if ( PcheckDecay(genpart,321, -211, 211, -211)==1){
return 48;}//anti-D0 decays to K+ pi- pi+ pi- 

if ( PcheckDecay(genpart,311, 211, -211, 111)==1){
return 49;}//anti-D0 decays to K0 pi+ pi- pi0 

if ( PcheckDecay(genpart,321, -211, 111, 111)==1){
return 50;}//anti-D0 decays to K+ pi- pi0 pi0 

if ( PcheckDecay(genpart,311, 111, 111, 111)==1){
return 51;}//anti-D0 decays to K0 pi0 pi0 pi0 

if ( PcheckDecay(genpart,321, -211, 211, -211, 111)==1){
return 52;}//anti-D0 decays to K+ pi- pi+ pi- pi0 

if ( PcheckDecay(genpart,321, -211, 111, 111, 111)==1){
return 53;}//anti-D0 decays to K+ pi- pi0 pi0 pi0 

if ( PcheckDecay(genpart,311, -211, 211, 111, 111)==1){
return 54;}//anti-D0 decays to K0 pi- pi+ pi0 pi0 

if ( PcheckDecay(genpart,311, -211, 211, 111, 111, 111)==1){
return 55;}//anti-D0 decays to K0 pi- pi+ pi0 pi0 pi0 

if ( PcheckDecay(genpart,333, 311)==1){
return 56;}//anti-D0 decays to phi K0 

if ( PcheckDecay(genpart,333, 313)==1){
return 57;}//anti-D0 decays to phi K*0 

if ( PcheckDecay(genpart,333, 321, -211)==1){
return 58;}//anti-D0 decays to phi K+ pi- 

if ( PcheckDecay(genpart,313, 321, -321)==1){
return 59;}//anti-D0 decays to K*0 K+ K- 

if ( PcheckDecay(genpart,311, 321, -321)==1){
return 60;}//anti-D0 decays to K0 K+ K- 

if ( PcheckDecay(genpart,-321, 321, 321, -211)==1){
return 61;}//anti-D0 decays to K- K+ K+ pi- 

//if ( PcheckDecay(genpart,return 62;}//anti-D0 decays to 

if ( PcheckDecay(genpart,310, 310, 310)==1){
return 63;}//anti-D0 decays to K_S0 K_S0 K_S0 

if ( PcheckDecay(genpart,130, 130, 130)==1){
return 64;}//anti-D0 decays to K_L0 K_L0 K_L0 

if ( PcheckDecay(genpart,321, -321)==1){
return 65;}//anti-D0 decays to K+ K- 

if ( PcheckDecay(genpart,310, 310)==1){
return 66;}//anti-D0 decays to K_S0 K_S0 

if ( PcheckDecay(genpart,130, 130)==1){
return 67;}//anti-D0 decays to K_L0 K_L0 

if ( PcheckDecay(genpart,323, -321)==1){
return 68;}//anti-D0 decays to K*+ K- 

if ( PcheckDecay(genpart,-323, 321)==1){
return 69;}//anti-D0 decays to K*- K+ 

if ( PcheckDecay(genpart,313, -313)==1){
return 70;}//anti-D0 decays to K*0 anti-K*0 

if ( PcheckDecay(genpart,333, 111)==1){
return 71;}//anti-D0 decays to phi pi0 

if ( PcheckDecay(genpart,333, 221)==1){
return 72;}//anti-D0 decays to phi eta 

if ( PcheckDecay(genpart,-311, 321, -211)==1){
return 73;}//anti-D0 decays to anti-K0 K+ pi- 

if ( PcheckDecay(genpart,311, -321, 211)==1){
return 74;}//anti-D0 decays to K0 K- pi+ 

if ( PcheckDecay(genpart,333, 113)==1){
return 75;}//anti-D0 decays to phi rho0 

if ( PcheckDecay(genpart,333, 211, -211)==1){
return 76;}//anti-D0 decays to phi pi+ pi- 

if ( PcheckDecay(genpart,321, -321, 113)==1){
return 77;}//anti-D0 decays to K+ K- rho0 

if ( PcheckDecay(genpart,321, -321, 111, 111)==1){
return 78;}//anti-D0 decays to K+ K- pi0 pi0 

if ( PcheckDecay(genpart,-313, 321, -211)==1){
return 79;}//anti-D0 decays to anti-K*0 K+ pi- 

if ( PcheckDecay(genpart,313, -321, 211)==1){
return 80;}//anti-D0 decays to K*0 K- pi+ 

if ( PcheckDecay(genpart,-311, 311, 211, -211)==1){
return 81;}//anti-D0 decays to anti-K0 K0 pi+ pi- 

if ( PcheckDecay(genpart,321, -321, 321, -211)==1){
return 82;}//anti-D0 decays to K+ K- K+ pi- 

if ( PcheckDecay(genpart,321, -321, 311, 111)==1){
return 83;}//anti-D0 decays to K+ K- K0 pi0 

if ( PcheckDecay(genpart,321, -321, 211, -211, 111)==1){
return 84;}//anti-D0 decays to K+ K- pi+ pi- pi0 

if ( PcheckDecay(genpart,211, -211)==1){
return 85;}//anti-D0 decays to pi+ pi- 

if ( PcheckDecay(genpart,111, 111)==1){
return 86;}//anti-D0 decays to pi0 pi0 

if ( PcheckDecay(genpart,221, 111)==1){
return 87;}//anti-D0 decays to eta pi0 

if ( PcheckDecay(genpart,331, 111)==1){
return 88;}//anti-D0 decays to eta' pi0 

if ( PcheckDecay(genpart,221, 221)==1){
return 89;}//anti-D0 decays to eta eta 

if ( PcheckDecay(genpart,213, -211)==1){
return 90;}//anti-D0 decays to rho+ pi- 

if ( PcheckDecay(genpart,-213, 211)==1){
return 91;}//anti-D0 decays to rho- pi+ 

if ( PcheckDecay(genpart,113, 111)==1){
return 92;}//anti-D0 decays to rho0 pi0 

if ( PcheckDecay(genpart,111, 111, 111)==1){
return 93;}//anti-D0 decays to pi0 pi0 pi0 

if ( PcheckDecay(genpart,211, 211, -211, -211)==1){
return 94;}//anti-D0 decays to pi+ pi+ pi- pi- 

if ( PcheckDecay(genpart,211, -211, 111, 111)==1){
return 95;}//anti-D0 decays to pi+ pi- pi0 pi0 

if ( PcheckDecay(genpart,211, -211, 211, -211, 111)==1){
return 96;}//anti-D0 decays to pi+ pi- pi+ pi- pi0 

if ( PcheckDecay(genpart,211, -211, 111, 111, 111)==1){
return 97;}//anti-D0 decays to pi+ pi- pi0 pi0 pi0 

if ( PcheckDecay(genpart,211, -211, 211, -211, 211, -211)==1){
return 98;}//anti-D0 decays to pi+ pi- pi+ pi- pi+ pi- 

if ( PcheckDecay(genpart,211, -321)==1){
return 99;}//anti-D0 decays to pi+ K- 

if ( PcheckDecay(genpart,-323, 211)==1){
return 100;}//anti-D0 decays to K*- pi+ 

if ( PcheckDecay(genpart,211, -321, 111)==1){
return 101;}//anti-D0 decays to pi+ K- pi0 

if ( PcheckDecay(genpart,211, -321, -211, 211)==1){
return 102;}//anti-D0 decays to pi+ K- pi- pi+ 

if ( PcheckDecay(genpart,-2212, -11)==1){
   return 3805;}//anti-D0 decays to p- e+  

 return -100;

 }
    if (AmId(genpart)==421){
  
  if ( PcheckDecay(genpart,-321, -11, 12)==1){
    return 103;}//D0 decays to K- e+ nu_e 
  
  if ( PcheckDecay(genpart,-323, -11, 12)==1){
    return 104;}//D0 decays to K*- e+ nu_e 
  
  if ( PcheckDecay(genpart,-325, -11, 12)==1){
    return 105;}//D0 decays to K_2*- e+ nu_e 
  
  if ( PcheckDecay(genpart,-211, -11, 12)==1){
    return 106;}//D0 decays to pi- e+ nu_e 
  
  if ( PcheckDecay(genpart,-213, -11, 12)==1){
    return 107;}//D0 decays to rho- e+ nu_e 
  
  if ( PcheckDecay(genpart,-10213, -11, 12)==1){
    return 108;}//D0 decays to b_1- e+ nu_e 
  
  if ( PcheckDecay(genpart,-321, 111, -11, 12)==1){
    return 109;}//D0 decays to K- pi0 e+ nu_e 
  
  if ( PcheckDecay(genpart,-311, -211, -11, 12)==1){
    return 110;}//D0 decays to anti-K0 pi- e+ nu_e 
  
  //if ( PcheckDecay(genpart,return 111;}//D0 decays to 
  
  if ( PcheckDecay(genpart,-321, -13, 14)==1){
    return 112;}//D0 decays to K- mu+ nu_mu 
  
  if ( PcheckDecay(genpart,-323, -13, 14)==1){
    return 113;}//D0 decays to K*- mu+ nu_mu 
  
  if ( PcheckDecay(genpart,-325, -13, 14)==1){
    return 114;}//D0 decays to K_2*- mu+ nu_mu 
  
  if ( PcheckDecay(genpart,-211, -13, 14)==1){
    return 115;}//D0 decays to pi- mu+ nu_mu 
  
  if ( PcheckDecay(genpart,-213, -13, 14)==1){
return 116;}//D0 decays to rho- mu+ nu_mu 

if ( PcheckDecay(genpart,-10213, -13, 14)==1){
return 117;}//D0 decays to b_1- mu+ nu_mu 

if ( PcheckDecay(genpart,-321, 111, -13, 14)==1){
return 118;}//D0 decays to K- pi0 mu+ nu_mu 

if ( PcheckDecay(genpart,-311, -211, -13, 14)==1){
return 119;}//D0 decays to anti-K0 pi- mu+ nu_mu 

if ( PcheckDecay(genpart,-321, 211)==1){
return 120;}//D0 decays to K- pi+ 

if ( PcheckDecay(genpart,-311, 111)==1){
return 121;}//D0 decays to anti-K0 pi0 

if ( PcheckDecay(genpart,-311, 221)==1){
return 122;}//D0 decays to anti-K0 eta 

if ( PcheckDecay(genpart,-311, 331)==1){
return 123;}//D0 decays to anti-K0 eta' 

if ( PcheckDecay(genpart,113, -311)==1){
return 124;}//D0 decays to rho0 anti-K0 

if ( PcheckDecay(genpart,213, -321)==1){
return 125;}//D0 decays to rho+ K- 

if ( PcheckDecay(genpart,223, -311)==1){
return 126;}//D0 decays to omega anti-K0 

if ( PcheckDecay(genpart,-313, 221)==1){
return 127;}//D0 decays to anti-K*0 eta 

if ( PcheckDecay(genpart,20213, -321)==1){
return 128;}//D0 decays to a_1+ K- 

if ( PcheckDecay(genpart,-323, 213)==1){
return 129;}//D0 decays to K*- rho+ 

if ( PcheckDecay(genpart,-313, 113)==1){
return 130;}//D0 decays to anti-K*0 rho0 

if ( PcheckDecay(genpart,-313, 223)==1){
return 131;}//D0 decays to anti-K*0 omega 

if ( PcheckDecay(genpart,-313, 111)==1){
return 132;}//D0 decays to anti-K*0 pi0 

if ( PcheckDecay(genpart,-323, 211)==1){
return 133;}//D0 decays to K*- pi+ 

if ( PcheckDecay(genpart,-10323, 211)==1){
return 134;}//D0 decays to K_1- pi+ 

if ( PcheckDecay(genpart,10221, -311)==1){
return 135;}//D0 decays to f_0 anti-K0 

if ( PcheckDecay(genpart,10331, -311)==1){
return 136;}//D0 decays to f'_0 anti-K0 

if ( PcheckDecay(genpart,-10321, 211)==1){
return 137;}//D0 decays to K_0*- pi+ 

if ( PcheckDecay(genpart,-325, 211)==1){
return 138;}//D0 decays to K_2*- pi+ 

if ( PcheckDecay(genpart,-10311, 111)==1){
return 139;}//D0 decays to anti-K_0*0 pi0 

if ( PcheckDecay(genpart,-321, 211, 111)==1){
return 140;}//D0 decays to K- pi+ pi0 

if ( PcheckDecay(genpart,-311, 111, 111)==1){
return 141;}//D0 decays to anti-K0 pi0 pi0 

if ( PcheckDecay(genpart,-313, 211, -211)==1){
return 142;}//D0 decays to anti-K*0 pi+ pi- 

if ( PcheckDecay(genpart,-313, 111, 111)==1){
return 143;}//D0 decays to anti-K*0 pi0 pi0 

if ( PcheckDecay(genpart,-321, 213, 111)==1){
return 144;}//D0 decays to K- rho+ pi0 

if ( PcheckDecay(genpart,-321, 211, 113)==1){
return 145;}//D0 decays to K- pi+ rho0 

if ( PcheckDecay(genpart,-321, 211, 223)==1){
return 146;}//D0 decays to K- pi+ omega 

if ( PcheckDecay(genpart,-321, 211, 221)==1){
return 147;}//D0 decays to K- pi+ eta 

if ( PcheckDecay(genpart,-321, 211, 331)==1){
return 148;}//D0 decays to K- pi+ eta' 

if ( PcheckDecay(genpart,-313, 211, -211, 111)==1){
return 149;}//D0 decays to anti-K*0 pi+ pi- pi0 

if ( PcheckDecay(genpart,-321, 211, 211, -211)==1){
return 150;}//D0 decays to K- pi+ pi+ pi- 

if ( PcheckDecay(genpart,-311, 211, -211, 111)==1){
return 151;}//D0 decays to anti-K0 pi+ pi- pi0 

if ( PcheckDecay(genpart,-321, 211, 111, 111)==1){
return 152;}//D0 decays to K- pi+ pi0 pi0 

if ( PcheckDecay(genpart,-311, 111, 111, 111)==1){
return 153;}//D0 decays to anti-K0 pi0 pi0 pi0 

if ( PcheckDecay(genpart,-321, 211, 211, -211, 111)==1){
return 154;}//D0 decays to K- pi+ pi+ pi- pi0 

if ( PcheckDecay(genpart,-321, 211, 111, 111, 111)==1){
return 155;}//D0 decays to K- pi+ pi0 pi0 pi0 

if ( PcheckDecay(genpart,-311, 211, -211, 111, 111)==1){
return 156;}//D0 decays to anti-K0 pi+ pi- pi0 pi0 

if ( PcheckDecay(genpart,-311, 211, -211, 111, 111, 111)==1){
return 157;}//D0 decays to anti-K0 pi+ pi- pi0 pi0 pi0 

if ( PcheckDecay(genpart,333, -311)==1){
return 158;}//D0 decays to phi anti-K0 

if ( PcheckDecay(genpart,333, -313)==1){
return 159;}//D0 decays to phi anti-K*0 

if ( PcheckDecay(genpart,333, -321, 211)==1){
return 160;}//D0 decays to phi K- pi+ 

if ( PcheckDecay(genpart,-313, 321, -321)==1){
return 161;}//D0 decays to anti-K*0 K+ K- 

if ( PcheckDecay(genpart,-311, 321, -321)==1){
return 162;}//D0 decays to anti-K0 K+ K- 

if ( PcheckDecay(genpart,321, -321, -321, 211)==1){
return 163;}//D0 decays to K+ K- K- pi+ 

     //if ( PcheckDecay(genpart,return 164;}//D0 decays to 

if ( PcheckDecay(genpart,310, 310, 310)==1){
return 165;}//D0 decays to K_S0 K_S0 K_S0 

if ( PcheckDecay(genpart,130, 130, 130)==1){
return 166;}//D0 decays to K_L0 K_L0 K_L0 

if ( PcheckDecay(genpart,321, -321)==1){
return 167;}//D0 decays to K+ K- 

if ( PcheckDecay(genpart,310, 310)==1){
return 168;}//D0 decays to K_S0 K_S0 

if ( PcheckDecay(genpart,130, 130)==1){
return 169;}//D0 decays to K_L0 K_L0 

if ( PcheckDecay(genpart,-323, 321)==1){
return 170;}//D0 decays to K*- K+ 

if ( PcheckDecay(genpart,323, -321)==1){
return 171;}//D0 decays to K*+ K- 

if ( PcheckDecay(genpart,-313, 313)==1){
return 172;}//D0 decays to anti-K*0 K*0 

if ( PcheckDecay(genpart,333, 111)==1){
return 173;}//D0 decays to phi pi0 

if ( PcheckDecay(genpart,333, 221)==1){
return 174;}//D0 decays to phi eta 

if ( PcheckDecay(genpart,311, -321, 211)==1){
return 175;}//D0 decays to K0 K- pi+ 

if ( PcheckDecay(genpart,-311, 321, -211)==1){
return 176;}//D0 decays to anti-K0 K+ pi- 

if ( PcheckDecay(genpart,333, 113)==1){
return 177;}//D0 decays to phi rho0 

if ( PcheckDecay(genpart,333, 211, -211)==1){
return 178;}//D0 decays to phi pi+ pi- 

if ( PcheckDecay(genpart,321, -321, 113)==1){
return 179;}//D0 decays to K+ K- rho0 

if ( PcheckDecay(genpart,321, -321, 111, 111)==1){
return 180;}//D0 decays to K+ K- pi0 pi0 

if ( PcheckDecay(genpart,313, -321, 211)==1){
return 181;}//D0 decays to K*0 K- pi+ 

if ( PcheckDecay(genpart,-313, 321, -211)==1){
return 182;}//D0 decays to anti-K*0 K+ pi- 

if ( PcheckDecay(genpart,-311, 311, 211, -211)==1){
return 183;}//D0 decays to anti-K0 K0 pi+ pi- 

if ( PcheckDecay(genpart,321, -321, -321, 211)==1){
return 184;}//D0 decays to K+ K- K- pi+ 

if ( PcheckDecay(genpart,321, -321, -311, 111)==1){
return 185;}//D0 decays to K+ K- anti-K0 pi0 

if ( PcheckDecay(genpart,321, -321, 211, -211, 111)==1){
return 186;}//D0 decays to K+ K- pi+ pi- pi0 

if ( PcheckDecay(genpart,211, -211)==1){
return 187;}//D0 decays to pi+ pi- 

if ( PcheckDecay(genpart,111, 111)==1){
return 188;}//D0 decays to pi0 pi0 

if ( PcheckDecay(genpart,221, 111)==1){
return 189;}//D0 decays to eta pi0 

if ( PcheckDecay(genpart,331, 111)==1){
return 190;}//D0 decays to eta' pi0 

if ( PcheckDecay(genpart,221, 221)==1){
return 191;}//D0 decays to eta eta 

if ( PcheckDecay(genpart,213, -211)==1){
return 192;}//D0 decays to rho+ pi- 

if ( PcheckDecay(genpart,-213, 211)==1){
return 193;}//D0 decays to rho- pi+ 

if ( PcheckDecay(genpart,113, 111)==1){
return 194;}//D0 decays to rho0 pi0 

if ( PcheckDecay(genpart,111, 111, 111)==1){
return 195;}//D0 decays to pi0 pi0 pi0 

if ( PcheckDecay(genpart,211, 211, -211, -211)==1){
return 196;}//D0 decays to pi+ pi+ pi- pi- 

if ( PcheckDecay(genpart,211, -211, 111, 111)==1){
return 197;}//D0 decays to pi+ pi- pi0 pi0 

if ( PcheckDecay(genpart,211, -211, 211, -211, 111)==1){
return 198;}//D0 decays to pi+ pi- pi+ pi- pi0 

if ( PcheckDecay(genpart,211, -211, 111, 111, 111)==1){
return 199;}//D0 decays to pi+ pi- pi0 pi0 pi0 

if ( PcheckDecay(genpart,211, -211, 211, -211, 211, -211)==1){
return 200;}//D0 decays to pi+ pi- pi+ pi- pi+ pi- 

if ( PcheckDecay(genpart,-211, 321)==1){
return 201;}//D0 decays to pi- K+ 

if ( PcheckDecay(genpart,323, -211)==1){
return 202;}//D0 decays to K*+ pi- 

if ( PcheckDecay(genpart,-211, 321, 111)==1){
return 203;}//D0 decays to pi- K+ pi0 

if ( PcheckDecay(genpart,-211, 321, 211, -211)==1){
return 204;}//D0 decays to pi- K+ pi+ pi- 

if ( PcheckDecay(genpart,310, 111, 111)==1){
  return 3804;}//D0 decays to ks pi0 pi0

 if ( PcheckDecay(genpart,2212, 11)==1){
   return 3806;}//D0 decays to p+ e-   
 
return -200;
 }

if (AmId(genpart)==411){
  if ( PcheckDecay(genpart,-313, -11, 12)==1){
    return 205;}//D+ decays to anti-K*0 e+ nu_e 
  
  if ( PcheckDecay(genpart,-311, -11, 12)==1){
    return 206;}//D+ decays to anti-K0 e+ nu_e 
  
if ( PcheckDecay(genpart,-10313, -11, 12)==1){
return 207;}//D+ decays to anti-K_10 e+ nu_e 

if ( PcheckDecay(genpart,-315, -11, 12)==1){
return 208;}//D+ decays to anti-K_2*0 e+ nu_e 

if ( PcheckDecay(genpart,111, -11, 12)==1){
return 209;}//D+ decays to pi0 e+ nu_e 

if ( PcheckDecay(genpart,221, -11, 12)==1){
return 210;}//D+ decays to eta e+ nu_e 

if ( PcheckDecay(genpart,331, -11, 12)==1){
return 211;}//D+ decays to eta' e+ nu_e 

if ( PcheckDecay(genpart,113, -11, 12)==1){
return 212;}//D+ decays to rho0 e+ nu_e 

if ( PcheckDecay(genpart,223, -11, 12)==1){
return 213;}//D+ decays to omega e+ nu_e 

if ( PcheckDecay(genpart,10113, -11, 12)==1){
return 214;}//D+ decays to b_10 e+ nu_e 

if ( PcheckDecay(genpart,-321, 211, -11, 12)==1){
return 215;}//D+ decays to K- pi+ e+ nu_e 

if ( PcheckDecay(genpart,-311, 111, -11, 12)==1){
return 216;}//D+ decays to anti-K0 pi0 e+ nu_e 

if ( PcheckDecay(genpart,-313, -13, 14)==1){
return 217;}//D+ decays to anti-K*0 mu+ nu_mu 

if ( PcheckDecay(genpart,-311, -13, 14)==1){
return 218;}//D+ decays to anti-K0 mu+ nu_mu 

if ( PcheckDecay(genpart,-10313, -13, 14)==1){
return 219;}//D+ decays to anti-K_10 mu+ nu_mu 

if ( PcheckDecay(genpart,-315, -13, 14)==1){
return 220;}//D+ decays to anti-K_2*0 mu+ nu_mu 

if ( PcheckDecay(genpart,111, -13, 14)==1){
return 221;}//D+ decays to pi0 mu+ nu_mu 

if ( PcheckDecay(genpart,221, -13, 14)==1){
return 222;}//D+ decays to eta mu+ nu_mu 

if ( PcheckDecay(genpart,331, -13, 14)==1){
return 223;}//D+ decays to eta' mu+ nu_mu 

if ( PcheckDecay(genpart,113, -13, 14)==1){
return 224;}//D+ decays to rho0 mu+ nu_mu 

if ( PcheckDecay(genpart,223, -13, 14)==1){
return 225;}//D+ decays to omega mu+ nu_mu 

if ( PcheckDecay(genpart,10113, -13, 14)==1){
return 226;}//D+ decays to b_10 mu+ nu_mu 

if ( PcheckDecay(genpart,-321, 211, -13, 14)==1){
return 227;}//D+ decays to K- pi+ mu+ nu_mu 

if ( PcheckDecay(genpart,-311, 111, -13, 14)==1){
return 228;}//D+ decays to anti-K0 pi0 mu+ nu_mu 

if ( PcheckDecay(genpart,-13, 14)==1){
return 229;}//D+ decays to mu+ nu_mu 

if ( PcheckDecay(genpart,-15, 16)==1){
return 230;}//D+ decays to tau+ nu_tau 

if ( PcheckDecay(genpart,-311, 211)==1){
return 231;}//D+ decays to anti-K0 pi+ 

if ( PcheckDecay(genpart,213, -311)==1){
return 232;}//D+ decays to rho+ anti-K0 

if ( PcheckDecay(genpart,20213, -311)==1){
return 233;}//D+ decays to a_1+ anti-K0 

if ( PcheckDecay(genpart,-313, 211)==1){
return 234;}//D+ decays to anti-K*0 pi+ 

if ( PcheckDecay(genpart,-20313, 211)==1){
return 235;}//D+ decays to anti-K'_10 pi+ 

if ( PcheckDecay(genpart,-10311, 211)==1){
return 236;}//D+ decays to anti-K_0*0 pi+ 

if ( PcheckDecay(genpart,-30313, 211)==1){
return 237;}//D+ decays to anti-K''*0 pi+ 

if ( PcheckDecay(genpart,-313, 213)==1){
return 238;}//D+ decays to anti-K*0 rho+ 

if ( PcheckDecay(genpart,-321, 211, 211)==1){
return 239;}//D+ decays to K- pi+ pi+ 

if ( PcheckDecay(genpart,-311, 211, 111)==1){
return 240;}//D+ decays to anti-K0 pi+ pi0 

if ( PcheckDecay(genpart,-311, 221, 211)==1){
return 241;}//D+ decays to anti-K0 eta pi+ 

if ( PcheckDecay(genpart,-311, 113, 211)==1){
return 242;}//D+ decays to anti-K0 rho0 pi+ 

if ( PcheckDecay(genpart,-311, 223, 211)==1){
return 243;}//D+ decays to anti-K0 omega pi+ 

if ( PcheckDecay(genpart,-321, 213, 211)==1){
return 244;}//D+ decays to K- rho+ pi+ 

if ( PcheckDecay(genpart,-323, 211, 211)==1){
return 245;}//D+ decays to K*- pi+ pi+ 

if ( PcheckDecay(genpart,-313, 111, 211)==1){
return 246;}//D+ decays to anti-K*0 pi0 pi+ 

if ( PcheckDecay(genpart,-313, 221, 211)==1){
return 247;}//D+ decays to anti-K*0 eta pi+ 

if ( PcheckDecay(genpart,-313, 113, 211)==1){
return 248;}//D+ decays to anti-K*0 rho0 pi+ 

if ( PcheckDecay(genpart,-313, 223, 211)==1){
return 249;}//D+ decays to anti-K*0 omega pi+ 

if ( PcheckDecay(genpart,-323, 213, 211)==1){
return 250;}//D+ decays to K*- rho+ pi+ 

if ( PcheckDecay(genpart,-321, 211, 211, 111)==1){
return 251;}//D+ decays to K- pi+ pi+ pi0 

if ( PcheckDecay(genpart,-311, 211, 111, 111)==1){
return 252;}//D+ decays to anti-K0 pi+ pi0 pi0 

if ( PcheckDecay(genpart,-321, 113, 211, 211)==1){
return 253;}//D+ decays to K- rho0 pi+ pi+ 

if ( PcheckDecay(genpart,-313, 211, 211, -211)==1){
return 254;}//D+ decays to anti-K*0 pi+ pi+ pi- 

if ( PcheckDecay(genpart,-321, 211, 211, 211, -211)==1){
return 255;}//D+ decays to K- pi+ pi+ pi+ pi- 

if ( PcheckDecay(genpart,-321, 211, 211, 111, 111)==1){
return 256;}//D+ decays to K- pi+ pi+ pi0 pi0 

if ( PcheckDecay(genpart,-311, 211, 211, -211, 111)==1){
return 257;}//D+ decays to anti-K0 pi+ pi+ pi- pi0 

if ( PcheckDecay(genpart,-311, 211, 111, 111, 111)==1){
return 258;}//D+ decays to anti-K0 pi+ pi0 pi0 pi0 

if ( PcheckDecay(genpart,-311, 211, 211, 211, -211, -211)==1){
return 259;}//D+ decays to anti-K0 pi+ pi+ pi+ pi- pi- 

if ( PcheckDecay(genpart,-311, -311, 321)==1){
return 260;}//D+ decays to anti-K0 anti-K0 K+ 

if ( PcheckDecay(genpart,-311, 321, -321, 211)==1){
return 261;}//D+ decays to anti-K0 K+ K- pi+ 

if ( PcheckDecay(genpart,333, 211)==1){
return 262;}//D+ decays to phi pi+ 

if ( PcheckDecay(genpart,333, 211, 111)==1){
return 263;}//D+ decays to phi pi+ pi0 

if ( PcheckDecay(genpart,-311, 321)==1){
return 264;}//D+ decays to anti-K0 K+ 

if ( PcheckDecay(genpart,-313, 321)==1){
return 265;}//D+ decays to anti-K*0 K+ 

if ( PcheckDecay(genpart,323, -311)==1){
return 266;}//D+ decays to K*+ anti-K0 

if ( PcheckDecay(genpart,-313, 323)==1){
return 267;}//D+ decays to anti-K*0 K*+ 

if ( PcheckDecay(genpart,321, -321, 211)==1){
return 268;}//D+ decays to K+ K- pi+ 

if ( PcheckDecay(genpart,321, -311, 111)==1){
return 269;}//D+ decays to K+ anti-K0 pi0 

if ( PcheckDecay(genpart,-311, 311, 211)==1){
return 270;}//D+ decays to anti-K0 K0 pi+ 

if ( PcheckDecay(genpart,323, -321, 211)==1){
return 271;}//D+ decays to K*+ K- pi+ 

if ( PcheckDecay(genpart,321, -323, 211)==1){
return 272;}//D+ decays to K+ K*- pi+ 

if ( PcheckDecay(genpart,323, -311, 111)==1){
return 273;}//D+ decays to K*+ anti-K0 pi0 

if ( PcheckDecay(genpart,321, -313, 111)==1){
return 274;}//D+ decays to K+ anti-K*0 pi0 

if ( PcheckDecay(genpart,-313, 311, 211)==1){
return 275;}//D+ decays to anti-K*0 K0 pi+ 

if ( PcheckDecay(genpart,-311, 313, 211)==1){
return 276;}//D+ decays to anti-K0 K*0 pi+ 

if ( PcheckDecay(genpart,111, 211)==1){
return 277;}//D+ decays to pi0 pi+ 

if ( PcheckDecay(genpart,113, 211)==1){
return 278;}//D+ decays to rho0 pi+ 

if ( PcheckDecay(genpart,211, 211, -211)==1){
return 279;}//D+ decays to pi+ pi+ pi- 

if ( PcheckDecay(genpart,211, 111, 111)==1){
return 280;}//D+ decays to pi+ pi0 pi0 

if ( PcheckDecay(genpart,211, 211, -211, 111)==1){
return 281;}//D+ decays to pi+ pi+ pi- pi0 

if ( PcheckDecay(genpart,211, 111, 111, 111)==1){
return 282;}//D+ decays to pi+ pi0 pi0 pi0 

if ( PcheckDecay(genpart,211, 211, 211, -211, -211)==1){
return 283;}//D+ decays to pi+ pi+ pi+ pi- pi- 

if ( PcheckDecay(genpart,211, 211, 211, -211, -211, 111)==1){
return 284;}//D+ decays to pi+ pi+ pi+ pi- pi- pi0 

if ( PcheckDecay(genpart,221, 211)==1){
return 285;}//D+ decays to eta pi+ 

if ( PcheckDecay(genpart,331, 211)==1){
return 286;}//D+ decays to eta' pi+ 

if ( PcheckDecay(genpart,225, 211)==1){
return 287;}//D+ decays to f_2 pi+ 

if ( PcheckDecay(genpart,221, 211, 111)==1){
return 288;}//D+ decays to eta pi+ pi0 

if ( PcheckDecay(genpart,221, 211, 211, -211)==1){
return 289;}//D+ decays to eta pi+ pi+ pi- 

if ( PcheckDecay(genpart,221, 211, 111, 111)==1){
return 290;}//D+ decays to eta pi+ pi0 pi0 

if ( PcheckDecay(genpart,313, 211)==1){
return 291;}//D+ decays to K*0 pi+ 

if ( PcheckDecay(genpart,113, 321)==1){
return 292;}//D+ decays to rho0 K+ 

if ( PcheckDecay(genpart,321, 211, -211)==1){
return 293;}//D+ decays to K+ pi+ pi- 

if ( PcheckDecay(genpart,321, 321, -321)==1){
return 294;}//D+ decays to K+ K+ K- 
 return -300;

 }


if (AmId(genpart)==-411){

  if ( PcheckDecay(genpart,313, 11, -12)==1){
    return 295;}//D- decays to K*0 e- anti-nu_e 
  
  if ( PcheckDecay(genpart,311, 11, -12)==1){
    return 296;}//D- decays to K0 e- anti-nu_e 
  
  if ( PcheckDecay(genpart,10313, 11, -12)==1){
return 297;}//D- decays to K_10 e- anti-nu_e 

if ( PcheckDecay(genpart,315, 11, -12)==1){
return 298;}//D- decays to K_2*0 e- anti-nu_e 

if ( PcheckDecay(genpart,111, 11, -12)==1){
return 299;}//D- decays to pi0 e- anti-nu_e 

if ( PcheckDecay(genpart,221, 11, -12)==1){
return 300;}//D- decays to eta e- anti-nu_e 

if ( PcheckDecay(genpart,331, 11, -12)==1){
return 301;}//D- decays to eta' e- anti-nu_e 

if ( PcheckDecay(genpart,113, 11, -12)==1){
return 302;}//D- decays to rho0 e- anti-nu_e 

if ( PcheckDecay(genpart,223, 11, -12)==1){
return 303;}//D- decays to omega e- anti-nu_e 
 
 if ( PcheckDecay(genpart,10113, 11, -12)==1){
   return 304;}//D- decays to b_10 e- anti-nu_e 
 
 if ( PcheckDecay(genpart,321, -211, 11, -12)==1){
   return 305;}//D- decays to K+ pi- e- anti-nu_e 
 
 if ( PcheckDecay(genpart,311, 111, 11, -12)==1){
return 306;}//D- decays to K0 pi0 e- anti-nu_e 

if ( PcheckDecay(genpart,313, 13, -14)==1){
return 307;}//D- decays to K*0 mu- anti-nu_mu 

if ( PcheckDecay(genpart,311, 13, -14)==1){
return 308;}//D- decays to K0 mu- anti-nu_mu 

if ( PcheckDecay(genpart,10313, 13, -14)==1){
return 309;}//D- decays to K_10 mu- anti-nu_mu 

if ( PcheckDecay(genpart,315, 13, -14)==1){
return 310;}//D- decays to K_2*0 mu- anti-nu_mu 

if ( PcheckDecay(genpart,111, 13, -14)==1){
return 311;}//D- decays to pi0 mu- anti-nu_mu 

if ( PcheckDecay(genpart,221, 13, -14)==1){
return 312;}//D- decays to eta mu- anti-nu_mu 

if ( PcheckDecay(genpart,331, 13, -14)==1){
return 313;}//D- decays to eta' mu- anti-nu_mu 

if ( PcheckDecay(genpart,113, 13, -14)==1){
return 314;}//D- decays to rho0 mu- anti-nu_mu 

if ( PcheckDecay(genpart,223, 13, -14)==1){
return 315;}//D- decays to omega mu- anti-nu_mu 

if ( PcheckDecay(genpart,10113, 13, -14)==1){
return 316;}//D- decays to b_10 mu- anti-nu_mu 

if ( PcheckDecay(genpart,321, -211, 13, -14)==1){
return 317;}//D- decays to K+ pi- mu- anti-nu_mu 

if ( PcheckDecay(genpart,311, 111, 13, -14)==1){
return 318;}//D- decays to K0 pi0 mu- anti-nu_mu 

if ( PcheckDecay(genpart,13, -14)==1){
return 319;}//D- decays to mu- anti-nu_mu 

if ( PcheckDecay(genpart,15, -16)==1){
return 320;}//D- decays to tau- anti-nu_tau 

if ( PcheckDecay(genpart,311, -211)==1){
return 321;}//D- decays to K0 pi- 

if ( PcheckDecay(genpart,-213, 311)==1){
return 322;}//D- decays to rho- K0 

if ( PcheckDecay(genpart,-20213, 311)==1){
return 323;}//D- decays to a_1- K0 

if ( PcheckDecay(genpart,313, -211)==1){
return 324;}//D- decays to K*0 pi- 

if ( PcheckDecay(genpart,20313, -211)==1){
return 325;}//D- decays to K'_10 pi- 

if ( PcheckDecay(genpart,10311, -211)==1){
return 326;}//D- decays to K_0*0 pi- 

if ( PcheckDecay(genpart,30313, -211)==1){
return 327;}//D- decays to K''*0 pi- 

if ( PcheckDecay(genpart,313, -213)==1){
return 328;}//D- decays to K*0 rho- 

if ( PcheckDecay(genpart,321, -211, -211)==1){
return 329;}//D- decays to K+ pi- pi- 

if ( PcheckDecay(genpart,311, -211, 111)==1){
return 330;}//D- decays to K0 pi- pi0 

if ( PcheckDecay(genpart,311, 221, -211)==1){
return 331;}//D- decays to K0 eta pi- 

if ( PcheckDecay(genpart,311, 113, -211)==1){
return 332;}//D- decays to K0 rho0 pi- 

if ( PcheckDecay(genpart,311, 223, -211)==1){
return 333;}//D- decays to K0 omega pi- 

if ( PcheckDecay(genpart,321, -213, -211)==1){
return 334;}//D- decays to K+ rho- pi- 

if ( PcheckDecay(genpart,323, -211, -211)==1){
return 335;}//D- decays to K*+ pi- pi- 

if ( PcheckDecay(genpart,313, 111, -211)==1){
return 336;}//D- decays to K*0 pi0 pi- 

if ( PcheckDecay(genpart,313, 221, -211)==1){
return 337;}//D- decays to K*0 eta pi- 

if ( PcheckDecay(genpart,313, 113, -211)==1){
return 338;}//D- decays to K*0 rho0 pi- 

if ( PcheckDecay(genpart,313, 223, -211)==1){
  return 339;}//D- decays to K*0 omega pi- 
 
 if ( PcheckDecay(genpart,323, -213, -211)==1){
  return 340;}//D- decays to K*+ rho- pi- 
 
 if ( PcheckDecay(genpart,321, -211, -211, 111)==1){
  return 341;}//D- decays to K+ pi- pi- pi0 
 
 if ( PcheckDecay(genpart,311, -211, 111, 111)==1){
  return 342;}//D- decays to K0 pi- pi0 pi0 
 
 if ( PcheckDecay(genpart,321, 113, -211, -211)==1){
  return 343;}//D- decays to K+ rho0 pi- pi- 
 
 if ( PcheckDecay(genpart,313, -211, -211, 211)==1){
  return 344;}//D- decays to K*0 pi- pi- pi+ 
 
 if ( PcheckDecay(genpart,321, -211, -211, -211, 211)==1){
  return 345;}//D- decays to K+ pi- pi- pi- pi+ 
 
 if ( PcheckDecay(genpart,321, -211, -211, 111, 111)==1){
  return 346;}//D- decays to K+ pi- pi- pi0 pi0 
 
 if ( PcheckDecay(genpart,311, -211, -211, 211, 111)==1){
  return 347;}//D- decays to K0 pi- pi- pi+ pi0 
 
 if ( PcheckDecay(genpart,311, -211, 111, 111, 111)==1){
return 348;}//D- decays to K0 pi- pi0 pi0 pi0 
 
 if ( PcheckDecay(genpart,311, -211, -211, -211, 211, 211)==1){
  return 349;}//D- decays to K0 pi- pi- pi- pi+ pi+ 
 
 if ( PcheckDecay(genpart,311, 311, -321)==1){
  return 350;}//D- decays to K0 K0 K- 
 
 if ( PcheckDecay(genpart,311, 321, -321, -211)==1){
  return 351;}//D- decays to K0 K+ K- pi- 
 
 if ( PcheckDecay(genpart,333, -211)==1){
  return 352;}//D- decays to phi pi- 
 
 if ( PcheckDecay(genpart,333, -211, 111)==1){
  return 353;}//D- decays to phi pi- pi0 
 
 if ( PcheckDecay(genpart,311, -321)==1){
  return 354;}//D- decays to K0 K- 
 
 if ( PcheckDecay(genpart,313, -321)==1){
  return 355;}//D- decays to K*0 K- 
 
 if ( PcheckDecay(genpart,-323, 311)==1){
  return 356;}//D- decays to K*- K0 
 
 if ( PcheckDecay(genpart,313, -323)==1){
  return 357;}//D- decays to K*0 K*- 
 
 if ( PcheckDecay(genpart,-321, 321, -211)==1){
  return 358;}//D- decays to K- K+ pi- 
 
 if ( PcheckDecay(genpart,-321, 311, 111)==1){
  return 359;}//D- decays to K- K0 pi0 
 
 if ( PcheckDecay(genpart,-311, 311, -211)==1){
  return 360;}//D- decays to anti-K0 K0 pi- 
 
if ( PcheckDecay(genpart,-323, 321, -211)==1){
  return 361;}//D- decays to K*- K+ pi- 

 if ( PcheckDecay(genpart,-321, 323, -211)==1){
  return 362;}//D- decays to K- K*+ pi- 
 
 if ( PcheckDecay(genpart,-323, 311, 111)==1){
  return 363;}//D- decays to K*- K0 pi0 
 
 if ( PcheckDecay(genpart,-321, 313, 111)==1){
  return 364;}//D- decays to K- K*0 pi0 
 
 if ( PcheckDecay(genpart,313, -311, -211)==1){
  return 365;}//D- decays to K*0 anti-K0 pi- 
 
 if ( PcheckDecay(genpart,311, -313, -211)==1){
  return 366;}//D- decays to K0 anti-K*0 pi- 
 
 if ( PcheckDecay(genpart,111, -211)==1){
  return 367;}//D- decays to pi0 pi- 
 
 if ( PcheckDecay(genpart,113, -211)==1){
  return 368;}//D- decays to rho0 pi- 
 
 if ( PcheckDecay(genpart,-211, 211, -211)==1){
  return 369;}//D- decays to pi- pi+ pi- 
 
 if ( PcheckDecay(genpart,-211, 111, 111)==1){
  return 370;}//D- decays to pi- pi0 pi0 
 
 if ( PcheckDecay(genpart,-211, 211, -211, 111)==1){
  return 371;}//D- decays to pi- pi+ pi- pi0 
 
 if ( PcheckDecay(genpart,-211, 111, 111, 111)==1){
  return 372;}//D- decays to pi- pi0 pi0 pi0 
 
 if ( PcheckDecay(genpart,-211, -211, -211, 211, 211)==1){
  return 373;}//D- decays to pi- pi- pi- pi+ pi+ 
 
 if ( PcheckDecay(genpart,-211, -211, -211, 211, 211, 111)==1){
  return 374;}//D- decays to pi- pi- pi- pi+ pi+ pi0 
 
 if ( PcheckDecay(genpart,221, -211)==1){
  return 375;}//D- decays to eta pi- 
 
 if ( PcheckDecay(genpart,331, -211)==1){
  return 376;}//D- decays to eta' pi- 
 
 if ( PcheckDecay(genpart,225, -211)==1){
  return 377;}//D- decays to f_2 pi- 
 
 if ( PcheckDecay(genpart,221, -211, 111)==1){
  return 378;}//D- decays to eta pi- pi0 
 
 if ( PcheckDecay(genpart,221, -211, 211, -211)==1){
  return 379;}//D- decays to eta pi- pi+ pi- 
 
 if ( PcheckDecay(genpart,221, -211, 111, 111)==1){
  return 380;}//D- decays to eta pi- pi0 pi0 
 
 if ( PcheckDecay(genpart,-313, -211)==1){
  return 381;}//D- decays to anti-K*0 pi- 
 
 if ( PcheckDecay(genpart,113, -321)==1){
  return 382;}//D- decays to rho0 K- 

 if ( PcheckDecay(genpart,-321, -211, 211)==1){
  return 383;}//D- decays to K- pi- pi+ 
 
 if ( PcheckDecay(genpart,-321, -321, 321)==1){
  return 384;}//D- decays to K- K- K+ 
 
 return -400;
 }

 
  }


  int  genT::Braremode(Gen_hepevt genpart)
{

  if(AmId(genpart)==-511){
  if ( PcheckDecay(genpart,-313, 22)==1){
return 1;}//anti-B0 decays to anti-K*0 gamma 

if ( PcheckDecay(genpart,-10313, 22)==1){
return 2;}//anti-B0 decays to anti-K_10 gamma 

if ( PcheckDecay(genpart,-315, 22)==1){
return 3;}//anti-B0 decays to anti-K_2*0 gamma 

if ( PcheckDecay(genpart,-30343, 22)==1){
return 4;}//anti-B0 decays to anti-Xsd gamma 

if ( PcheckDecay(genpart,113, 22)==1){
return 5;}//anti-B0 decays to rho0 gamma 

if ( PcheckDecay(genpart,223, 22)==1){
return 6;}//anti-B0 decays to omega gamma 

if ( PcheckDecay(genpart,30643, 22)==1){
return 7;}//anti-B0 decays to Xdd gamma 

if ( PcheckDecay(genpart,-311, -11, 11)==1){
return 8;}//anti-B0 decays to anti-K0 e+ e- 

if ( PcheckDecay(genpart,-313, -11, 11)==1){
return 9;}//anti-B0 decays to anti-K*0 e+ e- 

if ( PcheckDecay(genpart,-30343, -11, 11)==1){
return 10;}//anti-B0 decays to anti-Xsd e+ e- 

if ( PcheckDecay(genpart,-311, -13, 13)==1){
return 11;}//anti-B0 decays to anti-K0 mu+ mu- 

if ( PcheckDecay(genpart,-313, -13, 13)==1){
return 12;}//anti-B0 decays to anti-K*0 mu+ mu- 

if ( PcheckDecay(genpart,-30343, -13, 13)==1){
return 13;}//anti-B0 decays to anti-Xsd mu+ mu- 

if ( PcheckDecay(genpart,-311, -15, 15)==1){
return 14;}//anti-B0 decays to anti-K0 tau+ tau- 

if ( PcheckDecay(genpart,-313, -15, 15)==1){
return 15;}//anti-B0 decays to anti-K*0 tau+ tau- 

if ( PcheckDecay(genpart,-30343, -15, 15)==1){
return 16;}//anti-B0 decays to anti-Xsd tau+ tau- 

if ( PcheckDecay(genpart,-311, 12, -12)==1){
return 17;}//anti-B0 decays to anti-K0 nu_e anti-nu_e 

if ( PcheckDecay(genpart,-313, 12, -12)==1){
return 18;}//anti-B0 decays to anti-K*0 nu_e anti-nu_e 

if ( PcheckDecay(genpart,-30343, 12, -12)==1){
return 19;}//anti-B0 decays to anti-Xsd nu_e anti-nu_e 

if ( PcheckDecay(genpart,-311, 14, -14)==1){
return 20;}//anti-B0 decays to anti-K0 nu_mu anti-nu_mu 

if ( PcheckDecay(genpart,-313, 14, -14)==1){
return 21;}//anti-B0 decays to anti-K*0 nu_mu anti-nu_mu 

if ( PcheckDecay(genpart,-30343, 14, -14)==1){
return 22;}//anti-B0 decays to anti-Xsd nu_mu anti-nu_mu 

if ( PcheckDecay(genpart,-311, 16, -16)==1){
return 23;}//anti-B0 decays to anti-K0 nu_tau anti-nu_tau 

if ( PcheckDecay(genpart,-313, 16, -16)==1){
return 24;}//anti-B0 decays to anti-K*0 nu_tau anti-nu_tau 

if ( PcheckDecay(genpart,-30343, 16, -16)==1){
return 25;}//anti-B0 decays to anti-Xsd nu_tau anti-nu_tau 

if ( PcheckDecay(genpart,-321, 321)==1){
return 26;}//anti-B0 decays to K- K+ 

if ( PcheckDecay(genpart,-321, 211)==1){
return 27;}//anti-B0 decays to K- pi+ 

if ( PcheckDecay(genpart,-311, 111)==1){
return 28;}//anti-B0 decays to anti-K0 pi0 

if ( PcheckDecay(genpart,130, 130)==1){
return 29;}//anti-B0 decays to K_L0 K_L0 

if ( PcheckDecay(genpart,310, 310)==1){
return 30;}//anti-B0 decays to K_S0 K_S0 

if ( PcheckDecay(genpart,-211, 211)==1){
return 31;}//anti-B0 decays to pi- pi+ 

if ( PcheckDecay(genpart,111, 111)==1){
return 32;}//anti-B0 decays to pi0 pi0 

if ( PcheckDecay(genpart,-100323, 211)==1){
return 33;}//anti-B0 decays to K'*- pi+ 

if ( PcheckDecay(genpart,-100313, 111)==1){
return 34;}//anti-B0 decays to anti-K'*0 pi0 

if ( PcheckDecay(genpart,-10321, 211)==1){
return 35;}//anti-B0 decays to K_0*- pi+ 

if ( PcheckDecay(genpart,-10211, 211)==1){
return 36;}//anti-B0 decays to a_0- pi+ 

if ( PcheckDecay(genpart,10331, -311)==1){
return 37;}//anti-B0 decays to f'_0 anti-K0 

if ( PcheckDecay(genpart,10221, -311)==1){
return 38;}//anti-B0 decays to f_0 anti-K0 

if ( PcheckDecay(genpart,9030221, -311)==1){
return 39;}//anti-B0 decays to f_0(1500) anti-K0 

if ( PcheckDecay(genpart,10221, 10221)==1){
return 40;}//anti-B0 decays to f_0 f_0 

if ( PcheckDecay(genpart,10111, -311)==1){
return 41;}//anti-B0 decays to a_00 anti-K0 

if ( PcheckDecay(genpart,10211, -321)==1){
return 42;}//anti-B0 decays to a_0+ K- 

if ( PcheckDecay(genpart,225, -311)==1){
return 43;}//anti-B0 decays to f_2 anti-K0 

if ( PcheckDecay(genpart,-325, 211)==1){
return 44;}//anti-B0 decays to K_2*- pi+ 

if ( PcheckDecay(genpart,-315, 111)==1){
return 45;}//anti-B0 decays to anti-K_2*0 pi0 

if ( PcheckDecay(genpart,-323, 321)==1){
return 46;}//anti-B0 decays to K*- K+ 

if ( PcheckDecay(genpart,-323, 211)==1){
return 47;}//anti-B0 decays to K*- pi+ 

if ( PcheckDecay(genpart,323, -321)==1){
return 48;}//anti-B0 decays to K*+ K- 

if ( PcheckDecay(genpart,-313, 311)==1){
return 49;}//anti-B0 decays to anti-K*0 K0 

if ( PcheckDecay(genpart,313, -311)==1){
return 50;}//anti-B0 decays to K*0 anti-K0 

if ( PcheckDecay(genpart,-313, 111)==1){
return 51;}//anti-B0 decays to anti-K*0 pi0 

if ( PcheckDecay(genpart,-313, 10221)==1){
return 52;}//anti-B0 decays to anti-K*0 f_0 

if ( PcheckDecay(genpart,100113, -311)==1){
return 53;}//anti-B0 decays to rho(2S)0 anti-K0 

if ( PcheckDecay(genpart,30113, -311)==1){
return 54;}//anti-B0 decays to rho(3S)0 anti-K0 

if ( PcheckDecay(genpart,-213, 211)==1){
return 55;}//anti-B0 decays to rho- pi+ 

if ( PcheckDecay(genpart,213, -321)==1){
return 56;}//anti-B0 decays to rho+ K- 

if ( PcheckDecay(genpart,213, -211)==1){
return 57;}//anti-B0 decays to rho+ pi- 

if ( PcheckDecay(genpart,113, -311)==1){
return 58;}//anti-B0 decays to rho0 anti-K0 

if ( PcheckDecay(genpart,113, 111)==1){
return 59;}//anti-B0 decays to rho0 pi0 

if ( PcheckDecay(genpart,113, 10221)==1){
return 60;}//anti-B0 decays to rho0 f_0 

if ( PcheckDecay(genpart,-100213, 211)==1){
return 61;}//anti-B0 decays to rho(2S)- pi+ 

if ( PcheckDecay(genpart,-30213, 211)==1){
return 62;}//anti-B0 decays to rho(3S)- pi+ 

if ( PcheckDecay(genpart,100213, -211)==1){
return 63;}//anti-B0 decays to rho(2S)+ pi- 

if ( PcheckDecay(genpart,30213, -211)==1){
return 64;}//anti-B0 decays to rho(3S)+ pi- 

if ( PcheckDecay(genpart,100213, -321)==1){
return 65;}//anti-B0 decays to rho(2S)+ K- 

if ( PcheckDecay(genpart,30213, -321)==1){
return 66;}//anti-B0 decays to rho(3S)+ K- 

if ( PcheckDecay(genpart,-10323, 211)==1){
return 67;}//anti-B0 decays to K_1- pi+ 

if ( PcheckDecay(genpart,-20323, 211)==1){
return 68;}//anti-B0 decays to K'_1- pi+ 

if ( PcheckDecay(genpart,-30323, 211)==1){
return 69;}//anti-B0 decays to K''*- pi+ 

if ( PcheckDecay(genpart,-30313, 111)==1){
return 70;}//anti-B0 decays to anti-K''*0 pi0 

if ( PcheckDecay(genpart,10213, -321)==1){
return 71;}//anti-B0 decays to b_1+ K- 

if ( PcheckDecay(genpart,10213, -211)==1){
return 72;}//anti-B0 decays to b_1+ pi- 

if ( PcheckDecay(genpart,-10213, 211)==1){
return 73;}//anti-B0 decays to b_1- pi+ 

if ( PcheckDecay(genpart,-321, 321, -311)==1){
return 74;}//anti-B0 decays to K- K+ anti-K0 

if ( PcheckDecay(genpart,-321, 321, 111)==1){
return 75;}//anti-B0 decays to K- K+ pi0 

if ( PcheckDecay(genpart,-321, 311, 211)==1){
return 76;}//anti-B0 decays to K- K0 pi+ 

if ( PcheckDecay(genpart,-321, 211, 111)==1){
return 77;}//anti-B0 decays to K- pi+ pi0 

if ( PcheckDecay(genpart,-311, -211, 211)==1){
return 78;}//anti-B0 decays to anti-K0 pi- pi+ 

if ( PcheckDecay(genpart,310, 310, 310)==1){
return 79;}//anti-B0 decays to K_S0 K_S0 K_S0 

if ( PcheckDecay(genpart,130, 130, 130)==1){
return 80;}//anti-B0 decays to K_L0 K_L0 K_L0 

if ( PcheckDecay(genpart,310, 310, 130)==1){
return 81;}//anti-B0 decays to K_S0 K_S0 K_L0 

if ( PcheckDecay(genpart,310, 130, 130)==1){
return 82;}//anti-B0 decays to K_S0 K_L0 K_L0 

if ( PcheckDecay(genpart,-211, 211, 111, 211)==1){
return 83;}//anti-B0 decays to pi- pi+ pi0 pi+ 

if ( PcheckDecay(genpart,10221, 211, -211)==1){
return 84;}//anti-B0 decays to f_0 pi+ pi- 

if ( PcheckDecay(genpart,335, -313)==1){
return 85;}//anti-B0 decays to f'_2 anti-K*0 

if ( PcheckDecay(genpart,225, -313)==1){
return 86;}//anti-B0 decays to f_2 anti-K*0 

if ( PcheckDecay(genpart,213, -10321)==1){
return 87;}//anti-B0 decays to rho+ K_0*- 

if ( PcheckDecay(genpart,113, -10311)==1){
return 88;}//anti-B0 decays to rho0 anti-K_0*0 

if ( PcheckDecay(genpart,20213, -321)==1){
return 89;}//anti-B0 decays to a_1+ K- 

if ( PcheckDecay(genpart,20213, -211)==1){
return 90;}//anti-B0 decays to a_1+ pi- 

if ( PcheckDecay(genpart,-20213, 211)==1){
return 91;}//anti-B0 decays to a_1- pi+ 

if ( PcheckDecay(genpart,-323, 323)==1){
return 92;}//anti-B0 decays to K*- K*+ 

if ( PcheckDecay(genpart,-323, 213)==1){
return 93;}//anti-B0 decays to K*- rho+ 

if ( PcheckDecay(genpart,-313, 313)==1){
return 94;}//anti-B0 decays to anti-K*0 K*0 

if ( PcheckDecay(genpart,-313, -313)==1){
return 95;}//anti-B0 decays to anti-K*0 anti-K*0 

if ( PcheckDecay(genpart,-313, 113)==1){
return 96;}//anti-B0 decays to anti-K*0 rho0 

if ( PcheckDecay(genpart,-213, 213)==1){
return 97;}//anti-B0 decays to rho- rho+ 

if ( PcheckDecay(genpart,113, 113)==1){
return 98;}//anti-B0 decays to rho0 rho0 

if ( PcheckDecay(genpart,-213, 20213)==1){
return 99;}//anti-B0 decays to rho- a_1+ 

if ( PcheckDecay(genpart,213, -20213)==1){
return 100;}//anti-B0 decays to rho+ a_1- 

if ( PcheckDecay(genpart,-323, 211, 111)==1){
return 101;}//anti-B0 decays to K*- pi+ pi0 

if ( PcheckDecay(genpart,-313, -211, 211)==1){
return 102;}//anti-B0 decays to anti-K*0 pi- pi+ 

if ( PcheckDecay(genpart,-313, 111, 111)==1){
return 103;}//anti-B0 decays to anti-K*0 pi0 pi0 

if ( PcheckDecay(genpart,213, -321, 111)==1){
return 104;}//anti-B0 decays to rho+ K- pi0 

if ( PcheckDecay(genpart,213, -311, -211)==1){
return 105;}//anti-B0 decays to rho+ anti-K0 pi- 

if ( PcheckDecay(genpart,113, -321, 211)==1){
return 106;}//anti-B0 decays to rho0 K- pi+ 

if ( PcheckDecay(genpart,113, -311, 111)==1){
return 107;}//anti-B0 decays to rho0 anti-K0 pi0 

if ( PcheckDecay(genpart,113, 211, -211)==1){
return 108;}//anti-B0 decays to rho0 pi+ pi- 

if ( PcheckDecay(genpart,321, -321, -313)==1){
return 109;}//anti-B0 decays to K+ K- anti-K*0 

if ( PcheckDecay(genpart,-311, 311, -313)==1){
return 110;}//anti-B0 decays to anti-K0 K0 anti-K*0 

if ( PcheckDecay(genpart,-311, 313, -311)==1){
return 111;}//anti-B0 decays to anti-K0 K*0 anti-K0 

if ( PcheckDecay(genpart,321, -323, -311)==1){
return 112;}//anti-B0 decays to K+ K*- anti-K0 

if ( PcheckDecay(genpart,323, -321, -311)==1){
return 113;}//anti-B0 decays to K*+ K- anti-K0 

if ( PcheckDecay(genpart,-311, 221)==1){
return 114;}//anti-B0 decays to anti-K0 eta 

if ( PcheckDecay(genpart,-311, 331)==1){
return 115;}//anti-B0 decays to anti-K0 eta' 

if ( PcheckDecay(genpart,-311, 100221)==1){
return 116;}//anti-B0 decays to anti-K0 eta(2S) 

if ( PcheckDecay(genpart,-311, 9020221)==1){
return 117;}//anti-B0 decays to anti-K0 eta(1405) 

if ( PcheckDecay(genpart,-10311, 221)==1){
return 118;}//anti-B0 decays to anti-K_0*0 eta 

if ( PcheckDecay(genpart,111, 221)==1){
return 119;}//anti-B0 decays to pi0 eta 

if ( PcheckDecay(genpart,111, 331)==1){
return 120;}//anti-B0 decays to pi0 eta' 

if ( PcheckDecay(genpart,111, 100221)==1){
return 121;}//anti-B0 decays to pi0 eta(2S) 

if ( PcheckDecay(genpart,111, 9020221)==1){
return 122;}//anti-B0 decays to pi0 eta(1405) 

if ( PcheckDecay(genpart,221, 221)==1){
return 123;}//anti-B0 decays to eta eta 

if ( PcheckDecay(genpart,221, 331)==1){
return 124;}//anti-B0 decays to eta eta' 

if ( PcheckDecay(genpart,221, 100221)==1){
return 125;}//anti-B0 decays to eta eta(2S) 

if ( PcheckDecay(genpart,221, 9020221)==1){
return 126;}//anti-B0 decays to eta eta(1405) 

if ( PcheckDecay(genpart,331, 331)==1){
return 127;}//anti-B0 decays to eta' eta' 

if ( PcheckDecay(genpart,331, 100221)==1){
return 128;}//anti-B0 decays to eta' eta(2S) 

if ( PcheckDecay(genpart,331, 9020221)==1){
return 129;}//anti-B0 decays to eta' eta(1405) 

if ( PcheckDecay(genpart,100221, 100221)==1){
return 130;}//anti-B0 decays to eta(2S) eta(2S) 

if ( PcheckDecay(genpart,100221, 9020221)==1){
return 131;}//anti-B0 decays to eta(2S) eta(1405) 

if ( PcheckDecay(genpart,9020221, 9020221)==1){
return 132;}//anti-B0 decays to eta(1405) eta(1405) 

if ( PcheckDecay(genpart,10221, 221)==1){
return 133;}//anti-B0 decays to f_0 eta 

if ( PcheckDecay(genpart,10221, 331)==1){
return 134;}//anti-B0 decays to f_0 eta' 

if ( PcheckDecay(genpart,223, -311)==1){
return 135;}//anti-B0 decays to omega anti-K0 

if ( PcheckDecay(genpart,333, -311)==1){
return 136;}//anti-B0 decays to phi anti-K0 

if ( PcheckDecay(genpart,333, -10311)==1){
return 137;}//anti-B0 decays to phi anti-K_0*0 

if ( PcheckDecay(genpart,333, -30313)==1){
return 138;}//anti-B0 decays to phi anti-K''*0 

if ( PcheckDecay(genpart,333, -317)==1){
return 139;}//anti-B0 decays to phi anti-K_3*0 

if ( PcheckDecay(genpart,333, -319)==1){
return 140;}//anti-B0 decays to phi anti-K_4*0 

if ( PcheckDecay(genpart,223, 111)==1){
return 141;}//anti-B0 decays to omega pi0 

if ( PcheckDecay(genpart,333, 111)==1){
return 142;}//anti-B0 decays to phi pi0 

if ( PcheckDecay(genpart,-313, 221)==1){
return 143;}//anti-B0 decays to anti-K*0 eta 

if ( PcheckDecay(genpart,113, 221)==1){
return 144;}//anti-B0 decays to rho0 eta 

if ( PcheckDecay(genpart,223, 221)==1){
return 145;}//anti-B0 decays to omega eta 

if ( PcheckDecay(genpart,333, 221)==1){
return 146;}//anti-B0 decays to phi eta 

if ( PcheckDecay(genpart,-313, 331)==1){
return 147;}//anti-B0 decays to anti-K*0 eta' 

if ( PcheckDecay(genpart,113, 331)==1){
return 148;}//anti-B0 decays to rho0 eta' 

if ( PcheckDecay(genpart,223, 331)==1){
return 149;}//anti-B0 decays to omega eta' 

if ( PcheckDecay(genpart,333, 331)==1){
return 150;}//anti-B0 decays to phi eta' 

if ( PcheckDecay(genpart,-313, 100221)==1){
return 151;}//anti-B0 decays to anti-K*0 eta(2S) 

if ( PcheckDecay(genpart,113, 100221)==1){
return 152;}//anti-B0 decays to rho0 eta(2S) 

if ( PcheckDecay(genpart,223, 100221)==1){
return 153;}//anti-B0 decays to omega eta(2S) 

if ( PcheckDecay(genpart,333, 100221)==1){
return 154;}//anti-B0 decays to phi eta(2S) 

if ( PcheckDecay(genpart,-313, 9020221)==1){
return 155;}//anti-B0 decays to anti-K*0 eta(1405) 

if ( PcheckDecay(genpart,113, 9020221)==1){
return 156;}//anti-B0 decays to rho0 eta(1405) 

if ( PcheckDecay(genpart,223, 9020221)==1){
return 157;}//anti-B0 decays to omega eta(1405) 

if ( PcheckDecay(genpart,333, 9020221)==1){
return 158;}//anti-B0 decays to phi eta(1405) 

if ( PcheckDecay(genpart,223, 10221)==1){
return 159;}//anti-B0 decays to omega f_0 

if ( PcheckDecay(genpart,-313, 223)==1){
return 160;}//anti-B0 decays to anti-K*0 omega 

if ( PcheckDecay(genpart,-313, 333)==1){
return 161;}//anti-B0 decays to anti-K*0 phi 

if ( PcheckDecay(genpart,113, 223)==1){
return 162;}//anti-B0 decays to rho0 omega 

if ( PcheckDecay(genpart,113, 333)==1){
return 163;}//anti-B0 decays to rho0 phi 

if ( PcheckDecay(genpart,223, 223)==1){
return 164;}//anti-B0 decays to omega omega 

if ( PcheckDecay(genpart,223, 333)==1){
return 165;}//anti-B0 decays to omega phi 

if ( PcheckDecay(genpart,333, 333)==1){
return 166;}//anti-B0 decays to phi phi 

if ( PcheckDecay(genpart,-315, 333)==1){
return 167;}//anti-B0 decays to anti-K_2*0 phi 

if ( PcheckDecay(genpart,-315, 221)==1){
return 168;}//anti-B0 decays to anti-K_2*0 eta 

if ( PcheckDecay(genpart,-30343, 331)==1){
return 169;}//anti-B0 decays to anti-Xsd eta' 

if ( PcheckDecay(genpart,-30343, 221)==1){
return 170;}//anti-B0 decays to anti-Xsd eta 

if ( PcheckDecay(genpart,221, -321, 211)==1){
return 171;}//anti-B0 decays to eta K- pi+ 

if ( PcheckDecay(genpart,221, -211, 211)==1){
return 172;}//anti-B0 decays to eta pi- pi+ 

if ( PcheckDecay(genpart,221, -311, 111)==1){
return 173;}//anti-B0 decays to eta anti-K0 pi0 

if ( PcheckDecay(genpart,221, 111, 111)==1){
return 174;}//anti-B0 decays to eta pi0 pi0 

if ( PcheckDecay(genpart,331, -321, 211)==1){
return 175;}//anti-B0 decays to eta' K- pi+ 

if ( PcheckDecay(genpart,331, -211, 211)==1){
return 176;}//anti-B0 decays to eta' pi- pi+ 

if ( PcheckDecay(genpart,331, -311, 111)==1){
return 177;}//anti-B0 decays to eta' anti-K0 pi0 

if ( PcheckDecay(genpart,331, 111, 111)==1){
return 178;}//anti-B0 decays to eta' pi0 pi0 

if ( PcheckDecay(genpart,221, 221, -311)==1){
return 179;}//anti-B0 decays to eta eta anti-K0 

if ( PcheckDecay(genpart,221, 221, 111)==1){
return 180;}//anti-B0 decays to eta eta pi0 

if ( PcheckDecay(genpart,221, 331, -311)==1){
return 181;}//anti-B0 decays to eta eta' anti-K0 

if ( PcheckDecay(genpart,221, 331, 111)==1){
return 182;}//anti-B0 decays to eta eta' pi0 

if ( PcheckDecay(genpart,331, 331, -311)==1){
return 183;}//anti-B0 decays to eta' eta' anti-K0 

if ( PcheckDecay(genpart,331, 331, 111)==1){
return 184;}//anti-B0 decays to eta' eta' pi0 

if ( PcheckDecay(genpart,221, -323, 211)==1){
return 185;}//anti-B0 decays to eta K*- pi+ 

if ( PcheckDecay(genpart,221, -321, 213)==1){
return 186;}//anti-B0 decays to eta K- rho+ 

if ( PcheckDecay(genpart,221, -213, 211)==1){
return 187;}//anti-B0 decays to eta rho- pi+ 

if ( PcheckDecay(genpart,221, -211, 213)==1){
return 188;}//anti-B0 decays to eta pi- rho+ 

if ( PcheckDecay(genpart,221, -313, 111)==1){
return 189;}//anti-B0 decays to eta anti-K*0 pi0 

if ( PcheckDecay(genpart,221, -311, 113)==1){
return 190;}//anti-B0 decays to eta anti-K0 rho0 

if ( PcheckDecay(genpart,221, 113, 111)==1){
return 191;}//anti-B0 decays to eta rho0 pi0 

if ( PcheckDecay(genpart,331, -323, 211)==1){
return 192;}//anti-B0 decays to eta' K*- pi+ 

if ( PcheckDecay(genpart,331, -321, 213)==1){
return 193;}//anti-B0 decays to eta' K- rho+ 

if ( PcheckDecay(genpart,331, -213, 211)==1){
return 194;}//anti-B0 decays to eta' rho- pi+ 

if ( PcheckDecay(genpart,331, -211, 213)==1){
return 195;}//anti-B0 decays to eta' pi- rho+ 

if ( PcheckDecay(genpart,331, -313, 111)==1){
return 196;}//anti-B0 decays to eta' anti-K*0 pi0 

if ( PcheckDecay(genpart,331, -311, 113)==1){
return 197;}//anti-B0 decays to eta' anti-K0 rho0 

if ( PcheckDecay(genpart,331, 113, 111)==1){
return 198;}//anti-B0 decays to eta' rho0 pi0 

if ( PcheckDecay(genpart,221, 221, -313)==1){
return 199;}//anti-B0 decays to eta eta anti-K*0 

if ( PcheckDecay(genpart,221, 221, 113)==1){
return 200;}//anti-B0 decays to eta eta rho0 

if ( PcheckDecay(genpart,221, 331, -313)==1){
return 201;}//anti-B0 decays to eta eta' anti-K*0 

if ( PcheckDecay(genpart,221, 331, 113)==1){
return 202;}//anti-B0 decays to eta eta' rho0 

if ( PcheckDecay(genpart,331, 331, -313)==1){
return 203;}//anti-B0 decays to eta' eta' anti-K*0 

if ( PcheckDecay(genpart,331, 331, 113)==1){
return 204;}//anti-B0 decays to eta' eta' rho0 

if ( PcheckDecay(genpart,223, -321, 211)==1){
return 205;}//anti-B0 decays to omega K- pi+ 

if ( PcheckDecay(genpart,223, -211, 211)==1){
return 206;}//anti-B0 decays to omega pi- pi+ 

if ( PcheckDecay(genpart,223, -311, 111)==1){
return 207;}//anti-B0 decays to omega anti-K0 pi0 

if ( PcheckDecay(genpart,223, 111, 111)==1){
return 208;}//anti-B0 decays to omega pi0 pi0 

if ( PcheckDecay(genpart,333, -321, 211)==1){
return 209;}//anti-B0 decays to phi K- pi+ 

if ( PcheckDecay(genpart,333, -211, 211)==1){
return 210;}//anti-B0 decays to phi pi- pi+ 

if ( PcheckDecay(genpart,333, -311, 111)==1){
return 211;}//anti-B0 decays to phi anti-K0 pi0 

if ( PcheckDecay(genpart,333, 111, 111)==1){
return 212;}//anti-B0 decays to phi pi0 pi0 

if ( PcheckDecay(genpart,223, 221, -311)==1){
return 213;}//anti-B0 decays to omega eta anti-K0 

if ( PcheckDecay(genpart,223, 221, 111)==1){
return 214;}//anti-B0 decays to omega eta pi0 

if ( PcheckDecay(genpart,223, 331, -311)==1){
return 215;}//anti-B0 decays to omega eta' anti-K0 

if ( PcheckDecay(genpart,223, 331, 111)==1){
return 216;}//anti-B0 decays to omega eta' pi0 

if ( PcheckDecay(genpart,333, 221, -311)==1){
return 217;}//anti-B0 decays to phi eta anti-K0 

if ( PcheckDecay(genpart,333, 221, 111)==1){
return 218;}//anti-B0 decays to phi eta pi0 

if ( PcheckDecay(genpart,333, 331, -311)==1){
return 219;}//anti-B0 decays to phi eta' anti-K0 

if ( PcheckDecay(genpart,333, 331, 111)==1){
return 220;}//anti-B0 decays to phi eta' pi0 

if ( PcheckDecay(genpart,333, 333, -311)==1){
return 221;}//anti-B0 decays to phi phi anti-K0 

if ( PcheckDecay(genpart,10111, 111)==1){
return 222;}//anti-B0 decays to a_00 pi0 

if ( PcheckDecay(genpart,10111, 221)==1){
return 223;}//anti-B0 decays to a_00 eta 

if ( PcheckDecay(genpart,10111, 331)==1){
return 224;}//anti-B0 decays to a_00 eta' 

if ( PcheckDecay(genpart,10111, 100221)==1){
return 225;}//anti-B0 decays to a_00 eta(2S) 

if ( PcheckDecay(genpart,10111, 9020221)==1){
return 226;}//anti-B0 decays to a_00 eta(1405) 

if ( PcheckDecay(genpart,10111, 10111)==1){
return 227;}//anti-B0 decays to a_00 a_00 

if ( PcheckDecay(genpart,10211, -10211)==1){
return 228;}//anti-B0 decays to a_0+ a_0- 

if ( PcheckDecay(genpart,10111, 10221)==1){
return 229;}//anti-B0 decays to a_00 f_0 

if ( PcheckDecay(genpart,-313, 10111)==1){
return 230;}//anti-B0 decays to anti-K*0 a_00 

if ( PcheckDecay(genpart,-323, 10211)==1){
return 231;}//anti-B0 decays to K*- a_0+ 

if ( PcheckDecay(genpart,113, 10111)==1){
return 232;}//anti-B0 decays to rho0 a_00 

if ( PcheckDecay(genpart,-213, 10211)==1){
return 233;}//anti-B0 decays to rho- a_0+ 

if ( PcheckDecay(genpart,213, -10211)==1){
return 234;}//anti-B0 decays to rho+ a_0- 

if ( PcheckDecay(genpart,223, 10111)==1){
return 235;}//anti-B0 decays to omega a_00 

if ( PcheckDecay(genpart,333, 10111)==1){
return 236;}//anti-B0 decays to phi a_00 

if ( PcheckDecay(genpart,20113, 10111)==1){
return 237;}//anti-B0 decays to a_10 a_00 

if ( PcheckDecay(genpart,20213, -10211)==1){
return 238;}//anti-B0 decays to a_1+ a_0- 

if ( PcheckDecay(genpart,-20213, 10211)==1){
return 239;}//anti-B0 decays to a_1- a_0+ 

if ( PcheckDecay(genpart,20223, 10111)==1){
return 240;}//anti-B0 decays to f_1 a_00 

if ( PcheckDecay(genpart,10113, 10111)==1){
return 241;}//anti-B0 decays to b_10 a_00 

if ( PcheckDecay(genpart,10213, -10211)==1){
return 242;}//anti-B0 decays to b_1+ a_0- 

if ( PcheckDecay(genpart,-10213, 10211)==1){
return 243;}//anti-B0 decays to b_1- a_0+ 

if ( PcheckDecay(genpart,10223, 10111)==1){
return 244;}//anti-B0 decays to h_1 a_00 

if ( PcheckDecay(genpart,20113, -311)==1){
return 245;}//anti-B0 decays to a_10 anti-K0 

if ( PcheckDecay(genpart,20113, 111)==1){
return 246;}//anti-B0 decays to a_10 pi0 

if ( PcheckDecay(genpart,20113, 221)==1){
return 247;}//anti-B0 decays to a_10 eta 

if ( PcheckDecay(genpart,20113, 331)==1){
return 248;}//anti-B0 decays to a_10 eta' 

if ( PcheckDecay(genpart,20113, 100221)==1){
return 249;}//anti-B0 decays to a_10 eta(2S) 

if ( PcheckDecay(genpart,20113, 9020221)==1){
return 250;}//anti-B0 decays to a_10 eta(1405) 

if ( PcheckDecay(genpart,20113, 10221)==1){
return 251;}//anti-B0 decays to a_10 f_0 

if ( PcheckDecay(genpart,-313, 20113)==1){
return 252;}//anti-B0 decays to anti-K*0 a_10 

if ( PcheckDecay(genpart,-323, 20213)==1){
return 253;}//anti-B0 decays to K*- a_1+ 

if ( PcheckDecay(genpart,113, 20113)==1){
return 254;}//anti-B0 decays to rho0 a_10 

if ( PcheckDecay(genpart,223, 20113)==1){
return 255;}//anti-B0 decays to omega a_10 

if ( PcheckDecay(genpart,333, 20113)==1){
return 256;}//anti-B0 decays to phi a_10 

if ( PcheckDecay(genpart,20113, 20113)==1){
return 257;}//anti-B0 decays to a_10 a_10 

if ( PcheckDecay(genpart,20213, -20213)==1){
return 258;}//anti-B0 decays to a_1+ a_1- 

if ( PcheckDecay(genpart,20223, 20113)==1){
return 259;}//anti-B0 decays to f_1 a_10 

if ( PcheckDecay(genpart,10113, 20113)==1){
return 260;}//anti-B0 decays to b_10 a_10 

if ( PcheckDecay(genpart,10213, -20213)==1){
return 261;}//anti-B0 decays to b_1+ a_1- 

if ( PcheckDecay(genpart,-10213, 20213)==1){
return 262;}//anti-B0 decays to b_1- a_1+ 

if ( PcheckDecay(genpart,10223, 20113)==1){
return 263;}//anti-B0 decays to h_1 a_10 

if ( PcheckDecay(genpart,10221, 111)==1){
return 264;}//anti-B0 decays to f_0 pi0 

if ( PcheckDecay(genpart,10221, 100221)==1){
return 265;}//anti-B0 decays to f_0 eta(2S) 

if ( PcheckDecay(genpart,10221, 9020221)==1){
return 266;}//anti-B0 decays to f_0 eta(1405) 

if ( PcheckDecay(genpart,333, 10221)==1){
return 267;}//anti-B0 decays to phi f_0 

if ( PcheckDecay(genpart,20223, 10221)==1){
return 268;}//anti-B0 decays to f_1 f_0 

if ( PcheckDecay(genpart,10113, 10221)==1){
return 269;}//anti-B0 decays to b_10 f_0 

if ( PcheckDecay(genpart,10223, 10221)==1){
return 270;}//anti-B0 decays to h_1 f_0 

if ( PcheckDecay(genpart,20223, -311)==1){
return 271;}//anti-B0 decays to f_1 anti-K0 

if ( PcheckDecay(genpart,20223, 111)==1){
return 272;}//anti-B0 decays to f_1 pi0 

if ( PcheckDecay(genpart,20223, 221)==1){
return 273;}//anti-B0 decays to f_1 eta 

if ( PcheckDecay(genpart,20223, 331)==1){
return 274;}//anti-B0 decays to f_1 eta' 

if ( PcheckDecay(genpart,20223, 100221)==1){
return 275;}//anti-B0 decays to f_1 eta(2S) 

if ( PcheckDecay(genpart,20223, 9020221)==1){
return 276;}//anti-B0 decays to f_1 eta(1405) 

if ( PcheckDecay(genpart,-313, 20223)==1){
return 277;}//anti-B0 decays to anti-K*0 f_1 

if ( PcheckDecay(genpart,113, 20223)==1){
return 278;}//anti-B0 decays to rho0 f_1 

if ( PcheckDecay(genpart,223, 20223)==1){
return 279;}//anti-B0 decays to omega f_1 

if ( PcheckDecay(genpart,333, 20223)==1){
return 280;}//anti-B0 decays to phi f_1 

if ( PcheckDecay(genpart,20223, 20223)==1){
return 281;}//anti-B0 decays to f_1 f_1 

if ( PcheckDecay(genpart,10113, 20223)==1){
return 282;}//anti-B0 decays to b_10 f_1 

if ( PcheckDecay(genpart,10223, 20223)==1){
return 283;}//anti-B0 decays to h_1 f_1 

if ( PcheckDecay(genpart,10113, -311)==1){
return 284;}//anti-B0 decays to b_10 anti-K0 

if ( PcheckDecay(genpart,10113, 111)==1){
return 285;}//anti-B0 decays to b_10 pi0 

if ( PcheckDecay(genpart,10113, 221)==1){
return 286;}//anti-B0 decays to b_10 eta 

if ( PcheckDecay(genpart,10113, 331)==1){
return 287;}//anti-B0 decays to b_10 eta' 

if ( PcheckDecay(genpart,10113, 100221)==1){
return 288;}//anti-B0 decays to b_10 eta(2S) 

if ( PcheckDecay(genpart,10113, 9020221)==1){
return 289;}//anti-B0 decays to b_10 eta(1405) 

if ( PcheckDecay(genpart,-313, 10113)==1){
return 290;}//anti-B0 decays to anti-K*0 b_10 

if ( PcheckDecay(genpart,-323, 10213)==1){
return 291;}//anti-B0 decays to K*- b_1+ 

if ( PcheckDecay(genpart,113, 10113)==1){
return 292;}//anti-B0 decays to rho0 b_10 

if ( PcheckDecay(genpart,-213, 10213)==1){
return 293;}//anti-B0 decays to rho- b_1+ 

if ( PcheckDecay(genpart,213, -10213)==1){
return 294;}//anti-B0 decays to rho+ b_1- 

if ( PcheckDecay(genpart,223, 10113)==1){
return 295;}//anti-B0 decays to omega b_10 

if ( PcheckDecay(genpart,333, 10113)==1){
return 296;}//anti-B0 decays to phi b_10 

if ( PcheckDecay(genpart,10113, 10113)==1){
return 297;}//anti-B0 decays to b_10 b_10 

if ( PcheckDecay(genpart,10213, -10213)==1){
return 298;}//anti-B0 decays to b_1+ b_1- 

if ( PcheckDecay(genpart,10223, 10113)==1){
return 299;}//anti-B0 decays to h_1 b_10 

if ( PcheckDecay(genpart,10223, -311)==1){
return 300;}//anti-B0 decays to h_1 anti-K0 

if ( PcheckDecay(genpart,10223, 111)==1){
return 301;}//anti-B0 decays to h_1 pi0 

if ( PcheckDecay(genpart,10223, 221)==1){
return 302;}//anti-B0 decays to h_1 eta 

if ( PcheckDecay(genpart,10223, 331)==1){
return 303;}//anti-B0 decays to h_1 eta' 

if ( PcheckDecay(genpart,10223, 100221)==1){
return 304;}//anti-B0 decays to h_1 eta(2S) 

if ( PcheckDecay(genpart,10223, 9020221)==1){
return 305;}//anti-B0 decays to h_1 eta(1405) 

if ( PcheckDecay(genpart,-313, 10223)==1){
return 306;}//anti-B0 decays to anti-K*0 h_1 

if ( PcheckDecay(genpart,113, 10223)==1){
return 307;}//anti-B0 decays to rho0 h_1 

if ( PcheckDecay(genpart,223, 10223)==1){
return 308;}//anti-B0 decays to omega h_1 

if ( PcheckDecay(genpart,333, 10223)==1){
return 309;}//anti-B0 decays to phi h_1 

if ( PcheckDecay(genpart,10223, 10223)==1){
return 310;}//anti-B0 decays to h_1 h_1 

if ( PcheckDecay(genpart,215, -211)==1){
return 311;}//anti-B0 decays to a_2+ pi- 

if ( PcheckDecay(genpart,-215, 211)==1){
return 312;}//anti-B0 decays to a_2- pi+ 

if ( PcheckDecay(genpart,115, 111)==1){
return 313;}//anti-B0 decays to a_20 pi0 

if ( PcheckDecay(genpart,225, 111)==1){
return 314;}//anti-B0 decays to f_2 pi0 

if ( PcheckDecay(genpart,335, 111)==1){
return 315;}//anti-B0 decays to f'_2 pi0 

if ( PcheckDecay(genpart,115, 221)==1){
return 316;}//anti-B0 decays to a_20 eta 

if ( PcheckDecay(genpart,225, 221)==1){
return 317;}//anti-B0 decays to f_2 eta 

if ( PcheckDecay(genpart,335, 221)==1){
return 318;}//anti-B0 decays to f'_2 eta 

if ( PcheckDecay(genpart,115, 331)==1){
return 319;}//anti-B0 decays to a_20 eta' 

if ( PcheckDecay(genpart,225, 331)==1){
return 320;}//anti-B0 decays to f_2 eta' 

if ( PcheckDecay(genpart,335, 331)==1){
return 321;}//anti-B0 decays to f'_2 eta' 

if ( PcheckDecay(genpart,115, 100221)==1){
return 322;}//anti-B0 decays to a_20 eta(2S) 

if ( PcheckDecay(genpart,225, 100221)==1){
return 323;}//anti-B0 decays to f_2 eta(2S) 

if ( PcheckDecay(genpart,335, 100221)==1){
return 324;}//anti-B0 decays to f'_2 eta(2S) 

if ( PcheckDecay(genpart,115, 9020221)==1){
return 325;}//anti-B0 decays to a_20 eta(1405) 

if ( PcheckDecay(genpart,225, 9020221)==1){
return 326;}//anti-B0 decays to f_2 eta(1405) 

if ( PcheckDecay(genpart,335, 9020221)==1){
return 327;}//anti-B0 decays to f'_2 eta(1405) 

if ( PcheckDecay(genpart,315, -311)==1){
return 328;}//anti-B0 decays to K_2*0 anti-K0 

if ( PcheckDecay(genpart,-315, 311)==1){
return 329;}//anti-B0 decays to anti-K_2*0 K0 

if ( PcheckDecay(genpart,215, -321)==1){
return 330;}//anti-B0 decays to a_2+ K- 

if ( PcheckDecay(genpart,115, -311)==1){
return 331;}//anti-B0 decays to a_20 anti-K0 

if ( PcheckDecay(genpart,335, -311)==1){
return 332;}//anti-B0 decays to f'_2 anti-K0 

if ( PcheckDecay(genpart,-315, 331)==1){
return 333;}//anti-B0 decays to anti-K_2*0 eta' 

if ( PcheckDecay(genpart,-315, 100221)==1){
return 334;}//anti-B0 decays to anti-K_2*0 eta(2S) 

if ( PcheckDecay(genpart,-315, 9020221)==1){
return 335;}//anti-B0 decays to anti-K_2*0 eta(1405) 

if ( PcheckDecay(genpart,215, -213)==1){
return 336;}//anti-B0 decays to a_2+ rho- 

if ( PcheckDecay(genpart,-215, 213)==1){
return 337;}//anti-B0 decays to a_2- rho+ 

if ( PcheckDecay(genpart,115, 113)==1){
return 338;}//anti-B0 decays to a_20 rho0 

if ( PcheckDecay(genpart,225, 113)==1){
return 339;}//anti-B0 decays to f_2 rho0 

if ( PcheckDecay(genpart,335, 113)==1){
return 340;}//anti-B0 decays to f'_2 rho0 

if ( PcheckDecay(genpart,115, 223)==1){
return 341;}//anti-B0 decays to a_20 omega 

if ( PcheckDecay(genpart,225, 223)==1){
return 342;}//anti-B0 decays to f_2 omega 

if ( PcheckDecay(genpart,335, 223)==1){
return 343;}//anti-B0 decays to f'_2 omega 

if ( PcheckDecay(genpart,115, 333)==1){
return 344;}//anti-B0 decays to a_20 phi 

if ( PcheckDecay(genpart,225, 333)==1){
return 345;}//anti-B0 decays to f_2 phi 

if ( PcheckDecay(genpart,335, 333)==1){
return 346;}//anti-B0 decays to f'_2 phi 

if ( PcheckDecay(genpart,315, -313)==1){
return 347;}//anti-B0 decays to K_2*0 anti-K*0 

if ( PcheckDecay(genpart,-315, 313)==1){
return 348;}//anti-B0 decays to anti-K_2*0 K*0 

if ( PcheckDecay(genpart,215, -323)==1){
return 349;}//anti-B0 decays to a_2+ K*- 

if ( PcheckDecay(genpart,115, -313)==1){
return 350;}//anti-B0 decays to a_20 anti-K*0 

if ( PcheckDecay(genpart,-315, 113)==1){
return 351;}//anti-B0 decays to anti-K_2*0 rho0 

if ( PcheckDecay(genpart,-325, 213)==1){
return 352;}//anti-B0 decays to K_2*- rho+ 

if ( PcheckDecay(genpart,-315, 223)==1){
return 353;}//anti-B0 decays to anti-K_2*0 omega 

if ( PcheckDecay(genpart,3222, -2212, 111)==1){
return 354;}//anti-B0 decays to Sigma+ anti-p- pi0 

if ( PcheckDecay(genpart,3122, -2112, 111)==1){
return 355;}//anti-B0 decays to Lambda0 anti-n0 pi0 

if ( PcheckDecay(genpart,3212, -2112, 111)==1){
return 356;}//anti-B0 decays to Sigma0 anti-n0 pi0 

if ( PcheckDecay(genpart,3112, -1114, 111)==1){
return 357;}//anti-B0 decays to Sigma- anti-Delta+ pi0 

if ( PcheckDecay(genpart,3122, -2212, 211)==1){
return 358;}//anti-B0 decays to Lambda0 anti-p- pi+ 

if ( PcheckDecay(genpart,3212, -2212, 211)==1){
return 359;}//anti-B0 decays to Sigma0 anti-p- pi+ 

if ( PcheckDecay(genpart,3112, -2112, 211)==1){
return 360;}//anti-B0 decays to Sigma- anti-n0 pi+ 

if ( PcheckDecay(genpart,3222, -2112, -211)==1){
return 361;}//anti-B0 decays to Sigma+ anti-n0 pi- 

if ( PcheckDecay(genpart,3122, -1114, -211)==1){
return 362;}//anti-B0 decays to Lambda0 anti-Delta+ pi- 

if ( PcheckDecay(genpart,3212, -1114, -211)==1){
return 363;}//anti-B0 decays to Sigma0 anti-Delta+ pi- 

if ( PcheckDecay(genpart,2212, -2212, -311)==1){
return 364;}//anti-B0 decays to p+ anti-p- anti-K0 

if ( PcheckDecay(genpart,2112, -2112, -311)==1){
return 365;}//anti-B0 decays to n0 anti-n0 anti-K0 

if ( PcheckDecay(genpart,1114, -1114, -311)==1){
return 366;}//anti-B0 decays to Delta- anti-Delta+ anti-K0 

if ( PcheckDecay(genpart,2224, -2212, -321)==1){
return 367;}//anti-B0 decays to Delta++ anti-p- K- 

if ( PcheckDecay(genpart,2212, -2112, -321)==1){
return 368;}//anti-B0 decays to p+ anti-n0 K- 

if ( PcheckDecay(genpart,2112, -1114, -321)==1){
return 369;}//anti-B0 decays to n0 anti-Delta+ K- 

if ( PcheckDecay(genpart,3322, -3122, 111)==1){
return 370;}//anti-B0 decays to Xi0 anti-Lambda0 pi0 

if ( PcheckDecay(genpart,3322, -3212, 111)==1){
return 371;}//anti-B0 decays to Xi0 anti-Sigma0 pi0 

if ( PcheckDecay(genpart,3312, -3112, 111)==1){
return 372;}//anti-B0 decays to Xi- anti-Sigma+ pi0 

if ( PcheckDecay(genpart,3312, -3122, 211)==1){
return 373;}//anti-B0 decays to Xi- anti-Lambda0 pi+ 

if ( PcheckDecay(genpart,3312, -3212, 211)==1){
return 374;}//anti-B0 decays to Xi- anti-Sigma0 pi+ 

if ( PcheckDecay(genpart,3322, -3112, -211)==1){
return 375;}//anti-B0 decays to Xi0 anti-Sigma+ pi- 

if ( PcheckDecay(genpart,3322, -3222, 211)==1){
return 376;}//anti-B0 decays to Xi0 anti-Sigma- pi+ 

if ( PcheckDecay(genpart,3322, -2112, 311)==1){
return 377;}//anti-B0 decays to Xi0 anti-n0 K0 

if ( PcheckDecay(genpart,3322, -2212, 321)==1){
return 378;}//anti-B0 decays to Xi0 anti-p- K+ 

if ( PcheckDecay(genpart,3312, -1114, 311)==1){
return 379;}//anti-B0 decays to Xi- anti-Delta+ K0 

if ( PcheckDecay(genpart,3312, -2112, 321)==1){
return 380;}//anti-B0 decays to Xi- anti-n0 K+ 

if ( PcheckDecay(genpart,3122, -3122, -311)==1){
return 381;}//anti-B0 decays to Lambda0 anti-Lambda0 anti-K0 

if ( PcheckDecay(genpart,3122, -3212, -311)==1){
return 382;}//anti-B0 decays to Lambda0 anti-Sigma0 anti-K0 

if ( PcheckDecay(genpart,3212, -3122, -311)==1){
return 383;}//anti-B0 decays to Sigma0 anti-Lambda0 anti-K0 

if ( PcheckDecay(genpart,3212, -3212, -311)==1){
return 384;}//anti-B0 decays to Sigma0 anti-Sigma0 anti-K0 

if ( PcheckDecay(genpart,3222, -3122, -321)==1){
return 385;}//anti-B0 decays to Sigma+ anti-Lambda0 K- 

if ( PcheckDecay(genpart,3222, -3212, -321)==1){
return 386;}//anti-B0 decays to Sigma+ anti-Sigma0 K- 

if ( PcheckDecay(genpart,3112, -3112, -311)==1){
return 387;}//anti-B0 decays to Sigma- anti-Sigma+ anti-K0 

if ( PcheckDecay(genpart,3122, -3112, -321)==1){
return 388;}//anti-B0 decays to Lambda0 anti-Sigma+ K- 

if ( PcheckDecay(genpart,3212, -3112, -321)==1){
return 389;}//anti-B0 decays to Sigma0 anti-Sigma+ K- 

if ( PcheckDecay(genpart,3334, -3312, 111)==1){
return 390;}//anti-B0 decays to Omega- anti-Xi+ pi0 

if ( PcheckDecay(genpart,3334, -3322, 211)==1){
return 391;}//anti-B0 decays to Omega- anti-Xi0 pi+ 

if ( PcheckDecay(genpart,3334, -3112, 311)==1){
return 392;}//anti-B0 decays to Omega- anti-Sigma+ K0 

if ( PcheckDecay(genpart,3334, -3122, 321)==1){
return 393;}//anti-B0 decays to Omega- anti-Lambda0 K+ 

if ( PcheckDecay(genpart,3334, -3212, 321)==1){
return 394;}//anti-B0 decays to Omega- anti-Sigma0 K+ 

if ( PcheckDecay(genpart,3312, -3312, -311)==1){
return 395;}//anti-B0 decays to Xi- anti-Xi+ anti-K0 

if ( PcheckDecay(genpart,3322, -3312, -321)==1){
return 396;}//anti-B0 decays to Xi0 anti-Xi+ K- 

if ( PcheckDecay(genpart,3222, -2212, 113)==1){
return 397;}//anti-B0 decays to Sigma+ anti-p- rho0 

if ( PcheckDecay(genpart,3122, -2112, 113)==1){
return 398;}//anti-B0 decays to Lambda0 anti-n0 rho0 

if ( PcheckDecay(genpart,3212, -2112, 113)==1){
return 399;}//anti-B0 decays to Sigma0 anti-n0 rho0 

if ( PcheckDecay(genpart,3112, -1114, 113)==1){
return 400;}//anti-B0 decays to Sigma- anti-Delta+ rho0 

if ( PcheckDecay(genpart,3122, -2212, 213)==1){
return 401;}//anti-B0 decays to Lambda0 anti-p- rho+ 

if ( PcheckDecay(genpart,3212, -2212, 213)==1){
return 402;}//anti-B0 decays to Sigma0 anti-p- rho+ 

if ( PcheckDecay(genpart,3112, -2112, 213)==1){
return 403;}//anti-B0 decays to Sigma- anti-n0 rho+ 

if ( PcheckDecay(genpart,3222, -2112, -213)==1){
return 404;}//anti-B0 decays to Sigma+ anti-n0 rho- 

if ( PcheckDecay(genpart,3122, -1114, -213)==1){
return 405;}//anti-B0 decays to Lambda0 anti-Delta+ rho- 

if ( PcheckDecay(genpart,3212, -1114, -213)==1){
return 406;}//anti-B0 decays to Sigma0 anti-Delta+ rho- 

if ( PcheckDecay(genpart,2212, -2212, -313)==1){
return 407;}//anti-B0 decays to p+ anti-p- anti-K*0 

if ( PcheckDecay(genpart,2112, -2112, -313)==1){
return 408;}//anti-B0 decays to n0 anti-n0 anti-K*0 

if ( PcheckDecay(genpart,1114, -1114, -313)==1){
return 409;}//anti-B0 decays to Delta- anti-Delta+ anti-K*0 

if ( PcheckDecay(genpart,2224, -2212, -323)==1){
return 410;}//anti-B0 decays to Delta++ anti-p- K*- 

if ( PcheckDecay(genpart,2212, -2112, -323)==1){
return 411;}//anti-B0 decays to p+ anti-n0 K*- 

if ( PcheckDecay(genpart,2112, -1114, -323)==1){
return 412;}//anti-B0 decays to n0 anti-Delta+ K*- 

if ( PcheckDecay(genpart,3322, -3122, 113)==1){
return 413;}//anti-B0 decays to Xi0 anti-Lambda0 rho0 

if ( PcheckDecay(genpart,3322, -3212, 113)==1){
return 414;}//anti-B0 decays to Xi0 anti-Sigma0 rho0 

if ( PcheckDecay(genpart,3312, -3112, 113)==1){
return 415;}//anti-B0 decays to Xi- anti-Sigma+ rho0 

if ( PcheckDecay(genpart,3312, -3122, 213)==1){
return 416;}//anti-B0 decays to Xi- anti-Lambda0 rho+ 

if ( PcheckDecay(genpart,3312, -3212, 213)==1){
return 417;}//anti-B0 decays to Xi- anti-Sigma0 rho+ 

if ( PcheckDecay(genpart,3322, -3112, -213)==1){
return 418;}//anti-B0 decays to Xi0 anti-Sigma+ rho- 

if ( PcheckDecay(genpart,3322, -3222, 213)==1){
return 419;}//anti-B0 decays to Xi0 anti-Sigma- rho+ 

if ( PcheckDecay(genpart,3322, -2112, 313)==1){
return 420;}//anti-B0 decays to Xi0 anti-n0 K*0 

if ( PcheckDecay(genpart,3322, -2212, 323)==1){
return 421;}//anti-B0 decays to Xi0 anti-p- K*+ 

if ( PcheckDecay(genpart,3312, -1114, 313)==1){
return 422;}//anti-B0 decays to Xi- anti-Delta+ K*0 

if ( PcheckDecay(genpart,3312, -2112, 323)==1){
return 423;}//anti-B0 decays to Xi- anti-n0 K*+ 

if ( PcheckDecay(genpart,3122, -3122, -313)==1){
return 424;}//anti-B0 decays to Lambda0 anti-Lambda0 anti-K*0 

if ( PcheckDecay(genpart,3122, -3212, -313)==1){
return 425;}//anti-B0 decays to Lambda0 anti-Sigma0 anti-K*0 

if ( PcheckDecay(genpart,3212, -3122, -313)==1){
return 426;}//anti-B0 decays to Sigma0 anti-Lambda0 anti-K*0 

if ( PcheckDecay(genpart,3212, -3212, -313)==1){
return 427;}//anti-B0 decays to Sigma0 anti-Sigma0 anti-K*0 

if ( PcheckDecay(genpart,3222, -3122, -323)==1){
return 428;}//anti-B0 decays to Sigma+ anti-Lambda0 K*- 

if ( PcheckDecay(genpart,3222, -3212, -323)==1){
return 429;}//anti-B0 decays to Sigma+ anti-Sigma0 K*- 

if ( PcheckDecay(genpart,3112, -3112, -313)==1){
return 430;}//anti-B0 decays to Sigma- anti-Sigma+ anti-K*0 

if ( PcheckDecay(genpart,3122, -3112, -323)==1){
return 431;}//anti-B0 decays to Lambda0 anti-Sigma+ K*- 

if ( PcheckDecay(genpart,3212, -3112, -323)==1){
return 432;}//anti-B0 decays to Sigma0 anti-Sigma+ K*- 

if ( PcheckDecay(genpart,3334, -3312, 113)==1){
return 433;}//anti-B0 decays to Omega- anti-Xi+ rho0 

if ( PcheckDecay(genpart,3334, -3112, 313)==1){
return 434;}//anti-B0 decays to Omega- anti-Sigma+ K*0 

if ( PcheckDecay(genpart,3334, -3122, 323)==1){
return 435;}//anti-B0 decays to Omega- anti-Lambda0 K*+ 

if ( PcheckDecay(genpart,3334, -3212, 323)==1){
return 436;}//anti-B0 decays to Omega- anti-Sigma0 K*+ 

if ( PcheckDecay(genpart,3312, -3312, -313)==1){
return 437;}//anti-B0 decays to Xi- anti-Xi+ anti-K*0 

if ( PcheckDecay(genpart,3322, -3312, -323)==1){
return 438;}//anti-B0 decays to Xi0 anti-Xi+ K*- 

if ( PcheckDecay(genpart,2212, -2212, 111)==1){
return 439;}//anti-B0 decays to p+ anti-p- pi0 

if ( PcheckDecay(genpart,2112, -2112, 111)==1){
return 440;}//anti-B0 decays to n0 anti-n0 pi0 

if ( PcheckDecay(genpart,1114, -1114, 111)==1){
return 441;}//anti-B0 decays to Delta- anti-Delta+ pi0 

if ( PcheckDecay(genpart,2212, -2112, -211)==1){
return 442;}//anti-B0 decays to p+ anti-n0 pi- 

if ( PcheckDecay(genpart,2212, -2224, 211)==1){
return 443;}//anti-B0 decays to p+ anti-Delta-- pi+ 

if ( PcheckDecay(genpart,2224, -2212, -211)==1){
return 444;}//anti-B0 decays to Delta++ anti-p- pi- 

if ( PcheckDecay(genpart,2112, -2212, 211)==1){
return 445;}//anti-B0 decays to n0 anti-p- pi+ 

if ( PcheckDecay(genpart,2112, -1114, -211)==1){
return 446;}//anti-B0 decays to n0 anti-Delta+ pi- 

if ( PcheckDecay(genpart,1114, -2112, 211)==1){
return 447;}//anti-B0 decays to Delta- anti-n0 pi+ 

if ( PcheckDecay(genpart,3122, -3122, 111)==1){
return 448;}//anti-B0 decays to Lambda0 anti-Lambda0 pi0 

if ( PcheckDecay(genpart,3122, -3212, 111)==1){
return 449;}//anti-B0 decays to Lambda0 anti-Sigma0 pi0 

if ( PcheckDecay(genpart,3212, -3122, 111)==1){
return 450;}//anti-B0 decays to Sigma0 anti-Lambda0 pi0 

if ( PcheckDecay(genpart,3212, -3212, 111)==1){
return 451;}//anti-B0 decays to Sigma0 anti-Sigma0 pi0 

if ( PcheckDecay(genpart,3122, -3112, -211)==1){
return 452;}//anti-B0 decays to Lambda0 anti-Sigma+ pi- 

if ( PcheckDecay(genpart,3122, -3222, 211)==1){
return 453;}//anti-B0 decays to Lambda0 anti-Sigma- pi+ 

if ( PcheckDecay(genpart,3222, -3122, -211)==1){
return 454;}//anti-B0 decays to Sigma+ anti-Lambda0 pi- 

if ( PcheckDecay(genpart,3112, -3122, 211)==1){
return 455;}//anti-B0 decays to Sigma- anti-Lambda0 pi+ 

if ( PcheckDecay(genpart,3212, -3112, -211)==1){
return 456;}//anti-B0 decays to Sigma0 anti-Sigma+ pi- 

if ( PcheckDecay(genpart,3212, -3222, 211)==1){
return 457;}//anti-B0 decays to Sigma0 anti-Sigma- pi+ 

if ( PcheckDecay(genpart,3222, -3212, -211)==1){
return 458;}//anti-B0 decays to Sigma+ anti-Sigma0 pi- 

if ( PcheckDecay(genpart,3112, -3212, 211)==1){
return 459;}//anti-B0 decays to Sigma- anti-Sigma0 pi+ 

if ( PcheckDecay(genpart,3122, -2212, 321)==1){
return 460;}//anti-B0 decays to Lambda0 anti-p- K+ 

if ( PcheckDecay(genpart,3212, -2212, 321)==1){
return 461;}//anti-B0 decays to Sigma0 anti-p- K+ 

if ( PcheckDecay(genpart,3122, -2112, 311)==1){
return 462;}//anti-B0 decays to Lambda0 anti-n0 K0 

if ( PcheckDecay(genpart,3212, -2112, 311)==1){
return 463;}//anti-B0 decays to Sigma0 anti-n0 K0 

if ( PcheckDecay(genpart,2212, -3122, -321)==1){
return 464;}//anti-B0 decays to p+ anti-Lambda0 K- 

if ( PcheckDecay(genpart,2212, -3212, -321)==1){
return 465;}//anti-B0 decays to p+ anti-Sigma0 K- 

if ( PcheckDecay(genpart,2112, -3122, -311)==1){
return 466;}//anti-B0 decays to n0 anti-Lambda0 anti-K0 

if ( PcheckDecay(genpart,2112, -3212, -311)==1){
return 467;}//anti-B0 decays to n0 anti-Sigma0 anti-K0 

if ( PcheckDecay(genpart,3312, -3312, 111)==1){
return 468;}//anti-B0 decays to Xi- anti-Xi+ pi0 

if ( PcheckDecay(genpart,3322, -3312, -211)==1){
return 469;}//anti-B0 decays to Xi0 anti-Xi+ pi- 

if ( PcheckDecay(genpart,3312, -3322, 211)==1){
return 470;}//anti-B0 decays to Xi- anti-Xi0 pi+ 

if ( PcheckDecay(genpart,3112, -3312, -311)==1){
return 471;}//anti-B0 decays to Sigma- anti-Xi+ anti-K0 

if ( PcheckDecay(genpart,3122, -3312, -321)==1){
return 472;}//anti-B0 decays to Lambda0 anti-Xi+ K- 

if ( PcheckDecay(genpart,3212, -3312, -321)==1){
return 473;}//anti-B0 decays to Sigma0 anti-Xi+ K- 

if ( PcheckDecay(genpart,3312, -3112, 311)==1){
return 474;}//anti-B0 decays to Xi- anti-Sigma+ K0 

if ( PcheckDecay(genpart,3312, -3122, 321)==1){
return 475;}//anti-B0 decays to Xi- anti-Lambda0 K+ 

if ( PcheckDecay(genpart,3312, -3212, 321)==1){
return 476;}//anti-B0 decays to Xi- anti-Sigma0 K+ 

if ( PcheckDecay(genpart,2212, -2212, 113)==1){
return 477;}//anti-B0 decays to p+ anti-p- rho0 

if ( PcheckDecay(genpart,2112, -2112, 113)==1){
return 478;}//anti-B0 decays to n0 anti-n0 rho0 

if ( PcheckDecay(genpart,1114, -1114, 113)==1){
return 479;}//anti-B0 decays to Delta- anti-Delta+ rho0 

if ( PcheckDecay(genpart,2212, -2112, -213)==1){
return 480;}//anti-B0 decays to p+ anti-n0 rho- 

if ( PcheckDecay(genpart,2112, -2212, 213)==1){
return 481;}//anti-B0 decays to n0 anti-p- rho+ 

if ( PcheckDecay(genpart,1114, -2112, 213)==1){
return 482;}//anti-B0 decays to Delta- anti-n0 rho+ 

if ( PcheckDecay(genpart,2112, -1114, -213)==1){
return 483;}//anti-B0 decays to n0 anti-Delta+ rho- 

if ( PcheckDecay(genpart,3122, -3122, 113)==1){
return 484;}//anti-B0 decays to Lambda0 anti-Lambda0 rho0 

if ( PcheckDecay(genpart,3122, -3212, 113)==1){
return 485;}//anti-B0 decays to Lambda0 anti-Sigma0 rho0 

if ( PcheckDecay(genpart,3212, -3122, 113)==1){
return 486;}//anti-B0 decays to Sigma0 anti-Lambda0 rho0 

if ( PcheckDecay(genpart,3212, -3212, 113)==1){
return 487;}//anti-B0 decays to Sigma0 anti-Sigma0 rho0 

if ( PcheckDecay(genpart,3122, -3112, -213)==1){
return 488;}//anti-B0 decays to Lambda0 anti-Sigma+ rho- 

if ( PcheckDecay(genpart,3122, -3222, 213)==1){
return 489;}//anti-B0 decays to Lambda0 anti-Sigma- rho+ 

if ( PcheckDecay(genpart,3222, -3122, -213)==1){
return 490;}//anti-B0 decays to Sigma+ anti-Lambda0 rho- 

if ( PcheckDecay(genpart,3112, -3122, 213)==1){
return 491;}//anti-B0 decays to Sigma- anti-Lambda0 rho+ 

if ( PcheckDecay(genpart,3212, -3112, -213)==1){
return 492;}//anti-B0 decays to Sigma0 anti-Sigma+ rho- 

if ( PcheckDecay(genpart,3212, -3222, 213)==1){
return 493;}//anti-B0 decays to Sigma0 anti-Sigma- rho+ 

if ( PcheckDecay(genpart,3222, -3212, -213)==1){
return 494;}//anti-B0 decays to Sigma+ anti-Sigma0 rho- 

if ( PcheckDecay(genpart,3112, -3212, 213)==1){
return 495;}//anti-B0 decays to Sigma- anti-Sigma0 rho+ 

if ( PcheckDecay(genpart,3122, -2212, 323)==1){
return 496;}//anti-B0 decays to Lambda0 anti-p- K*+ 

if ( PcheckDecay(genpart,3212, -2212, 323)==1){
return 497;}//anti-B0 decays to Sigma0 anti-p- K*+ 

if ( PcheckDecay(genpart,3122, -2112, 313)==1){
return 498;}//anti-B0 decays to Lambda0 anti-n0 K*0 

if ( PcheckDecay(genpart,3212, -2112, 313)==1){
return 499;}//anti-B0 decays to Sigma0 anti-n0 K*0 

if ( PcheckDecay(genpart,2212, -3122, -323)==1){
return 500;}//anti-B0 decays to p+ anti-Lambda0 K*- 

if ( PcheckDecay(genpart,2212, -3212, -323)==1){
return 501;}//anti-B0 decays to p+ anti-Sigma0 K*- 

if ( PcheckDecay(genpart,2112, -3122, -313)==1){
return 502;}//anti-B0 decays to n0 anti-Lambda0 anti-K*0 

if ( PcheckDecay(genpart,2112, -3212, -313)==1){
return 503;}//anti-B0 decays to n0 anti-Sigma0 anti-K*0 

if ( PcheckDecay(genpart,3312, -3312, 113)==1){
return 504;}//anti-B0 decays to Xi- anti-Xi+ rho0 

if ( PcheckDecay(genpart,3322, -3312, -213)==1){
return 505;}//anti-B0 decays to Xi0 anti-Xi+ rho- 

if ( PcheckDecay(genpart,3312, -3322, 213)==1){
return 506;}//anti-B0 decays to Xi- anti-Xi0 rho+ 

if ( PcheckDecay(genpart,3112, -3312, -313)==1){
return 507;}//anti-B0 decays to Sigma- anti-Xi+ anti-K*0 

if ( PcheckDecay(genpart,3122, -3312, -323)==1){
return 508;}//anti-B0 decays to Lambda0 anti-Xi+ K*- 

if ( PcheckDecay(genpart,3212, -3312, -323)==1){
return 509;}//anti-B0 decays to Sigma0 anti-Xi+ K*- 

if ( PcheckDecay(genpart,3312, -3112, 313)==1){
return 510;}//anti-B0 decays to Xi- anti-Sigma+ K*0 

if ( PcheckDecay(genpart,3312, -3122, 323)==1){
return 511;}//anti-B0 decays to Xi- anti-Lambda0 K*+ 

if ( PcheckDecay(genpart,3312, -3212, 323)==1){
return 512;}//anti-B0 decays to Xi- anti-Sigma0 K*+ 

if ( PcheckDecay(genpart,3222, -2212)==1){
return 513;}//anti-B0 decays to Sigma+ anti-p- 

if ( PcheckDecay(genpart,3122, -2112)==1){
return 514;}//anti-B0 decays to Lambda0 anti-n0 

if ( PcheckDecay(genpart,3212, -2112)==1){
return 515;}//anti-B0 decays to Sigma0 anti-n0 

if ( PcheckDecay(genpart,3112, -1114)==1){
return 516;}//anti-B0 decays to Sigma- anti-Delta+ 

if ( PcheckDecay(genpart,3322, -3122)==1){
return 517;}//anti-B0 decays to Xi0 anti-Lambda0 

if ( PcheckDecay(genpart,3322, -3212)==1){
return 518;}//anti-B0 decays to Xi0 anti-Sigma0 

if ( PcheckDecay(genpart,3312, -3112)==1){
return 519;}//anti-B0 decays to Xi- anti-Sigma+ 

if ( PcheckDecay(genpart,3334, -3312)==1){
return 520;}//anti-B0 decays to Omega- anti-Xi+ 

if ( PcheckDecay(genpart,2212, -2212)==1){
return 521;}//anti-B0 decays to p+ anti-p- 

if ( PcheckDecay(genpart,2112, -2112)==1){
return 522;}//anti-B0 decays to n0 anti-n0 

if ( PcheckDecay(genpart,1114, -1114)==1){
return 523;}//anti-B0 decays to Delta- anti-Delta+ 

if ( PcheckDecay(genpart,3122, -3122)==1){
return 524;}//anti-B0 decays to Lambda0 anti-Lambda0 

if ( PcheckDecay(genpart,3122, -3212)==1){
return 525;}//anti-B0 decays to Lambda0 anti-Sigma0 

if ( PcheckDecay(genpart,3212, -3122)==1){
return 526;}//anti-B0 decays to Sigma0 anti-Lambda0 

if ( PcheckDecay(genpart,3212, -3212)==1){
return 527;}//anti-B0 decays to Sigma0 anti-Sigma0 

if ( PcheckDecay(genpart,3312, -3312)==1){
return 528;}//anti-B0 decays to Xi- anti-Xi+ 

if ( PcheckDecay(genpart,-431, 211)==1){
return 529;}//anti-B0 decays to D_s- pi+ 

if ( PcheckDecay(genpart,213, -431)==1){
return 530;}//anti-B0 decays to rho+ D_s- 

if ( PcheckDecay(genpart,-431, 211, 111)==1){
return 531;}//anti-B0 decays to D_s- pi+ pi0 

if ( PcheckDecay(genpart,20213, -431)==1){
return 532;}//anti-B0 decays to a_1+ D_s- 

if ( PcheckDecay(genpart,-433, 211)==1){
return 533;}//anti-B0 decays to D_s*- pi+ 

if ( PcheckDecay(genpart,-433, 213)==1){
return 534;}//anti-B0 decays to D_s*- rho+ 

if ( PcheckDecay(genpart,-433, 211, 111)==1){
return 535;}//anti-B0 decays to D_s*- pi+ pi0 

if ( PcheckDecay(genpart,-433, 20213)==1){
return 536;}//anti-B0 decays to D_s*- a_1+ 

if ( PcheckDecay(genpart,431, -321)==1){
return 537;}//anti-B0 decays to D_s+ K- 

if ( PcheckDecay(genpart,-323, 431)==1){
return 538;}//anti-B0 decays to K*- D_s+ 

if ( PcheckDecay(genpart,431, -321, 111)==1){
return 539;}//anti-B0 decays to D_s+ K- pi0 

if ( PcheckDecay(genpart,431, -311, -211)==1){
return 540;}//anti-B0 decays to D_s+ anti-K0 pi- 

if ( PcheckDecay(genpart,433, -321)==1){
return 541;}//anti-B0 decays to D_s*+ K- 

if ( PcheckDecay(genpart,433, -323)==1){
return 542;}//anti-B0 decays to D_s*+ K*- 

if ( PcheckDecay(genpart,433, -321, 111)==1){
return 543;}//anti-B0 decays to D_s*+ K- pi0 

if ( PcheckDecay(genpart,433, -311, -211)==1){
return 544;}//anti-B0 decays to D_s*+ anti-K0 pi- 

if ( PcheckDecay(genpart,443, 421)==1){
return 545;}//anti-B0 decays to psi D0 

if ( PcheckDecay(genpart,443, 423)==1){
return 546;}//anti-B0 decays to psi D*0 

if ( PcheckDecay(genpart,443, 411, -211)==1){
return 547;}//anti-B0 decays to psi D+ pi- 

if ( PcheckDecay(genpart,443, 421, 111)==1){
return 548;}//anti-B0 decays to psi D0 pi0 

if ( PcheckDecay(genpart,100443, 421)==1){
return 549;}//anti-B0 decays to psi(2S) D0 

if ( PcheckDecay(genpart,100443, 423)==1){
return 550;}//anti-B0 decays to psi(2S) D*0 

if ( PcheckDecay(genpart,100443, 411, -211)==1){
return 551;}//anti-B0 decays to psi(2S) D+ pi- 

if ( PcheckDecay(genpart,100443, 421, 111)==1){
return 552;}//anti-B0 decays to psi(2S) D0 pi0 

if ( PcheckDecay(genpart,423, 22)==1){
return 553;}//anti-B0 decays to D*0 gamma 

if ( PcheckDecay(genpart,443, 22)==1){
return 554;}//anti-B0 decays to psi gamma 

if ( PcheckDecay(genpart,100443, 22)==1){
return 555;}//anti-B0 decays to psi(2S) gamma 

if ( PcheckDecay(genpart,-15, 15)==1){
return 556;}//anti-B0 decays to tau+ tau- 

if ( PcheckDecay(genpart,22, 22)==1){
return 557;}//anti-B0 decays to gamma gamma 

if ( PcheckDecay(genpart,421, -421)==1){
return 558;}//anti-B0 decays to D0 anti-D0 

if ( PcheckDecay(genpart,423, -421)==1){
return 559;}//anti-B0 decays to D*0 anti-D0 

if ( PcheckDecay(genpart,-423, 421)==1){
return 560;}//anti-B0 decays to anti-D*0 D0 

if ( PcheckDecay(genpart,423, -423)==1){
return 561;}//anti-B0 decays to D*0 anti-D*0 

if ( PcheckDecay(genpart,431, -431)==1){
return 562;}//anti-B0 decays to D_s+ D_s- 

if ( PcheckDecay(genpart,433, -431)==1){
return 563;}//anti-B0 decays to D_s*+ D_s- 

if ( PcheckDecay(genpart,-433, 431)==1){
return 564;}//anti-B0 decays to D_s*- D_s+ 

if ( PcheckDecay(genpart,433, -433)==1){
return 565;}//anti-B0 decays to D_s*+ D_s*- 
 return -100;
      }
  else if (AmId(genpart)==511){

if ( PcheckDecay(genpart,313, 22)==1){
return 566;}//B0 decays to K*0 gamma 

if ( PcheckDecay(genpart,10313, 22)==1){
return 567;}//B0 decays to K_10 gamma 

if ( PcheckDecay(genpart,315, 22)==1){
return 568;}//B0 decays to K_2*0 gamma 

if ( PcheckDecay(genpart,30343, 22)==1){
return 569;}//B0 decays to Xsd gamma 

if ( PcheckDecay(genpart,113, 22)==1){
return 570;}//B0 decays to rho0 gamma 

if ( PcheckDecay(genpart,223, 22)==1){
return 571;}//B0 decays to omega gamma 

if ( PcheckDecay(genpart,30643, 22)==1){
return 572;}//B0 decays to Xdd gamma 

if ( PcheckDecay(genpart,311, -11, 11)==1){
return 573;}//B0 decays to K0 e+ e- 

if ( PcheckDecay(genpart,313, -11, 11)==1){
return 574;}//B0 decays to K*0 e+ e- 

if ( PcheckDecay(genpart,30343, -11, 11)==1){
return 575;}//B0 decays to Xsd e+ e- 

if ( PcheckDecay(genpart,311, -13, 13)==1){
return 576;}//B0 decays to K0 mu+ mu- 

if ( PcheckDecay(genpart,313, -13, 13)==1){
return 577;}//B0 decays to K*0 mu+ mu- 

if ( PcheckDecay(genpart,30343, -13, 13)==1){
return 578;}//B0 decays to Xsd mu+ mu- 

if ( PcheckDecay(genpart,311, -15, 15)==1){
return 579;}//B0 decays to K0 tau+ tau- 

if ( PcheckDecay(genpart,313, -15, 15)==1){
return 580;}//B0 decays to K*0 tau+ tau- 

if ( PcheckDecay(genpart,30343, -15, 15)==1){
return 581;}//B0 decays to Xsd tau+ tau- 

if ( PcheckDecay(genpart,311, 12, -12)==1){
return 582;}//B0 decays to K0 nu_e anti-nu_e 

if ( PcheckDecay(genpart,313, 12, -12)==1){
return 583;}//B0 decays to K*0 nu_e anti-nu_e 

if ( PcheckDecay(genpart,30343, 12, -12)==1){
return 584;}//B0 decays to Xsd nu_e anti-nu_e 

if ( PcheckDecay(genpart,311, 14, -14)==1){
return 585;}//B0 decays to K0 nu_mu anti-nu_mu 

if ( PcheckDecay(genpart,313, 14, -14)==1){
return 586;}//B0 decays to K*0 nu_mu anti-nu_mu 

if ( PcheckDecay(genpart,30343, 14, -14)==1){
return 587;}//B0 decays to Xsd nu_mu anti-nu_mu 

if ( PcheckDecay(genpart,311, 16, -16)==1){
return 588;}//B0 decays to K0 nu_tau anti-nu_tau 

if ( PcheckDecay(genpart,313, 16, -16)==1){
return 589;}//B0 decays to K*0 nu_tau anti-nu_tau 

if ( PcheckDecay(genpart,30343, 16, -16)==1){
return 590;}//B0 decays to Xsd nu_tau anti-nu_tau 

if ( PcheckDecay(genpart,321, -321)==1){
return 591;}//B0 decays to K+ K- 

if ( PcheckDecay(genpart,321, -211)==1){
return 592;}//B0 decays to K+ pi- 

if ( PcheckDecay(genpart,311, 111)==1){
return 593;}//B0 decays to K0 pi0 

if ( PcheckDecay(genpart,130, 130)==1){
return 594;}//B0 decays to K_L0 K_L0 

if ( PcheckDecay(genpart,310, 310)==1){
return 595;}//B0 decays to K_S0 K_S0 

if ( PcheckDecay(genpart,211, -211)==1){
return 596;}//B0 decays to pi+ pi- 

if ( PcheckDecay(genpart,111, 111)==1){
return 597;}//B0 decays to pi0 pi0 

if ( PcheckDecay(genpart,100323, -211)==1){
return 598;}//B0 decays to K'*+ pi- 

if ( PcheckDecay(genpart,100313, 111)==1){
return 599;}//B0 decays to K'*0 pi0 

if ( PcheckDecay(genpart,10321, -211)==1){
return 600;}//B0 decays to K_0*+ pi- 

if ( PcheckDecay(genpart,10211, -211)==1){
return 601;}//B0 decays to a_0+ pi- 

if ( PcheckDecay(genpart,10331, 311)==1){
return 602;}//B0 decays to f'_0 K0 

if ( PcheckDecay(genpart,10221, 311)==1){
return 603;}//B0 decays to f_0 K0 

if ( PcheckDecay(genpart,9030221, 311)==1){
return 604;}//B0 decays to f_0(1500) K0 

if ( PcheckDecay(genpart,10221, 10221)==1){
return 605;}//B0 decays to f_0 f_0 

if ( PcheckDecay(genpart,10111, 311)==1){
return 606;}//B0 decays to a_00 K0 

if ( PcheckDecay(genpart,-10211, 321)==1){
return 607;}//B0 decays to a_0- K+ 

if ( PcheckDecay(genpart,225, 311)==1){
return 608;}//B0 decays to f_2 K0 

if ( PcheckDecay(genpart,325, -211)==1){
return 609;}//B0 decays to K_2*+ pi- 

if ( PcheckDecay(genpart,315, 111)==1){
return 610;}//B0 decays to K_2*0 pi0 

if ( PcheckDecay(genpart,323, -321)==1){
return 611;}//B0 decays to K*+ K- 

if ( PcheckDecay(genpart,323, -211)==1){
return 612;}//B0 decays to K*+ pi- 

if ( PcheckDecay(genpart,-323, 321)==1){
return 613;}//B0 decays to K*- K+ 

if ( PcheckDecay(genpart,313, -311)==1){
return 614;}//B0 decays to K*0 anti-K0 

if ( PcheckDecay(genpart,-313, 311)==1){
return 615;}//B0 decays to anti-K*0 K0 

if ( PcheckDecay(genpart,313, 111)==1){
return 616;}//B0 decays to K*0 pi0 

if ( PcheckDecay(genpart,313, 10221)==1){
return 617;}//B0 decays to K*0 f_0 

if ( PcheckDecay(genpart,100113, 311)==1){
return 618;}//B0 decays to rho(2S)0 K0 

if ( PcheckDecay(genpart,30113, 311)==1){
return 619;}//B0 decays to rho(3S)0 K0 

if ( PcheckDecay(genpart,213, -211)==1){
return 620;}//B0 decays to rho+ pi- 

if ( PcheckDecay(genpart,-213, 321)==1){
return 621;}//B0 decays to rho- K+ 

if ( PcheckDecay(genpart,-213, 211)==1){
return 622;}//B0 decays to rho- pi+ 

if ( PcheckDecay(genpart,113, 311)==1){
return 623;}//B0 decays to rho0 K0 

if ( PcheckDecay(genpart,113, 111)==1){
return 624;}//B0 decays to rho0 pi0 

if ( PcheckDecay(genpart,113, 10221)==1){
return 625;}//B0 decays to rho0 f_0 

if ( PcheckDecay(genpart,100213, -211)==1){
return 626;}//B0 decays to rho(2S)+ pi- 

if ( PcheckDecay(genpart,30213, -211)==1){
return 627;}//B0 decays to rho(3S)+ pi- 

if ( PcheckDecay(genpart,-100213, 211)==1){
return 628;}//B0 decays to rho(2S)- pi+ 

if ( PcheckDecay(genpart,-30213, 211)==1){
return 629;}//B0 decays to rho(3S)- pi+ 

if ( PcheckDecay(genpart,-100213, 321)==1){
return 630;}//B0 decays to rho(2S)- K+ 

if ( PcheckDecay(genpart,-30213, 321)==1){
return 631;}//B0 decays to rho(3S)- K+ 

if ( PcheckDecay(genpart,10323, -211)==1){
return 632;}//B0 decays to K_1+ pi- 

if ( PcheckDecay(genpart,20323, -211)==1){
return 633;}//B0 decays to K'_1+ pi- 

if ( PcheckDecay(genpart,30323, -211)==1){
return 634;}//B0 decays to K''*+ pi- 

if ( PcheckDecay(genpart,30313, 111)==1){
return 635;}//B0 decays to K''*0 pi0 

if ( PcheckDecay(genpart,-10213, 321)==1){
return 636;}//B0 decays to b_1- K+ 

if ( PcheckDecay(genpart,-10213, 211)==1){
return 637;}//B0 decays to b_1- pi+ 

if ( PcheckDecay(genpart,10213, -211)==1){
return 638;}//B0 decays to b_1+ pi- 

if ( PcheckDecay(genpart,321, -321, 311)==1){
return 639;}//B0 decays to K+ K- K0 

if ( PcheckDecay(genpart,321, -321, 111)==1){
return 640;}//B0 decays to K+ K- pi0 

if ( PcheckDecay(genpart,321, -311, -211)==1){
return 641;}//B0 decays to K+ anti-K0 pi- 

if ( PcheckDecay(genpart,321, -211, 111)==1){
return 642;}//B0 decays to K+ pi- pi0 

if ( PcheckDecay(genpart,311, 211, -211)==1){
return 643;}//B0 decays to K0 pi+ pi- 

if ( PcheckDecay(genpart,310, 310, 310)==1){
return 644;}//B0 decays to K_S0 K_S0 K_S0 

if ( PcheckDecay(genpart,130, 130, 130)==1){
return 645;}//B0 decays to K_L0 K_L0 K_L0 

if ( PcheckDecay(genpart,310, 310, 130)==1){
return 646;}//B0 decays to K_S0 K_S0 K_L0 

if ( PcheckDecay(genpart,310, 130, 130)==1){
return 647;}//B0 decays to K_S0 K_L0 K_L0 

if ( PcheckDecay(genpart,211, -211, 111, 211)==1){
return 648;}//B0 decays to pi+ pi- pi0 pi+ 

if ( PcheckDecay(genpart,10221, -211, 211)==1){
return 649;}//B0 decays to f_0 pi- pi+ 

if ( PcheckDecay(genpart,335, 313)==1){
return 650;}//B0 decays to f'_2 K*0 

if ( PcheckDecay(genpart,225, 313)==1){
return 651;}//B0 decays to f_2 K*0 

if ( PcheckDecay(genpart,-213, 10321)==1){
return 652;}//B0 decays to rho- K_0*+ 

if ( PcheckDecay(genpart,113, 10311)==1){
return 653;}//B0 decays to rho0 K_0*0 

if ( PcheckDecay(genpart,-20213, 321)==1){
return 654;}//B0 decays to a_1- K+ 

if ( PcheckDecay(genpart,20213, -211)==1){
return 655;}//B0 decays to a_1+ pi- 

if ( PcheckDecay(genpart,-20213, 211)==1){
return 656;}//B0 decays to a_1- pi+ 

if ( PcheckDecay(genpart,323, -323)==1){
return 657;}//B0 decays to K*+ K*- 

if ( PcheckDecay(genpart,323, -213)==1){
return 658;}//B0 decays to K*+ rho- 

if ( PcheckDecay(genpart,313, -313)==1){
return 659;}//B0 decays to K*0 anti-K*0 

if ( PcheckDecay(genpart,313, 313)==1){
return 660;}//B0 decays to K*0 K*0 

if ( PcheckDecay(genpart,313, 113)==1){
return 661;}//B0 decays to K*0 rho0 

if ( PcheckDecay(genpart,213, -213)==1){
return 662;}//B0 decays to rho+ rho- 

if ( PcheckDecay(genpart,113, 113)==1){
return 663;}//B0 decays to rho0 rho0 

if ( PcheckDecay(genpart,213, -20213)==1){
return 664;}//B0 decays to rho+ a_1- 

if ( PcheckDecay(genpart,-213, 20213)==1){
return 665;}//B0 decays to rho- a_1+ 

if ( PcheckDecay(genpart,323, -211, 111)==1){
return 666;}//B0 decays to K*+ pi- pi0 

if ( PcheckDecay(genpart,313, 211, -211)==1){
return 667;}//B0 decays to K*0 pi+ pi- 

if ( PcheckDecay(genpart,313, 111, 111)==1){
return 668;}//B0 decays to K*0 pi0 pi0 

if ( PcheckDecay(genpart,-213, 321, 111)==1){
return 669;}//B0 decays to rho- K+ pi0 

if ( PcheckDecay(genpart,-213, 311, 211)==1){
return 670;}//B0 decays to rho- K0 pi+ 

if ( PcheckDecay(genpart,113, 321, -211)==1){
return 671;}//B0 decays to rho0 K+ pi- 

if ( PcheckDecay(genpart,113, 311, 111)==1){
return 672;}//B0 decays to rho0 K0 pi0 

if ( PcheckDecay(genpart,113, 211, -211)==1){
return 673;}//B0 decays to rho0 pi+ pi- 

if ( PcheckDecay(genpart,321, -321, 313)==1){
return 674;}//B0 decays to K+ K- K*0 

if ( PcheckDecay(genpart,311, -311, 313)==1){
return 675;}//B0 decays to K0 anti-K0 K*0 

if ( PcheckDecay(genpart,311, -313, 311)==1){
return 676;}//B0 decays to K0 anti-K*0 K0 

if ( PcheckDecay(genpart,-321, 323, 311)==1){
return 677;}//B0 decays to K- K*+ K0 

if ( PcheckDecay(genpart,-323, 321, 311)==1){
return 678;}//B0 decays to K*- K+ K0 

if ( PcheckDecay(genpart,311, 221)==1){
return 679;}//B0 decays to K0 eta 

if ( PcheckDecay(genpart,311, 331)==1){
return 680;}//B0 decays to K0 eta' 

if ( PcheckDecay(genpart,311, 100221)==1){
return 681;}//B0 decays to K0 eta(2S) 

if ( PcheckDecay(genpart,311, 9020221)==1){
return 682;}//B0 decays to K0 eta(1405) 

if ( PcheckDecay(genpart,10311, 221)==1){
return 683;}//B0 decays to K_0*0 eta 

if ( PcheckDecay(genpart,111, 221)==1){
return 684;}//B0 decays to pi0 eta 

if ( PcheckDecay(genpart,111, 331)==1){
return 685;}//B0 decays to pi0 eta' 

if ( PcheckDecay(genpart,111, 100221)==1){
return 686;}//B0 decays to pi0 eta(2S) 

if ( PcheckDecay(genpart,111, 9020221)==1){
return 687;}//B0 decays to pi0 eta(1405) 

if ( PcheckDecay(genpart,221, 221)==1){
return 688;}//B0 decays to eta eta 

if ( PcheckDecay(genpart,221, 331)==1){
return 689;}//B0 decays to eta eta' 

if ( PcheckDecay(genpart,221, 100221)==1){
return 690;}//B0 decays to eta eta(2S) 

if ( PcheckDecay(genpart,221, 9020221)==1){
return 691;}//B0 decays to eta eta(1405) 

if ( PcheckDecay(genpart,331, 331)==1){
return 692;}//B0 decays to eta' eta' 

if ( PcheckDecay(genpart,331, 100221)==1){
return 693;}//B0 decays to eta' eta(2S) 

if ( PcheckDecay(genpart,331, 9020221)==1){
return 694;}//B0 decays to eta' eta(1405) 

if ( PcheckDecay(genpart,100221, 100221)==1){
return 695;}//B0 decays to eta(2S) eta(2S) 

if ( PcheckDecay(genpart,100221, 9020221)==1){
return 696;}//B0 decays to eta(2S) eta(1405) 

if ( PcheckDecay(genpart,9020221, 9020221)==1){
return 697;}//B0 decays to eta(1405) eta(1405) 

if ( PcheckDecay(genpart,10221, 221)==1){
return 698;}//B0 decays to f_0 eta 

if ( PcheckDecay(genpart,10221, 331)==1){
return 699;}//B0 decays to f_0 eta' 

if ( PcheckDecay(genpart,223, 311)==1){
return 700;}//B0 decays to omega K0 

if ( PcheckDecay(genpart,333, 311)==1){
return 701;}//B0 decays to phi K0 

if ( PcheckDecay(genpart,333, 10311)==1){
return 702;}//B0 decays to phi K_0*0 

if ( PcheckDecay(genpart,333, 30313)==1){
return 703;}//B0 decays to phi K''*0 

if ( PcheckDecay(genpart,333, 317)==1){
return 704;}//B0 decays to phi K_3*0 

if ( PcheckDecay(genpart,333, 319)==1){
return 705;}//B0 decays to phi K_4*0 

if ( PcheckDecay(genpart,223, 111)==1){
return 706;}//B0 decays to omega pi0 

if ( PcheckDecay(genpart,333, 111)==1){
return 707;}//B0 decays to phi pi0 

if ( PcheckDecay(genpart,313, 221)==1){
return 708;}//B0 decays to K*0 eta 

if ( PcheckDecay(genpart,113, 221)==1){
return 709;}//B0 decays to rho0 eta 

if ( PcheckDecay(genpart,223, 221)==1){
return 710;}//B0 decays to omega eta 

if ( PcheckDecay(genpart,333, 221)==1){
return 711;}//B0 decays to phi eta 

if ( PcheckDecay(genpart,313, 331)==1){
return 712;}//B0 decays to K*0 eta' 

if ( PcheckDecay(genpart,113, 331)==1){
return 713;}//B0 decays to rho0 eta' 

if ( PcheckDecay(genpart,223, 331)==1){
return 714;}//B0 decays to omega eta' 

if ( PcheckDecay(genpart,333, 331)==1){
return 715;}//B0 decays to phi eta' 

if ( PcheckDecay(genpart,313, 100221)==1){
return 716;}//B0 decays to K*0 eta(2S) 

if ( PcheckDecay(genpart,113, 100221)==1){
return 717;}//B0 decays to rho0 eta(2S) 

if ( PcheckDecay(genpart,223, 100221)==1){
return 718;}//B0 decays to omega eta(2S) 

if ( PcheckDecay(genpart,333, 100221)==1){
return 719;}//B0 decays to phi eta(2S) 

if ( PcheckDecay(genpart,313, 9020221)==1){
return 720;}//B0 decays to K*0 eta(1405) 

if ( PcheckDecay(genpart,113, 9020221)==1){
return 721;}//B0 decays to rho0 eta(1405) 

if ( PcheckDecay(genpart,223, 9020221)==1){
return 722;}//B0 decays to omega eta(1405) 

if ( PcheckDecay(genpart,333, 9020221)==1){
return 723;}//B0 decays to phi eta(1405) 

if ( PcheckDecay(genpart,223, 10221)==1){
return 724;}//B0 decays to omega f_0 

if ( PcheckDecay(genpart,313, 223)==1){
return 725;}//B0 decays to K*0 omega 

if ( PcheckDecay(genpart,313, 333)==1){
return 726;}//B0 decays to K*0 phi 

if ( PcheckDecay(genpart,113, 223)==1){
return 727;}//B0 decays to rho0 omega 

if ( PcheckDecay(genpart,113, 333)==1){
return 728;}//B0 decays to rho0 phi 

if ( PcheckDecay(genpart,223, 223)==1){
return 729;}//B0 decays to omega omega 

if ( PcheckDecay(genpart,223, 333)==1){
return 730;}//B0 decays to omega phi 

if ( PcheckDecay(genpart,333, 333)==1){
return 731;}//B0 decays to phi phi 

if ( PcheckDecay(genpart,315, 333)==1){
return 732;}//B0 decays to K_2*0 phi 

if ( PcheckDecay(genpart,315, 221)==1){
return 733;}//B0 decays to K_2*0 eta 

if ( PcheckDecay(genpart,30343, 331)==1){
return 734;}//B0 decays to Xsd eta' 

if ( PcheckDecay(genpart,30343, 221)==1){
return 735;}//B0 decays to Xsd eta 

if ( PcheckDecay(genpart,221, 321, -211)==1){
return 736;}//B0 decays to eta K+ pi- 

if ( PcheckDecay(genpart,221, 211, -211)==1){
return 737;}//B0 decays to eta pi+ pi- 

if ( PcheckDecay(genpart,221, 311, 111)==1){
return 738;}//B0 decays to eta K0 pi0 

if ( PcheckDecay(genpart,221, 111, 111)==1){
return 739;}//B0 decays to eta pi0 pi0 

if ( PcheckDecay(genpart,331, 321, -211)==1){
return 740;}//B0 decays to eta' K+ pi- 

if ( PcheckDecay(genpart,331, -211, 211)==1){
return 741;}//B0 decays to eta' pi- pi+ 

if ( PcheckDecay(genpart,331, 311, 111)==1){
return 742;}//B0 decays to eta' K0 pi0 

if ( PcheckDecay(genpart,331, 111, 111)==1){
return 743;}//B0 decays to eta' pi0 pi0 

if ( PcheckDecay(genpart,221, 221, 311)==1){
return 744;}//B0 decays to eta eta K0 

if ( PcheckDecay(genpart,221, 221, 111)==1){
return 745;}//B0 decays to eta eta pi0 

if ( PcheckDecay(genpart,221, 331, 311)==1){
return 746;}//B0 decays to eta eta' K0 

if ( PcheckDecay(genpart,221, 331, 111)==1){
return 747;}//B0 decays to eta eta' pi0 

if ( PcheckDecay(genpart,331, 331, 311)==1){
return 748;}//B0 decays to eta' eta' K0 

if ( PcheckDecay(genpart,331, 331, 111)==1){
return 749;}//B0 decays to eta' eta' pi0 

if ( PcheckDecay(genpart,221, 323, -211)==1){
return 750;}//B0 decays to eta K*+ pi- 

if ( PcheckDecay(genpart,221, 321, -213)==1){
return 751;}//B0 decays to eta K+ rho- 

if ( PcheckDecay(genpart,221, -213, 211)==1){
return 752;}//B0 decays to eta rho- pi+ 

if ( PcheckDecay(genpart,221, -211, 213)==1){
return 753;}//B0 decays to eta pi- rho+ 

if ( PcheckDecay(genpart,221, 313, 111)==1){
return 754;}//B0 decays to eta K*0 pi0 

if ( PcheckDecay(genpart,221, 311, 113)==1){
return 755;}//B0 decays to eta K0 rho0 

if ( PcheckDecay(genpart,221, 113, 111)==1){
return 756;}//B0 decays to eta rho0 pi0 

if ( PcheckDecay(genpart,331, 323, -211)==1){
return 757;}//B0 decays to eta' K*+ pi- 

if ( PcheckDecay(genpart,331, 321, -213)==1){
return 758;}//B0 decays to eta' K+ rho- 

if ( PcheckDecay(genpart,331, -213, 211)==1){
return 759;}//B0 decays to eta' rho- pi+ 

if ( PcheckDecay(genpart,331, -211, 213)==1){
return 760;}//B0 decays to eta' pi- rho+ 

if ( PcheckDecay(genpart,331, 313, 111)==1){
return 761;}//B0 decays to eta' K*0 pi0 

if ( PcheckDecay(genpart,331, 311, 113)==1){
return 762;}//B0 decays to eta' K0 rho0 

if ( PcheckDecay(genpart,331, 113, 111)==1){
return 763;}//B0 decays to eta' rho0 pi0 

if ( PcheckDecay(genpart,221, 221, 313)==1){
return 764;}//B0 decays to eta eta K*0 

if ( PcheckDecay(genpart,221, 221, 113)==1){
return 765;}//B0 decays to eta eta rho0 

if ( PcheckDecay(genpart,221, 331, 313)==1){
return 766;}//B0 decays to eta eta' K*0 

if ( PcheckDecay(genpart,221, 331, 113)==1){
return 767;}//B0 decays to eta eta' rho0 

if ( PcheckDecay(genpart,331, 331, 313)==1){
return 768;}//B0 decays to eta' eta' K*0 

if ( PcheckDecay(genpart,331, 331, 113)==1){
return 769;}//B0 decays to eta' eta' rho0 

if ( PcheckDecay(genpart,223, 321, -211)==1){
return 770;}//B0 decays to omega K+ pi- 

if ( PcheckDecay(genpart,223, -211, 211)==1){
return 771;}//B0 decays to omega pi- pi+ 

if ( PcheckDecay(genpart,223, 311, 111)==1){
return 772;}//B0 decays to omega K0 pi0 

if ( PcheckDecay(genpart,223, 111, 111)==1){
return 773;}//B0 decays to omega pi0 pi0 

if ( PcheckDecay(genpart,333, 321, -211)==1){
return 774;}//B0 decays to phi K+ pi- 

if ( PcheckDecay(genpart,333, -211, 211)==1){
return 775;}//B0 decays to phi pi- pi+ 

if ( PcheckDecay(genpart,333, 311, 111)==1){
return 776;}//B0 decays to phi K0 pi0 

if ( PcheckDecay(genpart,333, 111, 111)==1){
return 777;}//B0 decays to phi pi0 pi0 

if ( PcheckDecay(genpart,223, 221, 311)==1){
return 778;}//B0 decays to omega eta K0 

if ( PcheckDecay(genpart,223, 221, 111)==1){
return 779;}//B0 decays to omega eta pi0 

if ( PcheckDecay(genpart,223, 331, 311)==1){
return 780;}//B0 decays to omega eta' K0 

if ( PcheckDecay(genpart,223, 331, 111)==1){
return 781;}//B0 decays to omega eta' pi0 

if ( PcheckDecay(genpart,333, 221, 311)==1){
return 782;}//B0 decays to phi eta K0 

if ( PcheckDecay(genpart,333, 221, 111)==1){
return 783;}//B0 decays to phi eta pi0 

if ( PcheckDecay(genpart,333, 331, 311)==1){
return 784;}//B0 decays to phi eta' K0 

if ( PcheckDecay(genpart,333, 331, 111)==1){
return 785;}//B0 decays to phi eta' pi0 

if ( PcheckDecay(genpart,333, 333, 311)==1){
return 786;}//B0 decays to phi phi K0 

if ( PcheckDecay(genpart,10111, 111)==1){
return 787;}//B0 decays to a_00 pi0 

if ( PcheckDecay(genpart,10111, 221)==1){
return 788;}//B0 decays to a_00 eta 

if ( PcheckDecay(genpart,10111, 331)==1){
return 789;}//B0 decays to a_00 eta' 

if ( PcheckDecay(genpart,10111, 100221)==1){
return 790;}//B0 decays to a_00 eta(2S) 

if ( PcheckDecay(genpart,10111, 9020221)==1){
return 791;}//B0 decays to a_00 eta(1405) 

if ( PcheckDecay(genpart,10111, 10111)==1){
return 792;}//B0 decays to a_00 a_00 

if ( PcheckDecay(genpart,10211, -10211)==1){
return 793;}//B0 decays to a_0+ a_0- 

if ( PcheckDecay(genpart,10111, 10221)==1){
return 794;}//B0 decays to a_00 f_0 

if ( PcheckDecay(genpart,313, 10111)==1){
return 795;}//B0 decays to K*0 a_00 

if ( PcheckDecay(genpart,323, -10211)==1){
return 796;}//B0 decays to K*+ a_0- 

if ( PcheckDecay(genpart,113, 10111)==1){
return 797;}//B0 decays to rho0 a_00 

if ( PcheckDecay(genpart,-213, 10211)==1){
return 798;}//B0 decays to rho- a_0+ 

if ( PcheckDecay(genpart,213, -10211)==1){
return 799;}//B0 decays to rho+ a_0- 

if ( PcheckDecay(genpart,223, 10111)==1){
return 800;}//B0 decays to omega a_00 

if ( PcheckDecay(genpart,333, 10111)==1){
return 801;}//B0 decays to phi a_00 

if ( PcheckDecay(genpart,20113, 10111)==1){
return 802;}//B0 decays to a_10 a_00 

if ( PcheckDecay(genpart,20213, -10211)==1){
return 803;}//B0 decays to a_1+ a_0- 

if ( PcheckDecay(genpart,-20213, 10211)==1){
return 804;}//B0 decays to a_1- a_0+ 

if ( PcheckDecay(genpart,20223, 10111)==1){
return 805;}//B0 decays to f_1 a_00 

if ( PcheckDecay(genpart,10113, 10111)==1){
return 806;}//B0 decays to b_10 a_00 

if ( PcheckDecay(genpart,10213, -10211)==1){
return 807;}//B0 decays to b_1+ a_0- 

if ( PcheckDecay(genpart,-10213, 10211)==1){
return 808;}//B0 decays to b_1- a_0+ 

if ( PcheckDecay(genpart,10223, 10111)==1){
return 809;}//B0 decays to h_1 a_00 

if ( PcheckDecay(genpart,20113, 311)==1){
return 810;}//B0 decays to a_10 K0 

if ( PcheckDecay(genpart,20113, 111)==1){
return 811;}//B0 decays to a_10 pi0 

if ( PcheckDecay(genpart,20113, 221)==1){
return 812;}//B0 decays to a_10 eta 

if ( PcheckDecay(genpart,20113, 331)==1){
return 813;}//B0 decays to a_10 eta' 

if ( PcheckDecay(genpart,20113, 100221)==1){
return 814;}//B0 decays to a_10 eta(2S) 

if ( PcheckDecay(genpart,20113, 9020221)==1){
return 815;}//B0 decays to a_10 eta(1405) 

if ( PcheckDecay(genpart,20113, 10221)==1){
return 816;}//B0 decays to a_10 f_0 

if ( PcheckDecay(genpart,313, 20113)==1){
return 817;}//B0 decays to K*0 a_10 

if ( PcheckDecay(genpart,323, -20213)==1){
return 818;}//B0 decays to K*+ a_1- 

if ( PcheckDecay(genpart,113, 20113)==1){
return 819;}//B0 decays to rho0 a_10 

if ( PcheckDecay(genpart,223, 20113)==1){
return 820;}//B0 decays to omega a_10 

if ( PcheckDecay(genpart,333, 20113)==1){
return 821;}//B0 decays to phi a_10 

if ( PcheckDecay(genpart,20113, 20113)==1){
return 822;}//B0 decays to a_10 a_10 

if ( PcheckDecay(genpart,20213, -20213)==1){
return 823;}//B0 decays to a_1+ a_1- 

if ( PcheckDecay(genpart,20223, 20113)==1){
return 824;}//B0 decays to f_1 a_10 

if ( PcheckDecay(genpart,10113, 20113)==1){
return 825;}//B0 decays to b_10 a_10 

if ( PcheckDecay(genpart,10213, -20213)==1){
return 826;}//B0 decays to b_1+ a_1- 

if ( PcheckDecay(genpart,-10213, 20213)==1){
return 827;}//B0 decays to b_1- a_1+ 

if ( PcheckDecay(genpart,10223, 20113)==1){
return 828;}//B0 decays to h_1 a_10 

if ( PcheckDecay(genpart,10221, 111)==1){
return 829;}//B0 decays to f_0 pi0 

if ( PcheckDecay(genpart,10221, 100221)==1){
return 830;}//B0 decays to f_0 eta(2S) 

if ( PcheckDecay(genpart,10221, 9020221)==1){
return 831;}//B0 decays to f_0 eta(1405) 

if ( PcheckDecay(genpart,333, 10221)==1){
return 832;}//B0 decays to phi f_0 

if ( PcheckDecay(genpart,20223, 10221)==1){
return 833;}//B0 decays to f_1 f_0 

if ( PcheckDecay(genpart,10113, 10221)==1){
return 834;}//B0 decays to b_10 f_0 

if ( PcheckDecay(genpart,10223, 10221)==1){
return 835;}//B0 decays to h_1 f_0 

if ( PcheckDecay(genpart,20223, 311)==1){
return 836;}//B0 decays to f_1 K0 

if ( PcheckDecay(genpart,20223, 111)==1){
return 837;}//B0 decays to f_1 pi0 

if ( PcheckDecay(genpart,20223, 221)==1){
return 838;}//B0 decays to f_1 eta 

if ( PcheckDecay(genpart,20223, 331)==1){
return 839;}//B0 decays to f_1 eta' 

if ( PcheckDecay(genpart,20223, 100221)==1){
return 840;}//B0 decays to f_1 eta(2S) 

if ( PcheckDecay(genpart,20223, 9020221)==1){
return 841;}//B0 decays to f_1 eta(1405) 

if ( PcheckDecay(genpart,313, 20223)==1){
return 842;}//B0 decays to K*0 f_1 

if ( PcheckDecay(genpart,113, 20223)==1){
return 843;}//B0 decays to rho0 f_1 

if ( PcheckDecay(genpart,223, 20223)==1){
return 844;}//B0 decays to omega f_1 

if ( PcheckDecay(genpart,333, 20223)==1){
return 845;}//B0 decays to phi f_1 

if ( PcheckDecay(genpart,20223, 20223)==1){
return 846;}//B0 decays to f_1 f_1 

if ( PcheckDecay(genpart,10113, 20223)==1){
return 847;}//B0 decays to b_10 f_1 

if ( PcheckDecay(genpart,10223, 20223)==1){
return 848;}//B0 decays to h_1 f_1 

if ( PcheckDecay(genpart,10113, 311)==1){
return 849;}//B0 decays to b_10 K0 

if ( PcheckDecay(genpart,10113, 111)==1){
return 850;}//B0 decays to b_10 pi0 

if ( PcheckDecay(genpart,10113, 221)==1){
return 851;}//B0 decays to b_10 eta 

if ( PcheckDecay(genpart,10113, 331)==1){
return 852;}//B0 decays to b_10 eta' 

if ( PcheckDecay(genpart,10113, 100221)==1){
return 853;}//B0 decays to b_10 eta(2S) 

if ( PcheckDecay(genpart,10113, 9020221)==1){
return 854;}//B0 decays to b_10 eta(1405) 

if ( PcheckDecay(genpart,313, 10113)==1){
return 855;}//B0 decays to K*0 b_10 

if ( PcheckDecay(genpart,323, -10213)==1){
return 856;}//B0 decays to K*+ b_1- 

if ( PcheckDecay(genpart,113, 10113)==1){
return 857;}//B0 decays to rho0 b_10 

if ( PcheckDecay(genpart,-213, 10213)==1){
return 858;}//B0 decays to rho- b_1+ 

if ( PcheckDecay(genpart,213, -10213)==1){
return 859;}//B0 decays to rho+ b_1- 

if ( PcheckDecay(genpart,223, 10113)==1){
return 860;}//B0 decays to omega b_10 

if ( PcheckDecay(genpart,333, 10113)==1){
return 861;}//B0 decays to phi b_10 

if ( PcheckDecay(genpart,10113, 10113)==1){
return 862;}//B0 decays to b_10 b_10 

if ( PcheckDecay(genpart,10213, -10213)==1){
return 863;}//B0 decays to b_1+ b_1- 

if ( PcheckDecay(genpart,10223, 10113)==1){
return 864;}//B0 decays to h_1 b_10 

if ( PcheckDecay(genpart,10223, 311)==1){
return 865;}//B0 decays to h_1 K0 

if ( PcheckDecay(genpart,10223, 111)==1){
return 866;}//B0 decays to h_1 pi0 

if ( PcheckDecay(genpart,10223, 221)==1){
return 867;}//B0 decays to h_1 eta 

if ( PcheckDecay(genpart,10223, 331)==1){
return 868;}//B0 decays to h_1 eta' 

if ( PcheckDecay(genpart,10223, 100221)==1){
return 869;}//B0 decays to h_1 eta(2S) 

if ( PcheckDecay(genpart,10223, 9020221)==1){
return 870;}//B0 decays to h_1 eta(1405) 

if ( PcheckDecay(genpart,313, 10223)==1){
return 871;}//B0 decays to K*0 h_1 

if ( PcheckDecay(genpart,113, 10223)==1){
return 872;}//B0 decays to rho0 h_1 

if ( PcheckDecay(genpart,223, 10223)==1){
return 873;}//B0 decays to omega h_1 

if ( PcheckDecay(genpart,333, 10223)==1){
return 874;}//B0 decays to phi h_1 

if ( PcheckDecay(genpart,10223, 10223)==1){
return 875;}//B0 decays to h_1 h_1 

if ( PcheckDecay(genpart,215, -211)==1){
return 876;}//B0 decays to a_2+ pi- 

if ( PcheckDecay(genpart,-215, 211)==1){
return 877;}//B0 decays to a_2- pi+ 

if ( PcheckDecay(genpart,115, 111)==1){
return 878;}//B0 decays to a_20 pi0 

if ( PcheckDecay(genpart,225, 111)==1){
return 879;}//B0 decays to f_2 pi0 

if ( PcheckDecay(genpart,335, 111)==1){
return 880;}//B0 decays to f'_2 pi0 

if ( PcheckDecay(genpart,115, 221)==1){
return 881;}//B0 decays to a_20 eta 

if ( PcheckDecay(genpart,225, 221)==1){
return 882;}//B0 decays to f_2 eta 

if ( PcheckDecay(genpart,335, 221)==1){
return 883;}//B0 decays to f'_2 eta 

if ( PcheckDecay(genpart,115, 331)==1){
return 884;}//B0 decays to a_20 eta' 

if ( PcheckDecay(genpart,225, 331)==1){
return 885;}//B0 decays to f_2 eta' 

if ( PcheckDecay(genpart,335, 331)==1){
return 886;}//B0 decays to f'_2 eta' 

if ( PcheckDecay(genpart,115, 100221)==1){
return 887;}//B0 decays to a_20 eta(2S) 

if ( PcheckDecay(genpart,225, 100221)==1){
return 888;}//B0 decays to f_2 eta(2S) 

if ( PcheckDecay(genpart,335, 100221)==1){
return 889;}//B0 decays to f'_2 eta(2S) 

if ( PcheckDecay(genpart,115, 9020221)==1){
return 890;}//B0 decays to a_20 eta(1405) 

if ( PcheckDecay(genpart,225, 9020221)==1){
return 891;}//B0 decays to f_2 eta(1405) 

if ( PcheckDecay(genpart,335, 9020221)==1){
return 892;}//B0 decays to f'_2 eta(1405) 

if ( PcheckDecay(genpart,315, -311)==1){
return 893;}//B0 decays to K_2*0 anti-K0 

if ( PcheckDecay(genpart,-315, 311)==1){
return 894;}//B0 decays to anti-K_2*0 K0 

if ( PcheckDecay(genpart,-215, 321)==1){
return 895;}//B0 decays to a_2- K+ 

if ( PcheckDecay(genpart,115, 311)==1){
return 896;}//B0 decays to a_20 K0 

if ( PcheckDecay(genpart,335, 311)==1){
return 897;}//B0 decays to f'_2 K0 

if ( PcheckDecay(genpart,315, 331)==1){
return 898;}//B0 decays to K_2*0 eta' 

if ( PcheckDecay(genpart,315, 100221)==1){
return 899;}//B0 decays to K_2*0 eta(2S) 

if ( PcheckDecay(genpart,315, 9020221)==1){
return 900;}//B0 decays to K_2*0 eta(1405) 

if ( PcheckDecay(genpart,215, -213)==1){
return 901;}//B0 decays to a_2+ rho- 

if ( PcheckDecay(genpart,-215, 213)==1){
return 902;}//B0 decays to a_2- rho+ 

if ( PcheckDecay(genpart,115, 113)==1){
return 903;}//B0 decays to a_20 rho0 

if ( PcheckDecay(genpart,225, 113)==1){
return 904;}//B0 decays to f_2 rho0 

if ( PcheckDecay(genpart,335, 113)==1){
return 905;}//B0 decays to f'_2 rho0 

if ( PcheckDecay(genpart,115, 223)==1){
return 906;}//B0 decays to a_20 omega 

if ( PcheckDecay(genpart,225, 223)==1){
return 907;}//B0 decays to f_2 omega 

if ( PcheckDecay(genpart,335, 223)==1){
return 908;}//B0 decays to f'_2 omega 

if ( PcheckDecay(genpart,115, 333)==1){
return 909;}//B0 decays to a_20 phi 

if ( PcheckDecay(genpart,225, 333)==1){
return 910;}//B0 decays to f_2 phi 

if ( PcheckDecay(genpart,335, 333)==1){
return 911;}//B0 decays to f'_2 phi 

if ( PcheckDecay(genpart,315, -313)==1){
return 912;}//B0 decays to K_2*0 anti-K*0 

if ( PcheckDecay(genpart,-315, 313)==1){
return 913;}//B0 decays to anti-K_2*0 K*0 

if ( PcheckDecay(genpart,-215, 323)==1){
return 914;}//B0 decays to a_2- K*+ 

if ( PcheckDecay(genpart,115, 313)==1){
return 915;}//B0 decays to a_20 K*0 

if ( PcheckDecay(genpart,315, 113)==1){
return 916;}//B0 decays to K_2*0 rho0 

if ( PcheckDecay(genpart,325, -213)==1){
return 917;}//B0 decays to K_2*+ rho- 

if ( PcheckDecay(genpart,315, 223)==1){
return 918;}//B0 decays to K_2*0 omega 

if ( PcheckDecay(genpart,-3222, 2212, 111)==1){
return 919;}//B0 decays to anti-Sigma- p+ pi0 

if ( PcheckDecay(genpart,-3122, 2112, 111)==1){
return 920;}//B0 decays to anti-Lambda0 n0 pi0 

if ( PcheckDecay(genpart,-3212, 2112, 111)==1){
return 921;}//B0 decays to anti-Sigma0 n0 pi0 

if ( PcheckDecay(genpart,-3112, 1114, 111)==1){
return 922;}//B0 decays to anti-Sigma+ Delta- pi0 

if ( PcheckDecay(genpart,-3122, 2212, -211)==1){
return 923;}//B0 decays to anti-Lambda0 p+ pi- 

if ( PcheckDecay(genpart,-3212, 2212, -211)==1){
return 924;}//B0 decays to anti-Sigma0 p+ pi- 

if ( PcheckDecay(genpart,-3112, 2112, -211)==1){
return 925;}//B0 decays to anti-Sigma+ n0 pi- 

if ( PcheckDecay(genpart,-3222, 2112, 211)==1){
return 926;}//B0 decays to anti-Sigma- n0 pi+ 

if ( PcheckDecay(genpart,-3122, 1114, 211)==1){
return 927;}//B0 decays to anti-Lambda0 Delta- pi+ 

if ( PcheckDecay(genpart,-3212, 1114, 211)==1){
return 928;}//B0 decays to anti-Sigma0 Delta- pi+ 

if ( PcheckDecay(genpart,-2212, 2212, 311)==1){
return 929;}//B0 decays to anti-p- p+ K0 

if ( PcheckDecay(genpart,-2112, 2112, 311)==1){
return 930;}//B0 decays to anti-n0 n0 K0 

if ( PcheckDecay(genpart,-1114, 1114, 311)==1){
return 931;}//B0 decays to anti-Delta+ Delta- K0 

if ( PcheckDecay(genpart,-2224, 2212, 321)==1){
return 932;}//B0 decays to anti-Delta-- p+ K+ 

if ( PcheckDecay(genpart,-2212, 2112, 321)==1){
return 933;}//B0 decays to anti-p- n0 K+ 

if ( PcheckDecay(genpart,-2112, 1114, 321)==1){
return 934;}//B0 decays to anti-n0 Delta- K+ 

if ( PcheckDecay(genpart,-3322, 3122, 111)==1){
return 935;}//B0 decays to anti-Xi0 Lambda0 pi0 

if ( PcheckDecay(genpart,-3322, 3212, 111)==1){
return 936;}//B0 decays to anti-Xi0 Sigma0 pi0 

if ( PcheckDecay(genpart,-3312, 3112, 111)==1){
return 937;}//B0 decays to anti-Xi+ Sigma- pi0 

if ( PcheckDecay(genpart,-3312, 3122, -211)==1){
return 938;}//B0 decays to anti-Xi+ Lambda0 pi- 

if ( PcheckDecay(genpart,-3312, 3212, -211)==1){
return 939;}//B0 decays to anti-Xi+ Sigma0 pi- 

if ( PcheckDecay(genpart,-3322, 3112, 211)==1){
return 940;}//B0 decays to anti-Xi0 Sigma- pi+ 

if ( PcheckDecay(genpart,-3322, 3222, -211)==1){
return 941;}//B0 decays to anti-Xi0 Sigma+ pi- 

if ( PcheckDecay(genpart,-3322, 2112, -311)==1){
return 942;}//B0 decays to anti-Xi0 n0 anti-K0 

if ( PcheckDecay(genpart,-3322, 2212, -321)==1){
return 943;}//B0 decays to anti-Xi0 p+ K- 

if ( PcheckDecay(genpart,-3312, 1114, -311)==1){
return 944;}//B0 decays to anti-Xi+ Delta- anti-K0 

if ( PcheckDecay(genpart,-3312, 2112, -321)==1){
return 945;}//B0 decays to anti-Xi+ n0 K- 

if ( PcheckDecay(genpart,-3122, 3122, 311)==1){
return 946;}//B0 decays to anti-Lambda0 Lambda0 K0 

if ( PcheckDecay(genpart,-3122, 3212, 311)==1){
return 947;}//B0 decays to anti-Lambda0 Sigma0 K0 

if ( PcheckDecay(genpart,-3212, 3122, 311)==1){
return 948;}//B0 decays to anti-Sigma0 Lambda0 K0 

if ( PcheckDecay(genpart,-3212, 3212, 311)==1){
return 949;}//B0 decays to anti-Sigma0 Sigma0 K0 

if ( PcheckDecay(genpart,-3222, 3122, 321)==1){
return 950;}//B0 decays to anti-Sigma- Lambda0 K+ 

if ( PcheckDecay(genpart,-3222, 3212, 321)==1){
return 951;}//B0 decays to anti-Sigma- Sigma0 K+ 

if ( PcheckDecay(genpart,-3112, 3112, 311)==1){
return 952;}//B0 decays to anti-Sigma+ Sigma- K0 

if ( PcheckDecay(genpart,-3122, 3112, 321)==1){
return 953;}//B0 decays to anti-Lambda0 Sigma- K+ 

if ( PcheckDecay(genpart,-3212, 3112, 321)==1){
return 954;}//B0 decays to anti-Sigma0 Sigma- K+ 

if ( PcheckDecay(genpart,-3334, 3312, 111)==1){
return 955;}//B0 decays to anti-Omega+ Xi- pi0 

if ( PcheckDecay(genpart,-3334, 3322, -211)==1){
return 956;}//B0 decays to anti-Omega+ Xi0 pi- 

if ( PcheckDecay(genpart,-3334, 3112, -311)==1){
return 957;}//B0 decays to anti-Omega+ Sigma- anti-K0 

if ( PcheckDecay(genpart,-3334, 3122, -321)==1){
return 958;}//B0 decays to anti-Omega+ Lambda0 K- 

if ( PcheckDecay(genpart,-3334, 3212, -321)==1){
return 959;}//B0 decays to anti-Omega+ Sigma0 K- 

if ( PcheckDecay(genpart,-3312, 3312, -311)==1){
return 960;}//B0 decays to anti-Xi+ Xi- anti-K0 

if ( PcheckDecay(genpart,-3322, 3312, 321)==1){
return 961;}//B0 decays to anti-Xi0 Xi- K+ 

if ( PcheckDecay(genpart,-3222, 2212, 113)==1){
return 962;}//B0 decays to anti-Sigma- p+ rho0 

if ( PcheckDecay(genpart,-3122, 2112, 113)==1){
return 963;}//B0 decays to anti-Lambda0 n0 rho0 

if ( PcheckDecay(genpart,-3212, 2112, 113)==1){
return 964;}//B0 decays to anti-Sigma0 n0 rho0 

if ( PcheckDecay(genpart,-3112, 1114, 113)==1){
return 965;}//B0 decays to anti-Sigma+ Delta- rho0 

if ( PcheckDecay(genpart,-3122, 2212, -213)==1){
return 966;}//B0 decays to anti-Lambda0 p+ rho- 

if ( PcheckDecay(genpart,-3212, 2212, -213)==1){
return 967;}//B0 decays to anti-Sigma0 p+ rho- 

if ( PcheckDecay(genpart,-3112, 2112, -213)==1){
return 968;}//B0 decays to anti-Sigma+ n0 rho- 

if ( PcheckDecay(genpart,-3222, 2112, 213)==1){
return 969;}//B0 decays to anti-Sigma- n0 rho+ 

if ( PcheckDecay(genpart,-3122, 1114, 213)==1){
return 970;}//B0 decays to anti-Lambda0 Delta- rho+ 

if ( PcheckDecay(genpart,-3212, 1114, 213)==1){
return 971;}//B0 decays to anti-Sigma0 Delta- rho+ 

if ( PcheckDecay(genpart,-2212, 2212, 313)==1){
return 972;}//B0 decays to anti-p- p+ K*0 

if ( PcheckDecay(genpart,-2112, 2112, 313)==1){
return 973;}//B0 decays to anti-n0 n0 K*0 

if ( PcheckDecay(genpart,-1114, 1114, 313)==1){
return 974;}//B0 decays to anti-Delta+ Delta- K*0 

if ( PcheckDecay(genpart,-2224, 2212, 323)==1){
return 975;}//B0 decays to anti-Delta-- p+ K*+ 

if ( PcheckDecay(genpart,-2212, 2112, 323)==1){
return 976;}//B0 decays to anti-p- n0 K*+ 

if ( PcheckDecay(genpart,-2112, 1114, 323)==1){
return 977;}//B0 decays to anti-n0 Delta- K*+ 

if ( PcheckDecay(genpart,-3322, 3122, 113)==1){
return 978;}//B0 decays to anti-Xi0 Lambda0 rho0 

if ( PcheckDecay(genpart,-3322, 3212, 113)==1){
return 979;}//B0 decays to anti-Xi0 Sigma0 rho0 

if ( PcheckDecay(genpart,-3312, 3112, 113)==1){
return 980;}//B0 decays to anti-Xi+ Sigma- rho0 

if ( PcheckDecay(genpart,-3312, 3122, -213)==1){
return 981;}//B0 decays to anti-Xi+ Lambda0 rho- 

if ( PcheckDecay(genpart,-3312, 3212, -213)==1){
return 982;}//B0 decays to anti-Xi+ Sigma0 rho- 

if ( PcheckDecay(genpart,-3322, 3112, 213)==1){
return 983;}//B0 decays to anti-Xi0 Sigma- rho+ 

if ( PcheckDecay(genpart,-3322, 3222, -213)==1){
return 984;}//B0 decays to anti-Xi0 Sigma+ rho- 

if ( PcheckDecay(genpart,-3322, 2112, -313)==1){
return 985;}//B0 decays to anti-Xi0 n0 anti-K*0 

if ( PcheckDecay(genpart,-3322, 2212, -323)==1){
return 986;}//B0 decays to anti-Xi0 p+ K*- 

if ( PcheckDecay(genpart,-3312, 1114, -313)==1){
return 987;}//B0 decays to anti-Xi+ Delta- anti-K*0 

if ( PcheckDecay(genpart,-3312, 2112, -323)==1){
return 988;}//B0 decays to anti-Xi+ n0 K*- 

if ( PcheckDecay(genpart,-3122, 3122, 313)==1){
return 989;}//B0 decays to anti-Lambda0 Lambda0 K*0 

if ( PcheckDecay(genpart,-3122, 3212, 313)==1){
return 990;}//B0 decays to anti-Lambda0 Sigma0 K*0 

if ( PcheckDecay(genpart,-3212, 3122, 313)==1){
return 991;}//B0 decays to anti-Sigma0 Lambda0 K*0 

if ( PcheckDecay(genpart,-3212, 3212, 313)==1){
return 992;}//B0 decays to anti-Sigma0 Sigma0 K*0 

if ( PcheckDecay(genpart,-3222, 3122, 323)==1){
return 993;}//B0 decays to anti-Sigma- Lambda0 K*+ 

if ( PcheckDecay(genpart,-3222, 3212, 323)==1){
return 994;}//B0 decays to anti-Sigma- Sigma0 K*+ 

if ( PcheckDecay(genpart,-3112, 3112, 313)==1){
return 995;}//B0 decays to anti-Sigma+ Sigma- K*0 

if ( PcheckDecay(genpart,-3122, 3112, 323)==1){
return 996;}//B0 decays to anti-Lambda0 Sigma- K*+ 

if ( PcheckDecay(genpart,-3212, 3112, 323)==1){
return 997;}//B0 decays to anti-Sigma0 Sigma- K*+ 

if ( PcheckDecay(genpart,-3334, 3312, 113)==1){
return 998;}//B0 decays to anti-Omega+ Xi- rho0 

if ( PcheckDecay(genpart,-3334, 3112, -313)==1){
return 999;}//B0 decays to anti-Omega+ Sigma- anti-K*0 

if ( PcheckDecay(genpart,-3334, 3122, -323)==1){
return 1000;}//B0 decays to anti-Omega+ Lambda0 K*- 

if ( PcheckDecay(genpart,-3334, 3212, -323)==1){
return 1001;}//B0 decays to anti-Omega+ Sigma0 K*- 

if ( PcheckDecay(genpart,-3312, 3312, -313)==1){
return 1002;}//B0 decays to anti-Xi+ Xi- anti-K*0 

if ( PcheckDecay(genpart,-3322, 3312, 323)==1){
return 1003;}//B0 decays to anti-Xi0 Xi- K*+ 

if ( PcheckDecay(genpart,-2212, 2212, 111)==1){
return 1004;}//B0 decays to anti-p- p+ pi0 

if ( PcheckDecay(genpart,-2112, 2112, 111)==1){
return 1005;}//B0 decays to anti-n0 n0 pi0 

if ( PcheckDecay(genpart,-1114, 1114, 111)==1){
return 1006;}//B0 decays to anti-Delta+ Delta- pi0 

if ( PcheckDecay(genpart,-2212, 2112, 211)==1){
return 1007;}//B0 decays to anti-p- n0 pi+ 

if ( PcheckDecay(genpart,-2212, 2224, -211)==1){
return 1008;}//B0 decays to anti-p- Delta++ pi- 

if ( PcheckDecay(genpart,-2224, 2212, 211)==1){
return 1009;}//B0 decays to anti-Delta-- p+ pi+ 

if ( PcheckDecay(genpart,-2112, 2212, -211)==1){
return 1010;}//B0 decays to anti-n0 p+ pi- 

if ( PcheckDecay(genpart,-1114, 2112, -211)==1){
return 1011;}//B0 decays to anti-Delta+ n0 pi- 

if ( PcheckDecay(genpart,-2112, 1114, 211)==1){
return 1012;}//B0 decays to anti-n0 Delta- pi+ 

if ( PcheckDecay(genpart,-3122, 3122, 111)==1){
return 1013;}//B0 decays to anti-Lambda0 Lambda0 pi0 

if ( PcheckDecay(genpart,-3122, 3212, 111)==1){
return 1014;}//B0 decays to anti-Lambda0 Sigma0 pi0 

if ( PcheckDecay(genpart,-3212, 3122, 111)==1){
return 1015;}//B0 decays to anti-Sigma0 Lambda0 pi0 

if ( PcheckDecay(genpart,-3212, 3212, 111)==1){
return 1016;}//B0 decays to anti-Sigma0 Sigma0 pi0 

if ( PcheckDecay(genpart,-3122, 3112, 211)==1){
return 1017;}//B0 decays to anti-Lambda0 Sigma- pi+ 

if ( PcheckDecay(genpart,-3122, 3222, -211)==1){
return 1018;}//B0 decays to anti-Lambda0 Sigma+ pi- 

if ( PcheckDecay(genpart,-3222, 3122, 211)==1){
return 1019;}//B0 decays to anti-Sigma- Lambda0 pi+ 

if ( PcheckDecay(genpart,-3112, 3122, -211)==1){
return 1020;}//B0 decays to anti-Sigma+ Lambda0 pi- 

if ( PcheckDecay(genpart,-3212, 3112, 211)==1){
return 1021;}//B0 decays to anti-Sigma0 Sigma- pi+ 

if ( PcheckDecay(genpart,-3212, 3222, -211)==1){
return 1022;}//B0 decays to anti-Sigma0 Sigma+ pi- 

if ( PcheckDecay(genpart,-3222, 3212, 211)==1){
return 1023;}//B0 decays to anti-Sigma- Sigma0 pi+ 

if ( PcheckDecay(genpart,-3112, 3212, -211)==1){
return 1024;}//B0 decays to anti-Sigma+ Sigma0 pi- 

if ( PcheckDecay(genpart,-3122, 2212, -321)==1){
return 1025;}//B0 decays to anti-Lambda0 p+ K- 

if ( PcheckDecay(genpart,-3212, 2212, -321)==1){
return 1026;}//B0 decays to anti-Sigma0 p+ K- 

if ( PcheckDecay(genpart,-3122, 2112, -311)==1){
return 1027;}//B0 decays to anti-Lambda0 n0 anti-K0 

if ( PcheckDecay(genpart,-3212, 2112, -311)==1){
return 1028;}//B0 decays to anti-Sigma0 n0 anti-K0 

if ( PcheckDecay(genpart,-2212, 3122, 321)==1){
return 1029;}//B0 decays to anti-p- Lambda0 K+ 

if ( PcheckDecay(genpart,-2212, 3212, 321)==1){
return 1030;}//B0 decays to anti-p- Sigma0 K+ 

if ( PcheckDecay(genpart,-2112, 3122, -311)==1){
return 1031;}//B0 decays to anti-n0 Lambda0 anti-K0 

if ( PcheckDecay(genpart,-2112, 3212, -311)==1){
return 1032;}//B0 decays to anti-n0 Sigma0 anti-K0 

if ( PcheckDecay(genpart,-3312, 3312, 111)==1){
return 1033;}//B0 decays to anti-Xi+ Xi- pi0 

if ( PcheckDecay(genpart,-3322, 3312, 211)==1){
return 1034;}//B0 decays to anti-Xi0 Xi- pi+ 

if ( PcheckDecay(genpart,-3312, 3322, -211)==1){
return 1035;}//B0 decays to anti-Xi+ Xi0 pi- 

if ( PcheckDecay(genpart,-3112, 3312, 311)==1){
return 1036;}//B0 decays to anti-Sigma+ Xi- K0 

if ( PcheckDecay(genpart,-3122, 3312, 321)==1){
return 1037;}//B0 decays to anti-Lambda0 Xi- K+ 

if ( PcheckDecay(genpart,-3212, 3312, 321)==1){
return 1038;}//B0 decays to anti-Sigma0 Xi- K+ 

if ( PcheckDecay(genpart,-3312, 3112, -311)==1){
return 1039;}//B0 decays to anti-Xi+ Sigma- anti-K0 

if ( PcheckDecay(genpart,-3312, 3122, -321)==1){
return 1040;}//B0 decays to anti-Xi+ Lambda0 K- 

if ( PcheckDecay(genpart,-3312, 3212, -321)==1){
return 1041;}//B0 decays to anti-Xi+ Sigma0 K- 

if ( PcheckDecay(genpart,-2212, 2212, 113)==1){
return 1042;}//B0 decays to anti-p- p+ rho0 

if ( PcheckDecay(genpart,-2112, 2112, 113)==1){
return 1043;}//B0 decays to anti-n0 n0 rho0 

if ( PcheckDecay(genpart,-1114, 1114, 113)==1){
return 1044;}//B0 decays to anti-Delta+ Delta- rho0 

if ( PcheckDecay(genpart,-2212, 2112, 213)==1){
return 1045;}//B0 decays to anti-p- n0 rho+ 

if ( PcheckDecay(genpart,-2112, 2212, -213)==1){
return 1046;}//B0 decays to anti-n0 p+ rho- 

if ( PcheckDecay(genpart,-1114, 2112, -213)==1){
return 1047;}//B0 decays to anti-Delta+ n0 rho- 

if ( PcheckDecay(genpart,-2112, 1114, 213)==1){
return 1048;}//B0 decays to anti-n0 Delta- rho+ 

if ( PcheckDecay(genpart,-3122, 3122, 113)==1){
return 1049;}//B0 decays to anti-Lambda0 Lambda0 rho0 

if ( PcheckDecay(genpart,-3122, 3212, 113)==1){
return 1050;}//B0 decays to anti-Lambda0 Sigma0 rho0 

if ( PcheckDecay(genpart,-3212, 3122, 113)==1){
return 1051;}//B0 decays to anti-Sigma0 Lambda0 rho0 

if ( PcheckDecay(genpart,-3212, 3212, 113)==1){
return 1052;}//B0 decays to anti-Sigma0 Sigma0 rho0 

if ( PcheckDecay(genpart,-3122, 3112, 213)==1){
return 1053;}//B0 decays to anti-Lambda0 Sigma- rho+ 

if ( PcheckDecay(genpart,-3122, 3222, -213)==1){
return 1054;}//B0 decays to anti-Lambda0 Sigma+ rho- 

if ( PcheckDecay(genpart,-3222, 3122, 213)==1){
return 1055;}//B0 decays to anti-Sigma- Lambda0 rho+ 

if ( PcheckDecay(genpart,-3112, 3122, -213)==1){
return 1056;}//B0 decays to anti-Sigma+ Lambda0 rho- 

if ( PcheckDecay(genpart,-3212, 3112, 213)==1){
return 1057;}//B0 decays to anti-Sigma0 Sigma- rho+ 

if ( PcheckDecay(genpart,-3212, 3222, -213)==1){
return 1058;}//B0 decays to anti-Sigma0 Sigma+ rho- 

if ( PcheckDecay(genpart,-3222, 3212, 213)==1){
return 1059;}//B0 decays to anti-Sigma- Sigma0 rho+ 

if ( PcheckDecay(genpart,-3112, 3212, -213)==1){
return 1060;}//B0 decays to anti-Sigma+ Sigma0 rho- 

if ( PcheckDecay(genpart,-3122, 2212, -323)==1){
return 1061;}//B0 decays to anti-Lambda0 p+ K*- 

if ( PcheckDecay(genpart,-3212, 2212, -323)==1){
return 1062;}//B0 decays to anti-Sigma0 p+ K*- 

if ( PcheckDecay(genpart,-3122, 2112, -313)==1){
return 1063;}//B0 decays to anti-Lambda0 n0 anti-K*0 

if ( PcheckDecay(genpart,-3212, 2112, -313)==1){
return 1064;}//B0 decays to anti-Sigma0 n0 anti-K*0 

if ( PcheckDecay(genpart,-2212, 3122, 323)==1){
return 1065;}//B0 decays to anti-p- Lambda0 K*+ 

if ( PcheckDecay(genpart,-2212, 3212, 323)==1){
return 1066;}//B0 decays to anti-p- Sigma0 K*+ 

if ( PcheckDecay(genpart,-2112, 3122, 313)==1){
return 1067;}//B0 decays to anti-n0 Lambda0 K*0 

if ( PcheckDecay(genpart,-2112, 3212, 313)==1){
return 1068;}//B0 decays to anti-n0 Sigma0 K*0 

if ( PcheckDecay(genpart,-3312, 3312, 113)==1){
return 1069;}//B0 decays to anti-Xi+ Xi- rho0 

if ( PcheckDecay(genpart,-3322, 3312, 213)==1){
return 1070;}//B0 decays to anti-Xi0 Xi- rho+ 

if ( PcheckDecay(genpart,-3312, 3322, -213)==1){
return 1071;}//B0 decays to anti-Xi+ Xi0 rho- 

if ( PcheckDecay(genpart,-3112, 3312, 313)==1){
return 1072;}//B0 decays to anti-Sigma+ Xi- K*0 

if ( PcheckDecay(genpart,-3122, 3312, 323)==1){
return 1073;}//B0 decays to anti-Lambda0 Xi- K*+ 

if ( PcheckDecay(genpart,-3212, 3312, 323)==1){
return 1074;}//B0 decays to anti-Sigma0 Xi- K*+ 

if ( PcheckDecay(genpart,-3312, 3112, -313)==1){
return 1075;}//B0 decays to anti-Xi+ Sigma- anti-K*0 

if ( PcheckDecay(genpart,-3312, 3122, -323)==1){
return 1076;}//B0 decays to anti-Xi+ Lambda0 K*- 

if ( PcheckDecay(genpart,-3312, 3212, -323)==1){
return 1077;}//B0 decays to anti-Xi+ Sigma0 K*- 

if ( PcheckDecay(genpart,-3222, 2212)==1){
return 1078;}//B0 decays to anti-Sigma- p+ 

if ( PcheckDecay(genpart,-3122, 2112)==1){
return 1079;}//B0 decays to anti-Lambda0 n0 

if ( PcheckDecay(genpart,-3212, 2112)==1){
return 1080;}//B0 decays to anti-Sigma0 n0 

if ( PcheckDecay(genpart,-3112, 1114)==1){
return 1081;}//B0 decays to anti-Sigma+ Delta- 

if ( PcheckDecay(genpart,-3322, 3122)==1){
return 1082;}//B0 decays to anti-Xi0 Lambda0 

if ( PcheckDecay(genpart,-3322, 3212)==1){
return 1083;}//B0 decays to anti-Xi0 Sigma0 

if ( PcheckDecay(genpart,-3312, 3112)==1){
return 1084;}//B0 decays to anti-Xi+ Sigma- 

if ( PcheckDecay(genpart,-3334, 3312)==1){
return 1085;}//B0 decays to anti-Omega+ Xi- 

if ( PcheckDecay(genpart,-2212, 2212)==1){
return 1086;}//B0 decays to anti-p- p+ 

if ( PcheckDecay(genpart,-2112, 2112)==1){
return 1087;}//B0 decays to anti-n0 n0 

if ( PcheckDecay(genpart,-1114, 1114)==1){
return 1088;}//B0 decays to anti-Delta+ Delta- 

if ( PcheckDecay(genpart,-3122, 3122)==1){
return 1089;}//B0 decays to anti-Lambda0 Lambda0 

if ( PcheckDecay(genpart,-3122, 3212)==1){
return 1090;}//B0 decays to anti-Lambda0 Sigma0 

if ( PcheckDecay(genpart,-3212, 3122)==1){
return 1091;}//B0 decays to anti-Sigma0 Lambda0 

if ( PcheckDecay(genpart,-3212, 3212)==1){
return 1092;}//B0 decays to anti-Sigma0 Sigma0 

if ( PcheckDecay(genpart,-3312, 3312)==1){
return 1093;}//B0 decays to anti-Xi+ Xi- 

if ( PcheckDecay(genpart,431, -211)==1){
return 1094;}//B0 decays to D_s+ pi- 

if ( PcheckDecay(genpart,-213, 431)==1){
return 1095;}//B0 decays to rho- D_s+ 

if ( PcheckDecay(genpart,431, -211, 111)==1){
return 1096;}//B0 decays to D_s+ pi- pi0 

if ( PcheckDecay(genpart,-20213, 431)==1){
return 1097;}//B0 decays to a_1- D_s+ 

if ( PcheckDecay(genpart,433, -211)==1){
return 1098;}//B0 decays to D_s*+ pi- 

if ( PcheckDecay(genpart,433, -213)==1){
return 1099;}//B0 decays to D_s*+ rho- 

if ( PcheckDecay(genpart,433, -211, 111)==1){
return 1100;}//B0 decays to D_s*+ pi- pi0 

if ( PcheckDecay(genpart,433, -20213)==1){
return 1101;}//B0 decays to D_s*+ a_1- 

if ( PcheckDecay(genpart,-431, 321)==1){
return 1102;}//B0 decays to D_s- K+ 

if ( PcheckDecay(genpart,323, -431)==1){
return 1103;}//B0 decays to K*+ D_s- 

if ( PcheckDecay(genpart,-431, 321, 111)==1){
return 1104;}//B0 decays to D_s- K+ pi0 

if ( PcheckDecay(genpart,-431, 311, 211)==1){
return 1105;}//B0 decays to D_s- K0 pi+ 

if ( PcheckDecay(genpart,-433, 321)==1){
return 1106;}//B0 decays to D_s*- K+ 

if ( PcheckDecay(genpart,-433, 323)==1){
return 1107;}//B0 decays to D_s*- K*+ 

if ( PcheckDecay(genpart,-433, 321, 111)==1){
return 1108;}//B0 decays to D_s*- K+ pi0 

if ( PcheckDecay(genpart,-433, 311, 211)==1){
return 1109;}//B0 decays to D_s*- K0 pi+ 

if ( PcheckDecay(genpart,443, -421)==1){
return 1110;}//B0 decays to psi anti-D0 

if ( PcheckDecay(genpart,443, -423)==1){
return 1111;}//B0 decays to psi anti-D*0 

if ( PcheckDecay(genpart,443, -411, 211)==1){
return 1112;}//B0 decays to psi D- pi+ 

if ( PcheckDecay(genpart,443, -421, 111)==1){
return 1113;}//B0 decays to psi anti-D0 pi0 

if ( PcheckDecay(genpart,100443, -421)==1){
return 1114;}//B0 decays to psi(2S) anti-D0 

if ( PcheckDecay(genpart,100443, -423)==1){
return 1115;}//B0 decays to psi(2S) anti-D*0 

if ( PcheckDecay(genpart,100443, -411, 211)==1){
return 1116;}//B0 decays to psi(2S) D- pi+ 

if ( PcheckDecay(genpart,100443, -421, 111)==1){
return 1117;}//B0 decays to psi(2S) anti-D0 pi0 

if ( PcheckDecay(genpart,-423, 22)==1){
return 1118;}//B0 decays to anti-D*0 gamma 

if ( PcheckDecay(genpart,443, 22)==1){
return 1119;}//B0 decays to psi gamma 

if ( PcheckDecay(genpart,100443, 22)==1){
return 1120;}//B0 decays to psi(2S) gamma 

if ( PcheckDecay(genpart,-15, 15)==1){
return 1121;}//B0 decays to tau+ tau- 

if ( PcheckDecay(genpart,22, 22)==1){
return 1122;}//B0 decays to gamma gamma 

if ( PcheckDecay(genpart,421, -421)==1){
return 1123;}//B0 decays to D0 anti-D0 

if ( PcheckDecay(genpart,423, -421)==1){
return 1124;}//B0 decays to D*0 anti-D0 

if ( PcheckDecay(genpart,-423, 421)==1){
return 1125;}//B0 decays to anti-D*0 D0 

if ( PcheckDecay(genpart,423, -423)==1){
return 1126;}//B0 decays to D*0 anti-D*0 

if ( PcheckDecay(genpart,431, -431)==1){
return 1127;}//B0 decays to D_s+ D_s- 

if ( PcheckDecay(genpart,433, -431)==1){
return 1128;}//B0 decays to D_s*+ D_s- 

if ( PcheckDecay(genpart,-433, 431)==1){
return 1129;}//B0 decays to D_s*- D_s+ 

if ( PcheckDecay(genpart,433, -433)==1){
return 1130;}//B0 decays to D_s*+ D_s*- 
    return -200;
  }
  else if(AmId(genpart)==521){
    if ( PcheckDecay(genpart,323, 22)==1){
      return 1131;}//B+ decays to K*+ gamma 
  
    if ( PcheckDecay(genpart,10323, 22)==1){
      return 1132;}//B+ decays to K_1+ gamma 
    
    if ( PcheckDecay(genpart,325, 22)==1){
      return 1133;}//B+ decays to K_2*+ gamma 
    
    if ( PcheckDecay(genpart,30353, 22)==1){
      return 1134;}//B+ decays to Xsu gamma 
    
    if ( PcheckDecay(genpart,213, 22)==1){
      return 1135;}//B+ decays to rho+ gamma 
    
    if ( PcheckDecay(genpart,30653, 22)==1){
      return 1136;}//B+ decays to Xdu+ gamma 
    
    if ( PcheckDecay(genpart,321, -11, 11)==1){
      return 1137;}//B+ decays to K+ e+ e- 
    
    if ( PcheckDecay(genpart,323, -11, 11)==1){
      return 1138;}//B+ decays to K*+ e+ e- 
    
    if ( PcheckDecay(genpart,30353, -11, 11)==1){
      return 1139;}//B+ decays to Xsu e+ e- 
    
    if ( PcheckDecay(genpart,321, -13, 13)==1){
      return 1140;}//B+ decays to K+ mu+ mu- 
    
    if ( PcheckDecay(genpart,323, -13, 13)==1){
      return 1141;}//B+ decays to K*+ mu+ mu- 
    
    if ( PcheckDecay(genpart,30353, -13, 13)==1){
      return 1142;}//B+ decays to Xsu mu+ mu- 
    
    if ( PcheckDecay(genpart,321, -15, 15)==1){
      return 1143;}//B+ decays to K+ tau+ tau- 
    
    if ( PcheckDecay(genpart,323, -15, 15)==1){
      return 1144;}//B+ decays to K*+ tau+ tau- 
    
    if ( PcheckDecay(genpart,30353, -15, 15)==1){
      return 1145;}//B+ decays to Xsu tau+ tau- 
    
    if ( PcheckDecay(genpart,321, 12, -12)==1){
      return 1146;}//B+ decays to K+ nu_e anti-nu_e 
    
    if ( PcheckDecay(genpart,323, 12, -12)==1){
      return 1147;}//B+ decays to K*+ nu_e anti-nu_e 
    
    if ( PcheckDecay(genpart,30353, 12, -12)==1){
      return 1148;}//B+ decays to Xsu nu_e anti-nu_e 
    
    if ( PcheckDecay(genpart,321, 14, -14)==1){
      return 1149;}//B+ decays to K+ nu_mu anti-nu_mu 
    
    if ( PcheckDecay(genpart,323, 14, -14)==1){
      return 1150;}//B+ decays to K*+ nu_mu anti-nu_mu 
    
    if ( PcheckDecay(genpart,30353, 14, -14)==1){
      return 1151;}//B+ decays to Xsu nu_mu anti-nu_mu 
    
    if ( PcheckDecay(genpart,321, 16, -16)==1){
      return 1152;}//B+ decays to K+ nu_tau anti-nu_tau 
    
    if ( PcheckDecay(genpart,323, 16, -16)==1){
      return 1153;}//B+ decays to K*+ nu_tau anti-nu_tau 
    
    if ( PcheckDecay(genpart,30353, 16, -16)==1){
      return 1154;}//B+ decays to Xsu nu_tau anti-nu_tau 
    
    if ( PcheckDecay(genpart,321, -311)==1){
      return 1155;}//B+ decays to K+ anti-K0 
    
    if ( PcheckDecay(genpart,321, 111)==1){
      return 1156;}//B+ decays to K+ pi0 
    
    if ( PcheckDecay(genpart,311, 211)==1){
      return 1157;}//B+ decays to K0 pi+ 
    
    if ( PcheckDecay(genpart,211, 111)==1){
      return 1158;}//B+ decays to pi+ pi0 
    
    if ( PcheckDecay(genpart,100323, 111)==1){
      return 1159;}//B+ decays to K'*+ pi0 
    
    if ( PcheckDecay(genpart,100313, 211)==1){
      return 1160;}//B+ decays to K'*0 pi+ 
    
    if ( PcheckDecay(genpart,10311, 211)==1){
      return 1161;}//B+ decays to K_0*0 pi+ 
    
    if ( PcheckDecay(genpart,-10311, 321)==1){
      return 1162;}//B+ decays to anti-K_0*0 K+ 
    
    if ( PcheckDecay(genpart,10331, 321)==1){
      return 1163;}//B+ decays to f'_0 K+ 
    
    if ( PcheckDecay(genpart,10221, 321)==1){
      return 1164;}//B+ decays to f_0 K+ 
    
    if ( PcheckDecay(genpart,9030221, 321)==1){
      return 1165;}//B+ decays to f_0(1500) K+ 
    
    if ( PcheckDecay(genpart,10331, 211)==1){
      return 1166;}//B+ decays to f'_0 pi+ 
  
    if ( PcheckDecay(genpart,10221, 211)==1){
      return 1167;}//B+ decays to f_0 pi+ 
    
    if ( PcheckDecay(genpart,9000221, 211)==1){
      return 1168;}//B+ decays to f_0(600) pi+ 
    
    if ( PcheckDecay(genpart,315, 211)==1){
      return 1169;}//B+ decays to K_2*0 pi+ 
    
    if ( PcheckDecay(genpart,335, 321)==1){
      return 1170;}//B+ decays to f'_2 K+ 
    
    if ( PcheckDecay(genpart,225, 321)==1){
      return 1171;}//B+ decays to f_2 K+ 
    
    if ( PcheckDecay(genpart,10111, 321)==1){
      return 1172;}//B+ decays to a_00 K+ 
    
    if ( PcheckDecay(genpart,10211, 311)==1){
      return 1173;}//B+ decays to a_0+ K0 
    
    if ( PcheckDecay(genpart,10111, 211)==1){
      return 1174;}//B+ decays to a_00 pi+ 
    
    if ( PcheckDecay(genpart,10211, 111)==1){
      return 1175;}//B+ decays to a_0+ pi0 
    
    if ( PcheckDecay(genpart,30313, 211)==1){
      return 1176;}//B+ decays to K''*0 pi+ 
    
    if ( PcheckDecay(genpart,323, -311)==1){
      return 1177;}//B+ decays to K*+ anti-K0 
    
    if ( PcheckDecay(genpart,323, 111)==1){
      return 1178;}//B+ decays to K*+ pi0 
    
    if ( PcheckDecay(genpart,323, 10221)==1){
      return 1179;}//B+ decays to K*+ f_0 
    
    if ( PcheckDecay(genpart,313, 211)==1){
      return 1180;}//B+ decays to K*0 pi+ 
    
    if ( PcheckDecay(genpart,-313, 321)==1){
      return 1181;}//B+ decays to anti-K*0 K+ 
    
    if ( PcheckDecay(genpart,113, 321)==1){
      return 1182;}//B+ decays to rho0 K+ 
    
    if ( PcheckDecay(genpart,100113, 321)==1){
      return 1183;}//B+ decays to rho(2S)0 K+ 
    
    if ( PcheckDecay(genpart,100113, 321)==1){
      return 1184;}//B+ decays to rho(2S)0 K+ 
    
    if ( PcheckDecay(genpart,213, 311)==1){
      return 1185;}//B+ decays to rho+ K0 
    
    if ( PcheckDecay(genpart,100213, 311)==1){
      return 1186;}//B+ decays to rho(2S)+ K0 
    
    if ( PcheckDecay(genpart,30213, 311)==1){
      return 1187;}//B+ decays to rho(3S)+ K0 
    
    if ( PcheckDecay(genpart,213, 111)==1){
      return 1188;}//B+ decays to rho+ pi0 
    
    if ( PcheckDecay(genpart,100213, 111)==1){
      return 1189;}//B+ decays to rho(2S)+ pi0 
    
    if ( PcheckDecay(genpart,30213, 111)==1){
      return 1190;}//B+ decays to rho(3S)+ pi0 
  
    if ( PcheckDecay(genpart,113, 211)==1){
      return 1191;}//B+ decays to rho0 pi+ 
    
    if ( PcheckDecay(genpart,100113, 211)==1){
      return 1192;}//B+ decays to rho(2S)0 pi+ 
    
    if ( PcheckDecay(genpart,30113, 211)==1){
      return 1193;}//B+ decays to rho(3S)0 pi+ 
    
    if ( PcheckDecay(genpart,20113, 321)==1){
      return 1194;}//B+ decays to a_10 K+ 
    
    if ( PcheckDecay(genpart,20213, 311)==1){
      return 1195;}//B+ decays to a_1+ K0 
    
    if ( PcheckDecay(genpart,20113, 211)==1){
      return 1196;}//B+ decays to a_10 pi+ 
    
    if ( PcheckDecay(genpart,20213, 111)==1){
      return 1197;}//B+ decays to a_1+ pi0 
    
    if ( PcheckDecay(genpart,213, 10221)==1){
      return 1198;}//B+ decays to rho+ f_0 
    
    if ( PcheckDecay(genpart,10113, 321)==1){
      return 1199;}//B+ decays to b_10 K+ 
    
    if ( PcheckDecay(genpart,10113, 211)==1){
      return 1200;}//B+ decays to b_10 pi+ 
    
    if ( PcheckDecay(genpart,225, 211)==1){
      return 1201;}//B+ decays to f_2 pi+ 
    
    if ( PcheckDecay(genpart,115, 321)==1){
      return 1202;}//B+ decays to a_20 K+ 
    
    if ( PcheckDecay(genpart,321, 321, -321)==1){
      return 1203;}//B+ decays to K+ K+ K- 
    
    if ( PcheckDecay(genpart,321, -321, 211)==1){
      return 1204;}//B+ decays to K+ K- pi+ 
    
    if ( PcheckDecay(genpart,321, 310, 310)==1){
      return 1205;}//B+ decays to K+ K_S0 K_S0 
    
    if ( PcheckDecay(genpart,321, -311, 111)==1){
      return 1206;}//B+ decays to K+ anti-K0 pi0 
    
    if ( PcheckDecay(genpart,321, 211, -211)==1){
      return 1207;}//B+ decays to K+ pi+ pi- 
    
    if ( PcheckDecay(genpart,311, 211, 111)==1){
      return 1208;}//B+ decays to K0 pi+ pi0 
    
    if ( PcheckDecay(genpart,211, 211, -211)==1){
      return 1209;}//B+ decays to pi+ pi+ pi- 
    
    if ( PcheckDecay(genpart,321, 321, -211)==1){
      return 1210;}//B+ decays to K+ K+ pi- 
    
    if ( PcheckDecay(genpart,-321, 211, 211)==1){
      return 1211;}//B+ decays to K- pi+ pi+ 
    
    if ( PcheckDecay(genpart,310, 310, 211)==1){
      return 1212;}//B+ decays to K_S0 K_S0 pi+ 
    
    if ( PcheckDecay(genpart,130, 130, 211)==1){
      return 1213;}//B+ decays to K_L0 K_L0 pi+ 
    
    if ( PcheckDecay(genpart,310, 130, 211)==1){
      return 1214;}//B+ decays to K_S0 K_L0 pi+ 
    
    if ( PcheckDecay(genpart,335, 323)==1){
      return 1215;}//B+ decays to f'_2 K*+ 
    
    if ( PcheckDecay(genpart,225, 323)==1){
      return 1216;}//B+ decays to f_2 K*+ 
    
    if ( PcheckDecay(genpart,213, 10311)==1){
      return 1217;}//B+ decays to rho+ K_0*0 
    
    if ( PcheckDecay(genpart,113, 10321)==1){
      return 1218;}//B+ decays to rho0 K_0*+ 
    
    if ( PcheckDecay(genpart,323, -313)==1){
      return 1219;}//B+ decays to K*+ anti-K*0 
    
    if ( PcheckDecay(genpart,323, 113)==1){
      return 1220;}//B+ decays to K*+ rho0 
    
    if ( PcheckDecay(genpart,313, 213)==1){
      return 1221;}//B+ decays to K*0 rho+ 
    
    if ( PcheckDecay(genpart,213, 113)==1){
      return 1222;}//B+ decays to rho+ rho0 
    
    if ( PcheckDecay(genpart,323, 211, -211)==1){
      return 1223;}//B+ decays to K*+ pi+ pi- 
    
    if ( PcheckDecay(genpart,323, 111, 111)==1){
      return 1224;}//B+ decays to K*+ pi0 pi0 
    
    if ( PcheckDecay(genpart,313, 211, 111)==1){
      return 1225;}//B+ decays to K*0 pi+ pi0 
    
    if ( PcheckDecay(genpart,213, 321, -211)==1){
      return 1226;}//B+ decays to rho+ K+ pi- 
    
    if ( PcheckDecay(genpart,213, 311, 111)==1){
      return 1227;}//B+ decays to rho+ K0 pi0 
    
    if ( PcheckDecay(genpart,113, 321, 111)==1){
      return 1228;}//B+ decays to rho0 K+ pi0 
    
    if ( PcheckDecay(genpart,113, 311, 211)==1){
      return 1229;}//B+ decays to rho0 K0 pi+ 
    
    if ( PcheckDecay(genpart,321, -321, 323)==1){
      return 1230;}//B+ decays to K+ K- K*+ 
    
    if ( PcheckDecay(genpart,311, -311, 323)==1){
      return 1231;}//B+ decays to K0 anti-K0 K*+ 
    
    if ( PcheckDecay(genpart,311, -313, 321)==1){
      return 1232;}//B+ decays to K0 anti-K*0 K+ 
    
    if ( PcheckDecay(genpart,313, -311, 321)==1){
      return 1233;}//B+ decays to K*0 anti-K0 K+ 
    
    if ( PcheckDecay(genpart,321, -323, 321)==1){
      return 1234;}//B+ decays to K+ K*- K+ 
    
    if ( PcheckDecay(genpart,323, -321, 211)==1){
      return 1235;}//B+ decays to K*+ K- pi+ 
    
    if ( PcheckDecay(genpart,323, 321, -211)==1){
      return 1236;}//B+ decays to K*+ K+ pi- 
    
    if ( PcheckDecay(genpart,-323, 321, 211)==1){
      return 1237;}//B+ decays to K*- K+ pi+ 
  
    if ( PcheckDecay(genpart,221, 321)==1){
      return 1238;}//B+ decays to eta K+ 
    
    if ( PcheckDecay(genpart,331, 321)==1){
      return 1239;}//B+ decays to eta' K+ 
    
    if ( PcheckDecay(genpart,10321, 221)==1){
      return 1240;}//B+ decays to K_0*+ eta 
    
    if ( PcheckDecay(genpart,100221, 321)==1){
      return 1241;}//B+ decays to eta(2S) K+ 
    
    if ( PcheckDecay(genpart,9020221, 321)==1){
      return 1242;}//B+ decays to eta(1405) K+ 
    
    if ( PcheckDecay(genpart,221, 211)==1){
      return 1243;}//B+ decays to eta pi+ 
    
    if ( PcheckDecay(genpart,331, 211)==1){
      return 1244;}//B+ decays to eta' pi+ 
    
    if ( PcheckDecay(genpart,100221, 211)==1){
      return 1245;}//B+ decays to eta(2S) pi+ 
    
    if ( PcheckDecay(genpart,9020221, 211)==1){
      return 1246;}//B+ decays to eta(1405) pi+ 
    
    if ( PcheckDecay(genpart,223, 321)==1){
      return 1247;}//B+ decays to omega K+ 
    
    if ( PcheckDecay(genpart,333, 321)==1){
      return 1248;}//B+ decays to phi K+ 
    
    if ( PcheckDecay(genpart,223, 211)==1){
      return 1249;}//B+ decays to omega pi+ 
    
    if ( PcheckDecay(genpart,333, 211)==1){
      return 1250;}//B+ decays to phi pi+ 
    
    if ( PcheckDecay(genpart,100333, 321)==1){
      return 1251;}//B+ decays to phi(1680) K+ 
    
    if ( PcheckDecay(genpart,323, 221)==1){
      return 1252;}//B+ decays to K*+ eta 
    
    if ( PcheckDecay(genpart,323, 331)==1){
      return 1253;}//B+ decays to K*+ eta' 
    
    if ( PcheckDecay(genpart,323, 100221)==1){
      return 1254;}//B+ decays to K*+ eta(2S) 
    
    if ( PcheckDecay(genpart,323, 9020221)==1){
      return 1255;}//B+ decays to K*+ eta(1405) 
    
    if ( PcheckDecay(genpart,213, 221)==1){
      return 1256;}//B+ decays to rho+ eta 
    
    if ( PcheckDecay(genpart,213, 331)==1){
      return 1257;}//B+ decays to rho+ eta' 
    
    if ( PcheckDecay(genpart,213, 100221)==1){
      return 1258;}//B+ decays to rho+ eta(2S) 
    
    if ( PcheckDecay(genpart,213, 9020221)==1){
      return 1259;}//B+ decays to rho+ eta(1405) 
    
    if ( PcheckDecay(genpart,223, 323)==1){
      return 1260;}//B+ decays to omega K*+ 

    if ( PcheckDecay(genpart,323, 333)==1){
      return 1261;}//B+ decays to K*+ phi 
    
    if ( PcheckDecay(genpart,223, 213)==1){
      return 1262;}//B+ decays to omega rho+ 
    
    if ( PcheckDecay(genpart,333, 213)==1){
      return 1263;}//B+ decays to phi rho+ 
    
    if ( PcheckDecay(genpart,325, 221)==1){
      return 1264;}//B+ decays to K_2*+ eta 
    
    if ( PcheckDecay(genpart,30353, 331)==1){
      return 1265;}//B+ decays to Xsu eta' 
    
    if ( PcheckDecay(genpart,30353, 221)==1){
      return 1266;}//B+ decays to Xsu eta 
    
    //Manual
    if ( PcheckDecay(genpart,221,311,211)==1){return 1267;}//B+ decays to eta K0 pi+
    
    if ( PcheckDecay(genpart,211,321,111)==1){return 1268;}//B+ decays to eta K+pi0
    
    if ( PcheckDecay(genpart,221, 111, 211)==1){
      return 1269;}//B+ decays to eta pi0 pi+ 
    
    if ( PcheckDecay(genpart,331, 311, 211)==1){
      return 1270;}//B+ decays to eta' K0 pi+ 
    
    if ( PcheckDecay(genpart,331, 321, 111)==1){
      return 1271;}//B+ decays to eta' K+ pi0 
    
    if ( PcheckDecay(genpart,331, 111, 211)==1){
      return 1272;}//B+ decays to eta' pi0 pi+ 
    
    if ( PcheckDecay(genpart,221, 221, 321)==1){
      return 1273;}//B+ decays to eta eta K+ 
    
    if ( PcheckDecay(genpart,221, 221, 211)==1){
      return 1274;}//B+ decays to eta eta pi+ 
    
    if ( PcheckDecay(genpart,221, 331, 321)==1){
      return 1275;}//B+ decays to eta eta' K+ 
    
    if ( PcheckDecay(genpart,221, 331, 211)==1){
      return 1276;}//B+ decays to eta eta' pi+ 
    
    if ( PcheckDecay(genpart,331, 331, 321)==1){
      return 1277;}//B+ decays to eta' eta' K+ 
    
    if ( PcheckDecay(genpart,331, 331, 211)==1){
      return 1278;}//B+ decays to eta' eta' pi+ 
    
    if ( PcheckDecay(genpart,221, 313, 211)==1){
      return 1279;}//B+ decays to eta K*0 pi+ 
    
    if ( PcheckDecay(genpart,221, 311, 213)==1){
      return 1280;}//B+ decays to eta K0 rho+ 
    
    if ( PcheckDecay(genpart,221, 323, 111)==1){
      return 1281;}//B+ decays to eta K*+ pi0 
    
    if ( PcheckDecay(genpart,221, 321, 113)==1){
      return 1282;}//B+ decays to eta K+ rho0 
    
    if ( PcheckDecay(genpart,221, 113, 211)==1){
      return 1283;}//B+ decays to eta rho0 pi+ 
    
    if ( PcheckDecay(genpart,221, 111, 213)==1){
      return 1284;}//B+ decays to eta pi0 rho+ 
    
    if ( PcheckDecay(genpart,331, 313, 211)==1){
      return 1285;}//B+ decays to eta' K*0 pi+ 
    
    if ( PcheckDecay(genpart,331, 311, 213)==1){
      return 1286;}//B+ decays to eta' K0 rho+ 
    
    if ( PcheckDecay(genpart,331, 323, 111)==1){
      return 1287;}//B+ decays to eta' K*+ pi0 
    
    if ( PcheckDecay(genpart,331, 321, 113)==1){
      return 1288;}//B+ decays to eta' K+ rho0 
    
    if ( PcheckDecay(genpart,331, 113, 211)==1){
      return 1289;}//B+ decays to eta' rho0 pi+ 
    
    if ( PcheckDecay(genpart,331, 111, 213)==1){
      return 1290;}//B+ decays to eta' pi0 rho+ 
    
    if ( PcheckDecay(genpart,221, 221, 323)==1){
      return 1291;}//B+ decays to eta eta K*+ 
    
    if ( PcheckDecay(genpart,221, 221, 213)==1){
      return 1292;}//B+ decays to eta eta rho+ 
    
    if ( PcheckDecay(genpart,221, 331, 323)==1){
      return 1293;}//B+ decays to eta eta' K*+ 
    
    if ( PcheckDecay(genpart,221, 331, 213)==1){
      return 1294;}//B+ decays to eta eta' rho+ 
    
    if ( PcheckDecay(genpart,331, 331, 323)==1){
      return 1295;}//B+ decays to eta' eta' K*+ 
    
    if ( PcheckDecay(genpart,331, 331, 213)==1){
      return 1296;}//B+ decays to eta' eta' rho+ 
    
    if ( PcheckDecay(genpart,223, 311, 211)==1){
      return 1297;}//B+ decays to omega K0 pi+ 
    
    if ( PcheckDecay(genpart,223, 321, 111)==1){
      return 1298;}//B+ decays to omega K+ pi0 
    
    if ( PcheckDecay(genpart,223, 111, 211)==1){
      return 1299;}//B+ decays to omega pi0 pi+ 
    
    if ( PcheckDecay(genpart,333, 311, 211)==1){
      return 1300;}//B+ decays to phi K0 pi+ 
    
    if ( PcheckDecay(genpart,333, 321, 111)==1){
      return 1301;}//B+ decays to phi K+ pi0 
    
    if ( PcheckDecay(genpart,333, 111, 211)==1){
      return 1302;}//B+ decays to phi pi0 pi+ 
    
    if ( PcheckDecay(genpart,223, 221, 321)==1){
      return 1303;}//B+ decays to omega eta K+ 
    
    if ( PcheckDecay(genpart,223, 221, 211)==1){
      return 1304;}//B+ decays to omega eta pi+ 
    
    if ( PcheckDecay(genpart,223, 331, 321)==1){
      return 1305;}//B+ decays to omega eta' K+ 
    
    if ( PcheckDecay(genpart,223, 331, 211)==1){
      return 1306;}//B+ decays to omega eta' pi+ 
    
    if ( PcheckDecay(genpart,333, 221, 321)==1){
      return 1307;}//B+ decays to phi eta K+ 
    
    if ( PcheckDecay(genpart,333, 221, 211)==1){
      return 1308;}//B+ decays to phi eta pi+ 
    
    if ( PcheckDecay(genpart,333, 331, 321)==1){
      return 1309;}//B+ decays to phi eta' K+ 
    
    if ( PcheckDecay(genpart,333, 331, 211)==1){
      return 1310;}//B+ decays to phi eta' pi+ 
    
    if ( PcheckDecay(genpart,333, 333, 321)==1){
      return 1311;}//B+ decays to phi phi K+ 
    
    if ( PcheckDecay(genpart,10211, 221)==1){
      return 1312;}//B+ decays to a_0+ eta 
    
    if ( PcheckDecay(genpart,10211, 331)==1){
      return 1313;}//B+ decays to a_0+ eta' 
    
    if ( PcheckDecay(genpart,10211, 100221)==1){
      return 1314;}//B+ decays to a_0+ eta(2S) 
    
    if ( PcheckDecay(genpart,10211, 9020221)==1){
      return 1315;}//B+ decays to a_0+ eta(1405) 
    
    if ( PcheckDecay(genpart,10211, 10111)==1){
      return 1316;}//B+ decays to a_0+ a_00 
    
    if ( PcheckDecay(genpart,10211, 10221)==1){
      return 1317;}//B+ decays to a_0+ f_0 
    
    if ( PcheckDecay(genpart,323, 10111)==1){
      return 1318;}//B+ decays to K*+ a_00 
    
    if ( PcheckDecay(genpart,313, 10211)==1){
      return 1319;}//B+ decays to K*0 a_0+ 
    
    if ( PcheckDecay(genpart,213, 10111)==1){
      return 1320;}//B+ decays to rho+ a_00 
    
    if ( PcheckDecay(genpart,113, 10211)==1){
      return 1321;}//B+ decays to rho0 a_0+ 
    
    if ( PcheckDecay(genpart,223, 10211)==1){
      return 1322;}//B+ decays to omega a_0+ 
    
    if ( PcheckDecay(genpart,333, 10211)==1){
      return 1323;}//B+ decays to phi a_0+ 
    
    if ( PcheckDecay(genpart,20113, 10211)==1){
      return 1324;}//B+ decays to a_10 a_0+ 
    
    if ( PcheckDecay(genpart,20213, 10111)==1){
      return 1325;}//B+ decays to a_1+ a_00 
    
    if ( PcheckDecay(genpart,20223, 10211)==1){
      return 1326;}//B+ decays to f_1 a_0+ 
    
    if ( PcheckDecay(genpart,10113, 10211)==1){
      return 1327;}//B+ decays to b_10 a_0+ 
    
    if ( PcheckDecay(genpart,10213, 10111)==1){
      return 1328;}//B+ decays to b_1+ a_00 
    
    if ( PcheckDecay(genpart,10223, 10211)==1){
      return 1329;}//B+ decays to h_1 a_0+ 
    
    if ( PcheckDecay(genpart,20213, 221)==1){
      return 1330;}//B+ decays to a_1+ eta 
    
    if ( PcheckDecay(genpart,20213, 331)==1){
      return 1331;}//B+ decays to a_1+ eta' 
    
    if ( PcheckDecay(genpart,20213, 100221)==1){
      return 1332;}//B+ decays to a_1+ eta(2S) 
    
    if ( PcheckDecay(genpart,20213, 9020221)==1){
      return 1333;}//B+ decays to a_1+ eta(1405) 
    
    if ( PcheckDecay(genpart,20213, 10221)==1){
      return 1334;}//B+ decays to a_1+ f_0 
    
    if ( PcheckDecay(genpart,323, 20113)==1){
      return 1335;}//B+ decays to K*+ a_10 
    
    if ( PcheckDecay(genpart,313, 20213)==1){
      return 1336;}//B+ decays to K*0 a_1+ 
    
    if ( PcheckDecay(genpart,213, 20113)==1){
      return 1337;}//B+ decays to rho+ a_10 
    
    if ( PcheckDecay(genpart,113, 20213)==1){
      return 1338;}//B+ decays to rho0 a_1+ 
    
    if ( PcheckDecay(genpart,223, 20213)==1){
      return 1339;}//B+ decays to omega a_1+ 
    
    if ( PcheckDecay(genpart,333, 20213)==1){
      return 1340;}//B+ decays to phi a_1+ 
    
    if ( PcheckDecay(genpart,20113, 20213)==1){
      return 1341;}//B+ decays to a_10 a_1+ 
    
    if ( PcheckDecay(genpart,20223, 20213)==1){
      return 1342;}//B+ decays to f_1 a_1+ 
    
    if ( PcheckDecay(genpart,10113, 20213)==1){
      return 1343;}//B+ decays to b_10 a_1+ 
    
    if ( PcheckDecay(genpart,10213, 20113)==1){
      return 1344;}//B+ decays to b_1+ a_10 
    
    if ( PcheckDecay(genpart,10223, 20213)==1){
      return 1345;}//B+ decays to h_1 a_1+ 

    if ( PcheckDecay(genpart,10213, 10221)==1){
      return 1346;}//B+ decays to b_1+ f_0 
    
    if ( PcheckDecay(genpart,20223, 321)==1){
      return 1347;}//B+ decays to f_1 K+ 
    
    if ( PcheckDecay(genpart,20223, 211)==1){
return 1348;}//B+ decays to f_1 pi+ 

if ( PcheckDecay(genpart,323, 20223)==1){
return 1349;}//B+ decays to K*+ f_1 

if ( PcheckDecay(genpart,213, 20223)==1){
return 1350;}//B+ decays to rho+ f_1 

if ( PcheckDecay(genpart,10213, 20223)==1){
return 1351;}//B+ decays to b_1+ f_1 

if ( PcheckDecay(genpart,10213, 311)==1){
return 1352;}//B+ decays to b_1+ K0 

if ( PcheckDecay(genpart,10213, 111)==1){
return 1353;}//B+ decays to b_1+ pi0 

if ( PcheckDecay(genpart,10213, 221)==1){
return 1354;}//B+ decays to b_1+ eta 

if ( PcheckDecay(genpart,10213, 331)==1){
return 1355;}//B+ decays to b_1+ eta' 

if ( PcheckDecay(genpart,10213, 100221)==1){
return 1356;}//B+ decays to b_1+ eta(2S) 

if ( PcheckDecay(genpart,10213, 9020221)==1){
return 1357;}//B+ decays to b_1+ eta(1405) 

if ( PcheckDecay(genpart,323, 10113)==1){
return 1358;}//B+ decays to K*+ b_10 

if ( PcheckDecay(genpart,313, 10213)==1){
return 1359;}//B+ decays to K*0 b_1+ 

if ( PcheckDecay(genpart,213, 10113)==1){
return 1360;}//B+ decays to rho+ b_10 

if ( PcheckDecay(genpart,113, 10213)==1){
return 1361;}//B+ decays to rho0 b_1+ 

if ( PcheckDecay(genpart,223, 10213)==1){
return 1362;}//B+ decays to omega b_1+ 

if ( PcheckDecay(genpart,333, 10213)==1){
return 1363;}//B+ decays to phi b_1+ 

if ( PcheckDecay(genpart,10113, 10213)==1){
return 1364;}//B+ decays to b_10 b_1+ 

if ( PcheckDecay(genpart,10223, 10213)==1){
return 1365;}//B+ decays to h_1 b_1+ 

if ( PcheckDecay(genpart,10223, 321)==1){
return 1366;}//B+ decays to h_1 K+ 

if ( PcheckDecay(genpart,10223, 211)==1){
return 1367;}//B+ decays to h_1 pi+ 

if ( PcheckDecay(genpart,323, 10223)==1){
return 1368;}//B+ decays to K*+ h_1 

if ( PcheckDecay(genpart,213, 10223)==1){
return 1369;}//B+ decays to rho+ h_1 

if ( PcheckDecay(genpart,115, 211)==1){
return 1370;}//B+ decays to a_20 pi+ 

if ( PcheckDecay(genpart,335, 211)==1){
return 1371;}//B+ decays to f'_2 pi+ 

if ( PcheckDecay(genpart,215, 111)==1){
return 1372;}//B+ decays to a_2+ pi0 

if ( PcheckDecay(genpart,215, 221)==1){
return 1373;}//B+ decays to a_2+ eta 

if ( PcheckDecay(genpart,215, 331)==1){
return 1374;}//B+ decays to a_2+ eta' 

if ( PcheckDecay(genpart,215, 100221)==1){
return 1375;}//B+ decays to a_2+ eta(2S) 

if ( PcheckDecay(genpart,215, 9020221)==1){
return 1376;}//B+ decays to a_2+ eta(1405) 

if ( PcheckDecay(genpart,325, -311)==1){
return 1377;}//B+ decays to K_2*+ anti-K0 

if ( PcheckDecay(genpart,-315, 321)==1){
return 1378;}//B+ decays to anti-K_2*0 K+ 

if ( PcheckDecay(genpart,215, 311)==1){
return 1379;}//B+ decays to a_2+ K0 

if ( PcheckDecay(genpart,325, 111)==1){
return 1380;}//B+ decays to K_2*+ pi0 

if ( PcheckDecay(genpart,325, 331)==1){
return 1381;}//B+ decays to K_2*+ eta' 

if ( PcheckDecay(genpart,325, 100221)==1){
return 1382;}//B+ decays to K_2*+ eta(2S) 

if ( PcheckDecay(genpart,325, 9020221)==1){
return 1383;}//B+ decays to K_2*+ eta(1405) 

if ( PcheckDecay(genpart,115, 213)==1){
return 1384;}//B+ decays to a_20 rho+ 

if ( PcheckDecay(genpart,225, 213)==1){
return 1385;}//B+ decays to f_2 rho+ 

if ( PcheckDecay(genpart,335, 213)==1){
return 1386;}//B+ decays to f'_2 rho+ 

if ( PcheckDecay(genpart,215, 113)==1){
return 1387;}//B+ decays to a_2+ rho0 

if ( PcheckDecay(genpart,215, 223)==1){
return 1388;}//B+ decays to a_2+ omega 

if ( PcheckDecay(genpart,215, 333)==1){
return 1389;}//B+ decays to a_2+ phi 

if ( PcheckDecay(genpart,325, -313)==1){
return 1390;}//B+ decays to K_2*+ anti-K*0 

if ( PcheckDecay(genpart,-315, 323)==1){
return 1391;}//B+ decays to anti-K_2*0 K*+ 

if ( PcheckDecay(genpart,115, 323)==1){
return 1392;}//B+ decays to a_20 K*+ 

if ( PcheckDecay(genpart,215, 313)==1){
return 1393;}//B+ decays to a_2+ K*0 

if ( PcheckDecay(genpart,325, 113)==1){
return 1394;}//B+ decays to K_2*+ rho0 

if ( PcheckDecay(genpart,325, 223)==1){
return 1395;}//B+ decays to K_2*+ omega 

if ( PcheckDecay(genpart,325, 333)==1){
return 1396;}//B+ decays to K_2*+ phi 

if ( PcheckDecay(genpart,315, 213)==1){
return 1397;}//B+ decays to K_2*0 rho+ 

if ( PcheckDecay(genpart,-3222, 2224, 111)==1){
return 1398;}//B+ decays to anti-Sigma- Delta++ pi0 

if ( PcheckDecay(genpart,-3122, 2212, 111)==1){
return 1399;}//B+ decays to anti-Lambda0 p+ pi0 

if ( PcheckDecay(genpart,-3212, 2212, 111)==1){
return 1400;}//B+ decays to anti-Sigma0 p+ pi0 

if ( PcheckDecay(genpart,-3112, 2112, 111)==1){
return 1401;}//B+ decays to anti-Sigma+ n0 pi0 

if ( PcheckDecay(genpart,-3122, 2224, -211)==1){
return 1402;}//B+ decays to anti-Lambda0 Delta++ pi- 

if ( PcheckDecay(genpart,-3212, 2224, -211)==1){
return 1403;}//B+ decays to anti-Sigma0 Delta++ pi- 

if ( PcheckDecay(genpart,-3112, 2212, -211)==1){
return 1404;}//B+ decays to anti-Sigma+ p+ pi- 

if ( PcheckDecay(genpart,-3222, 2212, 211)==1){
return 1405;}//B+ decays to anti-Sigma- p+ pi+ 

if ( PcheckDecay(genpart,-3122, 2112, 211)==1){
return 1406;}//B+ decays to anti-Lambda0 n0 pi+ 

if ( PcheckDecay(genpart,-3212, 2112, 211)==1){
return 1407;}//B+ decays to anti-Sigma0 n0 pi+ 

if ( PcheckDecay(genpart,-2212, 2224, 311)==1){
return 1408;}//B+ decays to anti-p- Delta++ K0 

if ( PcheckDecay(genpart,-2112, 2212, 311)==1){
return 1409;}//B+ decays to anti-n0 p+ K0 

if ( PcheckDecay(genpart,-1114, 2112, 311)==1){
return 1410;}//B+ decays to anti-Delta+ n0 K0 

if ( PcheckDecay(genpart,-2224, 2224, 321)==1){
return 1411;}//B+ decays to anti-Delta-- Delta++ K+ 

if ( PcheckDecay(genpart,-2212, 2212, 321)==1){
return 1412;}//B+ decays to anti-p- p+ K+ 

if ( PcheckDecay(genpart,-2112, 2112, 321)==1){
return 1413;}//B+ decays to anti-n0 n0 K+ 

if ( PcheckDecay(genpart,-3322, 3222, 111)==1){
return 1414;}//B+ decays to anti-Xi0 Sigma+ pi0 

if ( PcheckDecay(genpart,-3312, 3122, 111)==1){
return 1415;}//B+ decays to anti-Xi+ Lambda0 pi0 

if ( PcheckDecay(genpart,-3312, 3212, 111)==1){
return 1416;}//B+ decays to anti-Xi+ Sigma0 pi0 

if ( PcheckDecay(genpart,-3312, 3222, -211)==1){
return 1417;}//B+ decays to anti-Xi+ Sigma+ pi- 

if ( PcheckDecay(genpart,-3322, 3122, 211)==1){
return 1418;}//B+ decays to anti-Xi0 Lambda0 pi+ 

if ( PcheckDecay(genpart,-3322, 3212, 211)==1){
return 1419;}//B+ decays to anti-Xi0 Sigma0 pi+ 

if ( PcheckDecay(genpart,-3312, 3112, 211)==1){
return 1420;}//B+ decays to anti-Xi+ Sigma- pi+ 

if ( PcheckDecay(genpart,-3322, 2212, -311)==1){
return 1421;}//B+ decays to anti-Xi0 p+ anti-K0 

if ( PcheckDecay(genpart,-3322, 2224, -321)==1){
return 1422;}//B+ decays to anti-Xi0 Delta++ K- 

if ( PcheckDecay(genpart,-3312, 2112, -311)==1){
return 1423;}//B+ decays to anti-Xi+ n0 anti-K0 

if ( PcheckDecay(genpart,-3312, 2212, -321)==1){
return 1424;}//B+ decays to anti-Xi+ p+ K- 

if ( PcheckDecay(genpart,-3122, 3222, 311)==1){
return 1425;}//B+ decays to anti-Lambda0 Sigma+ K0 

if ( PcheckDecay(genpart,-3212, 3222, 311)==1){
return 1426;}//B+ decays to anti-Sigma0 Sigma+ K0 

if ( PcheckDecay(genpart,-3222, 3222, 321)==1){
return 1427;}//B+ decays to anti-Sigma- Sigma+ K+ 

if ( PcheckDecay(genpart,-3112, 3122, 311)==1){
return 1428;}//B+ decays to anti-Sigma+ Lambda0 K0 

if ( PcheckDecay(genpart,-3112, 3212, 311)==1){
return 1429;}//B+ decays to anti-Sigma+ Sigma0 K0 

if ( PcheckDecay(genpart,-3122, 3122, 321)==1){
return 1430;}//B+ decays to anti-Lambda0 Lambda0 K+ 

if ( PcheckDecay(genpart,-3122, 3212, 321)==1){
return 1431;}//B+ decays to anti-Lambda0 Sigma0 K+ 

if ( PcheckDecay(genpart,-3212, 3122, 321)==1){
return 1432;}//B+ decays to anti-Sigma0 Lambda0 K+ 

if ( PcheckDecay(genpart,-3212, 3212, 321)==1){
return 1433;}//B+ decays to anti-Sigma0 Sigma0 K+ 

if ( PcheckDecay(genpart,-3334, 3322, 111)==1){
return 1434;}//B+ decays to anti-Omega+ Xi0 pi0 

if ( PcheckDecay(genpart,-3334, 3312, 211)==1){
return 1435;}//B+ decays to anti-Omega+ Xi- pi+ 

if ( PcheckDecay(genpart,-3334, 3122, -311)==1){
return 1436;}//B+ decays to anti-Omega+ Lambda0 anti-K0 

if ( PcheckDecay(genpart,-3334, 3212, -311)==1){
return 1437;}//B+ decays to anti-Omega+ Sigma0 anti-K0 

if ( PcheckDecay(genpart,-3334, 3222, -321)==1){
return 1438;}//B+ decays to anti-Omega+ Sigma+ K- 

if ( PcheckDecay(genpart,-3312, 3322, 311)==1){
return 1439;}//B+ decays to anti-Xi+ Xi0 K0 

if ( PcheckDecay(genpart,-3322, 3322, 321)==1){
return 1440;}//B+ decays to anti-Xi0 Xi0 K+ 

if ( PcheckDecay(genpart,-3222, 2224, 113)==1){
return 1441;}//B+ decays to anti-Sigma- Delta++ rho0 

if ( PcheckDecay(genpart,-3122, 2212, 113)==1){
return 1442;}//B+ decays to anti-Lambda0 p+ rho0 

if ( PcheckDecay(genpart,-3212, 2212, 113)==1){
return 1443;}//B+ decays to anti-Sigma0 p+ rho0 

if ( PcheckDecay(genpart,-3112, 2112, 113)==1){
return 1444;}//B+ decays to anti-Sigma+ n0 rho0 

if ( PcheckDecay(genpart,-3122, 2224, -213)==1){
return 1445;}//B+ decays to anti-Lambda0 Delta++ rho- 

if ( PcheckDecay(genpart,-3212, 2224, -213)==1){
return 1446;}//B+ decays to anti-Sigma0 Delta++ rho- 

if ( PcheckDecay(genpart,-3112, 2212, -213)==1){
return 1447;}//B+ decays to anti-Sigma+ p+ rho- 

if ( PcheckDecay(genpart,-3222, 2212, 213)==1){
return 1448;}//B+ decays to anti-Sigma- p+ rho+ 

if ( PcheckDecay(genpart,-3122, 2112, 213)==1){
return 1449;}//B+ decays to anti-Lambda0 n0 rho+ 

if ( PcheckDecay(genpart,-3212, 2112, 213)==1){
return 1450;}//B+ decays to anti-Sigma0 n0 rho+ 

if ( PcheckDecay(genpart,-2212, 2224, 313)==1){
return 1451;}//B+ decays to anti-p- Delta++ K*0 

if ( PcheckDecay(genpart,-2112, 2212, 313)==1){
return 1452;}//B+ decays to anti-n0 p+ K*0 

if ( PcheckDecay(genpart,-1114, 2112, 313)==1){
return 1453;}//B+ decays to anti-Delta+ n0 K*0 

if ( PcheckDecay(genpart,-2224, 2224, 323)==1){
return 1454;}//B+ decays to anti-Delta-- Delta++ K*+ 

if ( PcheckDecay(genpart,-2212, 2212, 323)==1){
return 1455;}//B+ decays to anti-p- p+ K*+ 

if ( PcheckDecay(genpart,-2112, 2112, 323)==1){
return 1456;}//B+ decays to anti-n0 n0 K*+ 

if ( PcheckDecay(genpart,-3322, 3222, 113)==1){
return 1457;}//B+ decays to anti-Xi0 Sigma+ rho0 

if ( PcheckDecay(genpart,-3312, 3122, 113)==1){
return 1458;}//B+ decays to anti-Xi+ Lambda0 rho0 

if ( PcheckDecay(genpart,-3312, 3212, 113)==1){
return 1459;}//B+ decays to anti-Xi+ Sigma0 rho0 

if ( PcheckDecay(genpart,-3312, 3222, -213)==1){
return 1460;}//B+ decays to anti-Xi+ Sigma+ rho- 

if ( PcheckDecay(genpart,-3322, 3122, 213)==1){
return 1461;}//B+ decays to anti-Xi0 Lambda0 rho+ 

if ( PcheckDecay(genpart,-3322, 3212, 213)==1){
return 1462;}//B+ decays to anti-Xi0 Sigma0 rho+ 

if ( PcheckDecay(genpart,-3312, 3112, 213)==1){
return 1463;}//B+ decays to anti-Xi+ Sigma- rho+ 

if ( PcheckDecay(genpart,-3322, 2212, -313)==1){
return 1464;}//B+ decays to anti-Xi0 p+ anti-K*0 

if ( PcheckDecay(genpart,-3322, 2224, -323)==1){
return 1465;}//B+ decays to anti-Xi0 Delta++ K*- 

if ( PcheckDecay(genpart,-3312, 2112, -313)==1){
return 1466;}//B+ decays to anti-Xi+ n0 anti-K*0 

if ( PcheckDecay(genpart,-3312, 2212, -323)==1){
return 1467;}//B+ decays to anti-Xi+ p+ K*- 

if ( PcheckDecay(genpart,-3122, 3222, 313)==1){
return 1468;}//B+ decays to anti-Lambda0 Sigma+ K*0 

if ( PcheckDecay(genpart,-3212, 3222, 313)==1){
return 1469;}//B+ decays to anti-Sigma0 Sigma+ K*0 

if ( PcheckDecay(genpart,-3222, 3222, 323)==1){
return 1470;}//B+ decays to anti-Sigma- Sigma+ K*+ 

if ( PcheckDecay(genpart,-3112, 3122, 313)==1){
return 1471;}//B+ decays to anti-Sigma+ Lambda0 K*0 

if ( PcheckDecay(genpart,-3112, 3212, 313)==1){
return 1472;}//B+ decays to anti-Sigma+ Sigma0 K*0 

if ( PcheckDecay(genpart,-3122, 3122, 323)==1){
return 1473;}//B+ decays to anti-Lambda0 Lambda0 K*+ 

if ( PcheckDecay(genpart,-3122, 3212, 323)==1){
return 1474;}//B+ decays to anti-Lambda0 Sigma0 K*+ 

if ( PcheckDecay(genpart,-3212, 3122, 323)==1){
return 1475;}//B+ decays to anti-Sigma0 Lambda0 K*+ 

if ( PcheckDecay(genpart,-3212, 3212, 323)==1){
return 1476;}//B+ decays to anti-Sigma0 Sigma0 K*+ 

if ( PcheckDecay(genpart,-3334, 3322, 113)==1){
return 1477;}//B+ decays to anti-Omega+ Xi0 rho0 

if ( PcheckDecay(genpart,-3334, 3312, 213)==1){
return 1478;}//B+ decays to anti-Omega+ Xi- rho+ 

if ( PcheckDecay(genpart,-3334, 3122, -313)==1){
return 1479;}//B+ decays to anti-Omega+ Lambda0 anti-K*0 

if ( PcheckDecay(genpart,-3334, 3212, -313)==1){
return 1480;}//B+ decays to anti-Omega+ Sigma0 anti-K*0 

if ( PcheckDecay(genpart,-3334, 3222, -323)==1){
return 1481;}//B+ decays to anti-Omega+ Sigma+ K*- 

if ( PcheckDecay(genpart,-3312, 3322, 313)==1){
return 1482;}//B+ decays to anti-Xi+ Xi0 K*0 

if ( PcheckDecay(genpart,-3322, 3322, 323)==1){
return 1483;}//B+ decays to anti-Xi0 Xi0 K*+ 

if ( PcheckDecay(genpart,-2212, 2224, 111)==1){
return 1484;}//B+ decays to anti-p- Delta++ pi0 

if ( PcheckDecay(genpart,-2112, 2212, 111)==1){
return 1485;}//B+ decays to anti-n0 p+ pi0 

if ( PcheckDecay(genpart,-1114, 2112, 111)==1){
return 1486;}//B+ decays to anti-Delta+ n0 pi0 

if ( PcheckDecay(genpart,-2212, 2212, 211)==1){
return 1487;}//B+ decays to anti-p- p+ pi+ 

if ( PcheckDecay(genpart,-2112, 2224, -211)==1){
return 1488;}//B+ decays to anti-n0 Delta++ pi- 

if ( PcheckDecay(genpart,-2224, 2224, 211)==1){
return 1489;}//B+ decays to anti-Delta-- Delta++ pi+ 

if ( PcheckDecay(genpart,-2112, 2112, 211)==1){
return 1490;}//B+ decays to anti-n0 n0 pi+ 

if ( PcheckDecay(genpart,-1114, 2212, -211)==1){
return 1491;}//B+ decays to anti-Delta+ p+ pi- 

if ( PcheckDecay(genpart,-1114, 1114, 211)==1){
return 1492;}//B+ decays to anti-Delta+ Delta- pi+ 

if ( PcheckDecay(genpart,-3122, 3222, 111)==1){
return 1493;}//B+ decays to anti-Lambda0 Sigma+ pi0 

if ( PcheckDecay(genpart,-3212, 3222, 111)==1){
return 1494;}//B+ decays to anti-Sigma0 Sigma+ pi0 

if ( PcheckDecay(genpart,-3112, 3122, 111)==1){
return 1495;}//B+ decays to anti-Sigma+ Lambda0 pi0 

if ( PcheckDecay(genpart,-3112, 3212, 111)==1){
return 1496;}//B+ decays to anti-Sigma+ Sigma0 pi0 

if ( PcheckDecay(genpart,-3122, 3122, 211)==1){
return 1497;}//B+ decays to anti-Lambda0 Lambda0 pi+ 

if ( PcheckDecay(genpart,-3122, 3212, 211)==1){
return 1498;}//B+ decays to anti-Lambda0 Sigma0 pi+ 

if ( PcheckDecay(genpart,-3212, 3122, 211)==1){
return 1499;}//B+ decays to anti-Sigma0 Lambda0 pi+ 

if ( PcheckDecay(genpart,-3212, 3212, 211)==1){
return 1500;}//B+ decays to anti-Sigma0 Sigma0 pi+ 

if ( PcheckDecay(genpart,-3112, 3222, -211)==1){
return 1501;}//B+ decays to anti-Sigma+ Sigma+ pi- 

if ( PcheckDecay(genpart,-3222, 3222, 211)==1){
return 1502;}//B+ decays to anti-Sigma- Sigma+ pi+ 

if ( PcheckDecay(genpart,-3112, 3112, 211)==1){
return 1503;}//B+ decays to anti-Sigma+ Sigma- pi+ 

if ( PcheckDecay(genpart,-3122, 2212, -311)==1){
return 1504;}//B+ decays to anti-Lambda0 p+ anti-K0 

if ( PcheckDecay(genpart,-3212, 2212, -311)==1){
return 1505;}//B+ decays to anti-Sigma0 p+ anti-K0 

if ( PcheckDecay(genpart,-3122, 2224, -321)==1){
return 1506;}//B+ decays to anti-Lambda0 Delta++ K- 

if ( PcheckDecay(genpart,-3212, 2224, -321)==1){
return 1507;}//B+ decays to anti-Sigma0 Delta++ K- 

if ( PcheckDecay(genpart,-3112, 2112, -311)==1){
return 1508;}//B+ decays to anti-Sigma+ n0 anti-K0 

if ( PcheckDecay(genpart,-3112, 2214, -321)==1){
return 1509;}//B+ decays to anti-Sigma+ Delta+ K- 

if ( PcheckDecay(genpart,-2112, 3222, -311)==1){
return 1510;}//B+ decays to anti-n0 Sigma+ anti-K0 

if ( PcheckDecay(genpart,-2212, 3222, 321)==1){
return 1511;}//B+ decays to anti-p- Sigma+ K+ 

if ( PcheckDecay(genpart,-1114, 3122, 311)==1){
return 1512;}//B+ decays to anti-Delta+ Lambda0 K0 

if ( PcheckDecay(genpart,-1114, 3212, 311)==1){
return 1513;}//B+ decays to anti-Delta+ Sigma0 K0 

if ( PcheckDecay(genpart,-2112, 3122, 321)==1){
return 1514;}//B+ decays to anti-n0 Lambda0 K+ 

if ( PcheckDecay(genpart,-2112, 3212, 321)==1){
return 1515;}//B+ decays to anti-n0 Sigma0 K+ 

if ( PcheckDecay(genpart,-3312, 3322, 111)==1){
return 1516;}//B+ decays to anti-Xi+ Xi0 pi0 

if ( PcheckDecay(genpart,-3312, 3312, 211)==1){
return 1517;}//B+ decays to anti-Xi+ Xi- pi+ 

if ( PcheckDecay(genpart,-3322, 3322, 211)==1){
return 1518;}//B+ decays to anti-Xi0 Xi0 pi+ 

if ( PcheckDecay(genpart,-3312, 3122, -311)==1){
return 1519;}//B+ decays to anti-Xi+ Lambda0 anti-K0 

if ( PcheckDecay(genpart,-3312, 3212, -311)==1){
return 1520;}//B+ decays to anti-Xi+ Sigma0 anti-K0 

if ( PcheckDecay(genpart,-3312, 3222, -321)==1){
return 1521;}//B+ decays to anti-Xi+ Sigma+ K- 

if ( PcheckDecay(genpart,-3112, 3322, 311)==1){
return 1522;}//B+ decays to anti-Sigma+ Xi0 K0 

if ( PcheckDecay(genpart,-3122, 3322, 321)==1){
return 1523;}//B+ decays to anti-Lambda0 Xi0 K+ 

if ( PcheckDecay(genpart,-3212, 3322, 321)==1){
return 1524;}//B+ decays to anti-Sigma0 Xi0 K+ 

if ( PcheckDecay(genpart,-2212, 2224, 113)==1){
return 1525;}//B+ decays to anti-p- Delta++ rho0 

if ( PcheckDecay(genpart,-2112, 2212, 113)==1){
return 1526;}//B+ decays to anti-n0 p+ rho0 

if ( PcheckDecay(genpart,-1114, 2112, 113)==1){
return 1527;}//B+ decays to anti-Delta+ n0 rho0 

if ( PcheckDecay(genpart,-2212, 2212, 213)==1){
return 1528;}//B+ decays to anti-p- p+ rho+ 

if ( PcheckDecay(genpart,-2112, 2224, -213)==1){
return 1529;}//B+ decays to anti-n0 Delta++ rho- 

if ( PcheckDecay(genpart,-2224, 2224, 213)==1){
return 1530;}//B+ decays to anti-Delta-- Delta++ rho+ 

if ( PcheckDecay(genpart,-2112, 2112, 213)==1){
return 1531;}//B+ decays to anti-n0 n0 rho+ 

if ( PcheckDecay(genpart,-1114, 2212, -213)==1){
return 1532;}//B+ decays to anti-Delta+ p+ rho- 

if ( PcheckDecay(genpart,-1114, 1114, 213)==1){
return 1533;}//B+ decays to anti-Delta+ Delta- rho+ 

if ( PcheckDecay(genpart,-3122, 3222, 113)==1){
return 1534;}//B+ decays to anti-Lambda0 Sigma+ rho0 

if ( PcheckDecay(genpart,-3212, 3222, 113)==1){
return 1535;}//B+ decays to anti-Sigma0 Sigma+ rho0 

if ( PcheckDecay(genpart,-3112, 3122, 113)==1){
return 1536;}//B+ decays to anti-Sigma+ Lambda0 rho0 

if ( PcheckDecay(genpart,-3112, 3212, 113)==1){
return 1537;}//B+ decays to anti-Sigma+ Sigma0 rho0 

if ( PcheckDecay(genpart,-3122, 3122, 213)==1){
return 1538;}//B+ decays to anti-Lambda0 Lambda0 rho+ 

if ( PcheckDecay(genpart,-3122, 3212, 213)==1){
return 1539;}//B+ decays to anti-Lambda0 Sigma0 rho+ 

if ( PcheckDecay(genpart,-3212, 3122, 213)==1){
return 1540;}//B+ decays to anti-Sigma0 Lambda0 rho+ 

if ( PcheckDecay(genpart,-3212, 3212, 213)==1){
return 1541;}//B+ decays to anti-Sigma0 Sigma0 rho+ 

if ( PcheckDecay(genpart,-3112, 3222, -213)==1){
return 1542;}//B+ decays to anti-Sigma+ Sigma+ rho- 

if ( PcheckDecay(genpart,-3222, 3222, 213)==1){
return 1543;}//B+ decays to anti-Sigma- Sigma+ rho+ 

if ( PcheckDecay(genpart,-3112, 3112, 213)==1){
return 1544;}//B+ decays to anti-Sigma+ Sigma- rho+ 

if ( PcheckDecay(genpart,-3122, 2212, -313)==1){
return 1545;}//B+ decays to anti-Lambda0 p+ anti-K*0 

if ( PcheckDecay(genpart,-3212, 2212, -313)==1){
return 1546;}//B+ decays to anti-Sigma0 p+ anti-K*0 

if ( PcheckDecay(genpart,-3122, 2224, -323)==1){
return 1547;}//B+ decays to anti-Lambda0 Delta++ K*- 

if ( PcheckDecay(genpart,-3212, 2224, -323)==1){
return 1548;}//B+ decays to anti-Sigma0 Delta++ K*- 

if ( PcheckDecay(genpart,-3112, 2112, -313)==1){
return 1549;}//B+ decays to anti-Sigma+ n0 anti-K*0 

if ( PcheckDecay(genpart,-3112, 2214, -323)==1){
return 1550;}//B+ decays to anti-Sigma+ Delta+ K*- 

if ( PcheckDecay(genpart,-2112, 3222, 313)==1){
return 1551;}//B+ decays to anti-n0 Sigma+ K*0 

if ( PcheckDecay(genpart,-2212, 3222, 323)==1){
return 1552;}//B+ decays to anti-p- Sigma+ K*+ 

if ( PcheckDecay(genpart,-1114, 3122, 313)==1){
return 1553;}//B+ decays to anti-Delta+ Lambda0 K*0 

if ( PcheckDecay(genpart,-1114, 3212, 313)==1){
return 1554;}//B+ decays to anti-Delta+ Sigma0 K*0 

if ( PcheckDecay(genpart,-2112, 3122, 323)==1){
return 1555;}//B+ decays to anti-n0 Lambda0 K*+ 

if ( PcheckDecay(genpart,-2112, 3212, 323)==1){
return 1556;}//B+ decays to anti-n0 Sigma0 K*+ 

if ( PcheckDecay(genpart,-3312, 3322, 113)==1){
return 1557;}//B+ decays to anti-Xi+ Xi0 rho0 

if ( PcheckDecay(genpart,-3312, 3312, 213)==1){
return 1558;}//B+ decays to anti-Xi+ Xi- rho+ 

if ( PcheckDecay(genpart,-3322, 3322, 213)==1){
return 1559;}//B+ decays to anti-Xi0 Xi0 rho+ 

if ( PcheckDecay(genpart,-3312, 3122, -313)==1){
return 1560;}//B+ decays to anti-Xi+ Lambda0 anti-K*0 

if ( PcheckDecay(genpart,-3312, 3212, -313)==1){
return 1561;}//B+ decays to anti-Xi+ Sigma0 anti-K*0 

if ( PcheckDecay(genpart,-3312, 3222, -323)==1){
return 1562;}//B+ decays to anti-Xi+ Sigma+ K*- 

if ( PcheckDecay(genpart,-3112, 3322, 313)==1){
return 1563;}//B+ decays to anti-Sigma+ Xi0 K*0 

if ( PcheckDecay(genpart,-3122, 3322, 323)==1){
return 1564;}//B+ decays to anti-Lambda0 Xi0 K*+ 

if ( PcheckDecay(genpart,-3212, 3322, 323)==1){
return 1565;}//B+ decays to anti-Sigma0 Xi0 K*+ 

if ( PcheckDecay(genpart,-3222, 2224)==1){
return 1566;}//B+ decays to anti-Sigma- Delta++ 

if ( PcheckDecay(genpart,-3122, 2212)==1){
return 1567;}//B+ decays to anti-Lambda0 p+ 

if ( PcheckDecay(genpart,-3122, 2214)==1){
return 1568;}//B+ decays to anti-Lambda0 Delta+ 

if ( PcheckDecay(genpart,-3212, 2212)==1){
return 1569;}//B+ decays to anti-Sigma0 p+ 

if ( PcheckDecay(genpart,-3214, 2212)==1){
return 1570;}//B+ decays to anti-Sigma*0 p+ 

if ( PcheckDecay(genpart,-3112, 2112)==1){
return 1571;}//B+ decays to anti-Sigma+ n0 

if ( PcheckDecay(genpart,-3322, 3222)==1){
return 1572;}//B+ decays to anti-Xi0 Sigma+ 

if ( PcheckDecay(genpart,-3312, 3122)==1){
return 1573;}//B+ decays to anti-Xi+ Lambda0 

if ( PcheckDecay(genpart,-3312, 3212)==1){
return 1574;}//B+ decays to anti-Xi+ Sigma0 

if ( PcheckDecay(genpart,-3334, 3322)==1){
return 1575;}//B+ decays to anti-Omega+ Xi0 

if ( PcheckDecay(genpart,-2212, 2224)==1){
return 1576;}//B+ decays to anti-p- Delta++ 

if ( PcheckDecay(genpart,-2112, 2212)==1){
return 1577;}//B+ decays to anti-n0 p+ 

if ( PcheckDecay(genpart,-1114, 2112)==1){
return 1578;}//B+ decays to anti-Delta+ n0 

if ( PcheckDecay(genpart,-3122, 3222)==1){
return 1579;}//B+ decays to anti-Lambda0 Sigma+ 

if ( PcheckDecay(genpart,-3212, 3222)==1){
return 1580;}//B+ decays to anti-Sigma0 Sigma+ 

if ( PcheckDecay(genpart,-3112, 3122)==1){
return 1581;}//B+ decays to anti-Sigma+ Lambda0 

if ( PcheckDecay(genpart,-3112, 3212)==1){
return 1582;}//B+ decays to anti-Sigma+ Sigma0 

if ( PcheckDecay(genpart,-3312, 3322)==1){
return 1583;}//B+ decays to anti-Xi+ Xi0 

if ( PcheckDecay(genpart,3124, 2212)==1){
return 1584;}//B+ decays to Lambda(1520)0 p+ 

if ( PcheckDecay(genpart,431, 111)==1){
return 1585;}//B+ decays to D_s+ pi0 

if ( PcheckDecay(genpart,431, 221)==1){
return 1586;}//B+ decays to D_s+ eta 

if ( PcheckDecay(genpart,431, 331)==1){
return 1587;}//B+ decays to D_s+ eta' 

if ( PcheckDecay(genpart,113, 431)==1){
return 1588;}//B+ decays to rho0 D_s+ 

if ( PcheckDecay(genpart,431, 211, -211)==1){
return 1589;}//B+ decays to D_s+ pi+ pi- 

if ( PcheckDecay(genpart,431, 111, 111)==1){
return 1590;}//B+ decays to D_s+ pi0 pi0 

if ( PcheckDecay(genpart,20113, 431)==1){
return 1591;}//B+ decays to a_10 D_s+ 

if ( PcheckDecay(genpart,223, 431)==1){
return 1592;}//B+ decays to omega D_s+ 

if ( PcheckDecay(genpart,433, 111)==1){
return 1593;}//B+ decays to D_s*+ pi0 

if ( PcheckDecay(genpart,433, 221)==1){
return 1594;}//B+ decays to D_s*+ eta 

if ( PcheckDecay(genpart,433, 331)==1){
return 1595;}//B+ decays to D_s*+ eta' 

if ( PcheckDecay(genpart,433, 113)==1){
return 1596;}//B+ decays to D_s*+ rho0 

if ( PcheckDecay(genpart,433, 211, -211)==1){
return 1597;}//B+ decays to D_s*+ pi+ pi- 

if ( PcheckDecay(genpart,433, 111, 111)==1){
return 1598;}//B+ decays to D_s*+ pi0 pi0 

if ( PcheckDecay(genpart,433, 20113)==1){
return 1599;}//B+ decays to D_s*+ a_10 

if ( PcheckDecay(genpart,433, 223)==1){
return 1600;}//B+ decays to D_s*+ omega 

if ( PcheckDecay(genpart,411, 111)==1){
return 1601;}//B+ decays to D+ pi0 

if ( PcheckDecay(genpart,411, 221)==1){
return 1602;}//B+ decays to D+ eta 

if ( PcheckDecay(genpart,411, 331)==1){
return 1603;}//B+ decays to D+ eta' 

if ( PcheckDecay(genpart,113, 411)==1){
return 1604;}//B+ decays to rho0 D+ 

if ( PcheckDecay(genpart,411, 211, -211)==1){
return 1605;}//B+ decays to D+ pi+ pi- 

if ( PcheckDecay(genpart,411, 111, 111)==1){
return 1606;}//B+ decays to D+ pi0 pi0 

if ( PcheckDecay(genpart,20113, 411)==1){
return 1607;}//B+ decays to a_10 D+ 

if ( PcheckDecay(genpart,223, 411)==1){
return 1608;}//B+ decays to omega D+ 

if ( PcheckDecay(genpart,413, 111)==1){
return 1609;}//B+ decays to D*+ pi0 

if ( PcheckDecay(genpart,413, 221)==1){
return 1610;}//B+ decays to D*+ eta 

if ( PcheckDecay(genpart,413, 331)==1){
return 1611;}//B+ decays to D*+ eta' 

if ( PcheckDecay(genpart,413, 113)==1){
return 1612;}//B+ decays to D*+ rho0 

if ( PcheckDecay(genpart,413, 211, -211)==1){
return 1613;}//B+ decays to D*+ pi+ pi- 

if ( PcheckDecay(genpart,413, 111, 111)==1){
return 1614;}//B+ decays to D*+ pi0 pi0 

if ( PcheckDecay(genpart,413, 20113)==1){
return 1615;}//B+ decays to D*+ a_10 

if ( PcheckDecay(genpart,413, 223)==1){
return 1616;}//B+ decays to D*+ omega 

if ( PcheckDecay(genpart,411, 311)==1){
return 1617;}//B+ decays to D+ K0 

if ( PcheckDecay(genpart,313, 411)==1){
return 1618;}//B+ decays to K*0 D+ 

if ( PcheckDecay(genpart,411, 321, -211)==1){
return 1619;}//B+ decays to D+ K+ pi- 

if ( PcheckDecay(genpart,411, 311, 111)==1){
return 1620;}//B+ decays to D+ K0 pi0 

if ( PcheckDecay(genpart,413, 311)==1){
return 1621;}//B+ decays to D*+ K0 

if ( PcheckDecay(genpart,413, 313)==1){
return 1622;}//B+ decays to D*+ K*0 

if ( PcheckDecay(genpart,413, 321, -211)==1){
return 1623;}//B+ decays to D*+ K+ pi- 

if ( PcheckDecay(genpart,413, 311, 111)==1){
return 1624;}//B+ decays to D*+ K0 pi0 

if ( PcheckDecay(genpart,333, 431)==1){
return 1625;}//B+ decays to phi D_s+ 

if ( PcheckDecay(genpart,333, 433)==1){
return 1626;}//B+ decays to phi D_s*+ 

if ( PcheckDecay(genpart,431, -311)==1){
return 1627;}//B+ decays to D_s+ anti-K0 

if ( PcheckDecay(genpart,431, -321, 211)==1){
return 1628;}//B+ decays to D_s+ K- pi+ 

if ( PcheckDecay(genpart,431, -311, 111)==1){
return 1629;}//B+ decays to D_s+ anti-K0 pi0 

if ( PcheckDecay(genpart,-313, 431)==1){
return 1630;}//B+ decays to anti-K*0 D_s+ 

if ( PcheckDecay(genpart,433, -311)==1){
return 1631;}//B+ decays to D_s*+ anti-K0 

if ( PcheckDecay(genpart,433, -321, 211)==1){
return 1632;}//B+ decays to D_s*+ K- pi+ 

if ( PcheckDecay(genpart,433, -311, 111)==1){
return 1633;}//B+ decays to D_s*+ anti-K0 pi0 

if ( PcheckDecay(genpart,433, -313)==1){
return 1634;}//B+ decays to D_s*+ anti-K*0 

if ( PcheckDecay(genpart,433, 22)==1){
return 1635;}//B+ decays to D_s*+ gamma 

if ( PcheckDecay(genpart,413, 22)==1){
return 1636;}//B+ decays to D*+ gamma 

if ( PcheckDecay(genpart,-13, 14)==1){
return 1637;}//B+ decays to mu+ nu_mu 

if ( PcheckDecay(genpart,-15, 16)==1){
return 1638;}//B+ decays to tau+ nu_tau 

if ( PcheckDecay(genpart,-11, 12, 22)==1){
return 1639;}//B+ decays to e+ nu_e gamma 

if ( PcheckDecay(genpart,-13, 14, 22)==1){
  return 1640;}//B+ decays to mu+ nu_mu gamma 
    
    return -300;
  }
  else if (AmId(genpart)==-521){


if ( PcheckDecay(genpart,-323, 22)==1){
return 1641;}//B- decays to K*- gamma 

if ( PcheckDecay(genpart,-10323, 22)==1){
return 1642;}//B- decays to K_1- gamma 

if ( PcheckDecay(genpart,-325, 22)==1){
return 1643;}//B- decays to K_2*- gamma 

if ( PcheckDecay(genpart,-30353, 22)==1){
return 1644;}//B- decays to anti-Xsu gamma 

if ( PcheckDecay(genpart,-213, 22)==1){
return 1645;}//B- decays to rho- gamma 

if ( PcheckDecay(genpart,-30653, 22)==1){
return 1646;}//B- decays to anti-Xdu- gamma 

if ( PcheckDecay(genpart,-321, -11, 11)==1){
return 1647;}//B- decays to K- e+ e- 

if ( PcheckDecay(genpart,-323, -11, 11)==1){
return 1648;}//B- decays to K*- e+ e- 

if ( PcheckDecay(genpart,-30353, -11, 11)==1){
return 1649;}//B- decays to anti-Xsu e+ e- 

if ( PcheckDecay(genpart,-321, -13, 13)==1){
return 1650;}//B- decays to K- mu+ mu- 

if ( PcheckDecay(genpart,-323, -13, 13)==1){
return 1651;}//B- decays to K*- mu+ mu- 

if ( PcheckDecay(genpart,-30353, -13, 13)==1){
return 1652;}//B- decays to anti-Xsu mu+ mu- 

if ( PcheckDecay(genpart,-321, -15, 15)==1){
return 1653;}//B- decays to K- tau+ tau- 

if ( PcheckDecay(genpart,-323, -15, 15)==1){
return 1654;}//B- decays to K*- tau+ tau- 

if ( PcheckDecay(genpart,-30353, -15, 15)==1){
return 1655;}//B- decays to anti-Xsu tau+ tau- 

if ( PcheckDecay(genpart,-321, 12, -12)==1){
return 1656;}//B- decays to K- nu_e anti-nu_e 

if ( PcheckDecay(genpart,-323, 12, -12)==1){
return 1657;}//B- decays to K*- nu_e anti-nu_e 

if ( PcheckDecay(genpart,-30353, 12, -12)==1){
return 1658;}//B- decays to anti-Xsu nu_e anti-nu_e 

if ( PcheckDecay(genpart,-321, 14, -14)==1){
return 1659;}//B- decays to K- nu_mu anti-nu_mu 

if ( PcheckDecay(genpart,-323, 14, -14)==1){
return 1660;}//B- decays to K*- nu_mu anti-nu_mu 

if ( PcheckDecay(genpart,-30353, 14, -14)==1){
return 1661;}//B- decays to anti-Xsu nu_mu anti-nu_mu 

if ( PcheckDecay(genpart,-321, 16, -16)==1){
return 1662;}//B- decays to K- nu_tau anti-nu_tau 

if ( PcheckDecay(genpart,-323, 16, -16)==1){
return 1663;}//B- decays to K*- nu_tau anti-nu_tau 

if ( PcheckDecay(genpart,-30353, 16, -16)==1){
return 1664;}//B- decays to anti-Xsu nu_tau anti-nu_tau 

if ( PcheckDecay(genpart,-321, 311)==1){
return 1665;}//B- decays to K- K0 

if ( PcheckDecay(genpart,-321, 111)==1){
return 1666;}//B- decays to K- pi0 

if ( PcheckDecay(genpart,-311, -211)==1){
return 1667;}//B- decays to anti-K0 pi- 

if ( PcheckDecay(genpart,-211, 111)==1){
return 1668;}//B- decays to pi- pi0 

if ( PcheckDecay(genpart,-100323, 111)==1){
return 1669;}//B- decays to K'*- pi0 

if ( PcheckDecay(genpart,-100313, -211)==1){
return 1670;}//B- decays to anti-K'*0 pi- 

if ( PcheckDecay(genpart,-10311, -211)==1){
return 1671;}//B- decays to anti-K_0*0 pi- 

if ( PcheckDecay(genpart,10311, -321)==1){
return 1672;}//B- decays to K_0*0 K- 

if ( PcheckDecay(genpart,10331, -321)==1){
return 1673;}//B- decays to f'_0 K- 

if ( PcheckDecay(genpart,10221, -321)==1){
return 1674;}//B- decays to f_0 K- 

if ( PcheckDecay(genpart,9030221, -321)==1){
return 1675;}//B- decays to f_0(1500) K- 

if ( PcheckDecay(genpart,10331, -211)==1){
return 1676;}//B- decays to f'_0 pi- 

if ( PcheckDecay(genpart,10221, -211)==1){
return 1677;}//B- decays to f_0 pi- 

if ( PcheckDecay(genpart,9000221, -211)==1){
return 1678;}//B- decays to f_0(600) pi- 

if ( PcheckDecay(genpart,-315, -211)==1){
return 1679;}//B- decays to anti-K_2*0 pi- 

if ( PcheckDecay(genpart,335, -321)==1){
return 1680;}//B- decays to f'_2 K- 

if ( PcheckDecay(genpart,225, -321)==1){
return 1681;}//B- decays to f_2 K- 

if ( PcheckDecay(genpart,10111, -321)==1){
return 1682;}//B- decays to a_00 K- 

if ( PcheckDecay(genpart,-10211, -311)==1){
return 1683;}//B- decays to a_0- anti-K0 

if ( PcheckDecay(genpart,10111, -211)==1){
return 1684;}//B- decays to a_00 pi- 

if ( PcheckDecay(genpart,-10211, 111)==1){
return 1685;}//B- decays to a_0- pi0 

if ( PcheckDecay(genpart,-30313, -211)==1){
return 1686;}//B- decays to anti-K''*0 pi- 

if ( PcheckDecay(genpart,-323, 311)==1){
return 1687;}//B- decays to K*- K0 

if ( PcheckDecay(genpart,-323, 111)==1){
return 1688;}//B- decays to K*- pi0 

if ( PcheckDecay(genpart,-323, 10221)==1){
return 1689;}//B- decays to K*- f_0 

if ( PcheckDecay(genpart,-313, -211)==1){
return 1690;}//B- decays to anti-K*0 pi- 

if ( PcheckDecay(genpart,313, -321)==1){
return 1691;}//B- decays to K*0 K- 

if ( PcheckDecay(genpart,113, -321)==1){
return 1692;}//B- decays to rho0 K- 

if ( PcheckDecay(genpart,100113, -321)==1){
return 1693;}//B- decays to rho(2S)0 K- 

if ( PcheckDecay(genpart,30113, -321)==1){
return 1694;}//B- decays to rho(3S)0 K- 

if ( PcheckDecay(genpart,-213, -311)==1){
return 1695;}//B- decays to rho- anti-K0 

if ( PcheckDecay(genpart,-100213, -311)==1){
return 1696;}//B- decays to rho(2S)- anti-K0 

if ( PcheckDecay(genpart,-30213, -311)==1){
return 1697;}//B- decays to rho(3S)- anti-K0 

if ( PcheckDecay(genpart,-213, 111)==1){
return 1698;}//B- decays to rho- pi0 

if ( PcheckDecay(genpart,-100213, 111)==1){
return 1699;}//B- decays to rho(2S)- pi0 

if ( PcheckDecay(genpart,-30213, 111)==1){
return 1700;}//B- decays to rho(3S)- pi0 

if ( PcheckDecay(genpart,113, -211)==1){
return 1701;}//B- decays to rho0 pi- 

if ( PcheckDecay(genpart,100113, -211)==1){
return 1702;}//B- decays to rho(2S)0 pi- 

if ( PcheckDecay(genpart,30113, -211)==1){
return 1703;}//B- decays to rho(3S)0 pi- 

if ( PcheckDecay(genpart,20113, -321)==1){
return 1704;}//B- decays to a_10 K- 

if ( PcheckDecay(genpart,-20213, 311)==1){
return 1705;}//B- decays to a_1- K0 

if ( PcheckDecay(genpart,20113, -211)==1){
return 1706;}//B- decays to a_10 pi- 

if ( PcheckDecay(genpart,-20213, 111)==1){
return 1707;}//B- decays to a_1- pi0 

if ( PcheckDecay(genpart,-213, 10221)==1){
return 1708;}//B- decays to rho- f_0 

if ( PcheckDecay(genpart,10113, -321)==1){
return 1709;}//B- decays to b_10 K- 

if ( PcheckDecay(genpart,10113, -211)==1){
return 1710;}//B- decays to b_10 pi- 

if ( PcheckDecay(genpart,225, -211)==1){
return 1711;}//B- decays to f_2 pi- 

if ( PcheckDecay(genpart,115, -321)==1){
return 1712;}//B- decays to a_20 K- 

if ( PcheckDecay(genpart,-321, -321, 321)==1){
return 1713;}//B- decays to K- K- K+ 

if ( PcheckDecay(genpart,-321, 321, -211)==1){
return 1714;}//B- decays to K- K+ pi- 

if ( PcheckDecay(genpart,-321, 310, 310)==1){
return 1715;}//B- decays to K- K_S0 K_S0 

if ( PcheckDecay(genpart,-321, 311, 111)==1){
return 1716;}//B- decays to K- K0 pi0 

if ( PcheckDecay(genpart,-321, -211, 211)==1){
return 1717;}//B- decays to K- pi- pi+ 

if ( PcheckDecay(genpart,-311, -211, 111)==1){
return 1718;}//B- decays to anti-K0 pi- pi0 

if ( PcheckDecay(genpart,-211, -211, 211)==1){
return 1719;}//B- decays to pi- pi- pi+ 

if ( PcheckDecay(genpart,-321, -321, 211)==1){
return 1720;}//B- decays to K- K- pi+ 

if ( PcheckDecay(genpart,321, -211, -211)==1){
return 1721;}//B- decays to K+ pi- pi- 

if ( PcheckDecay(genpart,310, 310, -211)==1){
return 1722;}//B- decays to K_S0 K_S0 pi- 

if ( PcheckDecay(genpart,130, 130, -211)==1){
return 1723;}//B- decays to K_L0 K_L0 pi- 

if ( PcheckDecay(genpart,310, 130, -211)==1){
return 1724;}//B- decays to K_S0 K_L0 pi- 

if ( PcheckDecay(genpart,335, -323)==1){
return 1725;}//B- decays to f'_2 K*- 

if ( PcheckDecay(genpart,225, -323)==1){
return 1726;}//B- decays to f_2 K*- 

if ( PcheckDecay(genpart,-213, -10311)==1){
return 1727;}//B- decays to rho- anti-K_0*0 

if ( PcheckDecay(genpart,113, -10321)==1){
return 1728;}//B- decays to rho0 K_0*- 

if ( PcheckDecay(genpart,-323, 313)==1){
return 1729;}//B- decays to K*- K*0 

if ( PcheckDecay(genpart,-323, 113)==1){
return 1730;}//B- decays to K*- rho0 

if ( PcheckDecay(genpart,-313, -213)==1){
return 1731;}//B- decays to anti-K*0 rho- 

if ( PcheckDecay(genpart,-213, 113)==1){
return 1732;}//B- decays to rho- rho0 

if ( PcheckDecay(genpart,-323, -211, 211)==1){
return 1733;}//B- decays to K*- pi- pi+ 

if ( PcheckDecay(genpart,-323, 111, 111)==1){
return 1734;}//B- decays to K*- pi0 pi0 

if ( PcheckDecay(genpart,-313, -211, 111)==1){
return 1735;}//B- decays to anti-K*0 pi- pi0 

if ( PcheckDecay(genpart,-213, -321, 211)==1){
return 1736;}//B- decays to rho- K- pi+ 

if ( PcheckDecay(genpart,-213, -311, 111)==1){
return 1737;}//B- decays to rho- anti-K0 pi0 

if ( PcheckDecay(genpart,113, -321, 111)==1){
return 1738;}//B- decays to rho0 K- pi0 

if ( PcheckDecay(genpart,113, -311, -211)==1){
return 1739;}//B- decays to rho0 anti-K0 pi- 

if ( PcheckDecay(genpart,321, -321, -323)==1){
return 1740;}//B- decays to K+ K- K*- 

if ( PcheckDecay(genpart,311, -311, -323)==1){
return 1741;}//B- decays to K0 anti-K0 K*- 

if ( PcheckDecay(genpart,311, -313, -321)==1){
return 1742;}//B- decays to K0 anti-K*0 K- 

if ( PcheckDecay(genpart,313, -311, -321)==1){
return 1743;}//B- decays to K*0 anti-K0 K- 

if ( PcheckDecay(genpart,323, -321, -321)==1){
return 1744;}//B- decays to K*+ K- K- 

if ( PcheckDecay(genpart,-323, 321, -211)==1){
return 1745;}//B- decays to K*- K+ pi- 

if ( PcheckDecay(genpart,-323, -321, 211)==1){
return 1746;}//B- decays to K*- K- pi+ 

if ( PcheckDecay(genpart,323, -321, -211)==1){
return 1747;}//B- decays to K*+ K- pi- 

if ( PcheckDecay(genpart,221, -321)==1){
return 1748;}//B- decays to eta K- 

if ( PcheckDecay(genpart,331, -321)==1){
return 1749;}//B- decays to eta' K- 

if ( PcheckDecay(genpart,-10321, 221)==1){
return 1750;}//B- decays to K_0*- eta 

if ( PcheckDecay(genpart,100221, -321)==1){
return 1751;}//B- decays to eta(2S) K- 

if ( PcheckDecay(genpart,9020221, -321)==1){
return 1752;}//B- decays to eta(1405) K- 

if ( PcheckDecay(genpart,221, -211)==1){
return 1753;}//B- decays to eta pi- 

if ( PcheckDecay(genpart,331, -211)==1){
return 1754;}//B- decays to eta' pi- 

if ( PcheckDecay(genpart,100221, -211)==1){
return 1755;}//B- decays to eta(2S) pi- 

if ( PcheckDecay(genpart,9020221, -211)==1){
return 1756;}//B- decays to eta(1405) pi- 

if ( PcheckDecay(genpart,223, -321)==1){
return 1757;}//B- decays to omega K- 

if ( PcheckDecay(genpart,333, -321)==1){
return 1758;}//B- decays to phi K- 

if ( PcheckDecay(genpart,223, -211)==1){
return 1759;}//B- decays to omega pi- 

if ( PcheckDecay(genpart,333, -211)==1){
return 1760;}//B- decays to phi pi- 

if ( PcheckDecay(genpart,100333, -321)==1){
return 1761;}//B- decays to phi(1680) K- 

if ( PcheckDecay(genpart,-323, 221)==1){
return 1762;}//B- decays to K*- eta 

if ( PcheckDecay(genpart,-323, 331)==1){
return 1763;}//B- decays to K*- eta' 

if ( PcheckDecay(genpart,-323, 100221)==1){
return 1764;}//B- decays to K*- eta(2S) 

if ( PcheckDecay(genpart,-323, 9020221)==1){
return 1765;}//B- decays to K*- eta(1405) 

if ( PcheckDecay(genpart,-213, 221)==1){
return 1766;}//B- decays to rho- eta 

if ( PcheckDecay(genpart,-213, 331)==1){
return 1767;}//B- decays to rho- eta' 

if ( PcheckDecay(genpart,-213, 100221)==1){
return 1768;}//B- decays to rho- eta(2S) 

if ( PcheckDecay(genpart,-213, 9020221)==1){
return 1769;}//B- decays to rho- eta(1405) 

if ( PcheckDecay(genpart,223, -323)==1){
return 1770;}//B- decays to omega K*- 

if ( PcheckDecay(genpart,-323, 333)==1){
return 1771;}//B- decays to K*- phi 

if ( PcheckDecay(genpart,223, -213)==1){
return 1772;}//B- decays to omega rho- 

if ( PcheckDecay(genpart,333, -213)==1){
return 1773;}//B- decays to phi rho- 

if ( PcheckDecay(genpart,-325, 221)==1){
return 1774;}//B- decays to K_2*- eta 

if ( PcheckDecay(genpart,-30353, 331)==1){
return 1775;}//B- decays to anti-Xsu eta' 

if ( PcheckDecay(genpart,-30353, 221)==1){
return 1776;}//B- decays to anti-Xsu eta 

//Manual entry
 if ( PcheckDecay(genpart,221, -311,-211)==1){return 1777;}//B- decays to  eta anti-K0 pi-

 if ( PcheckDecay(genpart,221,-321,111)==1){return 1778;}//B- decays to eta K-  pi0

if ( PcheckDecay(genpart,221, 111, -211)==1){
return 1779;}//B- decays to eta pi0 pi- 

if ( PcheckDecay(genpart,331, -311, -211)==1){
return 1780;}//B- decays to eta' anti-K0 pi- 

if ( PcheckDecay(genpart,331, -321, 111)==1){
return 1781;}//B- decays to eta' K- pi0 

if ( PcheckDecay(genpart,331, 111, -211)==1){
return 1782;}//B- decays to eta' pi0 pi- 

if ( PcheckDecay(genpart,221, 221, -321)==1){
return 1783;}//B- decays to eta eta K- 

if ( PcheckDecay(genpart,221, 221, -211)==1){
return 1784;}//B- decays to eta eta pi- 

if ( PcheckDecay(genpart,221, 331, -321)==1){
return 1785;}//B- decays to eta eta' K- 

if ( PcheckDecay(genpart,221, 331, -211)==1){
return 1786;}//B- decays to eta eta' pi- 

if ( PcheckDecay(genpart,331, 331, -321)==1){
return 1787;}//B- decays to eta' eta' K- 

if ( PcheckDecay(genpart,331, 331, -211)==1){
return 1788;}//B- decays to eta' eta' pi- 

if ( PcheckDecay(genpart,221, -313, -211)==1){
return 1789;}//B- decays to eta anti-K*0 pi- 

if ( PcheckDecay(genpart,221, -311, -213)==1){
return 1790;}//B- decays to eta anti-K0 rho- 

if ( PcheckDecay(genpart,221, -323, 111)==1){
return 1791;}//B- decays to eta K*- pi0 

if ( PcheckDecay(genpart,221, -321, 113)==1){
return 1792;}//B- decays to eta K- rho0 

if ( PcheckDecay(genpart,221, 113, -211)==1){
return 1793;}//B- decays to eta rho0 pi- 

if ( PcheckDecay(genpart,221, 111, -213)==1){
return 1794;}//B- decays to eta pi0 rho- 

if ( PcheckDecay(genpart,331, -313, -211)==1){
return 1795;}//B- decays to eta' anti-K*0 pi- 

if ( PcheckDecay(genpart,331, -311, -213)==1){
return 1796;}//B- decays to eta' anti-K0 rho- 

if ( PcheckDecay(genpart,331, -323, 111)==1){
return 1797;}//B- decays to eta' K*- pi0 

if ( PcheckDecay(genpart,331, -321, 113)==1){
return 1798;}//B- decays to eta' K- rho0 

if ( PcheckDecay(genpart,331, 113, -211)==1){
return 1799;}//B- decays to eta' rho0 pi- 

if ( PcheckDecay(genpart,331, 111, -213)==1){
return 1800;}//B- decays to eta' pi0 rho- 

if ( PcheckDecay(genpart,221, 221, -323)==1){
return 1801;}//B- decays to eta eta K*- 

if ( PcheckDecay(genpart,221, 221, -213)==1){
return 1802;}//B- decays to eta eta rho- 

if ( PcheckDecay(genpart,221, 331, -323)==1){
return 1803;}//B- decays to eta eta' K*- 

if ( PcheckDecay(genpart,221, 331, -213)==1){
return 1804;}//B- decays to eta eta' rho- 

if ( PcheckDecay(genpart,331, 331, -323)==1){
return 1805;}//B- decays to eta' eta' K*- 

if ( PcheckDecay(genpart,331, 331, -213)==1){
return 1806;}//B- decays to eta' eta' rho- 

if ( PcheckDecay(genpart,223, -311, -211)==1){
return 1807;}//B- decays to omega anti-K0 pi- 

if ( PcheckDecay(genpart,223, -321, 111)==1){
return 1808;}//B- decays to omega K- pi0 

if ( PcheckDecay(genpart,223, 111, -211)==1){
return 1809;}//B- decays to omega pi0 pi- 

if ( PcheckDecay(genpart,333, -311, -211)==1){
return 1810;}//B- decays to phi anti-K0 pi- 

if ( PcheckDecay(genpart,333, -321, 111)==1){
return 1811;}//B- decays to phi K- pi0 

if ( PcheckDecay(genpart,333, 111, -211)==1){
return 1812;}//B- decays to phi pi0 pi- 

if ( PcheckDecay(genpart,223, 221, -321)==1){
return 1813;}//B- decays to omega eta K- 

if ( PcheckDecay(genpart,223, 221, -211)==1){
return 1814;}//B- decays to omega eta pi- 

if ( PcheckDecay(genpart,223, 331, -321)==1){
return 1815;}//B- decays to omega eta' K- 

if ( PcheckDecay(genpart,223, 331, -211)==1){
return 1816;}//B- decays to omega eta' pi- 

if ( PcheckDecay(genpart,333, 221, -321)==1){
return 1817;}//B- decays to phi eta K- 

if ( PcheckDecay(genpart,333, 221, -211)==1){
return 1818;}//B- decays to phi eta pi- 

if ( PcheckDecay(genpart,333, 331, -321)==1){
return 1819;}//B- decays to phi eta' K- 

if ( PcheckDecay(genpart,333, 331, -211)==1){
return 1820;}//B- decays to phi eta' pi- 

if ( PcheckDecay(genpart,333, 333, -321)==1){
return 1821;}//B- decays to phi phi K- 

if ( PcheckDecay(genpart,-10211, 221)==1){
return 1822;}//B- decays to a_0- eta 

if ( PcheckDecay(genpart,-10211, 331)==1){
return 1823;}//B- decays to a_0- eta' 

if ( PcheckDecay(genpart,-10211, 100221)==1){
return 1824;}//B- decays to a_0- eta(2S) 

if ( PcheckDecay(genpart,-10211, 9020221)==1){
return 1825;}//B- decays to a_0- eta(1405) 

if ( PcheckDecay(genpart,-10211, 10111)==1){
return 1826;}//B- decays to a_0- a_00 

if ( PcheckDecay(genpart,-10211, 10221)==1){
return 1827;}//B- decays to a_0- f_0 

if ( PcheckDecay(genpart,-323, 10111)==1){
return 1828;}//B- decays to K*- a_00 

if ( PcheckDecay(genpart,-313, -10211)==1){
return 1829;}//B- decays to anti-K*0 a_0- 

if ( PcheckDecay(genpart,-213, 10111)==1){
return 1830;}//B- decays to rho- a_00 

if ( PcheckDecay(genpart,113, -10211)==1){
return 1831;}//B- decays to rho0 a_0- 

if ( PcheckDecay(genpart,223, -10211)==1){
return 1832;}//B- decays to omega a_0- 

if ( PcheckDecay(genpart,333, -10211)==1){
return 1833;}//B- decays to phi a_0- 

if ( PcheckDecay(genpart,20113, -10211)==1){
return 1834;}//B- decays to a_10 a_0- 

if ( PcheckDecay(genpart,-20213, 10111)==1){
return 1835;}//B- decays to a_1- a_00 

if ( PcheckDecay(genpart,20223, -10211)==1){
return 1836;}//B- decays to f_1 a_0- 

if ( PcheckDecay(genpart,10113, -10211)==1){
return 1837;}//B- decays to b_10 a_0- 

if ( PcheckDecay(genpart,-10213, 10111)==1){
return 1838;}//B- decays to b_1- a_00 

if ( PcheckDecay(genpart,10223, -10211)==1){
return 1839;}//B- decays to h_1 a_0- 

if ( PcheckDecay(genpart,-20213, 221)==1){
return 1840;}//B- decays to a_1- eta 

if ( PcheckDecay(genpart,-20213, 331)==1){
return 1841;}//B- decays to a_1- eta' 

if ( PcheckDecay(genpart,-20213, 100221)==1){
return 1842;}//B- decays to a_1- eta(2S) 

if ( PcheckDecay(genpart,-20213, 9020221)==1){
return 1843;}//B- decays to a_1- eta(1405) 

if ( PcheckDecay(genpart,-20213, 10221)==1){
return 1844;}//B- decays to a_1- f_0 

if ( PcheckDecay(genpart,-323, 20113)==1){
return 1845;}//B- decays to K*- a_10 

if ( PcheckDecay(genpart,-313, -20213)==1){
return 1846;}//B- decays to anti-K*0 a_1- 

if ( PcheckDecay(genpart,-213, 20113)==1){
return 1847;}//B- decays to rho- a_10 

if ( PcheckDecay(genpart,113, -20213)==1){
return 1848;}//B- decays to rho0 a_1- 

if ( PcheckDecay(genpart,223, -20213)==1){
return 1849;}//B- decays to omega a_1- 

if ( PcheckDecay(genpart,333, -20213)==1){
return 1850;}//B- decays to phi a_1- 

if ( PcheckDecay(genpart,20113, -20213)==1){
return 1851;}//B- decays to a_10 a_1- 

if ( PcheckDecay(genpart,20223, -20213)==1){
return 1852;}//B- decays to f_1 a_1- 

if ( PcheckDecay(genpart,10113, -20213)==1){
return 1853;}//B- decays to b_10 a_1- 

if ( PcheckDecay(genpart,-10213, 20113)==1){
return 1854;}//B- decays to b_1- a_10 

if ( PcheckDecay(genpart,10223, -20213)==1){
return 1855;}//B- decays to h_1 a_1- 

if ( PcheckDecay(genpart,-10213, 10221)==1){
return 1856;}//B- decays to b_1- f_0 

if ( PcheckDecay(genpart,20223, -321)==1){
return 1857;}//B- decays to f_1 K- 
     
if ( PcheckDecay(genpart,20223, -211)==1){
return 1858;}//B- decays to f_1 pi- 

if ( PcheckDecay(genpart,-323, 20223)==1){
return 1859;}//B- decays to K*- f_1 

if ( PcheckDecay(genpart,-213, 20223)==1){
return 1860;}//B- decays to rho- f_1 

if ( PcheckDecay(genpart,-10213, 20223)==1){
return 1861;}//B- decays to b_1- f_1 

if ( PcheckDecay(genpart,-10213, -311)==1){
return 1862;}//B- decays to b_1- anti-K0 

if ( PcheckDecay(genpart,-10213, 111)==1){
return 1863;}//B- decays to b_1- pi0 

if ( PcheckDecay(genpart,-10213, 221)==1){
return 1864;}//B- decays to b_1- eta 

if ( PcheckDecay(genpart,-10213, 331)==1){
return 1865;}//B- decays to b_1- eta' 

if ( PcheckDecay(genpart,-10213, 100221)==1){
return 1866;}//B- decays to b_1- eta(2S) 

if ( PcheckDecay(genpart,-10213, 9020221)==1){
return 1867;}//B- decays to b_1- eta(1405) 

if ( PcheckDecay(genpart,-323, 10113)==1){
return 1868;}//B- decays to K*- b_10 

if ( PcheckDecay(genpart,-313, -10213)==1){
return 1869;}//B- decays to anti-K*0 b_1- 

if ( PcheckDecay(genpart,-213, 10113)==1){
return 1870;}//B- decays to rho- b_10 

if ( PcheckDecay(genpart,113, -10213)==1){
return 1871;}//B- decays to rho0 b_1- 

if ( PcheckDecay(genpart,223, -10213)==1){
return 1872;}//B- decays to omega b_1- 

if ( PcheckDecay(genpart,333, -10213)==1){
return 1873;}//B- decays to phi b_1- 

if ( PcheckDecay(genpart,10113, -10213)==1){
return 1874;}//B- decays to b_10 b_1- 

if ( PcheckDecay(genpart,10223, -10213)==1){
return 1875;}//B- decays to h_1 b_1- 

if ( PcheckDecay(genpart,10223, -321)==1){
return 1876;}//B- decays to h_1 K- 

if ( PcheckDecay(genpart,10223, -211)==1){
return 1877;}//B- decays to h_1 pi- 

if ( PcheckDecay(genpart,-323, 10223)==1){
return 1878;}//B- decays to K*- h_1 

if ( PcheckDecay(genpart,-213, 10223)==1){
return 1879;}//B- decays to rho- h_1 

if ( PcheckDecay(genpart,115, -211)==1){
return 1880;}//B- decays to a_20 pi- 

if ( PcheckDecay(genpart,335, -211)==1){
return 1881;}//B- decays to f'_2 pi- 

if ( PcheckDecay(genpart,-215, 111)==1){
return 1882;}//B- decays to a_2- pi0 

if ( PcheckDecay(genpart,-215, 221)==1){
return 1883;}//B- decays to a_2- eta 

if ( PcheckDecay(genpart,-215, 331)==1){
return 1884;}//B- decays to a_2- eta' 

if ( PcheckDecay(genpart,-215, 100221)==1){
return 1885;}//B- decays to a_2- eta(2S) 

if ( PcheckDecay(genpart,-215, 9020221)==1){
return 1886;}//B- decays to a_2- eta(1405) 

if ( PcheckDecay(genpart,-325, -311)==1){
return 1887;}//B- decays to K_2*- anti-K0 

if ( PcheckDecay(genpart,315, -321)==1){
return 1888;}//B- decays to K_2*0 K- 

if ( PcheckDecay(genpart,-215, -311)==1){
return 1889;}//B- decays to a_2- anti-K0 

if ( PcheckDecay(genpart,-325, 111)==1){
return 1890;}//B- decays to K_2*- pi0 

if ( PcheckDecay(genpart,-325, 331)==1){
return 1891;}//B- decays to K_2*- eta' 

if ( PcheckDecay(genpart,-325, 100221)==1){
return 1892;}//B- decays to K_2*- eta(2S) 

if ( PcheckDecay(genpart,-325, 9020221)==1){
return 1893;}//B- decays to K_2*- eta(1405) 

if ( PcheckDecay(genpart,115, -213)==1){
return 1894;}//B- decays to a_20 rho- 

if ( PcheckDecay(genpart,225, -213)==1){
return 1895;}//B- decays to f_2 rho- 

if ( PcheckDecay(genpart,335, -213)==1){
return 1896;}//B- decays to f'_2 rho- 

if ( PcheckDecay(genpart,-215, 113)==1){
return 1897;}//B- decays to a_2- rho0 

if ( PcheckDecay(genpart,-215, 223)==1){
return 1898;}//B- decays to a_2- omega 

if ( PcheckDecay(genpart,-215, 333)==1){
return 1899;}//B- decays to a_2- phi 

if ( PcheckDecay(genpart,-325, 313)==1){
return 1900;}//B- decays to K_2*- K*0 

if ( PcheckDecay(genpart,315, -323)==1){
return 1901;}//B- decays to K_2*0 K*- 

if ( PcheckDecay(genpart,115, -323)==1){
return 1902;}//B- decays to a_20 K*- 

if ( PcheckDecay(genpart,-215, -313)==1){
return 1903;}//B- decays to a_2- anti-K*0 

if ( PcheckDecay(genpart,-325, 113)==1){
return 1904;}//B- decays to K_2*- rho0 

if ( PcheckDecay(genpart,-325, 223)==1){
return 1905;}//B- decays to K_2*- omega 

if ( PcheckDecay(genpart,-325, 333)==1){
return 1906;}//B- decays to K_2*- phi 

if ( PcheckDecay(genpart,-315, -213)==1){
return 1907;}//B- decays to anti-K_2*0 rho- 

if ( PcheckDecay(genpart,3222, -2224, 111)==1){
return 1908;}//B- decays to Sigma+ anti-Delta-- pi0 

if ( PcheckDecay(genpart,3122, -2212, 111)==1){
return 1909;}//B- decays to Lambda0 anti-p- pi0 

if ( PcheckDecay(genpart,3212, -2212, 111)==1){
return 1910;}//B- decays to Sigma0 anti-p- pi0 

if ( PcheckDecay(genpart,3112, -2112, 111)==1){
return 1911;}//B- decays to Sigma- anti-n0 pi0 

if ( PcheckDecay(genpart,3122, -2224, 211)==1){
return 1912;}//B- decays to Lambda0 anti-Delta-- pi+ 

if ( PcheckDecay(genpart,3212, -2224, 211)==1){
return 1913;}//B- decays to Sigma0 anti-Delta-- pi+ 

if ( PcheckDecay(genpart,3112, -2212, 211)==1){
return 1914;}//B- decays to Sigma- anti-p- pi+ 

if ( PcheckDecay(genpart,3222, -2212, -211)==1){
return 1915;}//B- decays to Sigma+ anti-p- pi- 

if ( PcheckDecay(genpart,3122, -2112, -211)==1){
return 1916;}//B- decays to Lambda0 anti-n0 pi- 

if ( PcheckDecay(genpart,3212, -2112, -211)==1){
return 1917;}//B- decays to Sigma0 anti-n0 pi- 

if ( PcheckDecay(genpart,2212, -2224, -311)==1){
return 1918;}//B- decays to p+ anti-Delta-- anti-K0 

if ( PcheckDecay(genpart,2112, -2212, -311)==1){
return 1919;}//B- decays to n0 anti-p- anti-K0 

if ( PcheckDecay(genpart,1114, -2112, -311)==1){
return 1920;}//B- decays to Delta- anti-n0 anti-K0 

if ( PcheckDecay(genpart,2224, -2224, -321)==1){
return 1921;}//B- decays to Delta++ anti-Delta-- K- 

if ( PcheckDecay(genpart,2212, -2212, -321)==1){
return 1922;}//B- decays to p+ anti-p- K- 

if ( PcheckDecay(genpart,2112, -2112, -321)==1){
return 1923;}//B- decays to n0 anti-n0 K- 

if ( PcheckDecay(genpart,3322, -3222, 111)==1){
return 1924;}//B- decays to Xi0 anti-Sigma- pi0 

if ( PcheckDecay(genpart,3312, -3122, 111)==1){
return 1925;}//B- decays to Xi- anti-Lambda0 pi0 

if ( PcheckDecay(genpart,3312, -3212, 111)==1){
return 1926;}//B- decays to Xi- anti-Sigma0 pi0 

if ( PcheckDecay(genpart,3312, -3222, 211)==1){
return 1927;}//B- decays to Xi- anti-Sigma- pi+ 

if ( PcheckDecay(genpart,3322, -3122, -211)==1){
return 1928;}//B- decays to Xi0 anti-Lambda0 pi- 

if ( PcheckDecay(genpart,3322, -3212, -211)==1){
return 1929;}//B- decays to Xi0 anti-Sigma0 pi- 

if ( PcheckDecay(genpart,3312, -3112, -211)==1){
return 1930;}//B- decays to Xi- anti-Sigma+ pi- 

if ( PcheckDecay(genpart,3322, -2212, 311)==1){
return 1931;}//B- decays to Xi0 anti-p- K0 

if ( PcheckDecay(genpart,3322, -2224, 321)==1){
return 1932;}//B- decays to Xi0 anti-Delta-- K+ 

if ( PcheckDecay(genpart,3312, -2112, 311)==1){
return 1933;}//B- decays to Xi- anti-n0 K0 

if ( PcheckDecay(genpart,3312, -2212, 321)==1){
return 1934;}//B- decays to Xi- anti-p- K+ 

if ( PcheckDecay(genpart,3122, -3222, -311)==1){
return 1935;}//B- decays to Lambda0 anti-Sigma- anti-K0 

if ( PcheckDecay(genpart,3212, -3222, -311)==1){
return 1936;}//B- decays to Sigma0 anti-Sigma- anti-K0 

if ( PcheckDecay(genpart,3222, -3222, -321)==1){
return 1937;}//B- decays to Sigma+ anti-Sigma- K- 

if ( PcheckDecay(genpart,3112, -3122, -311)==1){
return 1938;}//B- decays to Sigma- anti-Lambda0 anti-K0 

if ( PcheckDecay(genpart,3112, -3212, -311)==1){
return 1939;}//B- decays to Sigma- anti-Sigma0 anti-K0 

if ( PcheckDecay(genpart,3122, -3122, -321)==1){
return 1940;}//B- decays to Lambda0 anti-Lambda0 K- 

if ( PcheckDecay(genpart,3122, -3212, -321)==1){
return 1941;}//B- decays to Lambda0 anti-Sigma0 K- 

if ( PcheckDecay(genpart,3212, -3122, -321)==1){
return 1942;}//B- decays to Sigma0 anti-Lambda0 K- 

if ( PcheckDecay(genpart,3212, -3212, -321)==1){
return 1943;}//B- decays to Sigma0 anti-Sigma0 K- 

if ( PcheckDecay(genpart,3334, -3322, 111)==1){
return 1944;}//B- decays to Omega- anti-Xi0 pi0 

if ( PcheckDecay(genpart,3334, -3312, -211)==1){
return 1945;}//B- decays to Omega- anti-Xi+ pi- 

if ( PcheckDecay(genpart,3334, -3122, 311)==1){
return 1946;}//B- decays to Omega- anti-Lambda0 K0 

if ( PcheckDecay(genpart,3334, -3212, 311)==1){
return 1947;}//B- decays to Omega- anti-Sigma0 K0 

if ( PcheckDecay(genpart,3334, -3222, 321)==1){
return 1948;}//B- decays to Omega- anti-Sigma- K+ 

if ( PcheckDecay(genpart,3312, -3322, -311)==1){
return 1949;}//B- decays to Xi- anti-Xi0 anti-K0 

if ( PcheckDecay(genpart,3322, -3322, -321)==1){
return 1950;}//B- decays to Xi0 anti-Xi0 K- 

if ( PcheckDecay(genpart,3222, -2224, 113)==1){
return 1951;}//B- decays to Sigma+ anti-Delta-- rho0 

if ( PcheckDecay(genpart,3122, -2212, 113)==1){
return 1952;}//B- decays to Lambda0 anti-p- rho0 

if ( PcheckDecay(genpart,3212, -2212, 113)==1){
return 1953;}//B- decays to Sigma0 anti-p- rho0 

if ( PcheckDecay(genpart,3112, -2112, 113)==1){
return 1954;}//B- decays to Sigma- anti-n0 rho0 

if ( PcheckDecay(genpart,3122, -2224, 213)==1){
return 1955;}//B- decays to Lambda0 anti-Delta-- rho+ 

if ( PcheckDecay(genpart,3212, -2224, 213)==1){
return 1956;}//B- decays to Sigma0 anti-Delta-- rho+ 

if ( PcheckDecay(genpart,3112, -2212, 213)==1){
return 1957;}//B- decays to Sigma- anti-p- rho+ 

if ( PcheckDecay(genpart,3222, -2212, -213)==1){
return 1958;}//B- decays to Sigma+ anti-p- rho- 

if ( PcheckDecay(genpart,3122, -2112, -213)==1){
return 1959;}//B- decays to Lambda0 anti-n0 rho- 

if ( PcheckDecay(genpart,3212, -2112, -213)==1){
return 1960;}//B- decays to Sigma0 anti-n0 rho- 

if ( PcheckDecay(genpart,2212, -2224, -313)==1){
return 1961;}//B- decays to p+ anti-Delta-- anti-K*0 

if ( PcheckDecay(genpart,2112, -2212, -313)==1){
return 1962;}//B- decays to n0 anti-p- anti-K*0 

if ( PcheckDecay(genpart,1114, -2112, -313)==1){
return 1963;}//B- decays to Delta- anti-n0 anti-K*0 

if ( PcheckDecay(genpart,2224, -2224, -323)==1){
return 1964;}//B- decays to Delta++ anti-Delta-- K*- 

if ( PcheckDecay(genpart,2212, -2212, -323)==1){
return 1965;}//B- decays to p+ anti-p- K*- 

if ( PcheckDecay(genpart,2112, -2112, -323)==1){
return 1966;}//B- decays to n0 anti-n0 K*- 

if ( PcheckDecay(genpart,3322, -3222, 113)==1){
return 1967;}//B- decays to Xi0 anti-Sigma- rho0 

if ( PcheckDecay(genpart,3312, -3122, 113)==1){
return 1968;}//B- decays to Xi- anti-Lambda0 rho0 

if ( PcheckDecay(genpart,3312, -3212, 113)==1){
return 1969;}//B- decays to Xi- anti-Sigma0 rho0 

if ( PcheckDecay(genpart,3312, -3222, 213)==1){
return 1970;}//B- decays to Xi- anti-Sigma- rho+ 

if ( PcheckDecay(genpart,3322, -3122, -213)==1){
return 1971;}//B- decays to Xi0 anti-Lambda0 rho- 

if ( PcheckDecay(genpart,3322, -3212, -213)==1){
return 1972;}//B- decays to Xi0 anti-Sigma0 rho- 

if ( PcheckDecay(genpart,3312, -3112, -213)==1){
return 1973;}//B- decays to Xi- anti-Sigma+ rho- 

if ( PcheckDecay(genpart,3322, -2212, 313)==1){
return 1974;}//B- decays to Xi0 anti-p- K*0 

if ( PcheckDecay(genpart,3322, -2224, 323)==1){
return 1975;}//B- decays to Xi0 anti-Delta-- K*+ 

if ( PcheckDecay(genpart,3312, -2112, 313)==1){
return 1976;}//B- decays to Xi- anti-n0 K*0 

if ( PcheckDecay(genpart,3312, -2212, 323)==1){
return 1977;}//B- decays to Xi- anti-p- K*+ 

if ( PcheckDecay(genpart,3122, -3222, -313)==1){
return 1978;}//B- decays to Lambda0 anti-Sigma- anti-K*0 

if ( PcheckDecay(genpart,3212, -3222, -313)==1){
return 1979;}//B- decays to Sigma0 anti-Sigma- anti-K*0 

if ( PcheckDecay(genpart,3222, -3222, -323)==1){
return 1980;}//B- decays to Sigma+ anti-Sigma- K*- 

if ( PcheckDecay(genpart,3112, -3122, -313)==1){
return 1981;}//B- decays to Sigma- anti-Lambda0 anti-K*0 

if ( PcheckDecay(genpart,3112, -3212, -313)==1){
return 1982;}//B- decays to Sigma- anti-Sigma0 anti-K*0 

if ( PcheckDecay(genpart,3122, -3122, -323)==1){
return 1983;}//B- decays to Lambda0 anti-Lambda0 K*- 

if ( PcheckDecay(genpart,3122, -3212, -323)==1){
return 1984;}//B- decays to Lambda0 anti-Sigma0 K*- 

if ( PcheckDecay(genpart,3212, -3122, -323)==1){
return 1985;}//B- decays to Sigma0 anti-Lambda0 K*- 

if ( PcheckDecay(genpart,3212, -3212, -323)==1){
return 1986;}//B- decays to Sigma0 anti-Sigma0 K*- 

if ( PcheckDecay(genpart,3334, -3322, 113)==1){
return 1987;}//B- decays to Omega- anti-Xi0 rho0 

if ( PcheckDecay(genpart,3334, -3312, -213)==1){
return 1988;}//B- decays to Omega- anti-Xi+ rho- 

if ( PcheckDecay(genpart,3334, -3122, 313)==1){
return 1989;}//B- decays to Omega- anti-Lambda0 K*0 

if ( PcheckDecay(genpart,3334, -3212, 313)==1){
return 1990;}//B- decays to Omega- anti-Sigma0 K*0 

if ( PcheckDecay(genpart,3334, -3222, 323)==1){
return 1991;}//B- decays to Omega- anti-Sigma- K*+ 

if ( PcheckDecay(genpart,3312, -3322, -313)==1){
return 1992;}//B- decays to Xi- anti-Xi0 anti-K*0 

if ( PcheckDecay(genpart,3322, -3322, -323)==1){
return 1993;}//B- decays to Xi0 anti-Xi0 K*- 

if ( PcheckDecay(genpart,2212, -2224, 111)==1){
return 1994;}//B- decays to p+ anti-Delta-- pi0 

if ( PcheckDecay(genpart,2112, -2212, 111)==1){
return 1995;}//B- decays to n0 anti-p- pi0 

if ( PcheckDecay(genpart,1114, -2112, 111)==1){
return 1996;}//B- decays to Delta- anti-n0 pi0 

if ( PcheckDecay(genpart,2212, -2212, -211)==1){
return 1997;}//B- decays to p+ anti-p- pi- 

if ( PcheckDecay(genpart,2112, -2224, 211)==1){
return 1998;}//B- decays to n0 anti-Delta-- pi+ 

if ( PcheckDecay(genpart,2224, -2224, -211)==1){
return 1999;}//B- decays to Delta++ anti-Delta-- pi- 

if ( PcheckDecay(genpart,2112, -2112, -211)==1){
return 2000;}//B- decays to n0 anti-n0 pi- 

if ( PcheckDecay(genpart,1114, -2212, 211)==1){
return 2001;}//B- decays to Delta- anti-p- pi+ 

if ( PcheckDecay(genpart,1114, -1114, -211)==1){
return 2002;}//B- decays to Delta- anti-Delta+ pi- 

if ( PcheckDecay(genpart,3122, -3222, 111)==1){
return 2003;}//B- decays to Lambda0 anti-Sigma- pi0 

if ( PcheckDecay(genpart,3212, -3222, 111)==1){
return 2004;}//B- decays to Sigma0 anti-Sigma- pi0 

if ( PcheckDecay(genpart,3112, -3122, 111)==1){
return 2005;}//B- decays to Sigma- anti-Lambda0 pi0 

if ( PcheckDecay(genpart,3112, -3212, 111)==1){
return 2006;}//B- decays to Sigma- anti-Sigma0 pi0 

if ( PcheckDecay(genpart,3122, -3122, -211)==1){
return 2007;}//B- decays to Lambda0 anti-Lambda0 pi- 

if ( PcheckDecay(genpart,3122, -3212, -211)==1){
return 2008;}//B- decays to Lambda0 anti-Sigma0 pi- 

if ( PcheckDecay(genpart,3212, -3122, -211)==1){
return 2009;}//B- decays to Sigma0 anti-Lambda0 pi- 

if ( PcheckDecay(genpart,3212, -3212, -211)==1){
return 2010;}//B- decays to Sigma0 anti-Sigma0 pi- 

if ( PcheckDecay(genpart,3112, -3222, 211)==1){
return 2011;}//B- decays to Sigma- anti-Sigma- pi+ 

if ( PcheckDecay(genpart,3222, -3222, -211)==1){
return 2012;}//B- decays to Sigma+ anti-Sigma- pi- 

if ( PcheckDecay(genpart,3112, -3112, -211)==1){
return 2013;}//B- decays to Sigma- anti-Sigma+ pi- 

if ( PcheckDecay(genpart,3122, -2212, 311)==1){
return 2014;}//B- decays to Lambda0 anti-p- K0 

if ( PcheckDecay(genpart,3212, -2212, 311)==1){
return 2015;}//B- decays to Sigma0 anti-p- K0 

if ( PcheckDecay(genpart,3122, -2224, 321)==1){
return 2016;}//B- decays to Lambda0 anti-Delta-- K+ 

if ( PcheckDecay(genpart,3212, -2224, 321)==1){
return 2017;}//B- decays to Sigma0 anti-Delta-- K+ 

if ( PcheckDecay(genpart,3112, -2112, 311)==1){
return 2018;}//B- decays to Sigma- anti-n0 K0 

if ( PcheckDecay(genpart,3112, -2214, 321)==1){
return 2019;}//B- decays to Sigma- anti-Delta- K+ 

if ( PcheckDecay(genpart,2112, -3222, -311)==1){
return 2020;}//B- decays to n0 anti-Sigma- anti-K0 

if ( PcheckDecay(genpart,2212, -3222, -321)==1){
return 2021;}//B- decays to p+ anti-Sigma- K- 

if ( PcheckDecay(genpart,1114, -3122, -311)==1){
return 2022;}//B- decays to Delta- anti-Lambda0 anti-K0 

if ( PcheckDecay(genpart,1114, -3212, -311)==1){
return 2023;}//B- decays to Delta- anti-Sigma0 anti-K0 

if ( PcheckDecay(genpart,2112, -3122, -321)==1){
return 2024;}//B- decays to n0 anti-Lambda0 K- 

if ( PcheckDecay(genpart,2112, -3212, -321)==1){
return 2025;}//B- decays to n0 anti-Sigma0 K- 

if ( PcheckDecay(genpart,3312, -3322, 111)==1){
return 2026;}//B- decays to Xi- anti-Xi0 pi0 

if ( PcheckDecay(genpart,3312, -3312, -211)==1){
return 2027;}//B- decays to Xi- anti-Xi+ pi- 

if ( PcheckDecay(genpart,3322, -3322, -211)==1){
return 2028;}//B- decays to Xi0 anti-Xi0 pi- 

if ( PcheckDecay(genpart,3312, -3122, 311)==1){
return 2029;}//B- decays to Xi- anti-Lambda0 K0 

if ( PcheckDecay(genpart,3312, -3212, 311)==1){
return 2030;}//B- decays to Xi- anti-Sigma0 K0 

if ( PcheckDecay(genpart,3312, -3222, 321)==1){
return 2031;}//B- decays to Xi- anti-Sigma- K+ 

if ( PcheckDecay(genpart,3112, -3322, -311)==1){
return 2032;}//B- decays to Sigma- anti-Xi0 anti-K0 

if ( PcheckDecay(genpart,3122, -3322, -321)==1){
return 2033;}//B- decays to Lambda0 anti-Xi0 K- 

if ( PcheckDecay(genpart,3212, -3322, -321)==1){
return 2034;}//B- decays to Sigma0 anti-Xi0 K- 

if ( PcheckDecay(genpart,2212, -2224, 113)==1){
return 2035;}//B- decays to p+ anti-Delta-- rho0 

if ( PcheckDecay(genpart,2112, -2212, 113)==1){
return 2036;}//B- decays to n0 anti-p- rho0 

if ( PcheckDecay(genpart,1114, -2112, 113)==1){
return 2037;}//B- decays to Delta- anti-n0 rho0 

if ( PcheckDecay(genpart,2212, -2212, -213)==1){
return 2038;}//B- decays to p+ anti-p- rho- 

if ( PcheckDecay(genpart,2112, -2224, 213)==1){
return 2039;}//B- decays to n0 anti-Delta-- rho+ 

if ( PcheckDecay(genpart,2224, -2224, -213)==1){
return 2040;}//B- decays to Delta++ anti-Delta-- rho- 

if ( PcheckDecay(genpart,2112, -2112, -213)==1){
return 2041;}//B- decays to n0 anti-n0 rho- 

if ( PcheckDecay(genpart,1114, -2212, 213)==1){
return 2042;}//B- decays to Delta- anti-p- rho+ 

if ( PcheckDecay(genpart,1114, -1114, -213)==1){
return 2043;}//B- decays to Delta- anti-Delta+ rho- 

if ( PcheckDecay(genpart,3122, -3222, 113)==1){
return 2044;}//B- decays to Lambda0 anti-Sigma- rho0 

if ( PcheckDecay(genpart,3212, -3222, 113)==1){
return 2045;}//B- decays to Sigma0 anti-Sigma- rho0 

if ( PcheckDecay(genpart,3112, -3122, 113)==1){
return 2046;}//B- decays to Sigma- anti-Lambda0 rho0 

if ( PcheckDecay(genpart,3112, -3212, 113)==1){
return 2047;}//B- decays to Sigma- anti-Sigma0 rho0 

if ( PcheckDecay(genpart,3122, -3122, -213)==1){
return 2048;}//B- decays to Lambda0 anti-Lambda0 rho- 

if ( PcheckDecay(genpart,3122, -3212, -213)==1){
return 2049;}//B- decays to Lambda0 anti-Sigma0 rho- 

if ( PcheckDecay(genpart,3212, -3122, -213)==1){
return 2050;}//B- decays to Sigma0 anti-Lambda0 rho- 

if ( PcheckDecay(genpart,3212, -3212, -213)==1){
return 2051;}//B- decays to Sigma0 anti-Sigma0 rho- 

if ( PcheckDecay(genpart,3112, -3222, 213)==1){
return 2052;}//B- decays to Sigma- anti-Sigma- rho+ 

if ( PcheckDecay(genpart,3222, -3222, -213)==1){
return 2053;}//B- decays to Sigma+ anti-Sigma- rho- 

if ( PcheckDecay(genpart,3112, -3112, -213)==1){
return 2054;}//B- decays to Sigma- anti-Sigma+ rho- 

if ( PcheckDecay(genpart,3122, -2212, 313)==1){
return 2055;}//B- decays to Lambda0 anti-p- K*0 

if ( PcheckDecay(genpart,3212, -2212, 313)==1){
return 2056;}//B- decays to Sigma0 anti-p- K*0 

if ( PcheckDecay(genpart,3122, -2224, 323)==1){
return 2057;}//B- decays to Lambda0 anti-Delta-- K*+ 

if ( PcheckDecay(genpart,3212, -2224, 323)==1){
return 2058;}//B- decays to Sigma0 anti-Delta-- K*+ 

if ( PcheckDecay(genpart,3112, -2112, 313)==1){
return 2059;}//B- decays to Sigma- anti-n0 K*0 

if ( PcheckDecay(genpart,3112, -2214, 323)==1){
return 2060;}//B- decays to Sigma- anti-Delta- K*+ 

if ( PcheckDecay(genpart,2112, -3222, -313)==1){
return 2061;}//B- decays to n0 anti-Sigma- anti-K*0 

if ( PcheckDecay(genpart,2212, -3222, -323)==1){
return 2062;}//B- decays to p+ anti-Sigma- K*- 

if ( PcheckDecay(genpart,1114, -3122, -313)==1){
return 2063;}//B- decays to Delta- anti-Lambda0 anti-K*0 

if ( PcheckDecay(genpart,1114, -3212, -313)==1){
return 2064;}//B- decays to Delta- anti-Sigma0 anti-K*0 

if ( PcheckDecay(genpart,2112, -3122, -323)==1){
return 2065;}//B- decays to n0 anti-Lambda0 K*- 

if ( PcheckDecay(genpart,2112, -3212, -323)==1){
return 2066;}//B- decays to n0 anti-Sigma0 K*- 

if ( PcheckDecay(genpart,3312, -3322, 113)==1){
return 2067;}//B- decays to Xi- anti-Xi0 rho0 

if ( PcheckDecay(genpart,3312, -3312, -213)==1){
return 2068;}//B- decays to Xi- anti-Xi+ rho- 

if ( PcheckDecay(genpart,3322, -3322, -213)==1){
return 2069;}//B- decays to Xi0 anti-Xi0 rho- 

if ( PcheckDecay(genpart,3312, -3122, 313)==1){
return 2070;}//B- decays to Xi- anti-Lambda0 K*0 

if ( PcheckDecay(genpart,3312, -3212, 313)==1){
return 2071;}//B- decays to Xi- anti-Sigma0 K*0 

if ( PcheckDecay(genpart,3312, -3222, 323)==1){
return 2072;}//B- decays to Xi- anti-Sigma- K*+ 

if ( PcheckDecay(genpart,3112, -3322, -313)==1){
return 2073;}//B- decays to Sigma- anti-Xi0 anti-K*0 

if ( PcheckDecay(genpart,3122, -3322, -323)==1){
return 2074;}//B- decays to Lambda0 anti-Xi0 K*- 

if ( PcheckDecay(genpart,3212, -3322, -323)==1){
return 2075;}//B- decays to Sigma0 anti-Xi0 K*- 

if ( PcheckDecay(genpart,3222, -2224)==1){
return 2076;}//B- decays to Sigma+ anti-Delta-- 

if ( PcheckDecay(genpart,3122, -2212)==1){
return 2077;}//B- decays to Lambda0 anti-p- 

if ( PcheckDecay(genpart,3122, -2214)==1){
return 2078;}//B- decays to Lambda0 anti-Delta- 

if ( PcheckDecay(genpart,3212, -2212)==1){
return 2079;}//B- decays to Sigma0 anti-p- 

if ( PcheckDecay(genpart,3214, -2212)==1){
return 2080;}//B- decays to Sigma*0 anti-p- 

if ( PcheckDecay(genpart,3112, -2112)==1){
return 2081;}//B- decays to Sigma- anti-n0 

if ( PcheckDecay(genpart,3322, -3222)==1){
return 2082;}//B- decays to Xi0 anti-Sigma- 

if ( PcheckDecay(genpart,3312, -3122)==1){
return 2083;}//B- decays to Xi- anti-Lambda0 

if ( PcheckDecay(genpart,3312, -3212)==1){
return 2084;}//B- decays to Xi- anti-Sigma0 

if ( PcheckDecay(genpart,3334, -3322)==1){
return 2085;}//B- decays to Omega- anti-Xi0 

if ( PcheckDecay(genpart,2212, -2224)==1){
return 2086;}//B- decays to p+ anti-Delta-- 

if ( PcheckDecay(genpart,2112, -2212)==1){
return 2087;}//B- decays to n0 anti-p- 

if ( PcheckDecay(genpart,1114, -2112)==1){
return 2088;}//B- decays to Delta- anti-n0 

if ( PcheckDecay(genpart,3122, -3222)==1){
return 2089;}//B- decays to Lambda0 anti-Sigma- 

if ( PcheckDecay(genpart,3212, -3222)==1){
return 2090;}//B- decays to Sigma0 anti-Sigma- 

if ( PcheckDecay(genpart,3112, -3122)==1){
return 2091;}//B- decays to Sigma- anti-Lambda0 

if ( PcheckDecay(genpart,3112, -3212)==1){
return 2092;}//B- decays to Sigma- anti-Sigma0 

if ( PcheckDecay(genpart,3312, -3322)==1){
return 2093;}//B- decays to Xi- anti-Xi0 

if ( PcheckDecay(genpart,3124, -2212)==1){
return 2094;}//B- decays to Lambda(1520)0 anti-p- 

if ( PcheckDecay(genpart,-431, 111)==1){
return 2095;}//B- decays to D_s- pi0 

if ( PcheckDecay(genpart,-431, 221)==1){
return 2096;}//B- decays to D_s- eta 

if ( PcheckDecay(genpart,-431, 331)==1){
return 2097;}//B- decays to D_s- eta' 

if ( PcheckDecay(genpart,113, -431)==1){
return 2098;}//B- decays to rho0 D_s- 

if ( PcheckDecay(genpart,-431, 211, -211)==1){
return 2099;}//B- decays to D_s- pi+ pi- 

if ( PcheckDecay(genpart,-431, 111, 111)==1){
return 2100;}//B- decays to D_s- pi0 pi0 

if ( PcheckDecay(genpart,20113, -431)==1){
return 2101;}//B- decays to a_10 D_s- 

if ( PcheckDecay(genpart,223, -431)==1){
return 2102;}//B- decays to omega D_s- 

if ( PcheckDecay(genpart,-433, 111)==1){
return 2103;}//B- decays to D_s*- pi0 

if ( PcheckDecay(genpart,-433, 221)==1){
return 2104;}//B- decays to D_s*- eta 

if ( PcheckDecay(genpart,-433, 331)==1){
return 2105;}//B- decays to D_s*- eta' 

if ( PcheckDecay(genpart,-433, 113)==1){
return 2106;}//B- decays to D_s*- rho0 

if ( PcheckDecay(genpart,-433, 211, -211)==1){
return 2107;}//B- decays to D_s*- pi+ pi- 

if ( PcheckDecay(genpart,-433, 111, 111)==1){
return 2108;}//B- decays to D_s*- pi0 pi0 

if ( PcheckDecay(genpart,-433, 20113)==1){
return 2109;}//B- decays to D_s*- a_10 

if ( PcheckDecay(genpart,-433, 223)==1){
return 2110;}//B- decays to D_s*- omega 

if ( PcheckDecay(genpart,-411, 111)==1){
return 2111;}//B- decays to D- pi0 

if ( PcheckDecay(genpart,-411, 221)==1){
return 2112;}//B- decays to D- eta 

if ( PcheckDecay(genpart,-411, 331)==1){
return 2113;}//B- decays to D- eta' 

if ( PcheckDecay(genpart,113, -411)==1){
return 2114;}//B- decays to rho0 D- 

if ( PcheckDecay(genpart,-411, 211, -211)==1){
return 2115;}//B- decays to D- pi+ pi- 

if ( PcheckDecay(genpart,-411, 111, 111)==1){
return 2116;}//B- decays to D- pi0 pi0 

if ( PcheckDecay(genpart,20113, -411)==1){
return 2117;}//B- decays to a_10 D- 

if ( PcheckDecay(genpart,223, -411)==1){
return 2118;}//B- decays to omega D- 

if ( PcheckDecay(genpart,-413, 111)==1){
return 2119;}//B- decays to D*- pi0 

if ( PcheckDecay(genpart,-413, 221)==1){
return 2120;}//B- decays to D*- eta 

if ( PcheckDecay(genpart,-413, 331)==1){
return 2121;}//B- decays to D*- eta' 

if ( PcheckDecay(genpart,-413, 113)==1){
return 2122;}//B- decays to D*- rho0 

if ( PcheckDecay(genpart,-413, 211, -211)==1){
return 2123;}//B- decays to D*- pi+ pi- 

if ( PcheckDecay(genpart,-413, 111, 111)==1){
return 2124;}//B- decays to D*- pi0 pi0 

if ( PcheckDecay(genpart,-413, 20113)==1){
return 2125;}//B- decays to D*- a_10 

if ( PcheckDecay(genpart,-413, 223)==1){
return 2126;}//B- decays to D*- omega 

if ( PcheckDecay(genpart,-411, -311)==1){
return 2127;}//B- decays to D- anti-K0 

if ( PcheckDecay(genpart,-313, -411)==1){
return 2128;}//B- decays to anti-K*0 D- 

if ( PcheckDecay(genpart,-411, -321, 211)==1){
return 2129;}//B- decays to D- K- pi+ 

if ( PcheckDecay(genpart,-411, -311, 111)==1){
return 2130;}//B- decays to D- anti-K0 pi0 

if ( PcheckDecay(genpart,-413, -311)==1){
return 2131;}//B- decays to D*- anti-K0 

if ( PcheckDecay(genpart,-413, -313)==1){
return 2132;}//B- decays to D*- anti-K*0 

if ( PcheckDecay(genpart,-413, -321, 211)==1){
return 2133;}//B- decays to D*- K- pi+ 

if ( PcheckDecay(genpart,-413, -311, 111)==1){
return 2134;}//B- decays to D*- anti-K0 pi0 

if ( PcheckDecay(genpart,333, -431)==1){
return 2135;}//B- decays to phi D_s- 

if ( PcheckDecay(genpart,333, -433)==1){
return 2136;}//B- decays to phi D_s*- 

if ( PcheckDecay(genpart,-431, 311)==1){
return 2137;}//B- decays to D_s- K0 

if ( PcheckDecay(genpart,-431, 321, -211)==1){
return 2138;}//B- decays to D_s- K+ pi- 

if ( PcheckDecay(genpart,-431, 311, 111)==1){
return 2139;}//B- decays to D_s- K0 pi0 

if ( PcheckDecay(genpart,313, -431)==1){
return 2140;}//B- decays to K*0 D_s- 

if ( PcheckDecay(genpart,-433, 311)==1){
return 2141;}//B- decays to D_s*- K0 

if ( PcheckDecay(genpart,-433, 321, -211)==1){
return 2142;}//B- decays to D_s*- K+ pi- 

if ( PcheckDecay(genpart,-433, 311, 111)==1){
return 2143;}//B- decays to D_s*- K0 pi0 

if ( PcheckDecay(genpart,-433, 313)==1){
return 2144;}//B- decays to D_s*- K*0 

if ( PcheckDecay(genpart,-433, 22)==1){
return 2145;}//B- decays to D_s*- gamma 

if ( PcheckDecay(genpart,-413, 22)==1){
return 2146;}//B- decays to D*- gamma 

if ( PcheckDecay(genpart,13, -14)==1){
return 2147;}//B- decays to mu- anti-nu_mu 

if ( PcheckDecay(genpart,15, -16)==1){
return 2148;}//B- decays to tau- anti-nu_tau 

if ( PcheckDecay(genpart,11, -12, 22)==1){
return 2149;}//B- decays to e- anti-nu_e gamma 

 if ( PcheckDecay(genpart,13, -14, 22)==1){
   return 2150;}//B- decays to mu- anti-nu_mu gamma 
 return -400;
  }

  else return -9876789;
}

  int  genT::MyownB(Gen_hepevt genpart)
  {
    
    if(abs(AmId(genpart))==521){
      if ( AbsPcheckDecay(genpart,10323,22)==1){
	return 1;}//anti-B0 decays to D*+ e- anti-nu_e 
      if ( AbsPcheckDecay(genpart,100323,22)==1){
	return 2;}//anti-B0 decays to D*+ e- anti-nu_e 
      if ( AbsPcheckDecay(genpart,30323,22)==1){
	return 3;}//anti-B0 decays to D*+ e- anti-nu_e 
      return -9;
    }

    if(abs(AmId(genpart))==511){
      if ( AbsPcheckDecay(genpart,10313,22)==1){
	return 11;}//anti-B0 decays to D*+ e- anti-nu_e 
      if ( AbsPcheckDecay(genpart,100313,22)==1){
	return 12;}//anti-B0 decays to D*+ e- anti-nu_e 
      if ( AbsPcheckDecay(genpart,30313,22)==1){
	return 13;}//anti-B0 decays to D*+ e- anti-nu_e 
      return -99;
    }
  }

  int  genT::Bmode(Gen_hepevt genpart)
  {
    
    if(AmId(genpart)==-511){
      if ( PcheckDecay(genpart,413, 11, -12)==1){
	return 1;}//anti-B0 decays to D*+ e- anti-nu_e 
      
      if ( PcheckDecay(genpart,411, 11, -12)==1){
	return 2;}//anti-B0 decays to D+ e- anti-nu_e 
      
      if ( PcheckDecay(genpart,10413, 11, -12)==1){
	return 3;}//anti-B0 decays to D_1+ e- anti-nu_e 
      
      if ( PcheckDecay(genpart,10411, 11, -12)==1){
	return 4;}//anti-B0 decays to D_0*+ e- anti-nu_e 
      
      if ( PcheckDecay(genpart,20413, 11, -12)==1){
	return 5;}//anti-B0 decays to D'_1+ e- anti-nu_e 

if ( PcheckDecay(genpart,415, 11, -12)==1){
return 6;}//anti-B0 decays to D_2*+ e- anti-nu_e 

if ( PcheckDecay(genpart,100411, 11, -12)==1){
return 7;}//anti-B0 decays to D(2S)+ e- anti-nu_e 

if ( PcheckDecay(genpart,100413, 11, -12)==1){
return 8;}//anti-B0 decays to D*(2S)+ e- anti-nu_e 

if ( PcheckDecay(genpart,413, 111, 11, -12)==1){
return 9;}//anti-B0 decays to D*+ pi0 e- anti-nu_e 

if ( PcheckDecay(genpart,423, 211, 11, -12)==1){
return 10;}//anti-B0 decays to D*0 pi+ e- anti-nu_e 

if ( PcheckDecay(genpart,411, 111, 11, -12)==1){
return 11;}//anti-B0 decays to D+ pi0 e- anti-nu_e 

if ( PcheckDecay(genpart,421, 211, 11, -12)==1){
return 12;}//anti-B0 decays to D0 pi+ e- anti-nu_e 

if ( PcheckDecay(genpart,413, 13, -14)==1){
return 13;}//anti-B0 decays to D*+ mu- anti-nu_mu 

if ( PcheckDecay(genpart,411, 13, -14)==1){
return 14;}//anti-B0 decays to D+ mu- anti-nu_mu 

if ( PcheckDecay(genpart,10413, 13, -14)==1){
return 15;}//anti-B0 decays to D_1+ mu- anti-nu_mu 

if ( PcheckDecay(genpart,10411, 13, -14)==1){
return 16;}//anti-B0 decays to D_0*+ mu- anti-nu_mu 

if ( PcheckDecay(genpart,20413, 13, -14)==1){
return 17;}//anti-B0 decays to D'_1+ mu- anti-nu_mu 

if ( PcheckDecay(genpart,415, 13, -14)==1){
return 18;}//anti-B0 decays to D_2*+ mu- anti-nu_mu 

if ( PcheckDecay(genpart,100411, 13, -14)==1){
return 19;}//anti-B0 decays to D(2S)+ mu- anti-nu_mu 

if ( PcheckDecay(genpart,100413, 13, -14)==1){
return 20;}//anti-B0 decays to D*(2S)+ mu- anti-nu_mu 

if ( PcheckDecay(genpart,413, 111, 13, -14)==1){
return 21;}//anti-B0 decays to D*+ pi0 mu- anti-nu_mu 

if ( PcheckDecay(genpart,423, 211, 13, -14)==1){
return 22;}//anti-B0 decays to D*0 pi+ mu- anti-nu_mu 

if ( PcheckDecay(genpart,411, 111, 13, -14)==1){
return 23;}//anti-B0 decays to D+ pi0 mu- anti-nu_mu 

if ( PcheckDecay(genpart,421, 211, 13, -14)==1){
return 24;}//anti-B0 decays to D0 pi+ mu- anti-nu_mu 

if ( PcheckDecay(genpart,413, 15, -16)==1){
return 25;}//anti-B0 decays to D*+ tau- anti-nu_tau 

if ( PcheckDecay(genpart,411, 15, -16)==1){
return 26;}//anti-B0 decays to D+ tau- anti-nu_tau 

if ( PcheckDecay(genpart,10413, 15, -16)==1){
return 27;}//anti-B0 decays to D_1+ tau- anti-nu_tau 

if ( PcheckDecay(genpart,10411, 15, -16)==1){
return 28;}//anti-B0 decays to D_0*+ tau- anti-nu_tau 

if ( PcheckDecay(genpart,20413, 15, -16)==1){
return 29;}//anti-B0 decays to D'_1+ tau- anti-nu_tau 

if ( PcheckDecay(genpart,415, 15, -16)==1){
return 30;}//anti-B0 decays to D_2*+ tau- anti-nu_tau 

if ( PcheckDecay(genpart,441, 310)==1){
return 31;}//anti-B0 decays to eta_c K_S0 

if ( PcheckDecay(genpart,441, 130)==1){
return 32;}//anti-B0 decays to eta_c K_L0 

if ( PcheckDecay(genpart,-313, 441)==1){
return 33;}//anti-B0 decays to anti-K*S eta_c 

if ( PcheckDecay(genpart,-313, 441)==1){
return 34;}//anti-B0 decays to anti-K*L eta_c 

if ( PcheckDecay(genpart,-313, 441)==1){
return 35;}//anti-B0 decays to anti-K*0T eta_c 

if ( PcheckDecay(genpart,441, -321, 211)==1){
return 36;}//anti-B0 decays to eta_c K- pi+ 

if ( PcheckDecay(genpart,441, -311, 111)==1){
return 37;}//anti-B0 decays to eta_c anti-K0 pi0 

if ( PcheckDecay(genpart,441, -311, 211, -211)==1){
return 38;}//anti-B0 decays to eta_c anti-K0 pi+ pi- 

if ( PcheckDecay(genpart,441, -311, 111, 111)==1){
return 39;}//anti-B0 decays to eta_c anti-K0 pi0 pi0 

if ( PcheckDecay(genpart,441, -321, 211, 111)==1){
return 40;}//anti-B0 decays to eta_c K- pi+ pi0 

if ( PcheckDecay(genpart,441, 3122, -2112)==1){
return 41;}//anti-B0 decays to eta_c Lambda0 anti-n0 

if ( PcheckDecay(genpart,441, 3222, -2212)==1){
return 42;}//anti-B0 decays to eta_c Sigma+ anti-p- 

if ( PcheckDecay(genpart,441, -30343)==1){
return 43;}//anti-B0 decays to eta_c anti-Xsd 

if ( PcheckDecay(genpart,100441, 310)==1){
return 44;}//anti-B0 decays to eta_c(2S) K_S0 

if ( PcheckDecay(genpart,100441, 130)==1){
return 45;}//anti-B0 decays to eta_c(2S) K_L0 

if ( PcheckDecay(genpart,-313, 100441)==1){
return 46;}//anti-B0 decays to anti-K*S eta_c(2S) 

if ( PcheckDecay(genpart,-313, 100441)==1){
return 47;}//anti-B0 decays to anti-K*L eta_c(2S) 

if ( PcheckDecay(genpart,-313, 100441)==1){
return 48;}//anti-B0 decays to anti-K*0T eta_c(2S) 

if ( PcheckDecay(genpart,100441, -321, 211)==1){
return 49;}//anti-B0 decays to eta_c(2S) K- pi+ 

if ( PcheckDecay(genpart,100441, -311, 111)==1){
return 50;}//anti-B0 decays to eta_c(2S) anti-K0 pi0 

if ( PcheckDecay(genpart,100441, -311, 211, -211)==1){
return 51;}//anti-B0 decays to eta_c(2S) anti-K0 pi+ pi- 

if ( PcheckDecay(genpart,100441, -311, 111, 111)==1){
return 52;}//anti-B0 decays to eta_c(2S) anti-K0 pi0 pi0 

if ( PcheckDecay(genpart,100441, -321, 211, 111)==1){
return 53;}//anti-B0 decays to eta_c(2S) K- pi+ pi0 

if ( PcheckDecay(genpart,100441, -30343)==1){
return 54;}//anti-B0 decays to eta_c(2S) anti-Xsd 

if ( PcheckDecay(genpart,443, 310)==1){
return 55;}//anti-B0 decays to psi K_S0 

if ( PcheckDecay(genpart,443, 130)==1){
return 56;}//anti-B0 decays to psi K_L0 

if ( PcheckDecay(genpart,443, -313)==1){
return 57;}//anti-B0 decays to psi anti-K*S 

if ( PcheckDecay(genpart,443, -313)==1){
return 58;}//anti-B0 decays to psi anti-K*L 

if ( PcheckDecay(genpart,443, -313)==1){
return 59;}//anti-B0 decays to psi anti-K*0T 

if ( PcheckDecay(genpart,443, -321, 211)==1){
return 60;}//anti-B0 decays to psi K- pi+ 

if ( PcheckDecay(genpart,443, -10313)==1){
return 61;}//anti-B0 decays to psi anti-K_10 

if ( PcheckDecay(genpart,443, -20313)==1){
return 62;}//anti-B0 decays to psi anti-K'_10 

if ( PcheckDecay(genpart,443, -311, 113)==1){
return 63;}//anti-B0 decays to psi anti-K0 rho0 

if ( PcheckDecay(genpart,443, -311, 211, -211)==1){
return 64;}//anti-B0 decays to psi anti-K0 pi+ pi- 

if ( PcheckDecay(genpart,443, -323, 211)==1){
return 65;}//anti-B0 decays to psi K*- pi+ 

if ( PcheckDecay(genpart,443, -313, 211, -211)==1){
return 66;}//anti-B0 decays to psi anti-K*0 pi+ pi- 

if ( PcheckDecay(genpart,443, 333, -311)==1){
return 67;}//anti-B0 decays to psi phi anti-K0 

if ( PcheckDecay(genpart,443, 111)==1){
return 68;}//anti-B0 decays to psi pi0 

if ( PcheckDecay(genpart,443, 221)==1){
return 69;}//anti-B0 decays to psi eta 

if ( PcheckDecay(genpart,443, 211, -211)==1){
return 70;}//anti-B0 decays to psi pi+ pi- 

if ( PcheckDecay(genpart,443, 113)==1){
return 71;}//anti-B0 decays to psi rho0 

if ( PcheckDecay(genpart,443, 223)==1){
return 72;}//anti-B0 decays to psi omega 

if ( PcheckDecay(genpart,443, -10311)==1){
return 73;}//anti-B0 decays to psi anti-K_0*0 

if ( PcheckDecay(genpart,443, -315)==1){
return 74;}//anti-B0 decays to psi anti-K_2*0 

if ( PcheckDecay(genpart,443, -30313)==1){
return 75;}//anti-B0 decays to psi anti-K''*0 

if ( PcheckDecay(genpart,443, 3122, -2112)==1){
return 76;}//anti-B0 decays to psi Lambda0 anti-n0 

if ( PcheckDecay(genpart,443, 3222, -2212)==1){
return 77;}//anti-B0 decays to psi Sigma+ anti-p- 

if ( PcheckDecay(genpart,443, -30343)==1){
return 78;}//anti-B0 decays to psi anti-Xsd 

if ( PcheckDecay(genpart,100443, 310)==1){
return 79;}//anti-B0 decays to psi(2S) K_S0 

if ( PcheckDecay(genpart,100443, 130)==1){
return 80;}//anti-B0 decays to psi(2S) K_L0 

if ( PcheckDecay(genpart,100443, -313)==1){
return 81;}//anti-B0 decays to psi(2S) anti-K*S 

if ( PcheckDecay(genpart,100443, -313)==1){
return 82;}//anti-B0 decays to psi(2S) anti-K*L 

if ( PcheckDecay(genpart,100443, -313)==1){
return 83;}//anti-B0 decays to psi(2S) anti-K*0T 

if ( PcheckDecay(genpart,100443, 111)==1){
return 84;}//anti-B0 decays to psi(2S) pi0 

if ( PcheckDecay(genpart,100443, -30343)==1){
return 85;}//anti-B0 decays to psi(2S) anti-Xsd 

if ( PcheckDecay(genpart,10441, 310)==1){
return 86;}//anti-B0 decays to chi_c0 K_S0 

if ( PcheckDecay(genpart,10441, 130)==1){
return 87;}//anti-B0 decays to chi_c0 K_L0 

if ( PcheckDecay(genpart,-313, 10441)==1){
return 88;}//anti-B0 decays to anti-K*S chi_c0 

if ( PcheckDecay(genpart,-313, 10441)==1){
return 89;}//anti-B0 decays to anti-K*L chi_c0 

if ( PcheckDecay(genpart,-313, 10441)==1){
return 90;}//anti-B0 decays to anti-K*0T chi_c0 

if ( PcheckDecay(genpart,10441, -30343)==1){
return 91;}//anti-B0 decays to chi_c0 anti-Xsd 

if ( PcheckDecay(genpart,20443, 310)==1){
return 92;}//anti-B0 decays to chi_c1 K_S0 

if ( PcheckDecay(genpart,20443, 130)==1){
return 93;}//anti-B0 decays to chi_c1 K_L0 

if ( PcheckDecay(genpart,20443, -313)==1){
return 94;}//anti-B0 decays to chi_c1 anti-K*S 

if ( PcheckDecay(genpart,20443, -313)==1){
return 95;}//anti-B0 decays to chi_c1 anti-K*L 

if ( PcheckDecay(genpart,20443, -313)==1){
return 96;}//anti-B0 decays to chi_c1 anti-K*0T 

if ( PcheckDecay(genpart,20443, 111)==1){
return 97;}//anti-B0 decays to chi_c1 pi0 

if ( PcheckDecay(genpart,20443, -321, 211)==1){
return 98;}//anti-B0 decays to chi_c1 K- pi+ 

if ( PcheckDecay(genpart,20443, -311, 111)==1){
return 99;}//anti-B0 decays to chi_c1 anti-K0 pi0 

if ( PcheckDecay(genpart,20443, -311, 211, -211)==1){
return 100;}//anti-B0 decays to chi_c1 anti-K0 pi+ pi- 

if ( PcheckDecay(genpart,20443, -311, 111, 111)==1){
return 101;}//anti-B0 decays to chi_c1 anti-K0 pi0 pi0 

if ( PcheckDecay(genpart,20443, -321, 211, 111)==1){
return 102;}//anti-B0 decays to chi_c1 K- pi+ pi0 

if ( PcheckDecay(genpart,20443, -30343)==1){
return 103;}//anti-B0 decays to chi_c1 anti-Xsd 

if ( PcheckDecay(genpart,445, 310)==1){
return 104;}//anti-B0 decays to chi_c2 K_S0 

if ( PcheckDecay(genpart,445, 130)==1){
return 105;}//anti-B0 decays to chi_c2 K_L0 

if ( PcheckDecay(genpart,445, -313)==1){
return 106;}//anti-B0 decays to chi_c2 anti-K*S 

if ( PcheckDecay(genpart,445, -313)==1){
return 107;}//anti-B0 decays to chi_c2 anti-K*L 

if ( PcheckDecay(genpart,445, -313)==1){
return 108;}//anti-B0 decays to chi_c2 anti-K*0T 

if ( PcheckDecay(genpart,445, -321, 211)==1){
return 109;}//anti-B0 decays to chi_c2 K- pi+ 

if ( PcheckDecay(genpart,445, -311, 111)==1){
return 110;}//anti-B0 decays to chi_c2 anti-K0 pi0 

if ( PcheckDecay(genpart,445, -311, 211, -211)==1){
return 111;}//anti-B0 decays to chi_c2 anti-K0 pi+ pi- 

if ( PcheckDecay(genpart,445, -311, 111, 111)==1){
return 112;}//anti-B0 decays to chi_c2 anti-K0 pi0 pi0 

if ( PcheckDecay(genpart,445, -321, 211, 111)==1){
return 113;}//anti-B0 decays to chi_c2 K- pi+ pi0 

if ( PcheckDecay(genpart,445, -30343)==1){
return 114;}//anti-B0 decays to chi_c2 anti-Xsd 

if ( PcheckDecay(genpart,10443, 310)==1){
return 115;}//anti-B0 decays to h_c K_S0 

if ( PcheckDecay(genpart,10443, 130)==1){
return 116;}//anti-B0 decays to h_c K_L0 

if ( PcheckDecay(genpart,10443, -313)==1){
return 117;}//anti-B0 decays to h_c anti-K*S 

if ( PcheckDecay(genpart,10443, -313)==1){
return 118;}//anti-B0 decays to h_c anti-K*L 

if ( PcheckDecay(genpart,10443, -313)==1){
return 119;}//anti-B0 decays to h_c anti-K*0T 

if ( PcheckDecay(genpart,10443, -30343)==1){
return 120;}//anti-B0 decays to h_c anti-Xsd 

if ( PcheckDecay(genpart,120443, 310)==1){
return 121;}//anti-B0 decays to X(3872) K_S0 

if ( PcheckDecay(genpart,120443, 130)==1){
return 122;}//anti-B0 decays to X(3872) K_L0 

if ( PcheckDecay(genpart,90000443, 310)==1){
return 123;}//anti-B0 decays to Y(3940) K_S0 

if ( PcheckDecay(genpart,90000443, 130)==1){
return 124;}//anti-B0 decays to Y(3940) K_L0 

if ( PcheckDecay(genpart,411, -431)==1){
return 125;}//anti-B0 decays to D+ D_s- 

if ( PcheckDecay(genpart,413, -431)==1){
return 126;}//anti-B0 decays to D*+ D_s- 

if ( PcheckDecay(genpart,-433, 411)==1){
return 127;}//anti-B0 decays to D_s*- D+ 

if ( PcheckDecay(genpart,413, -433)==1){
return 128;}//anti-B0 decays to D*+ D_s*- 

if ( PcheckDecay(genpart,-431, 10411)==1){
return 129;}//anti-B0 decays to D_s- D_0*+ 

if ( PcheckDecay(genpart,-433, 10411)==1){
return 130;}//anti-B0 decays to D_s*- D_0*+ 

if ( PcheckDecay(genpart,20413, -431)==1){
return 131;}//anti-B0 decays to D'_1+ D_s- 

if ( PcheckDecay(genpart,20413, -433)==1){
return 132;}//anti-B0 decays to D'_1+ D_s*- 

if ( PcheckDecay(genpart,10413, -431)==1){
return 133;}//anti-B0 decays to D_1+ D_s- 

if ( PcheckDecay(genpart,10413, -433)==1){
return 134;}//anti-B0 decays to D_1+ D_s*- 

if ( PcheckDecay(genpart,415, -431)==1){
return 135;}//anti-B0 decays to D_2*+ D_s- 

if ( PcheckDecay(genpart,415, -433)==1){
return 136;}//anti-B0 decays to D_2*+ D_s*- 

if ( PcheckDecay(genpart,411, -10431)==1){
return 137;}//anti-B0 decays to D+ D_s0*- 

if ( PcheckDecay(genpart,413, -10431)==1){
return 138;}//anti-B0 decays to D*+ D_s0*- 

if ( PcheckDecay(genpart,-20433, 411)==1){
return 139;}//anti-B0 decays to D'_s1- D+ 

if ( PcheckDecay(genpart,-20433, 413)==1){
return 140;}//anti-B0 decays to D'_s1- D*+ 

if ( PcheckDecay(genpart,-431, 411, 111)==1){
return 141;}//anti-B0 decays to D_s- D+ pi0 

if ( PcheckDecay(genpart,-431, 421, 211)==1){
return 142;}//anti-B0 decays to D_s- D0 pi+ 

if ( PcheckDecay(genpart,-433, 411, 111)==1){
return 143;}//anti-B0 decays to D_s*- D+ pi0 

if ( PcheckDecay(genpart,-433, 421, 211)==1){
return 144;}//anti-B0 decays to D_s*- D0 pi+ 

if ( PcheckDecay(genpart,-431, 411, 211, -211)==1){
return 145;}//anti-B0 decays to D_s- D+ pi+ pi- 

if ( PcheckDecay(genpart,-431, 411, 111, 111)==1){
return 146;}//anti-B0 decays to D_s- D+ pi0 pi0 

if ( PcheckDecay(genpart,-431, 421, 211, 111)==1){
return 147;}//anti-B0 decays to D_s- D0 pi+ pi0 

if ( PcheckDecay(genpart,-433, 411, 211, -211)==1){
return 148;}//anti-B0 decays to D_s*- D+ pi+ pi- 

if ( PcheckDecay(genpart,-433, 411, 111, 111)==1){
return 149;}//anti-B0 decays to D_s*- D+ pi0 pi0 

if ( PcheckDecay(genpart,-433, 421, 211, 111)==1){
return 150;}//anti-B0 decays to D_s*- D0 pi+ pi0 

if ( PcheckDecay(genpart,431, -211, -311)==1){
return 151;}//anti-B0 decays to D_s+ pi- anti-K0 

if ( PcheckDecay(genpart,433, -211, -311)==1){
return 152;}//anti-B0 decays to D_s*+ pi- anti-K0 

if ( PcheckDecay(genpart,431, -211, -311, 111)==1){
return 153;}//anti-B0 decays to D_s+ pi- anti-K0 pi0 

if ( PcheckDecay(genpart,433, -211, -311, 111)==1){
return 154;}//anti-B0 decays to D_s*+ pi- anti-K0 pi0 

if ( PcheckDecay(genpart,431, -211, -321, 211)==1){
return 155;}//anti-B0 decays to D_s+ pi- K- pi+ 

if ( PcheckDecay(genpart,433, -211, -321, 211)==1){
return 156;}//anti-B0 decays to D_s*+ pi- K- pi+ 

if ( PcheckDecay(genpart,431, -211, -311, 111, 111)==1){
return 157;}//anti-B0 decays to D_s+ pi- anti-K0 pi0 pi0 

if ( PcheckDecay(genpart,433, -211, -311, 111, 111)==1){
return 158;}//anti-B0 decays to D_s*+ pi- anti-K0 pi0 pi0 

if ( PcheckDecay(genpart,431, -211, -311, 211, -211)==1){
return 159;}//anti-B0 decays to D_s+ pi- anti-K0 pi+ pi- 

if ( PcheckDecay(genpart,433, -211, -311, 211, -211)==1){
return 160;}//anti-B0 decays to D_s*+ pi- anti-K0 pi+ pi- 

if ( PcheckDecay(genpart,431, -211, -321, 211, 111)==1){
return 161;}//anti-B0 decays to D_s+ pi- K- pi+ pi0 

if ( PcheckDecay(genpart,433, -211, -321, 211, 111)==1){
return 162;}//anti-B0 decays to D_s*+ pi- K- pi+ pi0 

if ( PcheckDecay(genpart,-411, 411)==1){
return 163;}//anti-B0 decays to D- D+ 

if ( PcheckDecay(genpart,-413, 411)==1){
return 164;}//anti-B0 decays to D*- D+ 

if ( PcheckDecay(genpart,413, -411)==1){
return 165;}//anti-B0 decays to D*+ D- 

if ( PcheckDecay(genpart,413, -413)==1){
return 166;}//anti-B0 decays to D*+ D*- 

if ( PcheckDecay(genpart,-20413, 411)==1){
return 167;}//anti-B0 decays to D'_1- D+ 

if ( PcheckDecay(genpart,20413, -411)==1){
return 168;}//anti-B0 decays to D'_1+ D- 

if ( PcheckDecay(genpart,-20413, 413)==1){
return 169;}//anti-B0 decays to D'_1- D*+ 

if ( PcheckDecay(genpart,20413, -413)==1){
return 170;}//anti-B0 decays to D'_1+ D*- 

if ( PcheckDecay(genpart,-10413, 411)==1){
return 171;}//anti-B0 decays to D_1- D+ 

if ( PcheckDecay(genpart,10413, -411)==1){
return 172;}//anti-B0 decays to D_1+ D- 

if ( PcheckDecay(genpart,-10413, 413)==1){
return 173;}//anti-B0 decays to D_1- D*+ 

if ( PcheckDecay(genpart,10413, -413)==1){
return 174;}//anti-B0 decays to D_1+ D*- 

if ( PcheckDecay(genpart,-415, 411)==1){
return 175;}//anti-B0 decays to D_2*- D+ 

if ( PcheckDecay(genpart,415, -411)==1){
return 176;}//anti-B0 decays to D_2*+ D- 

if ( PcheckDecay(genpart,-415, 413)==1){
return 177;}//anti-B0 decays to D_2*- D*+ 

if ( PcheckDecay(genpart,415, -413)==1){
return 178;}//anti-B0 decays to D_2*+ D*- 

if ( PcheckDecay(genpart,-10433, 411)==1){
return 179;}//anti-B0 decays to D_s1- D+ 

if ( PcheckDecay(genpart,-10433, 413)==1){
return 180;}//anti-B0 decays to D_s1- D*+ 

if ( PcheckDecay(genpart,-435, 411)==1){
return 181;}//anti-B0 decays to D_s2*- D+ 

if ( PcheckDecay(genpart,-435, 413)==1){
return 182;}//anti-B0 decays to D_s2*- D*+ 

if ( PcheckDecay(genpart,-9000433, 411)==1){
return 183;}//anti-B0 decays to D_sj(2700)- D+ 

if ( PcheckDecay(genpart,-9000433, 413)==1){
return 184;}//anti-B0 decays to D_sj(2700)- D*+ 

if ( PcheckDecay(genpart,30443, 310)==1){
return 185;}//anti-B0 decays to psi(3770) K_S0 

if ( PcheckDecay(genpart,30443, 130)==1){
return 186;}//anti-B0 decays to psi(3770) K_L0 

if ( PcheckDecay(genpart,30443, -313)==1){
return 187;}//anti-B0 decays to psi(3770) anti-K*S 

if ( PcheckDecay(genpart,30443, -313)==1){
return 188;}//anti-B0 decays to psi(3770) anti-K*L 

if ( PcheckDecay(genpart,30443, -313)==1){
return 189;}//anti-B0 decays to psi(3770) anti-K*0T 

if ( PcheckDecay(genpart,30443, -321, 211)==1){
return 190;}//anti-B0 decays to psi(3770) K- pi+ 

if ( PcheckDecay(genpart,30443, -311, 111)==1){
return 191;}//anti-B0 decays to psi(3770) anti-K0 pi0 

if ( PcheckDecay(genpart,30443, -311, 211, -211)==1){
return 192;}//anti-B0 decays to psi(3770) anti-K0 pi+ pi- 

if ( PcheckDecay(genpart,30443, -311, 111, 111)==1){
return 193;}//anti-B0 decays to psi(3770) anti-K0 pi0 pi0 

if ( PcheckDecay(genpart,30443, -321, 211, 111)==1){
return 194;}//anti-B0 decays to psi(3770) K- pi+ pi0 

if ( PcheckDecay(genpart,9000443, 310)==1){
return 195;}//anti-B0 decays to psi(4040) K_S0 

if ( PcheckDecay(genpart,9000443, 130)==1){
return 196;}//anti-B0 decays to psi(4040) K_L0 

if ( PcheckDecay(genpart,9000443, -313)==1){
return 197;}//anti-B0 decays to psi(4040) anti-K*S 

if ( PcheckDecay(genpart,9000443, -313)==1){
return 198;}//anti-B0 decays to psi(4040) anti-K*L 

if ( PcheckDecay(genpart,9000443, -313)==1){
return 199;}//anti-B0 decays to psi(4040) anti-K*0T 

if ( PcheckDecay(genpart,9000443, -321, 211)==1){
return 200;}//anti-B0 decays to psi(4040) K- pi+ 

if ( PcheckDecay(genpart,9000443, -311, 111)==1){
return 201;}//anti-B0 decays to psi(4040) anti-K0 pi0 

if ( PcheckDecay(genpart,9010443, 310)==1){
return 202;}//anti-B0 decays to psi(4160) K_S0 

if ( PcheckDecay(genpart,9010443, 130)==1){
return 203;}//anti-B0 decays to psi(4160) K_L0 

if ( PcheckDecay(genpart,9010443, -313)==1){
return 204;}//anti-B0 decays to psi(4160) anti-K*S 

if ( PcheckDecay(genpart,9010443, -313)==1){
return 205;}//anti-B0 decays to psi(4160) anti-K*L 

if ( PcheckDecay(genpart,9010443, -313)==1){
return 206;}//anti-B0 decays to psi(4160) anti-K*0T 

if ( PcheckDecay(genpart,9010443, -321, 211)==1){
return 207;}//anti-B0 decays to psi(4160) K- pi+ 

if ( PcheckDecay(genpart,9010443, -311, 111)==1){
return 208;}//anti-B0 decays to psi(4160) anti-K0 pi0 

//if ( PcheckDecay(genpart,return 209;}//anti-B0 decays to 

if ( PcheckDecay(genpart,-411, 411, 111)==1){
return 210;}//anti-B0 decays to D- D+ pi0 

if ( PcheckDecay(genpart,-411, 421, 211)==1){
return 211;}//anti-B0 decays to D- D0 pi+ 

if ( PcheckDecay(genpart,-421, 411, -211)==1){
return 212;}//anti-B0 decays to anti-D0 D+ pi- 

if ( PcheckDecay(genpart,-413, 411, 111)==1){
return 213;}//anti-B0 decays to D*- D+ pi0 

if ( PcheckDecay(genpart,-413, 421, 211)==1){
return 214;}//anti-B0 decays to D*- D0 pi+ 

if ( PcheckDecay(genpart,-423, 411, -211)==1){
return 215;}//anti-B0 decays to anti-D*0 D+ pi- 

if ( PcheckDecay(genpart,-411, 413, 111)==1){
return 216;}//anti-B0 decays to D- D*+ pi0 

if ( PcheckDecay(genpart,-411, 423, 211)==1){
return 217;}//anti-B0 decays to D- D*0 pi+ 

if ( PcheckDecay(genpart,-421, 413, -211)==1){
return 218;}//anti-B0 decays to anti-D0 D*+ pi- 

if ( PcheckDecay(genpart,-413, 413, 111)==1){
return 219;}//anti-B0 decays to D*- D*+ pi0 

if ( PcheckDecay(genpart,-413, 423, 211)==1){
return 220;}//anti-B0 decays to D*- D*0 pi+ 

if ( PcheckDecay(genpart,-423, 413, -211)==1){
return 221;}//anti-B0 decays to anti-D*0 D*+ pi- 

if ( PcheckDecay(genpart,-411, 411, 113)==1){
return 222;}//anti-B0 decays to D- D+ rho0 

if ( PcheckDecay(genpart,-411, 421, 213)==1){
return 223;}//anti-B0 decays to D- D0 rho+ 

if ( PcheckDecay(genpart,-421, 411, -213)==1){
return 224;}//anti-B0 decays to anti-D0 D+ rho- 

if ( PcheckDecay(genpart,-413, 411, 113)==1){
return 225;}//anti-B0 decays to D*- D+ rho0 

if ( PcheckDecay(genpart,-413, 421, 213)==1){
return 226;}//anti-B0 decays to D*- D0 rho+ 

if ( PcheckDecay(genpart,-423, 411, -213)==1){
return 227;}//anti-B0 decays to anti-D*0 D+ rho- 

if ( PcheckDecay(genpart,-411, 413, 113)==1){
return 228;}//anti-B0 decays to D- D*+ rho0 

if ( PcheckDecay(genpart,-411, 423, 213)==1){
return 229;}//anti-B0 decays to D- D*0 rho+ 

if ( PcheckDecay(genpart,-421, 413, -213)==1){
return 230;}//anti-B0 decays to anti-D0 D*+ rho- 

if ( PcheckDecay(genpart,-413, 413, 113)==1){
return 231;}//anti-B0 decays to D*- D*+ rho0 

if ( PcheckDecay(genpart,-413, 423, 213)==1){
return 232;}//anti-B0 decays to D*- D*0 rho+ 

if ( PcheckDecay(genpart,-423, 413, -213)==1){
return 233;}//anti-B0 decays to anti-D*0 D*+ rho- 

if ( PcheckDecay(genpart,-411, 411, 111, 111)==1){
return 234;}//anti-B0 decays to D- D+ pi0 pi0 

if ( PcheckDecay(genpart,-411, 411, 211, -211)==1){
return 235;}//anti-B0 decays to D- D+ pi+ pi- 

if ( PcheckDecay(genpart,-411, 421, 211, 111)==1){
return 236;}//anti-B0 decays to D- D0 pi+ pi0 

if ( PcheckDecay(genpart,-421, 411, -211, 111)==1){
return 237;}//anti-B0 decays to anti-D0 D+ pi- pi0 

if ( PcheckDecay(genpart,-421, 421, 211, -211)==1){
return 238;}//anti-B0 decays to anti-D0 D0 pi+ pi- 

if ( PcheckDecay(genpart,-413, 411, 111, 111)==1){
return 239;}//anti-B0 decays to D*- D+ pi0 pi0 

if ( PcheckDecay(genpart,-413, 411, 211, -211)==1){
return 240;}//anti-B0 decays to D*- D+ pi+ pi- 

if ( PcheckDecay(genpart,-413, 421, 211, 111)==1){
return 241;}//anti-B0 decays to D*- D0 pi+ pi0 

if ( PcheckDecay(genpart,-423, 411, -211, 111)==1){
return 242;}//anti-B0 decays to anti-D*0 D+ pi- pi0 

if ( PcheckDecay(genpart,-423, 421, 211, -211)==1){
return 243;}//anti-B0 decays to anti-D*0 D0 pi+ pi- 

if ( PcheckDecay(genpart,-411, 413, 111, 111)==1){
return 244;}//anti-B0 decays to D- D*+ pi0 pi0 

if ( PcheckDecay(genpart,-411, 413, 211, -211)==1){
return 245;}//anti-B0 decays to D- D*+ pi+ pi- 

if ( PcheckDecay(genpart,-411, 423, 211, 111)==1){
return 246;}//anti-B0 decays to D- D*0 pi+ pi0 

if ( PcheckDecay(genpart,-421, 413, -211, 111)==1){
return 247;}//anti-B0 decays to anti-D0 D*+ pi- pi0 

if ( PcheckDecay(genpart,-421, 423, 211, -211)==1){
return 248;}//anti-B0 decays to anti-D0 D*0 pi+ pi- 

if ( PcheckDecay(genpart,-413, 413, 111, 111)==1){
return 249;}//anti-B0 decays to D*- D*+ pi0 pi0 

if ( PcheckDecay(genpart,-413, 413, 211, -211)==1){
return 250;}//anti-B0 decays to D*- D*+ pi+ pi- 

if ( PcheckDecay(genpart,-413, 423, 211, 111)==1){
return 251;}//anti-B0 decays to D*- D*0 pi+ pi0 

if ( PcheckDecay(genpart,-423, 413, -211, 111)==1){
return 252;}//anti-B0 decays to anti-D*0 D*+ pi- pi0 

if ( PcheckDecay(genpart,-423, 423, 211, -211)==1){
return 253;}//anti-B0 decays to anti-D*0 D*0 pi+ pi- 

if ( PcheckDecay(genpart,411, -421, -321)==1){
return 254;}//anti-B0 decays to D+ anti-D0 K- 

if ( PcheckDecay(genpart,411, -411, -311)==1){
return 255;}//anti-B0 decays to D+ D- anti-K0 

if ( PcheckDecay(genpart,411, -421, -323)==1){
return 256;}//anti-B0 decays to D+ anti-D0 K*- 

if ( PcheckDecay(genpart,411, -411, -313)==1){
return 257;}//anti-B0 decays to D+ D- anti-K*0 

if ( PcheckDecay(genpart,413, -421, -321)==1){
return 258;}//anti-B0 decays to D*+ anti-D0 K- 

if ( PcheckDecay(genpart,413, -411, -311)==1){
return 259;}//anti-B0 decays to D*+ D- anti-K0 

if ( PcheckDecay(genpart,411, -423, -321)==1){
return 260;}//anti-B0 decays to D+ anti-D*0 K- 

if ( PcheckDecay(genpart,413, -421, -323)==1){
return 261;}//anti-B0 decays to D*+ anti-D0 K*- 

if ( PcheckDecay(genpart,413, -411, -313)==1){
return 262;}//anti-B0 decays to D*+ D- anti-K*0 

if ( PcheckDecay(genpart,411, -423, -323)==1){
return 263;}//anti-B0 decays to D+ anti-D*0 K*- 

if ( PcheckDecay(genpart,411, -413, -313)==1){
return 264;}//anti-B0 decays to D+ D*- anti-K*0 

if ( PcheckDecay(genpart,413, -423, -321)==1){
return 265;}//anti-B0 decays to D*+ anti-D*0 K- 

if ( PcheckDecay(genpart,413, -413, -311)==1){
return 266;}//anti-B0 decays to D*+ D*- anti-K0 

if ( PcheckDecay(genpart,413, -423, -323)==1){
return 267;}//anti-B0 decays to D*+ anti-D*0 K*- 

if ( PcheckDecay(genpart,413, -413, -313)==1){
return 268;}//anti-B0 decays to D*+ D*- anti-K*0 

if ( PcheckDecay(genpart,10413, -421, -321)==1){
return 269;}//anti-B0 decays to D_1+ anti-D0 K- 

if ( PcheckDecay(genpart,10413, -411, -311)==1){
return 270;}//anti-B0 decays to D_1+ D- anti-K0 

if ( PcheckDecay(genpart,411, -10423, -321)==1){
return 271;}//anti-B0 decays to D+ anti-D_10 K- 

if ( PcheckDecay(genpart,411, -10413, -311)==1){
return 272;}//anti-B0 decays to D+ D_1- anti-K0 

if ( PcheckDecay(genpart,10413, -423, -321)==1){
return 273;}//anti-B0 decays to D_1+ anti-D*0 K- 

if ( PcheckDecay(genpart,10413, -413, -311)==1){
return 274;}//anti-B0 decays to D_1+ D*- anti-K0 

if ( PcheckDecay(genpart,413, -10423, -321)==1){
return 275;}//anti-B0 decays to D*+ anti-D_10 K- 

if ( PcheckDecay(genpart,413, -10413, -311)==1){
return 276;}//anti-B0 decays to D*+ D_1- anti-K0 

if ( PcheckDecay(genpart,415, -421, -321)==1){
return 277;}//anti-B0 decays to D_2*+ anti-D0 K- 

if ( PcheckDecay(genpart,415, -411, -311)==1){
return 278;}//anti-B0 decays to D_2*+ D- anti-K0 

if ( PcheckDecay(genpart,411, -425, -321)==1){
return 279;}//anti-B0 decays to D+ anti-D_2*0 K- 

if ( PcheckDecay(genpart,411, -415, -311)==1){
return 280;}//anti-B0 decays to D+ D_2*- anti-K0 

if ( PcheckDecay(genpart,415, -423, -321)==1){
return 281;}//anti-B0 decays to D_2*+ anti-D*0 K- 

if ( PcheckDecay(genpart,415, -413, -311)==1){
return 282;}//anti-B0 decays to D_2*+ D*- anti-K0 

if ( PcheckDecay(genpart,413, -425, -321)==1){
return 283;}//anti-B0 decays to D*+ anti-D_2*0 K- 

if ( PcheckDecay(genpart,413, -415, -311)==1){
return 284;}//anti-B0 decays to D*+ D_2*- anti-K0 

if ( PcheckDecay(genpart,411, -421, -321, 111)==1){
return 285;}//anti-B0 decays to D+ anti-D0 K- pi0 

if ( PcheckDecay(genpart,421, -421, -321, 211)==1){
return 286;}//anti-B0 decays to D0 anti-D0 K- pi+ 

if ( PcheckDecay(genpart,411, -421, -311, -211)==1){
return 287;}//anti-B0 decays to D+ anti-D0 anti-K0 pi- 

if ( PcheckDecay(genpart,411, -411, -321, 211)==1){
return 288;}//anti-B0 decays to D+ D- K- pi+ 

if ( PcheckDecay(genpart,411, -411, -311, 111)==1){
return 289;}//anti-B0 decays to D+ D- anti-K0 pi0 

if ( PcheckDecay(genpart,421, -411, -311, 211)==1){
return 290;}//anti-B0 decays to D0 D- anti-K0 pi+ 

if ( PcheckDecay(genpart,413, -421, -321, 111)==1){
return 291;}//anti-B0 decays to D*+ anti-D0 K- pi0 

if ( PcheckDecay(genpart,423, -421, -321, 211)==1){
return 292;}//anti-B0 decays to D*0 anti-D0 K- pi+ 

if ( PcheckDecay(genpart,413, -421, -311, -211)==1){
return 293;}//anti-B0 decays to D*+ anti-D0 anti-K0 pi- 

if ( PcheckDecay(genpart,413, -411, -321, 211)==1){
return 294;}//anti-B0 decays to D*+ D- K- pi+ 

if ( PcheckDecay(genpart,413, -411, -311, 111)==1){
return 295;}//anti-B0 decays to D*+ D- anti-K0 pi0 

if ( PcheckDecay(genpart,423, -411, -311, 211)==1){
return 296;}//anti-B0 decays to D*0 D- anti-K0 pi+ 

if ( PcheckDecay(genpart,411, -423, -321, 111)==1){
return 297;}//anti-B0 decays to D+ anti-D*0 K- pi0 

if ( PcheckDecay(genpart,421, -423, -321, 211)==1){
return 298;}//anti-B0 decays to D0 anti-D*0 K- pi+ 

if ( PcheckDecay(genpart,411, -423, -311, -211)==1){
return 299;}//anti-B0 decays to D+ anti-D*0 anti-K0 pi- 

if ( PcheckDecay(genpart,411, -413, -321, 211)==1){
return 300;}//anti-B0 decays to D+ D*- K- pi+ 

if ( PcheckDecay(genpart,411, -413, -311, 111)==1){
return 301;}//anti-B0 decays to D+ D*- anti-K0 pi0 

if ( PcheckDecay(genpart,421, -413, -311, 211)==1){
return 302;}//anti-B0 decays to D0 D*- anti-K0 pi+ 

if ( PcheckDecay(genpart,413, -423, -321, 111)==1){
return 303;}//anti-B0 decays to D*+ anti-D*0 K- pi0 

if ( PcheckDecay(genpart,423, -423, -321, 211)==1){
return 304;}//anti-B0 decays to D*0 anti-D*0 K- pi+ 

if ( PcheckDecay(genpart,413, -423, -311, -211)==1){
return 305;}//anti-B0 decays to D*+ anti-D*0 anti-K0 pi- 

if ( PcheckDecay(genpart,413, -413, -321, 211)==1){
return 306;}//anti-B0 decays to D*+ D*- K- pi+ 

if ( PcheckDecay(genpart,413, -413, -311, 111)==1){
return 307;}//anti-B0 decays to D*+ D*- anti-K0 pi0 

if ( PcheckDecay(genpart,423, -413, -311, 211)==1){
return 308;}//anti-B0 decays to D*0 D*- anti-K0 pi+ 

if ( PcheckDecay(genpart,411, -421, -323, 111)==1){
return 309;}//anti-B0 decays to D+ anti-D0 K*- pi0 

if ( PcheckDecay(genpart,421, -421, -323, 211)==1){
return 310;}//anti-B0 decays to D0 anti-D0 K*- pi+ 

if ( PcheckDecay(genpart,411, -421, -313, -211)==1){
return 311;}//anti-B0 decays to D+ anti-D0 anti-K*0 pi- 

if ( PcheckDecay(genpart,411, -411, -323, 211)==1){
return 312;}//anti-B0 decays to D+ D- K*- pi+ 

if ( PcheckDecay(genpart,411, -411, -313, 111)==1){
return 313;}//anti-B0 decays to D+ D- anti-K*0 pi0 

if ( PcheckDecay(genpart,421, -411, -313, 211)==1){
return 314;}//anti-B0 decays to D0 D- anti-K*0 pi+ 

if ( PcheckDecay(genpart,411, -421, -321, 113)==1){
return 315;}//anti-B0 decays to D+ anti-D0 K- rho0 

if ( PcheckDecay(genpart,421, -421, -321, 213)==1){
return 316;}//anti-B0 decays to D0 anti-D0 K- rho+ 

if ( PcheckDecay(genpart,411, -421, -311, -213)==1){
return 317;}//anti-B0 decays to D+ anti-D0 anti-K0 rho- 

if ( PcheckDecay(genpart,411, -411, -321, 213)==1){
return 318;}//anti-B0 decays to D+ D- K- rho+ 

if ( PcheckDecay(genpart,411, -411, -311, 113)==1){
return 319;}//anti-B0 decays to D+ D- anti-K0 rho0 

if ( PcheckDecay(genpart,421, -411, -311, 213)==1){
return 320;}//anti-B0 decays to D0 D- anti-K0 rho+ 

if ( PcheckDecay(genpart,421, -421, -311)==1){
return 321;}//anti-B0 decays to D0 anti-D0 anti-K0 

if ( PcheckDecay(genpart,421, -421, -313)==1){
return 322;}//anti-B0 decays to D0 anti-D0 anti-K*0 

if ( PcheckDecay(genpart,421, -423, -313)==1){
return 323;}//anti-B0 decays to D0 anti-D*0 anti-K*0 

if ( PcheckDecay(genpart,423, -421, -311)==1){
return 324;}//anti-B0 decays to D*0 anti-D0 anti-K0 

if ( PcheckDecay(genpart,421, -423, -311)==1){
return 325;}//anti-B0 decays to D0 anti-D*0 anti-K0 

if ( PcheckDecay(genpart,423, -421, -313)==1){
return 326;}//anti-B0 decays to D*0 anti-D0 anti-K*0 

if ( PcheckDecay(genpart,423, -423, -311)==1){
return 327;}//anti-B0 decays to D*0 anti-D*0 anti-K0 

if ( PcheckDecay(genpart,423, -423, -313)==1){
return 328;}//anti-B0 decays to D*0 anti-D*0 anti-K*0 

if ( PcheckDecay(genpart,10423, -421, -311)==1){
return 329;}//anti-B0 decays to D_10 anti-D0 anti-K0 

if ( PcheckDecay(genpart,421, -10423, -311)==1){
return 330;}//anti-B0 decays to D0 anti-D_10 anti-K0 

if ( PcheckDecay(genpart,10423, -423, -311)==1){
return 331;}//anti-B0 decays to D_10 anti-D*0 anti-K0 

if ( PcheckDecay(genpart,423, -10423, -311)==1){
return 332;}//anti-B0 decays to D*0 anti-D_10 anti-K0 

if ( PcheckDecay(genpart,425, -421, -311)==1){
return 333;}//anti-B0 decays to D_2*0 anti-D0 anti-K0 

if ( PcheckDecay(genpart,421, -425, -311)==1){
return 334;}//anti-B0 decays to D0 anti-D_2*0 anti-K0 

if ( PcheckDecay(genpart,425, -423, -311)==1){
return 335;}//anti-B0 decays to D_2*0 anti-D*0 anti-K0 

if ( PcheckDecay(genpart,423, -425, -311)==1){
return 336;}//anti-B0 decays to D*0 anti-D_2*0 anti-K0 

if ( PcheckDecay(genpart,411, -211)==1){
return 337;}//anti-B0 decays to D+ pi- 

if ( PcheckDecay(genpart,411, -211, 111)==1){
return 338;}//anti-B0 decays to D+ pi- pi0 

if ( PcheckDecay(genpart,-213, 411)==1){
return 339;}//anti-B0 decays to rho- D+ 

if ( PcheckDecay(genpart,411, 211, -211, -211)==1){
return 340;}//anti-B0 decays to D+ pi+ pi- pi- 

if ( PcheckDecay(genpart,411, 111, -211, 111)==1){
return 341;}//anti-B0 decays to D+ pi0 pi- pi0 

if ( PcheckDecay(genpart,411, 113, -211)==1){
return 342;}//anti-B0 decays to D+ rho0 pi- 

if ( PcheckDecay(genpart,411, -213, 111)==1){
return 343;}//anti-B0 decays to D+ rho- pi0 

if ( PcheckDecay(genpart,-20213, 411)==1){
return 344;}//anti-B0 decays to a_1- D+ 

if ( PcheckDecay(genpart,-10213, 411)==1){
return 345;}//anti-B0 decays to b_1- D+ 

if ( PcheckDecay(genpart,411, 211, -211, -211, 111)==1){
return 346;}//anti-B0 decays to D+ pi+ pi- pi- pi0 

if ( PcheckDecay(genpart,411, 223, -211)==1){
return 347;}//anti-B0 decays to D+ omega pi- 

if ( PcheckDecay(genpart,411, 2212, -2212, -211)==1){
return 348;}//anti-B0 decays to D+ p+ anti-p- pi- 

if ( PcheckDecay(genpart,411, -2212, 2112)==1){
return 349;}//anti-B0 decays to D+ anti-p- n0 

if ( PcheckDecay(genpart,413, -211)==1){
return 350;}//anti-B0 decays to D*+ pi- 

if ( PcheckDecay(genpart,413, -211, 111)==1){
return 351;}//anti-B0 decays to D*+ pi- pi0 

if ( PcheckDecay(genpart,-213, 413)==1){
return 352;}//anti-B0 decays to rho- D*+ 

if ( PcheckDecay(genpart,413, 211, -211, -211)==1){
return 353;}//anti-B0 decays to D*+ pi+ pi- pi- 

if ( PcheckDecay(genpart,413, 111, -211, 111)==1){
return 354;}//anti-B0 decays to D*+ pi0 pi- pi0 

if ( PcheckDecay(genpart,413, 113, -211)==1){
return 355;}//anti-B0 decays to D*+ rho0 pi- 

if ( PcheckDecay(genpart,413, -213, 111)==1){
return 356;}//anti-B0 decays to D*+ rho- pi0 

if ( PcheckDecay(genpart,413, -20213)==1){
return 357;}//anti-B0 decays to D*+ a_1- 

if ( PcheckDecay(genpart,413, -10213)==1){
return 358;}//anti-B0 decays to D*+ b_1- 

if ( PcheckDecay(genpart,413, 211, -211, -211, 111)==1){
return 359;}//anti-B0 decays to D*+ pi+ pi- pi- pi0 

if ( PcheckDecay(genpart,413, 223, -211)==1){
return 360;}//anti-B0 decays to D*+ omega pi- 

if ( PcheckDecay(genpart,413, 2212, -2212, -211)==1){
return 361;}//anti-B0 decays to D*+ p+ anti-p- pi- 

if ( PcheckDecay(genpart,413, -2212, 2112)==1){
return 362;}//anti-B0 decays to D*+ anti-p- n0 

if ( PcheckDecay(genpart,10413, -211)==1){
return 363;}//anti-B0 decays to D_1+ pi- 

if ( PcheckDecay(genpart,415, -211)==1){
return 364;}//anti-B0 decays to D_2*+ pi- 

if ( PcheckDecay(genpart,10413, -213)==1){
return 365;}//anti-B0 decays to D_1+ rho- 

if ( PcheckDecay(genpart,415, -213)==1){
return 366;}//anti-B0 decays to D_2*+ rho- 

if ( PcheckDecay(genpart,425, 111)==1){
return 367;}//anti-B0 decays to D_2*0 pi0 

if ( PcheckDecay(genpart,10423, 111)==1){
return 368;}//anti-B0 decays to D_10 pi0 

if ( PcheckDecay(genpart,20423, 111)==1){
return 369;}//anti-B0 decays to D'_10 pi0 

if ( PcheckDecay(genpart,10421, 111)==1){
return 370;}//anti-B0 decays to D_0*0 pi0 

if ( PcheckDecay(genpart,425, 221)==1){
return 371;}//anti-B0 decays to D_2*0 eta 

if ( PcheckDecay(genpart,10423, 221)==1){
return 372;}//anti-B0 decays to D_10 eta 

if ( PcheckDecay(genpart,20423, 221)==1){
return 373;}//anti-B0 decays to D'_10 eta 

if ( PcheckDecay(genpart,10421, 221)==1){
return 374;}//anti-B0 decays to D_0*0 eta 

if ( PcheckDecay(genpart,425, 331)==1){
return 375;}//anti-B0 decays to D_2*0 eta' 

if ( PcheckDecay(genpart,10423, 331)==1){
return 376;}//anti-B0 decays to D_10 eta' 

if ( PcheckDecay(genpart,20423, 331)==1){
return 377;}//anti-B0 decays to D'_10 eta' 

if ( PcheckDecay(genpart,10421, 331)==1){
return 378;}//anti-B0 decays to D_0*0 eta' 

if ( PcheckDecay(genpart,425, 223)==1){
return 379;}//anti-B0 decays to D_2*0 omega 

if ( PcheckDecay(genpart,10423, 223)==1){
return 380;}//anti-B0 decays to D_10 omega 

if ( PcheckDecay(genpart,20423, 223)==1){
return 381;}//anti-B0 decays to D'_10 omega 

if ( PcheckDecay(genpart,10421, 223)==1){
return 382;}//anti-B0 decays to D_0*0 omega 

if ( PcheckDecay(genpart,425, 113)==1){
return 383;}//anti-B0 decays to D_2*0 rho0 

if ( PcheckDecay(genpart,10423, 113)==1){
return 384;}//anti-B0 decays to D_10 rho0 

if ( PcheckDecay(genpart,20423, 113)==1){
return 385;}//anti-B0 decays to D'_10 rho0 

if ( PcheckDecay(genpart,10421, 113)==1){
return 386;}//anti-B0 decays to D_0*0 rho0 

//if ( PcheckDecay(genpart,return 387;}//anti-B0 decays to 

if ( PcheckDecay(genpart,421, 111)==1){
return 388;}//anti-B0 decays to D0 pi0 

if ( PcheckDecay(genpart,421, 211, -211)==1){
return 389;}//anti-B0 decays to D0 pi+ pi- 

if ( PcheckDecay(genpart,113, 421)==1){
return 390;}//anti-B0 decays to rho0 D0 

if ( PcheckDecay(genpart,421, 111, 111)==1){
return 391;}//anti-B0 decays to D0 pi0 pi0 

if ( PcheckDecay(genpart,421, 221)==1){
return 392;}//anti-B0 decays to D0 eta 

if ( PcheckDecay(genpart,421, 331)==1){
return 393;}//anti-B0 decays to D0 eta' 

if ( PcheckDecay(genpart,223, 421)==1){
return 394;}//anti-B0 decays to omega D0 

if ( PcheckDecay(genpart,421, 211, -211, 111)==1){
return 395;}//anti-B0 decays to D0 pi+ pi- pi0 

if ( PcheckDecay(genpart,421, 111, 111, 111)==1){
return 396;}//anti-B0 decays to D0 pi0 pi0 pi0 

if ( PcheckDecay(genpart,421, 2212, -2212)==1){
return 397;}//anti-B0 decays to D0 p+ anti-p- 

if ( PcheckDecay(genpart,423, 111)==1){
return 398;}//anti-B0 decays to D*0 pi0 

if ( PcheckDecay(genpart,423, 211, -211)==1){
return 399;}//anti-B0 decays to D*0 pi+ pi- 

if ( PcheckDecay(genpart,423, 113)==1){
return 400;}//anti-B0 decays to D*0 rho0 

if ( PcheckDecay(genpart,423, 111, 111)==1){
return 401;}//anti-B0 decays to D*0 pi0 pi0 

if ( PcheckDecay(genpart,423, 221)==1){
return 402;}//anti-B0 decays to D*0 eta 

if ( PcheckDecay(genpart,423, 331)==1){
return 403;}//anti-B0 decays to D*0 eta' 

if ( PcheckDecay(genpart,423, 223)==1){
return 404;}//anti-B0 decays to D*0 omega 

if ( PcheckDecay(genpart,423, 211, -211, 111)==1){
return 405;}//anti-B0 decays to D*0 pi+ pi- pi0 

if ( PcheckDecay(genpart,423, 111, 111, 111)==1){
return 406;}//anti-B0 decays to D*0 pi0 pi0 pi0 

if ( PcheckDecay(genpart,423, 2212, -2212)==1){
return 407;}//anti-B0 decays to D*0 p+ anti-p- 

if ( PcheckDecay(genpart,411, -321)==1){
return 408;}//anti-B0 decays to D+ K- 

if ( PcheckDecay(genpart,-323, 411)==1){
return 409;}//anti-B0 decays to K*- D+ 

if ( PcheckDecay(genpart,413, -321)==1){
return 410;}//anti-B0 decays to D*+ K- 

if ( PcheckDecay(genpart,413, -323)==1){
return 411;}//anti-B0 decays to D*+ K*- 

if ( PcheckDecay(genpart,421, -311)==1){
return 412;}//anti-B0 decays to D0 anti-K0 

if ( PcheckDecay(genpart,-313, 421)==1){
return 413;}//anti-B0 decays to anti-K*0 D0 

if ( PcheckDecay(genpart,423, -311)==1){
return 414;}//anti-B0 decays to D*0 anti-K0 

if ( PcheckDecay(genpart,423, -313)==1){
return 415;}//anti-B0 decays to D*0 anti-K*0 

if ( PcheckDecay(genpart,-2, 3, 4, -1)==1){
return 416;}//anti-B0 decays to anti-u s c anti-d 

if ( PcheckDecay(genpart,-2, 1, 4, -1)==1){
return 417;}//anti-B0 decays to anti-u d c anti-d 

if ( PcheckDecay(genpart,4122, -2212)==1){
return 418;}//anti-B0 decays to Lambda_c+ anti-p- 

if ( PcheckDecay(genpart,4212, -2212)==1){
return 419;}//anti-B0 decays to Sigma_c+ anti-p- 

if ( PcheckDecay(genpart,4112, -2112)==1){
return 420;}//anti-B0 decays to Sigma_c0 anti-n0 

if ( PcheckDecay(genpart,4101, -2101)==1){
return 421;}//anti-B0 decays to cd_0 anti-ud_0 

if ( PcheckDecay(genpart,4103, -2103)==1){
return 422;}//anti-B0 decays to cd_1 anti-ud_1 

if ( PcheckDecay(genpart,4101, -2101)==1){
return 423;}//anti-B0 decays to cd_0 anti-ud_0 

if ( PcheckDecay(genpart,4103, -2103)==1){
return 424;}//anti-B0 decays to cd_1 anti-ud_1 

if ( PcheckDecay(genpart,4232, -4122)==1){
return 425;}//anti-B0 decays to Xi_c+ anti-Lambda_c- 

if ( PcheckDecay(genpart,4232, -4212)==1){
return 426;}//anti-B0 decays to Xi_c+ anti-Sigma_c- 

if ( PcheckDecay(genpart,4132, -4112)==1){
return 427;}//anti-B0 decays to Xi_c0 anti-Sigma_c0 

if ( PcheckDecay(genpart,4301, -4101)==1){
return 428;}//anti-B0 decays to cs_0 anti-cd_0 

if ( PcheckDecay(genpart,4232, -2212)==1){
return 429;}//anti-B0 decays to Xi_c+ anti-p- 

if ( PcheckDecay(genpart,4132, -2112)==1){
return 430;}//anti-B0 decays to Xi_c0 anti-n0 

if ( PcheckDecay(genpart,4301, -2101)==1){
return 431;}//anti-B0 decays to cs_0 anti-ud_0 

if ( PcheckDecay(genpart,4303, -2103)==1){
return 432;}//anti-B0 decays to cs_1 anti-ud_1 

if ( PcheckDecay(genpart,4301, -2101)==1){
return 433;}//anti-B0 decays to cs_0 anti-ud_0 

if ( PcheckDecay(genpart,4303, -2103)==1){
return 434;}//anti-B0 decays to cs_1 anti-ud_1 

if ( PcheckDecay(genpart,4122, -4122)==1){
return 435;}//anti-B0 decays to Lambda_c+ anti-Lambda_c- 

if ( PcheckDecay(genpart,4212, -4122)==1){
  return 436;}//anti-B0 decays to Sigma_c+ anti-Lambda_c- 
     
     if ( PcheckDecay(genpart,4122, -4212)==1){
  return 437;}//anti-B0 decays to Lambda_c+ anti-Sigma_c- 
     
     if ( PcheckDecay(genpart,4212, -4212)==1){
  return 438;}//anti-B0 decays to Sigma_c+ anti-Sigma_c- 
     
     if ( PcheckDecay(genpart,4112, -4112)==1){
  return 439;}//anti-B0 decays to Sigma_c0 anti-Sigma_c0 
     
     if ( PcheckDecay(genpart,4101, -4101)==1){
  return 440;}//anti-B0 decays to cd_0 anti-cd_0 

     return -100;
  }

 // End of anti-B0....
  
 else if (AmId(genpart)==511) {

if ( PcheckDecay(genpart,-413, -11, 12)==1){
return 441;}//B0 decays to D*- e+ nu_e 

if ( PcheckDecay(genpart,-411, -11, 12)==1){
return 442;}//B0 decays to D- e+ nu_e 

if ( PcheckDecay(genpart,-10413, -11, 12)==1){
return 443;}//B0 decays to D_1- e+ nu_e 

if ( PcheckDecay(genpart,-10411, -11, 12)==1){
return 444;}//B0 decays to D_0*- e+ nu_e 

if ( PcheckDecay(genpart,-20413, -11, 12)==1){
return 445;}//B0 decays to D'_1- e+ nu_e 

if ( PcheckDecay(genpart,-415, -11, 12)==1){
return 446;}//B0 decays to D_2*- e+ nu_e 

if ( PcheckDecay(genpart,-100411, -11, 12)==1){
return 447;}//B0 decays to D(2S)- e+ nu_e 

if ( PcheckDecay(genpart,-100413, -11, 12)==1){
return 448;}//B0 decays to D*(2S)- e+ nu_e 

if ( PcheckDecay(genpart,-413, 111, -11, 12)==1){
return 449;}//B0 decays to D*- pi0 e+ nu_e 

if ( PcheckDecay(genpart,-423, -211, -11, 12)==1){
return 450;}//B0 decays to anti-D*0 pi- e+ nu_e 

if ( PcheckDecay(genpart,-411, 111, -11, 12)==1){
return 451;}//B0 decays to D- pi0 e+ nu_e 

if ( PcheckDecay(genpart,-421, -211, -11, 12)==1){
return 452;}//B0 decays to anti-D0 pi- e+ nu_e 

if ( PcheckDecay(genpart,-413, -13, 14)==1){
return 453;}//B0 decays to D*- mu+ nu_mu 

if ( PcheckDecay(genpart,-411, -13, 14)==1){
return 454;}//B0 decays to D- mu+ nu_mu 

if ( PcheckDecay(genpart,-10413, -13, 14)==1){
return 455;}//B0 decays to D_1- mu+ nu_mu 

if ( PcheckDecay(genpart,-10411, -13, 14)==1){
return 456;}//B0 decays to D_0*- mu+ nu_mu 

if ( PcheckDecay(genpart,-20413, -13, 14)==1){
return 457;}//B0 decays to D'_1- mu+ nu_mu 

if ( PcheckDecay(genpart,-415, -13, 14)==1){
return 458;}//B0 decays to D_2*- mu+ nu_mu 

if ( PcheckDecay(genpart,-100411, -13, 14)==1){
return 459;}//B0 decays to D(2S)- mu+ nu_mu 

if ( PcheckDecay(genpart,-100413, -13, 14)==1){
return 460;}//B0 decays to D*(2S)- mu+ nu_mu 

if ( PcheckDecay(genpart,-413, 111, -13, 14)==1){
return 461;}//B0 decays to D*- pi0 mu+ nu_mu 

if ( PcheckDecay(genpart,-423, -211, -13, 14)==1){
return 462;}//B0 decays to anti-D*0 pi- mu+ nu_mu 

if ( PcheckDecay(genpart,-411, 111, -13, 14)==1){
return 463;}//B0 decays to D- pi0 mu+ nu_mu 

if ( PcheckDecay(genpart,-421, -211, -13, 14)==1){
return 464;}//B0 decays to anti-D0 pi- mu+ nu_mu 

//if ( PcheckDecay(genpart,return 465;}//B0 decays to 

if ( PcheckDecay(genpart,-413, -15, 16)==1){
return 466;}//B0 decays to D*- tau+ nu_tau 

if ( PcheckDecay(genpart,-411, -15, 16)==1){
return 467;}//B0 decays to D- tau+ nu_tau 

if ( PcheckDecay(genpart,-10413, -15, 16)==1){
return 468;}//B0 decays to D_1- tau+ nu_tau 

if ( PcheckDecay(genpart,-10411, -15, 16)==1){
return 469;}//B0 decays to D_0*- tau+ nu_tau 

if ( PcheckDecay(genpart,-20413, -15, 16)==1){
return 470;}//B0 decays to D'_1- tau+ nu_tau 

if ( PcheckDecay(genpart,-415, -15, 16)==1){
return 471;}//B0 decays to D_2*- tau+ nu_tau 

if ( PcheckDecay(genpart,441, 310)==1){
return 472;}//B0 decays to eta_c K_S0 

if ( PcheckDecay(genpart,441, 130)==1){
return 473;}//B0 decays to eta_c K_L0 

if ( PcheckDecay(genpart,313, 441)==1){
return 474;}//B0 decays to K*S eta_c 

if ( PcheckDecay(genpart,313, 441)==1){
return 475;}//B0 decays to K*L eta_c 

if ( PcheckDecay(genpart,313, 441)==1){
return 476;}//B0 decays to K*0T eta_c 

if ( PcheckDecay(genpart,441, 321, -211)==1){
return 477;}//B0 decays to eta_c K+ pi- 

if ( PcheckDecay(genpart,441, 311, 111)==1){
return 478;}//B0 decays to eta_c K0 pi0 

if ( PcheckDecay(genpart,441, 311, 211, -211)==1){
return 479;}//B0 decays to eta_c K0 pi+ pi- 

if ( PcheckDecay(genpart,441, 311, 111, 111)==1){
return 480;}//B0 decays to eta_c K0 pi0 pi0 

if ( PcheckDecay(genpart,441, 321, -211, 111)==1){
return 481;}//B0 decays to eta_c K+ pi- pi0 

if ( PcheckDecay(genpart,441, -3122, 2112)==1){
return 482;}//B0 decays to eta_c anti-Lambda0 n0 

if ( PcheckDecay(genpart,441, -3222, 2212)==1){
return 483;}//B0 decays to eta_c anti-Sigma- p+ 

if ( PcheckDecay(genpart,441, 30343)==1){
return 484;}//B0 decays to eta_c Xsd 

if ( PcheckDecay(genpart,100441, 310)==1){
return 485;}//B0 decays to eta_c(2S) K_S0 

if ( PcheckDecay(genpart,100441, 130)==1){
return 486;}//B0 decays to eta_c(2S) K_L0 

if ( PcheckDecay(genpart,313, 100441)==1){
return 487;}//B0 decays to K*S eta_c(2S) 

if ( PcheckDecay(genpart,313, 100441)==1){
return 488;}//B0 decays to K*L eta_c(2S) 

if ( PcheckDecay(genpart,313, 100441)==1){
return 489;}//B0 decays to K*0T eta_c(2S) 

if ( PcheckDecay(genpart,100441, 321, -211)==1){
return 490;}//B0 decays to eta_c(2S) K+ pi- 

if ( PcheckDecay(genpart,100441, 311, 111)==1){
return 491;}//B0 decays to eta_c(2S) K0 pi0 

if ( PcheckDecay(genpart,100441, 311, 211, -211)==1){
return 492;}//B0 decays to eta_c(2S) K0 pi+ pi- 

if ( PcheckDecay(genpart,100441, 311, 111, 111)==1){
return 493;}//B0 decays to eta_c(2S) K0 pi0 pi0 

if ( PcheckDecay(genpart,100441, 321, -211, 111)==1){
return 494;}//B0 decays to eta_c(2S) K+ pi- pi0 

if ( PcheckDecay(genpart,100441, 30343)==1){
return 495;}//B0 decays to eta_c(2S) Xsd 

if ( PcheckDecay(genpart,443, 310)==1){
return 496;}//B0 decays to psi K_S0 

if ( PcheckDecay(genpart,443, 130)==1){
return 497;}//B0 decays to psi K_L0 

if ( PcheckDecay(genpart,443, 313)==1){
return 498;}//B0 decays to psi K*S 

if ( PcheckDecay(genpart,443, 313)==1){
return 499;}//B0 decays to psi K*L 

if ( PcheckDecay(genpart,443, 313)==1){
return 500;}//B0 decays to psi K*0T 

if ( PcheckDecay(genpart,443, 321, -211)==1){
return 501;}//B0 decays to psi K+ pi- 

if ( PcheckDecay(genpart,443, 10313)==1){
return 502;}//B0 decays to psi K_10 

if ( PcheckDecay(genpart,443, 20313)==1){
return 503;}//B0 decays to psi K'_10 

if ( PcheckDecay(genpart,443, 311, 113)==1){
return 504;}//B0 decays to psi K0 rho0 

if ( PcheckDecay(genpart,443, 311, 211, -211)==1){
return 505;}//B0 decays to psi K0 pi+ pi- 

if ( PcheckDecay(genpart,443, 323, -211)==1){
return 506;}//B0 decays to psi K*+ pi- 

if ( PcheckDecay(genpart,443, 313, 211, -211)==1){
return 507;}//B0 decays to psi K*0 pi+ pi- 

if ( PcheckDecay(genpart,443, 333, 311)==1){
return 508;}//B0 decays to psi phi K0 

if ( PcheckDecay(genpart,443, 111)==1){
return 509;}//B0 decays to psi pi0 

if ( PcheckDecay(genpart,443, 221)==1){
return 510;}//B0 decays to psi eta 

if ( PcheckDecay(genpart,443, 211, -211)==1){
return 511;}//B0 decays to psi pi+ pi- 

if ( PcheckDecay(genpart,443, 113)==1){
return 512;}//B0 decays to psi rho0 

if ( PcheckDecay(genpart,443, 223)==1){
return 513;}//B0 decays to psi omega 

if ( PcheckDecay(genpart,443, 10311)==1){
return 514;}//B0 decays to psi K_0*0 

if ( PcheckDecay(genpart,443, 315)==1){
return 515;}//B0 decays to psi K_2*0 

if ( PcheckDecay(genpart,443, 30313)==1){
return 516;}//B0 decays to psi K''*0 

if ( PcheckDecay(genpart,443, -3122, 2112)==1){
return 517;}//B0 decays to psi anti-Lambda0 n0 

if ( PcheckDecay(genpart,443, -3222, 2212)==1){
return 518;}//B0 decays to psi anti-Sigma- p+ 

if ( PcheckDecay(genpart,443, 30343)==1){
return 519;}//B0 decays to psi Xsd 

if ( PcheckDecay(genpart,100443, 310)==1){
return 520;}//B0 decays to psi(2S) K_S0 

if ( PcheckDecay(genpart,100443, 130)==1){
return 521;}//B0 decays to psi(2S) K_L0 

if ( PcheckDecay(genpart,100443, 313)==1){
return 522;}//B0 decays to psi(2S) K*S 

if ( PcheckDecay(genpart,100443, 313)==1){
return 523;}//B0 decays to psi(2S) K*L 

if ( PcheckDecay(genpart,100443, 313)==1){
return 524;}//B0 decays to psi(2S) K*0T 

if ( PcheckDecay(genpart,100443, 111)==1){
return 525;}//B0 decays to psi(2S) pi0 

if ( PcheckDecay(genpart,100443, 30343)==1){
return 526;}//B0 decays to psi(2S) Xsd 

if ( PcheckDecay(genpart,10441, 310)==1){
return 527;}//B0 decays to chi_c0 K_S0 

if ( PcheckDecay(genpart,10441, 130)==1){
return 528;}//B0 decays to chi_c0 K_L0 

if ( PcheckDecay(genpart,313, 10441)==1){
return 529;}//B0 decays to K*S chi_c0 

if ( PcheckDecay(genpart,313, 10441)==1){
return 530;}//B0 decays to K*L chi_c0 

if ( PcheckDecay(genpart,313, 10441)==1){
return 531;}//B0 decays to K*0T chi_c0 

if ( PcheckDecay(genpart,10441, 30343)==1){
return 532;}//B0 decays to chi_c0 Xsd 

if ( PcheckDecay(genpart,20443, 310)==1){
return 533;}//B0 decays to chi_c1 K_S0 

if ( PcheckDecay(genpart,20443, 130)==1){
return 534;}//B0 decays to chi_c1 K_L0 

if ( PcheckDecay(genpart,20443, 313)==1){
return 535;}//B0 decays to chi_c1 K*S 

if ( PcheckDecay(genpart,20443, 313)==1){
return 536;}//B0 decays to chi_c1 K*L 

if ( PcheckDecay(genpart,20443, 313)==1){
return 537;}//B0 decays to chi_c1 K*0T 

if ( PcheckDecay(genpart,20443, 111)==1){
return 538;}//B0 decays to chi_c1 pi0 

if ( PcheckDecay(genpart,20443, 321, -211)==1){
return 539;}//B0 decays to chi_c1 K+ pi- 

if ( PcheckDecay(genpart,20443, 311, 111)==1){
return 540;}//B0 decays to chi_c1 K0 pi0 

if ( PcheckDecay(genpart,20443, 311, 211, -211)==1){
return 541;}//B0 decays to chi_c1 K0 pi+ pi- 

if ( PcheckDecay(genpart,20443, 311, 111, 111)==1){
return 542;}//B0 decays to chi_c1 K0 pi0 pi0 

if ( PcheckDecay(genpart,20443, 321, -211, 111)==1){
return 543;}//B0 decays to chi_c1 K+ pi- pi0 

if ( PcheckDecay(genpart,20443, 30343)==1){
return 544;}//B0 decays to chi_c1 Xsd 

if ( PcheckDecay(genpart,445, 310)==1){
return 545;}//B0 decays to chi_c2 K_S0 

if ( PcheckDecay(genpart,445, 130)==1){
return 546;}//B0 decays to chi_c2 K_L0 

if ( PcheckDecay(genpart,445, 313)==1){
return 547;}//B0 decays to chi_c2 K*S 

if ( PcheckDecay(genpart,445, 313)==1){
return 548;}//B0 decays to chi_c2 K*L 

if ( PcheckDecay(genpart,445, 313)==1){
return 549;}//B0 decays to chi_c2 K*0T 

if ( PcheckDecay(genpart,445, 321, -211)==1){
return 550;}//B0 decays to chi_c2 K+ pi- 

if ( PcheckDecay(genpart,445, 311, 111)==1){
return 551;}//B0 decays to chi_c2 K0 pi0 

if ( PcheckDecay(genpart,445, 311, 211, -211)==1){
return 552;}//B0 decays to chi_c2 K0 pi+ pi- 

if ( PcheckDecay(genpart,445, 311, 111, 111)==1){
return 553;}//B0 decays to chi_c2 K0 pi0 pi0 

if ( PcheckDecay(genpart,445, 321, -211, 111)==1){
return 554;}//B0 decays to chi_c2 K+ pi- pi0 

if ( PcheckDecay(genpart,445, 30343)==1){
return 555;}//B0 decays to chi_c2 Xsd 

if ( PcheckDecay(genpart,10443, 310)==1){
return 556;}//B0 decays to h_c K_S0 

if ( PcheckDecay(genpart,10443, 130)==1){
return 557;}//B0 decays to h_c K_L0 

if ( PcheckDecay(genpart,10443, 313)==1){
return 558;}//B0 decays to h_c K*S 

if ( PcheckDecay(genpart,10443, 313)==1){
return 559;}//B0 decays to h_c K*L 

if ( PcheckDecay(genpart,10443, 313)==1){
return 560;}//B0 decays to h_c K*0T 

if ( PcheckDecay(genpart,10443, 30343)==1){
return 561;}//B0 decays to h_c Xsd 

if ( PcheckDecay(genpart,120443, 310)==1){
return 562;}//B0 decays to X(3872) K_S0 

if ( PcheckDecay(genpart,120443, 130)==1){
return 563;}//B0 decays to X(3872) K_L0 

if ( PcheckDecay(genpart,90000443, 310)==1){
return 564;}//B0 decays to Y(3940) K_S0 

if ( PcheckDecay(genpart,90000443, 130)==1){
return 565;}//B0 decays to Y(3940) K_L0 

if ( PcheckDecay(genpart,-411, 431)==1){
return 566;}//B0 decays to D- D_s+ 

if ( PcheckDecay(genpart,-413, 431)==1){
return 567;}//B0 decays to D*- D_s+ 

if ( PcheckDecay(genpart,433, -411)==1){
return 568;}//B0 decays to D_s*+ D- 

if ( PcheckDecay(genpart,433, -413)==1){
return 569;}//B0 decays to D_s*+ D*- 

if ( PcheckDecay(genpart,431, -10411)==1){
return 570;}//B0 decays to D_s+ D_0*- 

if ( PcheckDecay(genpart,433, -10411)==1){
return 571;}//B0 decays to D_s*+ D_0*- 

if ( PcheckDecay(genpart,-20413, 431)==1){
return 572;}//B0 decays to D'_1- D_s+ 

if ( PcheckDecay(genpart,-20413, 433)==1){
return 573;}//B0 decays to D'_1- D_s*+ 

if ( PcheckDecay(genpart,-10413, 431)==1){
return 574;}//B0 decays to D_1- D_s+ 

if ( PcheckDecay(genpart,-10413, 433)==1){
return 575;}//B0 decays to D_1- D_s*+ 

if ( PcheckDecay(genpart,-415, 431)==1){
return 576;}//B0 decays to D_2*- D_s+ 

if ( PcheckDecay(genpart,-415, 433)==1){
return 577;}//B0 decays to D_2*- D_s*+ 

if ( PcheckDecay(genpart,-411, 10431)==1){
return 578;}//B0 decays to D- D_s0*+ 

if ( PcheckDecay(genpart,-413, 10431)==1){
return 579;}//B0 decays to D*- D_s0*+ 

if ( PcheckDecay(genpart,20433, -411)==1){
return 580;}//B0 decays to D'_s1+ D- 

if ( PcheckDecay(genpart,20433, -413)==1){
return 581;}//B0 decays to D'_s1+ D*- 

if ( PcheckDecay(genpart,431, -411, 111)==1){
return 582;}//B0 decays to D_s+ D- pi0 

if ( PcheckDecay(genpart,431, -421, -211)==1){
return 583;}//B0 decays to D_s+ anti-D0 pi- 

if ( PcheckDecay(genpart,433, -411, 111)==1){
return 584;}//B0 decays to D_s*+ D- pi0 

if ( PcheckDecay(genpart,433, -421, -211)==1){
return 585;}//B0 decays to D_s*+ anti-D0 pi- 

if ( PcheckDecay(genpart,431, -411, -211, 211)==1){
return 586;}//B0 decays to D_s+ D- pi- pi+ 

if ( PcheckDecay(genpart,431, -411, 111, 111)==1){
return 587;}//B0 decays to D_s+ D- pi0 pi0 

if ( PcheckDecay(genpart,431, -421, -211, 111)==1){
return 588;}//B0 decays to D_s+ anti-D0 pi- pi0 

if ( PcheckDecay(genpart,433, -411, -211, 211)==1){
return 589;}//B0 decays to D_s*+ D- pi- pi+ 

if ( PcheckDecay(genpart,433, -411, 111, 111)==1){
return 590;}//B0 decays to D_s*+ D- pi0 pi0 

if ( PcheckDecay(genpart,433, -421, -211, 111)==1){
return 591;}//B0 decays to D_s*+ anti-D0 pi- pi0 

if ( PcheckDecay(genpart,-431, 211, 311)==1){
return 592;}//B0 decays to D_s- pi+ K0 

if ( PcheckDecay(genpart,-433, 211, 311)==1){
return 593;}//B0 decays to D_s*- pi+ K0 

if ( PcheckDecay(genpart,-431, 211, 311, 111)==1){
return 594;}//B0 decays to D_s- pi+ K0 pi0 

if ( PcheckDecay(genpart,-433, 211, 311, 111)==1){
return 595;}//B0 decays to D_s*- pi+ K0 pi0 

if ( PcheckDecay(genpart,-431, 211, 321, -211)==1){
return 596;}//B0 decays to D_s- pi+ K+ pi- 

if ( PcheckDecay(genpart,-433, 211, 321, -211)==1){
return 597;}//B0 decays to D_s*- pi+ K+ pi- 

if ( PcheckDecay(genpart,-431, 211, 311, 111, 111)==1){
return 598;}//B0 decays to D_s- pi+ K0 pi0 pi0 

if ( PcheckDecay(genpart,-433, 211, 311, 111, 111)==1){
return 599;}//B0 decays to D_s*- pi+ K0 pi0 pi0 

if ( PcheckDecay(genpart,-431, 211, 311, 211, -211)==1){
return 600;}//B0 decays to D_s- pi+ K0 pi+ pi- 

if ( PcheckDecay(genpart,-433, 211, 311, 211, -211)==1){
return 601;}//B0 decays to D_s*- pi+ K0 pi+ pi- 

if ( PcheckDecay(genpart,-431, 211, 321, -211, 111)==1){
return 602;}//B0 decays to D_s- pi+ K+ pi- pi0 

if ( PcheckDecay(genpart,-433, 211, 321, -211, 111)==1){
return 603;}//B0 decays to D_s*- pi+ K+ pi- pi0 

if ( PcheckDecay(genpart,-411, 411)==1){
return 604;}//B0 decays to D- D+ 

if ( PcheckDecay(genpart,413, -411)==1){
return 605;}//B0 decays to D*+ D- 

if ( PcheckDecay(genpart,-413, 411)==1){
return 606;}//B0 decays to D*- D+ 

if ( PcheckDecay(genpart,413, -413)==1){
return 607;}//B0 decays to D*+ D*- 

if ( PcheckDecay(genpart,-20413, 411)==1){
return 608;}//B0 decays to D'_1- D+ 

if ( PcheckDecay(genpart,20413, -411)==1){
return 609;}//B0 decays to D'_1+ D- 

if ( PcheckDecay(genpart,-20413, 413)==1){
return 610;}//B0 decays to D'_1- D*+ 

if ( PcheckDecay(genpart,20413, -413)==1){
return 611;}//B0 decays to D'_1+ D*- 

if ( PcheckDecay(genpart,-10413, 411)==1){
return 612;}//B0 decays to D_1- D+ 

if ( PcheckDecay(genpart,10413, -411)==1){
return 613;}//B0 decays to D_1+ D- 

if ( PcheckDecay(genpart,-10413, 413)==1){
return 614;}//B0 decays to D_1- D*+ 

if ( PcheckDecay(genpart,10413, -413)==1){
return 615;}//B0 decays to D_1+ D*- 

if ( PcheckDecay(genpart,-415, 411)==1){
return 616;}//B0 decays to D_2*- D+ 

if ( PcheckDecay(genpart,415, -411)==1){
return 617;}//B0 decays to D_2*+ D- 

if ( PcheckDecay(genpart,-415, 413)==1){
return 618;}//B0 decays to D_2*- D*+ 

if ( PcheckDecay(genpart,415, -413)==1){
return 619;}//B0 decays to D_2*+ D*- 

if ( PcheckDecay(genpart,10433, -411)==1){
return 620;}//B0 decays to D_s1+ D- 

if ( PcheckDecay(genpart,10433, -413)==1){
return 621;}//B0 decays to D_s1+ D*- 

if ( PcheckDecay(genpart,435, -411)==1){
return 622;}//B0 decays to D_s2*+ D- 

if ( PcheckDecay(genpart,435, -413)==1){
return 623;}//B0 decays to D_s2*+ D*- 

if ( PcheckDecay(genpart,9000433, -411)==1){
return 624;}//B0 decays to D_sj(2700)+ D- 

if ( PcheckDecay(genpart,9000433, -413)==1){
return 625;}//B0 decays to D_sj(2700)+ D*- 

if ( PcheckDecay(genpart,30443, 310)==1){
return 626;}//B0 decays to psi(3770) K_S0 

if ( PcheckDecay(genpart,30443, 130)==1){
return 627;}//B0 decays to psi(3770) K_L0 

if ( PcheckDecay(genpart,30443, 313)==1){
return 628;}//B0 decays to psi(3770) K*S 

if ( PcheckDecay(genpart,30443, 313)==1){
return 629;}//B0 decays to psi(3770) K*L 

if ( PcheckDecay(genpart,30443, 313)==1){
return 630;}//B0 decays to psi(3770) K*0T 

if ( PcheckDecay(genpart,30443, 321, -211)==1){
return 631;}//B0 decays to psi(3770) K+ pi- 

if ( PcheckDecay(genpart,30443, 311, 111)==1){
return 632;}//B0 decays to psi(3770) K0 pi0 

if ( PcheckDecay(genpart,30443, 311, 211, -211)==1){
return 633;}//B0 decays to psi(3770) K0 pi+ pi- 

if ( PcheckDecay(genpart,30443, 311, 111, 111)==1){
return 634;}//B0 decays to psi(3770) K0 pi0 pi0 

if ( PcheckDecay(genpart,30443, 321, -211, 111)==1){
return 635;}//B0 decays to psi(3770) K+ pi- pi0 

if ( PcheckDecay(genpart,9000443, 310)==1){
return 636;}//B0 decays to psi(4040) K_S0 

if ( PcheckDecay(genpart,9000443, 130)==1){
return 637;}//B0 decays to psi(4040) K_L0 

if ( PcheckDecay(genpart,9000443, 313)==1){
return 638;}//B0 decays to psi(4040) K*S 

if ( PcheckDecay(genpart,9000443, 313)==1){
return 639;}//B0 decays to psi(4040) K*L 

if ( PcheckDecay(genpart,9000443, 313)==1){
return 640;}//B0 decays to psi(4040) K*0T 

if ( PcheckDecay(genpart,9000443, 321, -211)==1){
return 641;}//B0 decays to psi(4040) K+ pi- 

if ( PcheckDecay(genpart,9000443, 311, 111)==1){
return 642;}//B0 decays to psi(4040) K0 pi0 

if ( PcheckDecay(genpart,9010443, 310)==1){
return 643;}//B0 decays to psi(4160) K_S0 

if ( PcheckDecay(genpart,9010443, 130)==1){
return 644;}//B0 decays to psi(4160) K_L0 

if ( PcheckDecay(genpart,9010443, 313)==1){
return 645;}//B0 decays to psi(4160) K*S 

if ( PcheckDecay(genpart,9010443, 313)==1){
return 646;}//B0 decays to psi(4160) K*L 

if ( PcheckDecay(genpart,9010443, 313)==1){
return 647;}//B0 decays to psi(4160) K*0T 

if ( PcheckDecay(genpart,9010443, 321, -211)==1){
return 648;}//B0 decays to psi(4160) K+ pi- 

if ( PcheckDecay(genpart,9010443, 311, 111)==1){
return 649;}//B0 decays to psi(4160) K0 pi0 

//if ( PcheckDecay(genpart,return 650;}//B0 decays to 

if ( PcheckDecay(genpart,411, -411, 111)==1){
return 651;}//B0 decays to D+ D- pi0 

if ( PcheckDecay(genpart,411, -421, -211)==1){
return 652;}//B0 decays to D+ anti-D0 pi- 

if ( PcheckDecay(genpart,421, -411, 211)==1){
return 653;}//B0 decays to D0 D- pi+ 

if ( PcheckDecay(genpart,413, -411, 111)==1){
return 654;}//B0 decays to D*+ D- pi0 

if ( PcheckDecay(genpart,413, -421, -211)==1){
return 655;}//B0 decays to D*+ anti-D0 pi- 

if ( PcheckDecay(genpart,423, -411, 211)==1){
return 656;}//B0 decays to D*0 D- pi+ 

if ( PcheckDecay(genpart,411, -413, 111)==1){
return 657;}//B0 decays to D+ D*- pi0 

if ( PcheckDecay(genpart,411, -423, -211)==1){
return 658;}//B0 decays to D+ anti-D*0 pi- 

if ( PcheckDecay(genpart,421, -413, 211)==1){
return 659;}//B0 decays to D0 D*- pi+ 

if ( PcheckDecay(genpart,413, -413, 111)==1){
return 660;}//B0 decays to D*+ D*- pi0 

if ( PcheckDecay(genpart,413, -423, -211)==1){
return 661;}//B0 decays to D*+ anti-D*0 pi- 

if ( PcheckDecay(genpart,423, -413, 211)==1){
return 662;}//B0 decays to D*0 D*- pi+ 

if ( PcheckDecay(genpart,411, -411, 113)==1){
return 663;}//B0 decays to D+ D- rho0 

if ( PcheckDecay(genpart,411, -421, -213)==1){
return 664;}//B0 decays to D+ anti-D0 rho- 

if ( PcheckDecay(genpart,421, -411, 213)==1){
return 665;}//B0 decays to D0 D- rho+ 

if ( PcheckDecay(genpart,413, -411, 113)==1){
return 666;}//B0 decays to D*+ D- rho0 

if ( PcheckDecay(genpart,413, -421, -213)==1){
return 667;}//B0 decays to D*+ anti-D0 rho- 

if ( PcheckDecay(genpart,423, -411, 213)==1){
return 668;}//B0 decays to D*0 D- rho+ 

if ( PcheckDecay(genpart,411, -413, 113)==1){
return 669;}//B0 decays to D+ D*- rho0 

if ( PcheckDecay(genpart,411, -423, -213)==1){
return 670;}//B0 decays to D+ anti-D*0 rho- 

if ( PcheckDecay(genpart,421, -413, 213)==1){
return 671;}//B0 decays to D0 D*- rho+ 

if ( PcheckDecay(genpart,413, -413, 113)==1){
return 672;}//B0 decays to D*+ D*- rho0 

if ( PcheckDecay(genpart,413, -423, -213)==1){
return 673;}//B0 decays to D*+ anti-D*0 rho- 

if ( PcheckDecay(genpart,423, -413, 213)==1){
return 674;}//B0 decays to D*0 D*- rho+ 

if ( PcheckDecay(genpart,411, -411, 111, 111)==1){
return 675;}//B0 decays to D+ D- pi0 pi0 

if ( PcheckDecay(genpart,411, -411, -211, 211)==1){
return 676;}//B0 decays to D+ D- pi- pi+ 

if ( PcheckDecay(genpart,411, -421, -211, 111)==1){
return 677;}//B0 decays to D+ anti-D0 pi- pi0 

if ( PcheckDecay(genpart,421, -411, 211, 111)==1){
return 678;}//B0 decays to D0 D- pi+ pi0 

if ( PcheckDecay(genpart,421, -421, -211, 211)==1){
return 679;}//B0 decays to D0 anti-D0 pi- pi+ 

if ( PcheckDecay(genpart,413, -411, 111, 111)==1){
return 680;}//B0 decays to D*+ D- pi0 pi0 

if ( PcheckDecay(genpart,413, -411, -211, 211)==1){
return 681;}//B0 decays to D*+ D- pi- pi+ 

if ( PcheckDecay(genpart,413, -421, -211, 111)==1){
return 682;}//B0 decays to D*+ anti-D0 pi- pi0 

if ( PcheckDecay(genpart,423, -411, 211, 111)==1){
return 683;}//B0 decays to D*0 D- pi+ pi0 

if ( PcheckDecay(genpart,423, -421, -211, 211)==1){
return 684;}//B0 decays to D*0 anti-D0 pi- pi+ 

if ( PcheckDecay(genpart,411, -413, 111, 111)==1){
return 685;}//B0 decays to D+ D*- pi0 pi0 

if ( PcheckDecay(genpart,411, -413, -211, 211)==1){
return 686;}//B0 decays to D+ D*- pi- pi+ 

if ( PcheckDecay(genpart,411, -423, -211, 111)==1){
return 687;}//B0 decays to D+ anti-D*0 pi- pi0 

if ( PcheckDecay(genpart,421, -413, 211, 111)==1){
return 688;}//B0 decays to D0 D*- pi+ pi0 

if ( PcheckDecay(genpart,421, -423, -211, 211)==1){
return 689;}//B0 decays to D0 anti-D*0 pi- pi+ 

if ( PcheckDecay(genpart,413, -413, 111, 111)==1){
return 690;}//B0 decays to D*+ D*- pi0 pi0 

if ( PcheckDecay(genpart,413, -413, -211, 211)==1){
return 691;}//B0 decays to D*+ D*- pi- pi+ 

if ( PcheckDecay(genpart,413, -423, -211, 111)==1){
return 692;}//B0 decays to D*+ anti-D*0 pi- pi0 

if ( PcheckDecay(genpart,423, -413, 211, 111)==1){
return 693;}//B0 decays to D*0 D*- pi+ pi0 

if ( PcheckDecay(genpart,423, -423, -211, 211)==1){
return 694;}//B0 decays to D*0 anti-D*0 pi- pi+ 

if ( PcheckDecay(genpart,-411, 421, 321)==1){
return 695;}//B0 decays to D- D0 K+ 

if ( PcheckDecay(genpart,-411, 411, 311)==1){
return 696;}//B0 decays to D- D+ K0 

if ( PcheckDecay(genpart,-411, 421, 323)==1){
return 697;}//B0 decays to D- D0 K*+ 

if ( PcheckDecay(genpart,-411, 411, 313)==1){
return 698;}//B0 decays to D- D+ K*0 

if ( PcheckDecay(genpart,-413, 421, 321)==1){
return 699;}//B0 decays to D*- D0 K+ 

if ( PcheckDecay(genpart,-413, 411, 311)==1){
return 700;}//B0 decays to D*- D+ K0 

if ( PcheckDecay(genpart,-411, 423, 321)==1){
return 701;}//B0 decays to D- D*0 K+ 

if ( PcheckDecay(genpart,-413, 421, 323)==1){
return 702;}//B0 decays to D*- D0 K*+ 

if ( PcheckDecay(genpart,-413, 411, 313)==1){
return 703;}//B0 decays to D*- D+ K*0 

if ( PcheckDecay(genpart,-411, 423, 323)==1){
return 704;}//B0 decays to D- D*0 K*+ 

if ( PcheckDecay(genpart,-411, 413, 313)==1){
return 705;}//B0 decays to D- D*+ K*0 

if ( PcheckDecay(genpart,-413, 423, 321)==1){
return 706;}//B0 decays to D*- D*0 K+ 

if ( PcheckDecay(genpart,-413, 413, 311)==1){
return 707;}//B0 decays to D*- D*+ K0 

if ( PcheckDecay(genpart,-413, 423, 323)==1){
return 708;}//B0 decays to D*- D*0 K*+ 

if ( PcheckDecay(genpart,-413, 413, 313)==1){
return 709;}//B0 decays to D*- D*+ K*0 

if ( PcheckDecay(genpart,-10413, 421, 321)==1){
return 710;}//B0 decays to D_1- D0 K+ 

if ( PcheckDecay(genpart,-10413, 411, 311)==1){
return 711;}//B0 decays to D_1- D+ K0 

if ( PcheckDecay(genpart,-411, 10423, 321)==1){
return 712;}//B0 decays to D- D_10 K+ 

if ( PcheckDecay(genpart,-411, 10413, 311)==1){
return 713;}//B0 decays to D- D_1+ K0 

if ( PcheckDecay(genpart,-10413, 423, 321)==1){
return 714;}//B0 decays to D_1- D*0 K+ 

if ( PcheckDecay(genpart,-10413, 413, 311)==1){
return 715;}//B0 decays to D_1- D*+ K0 

if ( PcheckDecay(genpart,-413, 10423, 321)==1){
return 716;}//B0 decays to D*- D_10 K+ 

if ( PcheckDecay(genpart,-413, 10413, 311)==1){
return 717;}//B0 decays to D*- D_1+ K0 

if ( PcheckDecay(genpart,-415, 421, 321)==1){
return 718;}//B0 decays to D_2*- D0 K+ 

if ( PcheckDecay(genpart,-415, 411, 311)==1){
return 719;}//B0 decays to D_2*- D+ K0 

if ( PcheckDecay(genpart,-411, 425, 321)==1){
return 720;}//B0 decays to D- D_2*0 K+ 

if ( PcheckDecay(genpart,-411, 415, 311)==1){
return 721;}//B0 decays to D- D_2*+ K0 

if ( PcheckDecay(genpart,-415, 423, 321)==1){
return 722;}//B0 decays to D_2*- D*0 K+ 

if ( PcheckDecay(genpart,-415, 413, 311)==1){
return 723;}//B0 decays to D_2*- D*+ K0 

if ( PcheckDecay(genpart,-413, 425, 321)==1){
return 724;}//B0 decays to D*- D_2*0 K+ 

if ( PcheckDecay(genpart,-413, 415, 311)==1){
return 725;}//B0 decays to D*- D_2*+ K0 

if ( PcheckDecay(genpart,-411, 421, 321, 111)==1){
return 726;}//B0 decays to D- D0 K+ pi0 

if ( PcheckDecay(genpart,-421, 421, 321, -211)==1){
return 727;}//B0 decays to anti-D0 D0 K+ pi- 

if ( PcheckDecay(genpart,-411, 421, 311, 211)==1){
return 728;}//B0 decays to D- D0 K0 pi+ 

if ( PcheckDecay(genpart,-411, 411, 321, -211)==1){
return 729;}//B0 decays to D- D+ K+ pi- 

if ( PcheckDecay(genpart,-411, 411, 311, 111)==1){
return 730;}//B0 decays to D- D+ K0 pi0 

if ( PcheckDecay(genpart,-421, 411, 311, -211)==1){
return 731;}//B0 decays to anti-D0 D+ K0 pi- 

if ( PcheckDecay(genpart,-413, 421, 321, 111)==1){
return 732;}//B0 decays to D*- D0 K+ pi0 

if ( PcheckDecay(genpart,-423, 421, 321, -211)==1){
return 733;}//B0 decays to anti-D*0 D0 K+ pi- 

if ( PcheckDecay(genpart,-413, 421, 311, 211)==1){
return 734;}//B0 decays to D*- D0 K0 pi+ 

if ( PcheckDecay(genpart,-413, 411, 321, -211)==1){
return 735;}//B0 decays to D*- D+ K+ pi- 

if ( PcheckDecay(genpart,-413, 411, 311, 111)==1){
return 736;}//B0 decays to D*- D+ K0 pi0 

if ( PcheckDecay(genpart,-423, 411, 311, -211)==1){
return 737;}//B0 decays to anti-D*0 D+ K0 pi- 

if ( PcheckDecay(genpart,-411, 423, 321, 111)==1){
return 738;}//B0 decays to D- D*0 K+ pi0 

if ( PcheckDecay(genpart,-421, 423, 321, -211)==1){
return 739;}//B0 decays to anti-D0 D*0 K+ pi- 

if ( PcheckDecay(genpart,-411, 423, 311, 211)==1){
return 740;}//B0 decays to D- D*0 K0 pi+ 

if ( PcheckDecay(genpart,-411, 413, 321, -211)==1){
return 741;}//B0 decays to D- D*+ K+ pi- 

if ( PcheckDecay(genpart,-411, 413, 311, 111)==1){
return 742;}//B0 decays to D- D*+ K0 pi0 

if ( PcheckDecay(genpart,-421, 413, 311, -211)==1){
return 743;}//B0 decays to anti-D0 D*+ K0 pi- 

if ( PcheckDecay(genpart,-413, 423, 321, 111)==1){
return 744;}//B0 decays to D*- D*0 K+ pi0 

if ( PcheckDecay(genpart,-423, 423, 321, -211)==1){
return 745;}//B0 decays to anti-D*0 D*0 K+ pi- 

if ( PcheckDecay(genpart,-413, 423, 311, 211)==1){
return 746;}//B0 decays to D*- D*0 K0 pi+ 

if ( PcheckDecay(genpart,-413, 413, 321, -211)==1){
return 747;}//B0 decays to D*- D*+ K+ pi- 

if ( PcheckDecay(genpart,-413, 413, 311, 111)==1){
return 748;}//B0 decays to D*- D*+ K0 pi0 

if ( PcheckDecay(genpart,-423, 413, 311, -211)==1){
return 749;}//B0 decays to anti-D*0 D*+ K0 pi- 

if ( PcheckDecay(genpart,-411, 421, 323, 111)==1){
return 750;}//B0 decays to D- D0 K*+ pi0 

if ( PcheckDecay(genpart,-421, 421, 323, -211)==1){
return 751;}//B0 decays to anti-D0 D0 K*+ pi- 

if ( PcheckDecay(genpart,-411, 421, 313, 211)==1){
return 752;}//B0 decays to D- D0 K*0 pi+ 

if ( PcheckDecay(genpart,-411, 411, 323, -211)==1){
return 753;}//B0 decays to D- D+ K*+ pi- 

if ( PcheckDecay(genpart,-411, 411, 313, 111)==1){
return 754;}//B0 decays to D- D+ K*0 pi0 

if ( PcheckDecay(genpart,-421, 411, 313, -211)==1){
return 755;}//B0 decays to anti-D0 D+ K*0 pi- 

if ( PcheckDecay(genpart,-411, 421, 321, 113)==1){
return 756;}//B0 decays to D- D0 K+ rho0 

if ( PcheckDecay(genpart,-421, 421, 321, -213)==1){
return 757;}//B0 decays to anti-D0 D0 K+ rho- 

if ( PcheckDecay(genpart,-411, 421, 311, 213)==1){
return 758;}//B0 decays to D- D0 K0 rho+ 

if ( PcheckDecay(genpart,-411, 411, 321, -213)==1){
return 759;}//B0 decays to D- D+ K+ rho- 

if ( PcheckDecay(genpart,-411, 411, 311, 113)==1){
return 760;}//B0 decays to D- D+ K0 rho0 

if ( PcheckDecay(genpart,-421, 411, 311, -213)==1){
return 761;}//B0 decays to anti-D0 D+ K0 rho- 

if ( PcheckDecay(genpart,421, -421, 311)==1){
return 762;}//B0 decays to D0 anti-D0 K0 

if ( PcheckDecay(genpart,421, -421, 313)==1){
return 763;}//B0 decays to D0 anti-D0 K*0 

if ( PcheckDecay(genpart,421, -423, 311)==1){
return 764;}//B0 decays to D0 anti-D*0 K0 

if ( PcheckDecay(genpart,423, -421, 311)==1){
return 765;}//B0 decays to D*0 anti-D0 K0 

if ( PcheckDecay(genpart,421, -423, 313)==1){
return 766;}//B0 decays to D0 anti-D*0 K*0 

if ( PcheckDecay(genpart,423, -421, 313)==1){
return 767;}//B0 decays to D*0 anti-D0 K*0 

if ( PcheckDecay(genpart,423, -423, 311)==1){
return 768;}//B0 decays to D*0 anti-D*0 K0 

if ( PcheckDecay(genpart,423, -423, 313)==1){
return 769;}//B0 decays to D*0 anti-D*0 K*0 

if ( PcheckDecay(genpart,10423, -421, 311)==1){
return 770;}//B0 decays to D_10 anti-D0 K0 

if ( PcheckDecay(genpart,421, -10423, 311)==1){
return 771;}//B0 decays to D0 anti-D_10 K0 

if ( PcheckDecay(genpart,10423, -423, 311)==1){
return 772;}//B0 decays to D_10 anti-D*0 K0 

if ( PcheckDecay(genpart,423, -10423, 311)==1){
return 773;}//B0 decays to D*0 anti-D_10 K0 

if ( PcheckDecay(genpart,425, -421, 311)==1){
return 774;}//B0 decays to D_2*0 anti-D0 K0 

if ( PcheckDecay(genpart,421, -425, 311)==1){
return 775;}//B0 decays to D0 anti-D_2*0 K0 

if ( PcheckDecay(genpart,425, -423, 311)==1){
return 776;}//B0 decays to D_2*0 anti-D*0 K0 

if ( PcheckDecay(genpart,423, -425, 311)==1){
return 777;}//B0 decays to D*0 anti-D_2*0 K0 

if ( PcheckDecay(genpart,-411, 211)==1){
return 778;}//B0 decays to D- pi+ 

if ( PcheckDecay(genpart,-411, 211, 111)==1){
return 779;}//B0 decays to D- pi+ pi0 

if ( PcheckDecay(genpart,213, -411)==1){
return 780;}//B0 decays to rho+ D- 

if ( PcheckDecay(genpart,-411, -211, 211, 211)==1){
return 781;}//B0 decays to D- pi- pi+ pi+ 

if ( PcheckDecay(genpart,-411, 111, 211, 111)==1){
return 782;}//B0 decays to D- pi0 pi+ pi0 

if ( PcheckDecay(genpart,-411, 113, 211)==1){
return 783;}//B0 decays to D- rho0 pi+ 

if ( PcheckDecay(genpart,-411, 213, 111)==1){
return 784;}//B0 decays to D- rho+ pi0 

if ( PcheckDecay(genpart,20213, -411)==1){
return 785;}//B0 decays to a_1+ D- 

if ( PcheckDecay(genpart,10213, -411)==1){
return 786;}//B0 decays to b_1+ D- 

if ( PcheckDecay(genpart,-411, -211, 211, 211, 111)==1){
return 787;}//B0 decays to D- pi- pi+ pi+ pi0 

if ( PcheckDecay(genpart,-411, 223, 211)==1){
return 788;}//B0 decays to D- omega pi+ 

if ( PcheckDecay(genpart,-411, 2212, -2212, 211)==1){
return 789;}//B0 decays to D- p+ anti-p- pi+ 

if ( PcheckDecay(genpart,-411, 2212, -2112)==1){
return 790;}//B0 decays to D- p+ anti-n0 

if ( PcheckDecay(genpart,-413, 211)==1){
return 791;}//B0 decays to D*- pi+ 

if ( PcheckDecay(genpart,-413, 211, 111)==1){
return 792;}//B0 decays to D*- pi+ pi0 

if ( PcheckDecay(genpart,213, -413)==1){
return 793;}//B0 decays to rho+ D*- 

if ( PcheckDecay(genpart,-413, -211, 211, 211)==1){
return 794;}//B0 decays to D*- pi- pi+ pi+ 

if ( PcheckDecay(genpart,-413, 111, 211, 111)==1){
return 795;}//B0 decays to D*- pi0 pi+ pi0 

if ( PcheckDecay(genpart,-413, 113, 211)==1){
return 796;}//B0 decays to D*- rho0 pi+ 

if ( PcheckDecay(genpart,-413, 213, 111)==1){
return 797;}//B0 decays to D*- rho+ pi0 

if ( PcheckDecay(genpart,-413, 20213)==1){
return 798;}//B0 decays to D*- a_1+ 

if ( PcheckDecay(genpart,-413, 10213)==1){
return 799;}//B0 decays to D*- b_1+ 

if ( PcheckDecay(genpart,-413, -211, 211, 211, 111)==1){
return 800;}//B0 decays to D*- pi- pi+ pi+ pi0 

if ( PcheckDecay(genpart,-413, 223, 211)==1){
return 801;}//B0 decays to D*- omega pi+ 

if ( PcheckDecay(genpart,-413, 2212, -2212, 211)==1){
return 802;}//B0 decays to D*- p+ anti-p- pi+ 

if ( PcheckDecay(genpart,-413, 2212, -2112)==1){
return 803;}//B0 decays to D*- p+ anti-n0 

if ( PcheckDecay(genpart,-10413, 211)==1){
return 804;}//B0 decays to D_1- pi+ 

if ( PcheckDecay(genpart,-415, 211)==1){
return 805;}//B0 decays to D_2*- pi+ 

if ( PcheckDecay(genpart,-10413, 213)==1){
return 806;}//B0 decays to D_1- rho+ 

if ( PcheckDecay(genpart,-415, 213)==1){
return 807;}//B0 decays to D_2*- rho+ 

if ( PcheckDecay(genpart,-425, 111)==1){
return 808;}//B0 decays to anti-D_2*0 pi0 

if ( PcheckDecay(genpart,-10423, 111)==1){
return 809;}//B0 decays to anti-D_10 pi0 

if ( PcheckDecay(genpart,-20423, 111)==1){
return 810;}//B0 decays to anti-D'_10 pi0 

if ( PcheckDecay(genpart,-10421, 111)==1){
return 811;}//B0 decays to anti-D_0*0 pi0 

if ( PcheckDecay(genpart,-425, 221)==1){
return 812;}//B0 decays to anti-D_2*0 eta 

if ( PcheckDecay(genpart,-10423, 221)==1){
return 813;}//B0 decays to anti-D_10 eta 

if ( PcheckDecay(genpart,-20423, 221)==1){
return 814;}//B0 decays to anti-D'_10 eta 

if ( PcheckDecay(genpart,-10421, 221)==1){
return 815;}//B0 decays to anti-D_0*0 eta 

if ( PcheckDecay(genpart,-425, 331)==1){
return 816;}//B0 decays to anti-D_2*0 eta' 

if ( PcheckDecay(genpart,-10423, 331)==1){
return 817;}//B0 decays to anti-D_10 eta' 

if ( PcheckDecay(genpart,-20423, 331)==1){
return 818;}//B0 decays to anti-D'_10 eta' 

if ( PcheckDecay(genpart,-10421, 331)==1){
return 819;}//B0 decays to anti-D_0*0 eta' 

if ( PcheckDecay(genpart,-425, 223)==1){
return 820;}//B0 decays to anti-D_2*0 omega 

if ( PcheckDecay(genpart,-10423, 223)==1){
return 821;}//B0 decays to anti-D_10 omega 

if ( PcheckDecay(genpart,-20423, 223)==1){
return 822;}//B0 decays to anti-D'_10 omega 

if ( PcheckDecay(genpart,-10421, 223)==1){
return 823;}//B0 decays to anti-D_0*0 omega 

if ( PcheckDecay(genpart,-425, 113)==1){
return 824;}//B0 decays to anti-D_2*0 rho0 

if ( PcheckDecay(genpart,-10423, 113)==1){
return 825;}//B0 decays to anti-D_10 rho0 

if ( PcheckDecay(genpart,-20423, 113)==1){
return 826;}//B0 decays to anti-D'_10 rho0 

if ( PcheckDecay(genpart,-10421, 113)==1){
return 827;}//B0 decays to anti-D_0*0 rho0 

if ( PcheckDecay(genpart,-421, 111)==1){
return 828;}//B0 decays to anti-D0 pi0 

if ( PcheckDecay(genpart,-421, -211, 211)==1){
return 829;}//B0 decays to anti-D0 pi- pi+ 

if ( PcheckDecay(genpart,113, -421)==1){
return 830;}//B0 decays to rho0 anti-D0 

if ( PcheckDecay(genpart,-421, 111, 111)==1){
return 831;}//B0 decays to anti-D0 pi0 pi0 

if ( PcheckDecay(genpart,-421, 221)==1){
return 832;}//B0 decays to anti-D0 eta 

if ( PcheckDecay(genpart,-421, 331)==1){
return 833;}//B0 decays to anti-D0 eta' 

if ( PcheckDecay(genpart,223, -421)==1){
return 834;}//B0 decays to omega anti-D0 

if ( PcheckDecay(genpart,-421, -211, 211, 111)==1){
return 835;}//B0 decays to anti-D0 pi- pi+ pi0 

if ( PcheckDecay(genpart,-421, 111, 111, 111)==1){
return 836;}//B0 decays to anti-D0 pi0 pi0 pi0 

if ( PcheckDecay(genpart,-421, 2212, -2212)==1){
return 837;}//B0 decays to anti-D0 p+ anti-p- 

if ( PcheckDecay(genpart,-423, 111)==1){
return 838;}//B0 decays to anti-D*0 pi0 

if ( PcheckDecay(genpart,-423, -211, 211)==1){
return 839;}//B0 decays to anti-D*0 pi- pi+ 

if ( PcheckDecay(genpart,-423, 113)==1){
return 840;}//B0 decays to anti-D*0 rho0 

if ( PcheckDecay(genpart,-423, 111, 111)==1){
return 841;}//B0 decays to anti-D*0 pi0 pi0 

if ( PcheckDecay(genpart,-423, 221)==1){
return 842;}//B0 decays to anti-D*0 eta 

if ( PcheckDecay(genpart,-423, 331)==1){
return 843;}//B0 decays to anti-D*0 eta' 

if ( PcheckDecay(genpart,-423, 223)==1){
return 844;}//B0 decays to anti-D*0 omega 

if ( PcheckDecay(genpart,-423, -211, 211, 111)==1){
return 845;}//B0 decays to anti-D*0 pi- pi+ pi0 

if ( PcheckDecay(genpart,-423, 111, 111, 111)==1){
return 846;}//B0 decays to anti-D*0 pi0 pi0 pi0 

if ( PcheckDecay(genpart,-423, 2212, -2212)==1){
return 847;}//B0 decays to anti-D*0 p+ anti-p- 

if ( PcheckDecay(genpart,-411, 321)==1){
return 848;}//B0 decays to D- K+ 

if ( PcheckDecay(genpart,323, -411)==1){
return 849;}//B0 decays to K*+ D- 

if ( PcheckDecay(genpart,-413, 321)==1){
return 850;}//B0 decays to D*- K+ 

if ( PcheckDecay(genpart,-413, 323)==1){
return 851;}//B0 decays to D*- K*+ 

if ( PcheckDecay(genpart,-421, 311)==1){
return 852;}//B0 decays to anti-D0 K0 

if ( PcheckDecay(genpart,313, -421)==1){
return 853;}//B0 decays to K*0 anti-D0 

if ( PcheckDecay(genpart,-423, 311)==1){
return 854;}//B0 decays to anti-D*0 K0 

if ( PcheckDecay(genpart,-423, 313)==1){
return 855;}//B0 decays to anti-D*0 K*0 

if ( PcheckDecay(genpart,2, -3, -4, 1)==1){
return 856;}//B0 decays to u anti-s anti-c d 

if ( PcheckDecay(genpart,2, -1, -4, 1)==1){
return 857;}//B0 decays to u anti-d anti-c d 

if ( PcheckDecay(genpart,-4122, 2212)==1){
return 858;}//B0 decays to anti-Lambda_c- p+ 

if ( PcheckDecay(genpart,-4212, 2212)==1){
return 859;}//B0 decays to anti-Sigma_c- p+ 

if ( PcheckDecay(genpart,-4112, 2112)==1){
return 860;}//B0 decays to anti-Sigma_c0 n0 

if ( PcheckDecay(genpart,-4101, 2101)==1){
return 861;}//B0 decays to anti-cd_0 ud_0 

if ( PcheckDecay(genpart,-4103, 2103)==1){
return 862;}//B0 decays to anti-cd_1 ud_1 

if ( PcheckDecay(genpart,-4101, 2101)==1){
return 863;}//B0 decays to anti-cd_0 ud_0 

if ( PcheckDecay(genpart,-4103, 2103)==1){
return 864;}//B0 decays to anti-cd_1 ud_1 

if ( PcheckDecay(genpart,-4232, 4122)==1){
return 865;}//B0 decays to anti-Xi_c- Lambda_c+ 

if ( PcheckDecay(genpart,-4232, 4212)==1){
return 866;}//B0 decays to anti-Xi_c- Sigma_c+ 

if ( PcheckDecay(genpart,-4132, 4112)==1){
return 867;}//B0 decays to anti-Xi_c0 Sigma_c0 

if ( PcheckDecay(genpart,-4301, 4101)==1){
return 868;}//B0 decays to anti-cs_0 cd_0 

if ( PcheckDecay(genpart,-4232, 2212)==1){
return 869;}//B0 decays to anti-Xi_c- p+ 

if ( PcheckDecay(genpart,-4132, 2112)==1){
return 870;}//B0 decays to anti-Xi_c0 n0 

if ( PcheckDecay(genpart,-4301, 2101)==1){
return 871;}//B0 decays to anti-cs_0 ud_0 

if ( PcheckDecay(genpart,-4303, 2103)==1){
return 872;}//B0 decays to anti-cs_1 ud_1 

if ( PcheckDecay(genpart,-4301, 2101)==1){
return 873;}//B0 decays to anti-cs_0 ud_0 

if ( PcheckDecay(genpart,-4303, 2103)==1){
return 874;}//B0 decays to anti-cs_1 ud_1 

if ( PcheckDecay(genpart,-4122, 4122)==1){
return 875;}//B0 decays to anti-Lambda_c- Lambda_c+ 

if ( PcheckDecay(genpart,-4212, 4122)==1){
return 876;}//B0 decays to anti-Sigma_c- Lambda_c+ 

if ( PcheckDecay(genpart,-4122, 4212)==1){
return 877;}//B0 decays to anti-Lambda_c- Sigma_c+ 

if ( PcheckDecay(genpart,-4212, 4212)==1){
return 878;}//B0 decays to anti-Sigma_c- Sigma_c+ 

if ( PcheckDecay(genpart,-4112, 4112)==1){
return 879;}//B0 decays to anti-Sigma_c0 Sigma_c0 

if ( PcheckDecay(genpart,-4101, 4101)==1){
return 880;}//B0 decays to anti-cd_0 cd_0 
     return -200;
 }


 else if (AmId(genpart)==521){

if ( PcheckDecay(genpart,-423, -11, 12)==1){
return 881;}//B+ decays to anti-D*0 e+ nu_e 

if ( PcheckDecay(genpart,-421, -11, 12)==1){
return 882;}//B+ decays to anti-D0 e+ nu_e 

if ( PcheckDecay(genpart,-10423, -11, 12)==1){
return 883;}//B+ decays to anti-D_10 e+ nu_e 

if ( PcheckDecay(genpart,-10421, -11, 12)==1){
return 884;}//B+ decays to anti-D_0*0 e+ nu_e 

if ( PcheckDecay(genpart,-20423, -11, 12)==1){
return 885;}//B+ decays to anti-D'_10 e+ nu_e 

if ( PcheckDecay(genpart,-425, -11, 12)==1){
return 886;}//B+ decays to anti-D_2*0 e+ nu_e 

if ( PcheckDecay(genpart,-100421, -11, 12)==1){
return 887;}//B+ decays to anti-D(2S)0 e+ nu_e 

if ( PcheckDecay(genpart,-100423, -11, 12)==1){
return 888;}//B+ decays to anti-D*(2S)0 e+ nu_e 

if ( PcheckDecay(genpart,-423, 111, -11, 12)==1){
return 889;}//B+ decays to anti-D*0 pi0 e+ nu_e 

if ( PcheckDecay(genpart,-413, 211, -11, 12)==1){
return 890;}//B+ decays to D*- pi+ e+ nu_e 

if ( PcheckDecay(genpart,-421, 111, -11, 12)==1){
return 891;}//B+ decays to anti-D0 pi0 e+ nu_e 

if ( PcheckDecay(genpart,-411, 211, -11, 12)==1){
return 892;}//B+ decays to D- pi+ e+ nu_e 

if ( PcheckDecay(genpart,-423, -13, 14)==1){
return 893;}//B+ decays to anti-D*0 mu+ nu_mu 

if ( PcheckDecay(genpart,-421, -13, 14)==1){
return 894;}//B+ decays to anti-D0 mu+ nu_mu 

if ( PcheckDecay(genpart,-10423, -13, 14)==1){
return 895;}//B+ decays to anti-D_10 mu+ nu_mu 

if ( PcheckDecay(genpart,-10421, -13, 14)==1){
return 896;}//B+ decays to anti-D_0*0 mu+ nu_mu 

if ( PcheckDecay(genpart,-20423, -13, 14)==1){
return 897;}//B+ decays to anti-D'_10 mu+ nu_mu 

if ( PcheckDecay(genpart,-425, -13, 14)==1){
return 898;}//B+ decays to anti-D_2*0 mu+ nu_mu 

if ( PcheckDecay(genpart,-100421, -13, 14)==1){
return 899;}//B+ decays to anti-D(2S)0 mu+ nu_mu 

if ( PcheckDecay(genpart,-100423, -13, 14)==1){
return 900;}//B+ decays to anti-D*(2S)0 mu+ nu_mu 

if ( PcheckDecay(genpart,-423, 111, -13, 14)==1){
return 901;}//B+ decays to anti-D*0 pi0 mu+ nu_mu 

if ( PcheckDecay(genpart,-413, 211, -13, 14)==1){
return 902;}//B+ decays to D*- pi+ mu+ nu_mu 

if ( PcheckDecay(genpart,-421, 111, -13, 14)==1){
return 903;}//B+ decays to anti-D0 pi0 mu+ nu_mu 

if ( PcheckDecay(genpart,-411, 211, -13, 14)==1){
return 904;}//B+ decays to D- pi+ mu+ nu_mu 

if ( PcheckDecay(genpart,-423, -15, 16)==1){
return 905;}//B+ decays to anti-D*0 tau+ nu_tau 

if ( PcheckDecay(genpart,-421, -15, 16)==1){
return 906;}//B+ decays to anti-D0 tau+ nu_tau 

if ( PcheckDecay(genpart,-10423, -15, 16)==1){
return 907;}//B+ decays to anti-D_10 tau+ nu_tau 

if ( PcheckDecay(genpart,-10421, -15, 16)==1){
return 908;}//B+ decays to anti-D_0*0 tau+ nu_tau 

if ( PcheckDecay(genpart,-20423, -15, 16)==1){
return 909;}//B+ decays to anti-D'_10 tau+ nu_tau 

if ( PcheckDecay(genpart,-425, -15, 16)==1){
return 910;}//B+ decays to anti-D_2*0 tau+ nu_tau 

if ( PcheckDecay(genpart,441, 321)==1){
return 911;}//B+ decays to eta_c K+ 

if ( PcheckDecay(genpart,323, 441)==1){
return 912;}//B+ decays to K*+ eta_c 

if ( PcheckDecay(genpart,441, 311, 211)==1){
return 913;}//B+ decays to eta_c K0 pi+ 

if ( PcheckDecay(genpart,441, 321, 111)==1){
return 914;}//B+ decays to eta_c K+ pi0 

if ( PcheckDecay(genpart,441, 321, 211, -211)==1){
return 915;}//B+ decays to eta_c K+ pi+ pi- 

if ( PcheckDecay(genpart,441, 321, 111, 111)==1){
return 916;}//B+ decays to eta_c K+ pi0 pi0 

if ( PcheckDecay(genpart,441, 311, 211, 111)==1){
return 917;}//B+ decays to eta_c K0 pi+ pi0 

if ( PcheckDecay(genpart,441, -3122, 2212)==1){
return 918;}//B+ decays to eta_c anti-Lambda0 p+ 

if ( PcheckDecay(genpart,441, 3222, 2112)==1){
return 919;}//B+ decays to eta_c Sigma+ n0 

if ( PcheckDecay(genpart,441, 30353)==1){
return 920;}//B+ decays to eta_c Xsu 

if ( PcheckDecay(genpart,100441, 321)==1){
return 921;}//B+ decays to eta_c(2S) K+ 

if ( PcheckDecay(genpart,323, 100441)==1){
return 922;}//B+ decays to K*+ eta_c(2S) 

if ( PcheckDecay(genpart,100441, 311, 211)==1){
return 923;}//B+ decays to eta_c(2S) K0 pi+ 

if ( PcheckDecay(genpart,100441, 321, 111)==1){
return 924;}//B+ decays to eta_c(2S) K+ pi0 

if ( PcheckDecay(genpart,100441, 321, 211, -211)==1){
return 925;}//B+ decays to eta_c(2S) K+ pi+ pi- 

if ( PcheckDecay(genpart,100441, 321, 111, 111)==1){
return 926;}//B+ decays to eta_c(2S) K+ pi0 pi0 

if ( PcheckDecay(genpart,100441, 311, 211, 111)==1){
return 927;}//B+ decays to eta_c(2S) K0 pi+ pi0 

if ( PcheckDecay(genpart,100441, 30353)==1){
return 928;}//B+ decays to eta_c(2S) Xsu 

if ( PcheckDecay(genpart,443, 321)==1){
return 929;}//B+ decays to psi K+ 

if ( PcheckDecay(genpart,443, 323)==1){
return 930;}//B+ decays to psi K*+ 

if ( PcheckDecay(genpart,443, 10323)==1){
return 931;}//B+ decays to psi K_1+ 

if ( PcheckDecay(genpart,443, 20323)==1){
return 932;}//B+ decays to psi K'_1+ 

if ( PcheckDecay(genpart,443, 321, 211, -211)==1){
return 933;}//B+ decays to psi K+ pi+ pi- 

if ( PcheckDecay(genpart,443, 333, 321)==1){
return 934;}//B+ decays to psi phi K+ 

if ( PcheckDecay(genpart,443, 211)==1){
return 935;}//B+ decays to psi pi+ 

if ( PcheckDecay(genpart,443, 213)==1){
return 936;}//B+ decays to psi rho+ 

if ( PcheckDecay(genpart,443, 20213)==1){
return 937;}//B+ decays to psi a_1+ 

if ( PcheckDecay(genpart,443, 10321)==1){
return 938;}//B+ decays to psi K_0*+ 

if ( PcheckDecay(genpart,443, 325)==1){
return 939;}//B+ decays to psi K_2*+ 

if ( PcheckDecay(genpart,443, 30323)==1){
return 940;}//B+ decays to psi K''*+ 

if ( PcheckDecay(genpart,443, 311, 211)==1){
return 941;}//B+ decays to psi K0 pi+ 

if ( PcheckDecay(genpart,443, 321, 111)==1){
return 942;}//B+ decays to psi K+ pi0 

if ( PcheckDecay(genpart,443, 321, 111, 111)==1){
return 943;}//B+ decays to psi K+ pi0 pi0 

if ( PcheckDecay(genpart,443, 311, 211, 111)==1){
return 944;}//B+ decays to psi K0 pi+ pi0 

if ( PcheckDecay(genpart,443, -3122, 2212)==1){
return 945;}//B+ decays to psi anti-Lambda0 p+ 

if ( PcheckDecay(genpart,443, 3222, 2112)==1){
return 946;}//B+ decays to psi Sigma+ n0 

if ( PcheckDecay(genpart,443, 30353)==1){
return 947;}//B+ decays to psi Xsu 

if ( PcheckDecay(genpart,100443, 321)==1){
return 948;}//B+ decays to psi(2S) K+ 

if ( PcheckDecay(genpart,100443, 323)==1){
return 949;}//B+ decays to psi(2S) K*+ 

if ( PcheckDecay(genpart,100443, 321, 211, -211)==1){
return 950;}//B+ decays to psi(2S) K+ pi+ pi- 

if ( PcheckDecay(genpart,10441, 321)==1){
return 951;}//B+ decays to chi_c0 K+ 

if ( PcheckDecay(genpart,323, 10441)==1){
return 952;}//B+ decays to K*+ chi_c0 

if ( PcheckDecay(genpart,10441, 30353)==1){
return 953;}//B+ decays to chi_c0 Xsu 

if ( PcheckDecay(genpart,20443, 321)==1){
return 954;}//B+ decays to chi_c1 K+ 

if ( PcheckDecay(genpart,20443, 323)==1){
return 955;}//B+ decays to chi_c1 K*+ 

if ( PcheckDecay(genpart,20443, 211)==1){
return 956;}//B+ decays to chi_c1 pi+ 

if ( PcheckDecay(genpart,20443, 30353)==1){
return 957;}//B+ decays to chi_c1 Xsu 

if ( PcheckDecay(genpart,445, 321)==1){
return 958;}//B+ decays to chi_c2 K+ 

if ( PcheckDecay(genpart,445, 323)==1){
return 959;}//B+ decays to chi_c2 K*+ 

if ( PcheckDecay(genpart,445, 30353)==1){
return 960;}//B+ decays to chi_c2 Xsu 

if ( PcheckDecay(genpart,10443, 321)==1){
return 961;}//B+ decays to h_c K+ 

if ( PcheckDecay(genpart,10443, 323)==1){
return 962;}//B+ decays to h_c K*+ 

if ( PcheckDecay(genpart,10443, 30353)==1){
return 963;}//B+ decays to h_c Xsu 

if ( PcheckDecay(genpart,120443, 321)==1){
return 964;}//B+ decays to X(3872) K+ 

if ( PcheckDecay(genpart,90000443, 321)==1){
return 965;}//B+ decays to Y(3940) K+ 

if ( PcheckDecay(genpart,-421, 431)==1){
return 966;}//B+ decays to anti-D0 D_s+ 

if ( PcheckDecay(genpart,-423, 431)==1){
return 967;}//B+ decays to anti-D*0 D_s+ 

if ( PcheckDecay(genpart,433, -421)==1){
return 968;}//B+ decays to D_s*+ anti-D0 

if ( PcheckDecay(genpart,433, -423)==1){
return 969;}//B+ decays to D_s*+ anti-D*0 

if ( PcheckDecay(genpart,431, -10421)==1){
return 970;}//B+ decays to D_s+ anti-D_0*0 

if ( PcheckDecay(genpart,433, -10421)==1){
return 971;}//B+ decays to D_s*+ anti-D_0*0 

if ( PcheckDecay(genpart,-20423, 431)==1){
return 972;}//B+ decays to anti-D'_10 D_s+ 

if ( PcheckDecay(genpart,-20423, 433)==1){
return 973;}//B+ decays to anti-D'_10 D_s*+ 

if ( PcheckDecay(genpart,-10423, 431)==1){
return 974;}//B+ decays to anti-D_10 D_s+ 

if ( PcheckDecay(genpart,-10423, 433)==1){
return 975;}//B+ decays to anti-D_10 D_s*+ 

if ( PcheckDecay(genpart,-425, 431)==1){
return 976;}//B+ decays to anti-D_2*0 D_s+ 

if ( PcheckDecay(genpart,-425, 433)==1){
return 977;}//B+ decays to anti-D_2*0 D_s*+ 

if ( PcheckDecay(genpart,-421, 10431)==1){
return 978;}//B+ decays to anti-D0 D_s0*+ 

if ( PcheckDecay(genpart,-423, 10431)==1){
return 979;}//B+ decays to anti-D*0 D_s0*+ 

if ( PcheckDecay(genpart,20433, -421)==1){
return 980;}//B+ decays to D'_s1+ anti-D0 

if ( PcheckDecay(genpart,20433, -423)==1){
return 981;}//B+ decays to D'_s1+ anti-D*0 

if ( PcheckDecay(genpart,431, -411, 211)==1){
return 982;}//B+ decays to D_s+ D- pi+ 

if ( PcheckDecay(genpart,431, -421, 111)==1){
return 983;}//B+ decays to D_s+ anti-D0 pi0 

if ( PcheckDecay(genpart,433, -411, 211)==1){
return 984;}//B+ decays to D_s*+ D- pi+ 

if ( PcheckDecay(genpart,433, -421, 111)==1){
return 985;}//B+ decays to D_s*+ anti-D0 pi0 

if ( PcheckDecay(genpart,431, -411, 211, 111)==1){
return 986;}//B+ decays to D_s+ D- pi+ pi0 

if ( PcheckDecay(genpart,431, -421, 211, -211)==1){
return 987;}//B+ decays to D_s+ anti-D0 pi+ pi- 

if ( PcheckDecay(genpart,431, -421, 111, 111)==1){
return 988;}//B+ decays to D_s+ anti-D0 pi0 pi0 

if ( PcheckDecay(genpart,433, -411, 211, 111)==1){
return 989;}//B+ decays to D_s*+ D- pi+ pi0 

if ( PcheckDecay(genpart,433, -421, 211, -211)==1){
return 990;}//B+ decays to D_s*+ anti-D0 pi+ pi- 

if ( PcheckDecay(genpart,433, -421, 111, 111)==1){
return 991;}//B+ decays to D_s*+ anti-D0 pi0 pi0 

if ( PcheckDecay(genpart,-431, 211, 321)==1){
return 992;}//B+ decays to D_s- pi+ K+ 

if ( PcheckDecay(genpart,-433, 211, 321)==1){
return 993;}//B+ decays to D_s*- pi+ K+ 

if ( PcheckDecay(genpart,-431, 211, 321, 111)==1){
return 994;}//B+ decays to D_s- pi+ K+ pi0 

if ( PcheckDecay(genpart,-433, 211, 321, 111)==1){
return 995;}//B+ decays to D_s*- pi+ K+ pi0 

if ( PcheckDecay(genpart,-431, 211, 311, 211)==1){
return 996;}//B+ decays to D_s- pi+ K0 pi+ 

if ( PcheckDecay(genpart,-433, 211, 311, 211)==1){
return 997;}//B+ decays to D_s*- pi+ K0 pi+ 

if ( PcheckDecay(genpart,-431, 211, 321, 111, 111)==1){
return 998;}//B+ decays to D_s- pi+ K+ pi0 pi0 

if ( PcheckDecay(genpart,-433, 211, 321, 111, 111)==1){
return 999;}//B+ decays to D_s*- pi+ K+ pi0 pi0 

if ( PcheckDecay(genpart,-431, 211, 321, 211, -211)==1){
return 1000;}//B+ decays to D_s- pi+ K+ pi+ pi- 

if ( PcheckDecay(genpart,-433, 211, 321, 211, -211)==1){
return 1001;}//B+ decays to D_s*- pi+ K+ pi+ pi- 

if ( PcheckDecay(genpart,-431, 211, 311, 211, 111)==1){
return 1002;}//B+ decays to D_s- pi+ K0 pi+ pi0 

if ( PcheckDecay(genpart,-433, 211, 311, 211, 111)==1){
return 1003;}//B+ decays to D_s*- pi+ K0 pi+ pi0 

if ( PcheckDecay(genpart,411, -421)==1){
return 1004;}//B+ decays to D+ anti-D0 

if ( PcheckDecay(genpart,-423, 411)==1){
return 1005;}//B+ decays to anti-D*0 D+ 

if ( PcheckDecay(genpart,413, -421)==1){
return 1006;}//B+ decays to D*+ anti-D0 

if ( PcheckDecay(genpart,413, -423)==1){
return 1007;}//B+ decays to D*+ anti-D*0 

//if ( PcheckDecay(genpart,return 1008;}//B+ decays to 

if ( PcheckDecay(genpart,20413, -421)==1){
return 1009;}//B+ decays to D'_1+ anti-D0 

if ( PcheckDecay(genpart,-20423, 411)==1){
return 1010;}//B+ decays to anti-D'_10 D+ 

if ( PcheckDecay(genpart,20413, -423)==1){
return 1011;}//B+ decays to D'_1+ anti-D*0 

if ( PcheckDecay(genpart,-20423, 413)==1){
return 1012;}//B+ decays to anti-D'_10 D*+ 

if ( PcheckDecay(genpart,10413, -421)==1){
return 1013;}//B+ decays to D_1+ anti-D0 

if ( PcheckDecay(genpart,-10423, 411)==1){
return 1014;}//B+ decays to anti-D_10 D+ 

if ( PcheckDecay(genpart,10413, -423)==1){
return 1015;}//B+ decays to D_1+ anti-D*0 

if ( PcheckDecay(genpart,-10423, 413)==1){
return 1016;}//B+ decays to anti-D_10 D*+ 

if ( PcheckDecay(genpart,415, -421)==1){
return 1017;}//B+ decays to D_2*+ anti-D0 

if ( PcheckDecay(genpart,-425, 411)==1){
return 1018;}//B+ decays to anti-D_2*0 D+ 

if ( PcheckDecay(genpart,415, -423)==1){
return 1019;}//B+ decays to D_2*+ anti-D*0 

if ( PcheckDecay(genpart,-425, 413)==1){
return 1020;}//B+ decays to anti-D_2*0 D*+ 

//if ( PcheckDecay(genpart,return 1021;}//B+ decays to 

if ( PcheckDecay(genpart,10433, -421)==1){
return 1022;}//B+ decays to D_s1+ anti-D0 

if ( PcheckDecay(genpart,10433, -423)==1){
return 1023;}//B+ decays to D_s1+ anti-D*0 

if ( PcheckDecay(genpart,435, -421)==1){
return 1024;}//B+ decays to D_s2*+ anti-D0 

if ( PcheckDecay(genpart,435, -423)==1){
return 1025;}//B+ decays to D_s2*+ anti-D*0 

if ( PcheckDecay(genpart,9000433, -421)==1){
return 1026;}//B+ decays to D_sj(2700)+ anti-D0 

if ( PcheckDecay(genpart,9000433, -423)==1){
return 1027;}//B+ decays to D_sj(2700)+ anti-D*0 

if ( PcheckDecay(genpart,30443, 321)==1){
return 1028;}//B+ decays to psi(3770) K+ 

if ( PcheckDecay(genpart,30443, 323)==1){
return 1029;}//B+ decays to psi(3770) K*+ 

if ( PcheckDecay(genpart,30443, 311, 211)==1){
return 1030;}//B+ decays to psi(3770) K0 pi+ 

if ( PcheckDecay(genpart,30443, 321, 111)==1){
return 1031;}//B+ decays to psi(3770) K+ pi0 

if ( PcheckDecay(genpart,30443, 321, 111, 111)==1){
return 1032;}//B+ decays to psi(3770) K+ pi0 pi0 

if ( PcheckDecay(genpart,30443, 311, 211, 111)==1){
return 1033;}//B+ decays to psi(3770) K0 pi+ pi0 

if ( PcheckDecay(genpart,9000443, 321)==1){
return 1034;}//B+ decays to psi(4040) K+ 

if ( PcheckDecay(genpart,9000443, 323)==1){
return 1035;}//B+ decays to psi(4040) K*+ 

if ( PcheckDecay(genpart,9000443, 311, 211)==1){
return 1036;}//B+ decays to psi(4040) K0 pi+ 

if ( PcheckDecay(genpart,9000443, 321, 111)==1){
return 1037;}//B+ decays to psi(4040) K+ pi0 

if ( PcheckDecay(genpart,9010443, 321)==1){
return 1038;}//B+ decays to psi(4160) K+ 

if ( PcheckDecay(genpart,9010443, 323)==1){
return 1039;}//B+ decays to psi(4160) K*+ 

if ( PcheckDecay(genpart,9010443, 311, 211)==1){
return 1040;}//B+ decays to psi(4160) K0 pi+ 

if ( PcheckDecay(genpart,9010443, 321, 111)==1){
return 1041;}//B+ decays to psi(4160) K+ pi0 

//if ( PcheckDecay(genpart,return 1042;}//B+ decays to 

if ( PcheckDecay(genpart,411, -421, 111)==1){
return 1043;}//B+ decays to D+ anti-D0 pi0 

if ( PcheckDecay(genpart,411, -411, 211)==1){
return 1044;}//B+ decays to D+ D- pi+ 

if ( PcheckDecay(genpart,421, -421, 211)==1){
return 1045;}//B+ decays to D0 anti-D0 pi+ 

if ( PcheckDecay(genpart,413, -421, 111)==1){
return 1046;}//B+ decays to D*+ anti-D0 pi0 

if ( PcheckDecay(genpart,413, -411, 211)==1){
return 1047;}//B+ decays to D*+ D- pi+ 

if ( PcheckDecay(genpart,423, -421, 211)==1){
return 1048;}//B+ decays to D*0 anti-D0 pi+ 

if ( PcheckDecay(genpart,411, -423, 111)==1){
return 1049;}//B+ decays to D+ anti-D*0 pi0 

if ( PcheckDecay(genpart,411, -413, 211)==1){
return 1050;}//B+ decays to D+ D*- pi+ 

if ( PcheckDecay(genpart,421, -423, 211)==1){
return 1051;}//B+ decays to D0 anti-D*0 pi+ 

if ( PcheckDecay(genpart,413, -423, 111)==1){
return 1052;}//B+ decays to D*+ anti-D*0 pi0 

if ( PcheckDecay(genpart,413, -413, 211)==1){
return 1053;}//B+ decays to D*+ D*- pi+ 

if ( PcheckDecay(genpart,423, -423, 211)==1){
return 1054;}//B+ decays to D*0 anti-D*0 pi+ 

if ( PcheckDecay(genpart,411, -421, 113)==1){
return 1055;}//B+ decays to D+ anti-D0 rho0 

if ( PcheckDecay(genpart,411, -411, 213)==1){
return 1056;}//B+ decays to D+ D- rho+ 

if ( PcheckDecay(genpart,421, -421, 213)==1){
return 1057;}//B+ decays to D0 anti-D0 rho+ 

if ( PcheckDecay(genpart,413, -421, 113)==1){
return 1058;}//B+ decays to D*+ anti-D0 rho0 

if ( PcheckDecay(genpart,413, -411, 213)==1){
return 1059;}//B+ decays to D*+ D- rho+ 

if ( PcheckDecay(genpart,423, -421, 213)==1){
return 1060;}//B+ decays to D*0 anti-D0 rho+ 

if ( PcheckDecay(genpart,411, -423, 113)==1){
return 1061;}//B+ decays to D+ anti-D*0 rho0 

if ( PcheckDecay(genpart,411, -413, 213)==1){
return 1062;}//B+ decays to D+ D*- rho+ 

if ( PcheckDecay(genpart,421, -423, 213)==1){
return 1063;}//B+ decays to D0 anti-D*0 rho+ 

if ( PcheckDecay(genpart,413, -423, 113)==1){
return 1064;}//B+ decays to D*+ anti-D*0 rho0 

if ( PcheckDecay(genpart,413, -413, 213)==1){
return 1065;}//B+ decays to D*+ D*- rho+ 

if ( PcheckDecay(genpart,423, -423, 213)==1){
return 1066;}//B+ decays to D*0 anti-D*0 rho+ 

if ( PcheckDecay(genpart,411, -421, 111, 111)==1){
return 1067;}//B+ decays to D+ anti-D0 pi0 pi0 

if ( PcheckDecay(genpart,411, -421, 211, -211)==1){
return 1068;}//B+ decays to D+ anti-D0 pi+ pi- 

if ( PcheckDecay(genpart,411, -411, 211, 111)==1){
return 1069;}//B+ decays to D+ D- pi+ pi0 

if ( PcheckDecay(genpart,421, -421, 211, 111)==1){
return 1070;}//B+ decays to D0 anti-D0 pi+ pi0 

if ( PcheckDecay(genpart,421, -411, 211, 211)==1){
return 1071;}//B+ decays to D0 D- pi+ pi+ 

if ( PcheckDecay(genpart,413, -421, 111, 111)==1){
return 1072;}//B+ decays to D*+ anti-D0 pi0 pi0 

if ( PcheckDecay(genpart,413, -421, 211, -211)==1){
return 1073;}//B+ decays to D*+ anti-D0 pi+ pi- 

if ( PcheckDecay(genpart,413, -411, 211, 111)==1){
return 1074;}//B+ decays to D*+ D- pi+ pi0 

if ( PcheckDecay(genpart,423, -421, 211, 111)==1){
return 1075;}//B+ decays to D*0 anti-D0 pi+ pi0 

if ( PcheckDecay(genpart,423, -411, 211, 211)==1){
return 1076;}//B+ decays to D*0 D- pi+ pi+ 

if ( PcheckDecay(genpart,411, -423, 111, 111)==1){
return 1077;}//B+ decays to D+ anti-D*0 pi0 pi0 

if ( PcheckDecay(genpart,411, -423, 211, -211)==1){
return 1078;}//B+ decays to D+ anti-D*0 pi+ pi- 

if ( PcheckDecay(genpart,411, -413, 211, 111)==1){
return 1079;}//B+ decays to D+ D*- pi+ pi0 

if ( PcheckDecay(genpart,421, -423, 211, 111)==1){
return 1080;}//B+ decays to D0 anti-D*0 pi+ pi0 

if ( PcheckDecay(genpart,421, -413, 211, 211)==1){
return 1081;}//B+ decays to D0 D*- pi+ pi+ 

if ( PcheckDecay(genpart,413, -423, 111, 111)==1){
return 1082;}//B+ decays to D*+ anti-D*0 pi0 pi0 

if ( PcheckDecay(genpart,413, -423, 211, -211)==1){
return 1083;}//B+ decays to D*+ anti-D*0 pi+ pi- 

if ( PcheckDecay(genpart,413, -413, 211, 111)==1){
return 1084;}//B+ decays to D*+ D*- pi+ pi0 

if ( PcheckDecay(genpart,423, -423, 211, 111)==1){
return 1085;}//B+ decays to D*0 anti-D*0 pi+ pi0 

if ( PcheckDecay(genpart,423, -413, 211, 211)==1){
return 1086;}//B+ decays to D*0 D*- pi+ pi+ 

if ( PcheckDecay(genpart,-421, 421, 321)==1){
return 1087;}//B+ decays to anti-D0 D0 K+ 

if ( PcheckDecay(genpart,-421, 411, 311)==1){
return 1088;}//B+ decays to anti-D0 D+ K0 

if ( PcheckDecay(genpart,-421, 421, 323)==1){
return 1089;}//B+ decays to anti-D0 D0 K*+ 

if ( PcheckDecay(genpart,-421, 411, 313)==1){
return 1090;}//B+ decays to anti-D0 D+ K*0 

if ( PcheckDecay(genpart,-423, 421, 321)==1){
return 1091;}//B+ decays to anti-D*0 D0 K+ 

if ( PcheckDecay(genpart,-423, 411, 311)==1){
return 1092;}//B+ decays to anti-D*0 D+ K0 

if ( PcheckDecay(genpart,-421, 423, 321)==1){
return 1093;}//B+ decays to anti-D0 D*0 K+ 

if ( PcheckDecay(genpart,-421, 413, 311)==1){
return 1094;}//B+ decays to anti-D0 D*+ K0 

if ( PcheckDecay(genpart,-423, 421, 323)==1){
return 1095;}//B+ decays to anti-D*0 D0 K*+ 

if ( PcheckDecay(genpart,-423, 411, 313)==1){
return 1096;}//B+ decays to anti-D*0 D+ K*0 

if ( PcheckDecay(genpart,-421, 423, 323)==1){
return 1097;}//B+ decays to anti-D0 D*0 K*+ 

if ( PcheckDecay(genpart,-421, 413, 313)==1){
return 1098;}//B+ decays to anti-D0 D*+ K*0 

if ( PcheckDecay(genpart,-423, 423, 321)==1){
return 1099;}//B+ decays to anti-D*0 D*0 K+ 

if ( PcheckDecay(genpart,-423, 413, 311)==1){
return 1100;}//B+ decays to anti-D*0 D*+ K0 

if ( PcheckDecay(genpart,-423, 423, 323)==1){
return 1101;}//B+ decays to anti-D*0 D*0 K*+ 

if ( PcheckDecay(genpart,-423, 413, 313)==1){
return 1102;}//B+ decays to anti-D*0 D*+ K*0 

if ( PcheckDecay(genpart,-10423, 421, 321)==1){
return 1103;}//B+ decays to anti-D_10 D0 K+ 

if ( PcheckDecay(genpart,-10423, 411, 311)==1){
return 1104;}//B+ decays to anti-D_10 D+ K0 

if ( PcheckDecay(genpart,-421, 10423, 321)==1){
return 1105;}//B+ decays to anti-D0 D_10 K+ 

if ( PcheckDecay(genpart,-421, 10413, 311)==1){
return 1106;}//B+ decays to anti-D0 D_1+ K0 

if ( PcheckDecay(genpart,-10423, 423, 321)==1){
return 1107;}//B+ decays to anti-D_10 D*0 K+ 

if ( PcheckDecay(genpart,-10423, 413, 311)==1){
return 1108;}//B+ decays to anti-D_10 D*+ K0 

if ( PcheckDecay(genpart,-423, 10423, 321)==1){
return 1109;}//B+ decays to anti-D*0 D_10 K+ 

if ( PcheckDecay(genpart,-423, 10413, 311)==1){
return 1110;}//B+ decays to anti-D*0 D_1+ K0 

if ( PcheckDecay(genpart,-425, 421, 321)==1){
return 1111;}//B+ decays to anti-D_2*0 D0 K+ 

if ( PcheckDecay(genpart,-425, 411, 311)==1){
return 1112;}//B+ decays to anti-D_2*0 D+ K0 

if ( PcheckDecay(genpart,-421, 425, 321)==1){
return 1113;}//B+ decays to anti-D0 D_2*0 K+ 

if ( PcheckDecay(genpart,-421, 415, 311)==1){
return 1114;}//B+ decays to anti-D0 D_2*+ K0 

if ( PcheckDecay(genpart,-425, 423, 321)==1){
return 1115;}//B+ decays to anti-D_2*0 D*0 K+ 

if ( PcheckDecay(genpart,-425, 413, 311)==1){
return 1116;}//B+ decays to anti-D_2*0 D*+ K0 

if ( PcheckDecay(genpart,-423, 425, 321)==1){
return 1117;}//B+ decays to anti-D*0 D_2*0 K+ 

if ( PcheckDecay(genpart,-423, 415, 311)==1){
return 1118;}//B+ decays to anti-D*0 D_2*+ K0 

if ( PcheckDecay(genpart,-421, 421, 321, 111)==1){
return 1119;}//B+ decays to anti-D0 D0 K+ pi0 

if ( PcheckDecay(genpart,-411, 421, 321, 211)==1){
return 1120;}//B+ decays to D- D0 K+ pi+ 

if ( PcheckDecay(genpart,-421, 421, 311, 211)==1){
return 1121;}//B+ decays to anti-D0 D0 K0 pi+ 

if ( PcheckDecay(genpart,-421, 411, 321, -211)==1){
return 1122;}//B+ decays to anti-D0 D+ K+ pi- 

if ( PcheckDecay(genpart,-421, 411, 311, 111)==1){
return 1123;}//B+ decays to anti-D0 D+ K0 pi0 

if ( PcheckDecay(genpart,-411, 411, 311, 211)==1){
return 1124;}//B+ decays to D- D+ K0 pi+ 

if ( PcheckDecay(genpart,-423, 421, 321, 111)==1){
return 1125;}//B+ decays to anti-D*0 D0 K+ pi0 

if ( PcheckDecay(genpart,-413, 421, 321, 211)==1){
return 1126;}//B+ decays to D*- D0 K+ pi+ 

if ( PcheckDecay(genpart,-423, 421, 311, 211)==1){
return 1127;}//B+ decays to anti-D*0 D0 K0 pi+ 

if ( PcheckDecay(genpart,-423, 411, 321, -211)==1){
return 1128;}//B+ decays to anti-D*0 D+ K+ pi- 

if ( PcheckDecay(genpart,-423, 411, 311, 111)==1){
return 1129;}//B+ decays to anti-D*0 D+ K0 pi0 

if ( PcheckDecay(genpart,-413, 411, 311, 211)==1){
return 1130;}//B+ decays to D*- D+ K0 pi+ 

if ( PcheckDecay(genpart,-421, 423, 321, 111)==1){
return 1131;}//B+ decays to anti-D0 D*0 K+ pi0 

if ( PcheckDecay(genpart,-411, 423, 321, 211)==1){
return 1132;}//B+ decays to D- D*0 K+ pi+ 

if ( PcheckDecay(genpart,-421, 423, 311, 211)==1){
return 1133;}//B+ decays to anti-D0 D*0 K0 pi+ 

if ( PcheckDecay(genpart,-421, 413, 321, -211)==1){
return 1134;}//B+ decays to anti-D0 D*+ K+ pi- 

if ( PcheckDecay(genpart,-421, 413, 311, 111)==1){
return 1135;}//B+ decays to anti-D0 D*+ K0 pi0 

if ( PcheckDecay(genpart,-411, 413, 311, 211)==1){
return 1136;}//B+ decays to D- D*+ K0 pi+ 

if ( PcheckDecay(genpart,-423, 423, 321, 111)==1){
return 1137;}//B+ decays to anti-D*0 D*0 K+ pi0 

if ( PcheckDecay(genpart,-413, 423, 321, 211)==1){
return 1138;}//B+ decays to D*- D*0 K+ pi+ 

if ( PcheckDecay(genpart,-423, 423, 311, 211)==1){
return 1139;}//B+ decays to anti-D*0 D*0 K0 pi+ 

if ( PcheckDecay(genpart,-423, 413, 321, -211)==1){
return 1140;}//B+ decays to anti-D*0 D*+ K+ pi- 

if ( PcheckDecay(genpart,-423, 413, 311, 111)==1){
return 1141;}//B+ decays to anti-D*0 D*+ K0 pi0 

if ( PcheckDecay(genpart,-413, 413, 311, 211)==1){
return 1142;}//B+ decays to D*- D*+ K0 pi+ 

if ( PcheckDecay(genpart,-421, 421, 323, 111)==1){
return 1143;}//B+ decays to anti-D0 D0 K*+ pi0 

if ( PcheckDecay(genpart,-411, 421, 323, 211)==1){
return 1144;}//B+ decays to D- D0 K*+ pi+ 

if ( PcheckDecay(genpart,-421, 421, 313, 211)==1){
return 1145;}//B+ decays to anti-D0 D0 K*0 pi+ 

if ( PcheckDecay(genpart,-421, 411, 323, -211)==1){
return 1146;}//B+ decays to anti-D0 D+ K*+ pi- 

if ( PcheckDecay(genpart,-421, 411, 313, 111)==1){
return 1147;}//B+ decays to anti-D0 D+ K*0 pi0 

if ( PcheckDecay(genpart,-411, 411, 313, 211)==1){
return 1148;}//B+ decays to D- D+ K*0 pi+ 

if ( PcheckDecay(genpart,-421, 421, 321, 113)==1){
return 1149;}//B+ decays to anti-D0 D0 K+ rho0 

if ( PcheckDecay(genpart,-411, 421, 321, 213)==1){
return 1150;}//B+ decays to D- D0 K+ rho+ 

if ( PcheckDecay(genpart,-421, 421, 311, 213)==1){
return 1151;}//B+ decays to anti-D0 D0 K0 rho+ 

if ( PcheckDecay(genpart,-421, 411, 321, -213)==1){
return 1152;}//B+ decays to anti-D0 D+ K+ rho- 

if ( PcheckDecay(genpart,-421, 411, 311, 113)==1){
return 1153;}//B+ decays to anti-D0 D+ K0 rho0 

if ( PcheckDecay(genpart,-411, 411, 311, 213)==1){
return 1154;}//B+ decays to D- D+ K0 rho+ 

if ( PcheckDecay(genpart,-411, 411, 321)==1){
return 1155;}//B+ decays to D- D+ K+ 

if ( PcheckDecay(genpart,-411, 411, 323)==1){
return 1156;}//B+ decays to D- D+ K*+ 

if ( PcheckDecay(genpart,-413, 411, 321)==1){
return 1157;}//B+ decays to D*- D+ K+ 

if ( PcheckDecay(genpart,-413, 411, 323)==1){
return 1158;}//B+ decays to D*- D+ K*+ 

if ( PcheckDecay(genpart,-411, 413, 321)==1){
return 1159;}//B+ decays to D- D*+ K+ 

if ( PcheckDecay(genpart,-411, 413, 323)==1){
return 1160;}//B+ decays to D- D*+ K*+ 

if ( PcheckDecay(genpart,-413, 413, 321)==1){
return 1161;}//B+ decays to D*- D*+ K+ 

if ( PcheckDecay(genpart,-413, 413, 323)==1){
return 1162;}//B+ decays to D*- D*+ K*+ 

if ( PcheckDecay(genpart,-10413, 411, 321)==1){
return 1163;}//B+ decays to D_1- D+ K+ 

if ( PcheckDecay(genpart,-411, 10413, 321)==1){
return 1164;}//B+ decays to D- D_1+ K+ 

if ( PcheckDecay(genpart,-10413, 413, 321)==1){
return 1165;}//B+ decays to D_1- D*+ K+ 

if ( PcheckDecay(genpart,-413, 10413, 321)==1){
return 1166;}//B+ decays to D*- D_1+ K+ 

if ( PcheckDecay(genpart,-415, 411, 321)==1){
return 1167;}//B+ decays to D_2*- D+ K+ 

if ( PcheckDecay(genpart,-411, 415, 321)==1){
return 1168;}//B+ decays to D- D_2*+ K+ 

if ( PcheckDecay(genpart,-415, 413, 321)==1){
return 1169;}//B+ decays to D_2*- D*+ K+ 

if ( PcheckDecay(genpart,-413, 415, 321)==1){
return 1170;}//B+ decays to D*- D_2*+ K+ 

if ( PcheckDecay(genpart,-421, 211)==1){
return 1171;}//B+ decays to anti-D0 pi+ 

if ( PcheckDecay(genpart,-421, 111, 211)==1){
return 1172;}//B+ decays to anti-D0 pi0 pi+ 

if ( PcheckDecay(genpart,213, -421)==1){
return 1173;}//B+ decays to rho+ anti-D0 

if ( PcheckDecay(genpart,-421, -211, 211, 211)==1){
return 1174;}//B+ decays to anti-D0 pi- pi+ pi+ 

if ( PcheckDecay(genpart,-421, 113, 211)==1){
return 1175;}//B+ decays to anti-D0 rho0 pi+ 

if ( PcheckDecay(genpart,20213, -421)==1){
return 1176;}//B+ decays to a_1+ anti-D0 

if ( PcheckDecay(genpart,10213, -421)==1){
return 1177;}//B+ decays to b_1+ anti-D0 

if ( PcheckDecay(genpart,-421, -211, 211, 211, 111)==1){
return 1178;}//B+ decays to anti-D0 pi- pi+ pi+ pi0 

if ( PcheckDecay(genpart,-421, 223, 211)==1){
return 1179;}//B+ decays to anti-D0 omega pi+ 

if ( PcheckDecay(genpart,-423, 211)==1){
return 1180;}//B+ decays to anti-D*0 pi+ 

if ( PcheckDecay(genpart,-423, 111, 211)==1){
return 1181;}//B+ decays to anti-D*0 pi0 pi+ 

if ( PcheckDecay(genpart,-423, 213)==1){
return 1182;}//B+ decays to anti-D*0 rho+ 

if ( PcheckDecay(genpart,-423, -211, 211, 211)==1){
return 1183;}//B+ decays to anti-D*0 pi- pi+ pi+ 

if ( PcheckDecay(genpart,-423, 113, 211)==1){
return 1184;}//B+ decays to anti-D*0 rho0 pi+ 

if ( PcheckDecay(genpart,-423, 20213)==1){
return 1185;}//B+ decays to anti-D*0 a_1+ 

if ( PcheckDecay(genpart,-423, 10213)==1){
return 1186;}//B+ decays to anti-D*0 b_1+ 

if ( PcheckDecay(genpart,-423, 211, 111, 111)==1){
return 1187;}//B+ decays to anti-D*0 pi+ pi0 pi0 

if ( PcheckDecay(genpart,-423, 213, 111)==1){
return 1188;}//B+ decays to anti-D*0 rho+ pi0 

if ( PcheckDecay(genpart,-423, -211, 211, 211, 111)==1){
return 1189;}//B+ decays to anti-D*0 pi- pi+ pi+ pi0 

if ( PcheckDecay(genpart,-423, 223, 211)==1){
return 1190;}//B+ decays to anti-D*0 omega pi+ 

if ( PcheckDecay(genpart,-411, 211, 211)==1){
return 1191;}//B+ decays to D- pi+ pi+ 

if ( PcheckDecay(genpart,-411, 111, 211, 211)==1){
return 1192;}//B+ decays to D- pi0 pi+ pi+ 

if ( PcheckDecay(genpart,-411, 213, 211)==1){
return 1193;}//B+ decays to D- rho+ pi+ 

if ( PcheckDecay(genpart,-413, 211, 211)==1){
return 1194;}//B+ decays to D*- pi+ pi+ 

if ( PcheckDecay(genpart,-413, 111, 211, 211)==1){
return 1195;}//B+ decays to D*- pi0 pi+ pi+ 

if ( PcheckDecay(genpart,-413, 213, 211)==1){
return 1196;}//B+ decays to D*- rho+ pi+ 

if ( PcheckDecay(genpart,-10421, 211)==1){
return 1197;}//B+ decays to anti-D_0*0 pi+ 

if ( PcheckDecay(genpart,-10423, 211)==1){
return 1198;}//B+ decays to anti-D_10 pi+ 

if ( PcheckDecay(genpart,-20423, 211)==1){
return 1199;}//B+ decays to anti-D'_10 pi+ 

if ( PcheckDecay(genpart,-425, 211)==1){
return 1200;}//B+ decays to anti-D_2*0 pi+ 

if ( PcheckDecay(genpart,213, -10421)==1){
return 1201;}//B+ decays to rho+ anti-D_0*0 

if ( PcheckDecay(genpart,-10423, 213)==1){
return 1202;}//B+ decays to anti-D_10 rho+ 

if ( PcheckDecay(genpart,-20423, 213)==1){
return 1203;}//B+ decays to anti-D'_10 rho+ 

if ( PcheckDecay(genpart,-425, 213)==1){
return 1204;}//B+ decays to anti-D_2*0 rho+ 

if ( PcheckDecay(genpart,-421, 321)==1){
return 1205;}//B+ decays to anti-D0 K+ 

if ( PcheckDecay(genpart,323, -421)==1){
return 1206;}//B+ decays to K*+ anti-D0 

if ( PcheckDecay(genpart,-423, 321)==1){
return 1207;}//B+ decays to anti-D*0 K+ 

if ( PcheckDecay(genpart,-423, 323)==1){
return 1208;}//B+ decays to anti-D*0 K*+ 

if ( PcheckDecay(genpart,2, -3, -4, 2)==1){
return 1209;}//B+ decays to u anti-s anti-c u 

if ( PcheckDecay(genpart,2, -1, -4, 2)==1){
return 1210;}//B+ decays to u anti-d anti-c u 

if ( PcheckDecay(genpart,-4122, 2224)==1){
return 1211;}//B+ decays to anti-Lambda_c- Delta++ 

if ( PcheckDecay(genpart,-4212, 2224)==1){
return 1212;}//B+ decays to anti-Sigma_c- Delta++ 

if ( PcheckDecay(genpart,-4112, 2214)==1){
return 1213;}//B+ decays to anti-Sigma_c0 Delta+ 

if ( PcheckDecay(genpart,-4112, 2212)==1){
return 1214;}//B+ decays to anti-Sigma_c0 p+ 

if ( PcheckDecay(genpart,-4103, 2203)==1){
return 1215;}//B+ decays to anti-cd_1 uu_1 

if ( PcheckDecay(genpart,-4103, 2203)==1){
return 1216;}//B+ decays to anti-cd_1 uu_1 

if ( PcheckDecay(genpart,-4232, 4222)==1){
return 1217;}//B+ decays to anti-Xi_c- Sigma_c++ 

if ( PcheckDecay(genpart,-4132, 4122)==1){
return 1218;}//B+ decays to anti-Xi_c0 Lambda_c+ 

if ( PcheckDecay(genpart,-4132, 4212)==1){
return 1219;}//B+ decays to anti-Xi_c0 Sigma_c+ 

if ( PcheckDecay(genpart,-4301, 4201)==1){
return 1220;}//B+ decays to anti-cs_0 cu_0 

if ( PcheckDecay(genpart,-4232, 2224)==1){
return 1221;}//B+ decays to anti-Xi_c- Delta++ 

if ( PcheckDecay(genpart,-4132, 2212)==1){
return 1222;}//B+ decays to anti-Xi_c0 p+ 

if ( PcheckDecay(genpart,-4303, 2203)==1){
return 1223;}//B+ decays to anti-cs_1 uu_1 

if ( PcheckDecay(genpart,-4303, 2203)==1){
return 1224;}//B+ decays to anti-cs_1 uu_1 

if ( PcheckDecay(genpart,-4122, 4222)==1){
return 1225;}//B+ decays to anti-Lambda_c- Sigma_c++ 

if ( PcheckDecay(genpart,-4212, 4222)==1){
return 1226;}//B+ decays to anti-Sigma_c- Sigma_c++ 

if ( PcheckDecay(genpart,-4112, 4212)==1){
return 1227;}//B+ decays to anti-Sigma_c0 Sigma_c+ 

if ( PcheckDecay(genpart,-4101, 4201)==1){
return 1228;}//B+ decays to anti-cd_0 cu_0 
   return -300;
 }

  // ed of B+ loop
 else if (AmId(genpart)==-521){

if ( PcheckDecay(genpart,423, 11, -12)==1){
return 1229;}//B- decays to D*0 e- anti-nu_e 

if ( PcheckDecay(genpart,421, 11, -12)==1){
return 1230;}//B- decays to D0 e- anti-nu_e 

if ( PcheckDecay(genpart,10423, 11, -12)==1){
return 1231;}//B- decays to D_10 e- anti-nu_e 

if ( PcheckDecay(genpart,10421, 11, -12)==1){
return 1232;}//B- decays to D_0*0 e- anti-nu_e 

if ( PcheckDecay(genpart,20423, 11, -12)==1){
return 1233;}//B- decays to D'_10 e- anti-nu_e 

if ( PcheckDecay(genpart,425, 11, -12)==1){
return 1234;}//B- decays to D_2*0 e- anti-nu_e 

if ( PcheckDecay(genpart,100421, 11, -12)==1){
return 1235;}//B- decays to D(2S)0 e- anti-nu_e 

if ( PcheckDecay(genpart,100423, 11, -12)==1){
return 1236;}//B- decays to D*(2S)0 e- anti-nu_e 

if ( PcheckDecay(genpart,423, 111, 11, -12)==1){
return 1237;}//B- decays to D*0 pi0 e- anti-nu_e 

if ( PcheckDecay(genpart,413, -211, 11, -12)==1){
return 1238;}//B- decays to D*+ pi- e- anti-nu_e 

if ( PcheckDecay(genpart,421, 111, 11, -12)==1){
return 1239;}//B- decays to D0 pi0 e- anti-nu_e 

if ( PcheckDecay(genpart,411, -211, 11, -12)==1){
return 1240;}//B- decays to D+ pi- e- anti-nu_e 

if ( PcheckDecay(genpart,423, 13, -14)==1){
return 1241;}//B- decays to D*0 mu- anti-nu_mu 

if ( PcheckDecay(genpart,421, 13, -14)==1){
return 1242;}//B- decays to D0 mu- anti-nu_mu 

if ( PcheckDecay(genpart,10423, 13, -14)==1){
return 1243;}//B- decays to D_10 mu- anti-nu_mu 

if ( PcheckDecay(genpart,10421, 13, -14)==1){
return 1244;}//B- decays to D_0*0 mu- anti-nu_mu 

if ( PcheckDecay(genpart,20423, 13, -14)==1){
return 1245;}//B- decays to D'_10 mu- anti-nu_mu 

if ( PcheckDecay(genpart,425, 13, -14)==1){
return 1246;}//B- decays to D_2*0 mu- anti-nu_mu 

if ( PcheckDecay(genpart,100421, 13, -14)==1){
return 1247;}//B- decays to D(2S)0 mu- anti-nu_mu 

if ( PcheckDecay(genpart,100423, 13, -14)==1){
return 1248;}//B- decays to D*(2S)0 mu- anti-nu_mu 

if ( PcheckDecay(genpart,423, 111, 13, -14)==1){
return 1249;}//B- decays to D*0 pi0 mu- anti-nu_mu 

if ( PcheckDecay(genpart,413, -211, 13, -14)==1){
return 1250;}//B- decays to D*+ pi- mu- anti-nu_mu 

if ( PcheckDecay(genpart,421, 111, 13, -14)==1){
return 1251;}//B- decays to D0 pi0 mu- anti-nu_mu 

if ( PcheckDecay(genpart,411, -211, 13, -14)==1){
return 1252;}//B- decays to D+ pi- mu- anti-nu_mu 

if ( PcheckDecay(genpart,423, 15, -16)==1){
return 1253;}//B- decays to D*0 tau- anti-nu_tau 

if ( PcheckDecay(genpart,421, 15, -16)==1){
return 1254;}//B- decays to D0 tau- anti-nu_tau 

if ( PcheckDecay(genpart,10423, 15, -16)==1){
return 1255;}//B- decays to D_10 tau- anti-nu_tau 

if ( PcheckDecay(genpart,10421, 15, -16)==1){
return 1256;}//B- decays to D_0*0 tau- anti-nu_tau 

if ( PcheckDecay(genpart,20423, 15, -16)==1){
return 1257;}//B- decays to D'_10 tau- anti-nu_tau 

if ( PcheckDecay(genpart,425, 15, -16)==1){
return 1258;}//B- decays to D_2*0 tau- anti-nu_tau 

if ( PcheckDecay(genpart,441, -321)==1){
return 1259;}//B- decays to eta_c K- 

if ( PcheckDecay(genpart,-323, 441)==1){
return 1260;}//B- decays to K*- eta_c 

if ( PcheckDecay(genpart,441, -311, -211)==1){
return 1261;}//B- decays to eta_c anti-K0 pi- 

if ( PcheckDecay(genpart,441, -321, 111)==1){
return 1262;}//B- decays to eta_c K- pi0 

if ( PcheckDecay(genpart,441, -321, 211, -211)==1){
return 1263;}//B- decays to eta_c K- pi+ pi- 

if ( PcheckDecay(genpart,441, -321, 111, 111)==1){
return 1264;}//B- decays to eta_c K- pi0 pi0 

if ( PcheckDecay(genpart,441, -311, -211, 111)==1){
return 1265;}//B- decays to eta_c anti-K0 pi- pi0 

if ( PcheckDecay(genpart,441, 3122, -2212)==1){
return 1266;}//B- decays to eta_c Lambda0 anti-p- 

if ( PcheckDecay(genpart,441, 3112, -2112)==1){
return 1267;}//B- decays to eta_c Sigma- anti-n0 

if ( PcheckDecay(genpart,441, -30353)==1){
return 1268;}//B- decays to eta_c anti-Xsu 

if ( PcheckDecay(genpart,100441, -321)==1){
return 1269;}//B- decays to eta_c(2S) K- 

if ( PcheckDecay(genpart,-323, 100441)==1){
return 1270;}//B- decays to K*- eta_c(2S) 

if ( PcheckDecay(genpart,100441, -311, -211)==1){
return 1271;}//B- decays to eta_c(2S) anti-K0 pi- 

if ( PcheckDecay(genpart,100441, -321, 111)==1){
return 1272;}//B- decays to eta_c(2S) K- pi0 

if ( PcheckDecay(genpart,100441, -321, 211, -211)==1){
return 1273;}//B- decays to eta_c(2S) K- pi+ pi- 

if ( PcheckDecay(genpart,100441, -321, 111, 111)==1){
return 1274;}//B- decays to eta_c(2S) K- pi0 pi0 

if ( PcheckDecay(genpart,100441, -311, -211, 111)==1){
return 1275;}//B- decays to eta_c(2S) anti-K0 pi- pi0 

if ( PcheckDecay(genpart,100441, -30353)==1){
return 1276;}//B- decays to eta_c(2S) anti-Xsu 

if ( PcheckDecay(genpart,443, -321)==1){
return 1277;}//B- decays to psi K- 

if ( PcheckDecay(genpart,443, -323)==1){
return 1278;}//B- decays to psi K*- 

if ( PcheckDecay(genpart,443, -10323)==1){
return 1279;}//B- decays to psi K_1- 

if ( PcheckDecay(genpart,443, -20323)==1){
return 1280;}//B- decays to psi K'_1- 

if ( PcheckDecay(genpart,443, -321, 211, -211)==1){
return 1281;}//B- decays to psi K- pi+ pi- 

if ( PcheckDecay(genpart,443, 333, -321)==1){
return 1282;}//B- decays to psi phi K- 

if ( PcheckDecay(genpart,443, -211)==1){
return 1283;}//B- decays to psi pi- 

if ( PcheckDecay(genpart,443, -213)==1){
return 1284;}//B- decays to psi rho- 

if ( PcheckDecay(genpart,443, -20213)==1){
return 1285;}//B- decays to psi a_1- 

if ( PcheckDecay(genpart,443, -10321)==1){
return 1286;}//B- decays to psi K_0*- 

if ( PcheckDecay(genpart,443, -325)==1){
return 1287;}//B- decays to psi K_2*- 

if ( PcheckDecay(genpart,443, -30323)==1){
return 1288;}//B- decays to psi K''*- 

if ( PcheckDecay(genpart,443, -311, -211)==1){
return 1289;}//B- decays to psi anti-K0 pi- 

if ( PcheckDecay(genpart,443, -321, 111)==1){
return 1290;}//B- decays to psi K- pi0 

if ( PcheckDecay(genpart,443, -321, 111, 111)==1){
return 1291;}//B- decays to psi K- pi0 pi0 

if ( PcheckDecay(genpart,443, -311, -211, 111)==1){
return 1292;}//B- decays to psi anti-K0 pi- pi0 

if ( PcheckDecay(genpart,443, 3122, -2212)==1){
return 1293;}//B- decays to psi Lambda0 anti-p- 

if ( PcheckDecay(genpart,443, 3112, -2112)==1){
return 1294;}//B- decays to psi Sigma- anti-n0 

if ( PcheckDecay(genpart,443, -30353)==1){
return 1295;}//B- decays to psi anti-Xsu 

if ( PcheckDecay(genpart,100443, -321)==1){
return 1296;}//B- decays to psi(2S) K- 

if ( PcheckDecay(genpart,100443, -323)==1){
return 1297;}//B- decays to psi(2S) K*- 

if ( PcheckDecay(genpart,100443, -321, 211, -211)==1){
return 1298;}//B- decays to psi(2S) K- pi+ pi- 

if ( PcheckDecay(genpart,10441, -321)==1){
return 1299;}//B- decays to chi_c0 K- 

if ( PcheckDecay(genpart,-323, 10441)==1){
return 1300;}//B- decays to K*- chi_c0 

if ( PcheckDecay(genpart,10441, -30353)==1){
return 1301;}//B- decays to chi_c0 anti-Xsu 

if ( PcheckDecay(genpart,20443, -321)==1){
return 1302;}//B- decays to chi_c1 K- 

if ( PcheckDecay(genpart,20443, -323)==1){
return 1303;}//B- decays to chi_c1 K*- 

if ( PcheckDecay(genpart,20443, -211)==1){
return 1304;}//B- decays to chi_c1 pi- 

if ( PcheckDecay(genpart,20443, -30353)==1){
return 1305;}//B- decays to chi_c1 anti-Xsu 

if ( PcheckDecay(genpart,445, -321)==1){
return 1306;}//B- decays to chi_c2 K- 

if ( PcheckDecay(genpart,445, -323)==1){
return 1307;}//B- decays to chi_c2 K*- 

if ( PcheckDecay(genpart,445, -30353)==1){
return 1308;}//B- decays to chi_c2 anti-Xsu 

if ( PcheckDecay(genpart,10443, -321)==1){
return 1309;}//B- decays to h_c K- 

if ( PcheckDecay(genpart,10443, -323)==1){
return 1310;}//B- decays to h_c K*- 

if ( PcheckDecay(genpart,10443, -30353)==1){
return 1311;}//B- decays to h_c anti-Xsu 

if ( PcheckDecay(genpart,120443, -321)==1){
return 1312;}//B- decays to X(3872) K- 

if ( PcheckDecay(genpart,90000443, -321)==1){
return 1313;}//B- decays to Y(3940) K- 

if ( PcheckDecay(genpart,421, -431)==1){
return 1314;}//B- decays to D0 D_s- 

if ( PcheckDecay(genpart,423, -431)==1){
return 1315;}//B- decays to D*0 D_s- 

if ( PcheckDecay(genpart,-433, 421)==1){
return 1316;}//B- decays to D_s*- D0 

if ( PcheckDecay(genpart,-433, 423)==1){
return 1317;}//B- decays to D_s*- D*0 

if ( PcheckDecay(genpart,-431, 10421)==1){
return 1318;}//B- decays to D_s- D_0*0 

if ( PcheckDecay(genpart,-433, 10421)==1){
return 1319;}//B- decays to D_s*- D_0*0 

if ( PcheckDecay(genpart,20423, -431)==1){
return 1320;}//B- decays to D'_10 D_s- 

if ( PcheckDecay(genpart,20423, -433)==1){
return 1321;}//B- decays to D'_10 D_s*- 

if ( PcheckDecay(genpart,10423, -431)==1){
return 1322;}//B- decays to D_10 D_s- 

if ( PcheckDecay(genpart,10423, -433)==1){
return 1323;}//B- decays to D_10 D_s*- 

if ( PcheckDecay(genpart,425, -431)==1){
return 1324;}//B- decays to D_2*0 D_s- 

if ( PcheckDecay(genpart,425, -433)==1){
return 1325;}//B- decays to D_2*0 D_s*- 

if ( PcheckDecay(genpart,421, -10431)==1){
return 1326;}//B- decays to D0 D_s0*- 

if ( PcheckDecay(genpart,423, -10431)==1){
return 1327;}//B- decays to D*0 D_s0*- 

if ( PcheckDecay(genpart,-20433, 421)==1){
return 1328;}//B- decays to D'_s1- D0 

if ( PcheckDecay(genpart,-20433, 423)==1){
return 1329;}//B- decays to D'_s1- D*0 

if ( PcheckDecay(genpart,-431, 411, -211)==1){
return 1330;}//B- decays to D_s- D+ pi- 

if ( PcheckDecay(genpart,-431, 421, 111)==1){
return 1331;}//B- decays to D_s- D0 pi0 

if ( PcheckDecay(genpart,-433, 411, -211)==1){
return 1332;}//B- decays to D_s*- D+ pi- 

if ( PcheckDecay(genpart,-433, 421, 111)==1){
return 1333;}//B- decays to D_s*- D0 pi0 

if ( PcheckDecay(genpart,-431, 411, -211, 111)==1){
return 1334;}//B- decays to D_s- D+ pi- pi0 

if ( PcheckDecay(genpart,-431, 421, -211, 211)==1){
return 1335;}//B- decays to D_s- D0 pi- pi+ 

if ( PcheckDecay(genpart,-431, 421, 111, 111)==1){
return 1336;}//B- decays to D_s- D0 pi0 pi0 

if ( PcheckDecay(genpart,-433, 411, -211, 111)==1){
return 1337;}//B- decays to D_s*- D+ pi- pi0 

if ( PcheckDecay(genpart,-433, 421, -211, 211)==1){
return 1338;}//B- decays to D_s*- D0 pi- pi+ 

if ( PcheckDecay(genpart,-433, 421, 111, 111)==1){
return 1339;}//B- decays to D_s*- D0 pi0 pi0 

if ( PcheckDecay(genpart,431, -211, -321)==1){
return 1340;}//B- decays to D_s+ pi- K- 

if ( PcheckDecay(genpart,433, -211, -321)==1){
return 1341;}//B- decays to D_s*+ pi- K- 

if ( PcheckDecay(genpart,431, -211, -321, 111)==1){
return 1342;}//B- decays to D_s+ pi- K- pi0 

if ( PcheckDecay(genpart,433, -211, -321, 111)==1){
return 1343;}//B- decays to D_s*+ pi- K- pi0 

if ( PcheckDecay(genpart,431, -211, -311, -211)==1){
return 1344;}//B- decays to D_s+ pi- anti-K0 pi- 

if ( PcheckDecay(genpart,433, -211, -311, -211)==1){
return 1345;}//B- decays to D_s*+ pi- anti-K0 pi- 

if ( PcheckDecay(genpart,431, -211, -321, 111, 111)==1){
return 1346;}//B- decays to D_s+ pi- K- pi0 pi0 

if ( PcheckDecay(genpart,433, -211, -321, 111, 111)==1){
return 1347;}//B- decays to D_s*+ pi- K- pi0 pi0 

if ( PcheckDecay(genpart,431, -211, -321, 211, -211)==1){
return 1348;}//B- decays to D_s+ pi- K- pi+ pi- 

if ( PcheckDecay(genpart,433, -211, -321, 211, -211)==1){
return 1349;}//B- decays to D_s*+ pi- K- pi+ pi- 

if ( PcheckDecay(genpart,431, -211, -311, -211, 111)==1){
return 1350;}//B- decays to D_s+ pi- anti-K0 pi- pi0 

if ( PcheckDecay(genpart,433, -211, -311, -211, 111)==1){
return 1351;}//B- decays to D_s*+ pi- anti-K0 pi- pi0 

if ( PcheckDecay(genpart,-411, 421)==1){
return 1352;}//B- decays to D- D0 

if ( PcheckDecay(genpart,423, -411)==1){
return 1353;}//B- decays to D*0 D- 

if ( PcheckDecay(genpart,-413, 421)==1){
return 1354;}//B- decays to D*- D0 

if ( PcheckDecay(genpart,-413, 423)==1){
return 1355;}//B- decays to D*- D*0 

if ( PcheckDecay(genpart,-20413, 421)==1){
return 1356;}//B- decays to D'_1- D0 

if ( PcheckDecay(genpart,20423, -411)==1){
return 1357;}//B- decays to D'_10 D- 

if ( PcheckDecay(genpart,-20413, 423)==1){
return 1358;}//B- decays to D'_1- D*0 

if ( PcheckDecay(genpart,20423, -413)==1){
return 1359;}//B- decays to D'_10 D*- 

if ( PcheckDecay(genpart,-10413, 421)==1){
return 1360;}//B- decays to D_1- D0 

if ( PcheckDecay(genpart,10423, -411)==1){
return 1361;}//B- decays to D_10 D- 

if ( PcheckDecay(genpart,-10413, 423)==1){
return 1362;}//B- decays to D_1- D*0 

if ( PcheckDecay(genpart,10423, -413)==1){
return 1363;}//B- decays to D_10 D*- 

if ( PcheckDecay(genpart,-415, 421)==1){
return 1364;}//B- decays to D_2*- D0 

if ( PcheckDecay(genpart,425, -411)==1){
return 1365;}//B- decays to D_2*0 D- 

if ( PcheckDecay(genpart,-415, 423)==1){
return 1366;}//B- decays to D_2*- D*0 

if ( PcheckDecay(genpart,425, -413)==1){
return 1367;}//B- decays to D_2*0 D*- 

if ( PcheckDecay(genpart,-10433, 421)==1){
return 1368;}//B- decays to D_s1- D0 

if ( PcheckDecay(genpart,-10433, 423)==1){
return 1369;}//B- decays to D_s1- D*0 

if ( PcheckDecay(genpart,-435, 421)==1){
return 1370;}//B- decays to D_s2*- D0 

if ( PcheckDecay(genpart,-435, 423)==1){
return 1371;}//B- decays to D_s2*- D*0 

if ( PcheckDecay(genpart,-9000433, 421)==1){
return 1372;}//B- decays to D_sj(2700)- D0 

if ( PcheckDecay(genpart,-9000433, 423)==1){
return 1373;}//B- decays to D_sj(2700)- D*0 

if ( PcheckDecay(genpart,30443, -321)==1){
return 1374;}//B- decays to psi(3770) K- 

if ( PcheckDecay(genpart,30443, -323)==1){
return 1375;}//B- decays to psi(3770) K*- 

if ( PcheckDecay(genpart,30443, -311, -211)==1){
return 1376;}//B- decays to psi(3770) anti-K0 pi- 

if ( PcheckDecay(genpart,30443, -321, 111)==1){
return 1377;}//B- decays to psi(3770) K- pi0 

if ( PcheckDecay(genpart,30443, -321, 111, 111)==1){
return 1378;}//B- decays to psi(3770) K- pi0 pi0 

if ( PcheckDecay(genpart,30443, -311, -211, 111)==1){
return 1379;}//B- decays to psi(3770) anti-K0 pi- pi0 

if ( PcheckDecay(genpart,9000443, -321)==1){
return 1380;}//B- decays to psi(4040) K- 

if ( PcheckDecay(genpart,9000443, -323)==1){
return 1381;}//B- decays to psi(4040) K*- 

if ( PcheckDecay(genpart,9000443, -311, -211)==1){
return 1382;}//B- decays to psi(4040) anti-K0 pi- 

if ( PcheckDecay(genpart,9000443, -321, 111)==1){
return 1383;}//B- decays to psi(4040) K- pi0 

if ( PcheckDecay(genpart,9010443, -321)==1){
return 1384;}//B- decays to psi(4160) K- 

if ( PcheckDecay(genpart,9010443, -323)==1){
return 1385;}//B- decays to psi(4160) K*- 

if ( PcheckDecay(genpart,9010443, -311, -211)==1){
return 1386;}//B- decays to psi(4160) anti-K0 pi- 

if ( PcheckDecay(genpart,9010443, -321, 111)==1){
return 1387;}//B- decays to psi(4160) K- pi0 

//if ( PcheckDecay(genpart,return 1388;}//B- decays to 

if ( PcheckDecay(genpart,-411, 421, 111)==1){
return 1389;}//B- decays to D- D0 pi0 

if ( PcheckDecay(genpart,-411, 411, -211)==1){
return 1390;}//B- decays to D- D+ pi- 

if ( PcheckDecay(genpart,-421, 421, -211)==1){
return 1391;}//B- decays to anti-D0 D0 pi- 

if ( PcheckDecay(genpart,-413, 421, 111)==1){
return 1392;}//B- decays to D*- D0 pi0 

if ( PcheckDecay(genpart,-413, 411, -211)==1){
return 1393;}//B- decays to D*- D+ pi- 

if ( PcheckDecay(genpart,-423, 421, -211)==1){
return 1394;}//B- decays to anti-D*0 D0 pi- 

if ( PcheckDecay(genpart,-411, 423, 111)==1){
return 1395;}//B- decays to D- D*0 pi0 

if ( PcheckDecay(genpart,-411, 413, -211)==1){
return 1396;}//B- decays to D- D*+ pi- 

if ( PcheckDecay(genpart,-421, 423, -211)==1){
return 1397;}//B- decays to anti-D0 D*0 pi- 

if ( PcheckDecay(genpart,-413, 423, 111)==1){
return 1398;}//B- decays to D*- D*0 pi0 

if ( PcheckDecay(genpart,-413, 413, -211)==1){
return 1399;}//B- decays to D*- D*+ pi- 

if ( PcheckDecay(genpart,-423, 423, -211)==1){
return 1400;}//B- decays to anti-D*0 D*0 pi- 

if ( PcheckDecay(genpart,-411, 421, 113)==1){
return 1401;}//B- decays to D- D0 rho0 

if ( PcheckDecay(genpart,-411, 411, -213)==1){
return 1402;}//B- decays to D- D+ rho- 

if ( PcheckDecay(genpart,-421, 421, -213)==1){
return 1403;}//B- decays to anti-D0 D0 rho- 

if ( PcheckDecay(genpart,-413, 421, 113)==1){
return 1404;}//B- decays to D*- D0 rho0 

if ( PcheckDecay(genpart,-413, 411, -213)==1){
return 1405;}//B- decays to D*- D+ rho- 

if ( PcheckDecay(genpart,-423, 421, -213)==1){
return 1406;}//B- decays to anti-D*0 D0 rho- 

if ( PcheckDecay(genpart,-411, 423, 113)==1){
return 1407;}//B- decays to D- D*0 rho0 

if ( PcheckDecay(genpart,-411, 413, -213)==1){
return 1408;}//B- decays to D- D*+ rho- 

if ( PcheckDecay(genpart,-421, 423, -213)==1){
return 1409;}//B- decays to anti-D0 D*0 rho- 

if ( PcheckDecay(genpart,-413, 423, 113)==1){
return 1410;}//B- decays to D*- D*0 rho0 

if ( PcheckDecay(genpart,-413, 413, -213)==1){
return 1411;}//B- decays to D*- D*+ rho- 

if ( PcheckDecay(genpart,-423, 423, -213)==1){
return 1412;}//B- decays to anti-D*0 D*0 rho- 

if ( PcheckDecay(genpart,-411, 421, 111, 111)==1){
return 1413;}//B- decays to D- D0 pi0 pi0 

if ( PcheckDecay(genpart,-411, 421, 211, -211)==1){
return 1414;}//B- decays to D- D0 pi+ pi- 

if ( PcheckDecay(genpart,-411, 411, -211, 111)==1){
return 1415;}//B- decays to D- D+ pi- pi0 

if ( PcheckDecay(genpart,-421, 421, -211, 111)==1){
return 1416;}//B- decays to anti-D0 D0 pi- pi0 

if ( PcheckDecay(genpart,-421, 411, -211, -211)==1){
return 1417;}//B- decays to anti-D0 D+ pi- pi- 

if ( PcheckDecay(genpart,-413, 421, 111, 111)==1){
return 1418;}//B- decays to D*- D0 pi0 pi0 

if ( PcheckDecay(genpart,-413, 421, 211, -211)==1){
return 1419;}//B- decays to D*- D0 pi+ pi- 

if ( PcheckDecay(genpart,-413, 411, -211, 111)==1){
return 1420;}//B- decays to D*- D+ pi- pi0 

if ( PcheckDecay(genpart,-423, 421, -211, 111)==1){
return 1421;}//B- decays to anti-D*0 D0 pi- pi0 

if ( PcheckDecay(genpart,-423, 411, -211, -211)==1){
return 1422;}//B- decays to anti-D*0 D+ pi- pi- 

if ( PcheckDecay(genpart,-411, 423, 111, 111)==1){
return 1423;}//B- decays to D- D*0 pi0 pi0 

if ( PcheckDecay(genpart,-411, 423, 211, -211)==1){
return 1424;}//B- decays to D- D*0 pi+ pi- 

if ( PcheckDecay(genpart,-411, 413, -211, 111)==1){
return 1425;}//B- decays to D- D*+ pi- pi0 

if ( PcheckDecay(genpart,-421, 423, -211, 111)==1){
return 1426;}//B- decays to anti-D0 D*0 pi- pi0 

if ( PcheckDecay(genpart,-421, 413, -211, -211)==1){
return 1427;}//B- decays to anti-D0 D*+ pi- pi- 

if ( PcheckDecay(genpart,-413, 423, 111, 111)==1){
return 1428;}//B- decays to D*- D*0 pi0 pi0 

if ( PcheckDecay(genpart,-413, 423, 211, -211)==1){
return 1429;}//B- decays to D*- D*0 pi+ pi- 

if ( PcheckDecay(genpart,-413, 413, -211, 111)==1){
return 1430;}//B- decays to D*- D*+ pi- pi0 

if ( PcheckDecay(genpart,-423, 423, -211, 111)==1){
return 1431;}//B- decays to anti-D*0 D*0 pi- pi0 

if ( PcheckDecay(genpart,-423, 413, -211, -211)==1){
return 1432;}//B- decays to anti-D*0 D*+ pi- pi- 

if ( PcheckDecay(genpart,421, -421, -321)==1){
return 1433;}//B- decays to D0 anti-D0 K- 

if ( PcheckDecay(genpart,421, -411, -311)==1){
return 1434;}//B- decays to D0 D- anti-K0 

if ( PcheckDecay(genpart,421, -421, -323)==1){
return 1435;}//B- decays to D0 anti-D0 K*- 

if ( PcheckDecay(genpart,421, -411, -313)==1){
return 1436;}//B- decays to D0 D- anti-K*0 

if ( PcheckDecay(genpart,423, -421, -321)==1){
return 1437;}//B- decays to D*0 anti-D0 K- 

if ( PcheckDecay(genpart,423, -411, -311)==1){
return 1438;}//B- decays to D*0 D- anti-K0 

if ( PcheckDecay(genpart,421, -423, -321)==1){
return 1439;}//B- decays to D0 anti-D*0 K- 

if ( PcheckDecay(genpart,421, -413, -311)==1){
return 1440;}//B- decays to D0 D*- anti-K0 

if ( PcheckDecay(genpart,423, -421, -323)==1){
return 1441;}//B- decays to D*0 anti-D0 K*- 

if ( PcheckDecay(genpart,423, -411, -313)==1){
return 1442;}//B- decays to D*0 D- anti-K*0 

if ( PcheckDecay(genpart,421, -423, -323)==1){
return 1443;}//B- decays to D0 anti-D*0 K*- 

if ( PcheckDecay(genpart,421, -413, -313)==1){
return 1444;}//B- decays to D0 D*- anti-K*0 

if ( PcheckDecay(genpart,423, -423, -321)==1){
return 1445;}//B- decays to D*0 anti-D*0 K- 

if ( PcheckDecay(genpart,423, -413, -311)==1){
return 1446;}//B- decays to D*0 D*- anti-K0 

if ( PcheckDecay(genpart,423, -423, -323)==1){
return 1447;}//B- decays to D*0 anti-D*0 K*- 

if ( PcheckDecay(genpart,423, -413, -313)==1){
return 1448;}//B- decays to D*0 D*- anti-K*0 

if ( PcheckDecay(genpart,10423, -421, -321)==1){
return 1449;}//B- decays to D_10 anti-D0 K- 

if ( PcheckDecay(genpart,10423, -411, -311)==1){
return 1450;}//B- decays to D_10 D- anti-K0 

if ( PcheckDecay(genpart,421, -10423, -321)==1){
return 1451;}//B- decays to D0 anti-D_10 K- 

if ( PcheckDecay(genpart,421, -10413, -311)==1){
return 1452;}//B- decays to D0 D_1- anti-K0 

if ( PcheckDecay(genpart,10423, -423, -321)==1){
return 1453;}//B- decays to D_10 anti-D*0 K- 

if ( PcheckDecay(genpart,10423, -413, -311)==1){
return 1454;}//B- decays to D_10 D*- anti-K0 

if ( PcheckDecay(genpart,423, -10423, -321)==1){
return 1455;}//B- decays to D*0 anti-D_10 K- 

if ( PcheckDecay(genpart,423, -10413, -311)==1){
return 1456;}//B- decays to D*0 D_1- anti-K0 

if ( PcheckDecay(genpart,425, -421, -321)==1){
return 1457;}//B- decays to D_2*0 anti-D0 K- 

if ( PcheckDecay(genpart,425, -411, -311)==1){
return 1458;}//B- decays to D_2*0 D- anti-K0 

if ( PcheckDecay(genpart,421, -425, -321)==1){
return 1459;}//B- decays to D0 anti-D_2*0 K- 

if ( PcheckDecay(genpart,421, -415, -311)==1){
return 1460;}//B- decays to D0 D_2*- anti-K0 

if ( PcheckDecay(genpart,425, -423, -321)==1){
return 1461;}//B- decays to D_2*0 anti-D*0 K- 

if ( PcheckDecay(genpart,425, -413, -311)==1){
return 1462;}//B- decays to D_2*0 D*- anti-K0 

if ( PcheckDecay(genpart,423, -425, -321)==1){
return 1463;}//B- decays to D*0 anti-D_2*0 K- 

if ( PcheckDecay(genpart,423, -415, -311)==1){
return 1464;}//B- decays to D*0 D_2*- anti-K0 

if ( PcheckDecay(genpart,421, -421, -321, 111)==1){
return 1465;}//B- decays to D0 anti-D0 K- pi0 

if ( PcheckDecay(genpart,411, -421, -321, -211)==1){
return 1466;}//B- decays to D+ anti-D0 K- pi- 

if ( PcheckDecay(genpart,421, -421, -311, -211)==1){
return 1467;}//B- decays to D0 anti-D0 anti-K0 pi- 

if ( PcheckDecay(genpart,421, -411, -321, 211)==1){
return 1468;}//B- decays to D0 D- K- pi+ 

if ( PcheckDecay(genpart,421, -411, -311, 111)==1){
return 1469;}//B- decays to D0 D- anti-K0 pi0 

if ( PcheckDecay(genpart,411, -411, -311, -211)==1){
return 1470;}//B- decays to D+ D- anti-K0 pi- 

if ( PcheckDecay(genpart,423, -421, -321, 111)==1){
return 1471;}//B- decays to D*0 anti-D0 K- pi0 

if ( PcheckDecay(genpart,413, -421, -321, -211)==1){
return 1472;}//B- decays to D*+ anti-D0 K- pi- 

if ( PcheckDecay(genpart,423, -421, -311, -211)==1){
return 1473;}//B- decays to D*0 anti-D0 anti-K0 pi- 

if ( PcheckDecay(genpart,423, -411, -321, 211)==1){
return 1474;}//B- decays to D*0 D- K- pi+ 

if ( PcheckDecay(genpart,423, -411, -311, 111)==1){
return 1475;}//B- decays to D*0 D- anti-K0 pi0 

if ( PcheckDecay(genpart,413, -411, -311, -211)==1){
return 1476;}//B- decays to D*+ D- anti-K0 pi- 

if ( PcheckDecay(genpart,421, -423, -321, 111)==1){
return 1477;}//B- decays to D0 anti-D*0 K- pi0 

if ( PcheckDecay(genpart,411, -423, -321, -211)==1){
return 1478;}//B- decays to D+ anti-D*0 K- pi- 

if ( PcheckDecay(genpart,421, -423, -311, -211)==1){
return 1479;}//B- decays to D0 anti-D*0 anti-K0 pi- 

if ( PcheckDecay(genpart,421, -413, -321, 211)==1){
return 1480;}//B- decays to D0 D*- K- pi+ 

if ( PcheckDecay(genpart,421, -413, -311, 111)==1){
return 1481;}//B- decays to D0 D*- anti-K0 pi0 

if ( PcheckDecay(genpart,411, -413, -311, -211)==1){
return 1482;}//B- decays to D+ D*- anti-K0 pi- 

if ( PcheckDecay(genpart,423, -423, -321, 111)==1){
return 1483;}//B- decays to D*0 anti-D*0 K- pi0 

if ( PcheckDecay(genpart,413, -423, -321, -211)==1){
return 1484;}//B- decays to D*+ anti-D*0 K- pi- 

if ( PcheckDecay(genpart,423, -423, -311, -211)==1){
return 1485;}//B- decays to D*0 anti-D*0 anti-K0 pi- 

if ( PcheckDecay(genpart,423, -413, -321, 211)==1){
return 1486;}//B- decays to D*0 D*- K- pi+ 

if ( PcheckDecay(genpart,423, -413, -311, 111)==1){
return 1487;}//B- decays to D*0 D*- anti-K0 pi0 

if ( PcheckDecay(genpart,413, -413, -311, -211)==1){
return 1488;}//B- decays to D*+ D*- anti-K0 pi- 

if ( PcheckDecay(genpart,421, -421, -323, 111)==1){
return 1489;}//B- decays to D0 anti-D0 K*- pi0 

if ( PcheckDecay(genpart,411, -421, -323, -211)==1){
return 1490;}//B- decays to D+ anti-D0 K*- pi- 

if ( PcheckDecay(genpart,421, -421, -313, -211)==1){
return 1491;}//B- decays to D0 anti-D0 anti-K*0 pi- 

if ( PcheckDecay(genpart,421, -411, -323, 211)==1){
return 1492;}//B- decays to D0 D- K*- pi+ 

if ( PcheckDecay(genpart,421, -411, -313, 111)==1){
return 1493;}//B- decays to D0 D- anti-K*0 pi0 

if ( PcheckDecay(genpart,411, -411, -313, -211)==1){
return 1494;}//B- decays to D+ D- anti-K*0 pi- 

if ( PcheckDecay(genpart,421, -421, -321, 113)==1){
return 1495;}//B- decays to D0 anti-D0 K- rho0 

if ( PcheckDecay(genpart,411, -421, -321, -213)==1){
return 1496;}//B- decays to D+ anti-D0 K- rho- 

if ( PcheckDecay(genpart,421, -421, -311, -213)==1){
return 1497;}//B- decays to D0 anti-D0 anti-K0 rho- 

if ( PcheckDecay(genpart,421, -411, -321, 213)==1){
return 1498;}//B- decays to D0 D- K- rho+ 

if ( PcheckDecay(genpart,421, -411, -311, 113)==1){
return 1499;}//B- decays to D0 D- anti-K0 rho0 

if ( PcheckDecay(genpart,411, -411, -311, -213)==1){
return 1500;}//B- decays to D+ D- anti-K0 rho- 

if ( PcheckDecay(genpart,411, -411, -321)==1){
return 1501;}//B- decays to D+ D- K- 

if ( PcheckDecay(genpart,411, -411, -323)==1){
return 1502;}//B- decays to D+ D- K*- 

if ( PcheckDecay(genpart,413, -411, -321)==1){
return 1503;}//B- decays to D*+ D- K- 

if ( PcheckDecay(genpart,413, -411, -323)==1){
return 1504;}//B- decays to D*+ D- K*- 

if ( PcheckDecay(genpart,411, -413, -321)==1){
return 1505;}//B- decays to D+ D*- K- 

if ( PcheckDecay(genpart,411, -413, -323)==1){
return 1506;}//B- decays to D+ D*- K*- 

if ( PcheckDecay(genpart,413, -413, -321)==1){
return 1507;}//B- decays to D*+ D*- K- 

if ( PcheckDecay(genpart,413, -413, -323)==1){
return 1508;}//B- decays to D*+ D*- K*- 

if ( PcheckDecay(genpart,10413, -411, -321)==1){
return 1509;}//B- decays to D_1+ D- K- 

if ( PcheckDecay(genpart,411, -10413, -321)==1){
return 1510;}//B- decays to D+ D_1- K- 

if ( PcheckDecay(genpart,10413, -413, -321)==1){
return 1511;}//B- decays to D_1+ D*- K- 

if ( PcheckDecay(genpart,413, -10413, -321)==1){
return 1512;}//B- decays to D*+ D_1- K- 

if ( PcheckDecay(genpart,415, -411, -321)==1){
return 1513;}//B- decays to D_2*+ D- K- 

if ( PcheckDecay(genpart,411, -415, -321)==1){
return 1514;}//B- decays to D+ D_2*- K- 

if ( PcheckDecay(genpart,415, -413, -321)==1){
return 1515;}//B- decays to D_2*+ D*- K- 

if ( PcheckDecay(genpart,413, -415, -321)==1){
return 1516;}//B- decays to D*+ D_2*- K- 

if ( PcheckDecay(genpart,421, -211)==1){
return 1517;}//B- decays to D0 pi- 

if ( PcheckDecay(genpart,421, 111, -211)==1){
return 1518;}//B- decays to D0 pi0 pi- 

if ( PcheckDecay(genpart,-213, 421)==1){
return 1519;}//B- decays to rho- D0 

if ( PcheckDecay(genpart,421, 211, -211, -211)==1){
return 1520;}//B- decays to D0 pi+ pi- pi- 

if ( PcheckDecay(genpart,421, 113, -211)==1){
return 1521;}//B- decays to D0 rho0 pi- 

if ( PcheckDecay(genpart,-20213, 421)==1){
return 1522;}//B- decays to a_1- D0 

if ( PcheckDecay(genpart,-10213, 421)==1){
return 1523;}//B- decays to b_1- D0 

if ( PcheckDecay(genpart,421, 211, -211, -211, 111)==1){
return 1524;}//B- decays to D0 pi+ pi- pi- pi0 

if ( PcheckDecay(genpart,421, 223, -211)==1){
return 1525;}//B- decays to D0 omega pi- 

if ( PcheckDecay(genpart,423, -211)==1){
return 1526;}//B- decays to D*0 pi- 

if ( PcheckDecay(genpart,423, 111, -211)==1){
return 1527;}//B- decays to D*0 pi0 pi- 

if ( PcheckDecay(genpart,423, -213)==1){
return 1528;}//B- decays to D*0 rho- 

if ( PcheckDecay(genpart,423, 211, -211, -211)==1){
return 1529;}//B- decays to D*0 pi+ pi- pi- 

if ( PcheckDecay(genpart,423, 113, -211)==1){
return 1530;}//B- decays to D*0 rho0 pi- 

if ( PcheckDecay(genpart,423, -20213)==1){
return 1531;}//B- decays to D*0 a_1- 

if ( PcheckDecay(genpart,423, -10213)==1){
return 1532;}//B- decays to D*0 b_1- 

if ( PcheckDecay(genpart,423, -211, 111, 111)==1){
return 1533;}//B- decays to D*0 pi- pi0 pi0 

if ( PcheckDecay(genpart,423, -213, 111)==1){
return 1534;}//B- decays to D*0 rho- pi0 

if ( PcheckDecay(genpart,423, 211, -211, -211, 111)==1){
return 1535;}//B- decays to D*0 pi+ pi- pi- pi0 

if ( PcheckDecay(genpart,423, 223, -211)==1){
return 1536;}//B- decays to D*0 omega pi- 

if ( PcheckDecay(genpart,411, -211, -211)==1){
return 1537;}//B- decays to D+ pi- pi- 

if ( PcheckDecay(genpart,411, 111, -211, -211)==1){
return 1538;}//B- decays to D+ pi0 pi- pi- 

if ( PcheckDecay(genpart,411, -213, -211)==1){
return 1539;}//B- decays to D+ rho- pi- 

if ( PcheckDecay(genpart,413, -211, -211)==1){
return 1540;}//B- decays to D*+ pi- pi- 

if ( PcheckDecay(genpart,413, 111, -211, -211)==1){
return 1541;}//B- decays to D*+ pi0 pi- pi- 

if ( PcheckDecay(genpart,413, -213, -211)==1){
return 1542;}//B- decays to D*+ rho- pi- 

if ( PcheckDecay(genpart,10421, -211)==1){
return 1543;}//B- decays to D_0*0 pi- 

if ( PcheckDecay(genpart,10423, -211)==1){
return 1544;}//B- decays to D_10 pi- 

if ( PcheckDecay(genpart,20423, -211)==1){
return 1545;}//B- decays to D'_10 pi- 

if ( PcheckDecay(genpart,425, -211)==1){
return 1546;}//B- decays to D_2*0 pi- 

if ( PcheckDecay(genpart,-213, 10421)==1){
return 1547;}//B- decays to rho- D_0*0 

if ( PcheckDecay(genpart,10423, -213)==1){
return 1548;}//B- decays to D_10 rho- 

if ( PcheckDecay(genpart,20423, -213)==1){
return 1549;}//B- decays to D'_10 rho- 

if ( PcheckDecay(genpart,425, -213)==1){
return 1550;}//B- decays to D_2*0 rho- 

if ( PcheckDecay(genpart,421, -321)==1){
return 1551;}//B- decays to D0 K- 

if ( PcheckDecay(genpart,-323, 421)==1){
return 1552;}//B- decays to K*- D0 

if ( PcheckDecay(genpart,423, -321)==1){
return 1553;}//B- decays to D*0 K- 

if ( PcheckDecay(genpart,423, -323)==1){
return 1554;}//B- decays to D*0 K*- 

if ( PcheckDecay(genpart,-2, 3, 4, -2)==1){
return 1555;}//B- decays to anti-u s c anti-u 

if ( PcheckDecay(genpart,-2, 1, 4, -2)==1){
return 1556;}//B- decays to anti-u d c anti-u 

if ( PcheckDecay(genpart,4122, -2224)==1){
return 1557;}//B- decays to Lambda_c+ anti-Delta-- 

if ( PcheckDecay(genpart,4212, -2224)==1){
return 1558;}//B- decays to Sigma_c+ anti-Delta-- 

if ( PcheckDecay(genpart,4112, -2214)==1){
return 1559;}//B- decays to Sigma_c0 anti-Delta- 

if ( PcheckDecay(genpart,4112, -2212)==1){
return 1560;}//B- decays to Sigma_c0 anti-p- 

if ( PcheckDecay(genpart,4103, -2203)==1){
return 1561;}//B- decays to cd_1 anti-uu_1 

if ( PcheckDecay(genpart,4103, -2203)==1){
return 1562;}//B- decays to cd_1 anti-uu_1 

if ( PcheckDecay(genpart,4232, -4222)==1){
return 1563;}//B- decays to Xi_c+ anti-Sigma_c-- 

if ( PcheckDecay(genpart,4132, -4122)==1){
return 1564;}//B- decays to Xi_c0 anti-Lambda_c- 

if ( PcheckDecay(genpart,4132, -4212)==1){
return 1565;}//B- decays to Xi_c0 anti-Sigma_c- 

if ( PcheckDecay(genpart,4301, -4201)==1){
return 1566;}//B- decays to cs_0 anti-cu_0 

if ( PcheckDecay(genpart,4232, -2224)==1){
return 1567;}//B- decays to Xi_c+ anti-Delta-- 

if ( PcheckDecay(genpart,4132, -2212)==1){
return 1568;}//B- decays to Xi_c0 anti-p- 

if ( PcheckDecay(genpart,4303, -2203)==1){
return 1569;}//B- decays to cs_1 anti-uu_1 

if ( PcheckDecay(genpart,4303, -2203)==1){
return 1570;}//B- decays to cs_1 anti-uu_1 

if ( PcheckDecay(genpart,4122, -4222)==1){
return 1571;}//B- decays to Lambda_c+ anti-Sigma_c-- 

if ( PcheckDecay(genpart,4212, -4222)==1){
return 1572;}//B- decays to Sigma_c+ anti-Sigma_c-- 

if ( PcheckDecay(genpart,4112, -4212)==1){
return 1573;}//B- decays to Sigma_c0 anti-Sigma_c- 

if ( PcheckDecay(genpart,4101, -4201)==1){
return 1574;}//B- decays to cd_0 anti-cu_0 
   return -400;
 } // end of B-
  else return -9876789;
}


     
  
  








#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
