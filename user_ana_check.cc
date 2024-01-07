////////////////////////////////D+->ep///////////////////////////////////////
//  date jan 26.04.2020
// Souvik Maity
//////////////////////////////////////////////////////////////////////////////

// header files/////
#include "belle.h"
#include <math.h>
#include <time.h>
//#include <stream>
#include <iomanip>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"
#include "panther/panther.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Geometry/Point3D.h"
#include "eid/eid.h"        // to use eid class.
#include "mdst/Muid_mdst.h" // to use Muid_mdst class.
#include "mdst/mdst.h"      // to use Benergy func.
#include "kid/atc_pid.h"    // to use atc_pid class.
#include "helix/Helix.h"    // to use helix class.
#include "kfitter/kvertexfitter.h" // to use vertexfitter.
#include "kfitter/kmassfitter.h" // to use massfitter.
#include "particle/Particle.h"
#include "particle/Relation.h"
#include "particle/Momentum.h"
#include "particle/utility.h"
#include "ip/IpProfile.h"
#include "UserInfo.h"
#include "mdst/findKs.h"
#include "particle/combination.h"
#include "benergy/BeamEnergy.h"
#include "./geninfo.h"
//#include "./get_gen.h"
#include "./myutl.h"
#include "./get_gen_tag.h"

#include BELLETDF_H
#include BRECON_H
#include HEPEVT_H
#include MDST_H
#include EVTCLS_H



#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif
   class user_ana_check : public Module {
   public:
     user_ana_check (void);  // constructor
     ~user_ana_check (void){}; // destructor
     
     void init (int*);  //to be compatible with 20000430
     void term (){};
     void disp_stat (const char*) {}; // to be compatible with 20000430
     void hist_def (void);
     void event (BelleEvent*, int*);
     void begin_run ( BelleEvent*, int* );
     void compare(double x, int k, double &temp, int &tk);
     void end_run ( BelleEvent*, int * ){};
     void other (int*, BelleEvent*, int* ){};
     void vertexFit(Particle &p);//new addition
     void vertexDstrFit(Particle &p);
     bool refitSlowPion(Particle &p,const HepPoint3D &ip);
     int ievent,_everyNevt,nB;
   int iexpmc;  //get_gen   
   private:
     Hep3Vector CMBoost; HepLorentzVector cm, p_cm, p_beam;
     double E_LER,E_HER,E_beam,theta, oE_LER,oE_HER,er_E_beam,Ecm,oE_beam;
     BelleTuple* mnt_ds;
     BelleTuple* mnt_dst;
     motherid Get;
     myVfit Vf;
     genT Gen;
     //genT Get;
  }; //CLASS user_ana_check :public Module

   //____________________________________________________
   //          BASF  interface function 
   //____________________________________________________
     extern "C" Module_descr *mdcl_user_ana_check()
     {
       //Create the module
       user_ana_check *module = new user_ana_check;
       //Create the description of the module
       Module_descr *dscr = new Module_descr ( "user_ana_check", module);
       BeamEnergy::define_global(dscr);
       return dscr;
     }

     user_ana_check::user_ana_check(void) {     //construction
     ievent=0; nB=0; _everyNevt=1; }

     void user_ana_check::init ( int*){       //Initilalize function
     }

     void user_ana_check:: hist_def (void)     // histogram definition
     {

     extern BelleTupleManager* BASF_Histogram;
     BelleTupleManager* tm  = BASF_Histogram;
     mnt_ds  = tm->ntuple("Ds", "  mass mom d0mass q deltam dstrf d0f elf prf elec_mom pr_mom elec_mass pr_mass mult_ds mult_d0 kid dr pissvd elsvd prsvd masspis pismom chargepis dz elec_dz pr_dz elec_dr pr_dr pr_id elec_id elec_th pr_th elec_th1 pr_th1 pis_dedx el_dedx pr_dedx modedz modeadz modedp modedm d0chisq dschisq e1_l mu_l prk1_l prp1_l pis_cos elec_cos pr_cos mul0 mul1 mul2 mul3 mul4 mul5 mul6 mul7 mul8 mul9 bc0 bc1 bc2 bc3 bc4 bc5 bc6 bc7 bc8 bc9 etag lma1 lma2 prma1 prma2 pima1 pis_p pis_k k1svdx k1svdy k1svdz k2svdx k2svdy k2svdz k2_mom k1_mom ");
     mnt_dst= tm->ntuple("Dst","mass deltam etag");
  }
  //________________________________________________________________
  //                EVENT FUNCION LOOP START FROM HERE
  //________________________________________________________________
  void user_ana_check::event ( BelleEvent* evptr, int* status) {
    // ...... HAMLET OBJECT declaration...
    const double RadDeg=180.0/M_PI;     // Radians to Degrees
    const double DegRad=M_PI/180.0;     // Degrees to Radians
    E_LER     = BeamEnergy::E_LER();    //Beam Energy
    E_HER     = BeamEnergy::E_HER();
    E_beam    = BeamEnergy::E_beam_corr();
    oE_beam   = BeamEnergy::E_beam_orig();
    er_E_beam = BeamEnergy::E_beam_err();
    theta     = BeamEnergy::Cross_angle();
    CMBoost   = - BeamEnergy::CMBoost();
    oE_LER    = BeamEnergy::E_LER_orig();
    oE_HER    = BeamEnergy::E_HER_orig();
    p_beam    = BeamEnergy::p_beam();
    p_cm      = BeamEnergy::p_cm(  p_beam );
    Ecm       = BeamEnergy::Ecm();

    //...........................................
    // Set status code to select candidate event.
    //..........................................
    // *status = -1;     //   Discard if no candidate.
    ievent++;
    Belle_event_Manager& evmgr = Belle_event_Manager::get_manager();
    Belle_event& evthead = *evmgr.begin();
    int expno = evthead.ExpNo();
    int runno = evthead.RunNo();
    int evtno = (evthead.EvtNo() & 0x0fffffff);
    int farmno =  (evthead.EvtNo() >> 28 ); 
    
    Mdst_pi0_Manager&   Pi0Mgr  =    Mdst_pi0_Manager::get_manager();
    Mdst_vee2_Manager& Vee2Mgr= Mdst_vee2_Manager::get_manager(); // Ks
    Mdst_charged_Manager & ChgMgr= Mdst_charged_Manager::get_manager();
    Mdst_gamma_Manager & GamMgr = Mdst_gamma_Manager::get_manager();
    Mdst_ecl_aux_Manager & eclaux = Mdst_ecl_aux_Manager::get_manager();
    Mdst_ecl_Manager & EclMgr = Mdst_ecl_Manager::get_manager();
    Mdst_trk_Manager& mdsttrkmgr = Mdst_trk_Manager::get_manager();
    Mdst_trk_fit_Manager& mdsttrkfitmgr = Mdst_trk_fit_Manager::get_manager();
    Evtcls_hadron_info_Manager &HadMgr = Evtcls_hadron_info_Manager::get_manager();
    Belle_event_Manager& EveMgr = Belle_event_Manager::get_manager();
    Belle_event& Evt = *EveMgr.begin();
    Evtcls_hadron_info& Had = *HadMgr.begin();

    if(HadMgr.size() >0){
      int ntrk  = Had.Ntrk();
    }//HadMgr size loop
    
    std::vector<Particle> pion_pm;
    std::vector<Particle> pi_minus;
    std::vector<Particle> pi_plus;
    std::vector<Particle> kaon_pm;
    std::vector<Particle> K_plus;
    std::vector<Particle> Dpm;
    std::vector<Particle> Dstr;
    std::vector<Particle> Pi0; 
    std::vector<Particle> Gamma;
    std::vector<Particle> spion;
    std::vector<Particle> Ks;//Added for Ks
    std::vector<Particle> D0;//Added for Ks
    std::vector<Particle> Dp;
    std::vector<Particle> electron;
    std::vector<Particle> positron;	
    std::vector<Particle> proton;    
    std::vector<Particle> anti_proton;    
    std::vector<Particle> nD;
    std::vector<Particle> muon_minus;
    std::vector<Particle> muon_plus;
    

      Ptype ptype_positron("E+");
      Ptype ptype_electron("E-");
      Ptype ptype_proton("P+");
      Ptype ptype_antiproton("AP+");
      Ptype ptype_D0("D0");
      Ptype ptype_Dp("D+");
      Ptype ptype_Dm("D-");
      Ptype ptype_Gamma("GAMM");
      Ptype ptype_kaon_plus("K+");
      Ptype ptype_kaon_minus("K-");
      Ptype ptype_pi_plus("PI+");
      Ptype ptype_pi_minus("PI-");
      Ptype ptype_Pi0("PI0");
      Ptype ptype_Dsp("D*+");
      Ptype ptype_Dsm("D*-");
      Ptype ptype_muon_minus("MU-");
      Ptype ptype_muon_plus("MU+");

      Hep3Vector ipProf(IpProfile::position()); 
      atc_pid sel_k (3, 1, 5, 3, 2);  //Kaon vs Pion prob
      atc_pid sel_pi(3, 1, 5, 2, 3);  //pion vs kaon prob
      atc_pid selPK ( 3, 1, 5, 4, 3);  // p/K  separation (proton to kaon)           
      atc_pid selPPI( 3, 1, 5, 4, 2); //p/pi separation (proton to pion) 	

      double atc_PK, atc_PPI;
      atc_PK=-999, atc_PPI=-999;
      
      for(Mdst_charged_Manager::iterator p = ChgMgr.begin();
	  p!=ChgMgr.end(); p++){  //<<<<<<<< Charged particle 
	Mdst_charged& Chg = *p;
	
	eid sel_e(*p);    //electron ID
        Muid_mdst sel_mu(*p); //muon ID
	Mdst_trk& Trk = Chg.trk();
	
	Hep3Vector  Pch(Chg.px(),Chg.py(),Chg.pz());
	float ecl_energy = sel_e.eoverp()*Pch.mag();
	
	Mdst_trk_fit& TrkFit0 = Trk.mhyp(0) ;//for electron
	Mdst_trk_fit& TrkFit1 = Trk.mhyp(1) ;//for muon
	Mdst_trk_fit& TrkFit2 = Trk.mhyp(2) ;//for pion
	Mdst_trk_fit& TrkFit3 = Trk.mhyp(3) ;//for kaon
	Mdst_trk_fit& TrkFit4 = Trk.mhyp(4) ;//for proton

	int nhit3 = (*p).trk().mhyp(2).nhits(3);
	int nhit4= (*p).trk().mhyp(2).nhits(4);
	int svdhit = nhit3+nhit4 ;	
	double dxk = (*p).trk().mhyp(3).pivot_x( );
	double dyk = (*p).trk().mhyp(3).pivot_y( );
	double dzk = (*p).trk().mhyp(3).pivot_z( );
	double dz_k = (*p).trk().mhyp(3).helix(3);
	double dr_k = (*p).trk().mhyp(3).helix(0);
      
	if(fabs(dz_k) < 3.0 && fabs(dr_k) < 1.0){
	 if (sel_k.prob(*p) > 0.6) {
	    Particle kaon(*p, ((*p).charge()>0) ? ptype_kaon_plus : ptype_kaon_minus, ipProf);
	    kaon.userInfo( UserInfo() );
	    dynamic_cast<UserInfo&>(kaon.userInfo()).pidpi(sel_k.prob(*p));
	    dynamic_cast<UserInfo&>(kaon.userInfo()).deltz((*p).trk().mhyp(3).helix(3));
	    dynamic_cast<UserInfo&>(kaon.userInfo()).rphiv((*p).trk().mhyp(3).helix(0));
	    dynamic_cast<UserInfo&>(kaon.userInfo()).phi0((*p).trk().mhyp(3).helix(1));
	    dynamic_cast<UserInfo&>(kaon.userInfo()).kappa((*p).trk().mhyp(3).helix(2));
	    dynamic_cast<UserInfo&>(kaon.userInfo()).tanld((*p).trk().mhyp(3).helix(4));
	    dynamic_cast<UserInfo&>(kaon.userInfo()).nhit(svdhit);
	    dynamic_cast<UserInfo&>(kaon.userInfo()).pivotx(dxk);
	    dynamic_cast<UserInfo&>(kaon.userInfo()).pivoty(dyk);
	    dynamic_cast<UserInfo&>(kaon.userInfo()).pivotz(dzk);
	    kaon_pm.push_back(kaon);
	  } // kaon cut
	}//KID cut

	double dz2 = (*p).trk().mhyp(2).helix(3);
	double dr2 = (*p).trk().mhyp(2).helix(0);
	if(fabs(dz2) < 3.0 && fabs(dr2) < 1.0){
	if(sel_pi.prob(*p)>0.6){
	  Particle pion(*p,  ((*p).charge()>0) ? ptype_pi_plus : ptype_pi_minus, ipProf);
	  pion.userInfo( UserInfo() );
	  dynamic_cast<UserInfo&>(pion.userInfo()).pidpi(sel_pi.prob(*p));
	  dynamic_cast<UserInfo&>(pion.userInfo()).deltz((*p).trk().mhyp(2).helix(3));
	  dynamic_cast<UserInfo&>(pion.userInfo()).rphiv((*p).trk().mhyp(2).helix(0));
	  dynamic_cast<UserInfo&>(pion.userInfo()).phi0((*p).trk().mhyp(2).helix(1));
	  dynamic_cast<UserInfo&>(pion.userInfo()).kappa((*p).trk().mhyp(2).helix(2));
	  dynamic_cast<UserInfo&>(pion.userInfo()).tanld((*p).trk().mhyp(2).helix(4));
	  dynamic_cast<UserInfo&>(pion.userInfo()).nhit(svdhit);
	  pion_pm.push_back(pion);      
	}//pion cut
	}//pid cut
	//}//svdhit cut
	
	double dz_e = (*p).trk().mhyp(0).helix(3);
        double dr_e = (*p).trk().mhyp(0).helix(0);
	
	if (fabs(dz_e) < 3.0  && fabs(dr_e)<1.0) {
	if(sel_e.prob(3,-1,5)>0.9){	
	if((*p).charge() < 0.0){
	  Particle e_minus(*p, ptype_electron, ipProf);
	  e_minus.userInfo( UserInfo() );
	  dynamic_cast<UserInfo&>(e_minus.userInfo()).pidpi(sel_k.prob(*p));
	  dynamic_cast<UserInfo&>(e_minus.userInfo()).deltz((*p).trk().mhyp(2).helix(3));
	  dynamic_cast<UserInfo&>(e_minus.userInfo()).rphiv((*p).trk().mhyp(2).helix(0));
	  dynamic_cast<UserInfo&>(e_minus.userInfo()).phi0((*p).trk().mhyp(2).helix(1));
	  dynamic_cast<UserInfo&>(e_minus.userInfo()).kappa((*p).trk().mhyp(2).helix(2));
	  dynamic_cast<UserInfo&>(e_minus.userInfo()).tanld((*p).trk().mhyp(2).helix(4));
	  dynamic_cast<UserInfo&>(e_minus.userInfo()).nhit(svdhit);	  
	  electron.push_back(e_minus);
	  
	}
	else{
	  Particle e_plus(*p, ptype_positron, ipProf);
	  e_plus.userInfo( UserInfo() );
	  dynamic_cast<UserInfo&>(e_plus.userInfo()).pidpi(sel_k.prob(*p));
          dynamic_cast<UserInfo&>(e_plus.userInfo()).deltz((*p).trk().mhyp(2).helix(3));
          dynamic_cast<UserInfo&>(e_plus.userInfo()).rphiv((*p).trk().mhyp(2).helix(0));
          dynamic_cast<UserInfo&>(e_plus.userInfo()).phi0((*p).trk().mhyp(2).helix(1));
          dynamic_cast<UserInfo&>(e_plus.userInfo()).kappa((*p).trk().mhyp(2).helix(2));
          dynamic_cast<UserInfo&>(e_plus.userInfo()).tanld((*p).trk().mhyp(2).helix(4));
          dynamic_cast<UserInfo&>(e_plus.userInfo()).nhit(svdhit);	  
	  positron.push_back(e_plus);
	}	
	}
	} // electron & positron
	
	if (Evt.ExpMC()==2) {  setGenHepInfoF(electron);
	  setGenHepInfoF(positron);}
	
	//----------muon---------                                                                                                                     
	double dz_mu = (*p).trk().mhyp(1).helix(3);
	double dr_mu = (*p).trk().mhyp(1).helix(0);
        //std::cout<<"dz_mu "<<dz_mu<<std::endl;                                                                                                      
        if (fabs(dz_mu) <3.0 && fabs(dr_mu)<1.0) {
	  if(sel_mu.Muon_likelihood() > 0.9){
	    if((*p).charge() < 0.0){
	      Particle mu_minus(*p, ptype_muon_minus, ipProf);
	      mu_minus.userInfo( UserInfo() );
	      dynamic_cast<UserInfo&>(mu_minus.userInfo()).pidpi(sel_k.prob(*p));
	      dynamic_cast<UserInfo&>(mu_minus.userInfo()).deltz((*p).trk().mhyp(1).helix(3));
	      dynamic_cast<UserInfo&>(mu_minus.userInfo()).rphiv((*p).trk().mhyp(1).helix(0));
	      dynamic_cast<UserInfo&>(mu_minus.userInfo()).phi0((*p).trk().mhyp(1).helix(1));
	      dynamic_cast<UserInfo&>(mu_minus.userInfo()).kappa((*p).trk().mhyp(1).helix(2));
	      dynamic_cast<UserInfo&>(mu_minus.userInfo()).tanld((*p).trk().mhyp(1).helix(4));
	      dynamic_cast<UserInfo&>(mu_minus.userInfo()).nhit(svdhit);
	      muon_minus.push_back(mu_minus);
	    }else {
	      Particle mu_plus(*p, ptype_muon_plus, ipProf);
	      mu_plus.userInfo( UserInfo() );
	      dynamic_cast<UserInfo&>(mu_plus.userInfo()).pidpi(sel_k.prob(*p));
              dynamic_cast<UserInfo&>(mu_plus.userInfo()).deltz((*p).trk().mhyp(1).helix(3));
              dynamic_cast<UserInfo&>(mu_plus.userInfo()).rphiv((*p).trk().mhyp(1).helix(0));
              dynamic_cast<UserInfo&>(mu_plus.userInfo()).phi0((*p).trk().mhyp(1).helix(1));
              dynamic_cast<UserInfo&>(mu_plus.userInfo()).kappa((*p).trk().mhyp(1).helix(2));
              dynamic_cast<UserInfo&>(mu_plus.userInfo()).tanld((*p).trk().mhyp(1).helix(4));
              dynamic_cast<UserInfo&>(mu_plus.userInfo()).nhit(svdhit);
	      muon_plus.push_back(mu_plus);
	    }
	    if (Evt.ExpMC()==2) {       setGenHepInfoF(muon_minus);
              setGenHepInfoF(muon_plus);}
          }
        }
	
	
	atc_PK  = selPK.prob(Chg);
        atc_PPI = selPPI.prob(Chg);		
	double dz_p = (*p).trk().mhyp(4).helix(3);
        double dr_p = (*p).trk().mhyp(4).helix(0);
	
	if (fabs(dz_p) <3.0 && fabs(dr_p)<1.0){
	if((atc_PK > 0.6)&&(atc_PPI > 0.6)){
	  if((*p).charge()==+1){	
	    Particle p_plus(*p,ptype_proton,ipProf);
	    p_plus.userInfo( UserInfo() );
	    dynamic_cast<UserInfo&>(p_plus.userInfo()).pidpi(sel_k.prob(*p));
	    dynamic_cast<UserInfo&>(p_plus.userInfo()).deltz((*p).trk().mhyp(2).helix(3));
	    dynamic_cast<UserInfo&>(p_plus.userInfo()).rphiv((*p).trk().mhyp(2).helix(0));
	    dynamic_cast<UserInfo&>(p_plus.userInfo()).phi0((*p).trk().mhyp(2).helix(1));
	    dynamic_cast<UserInfo&>(p_plus.userInfo()).kappa((*p).trk().mhyp(2).helix(2));
	    dynamic_cast<UserInfo&>(p_plus.userInfo()).tanld((*p).trk().mhyp(2).helix(4));
	    dynamic_cast<UserInfo&>(p_plus.userInfo()).nhit(svdhit);
	    proton.push_back(p_plus);
	  } // proton
	  else{
	    Particle p_minus(*p,ptype_antiproton,ipProf);
            p_minus.userInfo( UserInfo() );
	    dynamic_cast<UserInfo&>(p_minus.userInfo()).pidpi(sel_k.prob(*p));
	    dynamic_cast<UserInfo&>(p_minus.userInfo()).deltz((*p).trk().mhyp(2).helix(3));
            dynamic_cast<UserInfo&>(p_minus.userInfo()).rphiv((*p).trk().mhyp(2).helix(0));
            dynamic_cast<UserInfo&>(p_minus.userInfo()).phi0((*p).trk().mhyp(2).helix(1));
            dynamic_cast<UserInfo&>(p_minus.userInfo()).kappa((*p).trk().mhyp(2).helix(2));
            dynamic_cast<UserInfo&>(p_minus.userInfo()).tanld((*p).trk().mhyp(2).helix(4));
            dynamic_cast<UserInfo&>(p_minus.userInfo()).nhit(svdhit);
	    anti_proton.push_back(p_minus);
	  }
	  if (Evt.ExpMC()==2) {  setGenHepInfoF(proton);
	    setGenHepInfoF(anti_proton);
	  } 
	}	
	}
      }//chgMgr ends here

      //-------------------------Tagging------------------------

      int mode_dp, mode_dm, mode_dz, mode_adz;
      mode_dp=mode_dm=mode_dz=mode_adz=-999;
      int tdp, tdm, tdz, tadz;
      tdp=tdm=tdz=tadz=0;
      if (Evt.ExpMC()==2){
	int abzc=0; int bzc=0;
	Gen_hepevt_Manager & hepevt_mag = Gen_hepevt_Manager::get_manager();
	for(std::vector<Gen_hepevt>::iterator it = hepevt_mag.begin();
	    it != hepevt_mag.end(); ++it){
	  Gen_hepevt& particle = *it;

	  if (particle.idhep()==-411){
	    tdp =tdp+1;
	    int modedp =Gen.Dmode(particle);
	    if (modedp<0) modedp=999;
	    if (tdp==1) mode_dp=modedp;
	    if (tdp ==2) mode_dp=modedp*1000+mode_dp;
	    if (tdp ==3) mode_dp=modedp*1000000+mode_dp;;
	  }

	  if (particle.idhep()==411){
	    tdm=tdm+1;
	    int modedm=Gen.Dmode(particle);
	    if (modedm<0) modedm=999;
	    if (tdm==1) mode_dm=modedm;
	    if (tdm==2) mode_dm=modedm*1000+mode_dm;
	    if (tdm==3) mode_dm=modedm*1000000+mode_dm;
	  }
	  if (particle.idhep()==-421){
	    tadz=tadz+1;
	    int modeadz=Gen.Dmode(particle);
	    if (modeadz<0) modeadz=999;
	    if (tadz==1) mode_adz=modeadz;
	    if (tadz==2) mode_adz=modeadz*1000+mode_adz;
	    if (tadz==3) mode_adz=modeadz*1000000+mode_adz;
	  }
	  if (particle.idhep()==421){
	    tdz=tdz+1;
	    int modedz=Gen.Dmode(particle);
	    //std::cout<<"modedz "<<modedz<<std::endl;
	    if (modedz<0) modedz=999;
	    if (tdz==1) mode_dz=modedz;
	    if (tdz==2) mode_dz=modedz*1000+mode_dz;
	    if (tdz==3) mode_dz=modedz*1000000+mode_dz;
	  }
	}
      }
      //std::cout<<"hello"<<std::endl;
	//---------------------GAMMA CORRECTION-------------------
	for(Mdst_gamma_Manager::iterator q = GamMgr.begin();
	    q!=GamMgr.end(); q++){
  	Mdst_gamma& gam = *q;
  	Mdst_ecl& ecl = gam.ecl();
 	 double e9ovr = eclaux[ecl.get_ID()-1].e9oe25();
  	int eclid = gam.ecl().get_ID();
  	int tinfo1 = Get.GTinfo(eclid,1);
  	int tinfo2 = Get.GTinfo(eclid,2);
  	int isthep = 55555;
  	if (Evt.ExpMC()==2) {isthep = Get.Gisthep(eclid);}
	//Here we check the unassociated track and their quality                                                                                                     
	if(ecl.match() == 0 && ecl.quality() == 0 && e9ovr > 0.85){
	  double E = ecl.energy(), Theta = ecl.theta(), Phi=ecl.phi();
	  HepLorentzVector Pgam(E*sin(Theta)*cos(Phi),
				E*sin(Theta)*sin(Phi),
				E*cos(Theta), E);
	  HepSymMatrix errEcl(3,0); //error in ecl matrix on E, theta, phi              
	  
	  errEcl[0][0] = ecl.error(0); //Energy                                                                                                                        
	  errEcl[1][0] = ecl.error(1);
	  errEcl[1][1] = ecl.error(2); //Phi                                                                                                                           
          errEcl[2][0] = ecl.error(3);
          errEcl[2][1] = ecl.error(4);
          errEcl[2][2] = ecl.error(5); //Theta         

	HepSymMatrix errCart(4,0);
          HepMatrix jacobian(4,3,0);

          jacobian[0][0] =         cos(Phi)*sin(Theta);
          jacobian[0][1] =      -E*sin(Phi)*sin(Theta);
          jacobian[0][2] =       E*cos(Phi)*cos(Theta);
          jacobian[1][0] =         sin(Phi)*sin(Theta);
          jacobian[1][1] =       E*cos(Phi)*sin(Theta);
          jacobian[1][2] =       E*sin(Phi)*cos(Theta);
          jacobian[2][0] =                  cos(Theta);
          jacobian[2][1] =                         0.0;
          jacobian[2][2] =               -E*sin(Theta);
          jacobian[3][0] =                         1.0;
          jacobian[3][1] =                         0.0;
          jacobian[3][2] =                         0.0;

	errCart = errEcl.similarity(jacobian);
          HepSymMatrix dX(IpProfile::position_err()); Hep3Vector X(ipProf);
          //position and error matrix    

	Particle TmpGam(*q); TmpGam.name("Gamma");
          TmpGam.momentum().momentum(Pgam,errCart);  //set momenta                                                                                                     
          TmpGam.momentum().position(X,dX);

          int gid = (*q).get_ID();


          TmpGam.userInfo( UserInfo() );
          dynamic_cast<UserInfo&>(TmpGam.userInfo()).e9ovr(e9ovr);
          dynamic_cast<UserInfo&>(TmpGam.userInfo()).gid(gid);
          dynamic_cast<UserInfo&>(TmpGam.userInfo()).tinfo1(tinfo1);
          dynamic_cast<UserInfo&>(TmpGam.userInfo()).tinfo2(tinfo2);
          if (Evt.ExpMC()==2) {
            int gma1 = Get.GMiD(eclid,1);
            int gma2 = Get.GMiD(eclid,2);
            int gma3 = Get.GMiD(eclid,3);
            int gma4 = Get.GMiD(eclid,4);
	    //std::cout<<gma1<<gma2<<std::endl;
            dynamic_cast<UserInfo&>(TmpGam.userInfo()).isthep(isthep);
            dynamic_cast<UserInfo&>(TmpGam.userInfo()).gma1(gma1);
            dynamic_cast<UserInfo&>(TmpGam.userInfo()).gma2(gma2);
            dynamic_cast<UserInfo&>(TmpGam.userInfo()).gma3(gma3);
            dynamic_cast<UserInfo&>(TmpGam.userInfo()).gma4(gma4);
	}
          Gamma.push_back(TmpGam);
        }
      }//mdst gamma     
	if (Evt.ExpMC()==2) { setGenHepInfoG(Gamma);}
	
	//std::cout<<"hello1"<<std::endl;
	//_____________________________________________________________                                                                              
	// reconstruction of  D0 from e- p+              
	//---------------------------------------------------------------------------                                      
	for( std::vector<Particle>::iterator i = electron.begin() ;i != electron.end(); i++ ) {
	  for( std::vector<Particle>::iterator j =  proton.begin(); j != proton.end(); j++ ) { 
	    if(checkSame(*i,*j))continue;
	    HepLorentzVector pEle((*i).momentum().p());
	    HepSymMatrix     pEleEr((*i).momentum().dp());
	    HepLorentzVector pEle_keep=((*i).momentum().p());
	    
	    HepLorentzVector prPos((*j).momentum().p());
	    HepSymMatrix     prPosEr((*j).momentum().dp());
	    HepLorentzVector pPro_keep=((*j).momentum().p());
	    
	    HepLorentzVector P_D = (pEle+prPos);
	    Particle d0_cand((*i).momentum().p()+(*j).momentum().p(), ptype_D0);
	    //std::cout<<"Hi "<<std::endl;
	    d0_cand.relation().append(*i);
	    d0_cand.relation().append(*j);

	    double confLevel = -9.9;
	    
	    int vtxflg,mcflg;
	    vtxflg = Vf.VXfit(d0_cand,confLevel,0);
	    double vt_x= Vf.VX().x();
	    double vt_y= Vf.VX().y();
	    double vt_z = Vf.VX().z();
	    double v_chisq = Vf.chisq();
	    
	    if (vtxflg !=0) { 	      
	      HepLorentzVector Pgam;
	      for(std::vector<Particle>::iterator it = Gamma.begin();
	      it!=Gamma.end(); it++){
		Pgam = HepLorentzVector((*it).momentum().p());
		double ang_em = pEle_keep.vect().angle(Pgam.vect());
		double momg = Pgam.vect().mag();
		if( pEle_keep.vect().angle(Pgam.vect())<0.05 || 
		    pPro_keep.vect().angle(Pgam.vect())< 0.05){
		  P_D+=Pgam;
		  double brem_px = Pgam.px();
		  d0_cand.relation().append(*it);
		  (*it).name((*it).name() + "Rad ");
		}
	      }// Gamma loop 
	      
	      d0_cand.momentum(P_D);
	      
	      double md0   = P_D.mag();
	      double d0mom = P_D.vect().mag();

	      if(md0>1.85 && md0<1.875)
		{ // 0.2 of PDG D0 mass
		  d0_cand.userInfo(UserInfo());
		  //vertexFit(d0_cand);	      
		dynamic_cast<UserInfo&>(d0_cand.userInfo()).mom(P_D.vect().mag());
		dynamic_cast<UserInfo&>(d0_cand.userInfo()).mass(P_D.mag());
		dynamic_cast<UserInfo&>(d0_cand.userInfo()).etag(1);
		dynamic_cast<UserInfo&>(d0_cand.userInfo()).chisq();
		double sum_chisq = dynamic_cast<UserInfo&>(d0_cand.userInfo()).chisq();
		D0.push_back(d0_cand);
		}
	    }//md0 cut
	  }
	}
      	
	//std::cout<<"hello2"<<std::endl;
	//_____________________________________________________________                                                                              
	// reconstruction of  D0 from e+ p-bar              
	//---------------------------------------------------------------------------                                      
	for( std::vector<Particle>::iterator i = positron.begin() ;i != positron.end(); i++ ) {
	  for( std::vector<Particle>::iterator j =  anti_proton.begin(); j != anti_proton.end(); j++ ) { 
	    if(checkSame(*i,*j))continue;
	    HepLorentzVector pEle((*i).momentum().p());
	    HepSymMatrix     pEleEr((*i).momentum().dp());
	    HepLorentzVector pEle_keep=((*i).momentum().p());
	    
	    HepLorentzVector prPos((*j).momentum().p());
	    HepSymMatrix     prPosEr((*j).momentum().dp());
	    HepLorentzVector pPro_keep=((*j).momentum().p());
	    
	    HepLorentzVector P_D = (pEle+prPos);
	    Particle d0_cand((*i).momentum().p()+(*j).momentum().p(), ptype_D0);
	    //std::cout<<"Hi "<<std::endl;
	    d0_cand.relation().append(*i);
	    d0_cand.relation().append(*j);

	    double confLevel = -9.9;
	    
	    int vtxflg,mcflg;
	    vtxflg = Vf.VXfit(d0_cand,confLevel,0);
	    double vt_x= Vf.VX().x();
	    double vt_y= Vf.VX().y();
	    double vt_z = Vf.VX().z();
	    double v_chisq = Vf.chisq();
	    
	    if (vtxflg !=0) { 
	      
	      HepLorentzVector Pgam;
	      
	      for(std::vector<Particle>::iterator it = Gamma.begin();
	      it!=Gamma.end(); it++){
		Pgam = HepLorentzVector((*it).momentum().p());
		double ang_em = pEle_keep.vect().angle(Pgam.vect());
		double momg = Pgam.vect().mag();
		if( pEle_keep.vect().angle(Pgam.vect())<0.05 || 
		    pPro_keep.vect().angle(Pgam.vect())< 0.05){
		  
		  P_D+=Pgam;
		  double brem_px = Pgam.px();
		  d0_cand.relation().append(*it);
		  (*it).name((*it).name() + "Rad ");
		}
	      }// Gamma loop 
	      
	      d0_cand.momentum(P_D);

	      double md0   = P_D.mag();
	      double d0mom = P_D.vect().mag();
	      
	      
	      if(md0>1.85 && md0<1.875){ // 0.2 of PDG D0 mass
		
		d0_cand.userInfo( UserInfo() );
		dynamic_cast<UserInfo&>(d0_cand.userInfo()).mom(P_D.vect().mag());
		dynamic_cast<UserInfo&>(d0_cand.userInfo()).mass(P_D.mag());
		dynamic_cast<UserInfo&>(d0_cand.userInfo()).etag(2);
		//vertexFit(d0_cand);
		dynamic_cast<UserInfo&>(d0_cand.userInfo()).chisq();
		double sum_chisq = dynamic_cast<UserInfo&>(d0_cand.userInfo()).chisq();
		
		D0.push_back(d0_cand);
	      }
	    }//md0 cut
	  }
	}
		
	
	//-----------reconstruction of D0 from mu- p+---------                                                                                        

        for(std::vector<Particle>::iterator i = muon_minus.begin(); i !=muon_minus.end(); i++) {
          for(std::vector<Particle>::iterator j = proton.begin(); j != proton.end(); j++){
	    
	    if(checkSame(*i,*j))continue;
            Particle d0_cand((*i).momentum().p()+(*j).momentum().p(), ptype_D0);
            d0_cand.relation().append(*i);
            d0_cand.relation().append(*j);
            
	    HepLorentzVector P_D = ((*i).momentum().p() + (*j).momentum().p());

            double md0   = P_D.mag();
            double d0mom = P_D.vect().mag();
            
	    if(md0>1.85 && md0<1.875){ // 0.2 of PDG D0 mass                                                                                             
	      d0_cand.userInfo( UserInfo() );
              dynamic_cast<UserInfo&>(d0_cand.userInfo()).mom(P_D.vect().mag());
              dynamic_cast<UserInfo&>(d0_cand.userInfo()).mass(P_D.mag());
	      dynamic_cast<UserInfo&>(d0_cand.userInfo()).etag(101);
              //vertexFit(d0_cand);
              dynamic_cast<UserInfo&>(d0_cand.userInfo()).chisq();
              double sum_chisq = dynamic_cast<UserInfo&>(d0_cand.userInfo()).chisq();
	      D0.push_back(d0_cand);
            }
          }
        }


	//-----------reconstruction of D0 from mu+ pbar ---------                                                                  

        for(std::vector<Particle>::iterator i = muon_plus.begin(); i !=muon_plus.end(); i++) {
          for(std::vector<Particle>::iterator j = anti_proton.begin(); j != anti_proton.end(); j++){
	    
	    if(checkSame(*i,*j))continue;
            Particle d0_cand((*i).momentum().p()+(*j).momentum().p(), ptype_D0);
            d0_cand.relation().append(*i);
            d0_cand.relation().append(*j);
            
	    HepLorentzVector P_D = ((*i).momentum().p() + (*j).momentum().p());
	    
            double md0   = P_D.mag();
            double d0mom = P_D.vect().mag();
            
	    if(md0>1.85 && md0<1.875){ // 0.2 of PDG D0 mass                                      

	      d0_cand.userInfo( UserInfo() );
              dynamic_cast<UserInfo&>(d0_cand.userInfo()).mom(P_D.vect().mag());
              dynamic_cast<UserInfo&>(d0_cand.userInfo()).mass(P_D.mag());
	      dynamic_cast<UserInfo&>(d0_cand.userInfo()).etag(102);
              //vertexFit(d0_cand);
              dynamic_cast<UserInfo&>(d0_cand.userInfo()).chisq();
              double sum_chisq = dynamic_cast<UserInfo&>(d0_cand.userInfo()).chisq();
	      D0.push_back(d0_cand);
            }
          }
        }
	

	//-----------reconstruction of D0 from KK---------- 
	for(std::vector<Particle>::iterator i = kaon_pm.begin(); i !=kaon_pm.end(); i++) {
          for(std::vector<Particle>::iterator j = kaon_pm.begin(); j != kaon_pm.end(); j++){

            if(j<=i)continue;
            if(checkSame(*i,*j))continue;
            Particle d0_cand((*i).momentum().p()+(*j).momentum().p(), ptype_D0);
            d0_cand.relation().append(*i);
            d0_cand.relation().append(*j);

            HepLorentzVector P_D = ((*i).momentum().p() + (*j).momentum().p());

            double md0   = P_D.mag();
            double d0mom = P_D.vect().mag();

            if(md0>1.85 && md0<1.875){ // 0.2 of PDG D0 mass 
	      d0_cand.userInfo( UserInfo() );

              dynamic_cast<UserInfo&>(d0_cand.userInfo()).mom(P_D.vect().mag());
              dynamic_cast<UserInfo&>(d0_cand.userInfo()).mass(P_D.mag());
              dynamic_cast<UserInfo&>(d0_cand.userInfo()).etag(103);
              ///vertexFit(d0_cand);
              dynamic_cast<UserInfo&>(d0_cand.userInfo()).chisq();
              double sum_chisq = dynamic_cast<UserInfo&>(d0_cand.userInfo()).chisq();
              D0.push_back(d0_cand);
            }
          }
        }
		
	if(Evt.ExpMC()==2){setGenHepInfoR(D0);}
	//std::cout<<D0.size()<<std::endl;
	//std::cout<<"hello5"<<std::endl;	

	//_____________________________________________________________                                                                           
	// reconstruction of  D*+ from D0 pi+                                                                                  
	//_____________________________________________________________                                       
                                   
	for(std::vector<Particle>::iterator i = D0.begin(); i !=D0.end(); i++) {
	  for(std::vector<Particle>::iterator j = pion_pm.begin(); j != pion_pm.end(); j++){
	    //if(j<=i)continue;
	    if (checkSame(*i,*j)) continue;
	    
	    Particle dsm_cand((*i).momentum().p()+(*j).momentum().p(), ((*j).charge()>0) ? ptype_Dsp : ptype_Dsm);
	    dsm_cand.relation().append(*i);
	    dsm_cand.relation().append(*j);

	    HepLorentzVector dsm = ((*i).momentum().p() + (*j).momentum().p());
	    double mass =dsm.mag();
	    HepLorentzVector dsm_cm = dsm;
	    //std::cout<<"Hey I am D*+"<<std::endl;
	    dsm_cm.boost(CMBoost);
	     
	    //std::cout<<"Hey I am vertex fit"<<std::endl;
	    // Cuts to reduce the background from direct B
	    if (dsm_cm.vect().mag() < 2.5) continue;

	    Particle &d0 = dsm_cand.child(0);
	    double deltam = mass -dynamic_cast<UserInfo&>((dsm_cand).child(0).userInfo()).mass(); 
	    double etag=dynamic_cast<UserInfo&>((dsm_cand).child(0).userInfo()).etag();
	    if(deltam>0.144 && deltam<0.147){
	    dsm_cand.userInfo( UserInfo() ); 
	    vertexFit(dsm_cand);
	    //dynamic_cast<UserInfo&>(dsm_cand.userInfo()).chisq(chi2);
	    //}
	    Dstr.push_back(dsm_cand);
	    }
	  mnt_dst->column("mass", mass);
	  mnt_dst->column("deltam",deltam);
	  mnt_dst->column("etag",etag);

	  mnt_dst->dumpData();
	     
	  }
	}
	if(Evt.ExpMC()==2){setGenHepInfoR(Dstr);}
	//std::cout<<"Dstr size is "<<Dstr.size()<<std::endl;
	//std::cout<<"hello4"<<std::endl;

	//___________________D* loop____________________
	
	int mul[10]={0,0,0,0,0,0,0,0,0,0} ;
	int v_BCS[10], d0_BCS[10] ;
	int order[10]={1,1,2,2,101,101,102,102,103,103};

	for (int j=0;j<10;j++){

	double v_temp=999999999;
	double v_d0 =999999999;
	v_BCS[j]=-100;
	//std::cout<<Dstr.size()<<std::endl;
	for(int k=0; k< int(Dstr.size()); k++)
          {
	   
            Particle Dee = Dstr[k];
	    int etag=int(dynamic_cast<UserInfo&>(Dee.child(0).userInfo()).etag());
	    //std::cout<<"tag number is "<<etag<<std::endl;
	    //if(deltam>0.135 && deltam<0.16){
	    HepLorentzVector dp_kp = (Dee.child(0).momentum().p()+
                                      Dee.child(1).momentum().p());
            double d0mass=dynamic_cast<UserInfo&>(Dee.child(0).userInfo()).mass();
            double mass = dp_kp.mag();
            double q=mass-d0mass-Dee.child(1).momentum().p().mag();
            double deltam= mass -d0mass;

	    //std::cout<<etag<<"    "<<order[j]<<std::endl;
	    if(etag != order[j]) continue;
	    double dstr_chisq = dynamic_cast<UserInfo&>(Dee.userInfo()).chisq();
	    double d0_chisq = dynamic_cast<UserInfo&>(Dee.child(0).userInfo()).chisq();
	    if(j==0 && Dee.child(1).charge()<0) continue;
	    if(j==1 && Dee.child(1).charge()>0) continue;
	    if(j==2 && Dee.child(1).charge()<0) continue;
	    if(j==3 && Dee.child(1).charge()>0) continue;
	    if(j==4 && Dee.child(1).charge()<0) continue;
	    if(j==5 && Dee.child(1).charge()>0) continue;
	    if(j==6 && Dee.child(1).charge()<0) continue;
	    if(j==7 && Dee.child(1).charge()>0) continue;
	    if(j==8 && Dee.child(1).charge()<0) continue;
	    if(j==9 && Dee.child(1).charge()>0) continue;

	    compare(dstr_chisq,k,v_temp,v_BCS[j]);
	    compare(d0_chisq,k,v_d0,d0_BCS[j]);
	    if(d0mass>1.85&&d0mass<1.875&&deltam>0.144&&deltam<0.147&&q<0.02)
	      mul[j]=mul[j]+1;
	  }

	}
	//std::cout<<"hello renu"<<std::endl;
	//std::cout<<"multiplicity "<<mul[0]<<std::endl;
	mnt_ds->column("mul0",mul[0]);
	mnt_ds->column("mul1",mul[1]);
	mnt_ds->column("mul2",mul[2]);
	mnt_ds->column("mul3",mul[3]);
	mnt_ds->column("mul4",mul[4]);
	mnt_ds->column("mul5",mul[5]);
	mnt_ds->column("mul6",mul[6]);
	mnt_ds->column("mul7",mul[7]);
	mnt_ds->column("mul8",mul[8]);
	mnt_ds->column("mul9",mul[9]);

	for(int k=0; k< int(Dstr.size()); k++)
	  {
	    Particle Dee = Dstr[k];
	    Particle d0 = Dee.child(0);
	    Particle pion_slow = Dee.child(1);
	
	    HepLorentzVector dp_kp = (Dee.child(0).momentum().p()+ 
				      Dee.child(1).momentum().p());
	    double d0mass=dynamic_cast<UserInfo&>(Dee.child(0).userInfo()).mass();
	    //double pimass=Dee.child(1).momentum().p().mag();//dynamic_cast<UserInfo&>(Dee.child(1).userInfo()).mass();
	    //std::cout<<Dee.child(1).momentum().p().mag()<<std::endl;
	    double mass = dp_kp.mag();
	    double d0mom=dynamic_cast<UserInfo&>(Dee.child(0).userInfo()).mom();
	    double q=mass-d0mass-Dee.child(1).momentum().p().mag();
	    double deltam= mass -d0mass;
	    double elec_mom = (Dee.child(0).child(0)).momentum().p().vect().mag();
	    double pr_mom = (Dee.child(0).child(1)).momentum().p().vect().mag();
	    double elec_mass=(Dee.child(0).child(0)).momentum().p().mag();	
	    double pr_mass=(Dee.child(0).child(1)).momentum().p().mag();	
	    double elec_th=(Dee.child(0).child(0)).momentum().p().theta(); // Angle between e- and beam axis.
	    double pr_th=(Dee.child(0).child(1)).momentum().p().theta();
	    double elec_th1=(Dee.child(0).child(0)).momentum().p().theta()* 180. * (7/22.); // Angle between e- and beam axis.
	    double pr_th1=(Dee.child(0).child(1)).momentum().p().theta()* 180. * (7/22.);
	    double el_dedx=(Dee.child(0).child(0).mdstCharged().trk().dEdx());
	    double pr_dedx=(Dee.child(0).child(1).mdstCharged().trk().dEdx());
	    double elec_theta1=dynamic_cast<UserInfo&>(Dee.child(0).child(0).userInfo()).theta();	
	    double d0_chisq = dynamic_cast<UserInfo&>(Dee.child(0).userInfo()).chisq();
	    double dstr_chisq = dynamic_cast<UserInfo&>(Dee.userInfo()).chisq();
	    int etag=int(dynamic_cast<UserInfo&>(Dee.child(0).userInfo()).etag());
	    

	    int charge_kaon1= (Dee.child(0).child(0)).charge();
	    int charge_kaon2= (Dee.child(0).child(1)).charge();
	    int bC[10], bd0[10];
	    
	    double k1_mom,k2_mom;
	    
	    if(etag==103){
	      if(charge_kaon1==+1)
		k1_mom = (Dee.child(0).child(0)).momentum().p().vect().mag();
	      if(charge_kaon2==+1)
		k1_mom = (Dee.child(0).child(1)).momentum().p().vect().mag();
	      if(charge_kaon1==-1)
		k2_mom = (Dee.child(0).child(0)).momentum().p().vect().mag();
	      if(charge_kaon2==-1)
		k2_mom = (Dee.child(0).child(1)).momentum().p().vect().mag();
	      
	    }
	    //std::cout<<" charge of the kaon "<<kaon1_mom<<std::endl; 
	    
	    for(int j=0;j<10;j++)
	      {
		bC[j]=bd0[j]=-100;
		if(etag==order[j])
		  {
		    if(v_BCS[j]==k) bC[j]=100;
		    if(d0_BCS[j]==k) bd0[j]=100;
		  }
	      }
	    //std::cout<<"filling"<<std::endl;
	    mnt_ds->column("k2_mom", k2_mom); //  saving variables                                                                                                      
            mnt_ds->column("k1_mom", k1_mom); //  saving variables                                                                  
	    
	    mnt_ds->column("mass", mass); //  saving variables 
	    mnt_ds->column("etag",etag);
	    mnt_ds->column("bc0",bC[0]);
	    mnt_ds->column("bc1",bC[1]);
	    mnt_ds->column("bc2",bC[2]);
	    mnt_ds->column("bc3",bC[3]);
	    mnt_ds->column("bc4",bC[4]);
	    mnt_ds->column("bc5",bC[5]);
	    mnt_ds->column("bc6",bC[6]);
	    mnt_ds->column("bc7",bC[7]);
	    mnt_ds->column("bc8",bC[8]);
	    mnt_ds->column("bc9",bC[9]);
	    mnt_ds->column("bd00",bd0[0]);
	    mnt_ds->column("bd01",bd0[1]);
	    mnt_ds->column("bd02",bd0[2]);
	    mnt_ds->column("bd03",bd0[3]);
	    mnt_ds->column("bd04",bd0[4]);
	    mnt_ds->column("bd05",bd0[5]);
	    mnt_ds->column("bd06",bd0[6]);
	    mnt_ds->column("bd07",bd0[7]);
	    mnt_ds->column("bd08",bd0[8]);
	    mnt_ds->column("bd09",bd0[9]);
	    mnt_ds->column("mom", d0mom);
	    mnt_ds->column("d0mass",d0mass);
	    mnt_ds->column("deltam",deltam);     
	    //mnt_ds->column("dschisq",dstr_chisq);
	    mnt_ds->column("q",q);
	    mnt_ds->column("elec_mom",elec_mom);
	    mnt_ds->column("pr_mom",pr_mom);
	    mnt_ds->column("elec_mass",elec_mass);
	    mnt_ds->column("pr_mass",pr_mass);
	    mnt_ds->column("elec_th",elec_th);
	    mnt_ds->column("pr_th",pr_th);
	    mnt_ds->column("elec_th1",elec_th1);
	    mnt_ds->column("pr_th1",pr_th1);
	    mnt_ds->column("pis_cos",Dee.child(1).momentum().p().cosTheta());
	    mnt_ds->column("elec_cos",Dee.child(0).child(0).momentum().p().cosTheta());
	    mnt_ds->column("pr_cos",Dee.child(0).child(1).momentum().p().cosTheta());
	    
	  
	  double dr = dynamic_cast<UserInfo&>(Dee.child(1).userInfo()).rphiv();
	  double dz = dynamic_cast<UserInfo&>(Dee.child(1).userInfo()).deltz();
	  double kid=dynamic_cast<UserInfo&>(Dee.child(1).userInfo()).pidpi();
	  double elec_id=dynamic_cast<UserInfo&>(Dee.child(0).child(0).userInfo()).pidpi();		
	  double pr_id=dynamic_cast<UserInfo&>(Dee.child(0).child(1).userInfo()).pidpi();
	  double elec_dz=dynamic_cast<UserInfo&>(Dee.child(0).child(0).userInfo()).deltz();
	  double pr_dz=dynamic_cast<UserInfo&>(Dee.child(0).child(1).userInfo()).deltz();
	  double elec_dr=dynamic_cast<UserInfo&>(Dee.child(0).child(0).userInfo()).rphiv();
	  double pr_dr=dynamic_cast<UserInfo&>(Dee.child(0).child(1).userInfo()).rphiv();
	  double pissvd = dynamic_cast<UserInfo&>(Dee.child(1).userInfo()).nhit();
	  double elsvd = dynamic_cast<UserInfo&>(Dee.child(0).child(0).userInfo()).nhit();
	  double prsvd = dynamic_cast<UserInfo&>(Dee.child(0).child(1).userInfo()).nhit();
	  double k1svdx = dynamic_cast<UserInfo&>(Dee.child(0).child(0).userInfo()).pivotx();
	  double k1svdy = dynamic_cast<UserInfo&>(Dee.child(0).child(0).userInfo()).pivoty();
	  double k1svdz = dynamic_cast<UserInfo&>(Dee.child(0).child(0).userInfo()).pivotz();
	  double k2svdx = dynamic_cast<UserInfo&>(Dee.child(0).child(1).userInfo()).pivotx();
          double k2svdy = dynamic_cast<UserInfo&>(Dee.child(0).child(1).userInfo()).pivoty();
          double k2svdz = dynamic_cast<UserInfo&>(Dee.child(0).child(1).userInfo()).pivotz();
	  
	  double prk1_l(selPK.prob(((Dee.child(0)).child(1)).mdstCharged()));
	  
	  double prp1_l(selPPI.prob(((Dee.child(0)).child(1)).mdstCharged()));

	  atc_pid pi( 3, 1, 5, 2, 3);
	  double pis_p(sel_pi.prob((Dee.child(1)).mdstCharged()));

	  atc_pid pi_k( 3, 1, 5, 3, 2);
	  double pis_k(sel_k.prob((Dee.child(1)).mdstCharged()));

	  
	  eid e1l((Dee.child(0)).child(0).mdstCharged());
	  double e1_l=e1l.prob();
	  
	  Muid_mdst mul((Dee.child(0)).child(0).mdstCharged());
	  double mu_l(mul.Muon_likelihood());
	  
	  mnt_ds->column("mu_l",mu_l);
	  mnt_ds->column("prk1_l",prk1_l);
	  mnt_ds->column("prp1_l",prp1_l);
	  mnt_ds->column("e1_l",e1_l);
	  mnt_ds->column("pis_p",pis_p);
          mnt_ds->column("pis_k",pis_k);

	  //std::cout<<"Hi "<<(((*i).child(0)).child(0)).momentum().p().theta()* 180. * (7/22.)<<std::endl;
	  mnt_ds->column("pissvd", pissvd);
	  mnt_ds->column("elsvd",elsvd);
	  mnt_ds->column("prsvd",prsvd);
	  mnt_ds->column("masspis",Dee.child(1).momentum().p().mag());
	  mnt_ds->column( "pismom", Dee.child(1).p().rho());
	  mnt_ds->column("chargepis",Dee.child(1).charge());
	  mnt_ds->column("dz", dz);
	  mnt_ds->column("dr", dr);
	  mnt_ds->column("kid",kid);
	  mnt_ds->column("pr_id", pr_id);
	  mnt_ds->column("elec_id",elec_id);
	  mnt_ds->column("elec_dz",elec_dz);
	  mnt_ds->column("pr_dz",pr_dz);
	  mnt_ds->column("elec_dr",elec_dr);
	  mnt_ds->column("pr_dr",pr_dr);
	  mnt_ds->column("d0chisq",dynamic_cast<UserInfo&>(Dee.child(0).userInfo()).chisq());	
	  mnt_ds->column("dschisq",dstr_chisq);
	  mnt_ds->column("elchisq",dynamic_cast<UserInfo&>(Dee.child(0).child(0).userInfo()).chisq());	
	  mnt_ds->column("pis_dedx",Dee.child(1).mdstCharged().trk().dEdx());
	  mnt_ds->column("el_dedx", el_dedx);
	  mnt_ds->column("pr_dedx",pr_dedx);
	  mnt_ds->column("modedz",mode_dz);
	  mnt_ds->column("modeadz",mode_adz);
	  mnt_ds->column("modedp",mode_dp);
          mnt_ds->column("modedm",mode_dm);
	  //std::cout<<"charge is "<<Dee.child(1).charge()<<std::endl;
	  
	  if(Evt.ExpMC()==2){
	    if(Dee){
	      setMCtruth(Dee);
	      mnt_ds->column("dstrf",getMCtruthFlag(Dee));
	      mnt_ds->column("d0f",getMCtruthFlag(Dee.child(0)));
	      mnt_ds->column("elf",getMCtruthFlag(Dee.child(0).child(0)));
	      mnt_ds->column("prf",getMCtruthFlag(Dee.child(0).child(1)));
	      mnt_ds->column("pisf",getMCtruthFlag(Dee.child(1)));
	      mnt_ds->column("dstrid",IDhep(Dee));
	      mnt_ds->column("d0id",IDhep(Dee.relation().child(0)));
	      mnt_ds->column("elid",IDhep(Dee.relation().child(0).child(0)));
	      mnt_ds->column("prid",IDhep(Dee.relation().child(0).child(1)));
	      mnt_ds->column("pisid",IDhep(Dee.relation().child(1)));
	      //std::cout<<"id is "<<IDhep(Dee.relation().child(0))<<std::endl;
	    }
	    //mnt_ds->column("lma1",Get.MiD(Dee.child(0).child(0),1));
	    //mnt_ds->column("lma2",Get.MiD(Dee.child(0).child(0),2));
	    //mnt_ds->column("prma1",Get.MiD(Dee.child(0).child(1),1));
	    //mnt_ds->column("prma2",Get.MiD(Dee.child(0).child(1),2));
	    //mnt_ds->column("pima1",Get.MiD(Dee.child(1),1));
	      
	  }
	  //std::cout<<"dumping"<<std::endl;
	  mnt_ds->dumpData();
	  

  }//Dstr loop
  
  
}//event shape

//----------vertex fitting-------------
void user_ana_check::vertexFit(Particle &p){
  kvertexfitter kvf; //always
  
  addTrack2fit(kvf, p.relation().child(0));  //for electron                                                                                      
  addTrack2fit(kvf, p.relation().child(1)); //for proton                                                                                       
  unsigned  err = kvf.fit();
  
  if(err == 0)
    {
      dynamic_cast<UserInfo&>(p.userInfo()).cl(kvf.cl());
      dynamic_cast<UserInfo&>(p.userInfo()).chisq(kvf.chisq());
      dynamic_cast<UserInfo&>(p.userInfo()).ndf(kvf.dgf());
      makeMother(kvf,p);
    }
  
  else
    { //fail                                                                                                                      
      dynamic_cast<UserInfo&>(p.userInfo()).cl(-1.0);
      HepPoint3D vtx(999.,999.,999.);
      HepSymMatrix errVtx(3,0);
      p.momentum().decayVertex(vtx,errVtx);
    }
  
}

void user_ana_check::vertexDstrFit(Particle &p){
  kvertexfitter kvf;
  addTrack2fit(kvf, p.relation().child(0));  //for pion1                                                                                      
  addTrack2fit(kvf, p.relation().child(1));
  unsigned  err = kvf.fit(); //always as its doing fitting                                                                                     
  if(err == 0)
    {// success                                                                                                                                 
      dynamic_cast<UserInfo&>(p.userInfo()).cl(kvf.cl());
      dynamic_cast<UserInfo&>(p.userInfo()).chisq(kvf.chisq());
      dynamic_cast<UserInfo&>(p.userInfo()).ndf(kvf.dgf());
      makeMother(kvf,p);
    }
  else
    { //fail                                                                                                                      
      dynamic_cast<UserInfo&>(p.userInfo()).cl(-1.0);
      HepPoint3D vtx(999.,999.,999.);
      HepSymMatrix errVtx(3,0);
      p.momentum().decayVertex(vtx,errVtx);
    }
}

// IP constraint fit

HepPoint3D
productionPoint(Particle &p,
		const HepPoint3D &ip,
		const HepSymMatrix &errIp,
		unsigned notRewrite = 1)
{
  kvertexfitter kvf;
  addTrack2fit(kvf,p); // add "D0" to kfitter.
  addBeam2fit(kvf,ip,errIp); // add known position with error to kfitter.
  unsigned err = kvf.fit(); // do "fitting".
  if(err == 0){ // success.
    if(notRewrite != 1){
      p.momentum().momentumPosition(kvf.momentum(0), // rewrite "D0" information.
				    kvf.position(0), // 4-momentum, position and error(7x7).
				    kvf.error(0));
    }
    p.momentum().vertex(kvf.vertex(),kvf.errVertex()); // set "D0" prodution point and error.
    return kvf.vertex(); // "D0" prodution point.
  }else{ // fail.
    return HepPoint3D(999.,999.,999.);
  }
}

bool user_ana_check::refitSlowPion(Particle &p,const HepPoint3D &ip)
//bool
//refitSlowPion(Particle &p,
//		const HepPoint3D &ip)
{
  kvertexfitter kvf;
  addTrack2fit(kvf,p); // add "slowPI+" to kfitter.
  kvf.initialVertex(ip);
  kvf.knownVertex();
  unsigned err = kvf.fit(); // do "fitting".
  if(err == 0){ // success.
    p.momentum().momentumPosition(kvf.momentum(0), // rewrite "slowPI+" information.
				  kvf.position(0), // 4-momentum, position and error(7x7).
				  kvf.error(0));
    p.momentum().vertex(kvf.vertex(),kvf.errVertex()); // set "slowPI+" prodution point and error.
    return true;
  }else{ // fail.
    return false;
  }
}

//================
// Begin_run function
//================

void user_ana_check::begin_run( BelleEvent*, int* )
{
   // Initialization for hamlet.
   // eid initialization is done in hamlet beginrun procedure
  eid::init_data(); 
  IpProfile::begin_run();  IpProfile::dump();
  BeamEnergy::begin_run(); //__ Initialization of beamenergy
  int nword = BsCouTab ( BELLE_RUNHEAD );
  if ( BsCouTab ( BELLE_RUNHEAD ) ) {
    Belle_runhead_Manager& runhead_mgr = Belle_runhead_Manager::get_manager();
    Belle_runhead& runhead(runhead_mgr[0]);
    iexpmc = runhead.ExpMC();
  }
}
void user_ana_check::compare(double de, int k, double &temp, int &tk){
  if (de < temp){
     temp = de;
     tk = k;
  }
}

#if defined(BELLE_NAMESPACE)
}
#endif

