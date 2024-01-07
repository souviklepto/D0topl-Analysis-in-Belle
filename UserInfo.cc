// Based on example from J. Tanaka 
//Prepared by : V. Bhardwaj  
#include"belle.h"
#include "UserInfo.h" 
#if defined(BELLE_NAMESPACE)
namespace Belle
{
#endif
UserInfo::UserInfo():
 m_user( 0.0 ),
 m_e9ovr( 0.0),
 m_gid( 0.0),
 m_etag( 0.0),
 m_tinfo1( 0.0),  
 m_tinfo2( 0.0),
 m_isthep( 0.0),
 m_gma1( 0.0),
 m_gma2( 0.0),  
 m_gma3( 0.0),
 m_gma4( 0.0),

 m_pidpi( 0.0),  
 m_deltz( 0.0),  
 m_rphiv( 0.0),  
 m_phi0( 0.0),  
 m_kappa( 0.0),  
 m_tanld( 0.0),  
 m_pivotx( 0.0),  
 m_pivoty( 0.0),  
 m_pivotz( 0.0),  
 m_nhit( 0.0),  
 m_mass( 0.0),  
 m_theta1( 0.0),  
 m_theta2( 0.0),  
 m_eg1( 0.0),  
 m_eg2( 0.0),  
 m_momp( 0.0),  
 m_theta( 0.0),  
 m_mom( 0.0),  
 m_chisq( 0.0), //Added for Ks
 m_ndf(0.0),
 m_fl( 0.0), //Added for Ks 
 m_zdist( 0.0), //Added for Ks 
 m_dphi( 0.0), //Added for Ks 
 m_nb_vlike(0.0),// added for nis ks me
 m_nb_nolam(0.0),// d0
////---8th july....
 m_pmag(0.0),
 m_drp(0.0),
 m_drn(0.0),
 m_decang(0.0),
 m_svdhit1(0.0),
 m_svdhit2(0.0),
 m_cdc_r1(0.0),
 m_cdc_z1(0.0),
 m_cdc_r2(0.0),
 m_cdc_z2(0.0),
 m_atc_1(0.0),
 m_atc_2(0.0),
 m_m_lambda(0.0),
 m_p_lab_p(0.0),
 m_p_lab_n(0.0),
 m_sin_th_p(0.0),
m_sin_th_n(0.0),

 m_test(0){}
UserInfo::~UserInfo()
{
  if( m_test ) delete m_test;
}
UserInfo::UserInfo(const UserInfo &x)
  :ParticleUserInfo(),
   m_user(        x.m_user   ),
 m_e9ovr( x.m_e9ovr),
 m_gid( x.m_gid),
 m_etag( x.m_etag),
 m_tinfo1( x.m_tinfo1), 
 m_tinfo2( x.m_tinfo2),
 m_isthep( x.m_isthep),
 m_gma1( x.m_gma1),
 m_gma2( x.m_gma2), 
 m_gma3( x.m_gma3),
 m_gma4( x.m_gma4),

 m_pidpi( x.m_pidpi), 
 m_deltz( x.m_deltz), 
 m_rphiv( x.m_rphiv), 
 m_phi0( x.m_phi0), 
 m_kappa( x.m_kappa), 
 m_tanld( x.m_tanld), 
 m_pivotx( x.m_pivotx), 
 m_pivoty( x.m_pivoty), 
 m_pivotz( x.m_pivotz), 
 m_nhit( x.m_nhit), 
 m_mass( x.m_mass), 
 m_theta1( x.m_theta1), 
 m_theta2( x.m_theta2), 
 m_eg1( x.m_eg1), 
 m_eg2( x.m_eg2), 
 m_momp( x.m_momp), 
 m_theta( x.m_theta), 
 m_mom( x.m_mom), 
 m_chisq( x.m_chisq),//Added for Ks
 m_ndf(x.m_ndf),
 m_fl( x.m_fl),//Added for Ks 
 m_zdist( x.m_zdist),//Added for Ks 
 m_dphi( x.m_dphi),//Added for Ks 
 m_nb_vlike(x.m_nb_vlike),//for nis me
 m_nb_nolam(x.m_nb_nolam),//do 
  
//8th july......
 m_pmag(x.m_pmag),
 m_drp(x.m_drp),
 m_drn(x.m_drn),
 m_decang(x.m_decang),
 m_svdhit1(x.m_svdhit1),
 m_svdhit2(x.m_svdhit2),
 m_cdc_r1(x.m_cdc_r1),
 m_cdc_z1(x.m_cdc_z1),
 m_cdc_r2(x.m_cdc_r2),
 m_cdc_z2(x.m_cdc_z2),
 m_atc_1(x.m_atc_1),
 m_atc_2(x.m_atc_2),
 m_m_lambda(x.m_m_lambda),
 m_p_lab_p(x.m_p_lab_p),
 m_p_lab_n(x.m_p_lab_n),
 m_sin_th_p(x.m_sin_th_p),
 m_sin_th_n(x.m_sin_th_n),






  
   m_test(0)
{
 if( x.m_test ) m_test = new double( *( x.m_test ) );
}
UserInfo* UserInfo::clone(void) const
{
  UserInfo *x = new UserInfo( *this );
  return x;
}
UserInfo & UserInfo::operator = (const UserInfo &x)
{
 m_user = x.m_user;
 m_e9ovr=  x.m_e9ovr;
 m_gid=  x.m_gid;
 m_tinfo1=  x.m_tinfo1; 
 m_tinfo2=  x.m_tinfo2;
 m_isthep=  x.m_isthep;
 m_gma1=  x.m_gma1;
 m_gma2=  x.m_gma2; 
 m_gma3=  x.m_gma3;
 m_gma4=  x.m_gma4;
 m_etag=  x.m_etag;

 m_pidpi=  x.m_pidpi; 
 m_deltz=  x.m_deltz; 
 m_rphiv=  x.m_rphiv; 
 m_phi0=  x.m_phi0; 
 m_kappa=  x.m_kappa; 
 m_tanld=  x.m_tanld; 
 m_pivotx=  x.m_pivotx; 
 m_pivoty=  x.m_pivoty; 
 m_pivotz=  x.m_pivotz; 
 m_nhit=  x.m_nhit; 
 m_mass=  x.m_mass; 
 m_theta1=  x.m_theta1; 
 m_theta2=  x.m_theta2; 
 m_eg1=  x.m_eg1; 
 m_eg2=  x.m_eg2; 
 m_momp=  x.m_momp; 
 m_theta=  x.m_theta; 
 m_mom=  x.m_mom; 
 m_chisq=  x.m_chisq;
 m_ndf=  x.m_ndf;


{
 if( x.m_test ) m_test = new double( *( x.m_test ) );
}




  if(   m_test ) delete m_test;
  if( x.m_test ){
    m_test = new double( *( x.m_test ) );
  }else{
    m_test = 0;
  }
  return *this;
}

void 
UserInfo::test(const double &v )
{
  if( m_test ){
    *m_test = v;
  }else{
    m_test = new double( v );
  }
}

#if defined(BELLE_NAMESPACE)
}//namespace Belle
#endif
 
