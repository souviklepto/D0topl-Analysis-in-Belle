#if !defined(USERINFO_H_INCLUDED)
#define USERINFO_H_INCLUDED 
  
#include "belle.h"
#include "particle/ParticleUserInfo.h"
///down change 
 
#if defined(BELLE_NAMESPACE)
namespace Belle{ 
#endif  
 
class UserInfo : public ParticleUserInfo 
{  
public: 
 
// Default constructor
UserInfo();
 
/// Copy constructor
UserInfo(const UserInfo &);
 
/// Destructor
 virtual ~UserInfo();
 
/// constructs self object.
 UserInfo * clone(void) const;
 
/// Copy operator
  UserInfo & operator = (const UserInfo &);
 
public : 
void user(const double &v) { m_user = v; }

 void e9ovr ( const  double &v) { m_e9ovr      = v; } 
 void pidpi ( const  double &v) { m_pidpi      = v; } 
 void etag ( const  double &v) { m_etag      = v; }
 void deltz ( const  double &v) { m_deltz      = v; } 
 void rphiv ( const  double &v) { m_rphiv      = v; } 
 void phi0 ( const  double &v) { m_phi0      = v; } 
 void kappa ( const  double &v) { m_kappa      = v; } 
 void tanld ( const  double &v) { m_tanld      = v; } 
 void pivotx ( const  double &v) { m_pivotx      = v; } 
 void pivoty ( const  double &v) { m_pivoty      = v; } 
 void pivotz ( const  double &v) { m_pivotz      = v; } 
 void nhit ( const  double &v) { m_nhit      = v; } 
 void mass ( const  double &v) { m_mass      = v; } 
 void theta1 ( const  double &v) { m_theta1      = v; } 
 void theta2 ( const  double &v) { m_theta2      = v; } 
 void eg1 ( const  double &v) { m_eg1      = v; } 
 void eg2 ( const  double &v) { m_eg2      = v; } 
 void momp ( const  double &v) { m_momp      = v; } 
 void theta ( const  double &v) { m_theta      = v; } 
 void mom ( const  double &v) { m_mom      = v; } 
 //Added for Ks
 void chisq ( const  double &v) { m_chisq      = v; }
 void fl ( const  double &v) { m_fl      = v; }
 void zdist ( const  double &v) { m_zdist      = v; }
 void dphi ( const  double &v) { m_dphi      = v; }
 void nb_vlike(const double &v) {m_nb_vlike   =v; }// for nis vlike me
 void nb_nolam(const double &v) {m_nb_nolam   =v; }
// 8th july 2014---------------------------------------------------------------me
 void pmag(const double &v) {m_pmag   =v; }
 void drp(const double &v) {m_drp   =v; }
 void drn(const double &v) {m_drn   =v; }
 void decang(const double &v) {m_decang   =v; }
 void svdhit1(const double &v) {m_svdhit1   =v; }
 void svdhit2(const double &v) {m_svdhit2   =v; }
 void cdc_r1(const double &v) {m_cdc_r1   =v; }
 void cdc_z1(const double &v) {m_cdc_z1   =v; }
 void cdc_r2(const double &v) {m_cdc_r2   =v; }
 void cdc_z2(const double &v) {m_cdc_z2   =v; }
 void atc_1(const double &v) {m_atc_1   =v; }
 void atc_2(const double &v) {m_atc_2   =v; }
 void m_lambda(const double &v) {m_m_lambda   =v; }
 void p_lab_p(const double &v) {m_p_lab_p   =v; }
 void p_lab_n(const double &v) {m_p_lab_n   =v; }
 void sin_th_p(const double &v) {m_sin_th_p   =v; }
 void sin_th_n(const double &v) {m_sin_th_n   =v; }
 void gid ( const  double &v) { m_gid      = v; }
 void tinfo1 ( const  double &v) { m_tinfo1      = v; } 
 void tinfo2 ( const  double &v) { m_tinfo2      = v; }
 void isthep ( const  double &v) { m_isthep      = v; }
 void gma1 ( const  double &v) { m_gma1      = v; }
 void gma2 ( const  double &v) { m_gma2      = v; } 
 void gma3 ( const  double &v) { m_gma3      = v; }
 void gma4 ( const  double &v) { m_gma4      = v; }





 void cl(const double &v) {m_cl =v;} //confidence level//me
 void ndf(const unsigned &v) {m_ndf =v;} //degree of freedom//me

 //  void (const double &v) {m_  =v;} // 
 
 
 void test(const double &);
 const double & user(void)  const { return m_user; } 
 const double & etag(void)     const {return m_etag; }
 const double & e9ovr(void)     const {return m_e9ovr; }
 const double & pidpi(void)     const {return m_pidpi; } 
 const double & deltz(void)     const {return m_deltz; } 
 const double & rphiv(void)     const {return m_rphiv; } 
 const double & phi0(void)     const {return m_phi0; } 
 const double & kappa(void)     const {return m_kappa; } 
 const double & tanld(void)     const {return m_tanld; } 
 const double & pivotx(void)     const {return m_pivotx; } 
 const double & pivoty(void)     const {return m_pivoty; } 
 const double & pivotz(void)     const {return m_pivotz; } 
 const double & nhit(void)     const {return m_nhit; } 
 const double & mass(void)     const {return m_mass; } 
 const double & theta1(void)     const {return m_theta1; } 
 const double & theta2(void)     const {return m_theta2; } 
 const double & eg1(void)     const {return m_eg1; } 
 const double & eg2(void)     const {return m_eg2; } 
 const double & momp(void)     const {return m_momp; } 
 const double & theta(void)     const {return m_theta; } 
 const double & mom(void)     const {return m_mom; } 
 //Added for Ks
 const double & chisq(void)     const {return m_chisq; }
 const double & fl(void)     const {return m_fl; }
 const double & zdist(void)     const {return m_zdist; }
 const double & dphi(void)     const {return m_dphi; }
 const double & nb_vlike(void)  const{return m_nb_vlike; }//for nis no lam me
 const double & nb_nolam(void)  const{return m_nb_nolam; }//for nis no lam me
//---------------8th july--
 const double & pmag(void)     const {return m_pmag; }
 const double & drp(void)     const {return m_drp; }
 const double & drn(void)     const {return m_drn; }
 const double & decang(void)     const {return m_decang; }
 const double & svdhit1(void)     const {return m_svdhit1; }
 const double & svdhit2(void)     const {return m_svdhit2; } 
 const double & cdc_r1(void)     const {return m_cdc_r1; }
 const double & cdc_z1(void)     const {return m_cdc_z1; }
 const double & cdc_r2(void)     const {return m_cdc_r2; }
 const double & cdc_z2(void)     const {return m_cdc_z2; }
 const double & atc_1(void)     const {return m_atc_1; }
 const double & atc_2(void)     const {return m_atc_2; }
 const double & m_lambda(void)     const {return m_m_lambda; }
 const double & p_lab_p(void)     const {return m_p_lab_p; }
 const double & p_lab_n(void)     const {return m_p_lab_n; }
 const double & sin_th_p(void)     const {return m_sin_th_p; }
 const double & sin_th_n(void)     const {return m_sin_th_n; }
 const double & gid(void)     const {return m_gid; } 
 const double & cl(void)       const { return m_cl; }//me 23rd may
 const unsigned & ndf(void)       const { return m_ndf; }//me
 const double & tinfo1(void)     const {return m_tinfo1; } 
 const double & tinfo2(void)     const {return m_tinfo2; }
 const double & isthep(void)     const {return m_isthep; }
 const double & gma1(void)     const {return m_gma1; }
 const double & gma2(void)     const {return m_gma2; } 
 const double & gma3(void)     const {return m_gma3; }
 const double & gma4(void)     const {return m_gma4; }
 

// const double & (void)          const { return m_; }
 const double & test(void)             const { return *m_test; }
private:
double m_user;
 double  m_etag;
 double  m_e9ovr;
 double  m_pidpi; 
 double  m_deltz; 
 double  m_rphiv; 
 double  m_phi0; 
 double  m_kappa; 
 double  m_tanld; 
 double  m_pivotx; 
 double  m_pivoty; 
 double  m_pivotz; 
 double  m_nhit; 
 double  m_mass; 
 double  m_theta1; 
 double  m_theta2; 
 double  m_eg1; 
 double  m_eg2; 
 double  m_momp; 
 double  m_theta; 
 double  m_mom; 
 //Added for Ks
 double  m_chisq;
 double  m_fl;
 double  m_zdist;
 double  m_dphi;
 double m_nb_vlike;//nis me
 double m_nb_nolam; // nis me
 double  m_gid;
  double m_pmag;
  double m_drp;
  double m_drn;
  double m_decang;
  double m_svdhit1;
  double m_svdhit2;
  double m_cdc_r1;
  double m_cdc_z1;
  double m_cdc_r2;
  double m_cdc_z2;
  double m_atc_1;
  double m_atc_2;
  double m_m_lambda;
  double m_p_lab_p;
  double m_p_lab_n;
  double m_sin_th_p;
  double m_sin_th_n;
  double  m_tinfo1; 
  double  m_tinfo2;
  double  m_isthep;
  double  m_gma1;
  double  m_gma2; 
  double  m_gma3; 
  double  m_gma4;




  unsigned m_ndf;//me
  double m_cl;//me

 
 double *m_test;
};
#if defined(BELLE_NAMESPACE)
}//namespce Bell
#endif
#endif
