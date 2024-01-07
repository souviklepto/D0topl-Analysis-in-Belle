void toygenrte(){
  using namespace RooFit;
  using namespace RooStats;
  
  gROOT->Reset();
  
  gStyle->SetOptFit(kTRUE);
  Color_t white=10;
  gStyle->SetCanvasColor(white);gSystem->Load("libRooFit");
  using namespace RooFit;
  
  
  RooRealVar *deltam = new RooRealVar("deltam","", 0.1395, 0.158);
  RooRealVar *md0 = new RooRealVar("md0","", 1.8, 1.9);
  
  // DeltaM fit                                                                                                                               
  RooRealVar *mean1_dm = new RooRealVar("mean1_dm", "MEAN of 1st gaussian",0.1453988);
  RooRealVar *sigma1_dm = new RooRealVar("sigma1_dm", "Sigma of 1st gaussian",0.0002429);
  RooGaussian *gauss_dm= new RooGaussian("gauss_dm", "1st gaussian PDF", *deltam, *mean1_dm, *sigma1_dm);
  RooRealVar *delm1_dm=new RooRealVar("delm1_dm", "Difference of Mean",0.000195);
  RooFormulaVar *mean2_dm= new RooFormulaVar("mean2_dm","Mean of 2nd gaussian", "@0+@1",RooArgList(*mean1_dm,*delm1_dm));
  RooRealVar *mean3_dm = new RooRealVar("mean3_dm", "MEAN of 1st gaussian",0.1454552);
  RooRealVar *sigma2_dm = new RooRealVar("sigma2_dm", "Sigma of 2nd gaussian",0.001053);
  RooRealVar *sigma3_dm = new RooRealVar("sigma3_dm", "Sigma of 2nd gaussian",0.0004650);
  RooRealVar *sigma4_dm = new RooRealVar("sigma4_dm", "Sigma of 2nd gaussian",0.00308);
  RooRealVar *mean4_dm = new RooRealVar("mean4_dm", "MEAN of 1st gaussian",0.146445);
  RooGaussian *gauss2_dm= new RooGaussian("gauss2_dm", "1st gaussian PDF", *deltam, *mean2_dm, *sigma2_dm);
  RooGaussian *gauss3_dm= new RooGaussian("gauss3_dm", "3rd gaussian PDF", *deltam, *mean3_dm, *sigma3_dm);
  RooGaussian *gauss4_dm= new RooGaussian("gauss4_dm", "4th gaussian PDF", *deltam, *mean4_dm, *sigma4_dm);
  RooRealVar *a2a_dm1=new RooRealVar("a2a_dm1", " Area2/Area1",0.727);
  RooRealVar *a2a_dm=new RooRealVar("a2a_dm", " Area2/Area1",0.415);
  RooRealVar *a2a_dm2=new RooRealVar("a2a_dm2", " Area2/Area1",0.9642);
  RooAddPdf *dG_dm3= new RooAddPdf("dG_dm3", " 1st GAuss + 2nd Gauss", RooArgList(*gauss_dm,*gauss2_dm), RooArgList(*a2a_dm1));
  RooAddPdf *dG_dm2= new RooAddPdf("dG_dm2", " 1st GAuss + 2nd Gauss", RooArgList(*dG_dm3,*gauss3_dm), RooArgList(*a2a_dm));
  RooAddPdf *dG_dm4= new RooAddPdf("dG_dm4", " 1st GAuss + 2nd Gauss", RooArgList(*dG_dm2,*gauss4_dm), RooArgList(*a2a_dm2));
  
  RooRealVar th("th","a",0.1395);
  RooRealVar a("a","a", 0.319);
  RooRealVar b("b","", -7.48);
  RooRealVar c("c","",0); 
  RooRealVar d("d","",0);  
  RooGenericPdf *t_bkg = new RooGenericPdf("t_bkg",
                                           "(deltam>=th)*pow((deltam-th),a)*exp(-b*(deltam-th))",  
                                           RooArgSet(*deltam,th,a,b));
  
  // M_D0 fit                                                                                                                                                                     
  RooRealVar *mean1_md = new RooRealVar("mean1_md", "MEAN of 1st gaussian",1.864580);
  RooRealVar *sigma1_md = new RooRealVar("sigma1_md", "Sigma of 1st gaussian",0.003151);
  RooGaussian *gauss_md= new RooGaussian("gauss_md", "1st gaussian PDF", *md0, *mean1_md, *sigma1_md);
  RooRealVar *mean_md = new RooRealVar("mean_md", "MEAN of 1st gaussian",1.8,1.9);
  RooRealVar *sigma_md = new RooRealVar("sigma_md", "Sigma of 1st gaussian",0.00782);
  RooGaussian *gauss_mds= new RooGaussian("gauss_mds", "1st gaussian PDF", *md0, *mean1_md, *sigma_md);
  RooRealVar *sls1_md=new RooRealVar("sls1_md","Ratio of Sigma",7.10);
  RooRealVar *srs1_md=new RooRealVar("srs1_md","Ratio of Sigma",0.780);
  RooRealVar *a2a_md=new RooRealVar("a2a_md", " Area2/Area1",0.1960);
  RooRealVar *delm1_md=new RooRealVar("delm1_md", "Difference of Mean",0.00428);
  RooFormulaVar *mean2_md= new RooFormulaVar("mean2_md","Mean of 2nd gaussian", "@0+@1",RooArgList(*mean1_md,*delm1_md));
  RooFormulaVar *sigmaL_md= new RooFormulaVar("sigmaL_md", "Sigma of 2nd gaussian", " @0*@1", RooArgList(*sigma1_md,*sls1_md));
  RooFormulaVar *sigmaR_md= new RooFormulaVar("sigmaR_md", "Sigma of 2nd gaussian", " @0*@1", RooArgList(*sigma1_md,*srs1_md));
  RooBifurGauss *gauss3_md= new RooBifurGauss("gauss3_md", "3nd gaussian PDF", *md0, *mean2_md, *sigmaL_md,*sigmaR_md);
  RooAddPdf *dG_md1= new RooAddPdf("dG_md1", " 1st GAuss + 3rd Gauss", RooArgList(*gauss3_md,*gauss_md), RooArgList(*a2a_md));
  RooRealVar *a2a_md1=new RooRealVar("a2a_md1", " Area2/Area1",0.8144);
  RooAddPdf *dG_md2= new RooAddPdf("dG_md2", " 2nd Gauss + Argus", RooArgList(*dG_md1,*gauss_mds), RooArgList(*a2a_md1));
  
  RooRealVar * slope2= new RooRealVar("slope2", "Slope of Polynomial", -0.0931);
  RooChebychev *chebpol1 = new RooChebychev("chebpol1","Chebshev Polynomial ", *md0, RooArgList(*slope2));
  
  RooProdPdf *pdf_signal_signal = new RooProdPdf("pdf_signal_signal","pdf_signal_signal",RooArgList(*dG_dm4,*dG_md2));
  RooProdPdf *pdf_signal = new RooProdPdf("pdf_signal","pdf_signal",RooArgList(*t_bkg,*chebpol1));
  
  RooRandom::randomGenerator()->SetSeed(100);
  //r ( int geni = 6; geni < 10; geni++){
  double SIG=15.0;
  int sig_2 =59;
  int BKG = 1424;
  
  for ( int i=1; i < 2001; i++) {
    //sig = 8.0
    
    Int_t gen_bkg_fff = 0;
    Int_t gen_sig_fff = 0; 
    
    TRandom3 rnd_1;
    rnd_1.SetSeed(0);
    gen_bkg_fff = rnd_1.Poisson(BKG);
    
    if (SIG!= 0) {gen_sig_fff = rnd_1.Poisson(SIG);}
    else {gen_sig_fff = 0;}
    
    
  string infile1="./incdat";
  std::ostringstream temp1;
  temp1 << infile1 << "/siginc-" << i;
  infile1 = temp1.str();
  infile1.append(".txt");
  std::cout  << infile1 << endl;
  RooDataSet *datainc1 = pdf_signal_signal->generate(RooArgSet(*deltam,*md0),gen_sig_fff);
  //RooDataSet *datafff = model_klm_fff_1->generate(RooArgSet(*x),gen_bkg_fff);
  if(datainc1 != NULL) datainc1->write(infile1.c_str());
  
  
  string infile11="./incbkg";
  std::ostringstream temp11;
  temp11 << infile11 << "/bkginc-" << i  ;
  infile11 = temp11.str();
  infile11.append(".txt");
  std::cout  << infile11 << endl;
  //RooDataSet *data_sigfff = smodel_fff->generate(RooArgSet(*x),gen_sig_fff);
  RooDataSet *datainc2 = pdf_signal->generate(RooArgSet(*deltam,*md0),gen_bkg_fff);
  if(datainc2 != NULL) datainc2->write(infile11.c_str()); 
  
  } 
}
