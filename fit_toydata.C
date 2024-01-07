#include "../Belle2Style.h"                                                                                                                                                   
#include "../Belle2Utils.h"                                                                                                                                               
#include "../Belle2Labels.h"

void SetBelle2Style() {                                                                                                                                                                                   
  static TStyle* belle2Style = 0;                                                                                                                                                                         
  cout << "\nApplying BELLE2 style settings...\n" << endl ;                                                                                                                                               
  if (!belle2Style)                                                                                                                                                                                       
    belle2Style = Belle2Style();                                                                                                                                                                          
  gROOT->SetStyle("BELLE2");                                                                                                                                                                              
  gROOT->ForceStyle();                                                                                                                                                                               
}

void fit_toydata()
{
  SetBelle2Style();
  gROOT->Reset();
  
  gStyle->SetOptFit(kTRUE);
  Color_t white=10;
  gStyle->SetCanvasColor(white);gSystem->Load("libRooFit");
  using namespace RooFit;
  using namespace RooStats;

  std::ostringstream name;
  std::ostringstream mode;
  
  // #ifdef __CINT__
  //gROOT->ProcessLineSync(".x nLGn.cxx+");
  //#endif

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
  RooRealVar a("a","a", 0.319,-10,10);
  RooRealVar b("b","", -7.48,-100,100);
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

  RooRealVar * slope2= new RooRealVar("slope2", "Slope of Polynomial", -0.0931,-1, 1);
  RooChebychev *chebpol1 = new RooChebychev("chebpol1","Chebshev Polynomial ", *md0, RooArgList(*slope2));
  
  RooRealVar *sig = new RooRealVar("sig","#signal events",-100,4000);
  RooRealVar *BKG = new RooRealVar("BKG","#signal events",-30,40000);
  
  RooProdPdf *pdf_signal_signal = new RooProdPdf("pdf_signal_signal","pdf_signal_signal",RooArgList(*dG_dm4,*dG_md2));
  RooProdPdf *pdf_signal = new RooProdPdf("pdf_signal","pdf_signal",RooArgList(*t_bkg,*chebpol1));
  RooAddPdf *total = new RooAddPdf("total","total",RooArgList(*pdf_signal_signal,*pdf_signal),RooArgList(*sig,*BKG));
  
  float tempe,temps,tempds,tempde,temp2,tempdsb,tempeb,tempsb,tempdas,tempdae,tempdab,tempdeb ;
  int Enty = 0;
  int Entx = 0;
  
  TCanvas c1("c1", "First canvas", 800, 800);
  c1.Divide(2,2);
  TCanvas c2("c2", "Second canvas", 800, 800);
  c2.Divide(2,2);
  
  TPostScript ps ("fit-toy1350.ps",112);
  FILE *f1, *f2;
  f1 = fopen("signal-1350-toy_2.txt","w");
  TPostScript ps1 ("fit-toy13501.ps",112);
  f2 = fopen("signal-13501-toy_2.txt","w");
  
  for ( int i=1; i<2001; i++){
    
    string infiles="./incbkg";
    std::ostringstream temp1;
    temp1 << infiles << "/bkginc-" << i;
    infiles = temp1.str();
    infiles.append(".txt");
    std::cout  << infiles << endl;
    RooDataSet *data = RooDataSet::read(infiles.c_str(), RooArgList(*deltam,*md0));
    
    string infile1="./incdat";
    std::ostringstream temp2;
    temp2 << infile1 << "/siginc-" << i;
    infile1 = temp2.str();
    infile1.append(".txt");
    std::cout  << infile1 << endl;
    RooDataSet *data1 = RooDataSet::read(infile1.c_str(), RooArgList(*deltam,*md0));
    data->append(*data1);
    
    RooFitResult *fitRes = total->fitTo(*data,Extended(true),Minos(true),Save());
    
    RooPlot* xframe = deltam->frame(Bins(40));
    data->plotOn(xframe);
    total->paramOn(xframe);
    total->plotOn(xframe);       
    total->plotOn(xframe,Components(RooArgSet(*total)),LineColor(kGreen),LineStyle(kDashed));

    name.str("");
    name << "#DeltaM [" << i <<"]";
    xframe->SetTitle(name.str().c_str());                                                                                               
    if (Entx==0) {ps.NewPage();}
    if (Entx==3 ) {c2.cd(4); Entx=0;
      xframe->Draw();
      c2.Update();}
    else {c2.cd(Entx+1); Entx++;
      xframe->Draw();}
    tempds = sig->getVal();
    tempde = sig->getError();
    tempdas = sig->getAsymErrorLo();
    tempdae = sig->getAsymErrorHi();
    //fprintf(f2,"%f \t %f \n",tempds,tempde);
    fprintf(f2,"%f \t %f \t %f \t %f \n",tempds,tempde,tempdas,tempdae);                                                                                                                  
    
    RooPlot* yframe = md0->frame(Bins(40));
    data->plotOn(yframe);
    total->paramOn(yframe);
    total->plotOn(yframe);
    total->plotOn(yframe,Components(RooArgSet(*total)),LineColor(kGreen),LineStyle(kDashed));
    
    name.str("");
    name << "M_{D^{0}} [" << i <<"]"; 
    yframe->SetTitle(name.str().c_str());
    //c3->SaveAs(finaF.c_str());
    if (Enty==0) {ps.NewPage();}
    if (Enty==3 ) {c1.cd(4); Enty=0;
      yframe->Draw();
      c1.Update();}
    else {c1.cd(Enty+1); Enty++;
      yframe->Draw();}
    
    tempsb = BKG->getVal();
    tempeb = BKG->getError();
    tempdab = BKG->getAsymErrorLo();
    tempdeb = BKG->getAsymErrorHi();
    fprintf(f1,"%f \t  %f \t  %f \t  %f \n",tempsb,tempeb,tempdab,tempdeb);
    //fprintf(f1,"%f \t  %f \n",tempsb,tempeb);
    //fprintf(f2,"%f \t %f \n",temp3,temp4);   
    
  }

  fclose (f1);
  fclose (f2);
  ps.Close();
  ps1.Close();
 
  return 0;
}


