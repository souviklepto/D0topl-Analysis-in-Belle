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

void fit_gauss()
  
{
  SetBelle2Style();
  gROOT->Reset();
  
  gStyle->SetOptFit(kTRUE);
  Color_t white=10;
  gStyle->SetCanvasColor(white);gSystem->Load("libRooFit");
  using namespace RooFit;
 
  RooRealVar *deltam = new RooRealVar("deltam","Pull", -6, 6);                                                                                                                                       
  //RooRealVar *deltam = new RooRealVar("deltam","Yield", -100, 100);
  
  RooDataSet *data = new RooDataSet("data","", RooArgSet(*deltam));
  int num = 2000;
  float p[3000][3],pull[3000];
  FILE *fp;
  
  fp = fopen("signal-13501-toy.txt","r");
  for (int i=0; i<num;i++){
    for (int j=0;j<4;j++){
      fscanf(fp,"%f",&p[i][j]);}
    pull[i]=p[i][0];
    pull[i]=(p[i][0]-5.0)*(1/p[i][1]);
    
    deltam->setVal(pull[i]);
    data->add(RooArgSet(*deltam));
  }
  
  
  RooRealVar *mean = new RooRealVar("mean", "MEAN of 1st gaussian",0.11,-1, 1);
  //RooRealVar *mean = new RooRealVar("mean", "MEAN of 1st gaussian", 5.0,-100, 100);                                                                                      
  
  //RooRealVar *sigma = new RooRealVar("sigma", "Sigma of 1st gaussian",0,100);
  RooRealVar *sigma = new RooRealVar("sigma", "Sigma of 1st gaussian",0,100);
  
  RooGaussian *dG= new RooGaussian("dG", "1st gaussian PDF", *deltam, *mean, *sigma);
  
  TCanvas *c1 =new TCanvas ("c1", "pull distribution plot", 700, 600);
  RooFitResult *fitresult = dG->fitTo(*data, Extended(false),Minos(true));
  RooPlot* deltaEplot = deltam->frame(Bins(100));
  data->plotOn(deltaEplot);
  dG->paramOn(deltaEplot);
  dG->plotOn(deltaEplot);
  deltaEplot->SetTitle("");
  deltaEplot->GetYaxis()->CenterTitle();
  deltaEplot->GetXaxis()->CenterTitle();
  deltaEplot->Draw();
  fclose(fp);
}
