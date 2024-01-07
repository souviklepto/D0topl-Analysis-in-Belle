#include "../../Belle2Style.h"
#include "../../Belle2Utils.h"
#include "../../Belle2Labels.h"

Double_t fitf(double *x, double *par)
{
  double val= par[0]*(1-pow(par[1],(-par[2] *x[0]))) + par[3];
  return val;
}

void SetBelle2Style() {
  static TStyle* belle2Style = 0;
  cout << "\nApplying BELLE2 style settings...\n" << endl ;
  if (!belle2Style)
    belle2Style = Belle2Style();
  gROOT->SetStyle("BELLE2");
  gROOT->ForceStyle();
}

void ul()  
{
  SetBelle2Style();
  
  int num = 2000;
  float p[3000][3],pull[3000];
  FILE *fp;
  int count;
  count=0;
  int q[5];
 
  for (int k=1; k<6; k++){

    //string infiles="/home/souvik/leptoquark_analysis/macros/toy_study/upper_limit/bc0";
    string infiles="/home/souvik/leptoquark_analysis/macros/toy_study/bc0/txt_files";
    std::ostringstream temp1;
    temp1 << infiles << "/out" << k;
    infiles = temp1.str();
    infiles.append(".txt");
    std::cout  << infiles << endl;

    //std::cout<<"Hi "<<std::endl;
    
    fp = fopen(infiles.c_str(),"r");
    for (int i=0; i<num;i++){
      for (int j=0;j<4;j++){
	fscanf(fp,"%f",&p[i][j]);}                                                                 
      //std::cout<<p[i][0]<< "    " << p[i][1] << std::endl;
      if(p[i][0]>5.33){
	count++ ;
      }                                                              
    //systematics calculations will be done here.
      //double max = p[i][0] + p[i][0]*5/100; //plus systematics
      //double min = p[i][0] - p[i][0]*2/100; //minus systematics
      //TRandom rnd ;
      
    }

    std::cout<<"Hi "<<std::endl;
    
    double p1, p2;
    p1=0;
    p1 = (count/2000.0)*100 ;
    q[k-1]=p1 ;
    
    std::cout<<"q "<<q[k]<<std::endl;                                            
    std::cout<<"p1 "<<p1<<std::endl;
    count=0;
  }

  int x[5]= {10,15,20,25,30};

  auto c1 = new TCanvas("c4","c4",200,0,600,400);
  TGraph *gr=new TGraph(5,x,q);
  gr->GetYaxis()->SetRangeUser(0,110);
  gr->GetXaxis()->SetRangeUser(-4,40);
  gr->GetXaxis()->SetTitle("Yield");
  gr->GetYaxis()->SetTitle("C.L.(%)");
  gr->GetYaxis()->CenterTitle();
  gr->GetXaxis()->CenterTitle();
  gr->Draw("AP");


  TF1 *f= new TF1("fitf",fitf,10,30,4);
  f->SetParNames("a","b","n","c");
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(kBlue);
  gr->SetFillColor(0);
  gr->SetLineColor(0);
  f->SetLineWidth(4);
  gr->SetTitle("C.L. without Systematics");
  f->SetLineColor(kBlue);
  f->SetParameters(0,1.2,1.8,-0.8);
  gr->Fit("fitf","R");
  gPad->Modified();
  c1->BuildLegend();
  //c1->Update();
  c1->GetFrame()->SetBorderSize(0);
  c1->Update();

    
  TLine *line = new TLine(10,90,13.0,90);
  line->SetLineColor(kCyan);
  line->SetLineWidth(2);
  line->SetLineStyle(kSolid);
  line->Draw();

  TLine *line1 = new TLine(13.0,8,13.0,90);
  line1->SetLineColor(kCyan);
  line1->SetLineWidth(2);
  line1->SetLineStyle(kSolid);
  line1->Draw();
  
  
}
