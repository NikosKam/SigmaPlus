//Author: Dr. P. Ganoti
//Edited: nik
//Analyse ALICE data, Pb-Pb collisions at sqrt(s_NN)=5.02 TeV


#include "TFitResultPtr.h"
//functions
Double_t Background(Double_t *x, Double_t *par);
Double_t Background4(Double_t *x, Double_t *par);
Double_t Background2(Double_t *x, Double_t *par);
Double_t BreitWignerCombFit(Double_t *x, Double_t *par);
Double_t BreitWigner(Double_t *x, Double_t *par);
Double_t gaussianCombined(Double_t *x, Double_t *par);
Double_t gaussianCombined4(Double_t *x, Double_t *par);
Double_t gaussian(Double_t *x, Double_t *par);

//main func 
void readFilePlus() {

  gStyle->SetOptTitle(0);

  Double_t xa[12]={1.,2.,2.5,3.,3.5,4.,4.5,5.5,7.0,9., 12., 15.};

  //opening the data file
  TFile *f1 = new TFile("outPbPbPionCut120.root");
  f1->ls();
    
  TList *l1 = (TList*) f1->Get("RsnOut_SigmaPlusMinus;1");
  // for MB
  //TList *l1 = (TList*) f1->Get("RsnOut_SigmaPMkINT7;1");

  // 
  TH1F *hAEventsVsMulti= (TH1F *) l1->FindObject("hAEventsVsMulti");// TH1F : 1-D histogram with 1 float per channel

  TH1F *hTotalEvents=new TH1F (*hAEventsVsMulti);

  //producing the 1st canvas
  TCanvas *c1 = new TCanvas("c1","",800,800);//TCanvas : an area mapped to a window directly under the control of the display manage
  hTotalEvents->Draw();// drawing the 1st histogram
    
  TH1D *mass=new TH1D("mass"," ",11 , xa);// TH1D : 1-D histogram with 1 double per channel
  TH1D *width=new TH1D("width"," ",11 , xa);
  TH1D *rawYield=new TH1D("rawYield"," ",11 , xa);
  TH1D *hchi2=new TH1D("hchi2"," ",11, xa);
  TH1D *hsignif=new TH1D("hsignif"," ",11, xa);

  TH2F *hVzVsCent1 = (TH2F*) l1->FindObject("hVzVsCent");// TH2F : 2-D histogram with 1 float per channel
    
  TH1D *vZ0_10_1= (TH1D *) hVzVsCent1->ProjectionY(Form("cent : %3.1f-%3.1f GeV/c",hVzVsCent1->GetXaxis()->GetBinCenter(hVzVsCent1->GetXaxis()->FindBin(0+0.001)),
              hVzVsCent1->GetXaxis()->GetBinCenter(hVzVsCent1->GetXaxis()->FindBin(10-0.001))), hVzVsCent1->GetXaxis()->FindBin(0+0.001), hVzVsCent1->GetXaxis()->FindBin(10-0.001), "e");
    
  TH1D *vZ0_10= new TH1D(*vZ0_10_1);

  TFitResultPtr fFitResultPtr=0;
    
  Double_t integral, integralFromHisto, integralBinCount;
  Double_t integralAll, significance;
  Double_t errorY;

  // producing 3 more canvases 
  TCanvas *c3=new TCanvas("c3","",800,800);
  c3->Divide(3,4);// will contain multiple graphs
  TCanvas *c4=new TCanvas("c4","",800,800);
  c4->Divide(3,4);
  TCanvas *c5=new TCanvas("c5","",800,500);
  TCanvas *c6=new TCanvas("c6","",800,500);
    
  Double_t binNorm=1.;
  Double_t lowEdge=0.;
  
  //defining fitting function
  TF1 *fau1=new TF1("fau1", BreitWignerCombFit, 1.2, 2., 7);// TF1: 1-D function

  for (Int_t j = 0; j<9; j++) {
    Double_t integralBinCounting=0.;
    binNorm=1.;
  

    THnSparse *h = 0;// THnSparse : histogram with only a small fraction of bins filled
    THnSparse *ha = 0;

    THnSparse *h2 = 0;
    THnSparse *ha2 = 0;
    
    THnSparse *h4 = 0;
    THnSparse *ha4 = 0;
    
    h = (THnSparse*) l1->FindObject("RsnTaskSigPM_SigPM_SigmaP");// for Sigma
    ha = (THnSparse*) l1->FindObject("RsnTaskSigPM_SigPM_ASigmaP");// for anti-Sigma
    
    TH3D *h3 = (TH3D*) h->Projection(0, 1, 2, "E");
    TH3D *h3a = (TH3D*) ha->Projection(0, 1, 2, "E");
    //TH1D *h1d = (TH1D*) h3->ProjectionX("", 20, 24, 5, 10, "e");    
    TH1D *h1d = (TH1D*) h3->ProjectionX(Form("p_{T} : %3.1f-%3.1f GeV/c",xa[j], xa[j+1]),// ProjectionX : project a 2-D histogram into a 1-D histogram along X-axis
                                                        h3->GetYaxis()->FindBin(xa[j]+0.001),
                                                        h3->GetYaxis()->FindBin(xa[j+1]-0.001), h3->GetZaxis()->FindBin(0+0.001), h3->GetZaxis()->FindBin(10-0.001), "e");

    TH1D *h1da = (TH1D*) h3a->ProjectionX(Form("ap_{T} : %3.1f-%3.1f GeV/c",xa[j], xa[j+1]),
                                                        h3->GetYaxis()->FindBin(xa[j]+0.001),
                                                        h3->GetYaxis()->FindBin(xa[j+1]-0.001), h3->GetZaxis()->FindBin(0+0.001), h3->GetZaxis()->FindBin(10-0.001), "e");
    
    // calculating the pair transverse energy
    TH1D *h1 = new TH1D(Form("pair p_{T} : %3.1f-%3.1f GeV/c",xa[j], xa[j+1]), "h1", 400, 1.2, 2.);
    TH1D *h1a = new TH1D(Form("apair p_{T} : %3.1f-%3.1f GeV/c",xa[j], xa[j+1]), "h1a", 400, 1.2, 2.); 
    for(Int_t i=1; i<(h1d->GetNbinsX() + 1); i++) {  
      h1->SetBinContent(i, h1d->GetBinContent(i));
      h1->SetBinError(i, h1d->GetBinError(i));
      
      h1a->SetBinContent(i, h1da->GetBinContent(i));
      h1a->SetBinError(i, h1da->GetBinError(i));  
    }
    
    delete h1d;
    delete h3;
    delete h1da;
    delete h3a;
    h=0;
    ha=0;
    
    //h1=0;
    //ha1=0;
    
    h1->Add(h1a);
    
    h = (THnSparse*) l1->FindObject("RsnTaskSigPM_SigPM_SigmaPmix");
    ha = (THnSparse*) l1->FindObject("RsnTaskSigPM_SigPM_ASigmaPmix");
    
    h3 = (TH3D*) h->Projection(0, 1, 2, "E");
    h3a = (TH3D*) ha->Projection(0, 1, 2, "E");
    
    h1d = (TH1D*) h3->ProjectionX(Form("p_{T} : %3.1f-%3.1f GeV/c",xa[j], xa[j+1]),
                                                h3->GetYaxis()->FindBin(xa[j]+0.001),
                                                h3->GetYaxis()->FindBin(xa[j+1]-0.001), h3->GetZaxis()->FindBin(0+0.001), h3->GetZaxis()->FindBin(10-0.001), "e");
    
    h1da = (TH1D*) h3a->ProjectionX(Form("ap_{T} : %3.1f-%3.1f GeV/c",xa[j], xa[j+1]),
                                                  h3->GetYaxis()->FindBin(xa[j]+0.001),
                                                  h3->GetYaxis()->FindBin(xa[j+1]-0.001), h3->GetZaxis()->FindBin(0+0.001), h3->GetZaxis()->FindBin(10-0.001), "e");
    
    TH1D *h1mix = new TH1D(Form("mix pair p_{T} : %3.1f-%3.1f GeV/c",xa[j], xa[j+1]), "h1mix", 400, 1.2, 2.);
    TH1D *h1amix = new TH1D(Form("mix apair p_{T} : %3.1f-%3.1f GeV/c",xa[j], xa[j+1]), "h1amix", 400, 1.2, 2.);
    for(Int_t i=1; i<(h1d->GetNbinsX() + 1); i++) {  
      h1mix->SetBinContent(i, h1d->GetBinContent(i));
      h1mix->SetBinError(i, h1d->GetBinError(i));
      
      h1amix->SetBinContent(i, h1da->GetBinContent(i));
      h1amix->SetBinError(i, h1da->GetBinError(i));    
    }
    
    delete h1d;
    delete h3;
    delete h1da;
    delete h3a;
    
    h1mix->Add(h1amix);
    

    Double_t h1_integral = h1->Integral(300,400);
    Double_t h1mix_integral = h1mix->Integral(300,400);
    Double_t ratio = h1_integral/h1mix_integral;

    h1mix->Scale(ratio);
      
    c3->cd(j+1);
    h1->SetMarkerStyle(20);
    h1->SetMarkerSize(0.3);
    h1->Draw();
    h1mix->SetMarkerStyle(20);
    h1mix->SetMarkerSize(0.3);
    h1mix->SetMarkerColor(2);
    h1mix->Draw("same");
    
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);

    c4->cd(j+1);
    TH1D *hSigma = new TH1D(*h1);
    TH1D *hSandB = new TH1D(*h1);
    hSigma->Add(h1mix,-1);

    /* if (j==1) { */
    /*   hSigma->Rebin(2); */
    /*   binNorm=2.; */
    /* } */

    TF1 *fu=new TF1("fu", BreitWignerCombFit, 1.2, 2., 7);
    
    fu->SetParameter(0,6.15758e-01);
    fu->SetParameter(1, 3.77086e+01);
    fu->SetParameter(2,-2.42891e+01);
    fu->SetParameter(3, 4.85593e+00);
    fu->SetParameter(4,2.01817e+03);
    fu->SetParameter(5,1.3828);
    fu->SetParameter(6,0.039);

    hSigma->SetMinimum(0.);
    fFitResultPtr = (hSigma->Fit("fu","S", "", 1.28, 1.8));
                
    TF1 *signal=new TF1("signal", BreitWigner, 1.28, 2., 3);
    TF1 *backg=new TF1("backg", Background, 1.28, 2., 4);
      
    //if(hSigma->GetFunction("fu")) {
    Double_t a1    = hSigma->GetFunction("fu")->GetParameter(0);
    Double_t a2 =    hSigma->GetFunction("fu")->GetParameter(1);
    Double_t a3 =    hSigma->GetFunction("fu")->GetParameter(2);
    Double_t a4   =  hSigma->GetFunction("fu")->GetParameter(3);
    Double_t N   =   hSigma->GetFunction("fu")->GetParameter(4);
    Double_t mean  = hSigma->GetFunction("fu")->GetParameter(5);
    Double_t gamma = hSigma->GetFunction("fu")->GetParameter(6);
      
    Double_t ea1   = hSigma->GetFunction("fu")->GetParError(0);
    Double_t ea2 =   hSigma->GetFunction("fu")->GetParError(1);
    Double_t ea3 =   hSigma->GetFunction("fu")->GetParError(2);
    Double_t ea4   = hSigma->GetFunction("fu")->GetParError(3);
    Double_t eN   =  hSigma->GetFunction("fu")->GetParError(4);
    Double_t emean = hSigma->GetFunction("fu")->GetParError(5);
    Double_t egamma =hSigma->GetFunction("fu")->GetParError(6);
    //}
      
      if (fu->GetNDF()!=0 )  {
        hchi2->SetBinContent(j+1, fu->GetChisquare()/fu->GetNDF());
      }
        
      signal->SetParameters(N, mean, gamma);
      backg->SetParameters(a1, a2, a3, a4);
      fau1->SetParameters(a1, a2, a3, a4, N, mean, gamma);
          
      signal->SetLineColor(4);
      backg->SetLineColor(30);
      signal->Draw("same");
      backg->Draw("same");
      

      integral = signal->Integral(mean-5.*(gamma), mean+5.*(gamma))/(0.002*binNorm);
      integralAll=hSandB->Integral(hSandB->FindBin(mean-5.*(gamma)), hSandB->FindBin(mean+5.*(gamma)));///(0.001*binNorm)

      mass->SetBinContent(j+1, mean);
      mass->SetBinError(j+1, emean);
      width->SetBinContent(j+1, gamma);
      width->SetBinError(j+1, egamma);
      rawYield->SetBinContent(j+1, integral);
      rawYield->SetBinError(j+1, TMath::Sqrt(integral));
        
      significance=integral/TMath::Sqrt(integralAll);
      hsignif->SetBinContent(j+1, significance);
    
  }

  c4->Print("sigmaPlusInvMasses0_10.pdf");
      
  TFile* outputFile = TFile::Open("sigmaPlusProperties0_10.root", "recreate");
  if (!outputFile || !outputFile->IsOpen()) {
    Printf("Cannot Open OUTPUT file");
    return;
  }
      
  mass->Write();
  width->Write();
  rawYield->Write();
  hchi2->Write();
  vZ0_10->SetName("zvertex0_10");
  vZ0_10->Write();
  hsignif->Write();

}

//----------------------------------------------------
Double_t Background4(Double_t *x, Double_t *par)
{
  return (par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0]);
}

//----------------------------------------------------
Double_t Background2(Double_t *x, Double_t *par)
{
  return (par[0]+par[1]*x[0]+par[2]*x[0]*x[0]);
}

//----------------------------------------------------
Double_t Background(Double_t *x, Double_t *par)
{
  return (TMath::Power((x[0]-(1.115+0.13957)), par[0]) * exp(par[1]+par[2]*x[0]+par[3]*x[0]*x[0]));
}

//----------------------------------------------------
Double_t BreitWignerCombFit(Double_t *x, Double_t *par)
{

  return (TMath::Power((x[0]-(1.115+0.13957)), par[0]) * exp(par[1]+par[2]*x[0]+par[3]*x[0]*x[0]))+par[4]*TMath::BreitWigner(x[0],par[5],par[6]);

}

//----------------------------------------------------
Double_t BreitWigner(Double_t *x, Double_t *par)
{

  return par[0]*TMath::BreitWigner(x[0],par[1],par[2]);

}

//----------------------------------------------------
Double_t gaussian(Double_t *x, Double_t *par)
{
  return par[0] * exp(-(TMath::Power((x[0] - par[1])/par[2], 2)/2.0));
}

//----------------------------------------------------
Double_t gaussianCombined(Double_t *x, Double_t *par)
{
 
  return (par[0] * exp(-(TMath::Power((x[0] - par[1])/par[2], 2)/2.0))) + (par[3]+par[4]*x[0]+par[5]*x[0]*x[0]);
}

//----------------------------------------------------
Double_t gaussianCombined4(Double_t *x, Double_t *par)
{
 
  return (par[0] * exp(-(TMath::Power((x[0] - par[1])/par[2], 2)/2.0))) + (par[3]+par[4]*x[0]+par[5]*x[0]*x[0]+par[6]*x[0]*x[0]*x[0]);
}

