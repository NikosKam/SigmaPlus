#include <TH2.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TFitResultPtr.h"

//main func 
void finale() {

    TLegend *leg = new TLegend(0.7,0.7,0.5,0.5);

    //opening the data file
    TFile *f1 = new TFile("sigmaPlusProperties0_10.root");
    TH1D *rawYield0_10= (TH1D *) f1->Get("rawYield");
    leg->AddEntry(rawYield0_10, "0-10", "p");

    TFile *f2 = new TFile("sigmaPlusProperties30_50.root");
    TH1D *rawYield30_50= (TH1D *) f2->Get("rawYield");
    leg->AddEntry(rawYield30_50, "30-50", "p");
    
    TFile *f3 = new TFile("sigmaPlusProperties50_90.root");
    TH1D *rawYield50_90= (TH1D *) f3->Get("rawYield");
    leg->AddEntry(rawYield50_90, "50-90", "p");
    
    TFile *f4 = new TFile("sigmaPlusProperties0_90.root");
    TH1D *rawYield0_90= (TH1D *) f4->Get("rawYield");
    leg->AddEntry(rawYield0_10, "0-90", "p");


    TCanvas *c1 = new TCanvas("c1","",800,800);//TCanvas : an area mapped to a window directly under the control of the display manage
    gPad->SetLogy();
    rawYield0_10->SetMarkerStyle(20);
    rawYield0_10->SetMarkerColor(1);
    rawYield0_10->SetLineColor(20);
    rawYield30_50->SetMarkerStyle(21);
    rawYield30_50->SetMarkerColor(2);
    rawYield30_50->SetLineColor(21);
    rawYield50_90->SetMarkerStyle(22);
    rawYield50_90->SetMarkerColor(3);
    rawYield50_90->SetLineColor(22);
    rawYield0_90->SetMarkerStyle(23);
    rawYield0_90->SetMarkerColor(4);
    rawYield0_90->SetLineColor(23);
    
    rawYield0_10->SetMaximum(0.1);
    //rawYield0_10->GetXaxis()->SetRangeUser(0.,9.);
    rawYield0_10->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    rawYield0_10->Draw("p");
    rawYield30_50->Draw("psame");
    rawYield50_90->Draw("psame");
    rawYield0_90->Draw("psame");
    leg->Draw();
    c1->Print("finale.pdf");
}