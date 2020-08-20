#include <Math/SpecFuncMathMore.h>
#include <TMath.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include "TFile.h"
#include <TH1.h>
#include <TPad.h>

const Char_t *hadrons[]={"Pion","Kaon","Proton"};
const Char_t *particle[]={"#pi^{+}","#pi^{-}","K^{+}","K^{-}","p","#bar{p}","#pi^{+}#pi^{-}","K^{+}K^{-}","p #bar{p}"};
const Int_t ngap=5;
const Color_t color[]={kMagenta, kRed, kBlue, kGreen, kBlack, kRed, kMagenta, kRed, kBlue, kGreen, kBlack, kRed};
const Color_t color_2[]={kMagenta, kRed, kBlue, kBlack, kRed, kMagenta, kRed, kBlue, kBlack, kRed};
const Int_t style[]={22, 8, 34, 29, 33, 23 };
const Int_t style_2[]={26, 24, 28, 30, 27, 32, 22, 8, 34, 29, 33, 23};
const Double_t n[ngap]={0.0, 0.05, 0.075, 0.2, 0.5};
//const Int_t col = 12;
//const Int_t col = 13;
//const Int_t col = 14;
//const Int_t col = 16;
const Int_t col = 15;

void makeplotstyle(){
  
  TStyle *mystyle = new TStyle("PlottingInStyle", "Style for Summary Plots");
  mystyle->SetPalette(1);
  mystyle->SetCanvasColor(10);
  mystyle->SetHistFillColor(10);
  mystyle->SetHistFillStyle(0);
  mystyle->SetOptTitle(0);
  mystyle->SetOptStat(0);
  mystyle->SetCanvasBorderMode(0);//removes the yellow frame around the canvas
  mystyle->SetPadLeftMargin(0.16);
  mystyle->SetPadBottomMargin(0.15);
  mystyle->SetPadTickX(1);
  mystyle->SetPadTickY(1);
  mystyle->SetAxisColor(1, "X");
  mystyle->SetAxisColor(1, "Y");
  mystyle->SetLabelColor(1, "X");
  mystyle->SetLabelColor(1, "Y");
  mystyle->SetTickLength(0.03, "X");
  mystyle->SetTickLength(0.03, "Y");
  mystyle->SetTitleXSize(0.05);
  mystyle->SetTitleYSize(0.05);
  mystyle->SetNdivisions(505, "X");
  mystyle->SetNdivisions(508, "Y");
  mystyle->SetTitleXOffset(1.2);
  mystyle->SetTitleYOffset(1.4);
  mystyle->SetLabelOffset(0.02, "X");
  mystyle->SetLabelOffset(0.02, "Y");
  mystyle->SetLabelSize(0.05, "X");
  mystyle->SetLabelSize(0.05, "Y");
  //mystyle->SetGridx();

  TFile f("style.root", "RECREATE");
  f.cd();
  mystyle->Write();
  f.Close();
}

void v_cent(Int_t harmonic,const Char_t *PID,  const Char_t *charge){

  makeplotstyle();

  TFile *f = new TFile("answer_flow.root","READ");
  
  Double_t x_centrality[9]  = {2.5, 7.5, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0,};
  Double_t x_centrality_err[9]  = {0.0};
  
  Double_t y[9]={0.0};
  Double_t y_err[9]={0.0};

  Double_t y_pos[9]={0.0};
  Double_t y_err_pos[9]={0.0};

  const Int_t n1=9;
  Double_t x_coord_1[n1], y_coord_1[n1];
  Double_t x_coord_1err[n1] = {0.0};
  Double_t y_coord_1err[n1] = {0.0};
  
  TFile *fstyle = new TFile("style.root");
  TStyle *tsty = (TStyle *)fstyle->Get("PlottingInStyle");
  tsty->cd();
  gStyle->SetOptStat(0);

  gROOT->ForceStyle();
  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1","Flow ",150,10,800,650);
  TPad *pad_up = new TPad("pad_up", "", 0.0, 0.35, 1.0, 1.0, 0, 0, 0);
  pad_up->SetFrameBorderMode(0);
  pad_up->SetBottomMargin(0);
  TPad *pad_rat = new TPad("pad_rat", "", 0.0, 0.0, 1.0, 0.35, 0, 0, 0);
  pad_rat->SetFrameBorderMode(0);
  pad_rat->SetTopMargin(0.005);
  pad_rat->SetBottomMargin(0.25);
 

  c1->GetFrame()->SetBorderSize(115);
  c1->cd();
  pad_up->Draw();
  pad_rat->Draw();
  gStyle->SetOptStat(0);

  c1->GetFrame()->SetFillColor(21);
  c1->GetFrame()->SetBorderSize(115);

  gStyle->SetOptStat(0);
  gStyle->SetPalette(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetHistFillColor(10);
  gStyle->SetHistFillStyle(0);
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  c1->SetBorderMode(0);

  c1->SetFillColor(0);

  char title1[800];
  sprintf(title1,"");

  float xmin1=0.0;
  float xmax1=100.0;
  float ymin1;
  float ymax1;

  float xmin1_rat=0.0;
  float xmax1_rat=100.0;
  float ymin1_rat=0.7;
  float ymax1_rat=1.3;

  float TLatex_y;

    Int_t nPID;

  if( strncmp(PID, "Proton",6)==0 ){
    if(strncmp(charge, "_pos",3)==0){
      nPID=4;
    }
    else{
      nPID=5;
    }
    if(strncmp(charge, "",1)==0){
      nPID=8;
    }
  }
  if( strncmp(PID, "Kaon",4)==0 ){
    if(strncmp(charge, "_pos",3)==0){
      nPID=2;
    }
    else{
      nPID=3;
    }
    if(strncmp(charge, "",1)==0){
      nPID=7;
    }
  }
  if( strncmp(PID, "Pion",4)==0 ){
    if(strncmp(charge, "_pos",3)==0){
      nPID=0;
    }
    else{
      nPID=1;
    }
    if(strncmp(charge, "",1)==0){
      nPID=6;
    }
  }

  if(harmonic == 3){
    ymin1=0.00;
    ymax1=0.04;
    TLatex_y=0.037;

  }
  if(harmonic == 2){
    ymin1=0.01;
    ymax1=0.11;
    TLatex_y=0.103;
  }

  TH2F *hr2 = new TH2F("hr2",title1, 2,xmin1,xmax1,2,ymin1,ymax1);
  //hr2->SetXTitle("Centrality [%] ");
  hr2->SetYTitle(Form("V_{%i}              ", harmonic));
  hr2->GetYaxis()->SetTitleOffset(1.5);

  TH2F *hr2_rat = new TH2F("hr2_rat",title1, 2,xmin1_rat,xmax1_rat,2,ymin1_rat,ymax1_rat);
  hr2_rat->SetXTitle("Centrality [%]");
  hr2_rat->SetYTitle("Ratio   ");
  hr2_rat->GetXaxis()->SetLabelSize(0.1);
  hr2_rat->GetYaxis()->SetLabelSize(0.1);
  hr2_rat->GetXaxis()->SetTitleSize(0.1);
  hr2_rat->GetYaxis()->SetTitleSize(0.1);
  hr2_rat->GetYaxis()->SetTitleOffset(0.5);
  hr2_rat->GetXaxis()->SetTitleOffset(1.0);
  
  TLine *line = new TLine(0,1,80,1);
  line->SetLineColor(kBlue);
  line->SetLineStyle(1);

  TLine *line_up = new TLine(0,1.1,80,1.1);
  line_up->SetLineColor(kBlack);
  line_up->SetLineStyle(2);

  TLine *line_up_2 = new TLine(0,1.2,80,1.2);
  line_up_2->SetLineColor(kBlack);
  line_up_2->SetLineStyle(2);

  TLine *line_down = new TLine(0,0.9,80,0.9);
  line_down->SetLineColor(kBlack);
  line_down->SetLineStyle(2);

  TLine *line_down_2 = new TLine(0,0.8,80,0.8);
  line_down_2->SetLineColor(kBlack);
  line_down_2->SetLineStyle(2);


  TLatex* harmMark = new TLatex(5, TLatex_y, Form("Au+Au #sqrt{s_{NN}}= 39 GeV, #Delta#eta-gap=0.1,  0.2<p_{T}<1.5 GeV/c"));
  harmMark->SetTextFont(47);
  harmMark->SetTextSize(21);

  TLatex* harmMark_2 = new TLatex(5, TLatex_y, Form("Au+Au #sqrt{s_{NN}}= 39 GeV,%s ,  0.2<p_{T}<1.5 GeV/c", particle[nPID]));
  harmMark_2->SetTextFont(47);
  harmMark_2->SetTextSize(21);

  //TCanvas *c1_1 = new TCanvas("c1_1","Flow analysis results",150,10,800,650);
  pad_up -> cd();
  hr2->Draw();
  harmMark->Draw("same");

  pad_rat ->cd();
  hr2_rat->Draw();
  line->Draw("same");
  line_up->Draw("same");
  line_down->Draw("same");
  line_up_2->Draw("same");
  line_down_2->Draw("same");

  TCanvas *c1_2 = new TCanvas("c1_2","Flow analysis results",150,10,800,650);
  hr2->Draw();
  harmMark_2->Draw("same");

  pad_up->cd();

  TProfile2D *p_v2_1_pos = (TProfile2D*)f -> Get(Form("v%i_cent_prof_n_2_%s_pos", harmonic, PID));
  TProfile2D *p_v2_1_neg = (TProfile2D*)f -> Get(Form("v%i_cent_prof_n_2_%s_neg", harmonic, PID));
  
  for(Int_t i=0; i<9; i++){
    y_pos[8-i] = p_v2_1_pos->GetBinContent(p_v2_1_pos->FindBin(i)) ;
    y_err_pos[8-i] = p_v2_1_pos->GetBinError(p_v2_1_pos->FindBin(i));
  }
  TGraphErrors *gr1_pos = new TGraphErrors(9, x_centrality, y_pos, x_centrality_err, y_err_pos);
  gr1_pos->SetMarkerColor(kBlue);
  gr1_pos->SetLineColor(kBlue);
  gr1_pos->SetMarkerStyle(27);
  gr1_pos->SetMarkerSize(1.4);
  gr1_pos->Draw("P");

  TF1 *f1 = new TF1("f1","pol6",0,80);
  f1 -> SetLineColor(kBlue);
  f1 -> SetLineWidth(1);
  //gr1_pos -> Fit(f1,"R");
  
  for(Int_t i=0; i<9; i++){
    y[8-i] = p_v2_1_neg->GetBinContent(p_v2_1_neg->FindBin(i)) ;
    y_err[8-i] = p_v2_1_neg->GetBinError(p_v2_1_neg->FindBin(i));

    x_coord_1[8-i] = x_centrality[8-i];
    y_coord_1[8-i] = (Double_t)(p_v2_1_neg->GetBinContent(p_v2_1_neg->FindBin(i)) / y_pos[8-i]);
    y_coord_1err[8-i] = TMath::Sqrt( pow( (1.0 / y_pos[8-i]) , 2 )*pow( y_err[8-i],2 ) 
                      + pow(( y[8-i] / pow(y_pos[8-i],2) ),2 ) * pow( y_err_pos[8-i],2 ) );


  }

  TGraphErrors *gr1_neg = new TGraphErrors(9, x_centrality, y, x_centrality_err, y_err);
  gr1_neg->SetMarkerColor(kRed);
  gr1_neg->SetLineColor(kRed);
  gr1_neg->SetMarkerStyle(33);
  gr1_neg->SetMarkerSize(1.4);
  gr1_neg->Draw("P");

  TLegend *legC2_1 = new TLegend(0.80,0.6,0.89,.79);
  legC2_1->AddEntry(gr1_neg,Form("%s ", particle[nPID+1]),"p");
  legC2_1->AddEntry(gr1_pos,Form("%s ", particle[nPID]),"p");
  legC2_1->SetFillColor(kWhite);
  legC2_1->SetBorderSize(0);
  legC2_1->Draw();

  pad_rat->cd();

  TGraphErrors *gr1_ratio = new TGraphErrors(9, x_coord_1, y_coord_1, x_centrality_err, y_coord_1err);
  gr1_ratio->SetMarkerColor(kRed);
  gr1_ratio->SetLineColor(kRed);
  gr1_ratio->SetMarkerStyle(33);
  gr1_ratio->SetMarkerSize(1.4);
  gr1_ratio->Draw("P");

  c1_2->cd();

  TLegend *legC2_2 = new TLegend(0.74,0.3,0.89,.8);
  legC2_2->SetFillColor(kWhite);
  legC2_2->SetBorderSize(0);

  for(Int_t j=0; j < ngap; j++){
    
    TProfile2D *p_v2_2 = (TProfile2D*)f -> Get(Form("v%i_cent_prof_n_%i_%s%s", harmonic, j+1, PID, charge));
    for(Int_t i=0; i<9; i++){
      y[8-i] = p_v2_2->GetBinContent(p_v2_2->FindBin(i));
      y_err[8-i] = p_v2_2->GetBinError(p_v2_2->FindBin(i));
    }

    TGraphErrors *gr1_neg = new TGraphErrors(9, x_centrality, y, x_centrality_err, y_err);
    gr1_neg->SetMarkerColor(color[j]);
    gr1_neg->SetLineColor(color[j]);
    gr1_neg->SetMarkerStyle(style[j]);
    gr1_neg->SetMarkerSize(1.4);
    gr1_neg->Draw("P");

    legC2_2->AddEntry(gr1_neg,Form("#Delta#eta-gap=%.2f", 2.0*n[j]),"p");

  }

  legC2_2->Draw();

  c1 -> SaveAs(Form("/home/demanov/NIRS/Diplom/identified_hadrons/pict/v%i_%s_centr.png" ,harmonic,PID));
  c1_2 -> SaveAs(Form("/home/demanov/NIRS/Diplom/identified_hadrons/pict/%s/v%i/v%i_%s_%s_centr_gap.png", charge, harmonic, harmonic,PID, charge));

}

void Res(){

  makeplotstyle();

  Double_t x_centrality[9]  = {2.5, 7.5, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0,};
  Double_t ex[9]={0.0};

  Double_t y_star[9]={0.324088, 0.420438, 0.505109, 0.531387, 0.487591, 0.408759, 0.309489, 0.210219, 0.145985};
  Double_t y_star_err[9]={0.0};

  Double_t Res_2[9] = { 0.0 };
  Double_t Res_2_err[9] = { 0.0 };
  Double_t Res_3[9] = { 0.0 };
  Double_t Res_3_err[9] = { 0.0 };

  TFile *f3 = new TFile("answer_all_flat.root","READ");
  
  TFile *fstyle = new TFile("style.root");
  TStyle *tsty = (TStyle *)fstyle->Get("PlottingInStyle");
  tsty->cd();

  gROOT->ForceStyle();
  gStyle->SetOptStat(0);

  char title1[800];
  sprintf(title1,"");
  float xmin1=0.0;
  float xmax1=90.0;
  float ymin1=0.0;
  float ymax1=0.7;
  float ymax2=0.4;

  TLatex* harmMark = new TLatex(5, 0.63, "Au+Au #sqrt{s_{NN}}= 39 GeV, 0.2 < p < 5.0 GeV/c");
        harmMark->SetTextFont(47);
        harmMark->SetTextSize(21);

  TLatex* harmMark1 = new TLatex(5, 0.35, "Au+Au #sqrt{s_{NN}}= 39 GeV, 0.2 < p < 5.0 GeV/c");
        harmMark1->SetTextFont(47);
        harmMark1->SetTextSize(21);

  TH2F *hr1 = new TH2F("hr1",title1, 9, xmin1, xmax1, 2, ymin1, ymax1);
  hr1->SetXTitle("Centrality  ");
  hr1->SetYTitle("Res     ");

  TH2F *hr2 = new TH2F("hr2",title1, 9, xmin1, xmax1, 2, ymin1, ymax1);
  hr2->SetXTitle("Centrality  ");
  hr2->SetYTitle("R_{2}     ");

  TH2F *hr3 = new TH2F("hr3",title1, 9, xmin1, xmax1, 2, ymin1, ymax2);
  hr3->SetXTitle("Centrality  ");
  hr3->SetYTitle("R_{3}     ");

  TCanvas *can11 = new TCanvas("can11","Flow analysis results",150,10,800,650);
  hr1->Draw();

  TCanvas *can22 = new TCanvas("can22","Flow analysis results",150,10,800,650);
  hr2->Draw();

  TCanvas *can33 = new TCanvas("can33","Flow analysis results",150,10,800,650);
  hr3->Draw();

  TProfile *p_Res_2[ngap];
  TProfile *p_Res_3[ngap];

  TLegend *legC2_2 = new TLegend(0.65,0.6,0.89,.8);
  legC2_2->SetFillColor(kWhite);
  legC2_2->SetBorderSize(0);

  for (Int_t i = 0; i < ngap; i++){
    p_Res_2[i] = (TProfile*)f3 -> Get(Form("p_Res_2_n_%i",i+1));
    p_Res_3[i] = (TProfile*)f3 -> Get(Form("p_Res_3_n_%i",i+1));
  }

  for(Int_t i=1; i<2; i++){
    for(Int_t j=0; j<9; j++){
      Res_2[8-j] = TMath::Sqrt( p_Res_2[i]->GetBinContent(p_Res_2[i] -> FindBin( j )));
      Res_2_err[8-j] = p_Res_2[i]->GetBinError(p_Res_2[i] -> FindBin( j )) / (2.0* TMath::Sqrt( p_Res_2[i]->GetBinContent(p_Res_2[i] -> FindBin( j ))));
      Res_3[8-j] = TMath::Sqrt( p_Res_3[i]->GetBinContent(p_Res_3[i] -> FindBin( j )));
      Res_3_err[8-j] = p_Res_3[i]->GetBinError(p_Res_3[i] -> FindBin( j )) / (2.0* TMath::Sqrt( p_Res_3[i]->GetBinContent(p_Res_3[i] -> FindBin( j ))));
    }

    TGraphErrors *gr_R2 = new TGraphErrors(9, x_centrality, Res_2, ex, Res_2_err);
    gr_R2 -> SetMarkerColor(color[i]);
    gr_R2 -> SetLineColor(color[i]);
    gr_R2 -> SetMarkerStyle(style_2[i]);
    gr_R2 -> SetMarkerSize(1.4);

    TGraphErrors *gr_star = new TGraphErrors(9, x_centrality, y_star, ex, y_star_err);
    gr_star -> SetMarkerColor(color[i+1]);
    gr_star -> SetLineColor(color[i+1]);
    gr_star -> SetMarkerStyle(style[i+2]);
    gr_star -> SetMarkerSize(1.4);
    
    TGraphErrors *gr_R3 = new TGraphErrors(9, x_centrality, Res_3, ex, Res_3_err);
    gr_R3 -> SetMarkerColor(color[i+3]);
    gr_R3 -> SetLineColor(color[i+3]);
    gr_R3 -> SetMarkerStyle(style[i]);
    gr_R3 -> SetMarkerSize(1.4);

    can11->cd();
    gr_star-> Draw("P");
    gr_R2 -> Draw("P");
    gr_R3 -> Draw("P");

    can22->cd();
    gr_star-> Draw("P");
    gr_R2 -> Draw("P");

    can33->cd();
    gr_R3 -> Draw("P");

    //legC2_2->AddEntry(gr_R2,Form("#Delta#eta-gap=%.2f", 2.0*n[i]),"p");
    legC2_2->AddEntry(gr_star,"R_{2} star (#eta-sub)","p");
    legC2_2->AddEntry(gr_R2,"R_{2}, #Delta#eta-gap=0.1","p");
    legC2_2->AddEntry(gr_R3,"R_{3}, #Delta#eta-gap=0.1","p");

  }
  can11->cd();
  legC2_2->Draw();
  harmMark->Draw("same");
  can22->cd();
  legC2_2->Draw();
  harmMark->Draw("same");
  can33->cd();
  legC2_2->Draw();
  harmMark1->Draw("same");

  can11 -> SaveAs("/home/demanov/NIRS/Diplom/identified_hadrons/pict/R23_2.png");
  can22 -> SaveAs("/home/demanov/NIRS/Diplom/identified_hadrons/pict/R2_2.png");
  can33 -> SaveAs("/home/demanov/NIRS/Diplom/identified_hadrons/pict/R3_2.png");
}

void v3_pt( Double_t pt39GeV[], Double_t bin[] , Int_t cent_bin_1 ,Int_t cent_bin_2,  const Char_t *PID,  const Char_t *charge){

  makeplotstyle();

  TFile *f = new TFile("answer_flow.root","READ");

  Int_t cent[10]={80,70,60,50,40,30,20,10,5,0};
  
  Double_t x[col]={0.0};
  Double_t y[col]={0.0};
  Double_t y_err[col]={0.0};
  Double_t y_rat[col]={0.0};
  Double_t y_rat_err[col]={0.0};
  
  
  Double_t pt39GeV_err[col] = {0};

  
  Double_t x_coord_1[col], y_coord_1[col];
  Double_t x_coord_1err[col] = {0.0};
  Double_t y_coord_1err[col] = {0.0};
  
  TFile *fstyle = new TFile("style.root");
  TStyle *tsty = (TStyle *)fstyle->Get("PlottingInStyle");
  tsty->cd();

  gROOT->ForceStyle();
  gStyle->SetOptStat(0);

  char title1[800];
  sprintf(title1,"");

  float xmin1=0.0;
  float xmax1=4.6;
  float ymin1=0.00;
  float ymax1=0.12;

  float xmin1_rat=0.0;
  float xmax1_rat=4.6;
  float ymin1_rat=0.8;
  float ymax1_rat=1.2;

  Int_t nPID;

  if( strncmp(PID, "Proton",6)==0 ){
    if(strncmp(charge, "_pos",3)==0){
      nPID=4;
    }
    else{
      nPID=5;
    }
    if(strncmp(charge, "",1)==0){
      nPID=8;
    }
  }
  if( strncmp(PID, "Kaon",4)==0 ){
    if(strncmp(charge, "_pos",3)==0){
      nPID=2;
    }
    else{
      nPID=3;
    }
    if(strncmp(charge, "",1)==0){
      nPID=7;
    }
  }
  if( strncmp(PID, "Pion",4)==0 ){
    if(strncmp(charge, "_pos",3)==0){
      nPID=0;
    }
    else{
      nPID=1;
    }
    if(strncmp(charge, "",1)==0){
      nPID=6;
    }
  }

  TH2F *hr2 = new TH2F("hr2",title1, 2,xmin1,xmax1,2,ymin1,ymax1);
  hr2->SetXTitle("P_{T} (GeV/c) ");
  hr2->SetYTitle("V_{3}          ");

  TH2F *hr2_rat = new TH2F("hr2_rat",title1, 2,xmin1_rat,xmax1_rat,2,ymin1_rat,ymax1_rat);
  hr2_rat->SetXTitle("P_{T} (GeV/c) ");
  hr2_rat->SetYTitle("Ratio #Delta#eta_{i} /  #Delta#eta=1.0         ");

  TLatex* harmMark = new TLatex(0.2, 0.1, Form("Au+Au #sqrt{s_{NN}}= 39 GeV, %i-%i %%, %s ", cent[cent_bin_2], cent[cent_bin_1-1], particle[nPID]));
  harmMark->SetTextFont(47);
  harmMark->SetTextSize(21);

  TLatex* harmMark_ratio = new TLatex(0.2, 1.15, Form("Au+Au #sqrt{s_{NN}}= 39 GeV, %i-%i %%, %s ", cent[cent_bin_2], cent[cent_bin_1-1], particle[nPID]));
  harmMark_ratio->SetTextFont(47);
  harmMark_ratio->SetTextSize(21);

  TLine *line = new TLine(0,1,4,1);
  line->SetLineColor(kBlack);
  line->SetLineStyle(2);

  TCanvas *c1_1 = new TCanvas("c1_1","Flow analysis results",150,10,800,650);
  hr2->Draw();
  harmMark->Draw("same");

  TCanvas *c1_2 = new TCanvas("c1_2","Flow analysis results",150,10,800,650);
  hr2->Draw();
  harmMark->Draw("same");

  TCanvas *c1_2_ratio = new TCanvas("c1_2_ratio","Flow analysis results",150,10,800,650);
  hr2_rat->Draw();
  harmMark_ratio->Draw("same");

  c1_1->cd();

  TProfile2D *p_v2_1 = (TProfile2D*)f -> Get(Form("v3_pt_prof_n_2_%s%s", PID, charge));
  TH1D *h1 = (TH1D*) p_v2_1->ProfileY("prof1",cent_bin_1,cent_bin_2)->Rebin(col,"h1",bin);
  for(Int_t i=0; i<col; i++){
    y[i] = h1->GetBinContent(h1->FindBin(pt39GeV[i]));
    x[i] = h1->GetBinCenter(h1->FindBin(pt39GeV[i]));
    y_err[i] = h1->GetBinError(h1->FindBin(pt39GeV[i]));
  }
  TGraphErrors *gr1 = new TGraphErrors(col, x, y, pt39GeV_err, y_err);
  gr1->SetMarkerColor(kRed);
  gr1->SetLineColor(kRed);
  gr1->SetMarkerStyle(33);
  gr1->SetMarkerSize(1.4);
  gr1->Draw("P");
  
  TLegend *legC3 = new TLegend(0.74,0.7,0.89,.89);
  legC3->AddEntry(gr1,"#Delta#eta-gap=0.1","p");
  legC3->SetFillColor(kWhite);
  legC3->SetBorderSize(0); legC3->Draw();


  TLegend *legC3_2 = new TLegend(0.74,0.3,0.89,.8);
  legC3_2->SetFillColor(kWhite);
  legC3_2->SetBorderSize(0);

  c1_2->cd();

  TProfile2D *p_v2_2 = (TProfile2D*)f -> Get(Form("v3_pt_prof_n_5_%s%s", PID, charge));
  TH1D *h2 = (TH1D*) p_v2_2->ProfileY("prof2",cent_bin_1,cent_bin_2)->Rebin(col,"h2",bin);

  for(Int_t i=0; i<col; i++){
    y_rat[i] = h2->GetBinContent(h2->FindBin(pt39GeV[i]));
    y_rat_err[i] = h2->GetBinError(h2->FindBin(pt39GeV[i]));
  }
  TGraphErrors *gr_r = new TGraphErrors(col, x, y_rat, pt39GeV_err, y_rat_err);
  gr_r->SetMarkerColor(kBlack);
  gr_r->SetLineColor(kBlack);
  gr_r->SetMarkerStyle(33);
  gr_r->SetMarkerSize(1.4);
  gr_r->Draw("P");

  TF1 *fit = new TF1("fit","pol9",0.3,1.5);
  fit -> SetLineWidth(1);
  gr_r -> Fit(fit,"R");

  Int_t flag=0.0;

  for(Int_t i=0; i < ngap; i++){
    flag++;
    TProfile2D *p_v2_test = (TProfile2D*)f -> Get(Form("v3_pt_prof_n_%i_%s%s", i+1, PID, charge));
    TH1D *h1 = (TH1D*) p_v2_test->ProfileY(Form("prof1_%i",i+1),cent_bin_1,cent_bin_2)->Rebin(col,Form("h1_%i",i+1),bin);
    for(Int_t j=0; j<col; j++){
      y[j] = h1->GetBinContent(h1->FindBin(pt39GeV[j]));
      x[i] = h1 ->GetBinCenter(h1->FindBin( pt39GeV[i] ));
      y_err[j] = h1->GetBinError(h1->FindBin(pt39GeV[j]));

      x_coord_1[j] = h1->GetBinCenter(h1->FindBin(pt39GeV[j]));
      y_coord_1[j] = (Double_t)(h1->GetBinContent(h1->FindBin(pt39GeV[j])) / (Double_t)gr_r->Eval(x_coord_1[j]));
      y_coord_1err[j] = TMath::Sqrt( pow( (1.0 / ( (Double_t)gr_r->Eval(x_coord_1[j]))) , 2 )*pow( y_err[j],2 ) 
                      + pow(( y[j] / pow((Double_t)gr_r->Eval(x_coord_1[j]),2) ),2 ) * pow( y_rat_err[j],2 ) );

      std::cout<< y[j] <<std::endl;

    }
  std::cout<< "\n\n";
    if(flag==5){
      for(int j=0; j<col; j++){
        y_coord_1err[j]=0.0;
      }
    }

    c1_2->cd();
    TGraphErrors *gr1 = new TGraphErrors(col, x, y, pt39GeV_err, y_err);
    gr1->SetMarkerColor(color[i]);
    gr1->SetLineColor(color[i]);
    gr1->SetMarkerStyle(style[i]);
    gr1->SetMarkerSize(1.4);
    gr1->Draw("P");

    c1_2_ratio->cd();
    TGraphErrors *gr1_ratio = new TGraphErrors(col, x_coord_1, y_coord_1, x_coord_1err, y_coord_1err);
    gr1_ratio->SetMarkerColor(color[i]);
    gr1_ratio->SetLineColor(color[i]);
    gr1_ratio->SetMarkerStyle(style[i]);
    gr1_ratio->SetMarkerSize(1.4);
    gr1_ratio->Draw("P");

    legC3_2->AddEntry(gr1,Form("#Delta#eta-gap=%.2f", 2.0*n[i]),"p");
  
  }

  c1_2->cd();
  legC3_2->Draw();
  c1_2_ratio->cd();
  legC3_2->Draw();
  line->Draw("same");

  c1_1 -> SaveAs(Form("/home/demanov/NIRS/Diplom/identified_hadrons/pict/%s/v3/v3_%s%s_%i_%i.png", charge, PID, charge, cent_bin_1, cent_bin_2));
  c1_2 -> SaveAs(Form("/home/demanov/NIRS/Diplom/identified_hadrons/pict/%s/v3/v3_%s%s_%i_%i_gap.png", charge, PID, charge, cent_bin_1, cent_bin_2));
  //c1_2_ratio -> SaveAs(Form("/home/demanov/NIRS/Diplom/identified_hadrons/pict/%s/v3/v3_%s_%i_%i_gap_ratio.png", charge, PID,cent_bin_1, cent_bin_2));

}


void v2_pt(Double_t v2_cent[], Double_t v2_cent_err[], Double_t pt39GeV[], Double_t bin[], Int_t cent_bin_1,Int_t cent_bin_2, const Char_t *PID,  const Char_t *charge  ){

  makeplotstyle();

  TFile *f = new TFile("answer_flow.root","READ");
  
  Int_t cent[10]={80,70,60,50,40,30,20,10,5,0};
  
  Double_t x[col]={0.0};
  Double_t y[col]={0.0};
  Double_t y_err[col]={0.0};

  Double_t x_coord_1[col], y_coord_1[col];
  Double_t x_coord_1err[col] = {0.0};
  Double_t y_coord_1err[col] = {0.0};
  
  
  Double_t pt39GeV_err[col] = {0};
  
  TFile *fstyle = new TFile("style.root");
  TStyle *tsty = (TStyle *)fstyle->Get("PlottingInStyle");
  tsty->cd();

  gROOT->ForceStyle();
  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1","Flow ",150,10,800,650);
  TPad *pad_up = new TPad("pad_up", "", 0.0, 0.33, 1.0, 1.0, 0, 0, 0);
  pad_up->SetFrameBorderMode(0);
  pad_up->SetBottomMargin(0);
  TPad *pad_rat = new TPad("pad_rat", "", 0.0, 0.0, 1.0, 0.33, 0, 0, 0);
  pad_rat->SetFrameBorderMode(0);
  pad_rat->SetTopMargin(0.005);
  pad_rat->SetBottomMargin(0.25);
 

  c1->GetFrame()->SetBorderSize(115);
  c1->cd();
  pad_up->Draw();
  pad_rat->Draw();
  gStyle->SetOptStat(0);

  c1->GetFrame()->SetFillColor(21);
  c1->GetFrame()->SetBorderSize(115);

  gStyle->SetOptStat(0);
  gStyle->SetPalette(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetHistFillColor(10);
  gStyle->SetHistFillStyle(0);
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  c1->SetBorderMode(0);

  c1->SetFillColor(0);

  char title1[800];
  sprintf(title1,"");

  float xmin1=0.0;
  float xmax1=4.6;
  float ymin1=0.00;
  float ymax1=0.28;
  float TLatex_y=0.26;

  if( strncmp(PID, "Proton",6)==0 ){
    ymax1 = 0.32;
    TLatex_y=0.30;
  }

  float xmin1_rat=0.0;
  float xmax1_rat=4.6;
  float ymin1_rat=0.8;
  float ymax1_rat=1.2;

  Int_t nPID;

  if( strncmp(PID, "Proton",6)==0 ){
    if(strncmp(charge, "_pos",3)==0){
      nPID=4;
    }
    else{
      nPID=5;
    }
    if(strncmp(charge, "",1)==0){
      nPID=8;
    }
  }
  if( strncmp(PID, "Kaon",4)==0 ){
    if(strncmp(charge, "_pos",3)==0){
      nPID=2;
    }
    else{
      nPID=3;
    }
    if(strncmp(charge, "",1)==0){
      nPID=7;
    }
  }
  if( strncmp(PID, "Pion",4)==0 ){
    if(strncmp(charge, "_pos",3)==0){
      nPID=0;
    }
    else{
      nPID=1;
    }
    if(strncmp(charge, "",1)==0){
      nPID=6;
    }
  }

  TH2F *hr2 = new TH2F("hr2",title1, 2,xmin1,xmax1,2,ymin1,ymax1);
  hr2->SetXTitle("P_{T} (GeV/c) ");
  hr2->SetYTitle("V_{2}          ");

  TH2F *hr2_rat = new TH2F("hr2_rat",title1, 2,xmin1_rat,xmax1_rat,2,ymin1_rat,ymax1_rat);
  hr2_rat->SetXTitle("p_{T} [GeV/c]");
  hr2_rat->SetYTitle("Ratio   ");
  hr2_rat->GetXaxis()->SetLabelSize(0.1);
  hr2_rat->GetYaxis()->SetLabelSize(0.1);
  hr2_rat->GetXaxis()->SetTitleSize(0.1);
  hr2_rat->GetYaxis()->SetTitleSize(0.1);
  hr2_rat->GetYaxis()->SetTitleOffset(0.5);
  hr2_rat->GetXaxis()->SetTitleOffset(1.0);

  TLatex* harmMark = new TLatex(0.2, TLatex_y, Form("Au+Au #sqrt{s_{NN}}= 39 GeV, %i-%i %%, %s", cent[cent_bin_2], cent[cent_bin_1-1], particle[nPID]));
  harmMark->SetTextFont(47);
  harmMark->SetTextSize(21);
  
  TLine *line = new TLine(0,1,4,1);
  line->SetLineColor(kBlue);
  line->SetLineStyle(1);

  TLine *line_up = new TLine(0,1.1,4,1.1);
  line_up->SetLineColor(kBlack);
  line_up->SetLineStyle(2);

  TLine *line_down = new TLine(0,0.9,4,0.9);
  line_down->SetLineColor(kBlack);
  line_down->SetLineStyle(2);

  pad_up->cd();
  hr2->Draw();
  harmMark->Draw("same");

  pad_rat->cd();
  hr2_rat->Draw();
  line->Draw("same");
  line_up->Draw("same");
  line_down->Draw("same");

  pad_up->cd();

  TGraphErrors *gr_star = new TGraphErrors(col, pt39GeV, v2_cent, pt39GeV_err, v2_cent_err);
  gr_star -> SetMarkerColor(kBlue);
  gr_star -> SetLineColor(kBlue);
  gr_star -> SetMarkerStyle(21);
  gr_star -> SetMarkerSize(1.4);
  gr_star -> Draw("P");

  TF1 *f1 = new TF1("f1","pol5",0.2,3.5);
  f1 -> SetLineColor(kBlue);
  f1 -> SetLineWidth(1);
  gr_star -> Fit(f1,"R");

  TProfile2D *p_v2_1 = (TProfile2D*)f -> Get(Form("v2_pt_prof_n_2_%s%s", PID, charge));
  TH1D *h1 = (TH1D*) p_v2_1->ProfileY("prof1",cent_bin_1,cent_bin_2)->Rebin(col,"h1",bin);
  for(Int_t i=0; i<col; i++){
    y[i] = h1->GetBinContent(h1->FindBin(pt39GeV[i]));
    x[i] =h1 ->GetBinCenter(h1->FindBin( pt39GeV[i] ));
    y_err[i] = h1->GetBinError(h1->FindBin(pt39GeV[i]));

    std::cout<< y[i]<<"\t\t"<< v2_cent[i] <<std::endl;

    x_coord_1[i] = h1->GetBinCenter(h1->FindBin(pt39GeV[i]));
    y_coord_1[i] = (Double_t)(h1->GetBinContent(h1->FindBin(pt39GeV[i])) / (Double_t)gr_star->Eval(x_coord_1[i]));
    y_coord_1err[i] = TMath::Sqrt( pow( (1.0 / ( (Double_t)gr_star->Eval(x_coord_1[i]))) , 2 )*pow( y_err[i],2 ) 
                      + pow(( y[i] / pow((Double_t)gr_star->Eval(x_coord_1[i]),2) ),2 ) * pow( v2_cent_err[i],2 ) );
  }
  TGraphErrors *gr1 = new TGraphErrors(col, x, y, pt39GeV_err, y_err);
  gr1->SetMarkerColor(kRed);
  gr1->SetLineColor(kRed);
  gr1->SetMarkerStyle(33);
  gr1->SetMarkerSize(1.4);
  gr1->Draw("P");

  TLegend *legC2_1 = new TLegend(0.74,0.7,0.89,.89);
  legC2_1->AddEntry(gr_star,"star","p");
  legC2_1->AddEntry(gr1,"#Delta#eta-gap=0.1","p");
  legC2_1->SetFillColor(kWhite);
  legC2_1->SetBorderSize(0); 
  legC2_1->Draw();

  pad_rat->cd();

  TGraphErrors *gr1_ratio = new TGraphErrors(col, x_coord_1, y_coord_1, x_coord_1err, y_coord_1err);
  gr1_ratio->SetMarkerColor(kRed);
  gr1_ratio->SetLineColor(kRed);
  gr1_ratio->SetMarkerStyle(33);
  gr1_ratio->SetMarkerSize(1.4);
  gr1_ratio->Draw("P");
  
  line->Draw("same");
  legC2_1->Draw();

  ///-------------------------------------------------
  TCanvas *c2 = new TCanvas("c2","Flow ",150,10,800,650);
  TPad *pad_up_2 = new TPad("pad_up_2", "", 0.0, 0.33, 1.0, 1.0, 0, 0, 0);
  pad_up_2->SetFrameBorderMode(0);
  pad_up_2->SetBottomMargin(0);
  TPad *pad_rat_2 = new TPad("pad_rat_2", "", 0.0, 0.0, 1.0, 0.33, 0, 0, 0);
  pad_rat_2->SetFrameBorderMode(0);
  pad_rat_2->SetTopMargin(0.005);
  pad_rat_2->SetBottomMargin(0.25);
 

  c2->GetFrame()->SetBorderSize(115);
  c2->cd();
  pad_up_2->Draw();
  pad_rat_2->Draw();
  gStyle->SetOptStat(0);

  c2->GetFrame()->SetFillColor(21);
  c2->GetFrame()->SetBorderSize(115);

  gStyle->SetOptStat(0);
  gStyle->SetPalette(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetHistFillColor(10);
  gStyle->SetHistFillStyle(0);
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  c2->SetBorderMode(0);

  c2->SetFillColor(0);

  pad_up_2->cd();
  hr2->Draw();
  harmMark->Draw("same");

  pad_rat_2->cd();
  hr2_rat->Draw();
  line->Draw("same");
  line_up->Draw("same");
  line_down->Draw("same");

  pad_up_2->cd();

  TLegend *legC2_2 = new TLegend(0.74,0.3,0.89,.8);
  legC2_2->SetFillColor(kWhite);
  legC2_2->SetBorderSize(0);
  legC2_2->AddEntry(gr_star,"star","p");

  pad_up_2 -> cd();
  gr_star -> Draw("P");

  for(Int_t j=0; j < ngap; j++){
    
    TProfile2D *p_v2_test = (TProfile2D*)f -> Get(Form("v2_pt_prof_n_%i_%s%s", j+1, PID, charge));
    TH1D *h1 = (TH1D*) p_v2_test->ProfileY(Form("prof1_%i",j+1),cent_bin_1,cent_bin_2)->Rebin(col,Form("h1_%i",j+1),bin);   
    
    for(Int_t i=0; i<col; i++){
      y[i] = h1->GetBinContent(h1->FindBin(pt39GeV[i]));
       x[i] =h1 ->GetBinCenter(h1->FindBin( pt39GeV[i] ));
      y_err[i] = h1->GetBinError(h1->FindBin(pt39GeV[i]));

      x_coord_1[i] = h1->GetBinCenter(h1->FindBin(pt39GeV[i]));
      y_coord_1[i] = (Double_t)(h1->GetBinContent(h1->FindBin(pt39GeV[i])) / (Double_t)gr_star->Eval(x_coord_1[i]));
      y_coord_1err[i] = TMath::Sqrt( pow( (1.0 / ( (Double_t)gr_star->Eval(x_coord_1[i]))) , 2 )*pow( y_err[i],2 ) 
                      + pow(( y[i] / pow((Double_t)gr_star->Eval(x_coord_1[i]),2) ),2 ) * pow( v2_cent_err[i],2 ) );
    }

    pad_up_2->cd();
    TGraphErrors *gr1 = new TGraphErrors(col, x, y, pt39GeV_err, y_err);
    gr1->SetMarkerColor(color[j]);
    gr1->SetLineColor(color[j]);
    gr1->SetMarkerStyle(style[j]);
    gr1->SetMarkerSize(1.4);
    gr1->Draw("P");

    legC2_2->AddEntry(gr1,Form("#Delta#eta-gap=%.2f", 2.0*n[j]),"p");

    pad_rat_2->cd();
    TGraphErrors *gr1_ratio = new TGraphErrors(col, x_coord_1, y_coord_1, x_coord_1err, y_coord_1err);
    gr1_ratio->SetMarkerColor(color[j]);
    gr1_ratio->SetLineColor(color[j]);
    gr1_ratio->SetMarkerStyle(style[j]);
    gr1_ratio->SetMarkerSize(1.4);
    gr1_ratio->Draw("P");
  
  }

  pad_up_2->cd();
  legC2_2->Draw();

  pad_rat_2->cd();
  legC2_2->Draw();
  line->Draw("same");

  c1 -> SaveAs(Form("/home/demanov/NIRS/Diplom/identified_hadrons/pict/%s/v2/v2_%s%s_%i_%i.png", charge, PID, charge, cent_bin_1, cent_bin_2));
  c2 -> SaveAs(Form("/home/demanov/NIRS/Diplom/identified_hadrons/pict/%s/v2/v2_%s%s_%i_%i_gap.png", charge, PID, charge, cent_bin_1, cent_bin_2));

}


void v_cent_rat(Int_t harmonic,const Char_t *PID){

  makeplotstyle();

  Double_t x_centrality[9]  = {2.5, 7.5, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0,};
  Double_t x_centrality_err[9]  = {0.0};
  
  Double_t y_39[9]={0.0};
  Double_t y_39_err[9]={0.0};

  Double_t y_19[9]={0.0};
  Double_t y_19_err[9]={0.0};

  const Int_t n1=9;
  Double_t x_coord_1[n1], y_coord_1[n1];
  Double_t x_coord_1err[n1] = {0.0};
  Double_t y_coord_1err[n1] = {0.0};
  
  TFile *fstyle = new TFile("style.root");
  TStyle *tsty = (TStyle *)fstyle->Get("PlottingInStyle");
  tsty->cd();
  gStyle->SetOptStat(0);

  gROOT->ForceStyle();
  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1","Flow ",150,10,800,650);
  TPad *pad_up = new TPad("pad_up", "", 0.0, 0.35, 1.0, 1.0, 0, 0, 0);
  pad_up->SetFrameBorderMode(0);
  pad_up->SetBottomMargin(0);
  TPad *pad_rat = new TPad("pad_rat", "", 0.0, 0.0, 1.0, 0.35, 0, 0, 0);
  pad_rat->SetFrameBorderMode(0);
  pad_rat->SetTopMargin(0.005);
  pad_rat->SetBottomMargin(0.25);
 

  c1->GetFrame()->SetBorderSize(115);
  c1->cd();
  pad_up->Draw();
  pad_rat->Draw();
  gStyle->SetOptStat(0);

  c1->GetFrame()->SetFillColor(21);
  c1->GetFrame()->SetBorderSize(115);

  gStyle->SetOptStat(0);
  gStyle->SetPalette(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetHistFillColor(10);
  gStyle->SetHistFillStyle(0);
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  c1->SetBorderMode(0);

  c1->SetFillColor(0);

  char title1[800];
  sprintf(title1,"");

  float xmin1=0.0;
  float xmax1=100.0;
  float ymin1;
  float ymax1;

  float xmin1_rat=0.0;
  float xmax1_rat=100.0;
  float ymin1_rat=0.7;
  float ymax1_rat=1.3;

  float TLatex_y;

  if(harmonic == 3){
    ymin1=0.0;
    ymax1=2.0;
    TLatex_y=1.85;
    ymin1_rat=0.95;
    ymax1_rat=1.6;

  }
  if(harmonic == 2){
    ymin1=0.8;
    ymax1=1.2;
    TLatex_y=1.15;
  }

  pad_up -> cd();
  
  TLatex *harmMark;
  //hr2->SetXTitle("Centrality [%] ");
  if( strncmp(PID, "Pion",4)==0 ){
    TH2F *hr2 = new TH2F("hr2",title1, 2,xmin1,xmax1,2,ymin1,ymax1);
    hr2->GetYaxis()->SetTitleOffset(1.5);
    hr2->SetYTitle(Form("V_{%i}(#pi^{+})/V_{%i}(#pi^{-})              ", harmonic, harmonic));
    harmMark = new TLatex(5, TLatex_y, Form("Au+Au, #pi^{+} #pi^{-}, #Delta#eta-gap=0.1,  0.2<P_{T}<1.5 GeV/c"));
    hr2->Draw();

  }
  if( strncmp(PID, "Kaon",4)==0 ){
    TH2F *hr2 = new TH2F("hr2",title1, 2,xmin1,xmax1,2,ymin1,ymax1);
    hr2->GetYaxis()->SetTitleOffset(1.5);
    hr2->SetYTitle(Form("V_{%i}(K^{+})/V_{%i}(K^{-})              ", harmonic, harmonic));
    harmMark = new TLatex(5, TLatex_y, Form("Au+Au, K^{+} K^{-}, #Delta#eta-gap=0.1,  0.2<P_{T}<1.5 GeV/c"));
    hr2->Draw();
    if(harmonic == 3){
      ymax1_rat=2.5;
    }

  }
  if( strncmp(PID, "Proton",6)==0 ){
    ymin1=1.0;
    ymax1=1.6;
    if(harmonic == 3){
      ymin1=0.0;
      ymax1=2.0;
      TLatex_y=1.9;
      ymin1_rat=0.95;
      ymax1_rat=2.0;
    }
    
    TH2F *hr2 = new TH2F("hr2",title1, 2,xmin1,xmax1,2,ymin1,ymax1);
    hr2->GetYaxis()->SetTitleOffset(1.5);
    hr2->SetYTitle(Form("V_{%i}(p)/V_{%i}(#bar{p})              ", harmonic, harmonic));
    harmMark = new TLatex(5, ymax1-0.1, Form("Au+Au, p #bar{p}, #Delta#eta-gap=0.1,  0.4<P_{T}<2.0 GeV/c"));
    hr2->Draw();
  }

  harmMark->SetTextFont(47);
  harmMark->SetTextSize(21);
  harmMark->Draw("same");

  TH2F *hr2_rat = new TH2F("hr2_rat",title1, 2,xmin1_rat,xmax1_rat,2,ymin1_rat,ymax1_rat);
  hr2_rat->SetXTitle("Centrality [%]");
  hr2_rat->SetYTitle("Ratio        ");
  hr2_rat->GetXaxis()->SetLabelSize(0.1);
  hr2_rat->GetYaxis()->SetLabelSize(0.1);
  hr2_rat->GetXaxis()->SetTitleSize(0.1);
  hr2_rat->GetYaxis()->SetTitleSize(0.1);
  hr2_rat->GetYaxis()->SetTitleOffset(0.5);
  hr2_rat->GetXaxis()->SetTitleOffset(1.0);
  
  TLine *line = new TLine(0,1,80,1);
  line->SetLineColor(kRed);
  line->SetLineStyle(1);

  TLine *line_up = new TLine(0,1.1,80,1.1);
  line_up->SetLineColor(kBlack);
  line_up->SetLineStyle(2);

  TLine *line_up_2 = new TLine(0,1.2,80,1.2);
  line_up_2->SetLineColor(kBlack);
  line_up_2->SetLineStyle(2);

  TLine *line_down = new TLine(0,0.9,80,0.9);
  line_down->SetLineColor(kBlack);
  line_down->SetLineStyle(2);

  TLine *line_down_2 = new TLine(0,0.8,80,0.8);
  line_down_2->SetLineColor(kBlack);
  line_down_2->SetLineStyle(2);

  //TCanvas *c1_1 = new TCanvas("c1_1","Flow analysis results",150,10,800,650);

  pad_rat ->cd();
  hr2_rat->Draw();
  line->Draw("same");
 // line_up->Draw("same");
 // line_down->Draw("same");
 // line_up_2->Draw("same");
 // line_down_2->Draw("same");


  pad_up->cd();

  TFile *f = new TFile("answer_flow.root","READ");

  TProfile2D *p_v2_1_39_pos = (TProfile2D*)f -> Get(Form("v%i_cent_prof_n_2_%s_pos", harmonic, PID));
  TProfile2D *p_v2_1_39_neg = (TProfile2D*)f -> Get(Form("v%i_cent_prof_n_2_%s_neg", harmonic, PID));

  TFile *f2 = new TFile("answer_flow(pasha).root","READ");

  TProfile2D *p_v2_1_19_pos = (TProfile2D*)f2 -> Get(Form("v%i_cent_prof_n_2_%s_pos", harmonic, PID));
  TProfile2D *p_v2_1_19_neg = (TProfile2D*)f2 -> Get(Form("v%i_cent_prof_n_2_%s_neg", harmonic, PID));
  
  for(Int_t i=0; i<9; i++){

    y_39[8-i] = p_v2_1_39_pos->GetBinContent(p_v2_1_39_pos->FindBin(i)) / p_v2_1_39_neg->GetBinContent(p_v2_1_39_neg->FindBin(i));
    y_39_err[8-i] = TMath::Sqrt( pow( p_v2_1_39_pos->GetBinError(p_v2_1_39_pos->FindBin(i)) / p_v2_1_39_neg->GetBinContent(p_v2_1_39_neg->FindBin(i)), 2) + 
                                 pow( p_v2_1_39_pos->GetBinContent(p_v2_1_39_pos->FindBin(i))*p_v2_1_39_neg->GetBinError(p_v2_1_39_neg->FindBin(i)) / pow(p_v2_1_39_neg->GetBinContent(p_v2_1_39_neg->FindBin(i)) , 2) , 2) );
    
    y_19[8-i] = p_v2_1_19_pos->GetBinContent(p_v2_1_19_pos->FindBin(i)) / p_v2_1_19_neg->GetBinContent(p_v2_1_19_neg->FindBin(i));
    y_19_err[8-i] = TMath::Sqrt( pow( p_v2_1_19_pos->GetBinError(p_v2_1_19_pos->FindBin(i)) / p_v2_1_19_neg->GetBinContent(p_v2_1_19_neg->FindBin(i)), 2) + 
                                 pow( p_v2_1_19_pos->GetBinContent(p_v2_1_19_pos->FindBin(i))*p_v2_1_19_neg->GetBinError(p_v2_1_19_neg->FindBin(i)) / pow(p_v2_1_19_neg->GetBinContent(p_v2_1_19_neg->FindBin(i)) , 2) , 2) );
 
    y_coord_1[8-i] = y_39[8-i] / y_19[8-i];
    y_coord_1err[8-i] = TMath::Sqrt( pow( y_39_err[8-i] / y_19[8-i] , 2) + 
                                     pow( y_39[8-i]*y_19_err[8-i] / pow( y_19[8-i] , 2) , 2) );

  }
  TGraphErrors *gr1_1 = new TGraphErrors(9, x_centrality, y_39, x_centrality_err, y_39_err);
  gr1_1->SetMarkerColor(kBlue);
  gr1_1->SetLineColor(kBlue);
  gr1_1->SetMarkerStyle(33);
  gr1_1->SetMarkerSize(1.4);
  gr1_1->Draw("P");

  TGraphErrors *gr1_2 = new TGraphErrors(9, x_centrality, y_19, x_centrality_err, y_19_err);
  gr1_2->SetMarkerColor(kRed);
  gr1_2->SetLineColor(kRed);
  gr1_2->SetMarkerStyle(33);
  gr1_2->SetMarkerSize(1.4);
  gr1_2->Draw("P");

  TLegend *legC2_1 = new TLegend(0.74,0.6,0.89,.79);
  legC2_1->AddEntry(gr1_1,"#sqrt{S_{NN}}=39 GeV","p");
  legC2_1->AddEntry(gr1_2,"#sqrt{S_{NN}}=19.6 GeV","p");
  legC2_1->SetFillColor(kWhite);
  legC2_1->SetBorderSize(0);
  legC2_1->Draw();

  pad_rat->cd();

  TGraphErrors *gr1_ratio = new TGraphErrors(9, x_centrality, y_coord_1, x_centrality_err, y_coord_1err);
  gr1_ratio->SetMarkerColor(kBlue);
  gr1_ratio->SetLineColor(kBlue);
  gr1_ratio->SetMarkerStyle(33);
  gr1_ratio->SetMarkerSize(1.4);
  gr1_ratio->Draw("P");

  TLatex *harmMark2 = new TLatex(60, ymax1_rat-0.1, "RESULT(39 GeV) / RESULT(19.6 GeV)");
  harmMark2 -> Draw("same");



  c1 -> SaveAs(Form("/home/demanov/NIRS/Diplom/identified_hadrons/pict/rat_v%i_%s_centr.png" ,harmonic,PID));

}


void v_cent_gev(Int_t harmonic,const Char_t *PID,  const Char_t *charge){

  makeplotstyle();
  
  Double_t x_centrality[9]  = {2.5, 7.5, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0,};
  Double_t x_centrality_err[9]  = {0.0};
  
  Double_t y[9]={0.0};
  Double_t y_err[9]={0.0};

  Double_t y_pos[9]={0.0};
  Double_t y_err_pos[9]={0.0};

  const Int_t n1=9;
  Double_t x_coord_1[n1], y_coord_1[n1];
  Double_t x_coord_1err[n1] = {0.0};
  Double_t y_coord_1err[n1] = {0.0};
  
  TFile *fstyle = new TFile("style.root");
  TStyle *tsty = (TStyle *)fstyle->Get("PlottingInStyle");
  tsty->cd();
  gStyle->SetOptStat(0);

  gROOT->ForceStyle();
  gStyle->SetOptStat(0);

  char title1[800];
  sprintf(title1,"");

  float xmin1=0.0;
  float xmax1=100.0;
  float ymin1;
  float ymax1;

  float xmin1_rat=0.0;
  float xmax1_rat=100.0;
  float ymin1_rat=0.7;
  float ymax1_rat=1.3;

  float TLatex_y;

  Int_t nPID;

  if( strncmp(PID, "Proton",6)==0 ){
    if(strncmp(charge, "_pos",3)==0){
      nPID=4;
    }
    else{
      nPID=5;
    }
    if(strncmp(charge, "",1)==0){
      nPID=8;
    }
  }
  if( strncmp(PID, "Kaon",4)==0 ){
    if(strncmp(charge, "_pos",3)==0){
      nPID=2;
    }
    else{
      nPID=3;
    }
    if(strncmp(charge, "",1)==0){
      nPID=7;
    }
  }
  if( strncmp(PID, "Pion",4)==0 ){
    if(strncmp(charge, "_pos",3)==0){
      nPID=0;
    }
    else{
      nPID=1;
    }
    if(strncmp(charge, "",1)==0){
      nPID=6;
    }
  }

  if(harmonic == 3){
    ymin1=0.00;
    ymax1=0.06;
    TLatex_y=0.057;

  }
  if(harmonic == 2){
    ymin1=0.01;
    ymax1=0.13;
    TLatex_y=0.123;
  }

  TH2F *hr2 = new TH2F("hr2",title1, 2,xmin1,xmax1,2,ymin1,ymax1);
  //hr2->SetXTitle("Centrality [%] ");
  hr2->SetYTitle(Form("V_{%i}              ", harmonic));
  hr2->GetYaxis()->SetTitleOffset(1.5);

  TLatex* harmMark3 = new TLatex(5, TLatex_y, Form("Au+Au  #Delta#eta-gap=0.1,  0.4<P_{T}<2.0 GeV/c"));
  harmMark3->SetTextFont(47);
  harmMark3->SetTextSize(21);


  TCanvas *ctar = new TCanvas("ctar","Flow ",150,10,800,650);

  ctar->cd();
  hr2->Draw();
  harmMark3->Draw("same");

  TFile *f = new TFile("answer_flow.root","READ");
  TProfile2D *p_v2_1_sasha = (TProfile2D*)f -> Get(Form("v%i_cent_prof_n_2_%s%s", harmonic, PID, charge));
  
  TFile *f2 = new TFile("answer_flow(pasha).root","READ");
  TProfile2D *p_v2_1_pasha = (TProfile2D*)f2 -> Get(Form("v%i_cent_prof_n_2_%s%s", harmonic, PID, charge));
  
  for(Int_t i=0; i<9; i++){
    y_pos[8-i] = p_v2_1_sasha->GetBinContent(p_v2_1_sasha->FindBin(i)) ;
    y_err_pos[8-i] = p_v2_1_sasha->GetBinError(p_v2_1_sasha->FindBin(i));

    y[8-i] = p_v2_1_pasha->GetBinContent(p_v2_1_pasha->FindBin(i)) ;
    y_err[8-i] = p_v2_1_pasha->GetBinError(p_v2_1_pasha->FindBin(i));
  }
  TGraphErrors *gr1_sasha = new TGraphErrors(9, x_centrality, y_pos, x_centrality_err, y_err_pos);
  gr1_sasha->SetMarkerColor(kBlue);
  gr1_sasha->SetLineColor(kBlue);
  gr1_sasha->SetMarkerStyle(33);
  gr1_sasha->SetMarkerSize(1.4);
  gr1_sasha->Draw("P");

  TGraphErrors *gr1_pasha = new TGraphErrors(9, x_centrality, y, x_centrality_err, y_err);
  gr1_pasha->SetMarkerColor(kRed);
  gr1_pasha->SetLineColor(kRed);
  gr1_pasha->SetMarkerStyle(33);
  gr1_pasha->SetMarkerSize(1.4);
  gr1_pasha->Draw("P");

  TLegend *legC2_1 = new TLegend(0.5,0.70,0.89,.79);
  legC2_1->AddEntry(gr1_pasha,Form("%s %s, #sqrt{s_{NN}}= 19.6 GeV", PID, charge),"p");
  legC2_1->AddEntry(gr1_sasha,Form("%s %s, #sqrt{s_{NN}}= 39 GeV", PID, charge),"p");
  legC2_1->SetFillColor(kWhite);
  legC2_1->SetBorderSize(0);
  legC2_1->Draw();

  ctar->SaveAs(Form("/home/demanov/NIRS/Diplom/identified_hadrons/39_19gev/v%i_%s%s.png", harmonic ,PID,charge));
}


void flow_vs_pt( Double_t pt39GeV[], Double_t bin[] ,  const Char_t *PID,  const Char_t *charge, Int_t harmonic){

  makeplotstyle();

  TFile *f = new TFile("answer_flow.root","READ");

  Int_t cent[10]={80,70,60,50,40,30,20,10,5,0};
  
  Double_t x[col]={0.0};
  Double_t y[col]={0.0};
  Double_t y_err[col]={0.0};
  Double_t y_rat[col]={0.0};
  Double_t y_rat_err[col]={0.0};
  
  Double_t pt39GeV_err[col] = {0};
  
  Double_t x_coord_1[col], y_coord_1[col];
  Double_t x_coord_1err[col] = {0.0};
  Double_t y_coord_1err[col] = {0.0};
  
  TFile *fstyle = new TFile("style.root");
  TStyle *tsty = (TStyle *)fstyle->Get("PlottingInStyle");
  tsty->cd();

  gROOT->ForceStyle();
  gStyle->SetOptStat(0);

  char title1[800];
  sprintf(title1,"");

  float xmin1=0.0;
  float xmax1=5.0;
  float ymin1=0.00;
  float ymax1;
  float leg;
  Int_t start;

  if(harmonic==3){
    ymax1=0.15;
    leg=0.137;
    start=3;
  }
  if(harmonic==2){
    ymax1=0.35;
    leg=0.32;
    start=1;
  }

  Int_t nPID;

  if( strncmp(PID, "Proton",6)==0 ){
    if(strncmp(charge, "_pos",3)==0){
      nPID=4;
    }
    else{
      nPID=5;
    }
    if(strncmp(charge, "",1)==0){
      nPID=8;
    }
  }
  if( strncmp(PID, "Kaon",4)==0 ){
    if(strncmp(charge, "_pos",3)==0){
      nPID=2;
    }
    else{
      nPID=3;
    }
    if(strncmp(charge, "",1)==0){
      nPID=7;
    }
  }
  if( strncmp(PID, "Pion",4)==0 ){
    if(strncmp(charge, "_pos",3)==0){
      nPID=0;
    }
    else{
      nPID=1;
    }
    if(strncmp(charge, "",1)==0){
      nPID=6;
    }
  }


  TLegend *legC3 = new TLegend(0.7,0.5,0.89,.85);
  legC3->SetHeader(" Centrality:");
  legC3->SetFillColor(kWhite);
  legC3->SetBorderSize(0); 

  TH2F *hr2 = new TH2F("hr2",title1, 2,xmin1,xmax1,2,ymin1,ymax1);
  hr2->SetXTitle("p_{T} (GeV/c) ");
  hr2->SetYTitle(Form("V_{%i}        ",harmonic));

  TLatex* harmMark = new TLatex(0.2, leg, Form("Au+Au #sqrt{s_{NN}}= 39 GeV, v_{%i} , %s ", harmonic,particle[nPID]));
  harmMark->SetTextFont(47);
  harmMark->SetTextSize(21);

  TCanvas *c1_1 = new TCanvas("c1_1","Flow analysis results",150,10,800,650);
  hr2->Draw();
  harmMark->Draw("same");

  TProfile2D *p_v2_1 = (TProfile2D*)f -> Get(Form("v%i_pt_prof_n_2_%s%s", harmonic, PID, charge));
  for(Int_t j=9; j>start; j--){
    TH1D *h1 = (TH1D*) p_v2_1->ProfileY("prof1",j,j)->Rebin(col,"h1",bin);
    for(Int_t i=0; i<col; i++){
      y[i] = h1->GetBinContent(h1->FindBin(pt39GeV[i]));
      x[i] = h1->GetBinCenter(h1->FindBin(pt39GeV[i]));
      y_err[i] = h1->GetBinError(h1->FindBin(pt39GeV[i]));
    }
  
    TGraphErrors *gr1 = new TGraphErrors(col, x, y, pt39GeV_err, y_err);
    gr1->SetMarkerColor(color_2[j]);
    gr1->SetLineColor(color_2[j]);
    gr1->SetMarkerStyle(style_2[j]);
    gr1->SetMarkerSize(1.4);
    gr1->Draw("P");

    legC3->AddEntry(gr1,Form("%i-%i %%",cent[j], cent[j-1] ),"p");
  }

  legC3->Draw();

  c1_1 -> SaveAs(Form("/home/demanov/NIRS/Diplom/identified_hadrons/pict/centr/v%i_%s%s.png", harmonic, PID, charge ));


}

void flow_par_antipar( Double_t pt39GeV[], Double_t bin[] , const Char_t *charge, Int_t harmonic, Int_t bin1, Int_t bin2){

  makeplotstyle();

  TFile *f = new TFile("answer_flow.root","READ");

  Int_t cent[10]={80,70,60,50,40,30,20,10,5,0};
  
  Double_t x[col]={0.0};
  Double_t y[col]={0.0};
  Double_t y_err[col]={0.0};

  Double_t y_2[col]={0.0};
  Double_t y_2_err[col]={0.0};

  Double_t y_3[col]={0.0};
  Double_t y_3_err[col]={0.0};
  
  Double_t pt39GeV_err[col] = {0};
  
  TFile *fstyle = new TFile("style.root");
  TStyle *tsty = (TStyle *)fstyle->Get("PlottingInStyle");
  tsty->cd();

  gROOT->ForceStyle();
  gStyle->SetOptStat(0);

  char title1[800];
  sprintf(title1,"");

  float xmin1=0.0;
  float xmax1=5.0;
  float ymin1=0.00;
  float ymax1;
  float leg;

  if(harmonic==3){
    ymax1=0.15;
    leg=0.137;
  }
  if(harmonic==2){
    ymax1=0.25;
    leg=0.23;
  }


  TLegend *legC3 = new TLegend(0.7,0.5,0.89,.85);
  legC3->SetHeader(" Particle:");
  legC3->SetFillColor(kWhite);
  legC3->SetBorderSize(0); 

  TH2F *hr2 = new TH2F("hr2",title1, 2,xmin1,xmax1,2,ymin1,ymax1);
  hr2->SetXTitle("P_{T} (GeV/c) ");
  hr2->SetYTitle(Form("V_{%i}        ",harmonic));

  TLatex* harmMark = new TLatex(0.2, leg, Form("Au+Au #sqrt{s_{NN}}= 39 GeV, v_{%i}, %i-%i %% ", harmonic, cent[bin2], cent[bin1-1]));
  harmMark->SetTextFont(47);
  harmMark->SetTextSize(21);

  TCanvas *c1_1 = new TCanvas("c1_1","Flow analysis results",150,10,800,650);
  hr2->Draw();
  harmMark->Draw("same");

  TProfile2D *p_v2_pion = (TProfile2D*)f -> Get(Form("v%i_pt_prof_n_2_Pion%s", harmonic,  charge));
  TH1D *h1 = (TH1D*) p_v2_pion->ProfileY("prof1",bin1, bin2)->Rebin(col,"h1",bin);
  
  TProfile2D *p_v2_kaon = (TProfile2D*)f -> Get(Form("v%i_pt_prof_n_2_Kaon%s", harmonic,  charge));
  TH1D *h2 = (TH1D*) p_v2_kaon->ProfileY("prof1",bin1, bin2)->Rebin(col,"h2",bin);
  
  TProfile2D *p_v2_proton = (TProfile2D*)f -> Get(Form("v%i_pt_prof_n_2_Proton%s", harmonic,  charge));
  TH1D *h3 = (TH1D*) p_v2_proton->ProfileY("prof1",bin1, bin2)->Rebin(col,"h3",bin);
  for(Int_t i=0; i<col; i++){
    y[i] = h1->GetBinContent(h1->FindBin(pt39GeV[i]));
    x[i] = h1->GetBinCenter(h1->FindBin(pt39GeV[i]));
    y_err[i] = h1->GetBinError(h1->FindBin(pt39GeV[i]));

    y_2[i] = h2->GetBinContent(h2->FindBin(pt39GeV[i]));
    y_2_err[i] = h2->GetBinError(h2->FindBin(pt39GeV[i]));

    y_3[i] = h3->GetBinContent(h3->FindBin(pt39GeV[i]));
    y_3_err[i] = h3->GetBinError(h3->FindBin(pt39GeV[i]));
  }

  TGraphErrors *gr1 = new TGraphErrors(col, x, y, pt39GeV_err, y_err);
  gr1->SetMarkerColor(kRed);
  gr1->SetLineColor(kRed);
  gr1->SetMarkerStyle(22);
  gr1->SetMarkerSize(1.4);
  gr1->Draw("P");

  TGraphErrors *gr2 = new TGraphErrors(col, x, y_2, pt39GeV_err, y_2_err);
  gr2->SetMarkerColor(kBlue);
  gr2->SetLineColor(kBlue);
  gr2->SetMarkerStyle(8);
  gr2->SetMarkerSize(1.4);
  gr2->Draw("P");

  TGraphErrors *gr3 = new TGraphErrors(col, x, y_3, pt39GeV_err, y_3_err);
  gr3->SetMarkerColor(kBlack);
  gr3->SetLineColor(kBlack);
  gr3->SetMarkerStyle(33);
  gr3->SetMarkerSize(1.4);
  gr3->Draw("P");

  if(strncmp(charge, "_pos",3)==0){
    legC3->AddEntry(gr1,"#pi^{+}","p");
    legC3->AddEntry(gr2,"K^{+}","p");
    legC3->AddEntry(gr3,"p","p");
  }
  else{
    legC3->AddEntry(gr1,"#pi^{-}","p");
    legC3->AddEntry(gr2,"K^{-}","p");
    legC3->AddEntry(gr3,"#bar{p}","p");
  }
  

  legC3->Draw();

  c1_1 -> SaveAs(Form("/home/demanov/NIRS/Diplom/identified_hadrons/pict/centr/particle_v%i%s.png", harmonic, charge));


}



void Plot_ind(){

  makeplotstyle();
  //Res_all();
  //Res();

  TFile *fstyle = new TFile("style.root");
  TStyle *tsty = (TStyle *)fstyle->Get("PlottingInStyle");
  tsty->cd();

  gROOT->ForceStyle();
  gStyle->SetOptStat(0);

  //Double_t pt39GeV[14]  = {0.3, 0.5, 0.7, 0.9, 1.1, 1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,3.0};
  Double_t pt39GeV_err[14] = {0};

  Double_t v2_cent_39GeV[9] = {0.02089, 0.03102, 0.04407, 0.05658, 0.06318, 0.06436, 0.06122, 0.05585, 0.05284};
  Double_t v2_cent_39GeV_err[9] = {0.00004, 0.00003, 0.00002, 0.00003, 0.00003, 0.00005, 0.00009, 0.00018, 0.00036};
  
  Double_t v2_10_20[14]  = {0.0210,0.0393,0.0558,0.0706,0.0846,0.0975,0.1099,0.1205,0.1302,0.1371,0.1444,0.1487,0.1524,0.1592};
  Double_t v2_10_20_err[14] = {.0000,.0000,.0001,.0001,.0001,.0001,.0002,.0002,.0003,.0004,.0006,.0008,.0011,.0011};

  Double_t v2_20_30[14]  = {0.0270,0.0509,0.0725,0.0919,0.1102,0.1270,0.1422,0.1556,0.1674,0.1770,0.1840,0.1913,0.1948,0.2023};
  Double_t v2_20_30_err[14] = {.0001,.0001,.0001,.0002,.0002,.0003,.0004,.0006,.0008,.0010,.0014,.0020,.0027,.0029};

  Double_t v2_30_40[14]  = {0.0304,0.0576,0.0824,0.1047,0.1251,0.1436,0.1604,0.1751,0.1868,0.1969,0.2032,0.2099,0.2182, 0.2168};
  Double_t v2_30_40_err[14] = {.0001,.0001,.0001,.0001,.0002,.0002,.0003,.0004,.0005,.0007,.0009,.0013,.0017,.0018};

  Double_t v2_50_60[14]  = {0.0298,0.0571,0.0839,0.1074,0.1286,0.1469,0.1639,0.1749,0.1885,0.1994,0.1999,0.2116,0.2119, 0.2163};
  Double_t v2_50_60_err[14] = {.0001,.0002,.0002,.0003,.0005,.0006,.0008,.0011,.0015,.0020,.0027,.0037,.0050,.0052};

 // ----PID

// ---------------- Particle species: Proton ----------------
Double_t pt_bin_center_Pr_0_80[15] = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
Double_t v2_values_Pr_0_80[15] = {0.00683153,0.0158907,0.0315793,0.0516096,0.0727302,0.0925249,0.110851,0.126435,0.140123,0.150104,0.160199,0.165746,0.172611,0.179436,0.184446};
Double_t v2_stat_error_Pr_0_80[15] = {0.000520759,0.000195399,0.000163827,0.000164381,0.00018273,0.000216225,0.000267061,0.000342627,0.000453176,0.000557974,0.000736763,0.00100473,0.00138938,0.00155757,0.00303019};  
// ---------------- Particle species: antiProton --------------
Double_t pt_bin_center_Ant_0_80[15] = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
Double_t v2_values_Ant_0_80[15] = {0.00494108,0.0114479,0.0261197,0.0449316,0.0663671,0.0860504,0.103827,0.120229,0.132341,0.142772,0.151479,0.165265,0.16601,0.17546,0.178903};
Double_t v2_stat_error_Ant_0_80[15] = {0.000897342,0.000371594,0.000306912,0.000307998,0.000345265,0.000415456,0.000522919,0.000689985,0.000821512,0.00110148,0.00148064,0.00206483,0.00291789,0.00336361,0.00684634};
// ---------------- Particle species: PiM ----------------
Double_t pt_bin_center_PiM_0_80[15] = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
Double_t v2_values_PiM_0_80[15] = {0.0223896,0.042742,0.0624241,0.0798173,0.0944707,0.106976,0.11663,0.125184,0.131548,0.135949,0.138445,0.14286,0.148789,0.143203,0.151072};
Double_t v2_stat_error_PiM_0_80[15] = {4.02884e-05,4.85174e-05,6.65846e-05,9.43872e-05,0.000134239,0.000189959,0.000271414,0.000385102,0.00053152,0.000717481,0.000968061,0.00129079,0.00172756,0.00185626,0.00332027};
// ---------------- Particle species: PiP ----------------
Double_t pt_bin_center_PiP_0_80[15] = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
Double_t v2_values_PiP_0_80[15] = {0.0208012,0.041395,0.0611387,0.0783936,0.0932306,0.105414,0.115887,0.123129,0.129545,0.134079,0.136608,0.141293,0.145453,0.144127,0.147837};
Double_t v2_stat_error_PiP_0_80[15] = {4.06899e-05,4.87531e-05,6.64991e-05,9.38308e-05,0.000132824,0.000187112,0.000265895,0.000375892,0.000515683,0.000692813,0.000918648,0.00122249,0.00162044,0.00172633,0.00303448};
// ---------------- Particle species: KM ----------------
Double_t pt_bin_center_KM_0_80[15] = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.265,2.475,2.685,2.955,3.345};
Double_t v2_values_KM_0_80[15] = {0.0129623,0.02816,0.0468495,0.0647472,0.0797767,0.0927023,0.102755,0.109907,0.117535,0.120576,0.121677,0.124428,0.126282,0.144784,0.139207};
Double_t v2_stat_error_KM_0_80[15] = {0.000275854,0.000175113,0.00017688,0.000208595,0.000262253,0.000341829,0.000461195,0.000635014,0.000904939,0.00131148,0.00196627,0.00289826,0.00425974,0.0049056,0.0110844};
// ---------------- Particle species: KP ----------------
Double_t pt_bin_center_KP_0_80[15] = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.265,2.475,2.685,2.955,3.345};
Double_t v2_values_KP_0_80[15] = {0.0131131,0.0280332,0.0473805,0.0650401,0.0797504,0.0919515,0.101452,0.108697,0.114705,0.120665,0.121731,0.128768,0.123222,0.127912,0.134311};
Double_t v2_stat_error_KP_0_80[15] = {0.000244437,0.000153924,0.000154119,0.000180181,0.000224963,0.000290565,0.000386755,0.000523833,0.000719747,0.0010333,0.00151358,0.002247,0.00334946,0.00396006,0.0063682};


// ---------------- Particle species: PiM ----------------
Double_t pt_bin_center_PiM_0_10[13] = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.355,2.805,3.465};
Double_t v2_values_PiM_0_10[13] = {0.0147436,0.0249519,0.0351931,0.0444924,0.0515255,0.0572386,0.0624405,0.0665682,0.0680159,0.0696394,0.0726135,0.0706268,0.0619339};
Double_t v2_stat_error_PiM_0_10[13] = {7.87931e-05,9.36376e-05,0.000127759,0.000180173,0.000255988,0.000362511,0.000519879,0.000740854,0.00102086,0.00138664,0.00150644,0.00248854,0.00578472};
// ---------------- Particle species: PiP ----------------
Double_t pt_bin_center_PiP_0_10[13] = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.355,2.805,3.465};
Double_t v2_values_PiP_0_10[13] = {0.0118138,0.0229824,0.0337497,0.0432688,0.0511433,0.0577776,0.0629835,0.0666356,0.0688837,0.0707009,0.0758873,0.0793376,0.0710317};
Double_t v2_stat_error_PiP_0_10[13] = {7.97649e-05,9.40754e-05,0.000127369,0.000178608,0.000252476,0.000355599,0.00050783,0.000718169,0.000987558,0.00132172,0.00141911,0.0022635,0.00528366};;
// ---------------- Particle species: KM ----------------
Double_t pt_bin_center_KM_0_10[13] = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.355,2.805,3.435};
Double_t v2_values_KM_0_10[13] = {0.00757911,0.0157954,0.026121,0.0363908,0.0452336,0.0522223,0.058466,0.0605987,0.0681485,0.0680222,0.0658601,0.0699299,0.10366};
Double_t v2_stat_error_KM_0_10[13] = {0.000547625,0.000341637,0.000339806,0.000396405,0.000494154,0.000642514,0.000865204,0.00118888,0.00169966,0.00245065,0.00300416,0.00624675,0.0188689};
// ---------------- Particle species: KP ----------------
Double_t pt_bin_center_KP_0_10[13] = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.355,2.805,3.435};
Double_t v2_values_KP_0_10[13] = {0.00555827,0.0135945,0.0257326,0.0364142,0.045515,0.0526201,0.0574833,0.0609207,0.0644966,0.0688554,0.0702958,0.0619705,0.0345815};
Double_t v2_stat_error_KP_0_10[13] = {0.000483909,0.000298923,0.00029486,0.000341234,0.000422957,0.000545761,0.00072527,0.000981056,0.00134703,0.00191111,0.00236173,0.00523082,0.0185284};
// ---------------- Particle species: Proton ----------------
Double_t pt_bin_center_Pr_0_10[15] = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
Double_t v2_values_Pr_0_10[15] = {0.000344971,0.0042966,0.0113989,0.0231444,0.0362371,0.0487195,0.0597448,0.0688017,0.0768114,0.0834786,0.0895703,0.0933591,0.0928367,0.0985281,0.0996354};
Double_t v2_stat_error_Pr_0_10[15] = {0.00107229,0.000396078,0.000324797,0.000318868,0.000347991,0.00040676,0.00049864,0.000637553,0.000842688,0.00103321,0.00136462,0.00186161,0.0025729,0.00289198,0.00563872};
// ---------------- Particle species: antiProton ----------------
Double_t pt_bin_center_Ant_0_10[15] = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
Double_t v2_values_Ant_0_10[15] = {0.00104822,0.00044334,0.00784047,0.0174989,0.0296773,0.042122,0.0521295,0.0608011,0.0658994,0.0712058,0.0747417,0.0849594,0.0807271,0.0906798,0.0848712};
Double_t v2_stat_error_Ant_0_10[15] = {0.00192614,0.000785758,0.000629482,0.000613868,0.000671617,0.000795191,0.000988842,0.00129467,0.00151703,0.00203264,0.00271891,0.00378032,0.00532838,0.00609857,0.0123108};




// ---------------- Particle species: Proton ----------------
Double_t pt_bin_center_Pr_10_40[15] = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
Double_t v2_values_Pr_10_40[15] = {0.00625568,0.0163721,0.0351739,0.0584888,0.0837165,0.107311,0.129351,0.148297,0.164633,0.176239,0.18757,0.195966,0.205624,0.212297,0.216848};
Double_t v2_stat_error_Pr_10_40[15] = {0.000544837,0.000204691,0.000172133,0.000173382,0.000193466,0.000229269,0.000283342,0.000363145,0.000479106,0.000593701,0.000783826,0.00106789,0.00147689,0.00165122,0.00321044};
// ---------------- Particle species: antiProton ----------------
Double_t pt_bin_center_Ant_10_40[15] = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
Double_t v2_values_Ant_10_40[15] = {0.00367436,0.0106725,0.0266682,0.0486995,0.0735285,0.097254,0.118446,0.1377,0.154429,0.16567,0.177735,0.194259,0.197949,0.205333,0.221208};
Double_t v2_stat_error_Ant_10_40[15] = {0.000925613,0.00038232,0.000316434,0.000318902,0.000358874,0.000432563,0.000545001,0.000718362,0.000865585,0.00115941,0.00155836,0.00217497,0.00308067,0.00356517,0.00727714};
// ---------------- Particle species: PiM ----------------
Double_t pt_bin_center_PiM_10_40[14] = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.295,2.475,2.835,3.465};
Double_t v2_values_PiM_10_40[14] = {0.024917,0.0485882,0.0711891,0.0911594,0.108139,0.122315,0.133424,0.142324,0.149065,0.15394,0.156085,0.159299,0.163019,0.16689};
Double_t v2_stat_error_PiM_10_40[14] = {4.20321e-05,5.07775e-05,6.97073e-05,9.88055e-05,0.000140512,0.000198458,0.000282845,0.000400946,0.000554477,0.000747931,0.00100679,0.0013418,0.00132684,0.00296483};
// ---------------- Particle species: PiP ----------------
Double_t pt_bin_center_PiP_10_40[14] = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.295,2.475,2.835,3.465};
Double_t v2_values_PiP_10_40[14] = {0.0236042,0.0474268,0.0701121,0.0900108,0.106989,0.120903,0.132197,0.140385,0.146647,0.152065,0.151277,0.157466,0.158328,0.155246};
Double_t v2_stat_error_PiP_10_40[14] = {4.24396e-05,5.10583e-05,6.97025e-05,9.83866e-05,0.000139269,0.000195935,0.000277684,0.000392379,0.000541797,0.00072417,0.000959308,0.00127277,0.00122297,0.00272259};
// ---------------- Particle species: KM ---------------
Double_t pt_bin_center_KM_10_40[14] = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.265,2.475,2.835,3.465};
Double_t v2_values_KM_10_40[14] = {0.0132178,0.0312904,0.0537999,0.0745947,0.0931274,0.108946,0.121914,0.131265,0.138464,0.142195,0.146256,0.150271,0.154251,0.144235};
Double_t v2_stat_error_KM_10_40[14] = {0.000290698,0.000185182,0.000187752,0.000222014,0.000279484,0.000362775,0.000485722,0.000669697,0.000946231,0.00137214,0.00202761,0.00303479,0.00317906,0.00722992};
// ---------------- Particle species: KP ----------------
Double_t pt_bin_center_KP_10_40[14] = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.265,2.475,2.805,3.435};
Double_t v2_values_KP_10_40[14] = {0.0143249,0.0324692,0.054843,0.075677,0.0935934,0.108188,0.120058,0.129584,0.137597,0.142903,0.145207,0.147236,0.143582,0.124331};
Double_t v2_stat_error_KP_10_40[14] = {0.000257904,0.00016312,0.000164087,0.000192338,0.000240327,0.000309352,0.00040976,0.000551598,0.00076111,0.00108504,0.00158232,0.00233087,0.00256974,0.00503966};



// ---------------- Particle species: Proton ----------------
Double_t pt_bin_center_Pr_40_80[15] = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
Double_t v2_values_Pr_40_80[15] = {0.0173274,0.0329489,0.0579302,0.0899572,0.122273,0.153005,0.181775,0.205587,0.22685,0.241718,0.259529,0.254186,0.274926,0.281757,0.299217};
Double_t v2_stat_error_Pr_40_80[15] = {0.00163998,0.000645123,0.000578342,0.000621392,0.000734358,0.000911515,0.00116596,0.00153552,0.00207608,0.00253886,0.00335845,0.00459636,0.00633159,0.00709883,0.013688};
// ---------------- Particle species: antiProton ----------------
Double_t pt_bin_center_Ant_40_80[12] = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.355,2.805};
Double_t v2_values_Ant_40_80[12] = {0.0110292,0.0256727,0.0503999,0.0802053,0.115991,0.143638,0.174016,0.205972,0.22776,0.253882,0.27293,0.298666};
Double_t v2_stat_error_Ant_40_80[12] = {0.00254538,0.00110941,0.000989204,0.00107262,0.00129246,0.00164622,0.00217035,0.00298469,0.00362504,0.00493413,0.00557032,0.0103886};
// ---------------- Particle species: PiM ----------------
Double_t pt_bin_center_PiM_40_80[12] = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.355,2.835};
Double_t v2_values_PiM_40_80[12] = {0.0289077,0.0593903,0.0895478,0.116196,0.140427,0.159371,0.172953,0.185262,0.197665,0.205783,0.21454,0.225593};
Double_t v2_stat_error_PiM_40_80[12] = {0.000142567,0.000177444,0.0002491,0.000359443,0.000519995,0.000745305,0.00106672,0.00149978,0.00204763,0.0027489,0.00295397,0.00481477};
// ---------------- Particle species: PiP ----------------
Double_t pt_bin_center_PiP_40_80[12] = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.355,2.835};
Double_t v2_values_PiP_40_80[12] = {0.028888,0.0589122,0.0884838,0.114564,0.138449,0.156339,0.173457,0.181754,0.19405,0.19931,0.210363,0.227665};
Double_t v2_stat_error_PiP_40_80[12] = {0.000143401,0.000178332,0.000249537,0.000358916,0.000517534,0.000736543,0.00105173,0.00148376,0.00201896,0.00269214,0.00284612,0.0045173};
// ---------------- Particle species: KM ----------------
Double_t pt_bin_center_KM_40_80[12] = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.355,2.835};
Double_t v2_values_KM_40_80[12] = {0.0225864,0.0437617,0.0699159,0.0982409,0.119146,0.139389,0.153179,0.167905,0.181927,0.188784,0.178602,0.210977};
Double_t v2_stat_error_KM_40_80[12] = {0.000963893,0.000642692,0.000679603,0.000830741,0.00107209,0.00142409,0.00193787,0.0026898,0.00384267,0.00561045,0.00660844,0.0116723};
// ---------------- Particle species: KP ----------------
Double_t pt_bin_center_KP_40_80[12] = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.355,2.835};
Double_t v2_values_KP_40_80[12] = {0.0239543,0.044623,0.0718769,0.0978395,0.117474,0.136554,0.15272,0.165516,0.17048,0.181233,0.184436,0.167675};
Double_t v2_stat_error_KP_40_80[12] = {0.0008657,0.000573558,0.000600633,0.000725502,0.000927351,0.00121505,0.00162533,0.00219898,0.00304011,0.00429892,0.00495433,0.00975087};


Double_t pt_bin_center_Pr_0_80_bin[16] = {0.315-0.1,0.495-0.1,0.705-0.1,0.885-0.1,1.095-0.1,1.275-0.1,1.485-0.1,1.695-0.1,1.875-0.1,2.085-0.1,2.295-0.1,2.475-0.1,2.685-0.1,2.955-0.1,3.345-0.1, 3.5};
Double_t pt_bin_center_Ant_0_80_bin[16] = {0.315-0.1,0.495-0.1,0.705-0.1,0.885-0.1,1.095-0.1,1.275-0.1,1.485-0.1,1.695-0.1,1.875-0.1,2.085-0.1,2.295-0.1,2.475-0.1,2.685-0.1,2.955-0.1,3.345-0.1, 3.5};
Double_t pt_bin_center_PiM_0_80_bin[16] = {0.285-0.1,0.495-0.1,0.675-0.1,0.885-0.1,1.065-0.1,1.275-0.1,1.485-0.1,1.665-0.1,1.875-0.1,2.085-0.1,2.295-0.1,2.475-0.1,2.685-0.1,2.955-0.1,3.345-0.1, 3.5};
Double_t pt_bin_center_PiP_0_80_bin[16] = {0.285-0.1,0.495-0.1,0.675-0.1,0.885-0.1,1.065-0.1,1.275-0.1,1.485-0.1,1.665-0.1,1.875-0.1,2.085-0.1,2.295-0.1,2.475-0.1,2.685-0.1,2.955-0.1,3.345-0.1, 3.5};
Double_t pt_bin_center_KM_0_80_bin[16] = {0.315-0.1,0.495-0.1,0.705-0.1,0.885-0.1,1.095-0.1,1.275-0.1,1.485-0.1,1.695-0.1,1.875-0.1,2.085-0.1,2.265-0.1,2.475-0.1,2.685-0.1,2.955-0.1,3.345-0.1, 3.5};
Double_t pt_bin_center_KP_0_80_bin[16] = {0.315-0.1,0.495-0.1,0.705-0.1,0.885-0.1,1.095-0.1,1.275-0.1,1.485-0.1,1.695-0.1,1.875-0.1,2.085-0.1,2.265-0.1,2.475-0.1,2.685-0.1,2.955-0.1,3.345-0.1, 3.5};

Double_t pt_bin_center_PiM_0_10_bin[14] = {0.285-0.1,0.495-0.1,0.675-0.1,0.885-0.1,1.065-0.1,1.275-0.1,1.485-0.1,1.665-0.1,1.875-0.1,2.085-0.1,2.355-0.1,2.805-0.1,3.465-0.1, 3.5};
Double_t pt_bin_center_PiP_0_10_bin[14] = {0.285-0.1,0.495-0.1,0.675-0.1,0.885-0.1,1.065-0.1,1.275-0.1,1.485-0.1,1.665-0.1,1.875-0.1,2.085-0.1,2.355-0.1,2.805-0.1,3.465-0.1, 3.5};
Double_t pt_bin_center_KM_0_10_bin[14] = {0.315-0.1,0.495-0.1,0.705-0.1,0.885-0.1,1.095-0.1,1.275-0.1,1.485-0.1,1.695-0.1,1.875-0.1,2.085-0.1,2.355-0.1,2.805-0.1,3.435-0.1, 3.5};
Double_t pt_bin_center_KP_0_10_bin[14] = {0.315-0.1,0.495-0.1,0.705-0.1,0.885-0.1,1.095-0.1,1.275-0.1,1.485-0.1,1.695-0.1,1.875-0.1,2.085-0.1,2.355-0.1,2.805-0.1,3.435-0.1, 3.5};
Double_t pt_bin_center_Pr_0_10_bin[16] = {0.315-0.1,0.495-0.1,0.705-0.1,0.885-0.1,1.095-0.1,1.275-0.1,1.485-0.1,1.695-0.1,1.875-0.1,2.085-0.1,2.295-0.1,2.475-0.1,2.685-0.1,2.955-0.1,3.345-0.1, 3.5};
Double_t pt_bin_center_Ant_0_10_bin[16] = {0.315-0.1,0.495-0.1,0.705-0.1,0.885-0.1,1.095-0.1,1.275-0.1,1.485-0.1,1.695-0.1,1.875-0.1,2.085-0.1,2.295-0.1,2.475-0.1,2.685-0.1,2.955-0.1,3.345-0.1, 3.5};


Double_t pt_bin_center_Pr_10_40_bin[16] = {0.315-0.1,0.495-0.1,0.705-0.1,0.885-0.1,1.095-0.1,1.275-0.1,1.485-0.1,1.695-0.1,1.875-0.1,2.085-0.1,2.295-0.1,2.475-0.1,2.685-0.1,2.955-0.1,3.345-0.1, 3.5};
Double_t pt_bin_center_Ant_10_40_bin[16] = {0.315-0.1,0.495-0.1,0.705-0.1,0.885-0.1,1.095-0.1,1.275-0.1,1.485-0.1,1.695-0.1,1.875-0.1,2.085-0.1,2.295-0.1,2.475-0.1,2.685-0.1,2.955-0.1,3.345-0.1, 3.5};
Double_t pt_bin_center_PiM_10_40_bin[15] = {0.285-0.1,0.495-0.1,0.675-0.1,0.885-0.1,1.065-0.1,1.275-0.1,1.485-0.1,1.665-0.1,1.875-0.1,2.085-0.1,2.295-0.1,2.475-0.1,2.835-0.1,3.465-0.1, 3.5};
Double_t pt_bin_center_PiP_10_40_bin[15] = {0.285-0.1,0.495-0.1,0.675-0.1,0.885-0.1,1.065-0.1,1.275-0.1,1.485-0.1,1.665-0.1,1.875-0.1,2.085-0.1,2.295-0.1,2.475-0.1,2.835-0.1,3.465-0.1, 3.5};
Double_t pt_bin_center_KM_10_40_bin[15] = {0.315-0.1,0.495-0.1,0.705-0.1,0.885-0.1,1.095-0.1,1.275-0.1,1.485-0.1,1.695-0.1,1.875-0.1,2.085-0.1,2.265-0.1,2.475-0.1,2.835-0.1,3.465-0.1, 3.5};
Double_t pt_bin_center_KP_10_40_bin[15] = {0.315-0.1,0.495-0.1,0.705-0.1,0.885-0.1,1.095-0.1,1.275-0.1,1.485-0.1,1.695-0.1,1.875-0.1,2.085-0.1,2.265-0.1,2.475-0.1,2.805-0.1,3.435-0.1, 3.5};

Double_t pt_bin_center_Pr_40_80_bin[16] = {0.315-0.1,0.495-0.1,0.705-0.1,0.885-0.1,1.095-0.1,1.275-0.1,1.485-0.1,1.695-0.1,1.875-0.1,2.085-0.1,2.295-0.1,2.475-0.1,2.685-0.1,2.955-0.1,3.345-0.1, 3.0};
Double_t pt_bin_center_Ant_40_80_bin[13] = {0.315-0.1,0.495-0.1,0.705-0.1,0.885-0.1,1.095-0.1,1.275-0.1,1.485-0.1,1.695-0.1,1.875-0.1,2.085-0.1,2.355-0.1,2.805-0.1, 3.0};
Double_t pt_bin_center_PiM_40_80_bin[13] = {0.285-0.1,0.495-0.1,0.675-0.1,0.885-0.1,1.065-0.1,1.275-0.1,1.485-0.1,1.665-0.1,1.875-0.1,2.085-0.1,2.355-0.1,2.835-0.1, 3.0};
Double_t pt_bin_center_PiP_40_80_bin[13] = {0.285-0.1,0.495-0.1,0.675-0.1,0.885-0.1,1.065-0.1,1.275-0.1,1.485-0.1,1.665-0.1,1.875-0.1,2.085-0.1,2.355-0.1,2.835-0.1, 3.0};
Double_t pt_bin_center_KM_40_80_bin[13] = {0.315-0.1,0.495-0.1,0.705-0.1,0.885-0.1,1.095-0.1,1.275-0.1,1.485-0.1,1.695-0.1,1.875-0.1,2.085-0.1,2.355-0.1,2.835-0.1, 3.0};
Double_t pt_bin_center_KP_40_80_bin[13] = {0.315-0.1,0.495-0.1,0.705-0.1,0.885-0.1,1.095-0.1,1.275-0.1,1.485-0.1,1.695-0.1,1.875-0.1,2.085-0.1,2.355-0.1,2.835-0.1, 3.0};

/*
  flow_vs_pt( pt_bin_center_PiM_0_80, pt_bin_center_PiM_0_80_bin, "Proton","", 2 );
  flow_vs_pt( pt_bin_center_PiM_0_80, pt_bin_center_PiM_0_80_bin, "Proton","_pos", 2 );
  flow_vs_pt( pt_bin_center_PiM_0_80, pt_bin_center_PiM_0_80_bin, "Proton","_neg", 2 );
  flow_vs_pt( pt_bin_center_PiM_0_80, pt_bin_center_PiM_0_80_bin, "Proton","", 3 );
  flow_vs_pt( pt_bin_center_PiM_0_80, pt_bin_center_PiM_0_80_bin, "Proton","_pos", 3 );
  flow_vs_pt( pt_bin_center_PiM_0_80, pt_bin_center_PiM_0_80_bin, "Proton","_neg", 3 );

  flow_vs_pt( pt_bin_center_PiM_0_80, pt_bin_center_PiM_0_80_bin, "Pion","", 2 );
  flow_vs_pt( pt_bin_center_PiM_0_80, pt_bin_center_PiM_0_80_bin, "Pion","_pos", 2 );
  flow_vs_pt( pt_bin_center_PiM_0_80, pt_bin_center_PiM_0_80_bin, "Pion","_neg", 2 );
  flow_vs_pt( pt_bin_center_PiM_0_80, pt_bin_center_PiM_0_80_bin, "Pion","", 3 );
  flow_vs_pt( pt_bin_center_PiM_0_80, pt_bin_center_PiM_0_80_bin, "Pion","_pos", 3 );
  flow_vs_pt( pt_bin_center_PiM_0_80, pt_bin_center_PiM_0_80_bin, "Pion","_neg", 3 );

  flow_vs_pt( pt_bin_center_PiM_0_80, pt_bin_center_PiM_0_80_bin, "Kaon","", 2 );
  flow_vs_pt( pt_bin_center_PiM_0_80, pt_bin_center_PiM_0_80_bin, "Kaon","_pos", 2 );
  flow_vs_pt( pt_bin_center_PiM_0_80, pt_bin_center_PiM_0_80_bin, "Kaon","_neg", 2 );
  flow_vs_pt( pt_bin_center_PiM_0_80, pt_bin_center_PiM_0_80_bin, "Kaon","", 3 );
  flow_vs_pt( pt_bin_center_PiM_0_80, pt_bin_center_PiM_0_80_bin, "Kaon","_pos", 3 );
  flow_vs_pt( pt_bin_center_PiM_0_80, pt_bin_center_PiM_0_80_bin, "Kaon","_neg", 3 );
*/
  flow_par_antipar( pt_bin_center_PiM_0_80, pt_bin_center_PiM_0_80_bin,"_pos", 2, 3,9 );
  flow_par_antipar( pt_bin_center_PiM_0_80, pt_bin_center_PiM_0_80_bin,"_neg", 2, 3,9 );
  flow_par_antipar( pt_bin_center_PiM_0_80, pt_bin_center_PiM_0_80_bin,"_pos", 3, 3,9 );
  flow_par_antipar( pt_bin_center_PiM_0_80, pt_bin_center_PiM_0_80_bin,"_neg", 3, 3,9 );

 // flow_vs_pt( pt_bin_center_PiM_0_80, pt_bin_center_PiM_0_80_bin, "Pion","_pos", 2 );

//----negative
  if(col==12){
  v2_pt( v2_values_PiM_40_80, v2_stat_error_PiM_40_80, pt_bin_center_PiM_40_80, pt_bin_center_PiM_40_80_bin , 1, 4, "Pion", "_neg");
  v3_pt(  pt_bin_center_PiM_40_80, pt_bin_center_PiM_40_80_bin , 1, 4, "Pion", "_neg");

  v2_pt( v2_values_KM_40_80, v2_stat_error_KM_40_80, pt_bin_center_KM_40_80, pt_bin_center_KM_40_80_bin , 1, 4, "Kaon", "_neg");
  v3_pt(  pt_bin_center_KM_40_80, pt_bin_center_KM_40_80_bin , 1, 4, "Kaon", "_neg");

  v2_pt( v2_values_Ant_40_80, v2_stat_error_Ant_40_80, pt_bin_center_Ant_40_80, pt_bin_center_Ant_40_80_bin , 1, 4, "Proton", "_neg");
  v3_pt(  pt_bin_center_Ant_40_80, pt_bin_center_Ant_40_80_bin , 1, 4, "Proton", "_neg");

  v2_pt( v2_values_PiP_40_80, v2_stat_error_PiP_40_80, pt_bin_center_PiP_40_80, pt_bin_center_PiP_40_80_bin , 1, 4, "Pion", "_pos");
  v3_pt(  pt_bin_center_PiP_40_80, pt_bin_center_PiP_40_80_bin , 1, 4, "Pion", "_pos");

  v2_pt( v2_values_KP_40_80, v2_stat_error_KP_40_80, pt_bin_center_KP_40_80, pt_bin_center_KP_40_80_bin , 1, 4, "Kaon", "_pos");
  v3_pt(  pt_bin_center_KP_40_80, pt_bin_center_KP_40_80_bin , 1, 4, "Kaon", "_pos");

  }

  if(col==13){
  v2_pt( v2_values_KM_0_10, v2_stat_error_KM_0_10, pt_bin_center_KM_0_10, pt_bin_center_KM_0_10_bin , 8, 9, "Kaon", "_neg");
  v3_pt(  pt_bin_center_KM_0_10, pt_bin_center_KM_0_10_bin , 8, 9, "Kaon", "_neg");
  v2_pt( v2_values_PiM_0_10, v2_stat_error_PiM_0_10, pt_bin_center_PiM_0_10, pt_bin_center_PiM_0_10_bin , 8, 9, "Pion", "_neg");
  v3_pt(  pt_bin_center_PiM_0_10, pt_bin_center_PiM_0_10_bin , 8, 9, "Pion", "_neg");

  v2_pt( v2_values_KP_0_10, v2_stat_error_KP_0_10, pt_bin_center_KP_0_10, pt_bin_center_KP_0_10_bin , 8, 9, "Kaon", "_pos");
  v3_pt(  pt_bin_center_KP_0_10, pt_bin_center_KP_0_10_bin , 8, 9, "Kaon", "_pos");
  v2_pt( v2_values_PiP_0_10, v2_stat_error_PiP_0_10, pt_bin_center_PiP_0_10, pt_bin_center_PiP_0_10_bin , 8, 9, "Pion", "_pos");
  v3_pt(  pt_bin_center_PiP_0_10, pt_bin_center_PiP_0_10_bin , 8, 9, "Pion", "_pos");
  }

  if(col==14){
  
  v2_pt( v2_values_PiM_10_40, v2_stat_error_PiM_10_40, pt_bin_center_PiM_10_40, pt_bin_center_PiM_10_40_bin , 5, 7, "Pion", "_neg");

  v3_pt(  pt_bin_center_PiM_10_40, pt_bin_center_PiM_10_40_bin , 5, 7, "Pion", "_neg");
    
  v2_pt( v2_values_KM_10_40, v2_stat_error_KM_10_40, pt_bin_center_KM_10_40, pt_bin_center_KM_10_40_bin , 5, 7, "Kaon", "_neg");
  v3_pt(  pt_bin_center_KM_10_40, pt_bin_center_KM_10_40_bin , 5, 7, "Kaon", "_neg");

  v2_pt( v2_values_PiP_10_40, v2_stat_error_PiP_10_40, pt_bin_center_PiP_10_40, pt_bin_center_PiP_10_40_bin , 5, 7, "Pion", "_pos");

  v3_pt(  pt_bin_center_PiP_10_40, pt_bin_center_PiP_10_40_bin , 5, 7, "Pion", "_pos");
    
  v2_pt( v2_values_KP_10_40, v2_stat_error_KP_10_40, pt_bin_center_KP_10_40, pt_bin_center_KP_10_40_bin , 5, 7, "Kaon", "_pos");
  v3_pt(  pt_bin_center_KP_10_40, pt_bin_center_KP_10_40_bin , 5, 7, "Kaon", "_pos");

  }

  if(col==155){
  v2_pt( v2_values_PiM_0_80, v2_stat_error_PiM_0_80, pt_bin_center_PiM_0_80, pt_bin_center_PiM_0_80_bin, 1, 9, "Pion", "_neg");
  v3_pt(  pt_bin_center_PiM_0_80, pt_bin_center_PiM_0_80_bin, 1, 9, "Pion", "_neg");
  
  v2_pt( v2_values_KM_0_80, v2_stat_error_KM_0_80, pt_bin_center_KM_0_80, pt_bin_center_KM_0_80_bin, 1, 9, "Kaon", "_neg");
  v3_pt(  pt_bin_center_KM_0_80, pt_bin_center_KM_0_80_bin, 1, 9, "Kaon", "_neg");
    

  v2_pt( v2_values_Ant_0_80, v2_stat_error_Ant_0_80, pt_bin_center_Ant_0_80, pt_bin_center_Ant_0_80_bin, 1, 9, "Proton", "_neg");
  v2_pt( v2_values_Ant_0_10, v2_stat_error_Ant_0_10, pt_bin_center_Ant_0_10, pt_bin_center_Ant_0_10_bin , 8, 9, "Proton", "_neg");
  v2_pt( v2_values_Ant_10_40, v2_stat_error_Ant_10_40, pt_bin_center_Ant_10_40, pt_bin_center_Ant_10_40_bin , 5, 7, "Proton", "_neg");
  //v2_pt( v2_values_Ant_40_80, v2_stat_error_Ant_40_80, pt_bin_center_Ant_40_80, pt_bin_center_Ant_40_80_bin , 1, 4, "Proton", "_neg");

  v3_pt(  pt_bin_center_Ant_0_80, pt_bin_center_Ant_0_80_bin, 1, 9, "Proton", "_neg");
  v3_pt(  pt_bin_center_Ant_0_10, pt_bin_center_Ant_0_10_bin , 8, 9, "Proton", "_neg");
  v3_pt(  pt_bin_center_Ant_10_40, pt_bin_center_Ant_10_40_bin , 5, 7, "Proton", "_neg");
  //v3_pt(  pt_bin_center_Ant_40_80, pt_bin_center_Ant_40_80_bin , 1, 4, "Proton", "_neg");

  v2_pt( v2_values_PiP_0_80, v2_stat_error_PiP_0_80, pt_bin_center_PiP_0_80, pt_bin_center_PiP_0_80_bin, 1, 9, "Pion", "_pos");
  v3_pt(  pt_bin_center_PiP_0_80, pt_bin_center_PiP_0_80_bin, 1, 9, "Pion", "_pos");
  
  v2_pt( v2_values_KP_0_80, v2_stat_error_KP_0_80, pt_bin_center_KP_0_80, pt_bin_center_KP_0_80_bin, 1, 9, "Kaon", "_pos");
  v3_pt(  pt_bin_center_KP_0_80, pt_bin_center_KP_0_80_bin, 1, 9, "Kaon", "_pos");
    

  v2_pt( v2_values_Pr_0_80, v2_stat_error_Pr_0_80, pt_bin_center_Pr_0_80, pt_bin_center_Pr_0_80_bin, 1, 9, "Proton", "_pos");
  v2_pt( v2_values_Pr_0_10, v2_stat_error_Pr_0_10, pt_bin_center_Pr_0_10, pt_bin_center_Pr_0_10_bin , 8, 9, "Proton", "_pos");
  v2_pt( v2_values_Pr_10_40, v2_stat_error_Pr_10_40, pt_bin_center_Pr_10_40, pt_bin_center_Pr_10_40_bin , 5, 7, "Proton", "_pos");
  v2_pt( v2_values_Pr_40_80, v2_stat_error_Pr_40_80, pt_bin_center_Pr_40_80, pt_bin_center_Pr_40_80_bin , 1, 4, "Proton", "_pos");

  v3_pt(  pt_bin_center_Pr_0_80, pt_bin_center_Pr_0_80_bin, 1, 9, "Proton", "_pos");
  v3_pt(  pt_bin_center_Pr_0_10, pt_bin_center_Pr_0_10_bin , 8, 9, "Proton", "_pos");
  v3_pt(  pt_bin_center_Pr_10_40, pt_bin_center_Pr_10_40_bin , 5, 7, "Proton", "_pos");
  v3_pt(  pt_bin_center_Pr_40_80, pt_bin_center_Pr_40_80_bin , 1, 4, "Proton", "_pos");

  }
  //KM
  //Ant
  //PiP
 if(col==16){
  v_cent(2,"Pion","_neg");
  v_cent(2,"Kaon","_neg");
  //v_cent(2,"Proton","_neg");
  v_cent(3,"Pion","_neg");
  v_cent(3,"Kaon","_neg");
  //v_cent(3,"Proton","_neg");

  v_cent(2,"Pion","_pos");
  v_cent(2,"Kaon","_pos");
  //v_cent(2,"Proton","_pos");
  v_cent(3,"Pion","_pos");
  v_cent(3,"Kaon","_pos");
  //v_cent(3,"Proton","_pos");
  }

  //Res();
  //v_cent_gev(2,"Pion","neg");
  //v_cent_gev(2,"Kaon","neg");
  //v_cent_gev(2,"Proton","neg");
  
  //v_cent_gev(3,"Pion","neg");
  //v_cent_gev(3,"Kaon","neg");
  //v_cent_gev(3,"Proton","neg");

  //v_cent_gev(2,"Pion","pos");
  //v_cent_gev(2,"Kaon","pos");
  //v_cent_gev(2,"Proton","pos");
  
  //v_cent_gev(3,"Pion","pos");
  //v_cent_gev(3,"Kaon","pos");
  //v_cent_gev(3,"Proton","pos");

/*   // для тараненко
  v_cent_rat(2, "Pion");
  v_cent_rat(2, "Kaon");
  v_cent_rat(2, "Proton");

  v_cent_rat(3, "Pion");
  v_cent_rat(3, "Kaon");
  v_cent_rat(3, "Proton");
*/

}