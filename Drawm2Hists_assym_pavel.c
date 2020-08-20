#include <TMath.h>
// particle: 1 - deuterons, 2 - tritons, 3 - He3, 4 - He4
Int_t Tparticle; Bool_t QMode = kFALSE, CFtracking = kTRUE;
Int_t GetStringInt = 0;

Double_t lorentzianPeak(Double_t *x, Double_t *par) {
  return (0.5*par[0]*par[1]/TMath::Pi()) /
        TMath::Max( 1.e-10,(x[0]-par[2])*(x[0]-par[2])
		+ .25*par[1]*par[1]);
		}

void Drawm2Hists_assym_pavel()
{
	std :: cout << "Enter particle" << std :: endl;
	std :: cin >> Tparticle;

	Double_t a, b, mom, c[5];

	TString inFile, inFile2;
	
	
	
	if (CFtracking) inFile = "/home/pepe-frog/Dubna/massSquare/m2Hists_asymmetric_pavel.root";
	else inFile = "/home/pepe-frog/Dubna/massSquare/m2Hists_asymmetric_pavel.root";

	gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C");
	

	//mpdloadlibs(kTRUE,kTRUE); // all libs
	
	TFile *iFile = new TFile(inFile);


	TH1D *De[60];
	TH1D *Tr[60];
	TH1D *He3[60];
	TH1D *He4[60];

	TH2D *m2_all; 
	TH2D *m2_de;
	TH2D *m2_tr; 
	TH2D *m2_he3;
	TH2D *m2_he4; 

	m2_de = (TH2D*)iFile->Get("m2_de");
	m2_tr = (TH2D*)iFile->Get("m2_tr");
	m2_he3 = (TH2D*)iFile->Get("m2_He3");
	m2_he4 = (TH2D*)iFile->Get("m2_He4");

	for (Int_t i = 18; i <= 57; i++){
		De[i] = (TH1D*)iFile->Get(Form("De_%i",i));
	}
	for (Int_t i = 1; i <= 60; i++){
		Tr[i] = (TH1D*)iFile->Get(Form("Tr_%i",i));
	}
	for (Int_t i = 1; i <= 60; i++){
		He3[i] = (TH1D*)iFile->Get(Form("He3_%i",i));
	}
	for (Int_t i = 1; i <= 60; i++){
		He4[i] = (TH1D*)iFile->Get(Form("He4_%i",i));
	}
	
	if (Tparticle == 1) /////////////////////////////////////// DEUTERONS
	{
		TF1 *Gaus[60];

		TCanvas *can[60];
		Double_t p = 1.6;

		for(Int_t i = 18; i <= 57; i++){

		//if (!QMode) can[i] = new TCanvas (Form("can%i",i), Form("deuterons, %1.2f < p < %1.2f, |eta| < 1.6",i*0.05,i*0.05+0.05), 1200, 800);


		De[i]->SetStats(kFALSE);
		De[i]->GetXaxis()->SetTitle("Energy deposit, ADC counts");
		De[i]->GetYaxis()->SetTitle("Counts");
		if (!QMode) De[i]->Draw();
		
		Gaus[i] = new TF1(Form("Gaus%i",i), "gaus", De[i]->GetXaxis()->GetXmin(), De[i]->GetXaxis()->GetXmax());
		De[i]->Fit(Gaus[i], "Q0R");
		if (!QMode) Gaus[i] -> Draw("same");
//		if(i<=31){
		//if (!QMode) can[i]->SaveAs(Form("/home/pepe-frog/Dubna/massSquare/pic/de/deuterons_%1.2f<p<%1.2f.png",i*0.05,i*0.05+0.05));
//		}
	}
	
	const Int_t nHists = 40;
	const Int_t k = 18;
	
	Double_t Xlow[nHists] = {0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2., 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85};
	Double_t Xhigh[nHists] = {0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2., 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9};
	Double_t X[nHists], Y[nHists], Y_err[nHists], Y_tail[nHists], X_err[nHists], minsigma[nHists], normerrminsigma[nHists], 
			 normminsigma[nHists], sigma1[nHists], errsigma1[nHists], sigma2[nHists], errsigma2[nHists], Y_tail_err[nHists];

	for (Int_t i=0; i<nHists; i++) {

			X[i] = (Xlow[i] + Xhigh[i]) / 2;

			Y[i] = Gaus[i+k]->GetParameter(1);
			Y_err[i] = Gaus[i+k]->GetParError(1);
			X_err[i] = 0.;


			sigma1[i] = Gaus[i+k]->GetParameter(2);

			errsigma1[i] = Gaus[i+k]->GetParError(2);
		}

	TGraphErrors *grDeuterons = new TGraphErrors(nHists,X,Y,X_err,Y_err);
	
	TCanvas *grcan1 = new TCanvas ("grcan1", "#m^{2} vs P, deuterons", 1200, 800);
	gPad->SetGrid();
	grDeuterons->SetTitle("");
	grDeuterons->GetYaxis()->SetTitle("#m^{2}");
	grDeuterons->GetYaxis()->SetTitleOffset(1.5);
	grDeuterons->GetXaxis()->SetTitle("Momentum, GeV/c");
	grDeuterons->SetMarkerStyle(21); grDeuterons->SetMarkerSize(2);
	grDeuterons->GetXaxis()->SetLabelSize(0.03); grDeuterons->GetYaxis()->SetLabelSize(0.03);
	grDeuterons->GetXaxis()->SetTitleSize(0.06); grDeuterons->GetYaxis()->SetTitleSize(0.06);
	grDeuterons->GetXaxis()->CenterTitle(); grDeuterons->GetYaxis()->CenterTitle();
	m2_de->Draw("colz");
	grDeuterons->Draw("same");


	grcan1->SaveAs("/home/pepe-frog/Dubna/massSquare/pic/de/fitting/m2_with_fit.png");
	
	TGraphErrors *grDeuterons_minsigmas = new TGraphErrors(nHists,X,sigma1,X_err,errsigma1);

	TCanvas *grcan4 = new TCanvas ("grcan4", "#m^{2} width vs P, deuterons", 1200, 800);
	gPad->SetGrid();
	grDeuterons_minsigmas->SetTitle("");
	grDeuterons_minsigmas->GetYaxis()->SetTitle("#sigma");
	grDeuterons_minsigmas->GetXaxis()->SetTitle("Momentum, GeV/c");
	grDeuterons_minsigmas->SetMarkerStyle(21);
	grDeuterons_minsigmas->SetMarkerSize(1);
	grDeuterons_minsigmas->GetXaxis()->SetLabelSize(0.03); grDeuterons_minsigmas->GetYaxis()->SetLabelSize(0.03);
	grDeuterons_minsigmas->GetXaxis()->SetTitleSize(0.05); grDeuterons_minsigmas->GetYaxis()->SetTitleSize(0.05);
	grDeuterons_minsigmas->GetXaxis()->CenterTitle(); grDeuterons_minsigmas->GetYaxis()->CenterTitle();
	grDeuterons_minsigmas->Draw("AP");


	grcan4->SaveAs("/home/pepe-frog/Dubna/massSquare/pic/de/fitting/sigma.png");
	}
	
	
	if (Tparticle == 2) /////////////////////////////////////////////// TRITONS
	{
		TF1 *Gaus[60];
		TCanvas *can[60];
		Double_t p = 1.6;

		for(Int_t i = 22; i <= 53; i++){

		if (!QMode) can[i] = new TCanvas (Form("can%i",i), Form("tritons, %1.2f < p < %1.2f, |eta| < 1.6",i*0.05,i*0.05+0.05), 1200, 800);

		Tr[i]->SetStats(kFALSE);
		Tr[i]->GetXaxis()->SetTitle("Energy deposit, ADC counts");
		Tr[i]->GetYaxis()->SetTitle("Counts");
		if (!QMode) Tr[i]->Draw();
		
		Gaus[i] = new TF1(Form("Gaus%i",i), "gaus", Tr[i]->GetXaxis()->GetXmin(), Tr[i]->GetXaxis()->GetXmax());
		Tr[i]->Fit(Gaus[i], "Q0R");

		if (!QMode) Gaus[i]->Draw("same");

		if (!QMode) can[i]->SaveAs(Form("/home/pepe-frog/Dubna/massSquare/pic/tr/tritons_%1.2f<p<%1.2f.png",i*0.05,i*0.05+0.05));


	}
	
	const Int_t nHists = 31;
	const Int_t k = 22;

	Double_t Xlow[nHists] = { 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2., 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6};
	Double_t Xhigh[nHists] = { 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2., 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65};
	Double_t X[nHists], Y[nHists], Y_err[nHists], Y_tail[nHists], X_err[nHists], minsigma[nHists], normerrminsigma[nHists], 
			 normminsigma[nHists], sigma1[nHists], errsigma1[nHists], sigma2[nHists], errsigma2[nHists], Y_tail_err[nHists];

	for (Int_t i=0; i<nHists; i++) {

			X[i] = (Xlow[i] + Xhigh[i]) / 2;

			Y[i] = Gaus[i+k]->GetParameter(1);
			Y_err[i] = Gaus[i+k]->GetParError(1);
			X_err[i] = 0.;


			sigma1[i] = Gaus[i+k]->GetParameter(2);

			errsigma1[i] = Gaus[i+k]->GetParError(2);
		}
	
	TGraphErrors *grTritons = new TGraphErrors(nHists,X,Y,X_err,Y_err);
	
	TCanvas *grcan1 = new TCanvas ("grcan1", "#m^{2} vs P, tritons", 1200, 800);
	gPad->SetGrid();
	grTritons->SetTitle("");
	grTritons->GetYaxis()->SetTitle("#m^{2}");
	grTritons->GetXaxis()->SetTitle("Momentum, GeV/c");
	grTritons->SetMarkerStyle(21); grTritons->SetMarkerSize(2);
	grTritons->GetXaxis()->SetLabelSize(0.03); grTritons->GetYaxis()->SetLabelSize(0.03);
	grTritons->GetXaxis()->SetTitleSize(0.05); grTritons->GetYaxis()->SetTitleSize(0.05);
	grTritons->GetXaxis()->CenterTitle(); grTritons->GetYaxis()->CenterTitle();
	m2_tr->Draw("colz");
	grTritons->Draw("same");
	
	grcan1->SaveAs("/home/pepe-frog/Dubna/massSquare/pic/tr/fiting/m2_with_fit.png");


	TGraphErrors *grTritons_minsigmas = new TGraphErrors(nHists,X,sigma1,X_err,errsigma1);

	TCanvas *grcan4 = new TCanvas ("grcan4", "#m^{2} width vs P, tritons", 1200, 800);
	gPad->SetGrid();
	grTritons_minsigmas->SetTitle("");
	grTritons_minsigmas->GetYaxis()->SetTitle("#sigma");
	grTritons_minsigmas->GetXaxis()->SetTitle("Momentum, GeV/c");
	grTritons_minsigmas->SetMarkerStyle(21); grTritons_minsigmas->SetMarkerSize(1);
	grTritons_minsigmas->GetXaxis()->SetLabelSize(0.03); grTritons_minsigmas->GetYaxis()->SetLabelSize(0.03);
	grTritons_minsigmas->GetXaxis()->SetTitleSize(0.05); grTritons_minsigmas->GetYaxis()->SetTitleSize(0.05);
	grTritons_minsigmas->GetXaxis()->CenterTitle(); grTritons_minsigmas->GetYaxis()->CenterTitle();
	grTritons_minsigmas->Draw("AP");
	

	grcan4->SaveAs("/home/pepe-frog/Dubna/massSquare/pic/tr/fiting/sigma.png");

}
	
	if (Tparticle == 3) //////////////////////////////////////// He^3
	{
		TF1 *Gaus[33];
		TF1 *AGaus[33];
		TF1 *Novosib[33];
		TCanvas *can[33];
		Double_t p =1.6;

		for(Int_t i = 3; i <= 33; i++){
		if(i<=31){
		if (!QMode) can[i] = new TCanvas (Form("can%i",i), Form("he3, %1.2f < p < %1.2f, |eta| < 1.6",i*0.05,i*0.05+0.05), 1200, 800);
		}
		if(i>31){
		if (!QMode) can[i] = new TCanvas (Form("can%i",i), Form("he3, %1.2f < p < %1.2f, |eta| < 1.6",p,p+0.1), 1200, 800);
		}
		He3[i]->SetStats(kFALSE);
		He3[i]->GetXaxis()->SetTitle("Energy deposit, ADC counts");
		He3[i]->GetYaxis()->SetTitle("Counts");
		if (!QMode) He3[i]->Draw();
		
		Gaus[i] = new TF1(Form("Gaus%i",i), "gaus", He3[i]->GetXaxis()->GetXmin(), He3[i]->GetXaxis()->GetXmax());
		He3[i]->Fit(Gaus[i], "Q0R");

		Novosib[i] = new TF1(Form("Novosib%i",i), Novosibirsk, He3[i]->GetXaxis()->GetXmin(), He3[i]->GetXaxis()->GetXmax(), 4);
		Novosib[i]->SetParameters(Gaus[i]->GetParameter(0), 0.01, Gaus[i]->GetParameter(2), Gaus[i]->GetParameter(1)); Novosib[i]->SetLineColor(kBlue);
		He3[i]->Fit(Form("Novosib%i",i),"Q0RW");
		
		if (!QMode) Novosib[i]->Draw("same");
		AGaus[i] = new TF1(Form("AGaus%i",i), AsymGaus, He3[i]->GetXaxis()->GetXmin(), He3[i]->GetXaxis()->GetXmax(), 4);
		AGaus[i]->SetParameters(Gaus[i]->GetParameter(0), Gaus[i]->GetParameter(1), Gaus[i]->GetParameter(2), 0.01);
		He3[i]->Fit(Form("AGaus%i",i),"Q0RW");
		if (!QMode) AGaus[i]->Draw("same");
/*		if(i<=31){
		if (!QMode) can[i]->SaveAs(Form("/home/pepe-frog/Dubna/dEdx/pic/he3/he3_%1.2f<p<%1.2f.png",i*0.05,i*0.05+0.05));
		}
		if(i>31){
		if (!QMode) can[i]->SaveAs(Form("/home/pepe-frog/Dubna/dEdx/pic/he3/he3_%1.2f<p<%1.2f.png",p,p+0.1));
		p = p + 0.1;
		}*/
	}
	
	
	const Int_t nHists = 25;
	const Int_t k = 6;

	Double_t Xlow[nHists] = { 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5};
	Double_t Xhigh[nHists] = { 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55};
	Double_t X[nHists], Y[nHists], Y_err[nHists], Y_tail[nHists], X_err[nHists], minsigma[nHists], normerrminsigma[nHists], 
			 normminsigma[nHists], sigma1[nHists], errsigma1[nHists], sigma2[nHists], errsigma2[nHists], Y_tail_err[nHists];

	FillMassive(AGaus, Novosib, nHists, Xlow, Xhigh, X, Y, Y_err, Y_tail, X_err, minsigma, 
		normerrminsigma, normminsigma, sigma1, errsigma1, sigma2, errsigma2, Y_tail_err, k);

	
	TGraphErrors *grHe3 = new TGraphErrors(nHists,X,Y,X_err,minsigma);
	
	TCanvas *grcan1 = new TCanvas ("grcan1", "dE/dx vs P, he3", 1200, 800);
	gPad->SetGrid();
	grHe3->SetTitle("");
	grHe3->GetYaxis()->SetTitle("dE/dx, ADC counts");
	grHe3->GetXaxis()->SetTitle("Momentum, GeV/c");
	grHe3->SetMarkerStyle(21); grHe3->SetMarkerSize(2);
	grHe3->GetXaxis()->SetLabelSize(0.03); grHe3->GetYaxis()->SetLabelSize(0.03);
	grHe3->GetXaxis()->SetTitleSize(0.05); grHe3->GetYaxis()->SetTitleSize(0.05);
	grHe3->GetXaxis()->CenterTitle(); grHe3->GetYaxis()->CenterTitle();
	m2_he3->Draw("colz");
	grHe3->Draw("same");
	
	TF1* fitparHe3BB1 = new TF1("fitparHe3BB1", "pol1(0)", 0., 0.35);
	grHe3->Fit(fitparHe3BB1, "Q0", "", 0.15, 0.35);
	TF1* fitparHe3BB2 = new TF1("fitparHe3BB2", "pol1(0)", 0.35, 0.55);
	grHe3->Fit(fitparHe3BB2, "Q0", "", 0.35, 0.55);
	TF1* fitparHe3BB3 = new TF1("fitparHe3BB3",fitparHe3BB->GetExpFormula(), 0.55, 1.8);
	grHe3->Fit(fitparHe3BB3, "Q0", "", 0.55, 1.8);
	
	TF1* fitparHe3BB4 = new TF1("fitparHe3BB4", "pol1(0)", 0., 0.325); fitparHe3BB4->SetParameters(fitparHe3BB1->GetParameters());
	TF1* fitparHe3BB5 = new TF1("fitparHe3BB5", "pol1(0)", 0.325, 0.55); fitparHe3BB5->SetParameters(fitparHe3BB2->GetParameters());
	TF1* fitparHe3BB6 = new TF1("fitparHe3BB6", fitparHe3BB->GetExpFormula(), 0.55, 1.8); fitparHe3BB6->SetParameters(fitparHe3BB3->GetParameters());
	fitparHe3BB4->Draw("same"); fitparHe3BB5->Draw("same"); fitparHe3BB6->Draw("same");
	grcan1->SaveAs("/home/pepe-frog/Dubna/dEdx/pic/he3/fiting/dedx_with_fit.png");

	for (Int_t i = 0; i < fitparHe3BB4->GetNpar(); i++) cout << "c[" << i << "] = " << fitparHe3BB4->GetParameter(i) << "    ";
	cout << endl;
	for (Int_t i = 0; i < fitparHe3BB5->GetNpar(); i++) cout << "c[" << i << "] = " << fitparHe3BB5->GetParameter(i) << "    ";
	cout << endl;
	for (Int_t i = 0; i < fitparHe3BB6->GetNpar(); i++) cout << "c[" << i << "] = " << fitparHe3BB6->GetParameter(i) << "    ";
	cout << endl;
	
	TGraphErrors *grHe3_minsigmas = new TGraphErrors(nHists,X,normminsigma,X_err,normerrminsigma);
	TGraphErrors *grHe3_tails = new TGraphErrors(nHists,X,Y_tail,X_err,Y_tail_err);
	
	TCanvas *grcan3 = new TCanvas ("grcan3", "Tails vs P, he3", 1200, 800);
	gPad->SetGrid();
	grHe3_tails->SetTitle("");
	grHe3_tails->GetYaxis()->SetTitle("#delta");
	grHe3_tails->GetXaxis()->SetTitle("Momentum, GeV/c");
	grHe3_tails->SetMarkerStyle(21); grHe3_tails->SetMarkerSize(2);
	grHe3_tails->GetXaxis()->SetLabelSize(0.03); grHe3_tails->GetYaxis()->SetLabelSize(0.03);
	grHe3_tails->GetXaxis()->SetTitleSize(0.05); grHe3_tails->GetYaxis()->SetTitleSize(0.05);
	grHe3_tails->GetXaxis()->CenterTitle(); grHe3_tails->GetYaxis()->CenterTitle();
	grHe3_tails->Draw("AP");
	TF1 *deltaall = new TF1("deltaall","pol3(0)",0.75,1.5);grHe3_tails->Fit(deltaall,"Q0R");
	TF1 *deltaHe3LowP = new TF1("deltaHe3LowP", "pol1(0)", 0.35, 0.7); grHe3_tails->Fit(deltaHe3LowP, "Q0R");
	//TF1 *deltaHe3MidP = new TF1("deltaHe3MidP", "pol1(0)", 0.685, 0.99); deltaHe3MidP->SetParameters(11.275, -11.666);
	//TF1 *deltaHe3HighP = new TF1("deltaHe3HighP", "pol1(0)", 0.99, 1.55); grHe3_tails->Fit(deltaHe3HighP, "Q0R");
	/*cout << "deltaHe3LowP[0] = " << deltaHe3LowP->GetParameter(0) << ", deltaHe3LowP[1] = " << deltaHe3LowP->GetParameter(1) << endl;
	cout << "deltaHe3HighP[0] = " << deltaHe3HighP->GetParameter(0) << ", deltaHe3HighP[1] = " << deltaHe3HighP->GetParameter(1) << endl;
	*/deltaHe3LowP->Draw("same"); /*deltaHe3MidP->Draw("same"); deltaHe3HighP->Draw("same");*/deltaall->Draw("same");
	grcan3->SaveAs("/home/pepe-frog/Dubna/dEdx/pic/he3/fiting/delta.png");

	TCanvas *grcan4 = new TCanvas ("grcan4", "dE/dx width vs P, he3", 1200, 800);
	gPad->SetGrid();
	grHe3_minsigmas->SetTitle("");
	grHe3_minsigmas->GetYaxis()->SetTitle("#sigma");
	grHe3_minsigmas->GetXaxis()->SetTitle("Momentum, GeV/c");
	grHe3_minsigmas->SetMarkerStyle(21); grHe3_minsigmas->SetMarkerSize(2);
	grHe3_minsigmas->GetXaxis()->SetLabelSize(0.03); grHe3_minsigmas->GetYaxis()->SetLabelSize(0.03);
	grHe3_minsigmas->GetXaxis()->SetTitleSize(0.05); grHe3_minsigmas->GetYaxis()->SetTitleSize(0.05);
	grHe3_minsigmas->GetXaxis()->CenterTitle(); grHe3_minsigmas->GetYaxis()->CenterTitle();
	grHe3_minsigmas->Draw("AP");
	
	TF1 *polHe3LowP = new TF1("polHe3LowP", "pol1(0)", 0., 0.675); polHe3LowP->SetParameters(0.44225, -0.47);
	TF1 *polHe3MidP = new TF1("polHe3MidP", "pol1(0)", 0.675, 1.475); polHe3MidP->SetParameters(0.17984375, -0.08125);
	TF1 *polHe3HighP = new TF1("polHe3HighP", "pol1(0)", 1.475, 3.); polHe3HighP->SetParameters(-0.261818, 0.21818);
	cout << "polHe3HighP[0] = " << polHe3HighP->GetParameter(0) << ", polHe3HighP[1] = " << polHe3HighP->GetParameter(1) << endl;
	polHe3LowP->Draw("same"); polHe3MidP->Draw("same"); polHe3HighP->Draw("same");
	grcan4->SaveAs("/home/pepe-frog/Dubna/dEdx/pic/he3/fiting/sigma.png");

	}
	
	if (Tparticle == 4) //////////////////////////////////////// He^4
	{
	
		TF1 *Gaus[31];
		TF1 *AGaus[31];
		TF1 *Novosib[31];
		TCanvas *can[31];
		Double_t p = 1.6;

		for(Int_t i = 4; i <= 31; i++){	
		if(i<=31){
		if (!QMode) can[i] = new TCanvas (Form("can%i",i), Form("he4, %1.2f < p < %1.2f, |eta| < 1.6",i*0.05,i*0.05+0.05), 1200, 800);
		}
		if(i>31){
		if (!QMode) can[i] = new TCanvas (Form("can%i",i), Form("he4, %1.2f < p < %1.2f, |eta| < 1.6",p,p+0.1), 1200, 800);
		}
		He4[i]->SetStats(kFALSE);
		He4[i]->GetXaxis()->SetTitle("Energy deposit, ADC counts");
		He4[i]->GetYaxis()->SetTitle("Counts");
		if (!QMode) He4[i]->Draw();
		
		Gaus[i] = new TF1(Form("Gaus%i",i), "gaus", He4[i]->GetXaxis()->GetXmin(), He4[i]->GetXaxis()->GetXmax());
		He4[i]->Fit(Gaus[i], "Q0R");

		Novosib[i] = new TF1(Form("Novosib%i",i), Novosibirsk, He4[i]->GetXaxis()->GetXmin(), He4[i]->GetXaxis()->GetXmax(), 4);
		Novosib[i]->SetParameters(Gaus[i]->GetParameter(0), 0.01, Gaus[i]->GetParameter(2), Gaus[i]->GetParameter(1)); Novosib[i]->SetLineColor(kBlue);
		He4[i]->Fit(Form("Novosib%i",i),"Q0RW");
		
		if (!QMode) Novosib[i]->Draw("same");
		AGaus[i] = new TF1(Form("AGaus%i",i), AsymGaus, He4[i]->GetXaxis()->GetXmin(), He4[i]->GetXaxis()->GetXmax(), 4);
		AGaus[i]->SetParameters(Gaus[i]->GetParameter(0), Gaus[i]->GetParameter(1), Gaus[i]->GetParameter(2), 0.01);
		He4[i]->Fit(Form("AGaus%i",i),"Q0RW");
		if (!QMode) AGaus[i]->Draw("same");
/*		if(i<=31){
		if (!QMode) can[i]->SaveAs(Form("/home/pepe-frog/Dubna/dEdx/pic/he4/he4_%1.2f<p<%1.2f.png",i*0.05,i*0.05+0.05));
		}
		if(i>31){
		if (!QMode) can[i]->SaveAs(Form("/home/pepe-frog/Dubna/dEdx/pic/he4/he4_%1.2f<p<%1.2f.png",p,p+0.1));
		p = p + 0.1;
		}*/
	}
	
	const Int_t nHists = 23;
	const Int_t k =7;

	Double_t Xlow[nHists] = { 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45};
	Double_t Xhigh[nHists] = { 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5};
		Double_t X[nHists], Y[nHists], Y_err[nHists], Y_tail[nHists], X_err[nHists], minsigma[nHists], normerrminsigma[nHists], 
			 normminsigma[nHists], sigma1[nHists], errsigma1[nHists], sigma2[nHists], errsigma2[nHists], Y_tail_err[nHists];

	FillMassive(AGaus, Novosib, nHists, Xlow, Xhigh, X, Y, Y_err, Y_tail, X_err, minsigma, 
		normerrminsigma, normminsigma, sigma1, errsigma1, sigma2, errsigma2, Y_tail_err, k);

	TGraphErrors *grHe4 = new TGraphErrors(nHists,X,Y,X_err,minsigma);
	
	TCanvas *grcan1 = new TCanvas ("grcan1", "dE/dx vs P, he4", 1200, 800);
	gPad->SetGrid();
	grHe4->SetTitle("");
	grHe4->GetYaxis()->SetTitle("dE/dx, ADC counts");
	grHe4->GetXaxis()->SetTitle("Momentum, GeV/c");
	grHe4->SetMarkerStyle(21); grHe4->SetMarkerSize(2);
	grHe4->GetXaxis()->SetLabelSize(0.03); grHe4->GetYaxis()->SetLabelSize(0.03);
	grHe4->GetXaxis()->SetTitleSize(0.05); grHe4->GetYaxis()->SetTitleSize(0.05);
	grHe4->GetXaxis()->CenterTitle(); grHe4->GetYaxis()->CenterTitle();
	m2_he4->Draw("colz");
	grHe4->Draw("same");
	
	TF1* fitparHe4BB1 = new TF1("fitparHe4BB1", "pol1(0)", 0., 0.35);
	grHe4->Fit(fitparHe4BB1, "Q0", "", 0.2, 0.35);
	TF1* fitparHe4BB2 = new TF1("fitparHe4BB2", "pol1(0)", 0.35, 0.5);
	grHe4->Fit(fitparHe4BB2, "Q0", "", 0.35, 0.5);
	TF1* fitparHe4BB3 = new TF1("fitparHe4BB3",fitparHe4BB->GetExpFormula(), 0.45, 1.6);
	grHe4->Fit(fitparHe4BB3, "Q0", "", 0.45, 1.6);
	
	TF1* fitparHe4BB4 = new TF1("fitparHe4BB4", "pol1(0)", 0., 0.365); fitparHe4BB4->SetParameters(fitparHe4BB1->GetParameters());
	TF1* fitparHe4BB5 = new TF1("fitparHe4BB5", "pol1(0)", 0.365, 0.475); fitparHe4BB5->SetParameters(fitparHe4BB2->GetParameters());
	TF1* fitparHe4BB6 = new TF1("fitparHe4BB6", fitparHe4BB->GetExpFormula(), 0.475, 1.6); fitparHe4BB6->SetParameters(fitparHe4BB3->GetParameters());
	fitparHe4BB4->Draw("same"); fitparHe4BB5->Draw("same"); fitparHe4BB6->Draw("same");
		
	grcan1->SaveAs("/home/pepe-frog/Dubna/dEdx/pic/he4/fiting/dedx_with_fit.png");

	for (Int_t i = 0; i < fitparHe4BB4->GetNpar(); i++) cout << "c[" << i << "] = " << fitparHe4BB4->GetParameter(i) << "    ";
	cout << endl;
	for (Int_t i = 0; i < fitparHe4BB5->GetNpar(); i++) cout << "c[" << i << "] = " << fitparHe4BB5->GetParameter(i) << "    ";
	cout << endl;
	for (Int_t i = 0; i < fitparHe4BB6->GetNpar(); i++) cout << "c[" << i << "] = " << fitparHe4BB6->GetParameter(i) << "    ";
	cout << endl;
	
	TGraphErrors *grHe4_minsigmas = new TGraphErrors(nHists,X,normminsigma,X_err,normerrminsigma);
	TGraphErrors *grHe4_tails = new TGraphErrors(nHists,X,Y_tail,X_err,Y_tail_err);
	
	TCanvas *grcan3 = new TCanvas ("grcan3", "Tails vs P, he4", 1200, 800);
	gPad->SetGrid();
	grHe4_tails->SetTitle("");
	grHe4_tails->GetYaxis()->SetTitle("#delta");
	grHe4_tails->GetXaxis()->SetTitle("Momentum, GeV/c");
	grHe4_tails->SetMarkerStyle(21); grHe4_tails->SetMarkerSize(2);
	grHe4_tails->GetXaxis()->SetLabelSize(0.03); grHe4_tails->GetYaxis()->SetLabelSize(0.03);
	grHe4_tails->GetXaxis()->SetTitleSize(0.05); grHe4_tails->GetYaxis()->SetTitleSize(0.05);
	grHe4_tails->GetXaxis()->CenterTitle(); grHe4_tails->GetYaxis()->CenterTitle();
	grHe4_tails->Draw("AP");
	TF1 *deltaHe4LowP = new TF1("deltaHe4LowP", "pol1(0)", 0., 0.575); deltaHe4LowP->SetParameters(1.45, -2.);
	TF1 *deltaHe4MidP = new TF1("deltaHe4MidP", "pol1(0)", 0.575, 0.875); grHe4_tails->Fit(deltaHe4MidP, "Q0R");
	TF1 *deltaHe4HighP = new TF1("deltaHe4HighP", "pol2(0)", 0.875, 1.6); grHe4_tails->Fit(deltaHe4HighP, "Q0R");
	cout << "deltaHe4MidP[0] = " << deltaHe4MidP->GetParameter(0) << ", deltaHe4MidP[1] = " << deltaHe4MidP->GetParameter(1) << endl;
	cout << "deltaHe4HighP[0] = " << deltaHe4HighP->GetParameter(0) << ", deltaHe4HighP[1] = " << deltaHe4HighP->GetParameter(1) << endl;
	deltaHe4LowP->Draw("same"); deltaHe4MidP->Draw("same"); deltaHe4HighP->Draw("same");
		grcan3->SaveAs("/home/pepe-frog/Dubna/dEdx/pic/he4/fiting/delta.png");

	TCanvas *grcan4 = new TCanvas ("grcan4", "dE/dx width vs P, he4", 1200, 800);
	gPad->SetGrid();
	grHe4_minsigmas->SetTitle("");
	grHe4_minsigmas->GetYaxis()->SetTitle("#sigma");
	grHe4_minsigmas->GetXaxis()->SetTitle("Momentum, GeV/c");
	grHe4_minsigmas->SetMarkerStyle(21); grHe4_minsigmas->SetMarkerSize(2);
	grHe4_minsigmas->GetXaxis()->SetLabelSize(0.03); grHe4_minsigmas->GetYaxis()->SetLabelSize(0.03);
	grHe4_minsigmas->GetXaxis()->SetTitleSize(0.05); grHe4_minsigmas->GetYaxis()->SetTitleSize(0.05);
	grHe4_minsigmas->GetXaxis()->CenterTitle(); grHe4_minsigmas->GetYaxis()->CenterTitle();
	grHe4_minsigmas->Draw("AP");
	

	TF1 *polHe4LowP = new TF1("polHe4LowP", "pol1(0)", 0., 0.475); polHe4LowP->SetParameters(0.37188, -0.225);
	TF1 *polHe4Mid1P = new TF1("polHe4Mid1P", "pol1(0)", 0.475, 0.7); polHe4Mid1P->SetParameters(0.19375, 0.15);
	TF1 *polHe4Mid2P = new TF1("polHe4Mid2P", "pol1(0)", 0.7, 0.925); polHe4Mid2P->SetParameters(0.72231, -0.6025);
	TF1 *polHe4HighP = new TF1("polHe4HighP", "pol1(0)", 0.925, 3.); polHe4HighP->SetParameters(0.32477, -0.17273);
	polHe4LowP->Draw("same"); polHe4Mid1P->Draw("same"); polHe4Mid2P->Draw("same"); polHe4HighP->Draw("same");
		grcan4->SaveAs("/home/pepe-frog/Dubna/dEdx/pic/he4/fiting/sigma.png");

	}
}
