#include <TMath.h>
// particle: 1 - deuterons, 2 - tritons, 3 - He3, 4 - He4
Int_t Tparticle; Bool_t QMode = kFALSE, CFtracking = kTRUE;
Int_t GetStringInt = 0;

Double_t GetDedxDeParam(TF1 *fitparDeBB, Double_t p){
   Double_t dedx=fitparDeBB->Eval(p);
   return dedx;
}

Double_t GetDedxTrParam(TF1 *fitparTrBB,Double_t p){
   Double_t dedx=fitparTrBB->Eval(p);
   return dedx;
}

Double_t GetDedxHe3Param(TF1 *fitparHe3BB,Double_t p){
   Double_t dedx=fitparHe3BB->Eval(p);
   return dedx;
}

Double_t GetDedxHe4Param(TF1 *fitparHe4BB,Double_t p){
   Double_t dedx=fitparHe4BB->Eval(p);
   return dedx;
}

Double_t AsymGaus(Double_t *x, Double_t *par){
	Double_t peak = par[1];
	if (x[0] < peak)
	return par[0]*TMath::Exp( -0.5 * TMath::Power( ( (x[0] - par[1]) / par[2] ), 2 ));
	else
	return par[0]*TMath::Exp( -0.5 * TMath::Power( ( (x[0] - par[1]) / ((1. + par[3]) * par[2]) ), 2 ));
}

Double_t Novosibirsk(Double_t *x, Double_t *par){
	Double_t tail = par[1];
	Double_t width = par[2];
	Double_t peak = par[3];
	if (TMath::Abs(tail) < 1.e-7) 
	return par[0]*TMath::Exp( -0.5 * TMath::Power( ( (x[0] - peak) / width ), 2 ));
	Double_t arg = 1.0 - ( x[0] - peak ) * tail / width;
	if (arg < 1.e-6) return 0.0;     //Argument of logaritem negative. Real continuation -> function equals zero
	
	Double_t log = TMath::Log(arg);   
	static const Double_t xi = 2.3548200450309494; // 2 Sqrt( Ln(4) )
	Double_t width_zero = ( 2.0 / xi ) * TMath::ASinH( tail * xi * 0.5 );
	Double_t width_zero2 = width_zero * width_zero;
	Double_t exponent = ( -0.5 / (width_zero2) * log * log ) - ( width_zero2 * 0.5 );
	
	return par[0]*TMath::Exp(exponent);
}
Double_t lorentzianPeak(Double_t *x, Double_t *par) {
  return (0.5*par[0]*par[1]/TMath::Pi()) /
        TMath::Max( 1.e-10,(x[0]-par[2])*(x[0]-par[2])
		+ .25*par[1]*par[1]);
		}

void FillMassive(TF1* AGaus[], TF1* Novosib[], Int_t nHists, Double_t Xlow[], Double_t Xhigh[], Double_t X[], Double_t Y[], Double_t Y_err[], Double_t Y_tail[], Double_t X_err[],
				 Double_t minsigma[], Double_t normerrminsigma[], Double_t normminsigma[], Double_t sigma1[], 
				 Double_t errsigma1[], Double_t sigma2[], Double_t errsigma2[], Double_t Y_tail_err[], Int_t k){
	for (Int_t i=0; i<nHists; i++) {

		X[i] = (Xlow[i] + Xhigh[i]) / 2;

		Y[i] = Novosib[i+k]->GetParameter(3);

		Y_err[i] = Novosib[i+k]->GetParError(3);

		Y_tail[i] = AGaus[i+k]->GetParameter(3);

		Y_tail_err[i] = AGaus[i+k]->GetParError(3);

		X_err[i] = 0.;

		sigma1[i] = AGaus[i+k]->GetParameter(2);

		errsigma1[i] = AGaus[i+k]->GetParError(2);

		sigma2[i] = (1. + AGaus[i+k]->GetParameter(3)) * AGaus[i+k]->GetParameter(2);

		errsigma2[i] = (1. + AGaus[i+k]->GetParameter(3)) * AGaus[i+k]->GetParError(2) + AGaus[i+k]->GetParameter(2) * AGaus[i+k]->GetParError(3);

		errsigma1[i] = errsigma1[i] / Y[i] + Y_err[i] * sigma1[i] / (Y[i]*Y[i]); errsigma2[i] = errsigma2[i] / Y[i] + Y_err[i] * sigma2[i] / (Y[i]*Y[i]); 

		if (sigma1[i] <= sigma2[i]) {
		
			minsigma[i] = sigma1[i]; 
			normerrminsigma[i] = errsigma1[i];
		} 
		else {
		
			minsigma[i] = sigma2[i]; 
			normerrminsigma[i] = errsigma2[i];
		}
		
		normminsigma[i] = minsigma[i] / Y[i];
	}

}

void DrawdedxHists_asymmetric_pavel()
{
	std :: cout << "Enter particle" << std :: endl;
	std :: cin >> Tparticle;

	Double_t a, b, mom, c[5];

	TString inFile, inFile2;
	
	TF1* fitparDeBB = new TF1("fitparDeBB","[0]/pow(x/sqrt(x*x+3.52),[3])*([1]-pow(x/sqrt(x*x+3.52),[3])-log([2]+pow(1./(x/1.876),[4])) )", 0., 3.);
	if (CFtracking) fitparDeBB->SetParameters(-2564.17,0.0987747,1.32963,0.481253,3.51363);
	else fitparDeBB->SetParameters(4.89e-08,22.953,1.876e-05,2.59,4.153);
	
	TF1* fitparTrBB = new TF1("fitparTrBB","[0]/pow(x/sqrt(x*x+7.89),[3])*([1]-pow(x/sqrt(x*x+7.89),[3])-log([2]+pow(1./(x/2.81),[4])) )", 0., 3.);
	if (CFtracking) fitparTrBB->SetParameters(-7339.89,1.76536,3.22294,0.0533274,2.90796);
	else fitparTrBB->SetParameters(3.428e-07,3.768,-1.74e-02,2.284,0.963);
	
	TF1* fitparHe3BB = new TF1("fitparHe3BB","[0]*((1+(x/1.4047)**2)/pow(x/1.4047,[3])*([1]+[2]*log(1+pow(x/1.4047,2)))-1.)", 0., 3.);
	if (CFtracking) fitparHe3BB->SetParameters(2.86201,2.10168,0.274807,1.86774);
	else fitparHe3BB->SetParameters(2.86201e-06,2.10168,0.274807,1.86774);
	
	TF1* fitparHe4BB = new TF1("fitparHe4BB","[0]*((1+(x/1.863)**2)/pow(x/1.863,[3])*([1]+[2]*log(1+pow(x/1.863,2)))-1.)", 0., 3.);
	if (CFtracking) fitparHe4BB->SetParameters(2.96,2.085,0.256,1.85);
	else fitparHe4BB->SetParameters(2.96e-06,2.085,0.256,1.85);
	
	if (CFtracking) inFile = "/home/pepe-frog/Dubna/dEdx/dedxHists_asymmetric_pavel.root";
	else inFile = "/home/pepe-frog/Dubna/dEdx/dedxHists_asymmetric_pavel.root";

	gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C");
	

	//mpdloadlibs(kTRUE,kTRUE); // all libs
	
	TFile *iFile = new TFile(inFile);


	TH1D *De[41];
	TH1D *Tr[41];
	TH1D *He3[32];
	TH1D *He4[34];

	TH2D *banana_all; 
	TH2D *banana_de;
	TH2D *banana_tr; 
	TH2D *banana_he3;
	TH2D *banana_he4; 

	banana_de = (TH2D*)iFile->Get("dedxVSp_De");
	banana_tr = (TH2D*)iFile->Get("dedxVSp_Tr");
	banana_he3 = (TH2D*)iFile->Get("dedxVSp_He3");
	banana_he4 = (TH2D*)iFile->Get("dedxVSp_He4");

	for (Int_t i = 3; i <= 40; i++){
		De[i] = (TH1D*)iFile->Get(Form("De%i",i));
	}
	for (Int_t i = 4; i <= 40; i++){
		Tr[i] = (TH1D*)iFile->Get(Form("Tr%i",i));
	}
	for (Int_t i = 3; i <= 33; i++){
		He3[i] = (TH1D*)iFile->Get(Form("He3_%i",i));
	}
	for (Int_t i = 4; i <= 31; i++){
		He4[i] = (TH1D*)iFile->Get(Form("He4_%i",i));
	}
	
	if (Tparticle == 1) /////////////////////////////////////// DEUTERONS
	{
		TF1 *Gaus[41];
		TF1 *AGaus[41];
		TF1 *Novosib[41];
		TCanvas *can[41];
		Double_t p = 1.6;

		for(Int_t i = 3; i <= 40; i++){
/*		if(i<=31){
		if (!QMode) can[i] = new TCanvas (Form("can%i",i), Form("deuterons, %1.2f < p < %1.2f, |eta| < 1.6",i*0.05,i*0.05+0.05), 1200, 800);
		}
		if(i>31){
		if (!QMode) can[i] = new TCanvas (Form("can%i",i), Form("deuterons, %1.2f < p < %1.2f, |eta| < 1.6",p,p+0.1), 1200, 800);
		}*/
		De[i]->SetStats(kFALSE);
		De[i]->GetXaxis()->SetTitle("Energy deposit, ADC counts");
		De[i]->GetYaxis()->SetTitle("Counts");
		if (!QMode) De[i]->Draw();
		
		Gaus[i] = new TF1(Form("Gaus%i",i), "gaus", De[i]->GetXaxis()->GetXmin(), De[i]->GetXaxis()->GetXmax());
		De[i]->Fit(Gaus[i], "Q0R");

		Novosib[i] = new TF1(Form("Novosib%i",i), Novosibirsk, De[i]->GetXaxis()->GetXmin(), De[i]->GetXaxis()->GetXmax(), 4);
		Novosib[i]->SetParameters(Gaus[i]->GetParameter(0), 0.01, Gaus[i]->GetParameter(2), Gaus[i]->GetParameter(1)); Novosib[i]->SetLineColor(kBlue);
		De[i]->Fit(Form("Novosib%i",i),"Q0RW");
		
		if (!QMode) Novosib[i]->Draw("same");
		AGaus[i] = new TF1(Form("AGaus%i",i), AsymGaus, De[i]->GetXaxis()->GetXmin(), De[i]->GetXaxis()->GetXmax(), 4);
		AGaus[i]->SetParameters(Gaus[i]->GetParameter(0), Gaus[i]->GetParameter(1), Gaus[i]->GetParameter(2), 0.01);
		De[i]->Fit(Form("AGaus%i",i),"Q0RW");
		if (!QMode) AGaus[i]->Draw("same");
	/*	if(i<=31){
		if (!QMode) can[i]->SaveAs(Form("/home/pepe-frog/Dubna/dEdx/pic/deutrons/deuterons_%1.2f<p<%1.2f.png",i*0.05,i*0.05+0.05));
		}
		if(i>31){
		if (!QMode) can[i]->SaveAs(Form("/home/pepe-frog/Dubna/dEdx/pic/deutrons/deuterons_%1.2f<p<%1.2f.png",p,p+0.1));
		p = p + 0.1;
		}*/
	}
	
	const Int_t nHists = 38;
	const Int_t k = 3;
	
	Double_t Xlow[nHists] = {0.15,0.20,0.25,0.30,0.35,0.40,0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.4, 2.6 };
	Double_t Xhigh[nHists] = {0.20,0.25,0.30,0.35,0.40,0.45,0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.4, 2.6, 3. };
	Double_t X[nHists], Y[nHists], Y_err[nHists], Y_tail[nHists], X_err[nHists], minsigma[nHists], normerrminsigma[nHists], 
			 normminsigma[nHists], sigma1[nHists], errsigma1[nHists], sigma2[nHists], errsigma2[nHists], Y_tail_err[nHists];

	FillMassive(AGaus, Novosib, nHists, Xlow, Xhigh, X, Y, Y_err, Y_tail, X_err, minsigma, 
		normerrminsigma, normminsigma, sigma1, errsigma1, sigma2, errsigma2, Y_tail_err, k);
	

	TGraphErrors *grDeuterons = new TGraphErrors(nHists,X,Y,X_err,minsigma);
	
	TCanvas *grcan1 = new TCanvas ("grcan1", "dE/dx vs P, deuterons", 1200, 800);
	gPad->SetGrid();
	grDeuterons->SetTitle("");
	grDeuterons->GetYaxis()->SetTitle("dE/dx, ADC counts");
	grDeuterons->GetYaxis()->SetTitleOffset(1.5);
	grDeuterons->GetXaxis()->SetTitle("Momentum, GeV/c");
	grDeuterons->SetMarkerStyle(21); grDeuterons->SetMarkerSize(2);
	grDeuterons->GetXaxis()->SetLabelSize(0.03); grDeuterons->GetYaxis()->SetLabelSize(0.03);
	grDeuterons->GetXaxis()->SetTitleSize(0.06); grDeuterons->GetYaxis()->SetTitleSize(0.06);
	grDeuterons->GetXaxis()->CenterTitle(); grDeuterons->GetYaxis()->CenterTitle();
	banana_de->Draw("colz");
	grDeuterons->Draw("same");
	
	TF1* fitparDeBB1 = new TF1("fitparDeBB1",fitparDeBB->GetExpFormula(), 0., 3.);
	
	grDeuterons->Fit(fitparDeBB1, "Q0", "", Xlow[1], Xlow[1] + (Xhigh[nHists-1] - Xlow[1]) / 3.);
	
	for (Int_t i = 0; i < fitparDeBB->GetNpar(); i++) {
		c[i] = fitparDeBB1->GetParameter(i);
	}
	
	TF1* fitparDeBB2 = new TF1("fitparDeBB2",fitparDeBB->GetExpFormula(), 0., 3.);
	
	for (Int_t i = 0; i < fitparDeBB->GetNpar(); i++) {
		fitparDeBB2->SetParameter(i, c[i]);
	}
	
	grDeuterons->Fit(fitparDeBB2, "Q0", "", Xlow[1], Xlow[1] + 2. * (Xhigh[nHists-1] - Xlow[1]) / 3.);
	
	for (Int_t i = 0; i < fitparDeBB->GetNpar(); i++){

	c[i] = fitparDeBB2->GetParameter(i);
	}
	
	TF1* fitparDeBB3 = new TF1("fitparDeBB3",fitparDeBB->GetExpFormula(), 0., 3.);
	
	for (Int_t i = 0; i < fitparDeBB->GetNpar(); i++) {
		fitparDeBB3->SetParameter(i, c[i]);
	}
	
	grDeuterons->Fit(fitparDeBB3, "Q0", "", Xlow[1], Xhigh[nHists-1]);
	
	for (Int_t i = 0; i < fitparDeBB->GetNpar(); i++){ 
		c[i] = fitparDeBB3->GetParameter(i);
	}
	
	for (Int_t i = 0; i < fitparDeBB->GetNpar(); i++) {
		fitparDeBB->SetParameter(i, c[i]);
	}
	
	TF1* fitparDeBB6 = new TF1("fitparDeBB6", "pol1(0)", 0.15, 0.25);
	
	grDeuterons->Fit(fitparDeBB6, "Q0R");
	
	TF1* fitparDeBB4 = new TF1("fitparDeBB4", "pol1(0)", 0., 0.4);
	
	grDeuterons->Fit(fitparDeBB4, "Q0", "", 0.2, 0.4);
	
	TF1* fitparDeBB5 = new TF1("fitparDeBB5",fitparDeBB->GetExpFormula(), 0.365, 3.);
	
	grDeuterons->Fit(fitparDeBB5, "Q0", "", 0.35, 3.);
	
	TF1* fitparDeBB7 = new TF1("fitparDeBB7", "pol1(0)", 0., 0.245); 
	fitparDeBB7->SetParameters(fitparDeBB6->GetParameters());
	
	TF1* fitparDeBB8 = new TF1("fitparDeBB8", "pol1(0)", 0.245, 0.385);
	 fitparDeBB8->SetParameters(fitparDeBB4->GetParameters());
	
	TF1* fitparDeBB9 = new TF1("fitparDeBB9", fitparDeBB->GetExpFormula(), 0.374, 3.); 
	fitparDeBB9->SetParameters(fitparDeBB5->GetParameters());
	
	fitparDeBB7->Draw("same"); 
	fitparDeBB8->Draw("same"); 
	fitparDeBB9->Draw("same");



	grcan1->SaveAs("/home/pepe-frog/Dubna/dEdx/pic/deutrons/fiting/dedx_with_fit.png");
	for (Int_t i = 0; i < fitparDeBB7->GetNpar(); i++) cout << "c[" << i << "] = " << fitparDeBB7->GetParameter(i) << "    ";
	cout << endl;
	for (Int_t i = 0; i < fitparDeBB8->GetNpar(); i++) cout << "c[" << i << "] = " << fitparDeBB8->GetParameter(i) << "    ";
	cout << endl;
	for (Int_t i = 0; i < fitparDeBB9->GetNpar(); i++) cout << "c[" << i << "] = " << fitparDeBB9->GetParameter(i) << "    ";
	cout << endl;
	
	TF1* parDeBBMultX = new TF1("parDeBBMultX","x*("+fitparDeBB->GetExpFormula()+")", 0., 3.);
	
	for (Int_t k=0; k<fitparDeBB->GetNpar(); k++) {
		parDeBBMultX->SetParameter(k, fitparDeBB9->GetParameter(k));
	}
	
	TF1* parDeBBMultX1 = new TF1("parDeBBMultX1","x*([0] + [1] * x)", 0., 3.); parDeBBMultX1->SetParameters(fitparDeBB6->GetParameters());
	
	TF1* parDeBBMultX2 = new TF1("parDeBBMultX2","x*([0] + [1] * x)", 0., 3.); parDeBBMultX2->SetParameters(fitparDeBB4->GetParameters());
	
	Double_t Y_ratio[nHists], Y_ratio_err[nHists];
	
	for (Int_t i = 0; i < nHists; i++)
	{
		a = Xlow[i]; b = Xhigh[i];
		if (b <= 0.25) { mom = parDeBBMultX1->Integral(a,b)/(fitparDeBB7->Integral(a,b)); Y_ratio[i]=Y[i]/(GetDedxDeParam(fitparDeBB7, mom)); }
		else if (b <= 0.4) { mom = parDeBBMultX2->Integral(a,b)/(fitparDeBB8->Integral(a,b)); Y_ratio[i]=Y[i]/(GetDedxDeParam(fitparDeBB8, mom)); }
		else { mom = parDeBBMultX->Integral(a,b)/(fitparDeBB9->Integral(a,b)); Y_ratio[i]=Y[i]/(GetDedxDeParam(fitparDeBB9, mom)); }
		Y_ratio_err[i]=normminsigma[i];
	}
	
	TGraphErrors *grDeuterons_minsigmas = new TGraphErrors(nHists,X,normminsigma,X_err,normerrminsigma);
	TGraphErrors *grDeuterons_tails = new TGraphErrors(nHists,X,Y_tail,X_err,Y_tail_err);
	TGraphErrors *grDeuterons_ratio = new TGraphErrors(nHists,X,Y_ratio,X_err,Y_ratio_err);
	
	TCanvas *grcan2 = new TCanvas ("grcan2", "Ratio, deuterons", 1200, 800);
	gPad->SetGrid();
	grDeuterons_ratio->SetTitle("");
	grDeuterons_ratio->GetYaxis()->SetTitle("<dE/dx> / <Pid Bethe-Bloch>");
	grDeuterons_ratio->GetXaxis()->SetTitle("Momentum, GeV/c");
	grDeuterons_ratio->SetMarkerStyle(21); grDeuterons_ratio->SetMarkerSize(2);
	grDeuterons_ratio->GetXaxis()->SetLabelSize(0.03); grDeuterons_ratio->GetYaxis()->SetLabelSize(0.03);
	grDeuterons_ratio->GetXaxis()->SetTitleSize(0.05); grDeuterons_ratio->GetYaxis()->SetTitleSize(0.05);
	grDeuterons_ratio->GetXaxis()->CenterTitle(); grDeuterons_ratio->GetYaxis()->CenterTitle();
	grDeuterons_ratio->GetYaxis()->SetRangeUser(0.8,1.2);
	grDeuterons_ratio->Draw("AP");
	
	TCanvas *grcan3 = new TCanvas ("grcan3", "Tails vs P, deuterons", 1200, 800);
	gPad->SetGrid();
	grDeuterons_tails->SetTitle("");
	grDeuterons_tails->GetYaxis()->SetTitle("#delta");
	grDeuterons_tails->GetXaxis()->SetTitle("Momentum, GeV/c");
	grDeuterons_tails->SetMarkerStyle(21); grDeuterons_tails->SetMarkerSize(1);
	grDeuterons_tails->GetXaxis()->SetLabelSize(0.03); grDeuterons_tails->GetYaxis()->SetLabelSize(0.03);
	grDeuterons_tails->GetXaxis()->SetTitleSize(0.05); grDeuterons_tails->GetYaxis()->SetTitleSize(0.05);
	grDeuterons_tails->GetXaxis()->CenterTitle(); grDeuterons_tails->GetYaxis()->CenterTitle();
	grDeuterons_tails->Draw("AP");
	

	TF1 *deltaDe1 = new TF1("deltaDe1","pol1",1.65,3.06);
	grDeuterons_tails -> Fit(deltaDe1,"Q0R");
	TF1 *deltaDe2 = new TF1("deltaDe2","pol1",0.98,1.68);
	grDeuterons_tails -> Fit(deltaDe2,"Q0R");
	TF1 *deltaDe3 = new TF1("deltaDe3","expo",0.604,1.01);
	grDeuterons_tails -> Fit(deltaDe3,"Q0R");
	TF1 *deltaDe4 = new TF1("deltaDe4",lorentzianPeak,0.28,0.608,3);
	deltaDe4 -> SetParameters(1.4,0.08,0.4);
	grDeuterons_tails -> Fit(deltaDe4,"Q0R");
	TF1 *deltaDe6 = new TF1("deltaDe6","pol1",0.0,0.285);
	grDeuterons_tails -> Fit(deltaDe6,"Q0R");
	

	deltaDe1 -> Draw("same");
	deltaDe2 -> Draw("same");
	deltaDe3 -> Draw("same");
	deltaDe4 -> Draw("same");
	deltaDe6 -> Draw("same");



	grcan3->SaveAs("/home/pepe-frog/Dubna/dEdx/pic/deutrons/fiting/delta.png");
	TCanvas *grcan4 = new TCanvas ("grcan4", "dE/dx width vs P, deuterons", 1200, 800);
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
	
	TF1 *sigmaDe1 = new TF1("sigmaDe1","pol1",1.74,3.06);
	grDeuterons_minsigmas -> Fit(sigmaDe1,"Q0R");
	TF1 *sigmaDe2 = new TF1("sigmaDe2","pol1",0.945,1.76);
	grDeuterons_minsigmas -> Fit(sigmaDe2,"Q0R");
	TF1 *sigmaDe3 = new TF1("sigmaDe3","pol1",0.64,0.955);
	grDeuterons_minsigmas -> Fit(sigmaDe3,"Q0R");
	TF1 *sigmaDe7 = new TF1("sigmaDe7","pol1",0.555,0.642);
	grDeuterons_minsigmas -> Fit(sigmaDe7,"Q0R");
	TF1 *sigmaDe4 = new TF1("sigmaDe4","pol1",0.428,0.557);
	grDeuterons_minsigmas -> Fit(sigmaDe4,"Q0R");
	TF1 *sigmaDe5 = new TF1("sigmaDe5","pol1",0.227,0.428);
	grDeuterons_minsigmas -> Fit(sigmaDe5,"Q0R");
	TF1 *sigmaDe6 = new TF1("sigmaDe6","pol1",0.0,0.228);
	grDeuterons_minsigmas -> Fit(sigmaDe6,"Q0R");
	

	sigmaDe1 -> Draw("same");
	sigmaDe2 -> Draw("same");
	sigmaDe3 -> Draw("same");
	sigmaDe4 -> Draw("same");
	sigmaDe5 -> Draw("same");
	sigmaDe6 -> Draw("same");
	sigmaDe7 -> Draw("same");


	grcan4->SaveAs("/home/pepe-frog/Dubna/dEdx/pic/deutrons/fiting/sigma.png");
	}
	
	
	if (Tparticle == 2) /////////////////////////////////////////////// TRITONS
	{
		TF1 *Gaus[40];
		TF1 *AGaus[40];
		TF1 *Novosib[40];
		TCanvas *can[40];
		Double_t p = 1.6;

		for(Int_t i = 4; i <= 40; i++){
		/*if(i<=31){
		if (!QMode) can[i] = new TCanvas (Form("can%i",i), Form("tritons, %1.2f < p < %1.2f, |eta| < 1.6",i*0.05,i*0.05+0.05), 1200, 800);
		}
		if(i>31){
		if (!QMode) can[i] = new TCanvas (Form("can%i",i), Form("tritons, %1.2f < p < %1.2f, |eta| < 1.6",p,p+0.1), 1200, 800);
		}*/
		Tr[i]->SetStats(kFALSE);
		Tr[i]->GetXaxis()->SetTitle("Energy deposit, ADC counts");
		Tr[i]->GetYaxis()->SetTitle("Counts");
		if (!QMode) Tr[i]->Draw();
		
		Gaus[i] = new TF1(Form("Gaus%i",i), "gaus", Tr[i]->GetXaxis()->GetXmin(), Tr[i]->GetXaxis()->GetXmax());
		Tr[i]->Fit(Gaus[i], "Q0R");

		Novosib[i] = new TF1(Form("Novosib%i",i), Novosibirsk, Tr[i]->GetXaxis()->GetXmin(), Tr[i]->GetXaxis()->GetXmax(), 4);
		Novosib[i]->SetParameters(Gaus[i]->GetParameter(0), 0.01, Gaus[i]->GetParameter(2), Gaus[i]->GetParameter(1)); Novosib[i]->SetLineColor(kBlue);
		Tr[i]->Fit(Form("Novosib%i",i),"Q0RW");
		
		if (!QMode) Novosib[i]->Draw("same");
		AGaus[i] = new TF1(Form("AGaus%i",i), AsymGaus, Tr[i]->GetXaxis()->GetXmin(), Tr[i]->GetXaxis()->GetXmax(), 4);
		AGaus[i]->SetParameters(Gaus[i]->GetParameter(0), Gaus[i]->GetParameter(1), Gaus[i]->GetParameter(2), 0.01);
		Tr[i]->Fit(Form("AGaus%i",i),"Q0RW");
		if (!QMode) AGaus[i]->Draw("same");
/*		if(i<=31){
		if (!QMode) can[i]->SaveAs(Form("/home/pepe-frog/Dubna/dEdx/pic/tritons/tritons_%1.2f<p<%1.2f.png",i*0.05,i*0.05+0.05));
		}
		if(i>31){
		if (!QMode) can[i]->SaveAs(Form("/home/pepe-frog/Dubna/dEdx/pic/tritons/tritons_%1.2f<p<%1.2f.png",p,p+0.1));
		p = p + 0.1;
		}*/
	}
	
	const Int_t nHists = 37;
	const Int_t k = 4;

	Double_t Xlow[nHists] = { 0.20,0.25,0.30,0.35,0.40,0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.4, 2.6 };
	Double_t Xhigh[nHists] = { 0.25,0.30,0.35,0.40,0.45,0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.4, 2.6, 3. };
	Double_t X[nHists], Y[nHists], Y_err[nHists], Y_tail[nHists], X_err[nHists], minsigma[nHists], normerrminsigma[nHists], 
			 normminsigma[nHists], sigma1[nHists], errsigma1[nHists], sigma2[nHists], errsigma2[nHists], Y_tail_err[nHists];

	FillMassive(AGaus, Novosib, nHists, Xlow, Xhigh, X, Y, Y_err, Y_tail, X_err, minsigma, 
		normerrminsigma, normminsigma, sigma1, errsigma1, sigma2, errsigma2, Y_tail_err, k);
	
	TGraphErrors *grTritons = new TGraphErrors(nHists,X,Y,X_err,minsigma);
	
	TCanvas *grcan1 = new TCanvas ("grcan1", "dE/dx vs P, tritons", 1200, 800);
	gPad->SetGrid();
	grTritons->SetTitle("");
	grTritons->GetYaxis()->SetTitle("dE/dx, ADC counts");
	grTritons->GetXaxis()->SetTitle("Momentum, GeV/c");
	grTritons->SetMarkerStyle(21); grTritons->SetMarkerSize(2);
	grTritons->GetXaxis()->SetLabelSize(0.03); grTritons->GetYaxis()->SetLabelSize(0.03);
	grTritons->GetXaxis()->SetTitleSize(0.05); grTritons->GetYaxis()->SetTitleSize(0.05);
	grTritons->GetXaxis()->CenterTitle(); grTritons->GetYaxis()->CenterTitle();
	banana_tr->Draw("colz");
	grTritons->Draw("same");
	
	TF1* fitparTrBB1 = new TF1("fitparTrBB1", "pol1(0)", 0., 0.3);
	grTritons->Fit(fitparTrBB1, "Q0", "", 0.2, 0.3);
	TF1* fitparTrBB2 = new TF1("fitparTrBB2", "pol1(0)", 0.3, 0.55);
	grTritons->Fit(fitparTrBB2, "Q0", "", 0.3, 0.55);
	TF1* fitparTrBB3 = new TF1("fitparTrBB3",fitparTrBB->GetExpFormula(), 0.55, 3.);
	grTritons->Fit(fitparTrBB3, "Q0", "", 0.55, 3.);
	
	TF1* fitparTrBB4 = new TF1("fitparTrBB4", "pol1(0)", 0., 0.325); 
	fitparTrBB4->SetParameters(fitparTrBB1->GetParameters());
	TF1* fitparTrBB5 = new TF1("fitparTrBB5", "pol1(0)", 0.312, 0.547); 
	fitparTrBB5->SetParameters(fitparTrBB2->GetParameters());
	TF1* fitparTrBB6 = new TF1("fitparTrBB6", fitparTrBB->GetExpFormula(), 0.535, 3.);
	 fitparTrBB6->SetParameters(fitparTrBB3->GetParameters());
	fitparTrBB4->Draw("same"); 
	fitparTrBB5->Draw("same"); 
	fitparTrBB6->Draw("same");
	grcan1->SaveAs("/home/pepe-frog/Dubna/dEdx/pic/tritons/fiting/dedx_with_fit.png");

	for (Int_t i = 0; i < fitparTrBB4->GetNpar(); i++) cout << "c[" << i << "] = " << fitparTrBB4->GetParameter(i) << "    ";
	cout << endl;
	for (Int_t i = 0; i < fitparTrBB5->GetNpar(); i++) cout << "c[" << i << "] = " << fitparTrBB5->GetParameter(i) << "    ";
	cout << endl;
	for (Int_t i = 0; i < fitparTrBB6->GetNpar(); i++) cout << "c[" << i << "] = " << fitparTrBB6->GetParameter(i) << "    ";
	cout << endl;
	
	TGraphErrors *grTritons_minsigmas = new TGraphErrors(nHists,X,normminsigma,X_err,normerrminsigma);
	TGraphErrors *grTritons_tails = new TGraphErrors(nHists,X,Y_tail,X_err,Y_tail_err);
	
	TCanvas *grcan3 = new TCanvas ("grcan3", "Tails vs P, tritons", 1200, 800);
	gPad->SetGrid();
	grTritons_tails->SetTitle("");
	grTritons_tails->GetYaxis()->SetTitle("#delta");
	grTritons_tails->GetXaxis()->SetTitle("Momentum, GeV/c");
	grTritons_tails->SetMarkerStyle(21); grTritons_tails->SetMarkerSize(1);
	grTritons_tails->GetXaxis()->SetLabelSize(0.03); grTritons_tails->GetYaxis()->SetLabelSize(0.03);
	grTritons_tails->GetXaxis()->SetTitleSize(0.05); grTritons_tails->GetYaxis()->SetTitleSize(0.05);
	grTritons_tails->GetXaxis()->CenterTitle(); grTritons_tails->GetYaxis()->CenterTitle();
	grTritons_tails->Draw("AP");
	
	TF1 *deltaTr1 = new TF1("deltaTr1","pol2",1.64,3.06);
	grTritons_tails -> Fit(deltaTr1,"Q0R");
	TF1 *deltaTr2 = new TF1("deltaTr2","pol5",0.83,1.68);
	grTritons_tails -> Fit(deltaTr2,"Q0R");
	TF1 *deltaTr3 = new TF1("deltaTr3","expo",0.738,0.84);
	grTritons_tails -> Fit(deltaTr3,"Q0R");
	TF1 *deltaTr4 = new TF1("deltaTr1",lorentzianPeak,0.37,0.745,3);
	deltaTr4 -> SetParameters(1.17963,0.147218,0.57205);
	grTritons_tails->Fit(deltaTr4,"Q0R");
	TF1 *deltaTr6 = new TF1("deltaTr6","pol1",0.315,0.375);
	grTritons_tails -> Fit(deltaTr6,"Q0R");
	TF1 *deltaTr7 = new TF1("deltaTr7","pol1",0.,0.32);
	grTritons_tails -> Fit(deltaTr7,"Q0R");

	deltaTr1 -> Draw("same");
	deltaTr2 -> Draw("same");
	deltaTr3 -> Draw("same");
	deltaTr4 -> Draw("same");
	deltaTr6 -> Draw("same");
	deltaTr7 -> Draw("same");


	
	grcan3->SaveAs("/home/pepe-frog/Dubna/dEdx/pic/tritons/fiting/delta.png");

	TCanvas *grcan4 = new TCanvas ("grcan4", "dE/dx width vs P, tritons", 1200, 800);
	gPad->SetGrid();
	grTritons_minsigmas->SetTitle("");
	grTritons_minsigmas->GetYaxis()->SetTitle("#sigma");
	grTritons_minsigmas->GetXaxis()->SetTitle("Momentum, GeV/c");
	grTritons_minsigmas->SetMarkerStyle(21); grTritons_minsigmas->SetMarkerSize(1);
	grTritons_minsigmas->GetXaxis()->SetLabelSize(0.03); grTritons_minsigmas->GetYaxis()->SetLabelSize(0.03);
	grTritons_minsigmas->GetXaxis()->SetTitleSize(0.05); grTritons_minsigmas->GetYaxis()->SetTitleSize(0.05);
	grTritons_minsigmas->GetXaxis()->CenterTitle(); grTritons_minsigmas->GetYaxis()->CenterTitle();
	grTritons_minsigmas->Draw("AP");
	
	TF1 *sigmaTr1 = new TF1("sigmaTr1","pol1",2.05,3.06);
	grTritons_minsigmas -> Fit(sigmaTr1,"Q0R");
	TF1 *sigmaTr2 = new TF1("sigmaTr2","pol1",1.65,2.07);
	grTritons_minsigmas -> Fit(sigmaTr2,"Q0R");
	TF1 *sigmaTr3 = new TF1("sigmaTr3","pol2",1.143,1.658);
	grTritons_minsigmas -> Fit(sigmaTr3,"Q0R");
	TF1 *sigmaTr7 = new TF1("sigmaTr7","pol2",0.845,1.15);
	grTritons_minsigmas -> Fit(sigmaTr7,"Q0R");
	TF1 *sigmaTr4 = new TF1("sigmaTr4","pol1",0.5795,0.855);
	grTritons_minsigmas -> Fit(sigmaTr4,"Q0R");
	TF1 *sigmaTr5 = new TF1("sigmaTr5","pol1",0.385,0.583);
	grTritons_minsigmas -> Fit(sigmaTr5,"Q0R");
	TF1 *sigmaTr6 = new TF1("sigmaTr6","pol1",0.318,0.389);
	grTritons_minsigmas -> Fit(sigmaTr6,"Q0R");
	TF1 *sigmaTr8 = new TF1("sigmaTr8","pol1",0.22,0.32);
	grTritons_minsigmas -> Fit(sigmaTr8,"Q0R");
	TF1 *sigmaTr9 = new TF1("sigmaTr9","pol0",0.,0.235);
	grTritons_minsigmas -> Fit(sigmaTr9,"Q0R");

	sigmaTr1 -> Draw("same");
	sigmaTr2 -> Draw("same");
	sigmaTr3 -> Draw("same");
	sigmaTr4 -> Draw("same");
	sigmaTr5 -> Draw("same");
	sigmaTr6 -> Draw("same");
	sigmaTr7 -> Draw("same");
	sigmaTr8 -> Draw("same");
	sigmaTr9 -> Draw("same");

	grcan4->SaveAs("/home/pepe-frog/Dubna/dEdx/pic/tritons/fiting/sigma.png");

}
	
	if (Tparticle == 3) //////////////////////////////////////// He^3
	{
		TF1 *Gaus[33];
		TF1 *AGaus[33];
		TF1 *Novosib[33];
		TCanvas *can[33];
		Double_t p =1.6;

		for(Int_t i = 3; i <= 33; i++){
/*		if(i<=31){
		if (!QMode) can[i] = new TCanvas (Form("can%i",i), Form("he3, %1.2f < p < %1.2f, |eta| < 1.6",i*0.05,i*0.05+0.05), 1200, 800);
		}
		if(i>31){
		if (!QMode) can[i] = new TCanvas (Form("can%i",i), Form("he3, %1.2f < p < %1.2f, |eta| < 1.6",p,p+0.1), 1200, 800);
		}*/
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
	
	
	const Int_t nHists = 31;
	const Int_t k = 3;

	Double_t Xlow[nHists] = { 0.15,0.20,0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.7,};
	Double_t Xhigh[nHists] = { 0.20,0.25,0.30,0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.7, 1.8};
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
	banana_he3->Draw("colz");
	grHe3->Draw("same");
	
	TF1* fitparHe3BB1 = new TF1("fitparHe3BB1", "pol1(0)", 0., 0.35);
	grHe3->Fit(fitparHe3BB1, "Q0", "", 0.15, 0.35);
	TF1* fitparHe3BB2 = new TF1("fitparHe3BB2", "pol1(0)", 0.35, 0.55);
	grHe3->Fit(fitparHe3BB2, "Q0", "", 0.35, 0.55);
	TF1* fitparHe3BB3 = new TF1("fitparHe3BB3",fitparHe3BB->GetExpFormula(), 0.55, 1.8);
	grHe3->Fit(fitparHe3BB3, "Q0", "", 0.55, 1.8);
	
	TF1* fitparHe3BB4 = new TF1("fitparHe3BB4", "pol1(0)", 0., 0.339); fitparHe3BB4->SetParameters(fitparHe3BB1->GetParameters());
	TF1* fitparHe3BB5 = new TF1("fitparHe3BB5", "pol1(0)", 0.338, 0.555); fitparHe3BB5->SetParameters(fitparHe3BB2->GetParameters());
	TF1* fitparHe3BB6 = new TF1("fitparHe3BB6", fitparHe3BB->GetExpFormula(), 0.545, 3.0); fitparHe3BB6->SetParameters(fitparHe3BB3->GetParameters());
	
	fitparHe3BB4->Draw("same"); 
	fitparHe3BB5->Draw("same"); 
	fitparHe3BB6->Draw("same");

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

	TF1 *deltaHe3_1 = new TF1("deltaHe3_1","pol0",0.02,0.175);
	grHe3_tails -> Fit(deltaHe3_1,"Q0R");
	TF1 *deltaHe3_2 = new TF1("deltaHe3_2","pol1",0.1745,0.227);
	grHe3_tails -> Fit(deltaHe3_2,"Q0R");
	TF1 *deltaHe3_3 = new TF1("deltaHe3_3","pol4",0.225,0.43);
	grHe3_tails -> Fit(deltaHe3_3,"Q0R");
	TF1 *deltaHe3_4_1 = new TF1("deltaHe3_4_1","pol1",0.425,0.625);
	grHe3_tails->Fit(deltaHe3_4_1,"Q0R");
	TF1 *deltaHe3_4_2 = new TF1("deltaHe3_4_2","pol1",0.623,0.717);
	grHe3_tails->Fit(deltaHe3_4_2,"Q0R");
	TF1 *deltaHe3_4_3 = new TF1("deltaHe3_4_3","pol1",0.716,0.902);
	grHe3_tails->Fit(deltaHe3_4_3,"Q0R");
	TF1 *deltaHe3_4_4 = new TF1("deltaHe3_4_4","pol1",0.9,1.0);
	grHe3_tails->Fit(deltaHe3_4_4,"Q0R");
	TF1 *deltaHe3_6 = new TF1("deltaHe3_6","pol2",0.99,1.28);
	grHe3_tails -> Fit(deltaHe3_6,"Q0R");
	TF1 *deltaHe3_7 = new TF1("deltaHe3_7","pol1",1.273,1.4199);
	grHe3_tails -> Fit(deltaHe3_7,"Q0R");
	TF1 *deltaHe3_8 = new TF1("deltaHe3_8","pol1",1.405,1.51);
	grHe3_tails -> Fit(deltaHe3_8,"Q0R");
	TF1 *deltaHe3_9 = new TF1("deltaHe3_9","pol1",1.525,1.576);
	grHe3_tails -> Fit(deltaHe3_9,"Q0R");
	TF1 *deltaHe3_10 = new TF1("deltaHe3_10","pol1",1.5745,1.65);
	grHe3_tails -> Fit(deltaHe3_10,"Q0R");
	TF1 *deltaHe3_11 = new TF1("deltaHe3_11","pol1",1.649,1.7505);
	grHe3_tails -> Fit(deltaHe3_11,"Q0R");
	TF1 *deltaHe3_12 = new TF1("deltaHe3_12","pol0",1.745,1.91);
	grHe3_tails -> Fit(deltaHe3_12,"Q0R");
	TF1 *deltaHe3_13 = new TF1("deltaHe3_13","pol1",1.418,1.5265);
	grHe3_tails -> Fit(deltaHe3_13,"Q0R");


	deltaHe3_1 -> Draw("same");
	deltaHe3_2 -> Draw("same");
	deltaHe3_3 -> Draw("same");
	deltaHe3_4_1 -> Draw("same");
	deltaHe3_4_2 -> Draw("same");
	deltaHe3_4_3 -> Draw("same");
	deltaHe3_4_4 -> Draw("same");
	deltaHe3_6 -> Draw("same");
	deltaHe3_7 -> Draw("same");
	deltaHe3_9 -> Draw("same");
	deltaHe3_10 -> Draw("same");
	deltaHe3_11 -> Draw("same");
	deltaHe3_12 -> Draw("same");
	deltaHe3_13 -> Draw("same");
	
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

	TF1 *sigmaHe3_1 = new TF1("sigmaHe3_1","pol0",0.02,0.175);
	grHe3_minsigmas -> Fit(sigmaHe3_1,"Q0R");
	TF1 *sigmaHe3_2 = new TF1("sigmaHe3_2","pol4",0.174,0.4255);
	grHe3_minsigmas -> Fit(sigmaHe3_2,"Q0R");
	TF1 *sigmaHe3_3 = new TF1("sigmaHe3_3","pol3",0.415,0.725);
	grHe3_minsigmas -> Fit(sigmaHe3_3,"Q0R");
	TF1 *sigmaHe3_7 = new TF1("sigmaHe3_7","pol3",0.72,1.239);
	grHe3_minsigmas -> Fit(sigmaHe3_7,"Q0R");
	TF1 *sigmaHe3_4 = new TF1("sigmaHe3_4","pol1",1.237,1.327);
	grHe3_minsigmas -> Fit(sigmaHe3_4,"Q0R");
	TF1 *sigmaHe3_5 = new TF1("sigmaHe3_5","pol1",1.323,1.375);
	grHe3_minsigmas -> Fit(sigmaHe3_5,"Q0R");
	TF1 *sigmaHe3_6 = new TF1("sigmaHe3_6","pol1",1.374,1.647);
	grHe3_minsigmas -> Fit(sigmaHe3_6,"Q0R");
	TF1 *sigmaHe3_10 = new TF1("sigmaHe3_10","pol0",1.645,1.91);
	grHe3_minsigmas -> Fit(sigmaHe3_10,"Q0R");

	sigmaHe3_1 -> Draw("same");
	sigmaHe3_2 -> Draw("same");
	sigmaHe3_3 -> Draw("same");
	sigmaHe3_4 -> Draw("same");
	sigmaHe3_5 -> Draw("same");
	sigmaHe3_6 -> Draw("same");
	sigmaHe3_7 -> Draw("same");
	sigmaHe3_10 -> Draw("same");
	
	
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
/*		if(i<=31){
		if (!QMode) can[i] = new TCanvas (Form("can%i",i), Form("he4, %1.2f < p < %1.2f, |eta| < 1.6",i*0.05,i*0.05+0.05), 1200, 800);
		}
		if(i>31){
		if (!QMode) can[i] = new TCanvas (Form("can%i",i), Form("he4, %1.2f < p < %1.2f, |eta| < 1.6",p,p+0.1), 1200, 800);
		}*/
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
	
	const Int_t nHists = 28;
	const Int_t k = 4;

	Double_t Xlow[nHists] = { 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55};
	Double_t Xhigh[nHists] = { 0.25, 0.3, 0.35,0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6};
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
	banana_he4->Draw("colz");
	grHe4->Draw("same");
	
	TF1* fitparHe4BB1 = new TF1("fitparHe4BB1", "pol1(0)", 0., 0.35);
	grHe4->Fit(fitparHe4BB1, "Q0", "", 0.2, 0.35);
	TF1* fitparHe4BB2 = new TF1("fitparHe4BB2", "pol1(0)", 0.35, 0.5);
	grHe4->Fit(fitparHe4BB2, "Q0", "", 0.35, 0.529);
	TF1* fitparHe4BB3 = new TF1("fitparHe4BB3",fitparHe4BB->GetExpFormula(), 0.45, 1.6);
	grHe4->Fit(fitparHe4BB3, "Q0", "", 0.45, 1.6);
	
	TF1* fitparHe4BB4 = new TF1("fitparHe4BB4", "pol1(0)", 0., 0.36); 
	fitparHe4BB4->SetParameters(fitparHe4BB1->GetParameters());
	TF1* fitparHe4BB5 = new TF1("fitparHe4BB5", "pol1(0)", 0.357, 0.529);
	 fitparHe4BB5->SetParameters(fitparHe4BB2->GetParameters());
	TF1* fitparHe4BB6 = new TF1("fitparHe4BB6", fitparHe4BB->GetExpFormula(), 0.515, 3.0); 
	fitparHe4BB6->SetParameters(fitparHe4BB3->GetParameters());
	fitparHe4BB4->Draw("same");
	 fitparHe4BB5->Draw("same");
	  fitparHe4BB6->Draw("same");
		
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

	TF1 *deltaHe4_1 = new TF1("deltaHe4_1","pol0",0.0,0.23);
	grHe4_tails -> Fit(deltaHe4_1,"Q0R");
	TF1 *deltaHe4_2 = new TF1("deltaHe4_2","pol3",0.225,0.383);
	grHe4_tails -> Fit(deltaHe4_2,"Q0R");
	TF1 *deltaHe4_3 = new TF1("deltaHe4_3","pol1",0.3795,0.53);
	grHe4_tails -> Fit(deltaHe4_3,"Q0R");
	TF1 *deltaHe4_4 = new TF1("deltaHe4_4","pol3",0.52,0.945);
	grHe4_tails->Fit(deltaHe4_4,"Q0R");
	TF1 *deltaHe4_6 = new TF1("deltaHe4_6","pol1",0.94,1.221);
	grHe4_tails -> Fit(deltaHe4_6,"Q0R");
	TF1 *deltaHe4_7 = new TF1("deltaHe4_7","pol1",1.21,1.221);
	grHe4_tails -> Fit(deltaHe4_7,"Q0R");
	TF1 *deltaHe4_8_1 = new TF1("deltaHe4_8_1","pol1",1.22,1.275);
	grHe4_tails -> Fit(deltaHe4_8_1,"Q0R");
	TF1 *deltaHe4_8_2 = new TF1("deltaHe4_8_2","pol1",1.2745,1.326);
	grHe4_tails -> Fit(deltaHe4_8_2,"Q0R");
	TF1 *deltaHe4_8_3 = new TF1("deltaHe4_8_3","pol1",1.325,1.3775);
	grHe4_tails -> Fit(deltaHe4_8_3,"Q0R");
	TF1 *deltaHe4_9 = new TF1("deltaHe4_9","pol1",1.374,1.477);
	grHe4_tails -> Fit(deltaHe4_9,"Q0R");
	TF1 *deltaHe4_10 = new TF1("deltaHe4_10","pol1",1.476,1.579);
	grHe4_tails -> Fit(deltaHe4_10,"Q0R");
	TF1 *deltaHe4_11 = new TF1("deltaHe4_11","pol0",1.575,1.71);
	grHe4_tails -> Fit(deltaHe4_11,"Q0R");

	deltaHe4_1 -> Draw("same");
	deltaHe4_2 -> Draw("same");
	deltaHe4_3 -> Draw("same");
	deltaHe4_4 -> Draw("same");
	deltaHe4_6 -> Draw("same");
	deltaHe4_7 -> Draw("same");
	deltaHe4_8_1 -> Draw("same");
	deltaHe4_8_2 -> Draw("same");
	deltaHe4_8_3 -> Draw("same");
	deltaHe4_9 -> Draw("same");
	deltaHe4_10 -> Draw("same");
	deltaHe4_11 -> Draw("same");

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
	
	TF1 *sigmaHe4_1 = new TF1("sigmaHe4_1","pol0",0.02,0.225);
	grHe4_minsigmas -> Fit(sigmaHe4_1,"Q0R");
	TF1 *sigmaHe4_2 = new TF1("sigmaHe4_2","pol4",0.224,0.428);
	grHe4_minsigmas -> Fit(sigmaHe4_2,"Q0R");
	TF1 *sigmaHe4_3_1 = new TF1("sigmaHe4_3_1","pol1",0.424,0.475);
	grHe4_minsigmas -> Fit(sigmaHe4_3_1,"Q0R");
	TF1 *sigmaHe4_3_2 = new TF1("sigmaHe4_3_2","pol3",0.472,0.768);
	grHe4_minsigmas -> Fit(sigmaHe4_3_2,"Q0R");
	TF1 *sigmaHe4_7 = new TF1("sigmaHe4_7","pol3",0.7635,1.217);
	grHe4_minsigmas -> Fit(sigmaHe4_7,"Q0R");
	TF1 *sigmaHe4_4 = new TF1("sigmaHe4_4","pol1",1.237,1.327);
	grHe4_minsigmas -> Fit(sigmaHe4_4,"Q0R");
	TF1 *sigmaHe4_5 = new TF1("sigmaHe4_5","pol1",1.214,1.2765);
	grHe4_minsigmas -> Fit(sigmaHe4_5,"Q0R");
	TF1 *sigmaHe4_6 = new TF1("sigmaHe4_6","pol1",1.275,1.327);
	grHe4_minsigmas -> Fit(sigmaHe4_6,"Q0R");
	TF1 *sigmaHe4_6_1 = new TF1("sigmaHe4_6_1","pol1",1.32,1.425);
	grHe4_minsigmas -> Fit(sigmaHe4_6_1,"Q0R");
	TF1 *sigmaHe4_6_2 = new TF1("sigmaHe4_6_2","pol1",1.42,1.475);
	grHe4_minsigmas -> Fit(sigmaHe4_6_2,"Q0R");
	TF1 *sigmaHe4_6_3 = new TF1("sigmaHe4_6_3","pol1",1.474,1.525);
	grHe4_minsigmas -> Fit(sigmaHe4_6_3,"Q0R");
	TF1 *sigmaHe4_10 = new TF1("sigmaHe4_10","pol0",1.523,1.71);
	grHe4_minsigmas -> Fit(sigmaHe4_10,"Q0R");

	sigmaHe4_1 -> Draw("same");
	sigmaHe4_2 -> Draw("same");
	sigmaHe4_3_1 -> Draw("same");
	sigmaHe4_3_2 -> Draw("same");
	sigmaHe4_5 -> Draw("same");
	sigmaHe4_6 -> Draw("same");
	sigmaHe4_6_1 -> Draw("same");
	sigmaHe4_6_2 -> Draw("same");
	sigmaHe4_6_3 -> Draw("same");
	sigmaHe4_7 -> Draw("same");
	sigmaHe4_10 -> Draw("same");

		
	grcan4->SaveAs("/home/pepe-frog/Dubna/dEdx/pic/he4/fiting/sigma.png");

	}
}
