#include <iostream>
#include "TCanvas.h"
#include "TH1.h"
#include "TFitResultPtr.h"
#include "TLine.h"
#include "TF1.h"

void task4_lntoymc() {
    double N = 1000; // кол-во ген событий
    double Nbins = 100;

    auto c1 = new TCanvas("c1", "", 1200, 900);
    c1->Divide(2, 1);

    auto h = new TH1D("h", "h", Nbins, -2., 2.); // гистограмма с данными
    auto h2 = new TH1D("h", "h", Nbins, -2., 2.); // гистограмма с данными
    auto hfcn_true = new TH1D("hfcn_true", "", 100, 30, 130);
    auto hfcn_false = new TH1D("hfcn_false", "", 100, 30, 130);
    //auto gausfunc = new TF1("gausfunc", "[0] * TMath::Gaus(x, 0, 1, 1)", -2., 2.);
    auto fitfunc = new TF1("fitfunc", "[0] + [1]*x + [2]*x*x", -2., 2.);
    
    TFitResultPtr rold;
    TFitResultPtr r;
    int j;

    c1->cd(1);
    h->FillRandom("gaus", N);
    h->Draw();
    rold = h->Fit(fitfunc, "LS", "ep");\
    double fcn = rold->MinFcnValue();

    c1->cd(2);
    hfcn_true->Reset();
    hfcn_false->Reset();
    for(int i = 0; i < N * 10; i++) {
        h2->Reset();
        h2->FillRandom("gaus", gRandom->Poisson(N));
        r = h2->Fit("pol2", "LSQN");
        hfcn_false->Fill(r->MinFcnValue());
    }

    int nfcn = 0;
    for(int i = 0; i < N * 10; i++) {
        h2->Reset();
        h2->FillRandom("fitfunc", gRandom->Poisson(N));
        r = h2->Fit("pol2", "LSQN");
        hfcn_true->Fill(r->MinFcnValue());
        if(r->MinFcnValue() > fcn) nfcn++;
    }

    for(j = 50; j > 0; j--) {
     if(hfcn_true->Integral(j, 100) > 500) break;
    }

    hfcn_false->SetLineColor(kRed);
    hfcn_true->Draw();
    hfcn_false->Draw("same");
    c1->Draw();
    std::cout << "p = " << nfcn/10000. << "\tfcn = " << fcn << "\tfcn_max = " << hfcn_true->GetBinCenter(j)
        << "\tp_false = " << hfcn_false->Integral(1, j)/10000. << std::endl;


}
