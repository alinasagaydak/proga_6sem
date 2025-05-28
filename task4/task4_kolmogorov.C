#include <iostream>
#include "TCanvas.h"
#include "TH1.h"
#include "TFitResultPtr.h"
#include "TLine.h"
#include "TF1.h"

void task4_kolmogorov() {
    double N = 1000; // кол-во ген событий
    double Nbins = 100;

    auto c1 = new TCanvas("c1", "", 1200, 900);
    c1->Divide(2, 1);

    auto h = new TH1D("h", "h", Nbins, -2., 2.); // гистограмма с данными
    auto h2 = new TH1D("h", "h", Nbins, -2., 2.); // гистограмма с данными
    auto fitfunc = new TF1("fitfunc", "[0] + [1]*x + [2]*x*x", -2., 2.);

    auto h4 = new TH1D("h4", "", Nbins, -2., 2.);
    auto htest_true = new TH1D("htest_true", "", Nbins, 0., 1.);
    auto htest_false = new TH1D("htest_false", "", Nbins, 0., 1.);
    h4->FillRandom("gaus", N * 100);
    double p;

    TFitResultPtr r;

    c1->cd(1);
    h->FillRandom("gaus", N);
    h->Draw();
    r = h->Fit(fitfunc, "LS", "ep");\

    c1->cd(2);
    htest_false->Reset();
    htest_true->Reset();
    
    for(int i = 0; i < N * 10; i++) {
        h2->Reset();
        h2->FillRandom("gaus", N);
        p = h2->KolmogorovTest(h4);
        htest_true->Fill(p);
    }

    for(int i = 0; i < N * 10; i++) {
        h2->Reset();
        h2->FillRandom("fitfunc", N);
        p = h2->KolmogorovTest(h4);
        htest_false->Fill(p);
    }

    htest_false->SetLineColor(kRed);
    htest_false->Draw();
    htest_true->Draw("same");
    c1->Draw();

    std::cout << "p_false = " << htest_false->Integral(6, 100)/(N * 10.) << std::endl;
}
