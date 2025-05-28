#include <iostream>
#include "TCanvas.h"
#include "TH1.h"
#include "TFitResultPtr.h"
#include "TLine.h"
#include "TF1.h"

void task4() {
    double N = 1000; // кол-во ген событий
    double Nbins = 100;

    auto c1 = new TCanvas("c1", "", 1200, 900);
    c1->Divide(2, 1);

    auto h = new TH1D("h", "h", Nbins, -2., 2.); // гистограмма с данными
    auto h1 = new TH1D("h1", "h", Nbins, -2., 2.); // гистограмма с данными
    auto h2 = new TH1D("h2", "", 50, 0., 50.); // гистограмма с тестовой статистикой
    TFitResultPtr r;

    c1->cd(1);
    h->FillRandom("gaus", N);
    h->Draw();
    r = h->Fit("pol2", "LS", "ep");

    c1->cd(2);
    auto fchi2 = new TF1("fchi2", "[0] * ROOT::Math::chisquared_pdf(x, [1])", 0., 300.);
    fchi2->SetParameters(N*10, Nbins-3);
    TLine line;
    h2->Reset();
    for(int i = 0; i < N * 10; i++) {
        h1->Reset();
        h1->FillRandom("gaus", N);
        r = h1->Fit("pol2", "LSQN");
        h2->Fill(r->Chi2());
    }

    h2->Draw();
    fchi2->SetLineColor(kRed);
    fchi2->Draw("same");
    line.DrawLine(TMath::ChisquareQuantile(0.95, Nbins-3), 0, TMath::ChisquareQuantile(0.95, Nbins-3), 800);

    c1->Draw();

    std::cout << "p = " << r->Prob() << "\tchi2 = " << r->Chi2() << "\tndf = " << r->Ndf() << std::endl;
    std::cout << "chi2_max = " << TMath::ChisquareQuantile(0.95, Nbins-3) << std::endl;
}
