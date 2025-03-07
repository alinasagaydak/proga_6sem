#include "TLegend.h"

TF1* initGausDist(double N, double min, double max) {
    auto f = new TF1("f", "TMath::Gaus(x, [0], sqrt([0]), kTRUE)", min, max);
    f->SetLineColor(kRed);
    f->SetTitle(" ; ;");
    f->SetParameter(0, N);
    f->SetNpx(10000); 
    return f;
}

TF1* initPoissonDist(double N, double min, double max) {
    auto f = new TF1("f", "TMath::Poisson(x, [0])", min, max);
    f->SetLineColor(kBlue);
    f->SetTitle(" ");
    f->SetParameter(0, N);
    f->SetNpx(10000); 
    return f;
}

void task1_2() {
    double min = 0.;
    double max = 120.;
    double N = 3.;

    for (int i = 0; i < 10; i++) {
        N = 10 * i + 3.;
        auto f_gaus = initGausDist(N, min, max);
        auto f_pois = initPoissonDist(N, min, max);
        
        if (i == 0) {
            f_pois->Draw();
            f_gaus->Draw("same");
        }

        if (i % 2 != 0) {
            f_pois->Draw("same");
            f_gaus->Draw("same");
        }
    }
     
    auto f_gaus = initGausDist(N, min, max); // для отрисовки легенды
    auto f_pois = initPoissonDist(N, min, max);
    auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(f_gaus, "Gaus");
    legend->AddEntry(f_pois, "Poisson");
    legend->Draw();
}
