#include "TF1.h"
#include "TLegend.h"

TF1* initPoissonDist(double mu, double N) {
    auto f = new TF1("f", "TMath::PoissonI(x, [0])", 0, N);
    f->SetParameter(0, mu);
    f->SetTitle(" ");
    f->SetLineColor(kRed);
    f->SetLineWidth(2);
    f->SetNpx(100);
    return f;
}

TF1* initBinomialDist(double p, unsigned int N) {
    auto f = new TF1("f", "ROOT::Math::binomial_pdf(x, [0], [1])", 0, N);
    f->SetParameters(p, N);
    f->SetLineColor(kBlue);
    f->SetTitle(" ");
    f->SetLineWidth(2);
    f->SetNpx(100);
    return f;
}

void task1_3() {
    unsigned int N = 100;
    double p = 0.1;

    for (int i = 0; i < 9; i++) {
        p = 0.1 * (i + 1);
        auto f_p = initPoissonDist(N*p, N);
        auto f_b = initBinomialDist(p, N);

        if (i == 0) {
            f_b->Draw("HIST");
            f_p->Draw("HIST same");
        }
        if (i % 2 != 0) {
            f_b->Draw("HIST same");
            f_p->Draw("HIST same");
        }
    }

    auto f_p = initPoissonDist(N*p, N); // для отрисовки легенды
    auto f_b = initBinomialDist(p, N);
    auto legend = new TLegend(0.4, 0.7, 0.6, 0.8);
    legend->AddEntry(f_b, "Binomial");
    legend->AddEntry(f_p, "Poisson");
    legend->Draw();
}
