#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"

TF1* initGausDist(double N, double min, double max) {
    auto f = new TF1("f", "TMath::Gaus(x, [0], sqrt(2 * [0]), kTRUE)", min, max);
    f->SetLineColor(kRed);
    f->SetTitle(" ; ;");
    f->SetParameter(0, N);
    f->SetNpx(10000);
    return f;
}

TF1* initChiDist(double ndf, double min, double max) {
    auto f = new TF1("f", "ROOT::Math::gamma_pdf(x, [0], [1])", min, max);
    f->SetLineColor(kBlue);
    f->SetTitle(" ; ;");
    f->SetParameter(0, ndf / 2.);
    f->SetParameter(1, 2);
    f->SetNpx(10000);
    return f;

}

void task1_4() {
    double N = 3.;
    double min = 0.;
    double max = 160.;

    for (int i = 0; i < 10; i++) {
        N = 15 * i + 3;
        if (i == 0) {
            auto f_gaus = initGausDist(N, min, max);
            auto f_chi = initChiDist(N, min, max);
            f_gaus->Draw();
            f_chi->Draw();
        }
        if (i % 2 == 0) {
            auto f_gaus = initGausDist(N, min, max);
            auto f_chi = initChiDist(N, min, max);
            f_gaus->Draw("same");
            f_chi->Draw("same");
        } 
    }

    auto f_gaus = initGausDist(N, min, max); // для отрисовки легенды
    auto f_chi = initChiDist(N, min, max);
    auto legend = new TLegend(0.6, 0.7, 0.8, 0.8);
    legend->AddEntry(f_gaus, "Gaussian");
    legend->AddEntry(f_chi, "#chi^{2} distribution");
    legend->Draw();
}
