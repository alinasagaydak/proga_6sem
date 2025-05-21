#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>

void add_quest() {
    int N = 1000;
    auto c1 = new TCanvas("c1", "", 1200, 800);

    auto h1 = new TH1D("h1", Form("N = %i", N), 40, -100., 100.);
    auto h_d = new TH1D("h_d", "", 300, 0., 50.);

    TF1* func = new TF1("func", "[0]+[1]*(x/100.)+[2]*(x/100.)**2+[3]*(x/100.)**3+[4]*(x/100.)**4", -100, 100);
    func->SetParameters(10,-6,3,7,-3);

    auto fitfunc5 = new TF1("fitfunc5", "[0] + [1]*(x/100.) + [2]*(x/100.)**2 + [3]*(x/100.)**3 + [4]*(x/100.)**4 + [5]*(x/100.)**5", -100., 100.);
    auto fitfunc9 = new TF1("fitfunc9",
            "[0] + [1]*(x/100.) + [2]*(x/100.)**2 + [3]*(x/100.)**3 + [4]*(x/100.)**4 + [5]*(x/100.)**5 + [6]*(x/100.)**6 + [7]*(x/100.)**7 + [8]*(x/100.)**8 + [9]*(x/100.)**9",
            -100., 100.);

    TFitResultPtr F5;
    TFitResultPtr F9;

    h_d->Reset(); 
    for (int i = 0; i < N; i++) {
        h1->Reset();
        h1->FillRandom("func", N);
        F5 = h1->Fit(fitfunc5, "LSQ0");
        F9 = h1->Fit(fitfunc9, "LSQ0");
        double d = -2 * (F9->MinFcnValue() - F5->MinFcnValue());
        h_d->Fill(d);
    }

    TF1* chifunc = new TF1("chifunc", "[0] * ROOT::Math::chisquared_pdf(x, 4)", 0., 50.);
    chifunc->SetParameter(0, N*0.25);
    h_d->Draw();
    h_d->Fit(chifunc, "L");
}
