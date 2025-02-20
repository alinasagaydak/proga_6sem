#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>

void task2() {
    int N = 100000;
    auto c1 = new TCanvas("c1", "", 1200, 800);
    auto h1 = new TH1D("h1", Form("N = %i", N), 40, -100., 100.);
    TF1* func = new TF1("func", "[0]+[1]*(x/100.)+[2]*(x/100.)**2+[3]*(x/100.)**3+[4]*(x/100.)**4", -100, 100);
    func->SetParameters(10,-6,3,7,-3);
    h1->FillRandom("func", N);
    h1->SetLineWidth(2);
    h1->Draw("e");

    auto fitfunc0 = new TF1("fitfunc0", "[0]", -100., 100.);
    auto fitfunc1 = new TF1("fitfunc1", "[0] + [1]*(x/100.)", -100., 100.);
    auto fitfunc2 = new TF1("fitfunc2", "[0] + [1]*(x/100.) + [2]*(x/100.)**2", -100., 100.);
    auto fitfunc3 = new TF1("fitfunc3", "[0] + [1]*(x/100.) + [2]*(x/100.)**2 + [3]*(x/100.)**3", -100., 100.);
    auto fitfunc4 = new TF1("fitfunc4", "[0] + [1]*(x/100.) + [2]*(x/100.)**2 + [3]*(x/100.)**3 + [4]*(x/100.)**4", -100., 100.);
    auto fitfunc5 = new TF1("fitfunc5", "[0] + [1]*(x/100.) + [2]*(x/100.)**2 + [3]*(x/100.)**3 + [4]*(x/100.)**4 + [5]*(x/100.)**5", -100., 100.);
    auto fitfunc6 = new TF1("fitfunc6", "[0] + [1]*(x/100.) + [2]*(x/100.)**2 + [3]*(x/100.)**3 + [4]*(x/100.)**4 + [5]*(x/100.)**5 + [6]*(x/100.)**6", -100., 100.);
    auto fitfunc7 = new TF1("fitfunc7", "[0] + [1]*(x/100.) + [2]*(x/100.)**2 + [3]*(x/100.)**3 + [4]*(x/100.)**4 + [5]*(x/100.)**5 + [6]*(x/100.)**6 + [7]*(x/100.)**7", -100., 100.);
    auto fitfunc8 = new TF1("fitfunc8",
            "[0] + [1]*(x/100.) + [2]*(x/100.)**2 + [3]*(x/100.)**3 + [4]*(x/100.)**4 + [5]*(x/100.)**5 + [6]*(x/100.)**6 + [7]*(x/100.)**7 + [8]*(x/100.)**8", -100., 100.);
    auto fitfunc9 = new TF1("fitfunc9",
            "[0] + [1]*(x/100.) + [2]*(x/100.)**2 + [3]*(x/100.)**3 + [4]*(x/100.)**4 + [5]*(x/100.)**5 + [6]*(x/100.)**6 + [7]*(x/100.)**7 + [8]*(x/100.)**8 + [9]*(x/100.)**9",
            -100., 100.);

    fitfunc0->SetLineColor(2); fitfunc0->SetLineWidth(3);
    fitfunc1->SetLineColor(3); fitfunc1->SetLineWidth(3);
    fitfunc2->SetLineColor(4); fitfunc2->SetLineWidth(3);
    fitfunc3->SetLineColor(5); fitfunc3->SetLineWidth(3);
    fitfunc4->SetLineColor(6); fitfunc4->SetLineWidth(3);
    fitfunc5->SetLineColor(7); fitfunc5->SetLineWidth(3);
    fitfunc6->SetLineColor(8); fitfunc6->SetLineWidth(3);
    fitfunc7->SetLineColor(9); fitfunc7->SetLineWidth(3);
    fitfunc8->SetLineColor(41); fitfunc8->SetLineWidth(3);
    fitfunc9->SetLineColor(46); fitfunc9->SetLineWidth(3);

    TFitResultPtr F0 = h1->Fit(fitfunc0, "LS+", "same");
    TFitResultPtr F1 = h1->Fit(fitfunc1, "LS+", "same");
    TFitResultPtr F2 = h1->Fit(fitfunc2, "LS+", "same");
    TFitResultPtr F3 = h1->Fit(fitfunc3, "LS+", "same");
    TFitResultPtr F4 = h1->Fit(fitfunc4, "LS+", "same");
    TFitResultPtr F5 = h1->Fit(fitfunc5, "LS+", "same");
    TFitResultPtr F6 = h1->Fit(fitfunc6, "LS+", "same");
    TFitResultPtr F7 = h1->Fit(fitfunc7, "LS+", "same");
    TFitResultPtr F8 = h1->Fit(fitfunc8, "LS+", "same");
    TFitResultPtr F9 = h1->Fit(fitfunc9, "LS+", "same");

    auto legend = new TLegend(0.1, 0.1, 0.2, 0.3);
    legend->AddEntry(fitfunc0, "pol0");
    legend->AddEntry(fitfunc1, "pol1");
    legend->AddEntry(fitfunc2, "pol2");
    legend->AddEntry(fitfunc3, "pol3");
    legend->AddEntry(fitfunc4, "pol4");
    legend->AddEntry(fitfunc5, "pol5");
    legend->AddEntry(fitfunc6, "pol6");
    legend->AddEntry(fitfunc7, "pol7");
    legend->AddEntry(fitfunc8, "pol8");
    legend->AddEntry(fitfunc9, "pol9");
    legend->Draw();

    double D01 = -2 * (F1->MinFcnValue() - F0->MinFcnValue());
    double D12 = -2 * (F2->MinFcnValue() - F1->MinFcnValue());
    double D23 = -2 * (F3->MinFcnValue() - F2->MinFcnValue());
    double D34 = -2 * (F4->MinFcnValue() - F3->MinFcnValue());
    double D45 = -2 * (F5->MinFcnValue() - F4->MinFcnValue());
    double D56 = -2 * (F6->MinFcnValue() - F5->MinFcnValue());
    double D67 = -2 * (F7->MinFcnValue() - F6->MinFcnValue());
    double D78 = -2 * (F8->MinFcnValue() - F7->MinFcnValue());
    double D89 = -2 * (F9->MinFcnValue() - F8->MinFcnValue());

    std::cout << "D01:\t" << D01 << std::endl;
    std::cout << "D12:\t" << D12 << std::endl;
    std::cout << "D23:\t" << D23 << std::endl;
    std::cout << "D34:\t" << D34 << std::endl;
    std::cout << "D45:\t" << D45 << std::endl;
    std::cout << "D56:\t" << D56 << std::endl;
    std::cout << "D67:\t" << D67 << std::endl;
    std::cout << "D78:\t" << D78 << std::endl;
    std::cout << "D89:\t" << D89 << std::endl;
}
