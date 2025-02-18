#include "TLegend.h"

TF1* initGausFunc(double min, double max, double N) {
    auto f = new TF1("f", "TMath::Gaus(x, [0], sqrt([0]), kTRUE)", min, max);
    f->SetLineColor(kRed);
    f->SetTitle(" ; ;");
    f->SetParameter(0, N);
    f->SetNpx(10000); 
    return f;
}

TF1* initPoissFunc(double min, double max, double N) {
    auto f = new TF1("f", "TMath::Poisson(x, [0])", min, max);
    f->SetLineColor(kBlue);
    f->SetParameter(0, N);
    f->SetNpx(10000); 
    return f;
}

void task1_2() {
    double min = 0.;
    double max = 120.;

    double N1 = 3.;
    double N2 = 10.;
    double N3 = 15.;
    double N4 = 20.;
    double N5 = 25.;
    double N6 = 30.;
    double N7 = 40.;
    double N8 = 50.;
    double N9 = 60.;
    double N10 = 70.;
    double N11 = 80.;
    double N12 = 90.;
    double N13 = 100.;

    auto f1_1 = initGausFunc(min, max, N1);
    auto f2_1 = initPoissFunc(min, max, N1);
    f1_1->Draw();
    f2_1->Draw("same");

    auto f1_2 = initGausFunc(min, max, N2);
    auto f2_2 = initPoissFunc(min, max, N2);
    f1_2->Draw("same");
    f2_2->Draw("same");

    auto f1_3 = initGausFunc(min, max, N3);
    auto f2_3 = initPoissFunc(min, max, N3);
    f1_3->Draw("same");
    f2_3->Draw("same");

    auto f1_4 = initGausFunc(min, max, N4);
    auto f2_4 = initPoissFunc(min, max, N4);
    f1_4->Draw("same");
    f2_4->Draw("same");

    auto f1_5 = initGausFunc(min, max, N5);
    auto f2_5 = initPoissFunc(min, max, N5);
    f1_5->Draw("same");
    f2_5->Draw("same");

    auto f1_6 = initGausFunc(min, max, N6);
    auto f2_6 = initPoissFunc(min, max, N6);
    f1_6->Draw("same");
    f2_6->Draw("same");

    auto f1_7 = initGausFunc(min, max, N7);
    auto f2_7 = initPoissFunc(min, max, N7);
    f1_7->Draw("same");
    f2_7->Draw("same");

    auto f1_8 = initGausFunc(min, max, N8);
    auto f2_8 = initPoissFunc(min, max, N8);
    f1_8->Draw("same");
    f2_8->Draw("same");

    auto f1_9 = initGausFunc(min, max, N9);
    auto f2_9 = initPoissFunc(min, max, N9);
    f1_9->Draw("same");
    f2_9->Draw("same");

    auto f1_10 = initGausFunc(min, max, N10);
    auto f2_10 = initPoissFunc(min, max, N10);
    f1_10->Draw("same");
    f2_10->Draw("same");

    auto f1_11 = initGausFunc(min, max, N11);
    auto f2_11 = initPoissFunc(min, max, N11);
    f1_11->Draw("same");
    f2_11->Draw("same");

    auto f1_12 = initGausFunc(min, max, N12);
    auto f2_12 = initPoissFunc(min, max, N12);
    f1_12->Draw("same");
    f2_12->Draw("same");

    auto f1_13 = initGausFunc(min, max, N13);
    auto f2_13 = initPoissFunc(min, max, N13);
    f1_13->Draw("same");
    f2_13->Draw("same");

    auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(f1_13, "Gaus");
    legend->AddEntry(f2_13, "Poisson");
    legend->Draw();
}
