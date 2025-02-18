void task1_1() {
    auto hist = new TH1D("hist", "", 100, -1., 1.);
    auto rand = new TRandom3(0);
    for (int i = 0; i < 25000; i++) {double x = 2 * (rand->Rndm(0)) - 1; hist->Fill(x); }
    hist->SetLineWidth(2);
    hist->Draw("e");
    hist->Fit("pol0");

    int out1sigma = 0;
    int out2sigma = 0;
    int out3sigma = 0;

    for (int i = 1; i < 100; i++) {
        double val = hist->GetBinContent(i);
        if (fabs(val - 250) >= sqrt(250)) out1sigma++;
        if (fabs(val - 250) >= 2 * sqrt(250)) out2sigma++;
        if (fabs(val - 250) >= 3 * sqrt(250)) out3sigma++;
    }
    std::cout << "Beyond 1 sigma:\t" << out1sigma << std::endl;
    std::cout << "Beyond 2 sigma:\t" << out2sigma << std::endl;
    std::cout << "Beyond 3 sigma:\t" << out3sigma << std::endl;

}
