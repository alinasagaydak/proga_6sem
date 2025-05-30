#include <iostream>
#include "TFeldmanCousins.h"
#include "TRolke.h"
#include "TLimit.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TLegend.h"

void task8() {
    auto gr_fc = new TGraph();
    auto gr_rl = new TGraph();
    auto gr_l = new TGraph();
    
    double Nbkg = 10;
    double Nsignal = 2;

    TFeldmanCousins *fc = new TFeldmanCousins(0.95);
    TRolke *rl = new TRolke(0.95);
    
    TConfidenceLevel *myconfidence;
    double CLs = -100, s_ = 0.1, s = 0.1;
    double ul = 0;
    
    for (int i = 0; i < 7; i++) {
        gr_fc->SetPoint(i, Nsignal, fc->CalculateUpperLimit(Nsignal, Nbkg));
        rl->SetKnownBkgGaussEff(Nsignal, 1., 0, Nbkg);
        gr_rl->SetPoint(i, Nsignal, rl->GetUpperLimit());
                
        do {
            myconfidence = TLimit::ComputeLimit(s, Nbkg, Nsignal);
            if(myconfidence->CLs() <= 0.1 && CLs > 0.1) {
                ul = s_ + (0.1 - CLs)*(s - s_)/(myconfidence->CLs() - CLs);
                break;
            }
            CLs = myconfidence->CLs();
            s_=s;
            s += 0.1;
        } while(myconfidence->CLs() > 0.05);
        gr_l->SetPoint(i, Nsignal, ul);

        Nsignal+=2;
    }

    gr_fc->SetMarkerStyle(20); gr_fc->SetMarkerColor(kRed);
    gr_rl->SetMarkerStyle(21); gr_rl->SetMarkerColor(kBlue);
    gr_l->SetMarkerStyle(22); gr_l->SetMarkerColor(kGreen);

    auto gr = new TMultiGraph();
    gr->SetTitle(";Nsignal;Upper limit");
    gr->Add(gr_fc);
    gr->Add(gr_rl);
    gr->Add(gr_l);
    gr->Draw("ALP");

    auto legend = new TLegend(0.2, 0.7, 0.4, 0.8);
    legend->AddEntry(gr_fc, "TFeldmanCousins");
    legend->AddEntry(gr_rl, "TRolke");
    legend->AddEntry(gr_l, "TLimit");
    legend->Draw();


    Nsignal = 9;
    std::cout << "Nsignal = " << Nsignal << std::endl;
    std::cout << "/// TFeldmanCousins ///" << std::endl;
    std::cout << "FC upper limit = " << fc -> CalculateUpperLimit(Nsignal, Nbkg) << std::endl;
    std::cout << "FC lower limit = " << fc->CalculateLowerLimit(Nsignal, Nbkg) << std::endl;

    std::cout << "/// TRolke ///" << std::endl;
    rl->SetKnownBkgGaussEff(Nsignal, 1., 0, Nbkg);
    std::cout << "Rolke upper limit = " << rl->GetUpperLimit() << std::endl;
    std::cout << "Rolke lower limit = " << rl->GetLowerLimit() << std::endl;


    do {
        myconfidence = TLimit::ComputeLimit(s, Nbkg, Nsignal);
        if(myconfidence->CLs() <= 0.1 && CLs > 0.1) {
            ul = s_ + (0.1 - CLs)*(s - s_)/(myconfidence->CLs() - CLs);
            break;
        }
        CLs = myconfidence->CLs();
        s_=s;
        s += 0.1;
    } while(myconfidence->CLs() > 0.05);
    std::cout << "/// TLimit ///" << std::endl;
    std::cout << "CLs upper limit = " << ul << std::endl;
}
