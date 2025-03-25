#include "TH1.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TCanvas.h"

struct GlobalLL2 {
   GlobalLL2(  TH1D* hist_, TF1* func_) :
      hist(hist_),func(func_) {}

   double operator() (const double *par) const {
       func->SetParameters(par);
       double NLL=0;
       double s = 5.;
       double error_s = 0.1;
       for(int i=1;i<=hist->GetNbinsX();i++){
           double mu = func->Eval( hist->GetBinCenter(i) );
           mu = TMath::Max(mu,0.);
           double n = hist->GetBinContent(i);
           NLL += ( mu - n + n*(ROOT::Math::Util::EvalLog(n)- ROOT::Math::Util::EvalLog(mu)) );
       }
       NLL += 0.5 * pow((par[2]-s)/error_s, 2);
      return NLL;
   }
    TH1D* hist; 
    TF1* func; 
};

double fit_func(double* x, double* par) {
    return par[0] + par[1] * TMath::Gaus(x[0], par[2], par[3], kTRUE);
}

void task3() {
    TVirtualFitter::SetDefaultFitter("Minuit2");

/// Create data ///
    auto gaus_func = new TF1("gaus_func", "TMath::Gaus(x, [0], [1], kTRUE)", -50., 50.);
    gaus_func->SetParameters(0, 5);
    auto h = new TH1D("h", "Gaus + const", 100, -50., 50.);
    h->FillRandom("pol0", 100000);
    h->FillRandom("gaus_func", 700);
    
    auto c1 = new TCanvas("c1", "", 900, 600);
    c1->SetGrid();
    h->Draw("e");

/// Fitting ///
    ROOT::Fit::DataOptions opt;
    opt.fUseEmpty = true;
    ROOT::Fit::DataRange range;

    range.SetRange(-50., 50.);
    ROOT::Fit::BinData data(opt,range);
    ROOT::Fit::FillData(data, h);

    auto fitfunc = new TF1("fitfunc", fit_func, -50., 50., 4);
    fitfunc->SetParNames("const", "ampl", "mean", "sigma");
    GlobalLL2 globalLL2(h,fitfunc);

    ROOT::Fit::Fitter fitter;
    fitter.Config().MinimizerOptions().SetErrorDef(0.5);
    fitter.Config().MinimizerOptions().SetPrintLevel(0);
    double par[4]={1000., 1., 0., 1.};

    fitter.FitFCN(4, globalLL2, par, data.Size(), false);
    ROOT::Fit::FitResult result = fitter.Result();
    result.Print(std::cout);

/// Draw fit function ///
    fitfunc->SetParameters(result.GetParams());
    fitfunc->SetLineColor(kRed);
    fitfunc->Draw("same");
}
