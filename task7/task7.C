#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TLegend.h"

class Event {
public:
    Event():
    E(1.02),
    m_eta(0.547862),
    sigmaE(0.04),
    sigmaT(0.005),
    alfa(9),alfa0(9),Valfa(9,9),d(5),D(5,9)
    {
        rnd = gRandom;
        m[0]=0.; m[1]=0.; m[2]=0.;
    }
    ~Event(){
    }
    
    void GenEvent(){ 
        
        double s1,s2,s3;
        double s3min,s3max,s1min,s1max;
        double wei;                
        
        double weimax = ((E-m[0])*(E-m[0]) - (m[1]+m[2])*(m[1]+m[2]))*
            ((E-m[2])*(E-m[2]) - (m[0]+m[1])*(m[0]+m[1]))*8*TMath::Pi()*TMath::Pi();
        
        do{
            // s1 = (p2+p3)^2
            // s3 = (p1+p2)^2
            s1min = (m[1]+m[2])*(m[1]+m[2]);
            s1max = (E-m[0])*(E-m[0]);
            s1 = s1min + rnd->Rndm()*(s1max - s1min);
            
            do{
                s1 = rnd->BreitWigner(m_eta*m_eta,1.3e-4); // 1.3e-6 for eta
            }while(s1>s1max || s1<s1min);
            
            s3min=m[0]*m[0]+m[1]*m[1]-1/2./s1*
            (   (s1-E*E+m[0]*m[0])*(s1+m[1]*m[1]-m[2]*m[2])+
            sqrt( Lambda(s1,E*E,m[0]*m[0])*
            Lambda(s1,m[1]*m[1],m[2]*m[2]) )   );
        
            s3max=m[0]*m[0]+m[1]*m[1]-1/2./s1*
            (   (s1-E*E+m[0]*m[0])*(s1+m[1]*m[1]-m[2]*m[2])-
            sqrt( Lambda(s1,E*E,m[0]*m[0])*
            Lambda(s1,m[1]*m[1],m[2]*m[2]) )   );
            s3 = s3min + rnd->Rndm()*(s3max - s3min);
                
            wei = (s1max-s1min)*(s3max-s3min)*8*TMath::Pi()*TMath::Pi();
        }while( wei < weimax * rnd->Rndm() );
        
        s2=E*E+m[0]*m[0]+m[1]*m[1]+m[2]*m[2]-s1-s3;
        
        double E1 =(E*E+m[0]*m[0]-s1)/2./E;
        double P_1_=sqrt(E1*E1-m[0]*m[0]);
        double E2 =(E*E+m[1]*m[1]-s2)/2./E;
        double P_2_=sqrt(E2*E2-m[1]*m[1]);
        double E3 =(E*E+m[2]*m[2]-s3)/2./E;
        double P_3_=sqrt(E3*E3-m[2]*m[2]);
        
        double cos2=(P_1_*P_1_+P_2_*P_2_-P_3_*P_3_)/2./P_1_/P_2_;
        double cos3=(P_1_*P_1_+P_3_*P_3_-P_2_*P_2_)/2./P_1_/P_3_;
        double sin2=sqrt(1-cos2*cos2);
        double sin3=sqrt(1-cos3*cos3);                       
        
        double phi3 = 2*TMath::Pi()*rnd->Rndm();
        P[0].SetXYZM(0,0,P_1_,m[0]);
        P[1].SetXYZM(P_2_*sin2*cos(phi3),P_2_*sin2*sin(phi3),-P_2_*cos2,m[1]);
        P[2].SetXYZM(-P_3_*sin3*cos(phi3),-P_3_*sin3*sin(phi3),-P_3_*cos3,m[2]);
        
        double cos1 = -1 + 2*rnd->Rndm();
        double phi1 = 2*TMath::Pi()*rnd->Rndm();
        TRotation Oy; Oy.RotateY(acos(cos1));
        TRotation Oz; Oz.RotateZ(phi1);
                                 
        P[0].Transform(Oz*Oy);
        P[1].Transform(Oz*Oy);
        P[2].Transform(Oz*Oy); 
    }
    double Lambda(double x, double y, double z){
        return (x-y-z)*(x-y-z)-4.*y*z;        
    }
    void DetectorRespond(){
        double Enew = P[0].E()+ gRandom->Gaus(0,sigmaE*P[0].E());
        Enew = Enew > 0 ? Enew : 0.;
        P[0].SetRho( sqrt(Enew*Enew-m[0]*m[0]));
        P[0].SetE(Enew);
        Enew = P[1].E()+ gRandom->Gaus(0,sigmaE*P[1].E());
        Enew = Enew > 0 ? Enew : 0.;
        P[1].SetRho( sqrt(Enew*Enew-m[1]*m[1]));
        P[1].SetE(Enew);
        Enew = P[2].E()+ gRandom->Gaus(0,sigmaE*P[2].E());
        Enew = Enew > 0 ? Enew : 0.;
        P[2].SetRho( sqrt(Enew*Enew-m[2]*m[2]));
        P[2].SetE(Enew);
        
        P[0].SetTheta( P[0].Theta()+ gRandom->Gaus(0,sigmaT/sqrt(2)*TMath::Pi()));
        P[1].SetTheta( P[1].Theta()+ gRandom->Gaus(0,sigmaT/sqrt(2)*TMath::Pi()));
        P[2].SetTheta( P[2].Theta()+ gRandom->Gaus(0,sigmaT/sqrt(2)*TMath::Pi()));
        
        double sigmaT_ = sigmaT/sqrt(2)/TMath::Sin(P[0].Theta());
        P[0].SetPhi( P[0].Phi()+ gRandom->Gaus(0,sigmaT_*TMath::Pi()));
        sigmaT_ = sigmaT/sqrt(2)/TMath::Sin(P[1].Theta());
        P[1].SetPhi( P[1].Phi()+ gRandom->Gaus(0,sigmaT_*TMath::Pi()));
        sigmaT_ = sigmaT/sqrt(2)/TMath::Sin(P[2].Theta());
        P[2].SetPhi( P[2].Phi()+ gRandom->Gaus(0,sigmaT_*TMath::Pi()));
    }    
    void KinFitInit(){
        alfa*=0; alfa0*=0; Valfa*=0;
        for(int i=0;i<3;i++){
            alfa[3*i]   = P[i].Rho();
            alfa[3*i+1] = P[i].Theta();
            alfa[3*i+2] = P[i].Phi(); 
            
            Valfa[3*i][3*i]     = sigmaE*alfa[3*i] * sigmaE*alfa[3*i];
            Valfa[3*i+1][3*i+1] = sigmaT*TMath::Pi()/sqrt(2) * sigmaT*TMath::Pi()/sqrt(2);
            Valfa[3*i+2][3*i+2] = sigmaT*TMath::Pi()/sqrt(2)/TMath::Sin(alfa[3*i+1]) * sigmaT*TMath::Pi()/sqrt(2)/TMath::Sin(alfa[3*i+1]); 
        }
        alfa0 = alfa;
    }
    
    double KinFitIter(){
        D*=0; d*=0;
        for(int i=0;i<3;i++){            
            D[0][3*i]   = alfa[3*i]/sqrt(alfa[3*i]*alfa[3*i]+m[i]*m[i]);
            
            D[1][3*i]   = sin(alfa[3*i+1])*cos(alfa[3*i+2]);
            D[1][3*i+1] = alfa[3*i]*cos(alfa[3*i+1])*cos(alfa[3*i+2]);
            D[1][3*i+2] = -alfa[3*i]*sin(alfa[3*i+1])*sin(alfa[3*i+2]);
            
            D[2][3*i]   = sin(alfa[3*i+1])*sin(alfa[3*i+2]);
            D[2][3*i+1] = alfa[3*i]*cos(alfa[3*i+1])*sin(alfa[3*i+2]);
            D[2][3*i+2] = alfa[3*i]*sin(alfa[3*i+1])*cos(alfa[3*i+2]);
            
            D[3][3*i]   = cos(alfa[3*i+1]);
            D[3][3*i+1] = -alfa[3*i]*sin(alfa[3*i+1]);

            d[0] = d[0] + sqrt(alfa[3*i]*alfa[3*i]+m[i]*m[i]);
            d[1] = d[1] + alfa[3*i]*sin(alfa[3*i+1])*cos(alfa[3*i+2]);
            d[2] = d[2] + alfa[3*i]*sin(alfa[3*i+1])*sin(alfa[3*i+2]);
            d[3] = d[3] + alfa[3*i]*cos(alfa[3*i+1]);            
        }
        d[0] = d[0] - E;

 ////// Ограничение на инв. массу. Фотоны 1 и 2 из эта-мезона? (нумерация с 0) //////       
        double E_1 = sqrt(alfa[3*1]*alfa[3*1]+m[1]*m[1]);
        double E_2 = sqrt(alfa[3*2]*alfa[3*2]+m[2]*m[2]);

        double p_x1 = alfa[3*1]*sin(alfa[3*1+1])*cos(alfa[3*1+2]);
        double p_y1 = alfa[3*1]*sin(alfa[3*1+1])*sin(alfa[3*1+2]);
        double p_z1 = alfa[3*1]*cos(alfa[3*1+1]);

        double p_x2 = alfa[3*2]*sin(alfa[3*2+1])*cos(alfa[3*2+2]);
        double p_y2 = alfa[3*2]*sin(alfa[3*2+1])*sin(alfa[3*2+2]);
        double p_z2 = alfa[3*2]*cos(alfa[3*2+1]);
        
        double p_sum_x = p_x1 + p_x2; 
        double p_sum_y = p_y1 + p_y2; 
        double p_sum_z = p_z1 + p_z2; 

        for (int i = 1; i < 3; i++) {
            double tmp_p1 = 2 * (E_1 + E_2) * alfa[3*i]/sqrt(alfa[3*i]*alfa[3*i]+m[i]*m[i]);
            double tmp_p2 = 2 * p_sum_x * sin(alfa[3*i+1])*cos(alfa[3*i+2]);
            double tmp_p3 = 2 * p_sum_y * sin(alfa[3*i+1])*sin(alfa[3*i+2]);
            double tmp_p4 = 2 * p_sum_z * cos(alfa[3*i+1]);

            double tmp_th1 = 2 * p_sum_x * alfa[3*i] * cos(alfa[3*i+1])*cos(alfa[3*i+2]);
            double tmp_th2 = 2 * p_sum_y * alfa[3*i] * cos(alfa[3*i+1])*sin(alfa[3*i+2]);
            double tmp_th3 = 2 * p_sum_z * alfa[3*i] * sin(alfa[3*i+1]);

            double tmp_phi1 = 2 * p_sum_x * alfa[3*i] * sin(alfa[3*i+1])*sin(alfa[3*i+2]);
            double tmp_phi2 = 2 * p_sum_y * alfa[3*i] * sin(alfa[3*i+1])*cos(alfa[3*i+2]);

            D[4][3*i] = tmp_p1 - tmp_p2 + tmp_p3 + tmp_p4;
            D[4][3*i+1] = - tmp_th1 - tmp_th2 + tmp_th3;
            D[4][3*i+2] = tmp_phi1 - tmp_phi2;
        }
        d[4] = sqrt( (E_1 + E_2)*(E_1 + E_2) - (p_sum_x*p_sum_x + p_sum_y*p_sum_y + p_sum_z*p_sum_z)) - m_eta;
        
        

        TMatrixD VDinv = TMatrixD( TMatrixD(D,TMatrixD::kMult,Valfa),TMatrixD::kMultTranspose,D);
        TMatrixD VD = TMatrixD(TMatrixD::kInverted,VDinv);
        TVectorD lambda = VD* ( D*(alfa0 - alfa) + d);
        alfa = alfa0 - TMatrixD(Valfa,TMatrixD::kMultTranspose,D)*lambda;
        double chi2_ = lambda*(VDinv*lambda);
        return chi2_;
    }
    
    double KinFit(){
        KinFitInit();
        double chi2_=1e10, chi2=1e12;
        for(int i=0;i<100;i++){
            chi2 = KinFitIter(); 
            if(abs(chi2-chi2_)<1e-5) break;
            chi2_ = chi2;
        }
        for(int i=0;i<3;i++){
            double Px = alfa[3*i]*sin(alfa[3*i+1])*cos(alfa[3*i+2]);
            double Py = alfa[3*i]*sin(alfa[3*i+1])*sin(alfa[3*i+2]);
            double Pz = alfa[3*i]*cos(alfa[3*i+1]);
            P[i].SetXYZM(Px,Py,Pz,m[i]);
        }
        return chi2;
    }
    
    void SetE(double E_=1.02){E = E_;}
    void SetResolution(double sigmaE_, double sigmaT_){sigmaE = sigmaE_;  sigmaT = sigmaT_;}
    double GetM12(){ return (P[0]+P[1]).M2(); }
    double GetM23(){ return (P[1]+P[2]).M2(); }
    double GetM13(){ return (P[0]+P[2]).M2(); }
    double GetE1(){ return P[0].E(); }
    double GetE2(){ return P[1].E(); }
    
    TLorentzVector* GetP_1(){return &P[0];}
    TLorentzVector* GetP_2(){return &P[1];}
    TLorentzVector* GetP_3(){return &P[2];}
    TLorentzVector* GetP_kin1(){return &(P_kin[0]);}
    TLorentzVector* GetP_kin2(){return &(P_kin[1]);}
    TLorentzVector* GetP_kin3(){return &(P_kin[2]);}
private:
    TLorentzVector P[3];    
    TLorentzVector P_kin[3];
    TVectorD alfa,alfa0,d;
    double vec[9];
    TMatrixD Valfa,D;
//    double m[0],m[1],m[2];
    double m[3];
    double E; // GeV
    double m_eta;
    double sigmaE,sigmaT;
    TRandom* rnd;
    double arglist[10];
    int ierflg;
};

void task7() {
    Event evt;

    auto c1 = new TCanvas("c1", "", 900, 1200);
    c1->Divide(1,2);

    auto gr = new TGraph(10000);
    auto h = new TH1D("h","",100,0.2,0.4);
    auto h2 = new TH1D("h2","",100,0.2,0.4);
    auto h3 = new TH1D("h3","",100,0.2,0.4);
    auto hchi2 = new TH1D("hchi2","",100,0.,15);
    auto fchi2 = new TF1("fchi2","[0]*ROOT::Math::chisquared_pdf(x,[1])",0,15);
    fchi2->SetParameters(10000*15/100,5);

    for(int i = 0; i < 10000; i++) {
        evt.GenEvent();
        h->Fill(evt.GetM23());
        evt.DetectorRespond();
    //    gr.SetPoint(i, evt.GetE1(),evt.GetE2());
        gr->SetPoint(i, evt.GetM12(),evt.GetM23());
        h2->Fill(evt.GetM23());
        double chi2 = evt.KinFit();
        h3->Fill(evt.GetM23());
        hchi2->Fill(chi2);
    }

    c1->cd(1);
    h->Draw();
    h2->SetLineColor(kRed);
    h3->SetLineColor(kGreen);
    h2->Draw("same");
    h3->Draw("same");

    auto legend = new TLegend(0.2, 0.6, 0.4, 0.8);
    legend->AddEntry(h, "Gen events");
    legend->AddEntry(h2, "Row data");
    legend->AddEntry(h3, "Kinematic fit");
    legend->Draw();

    c1->cd(2);
    hchi2->Fit("fchi2","L","ep");


}
