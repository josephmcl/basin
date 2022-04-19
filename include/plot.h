#pragma once

#include <vector>

#include "definitions.h"

#include "TCanvas.h"
#include "TGraph2D.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TH1.h"
#include "TF1.h"
#include "TF2.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TAxis.h"

// Create, Draw and fit a TGraph2DErrors

namespace plot {

void matrix(std::vector<type::real_t> const &x, 
            std::vector<type::real_t> const &y, 
            std::vector<type::real_t> const &z){
   
    gStyle->SetPalette(kBird);
    const double e = 0.3;
    const int nd = 500;

    TCanvas *c = new TCanvas("c","Graph2D example",0,0,600,400);
    Int_t np = 200;
    TGraph2D *dt = new TGraph2D();
    dt->SetTitle("Graph title; X axis title; Y axis title; Z axis title");
    TRandom *r = new TRandom();

    for (auto i = 0; i < x.size(); ++i) {
        for (auto j = 0; j < y.size(); ++j) {
            auto k = (i * j) + j;
            dt->SetPoint(k, x[i], y[j], z[k]);
        }
    }
    // gStyle->SetPalette(1);
    c->SetTheta(85);
    c->SetPhi(280);
    dt->Draw("surf1");
    c->Print("example.pdf");
}

};