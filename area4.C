#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"

#include <string>
#include <vector>

#include "git/config/include/configurer.h"

#include "include/cosmetics.h"
#include "include/electrontree.h"

float area_for_abs_sceta(float eta) {
    static constexpr float areas[7] = {
        0.1440, 0.1562, 0.1032, 0.0859, 0.1116, 0.1321, 0.1654
    };

    static constexpr float etas[7] = {
        1.000, 1.479, 2.000, 2.200, 2.300, 2.400, 2.500
    };

    int ieta = 0;
    for (; eta > etas[ieta] && ieta < 7; ++ieta);

    return areas[ieta];
}

void generate_bins_from(std::vector<float>& edges, int& count,
                        std::string var) {
    if (count && edges.size() != 2)
        printf("invalid bin configuration for %s\n", var.data());

    if (count) {
        float base = edges[0];
        float interval = (edges[1] - edges[0]) / count;

        edges = std::vector<float>(count + 1);
        std::iota(edges.begin(), edges.end(), 0);
        std::for_each(edges.begin(), edges.end(), [&] (float& f) {
            f = base + f * interval; });
    } else {
        count = edges.size() - 1;
    }
}

double get_cutoff_point(TH1D* hist, double cutoff) {
    double quantile = 9999.;
    hist->GetQuantiles(1, &quantile, &cutoff);

    return quantile;
}

float evaluate_cutoff_error(TH1D* slice, TH1D* h_isooverrho, int ieta,
                            int irho, float rho_average, double cutoff,
                            int nisos, std::vector<float>& isos) {
    TH1D* h_trial_cutoffs = new TH1D(
        Form("h_trial_cutoffs_eta%i_rho%i", ieta, irho), "", nisos, &isos[0]);

    int zero_size = slice->GetBinContent(0);
    int sample_size = slice->GetSumOfWeights();

    constexpr int trials = 100;
    for (int i = 0; i < trials; ++i) {
        int trial_size_zero = gRandom->Poisson(zero_size);
        int trial_size_nonzero = gRandom->Poisson(sample_size - zero_size);

        TH1D* h_trial = new TH1D(Form("trial%i_eta%i_rho%i", i, ieta, irho),
                                 "", nisos, &isos[0]);
        h_trial->SetBinContent(0, trial_size_zero);
        for (int j = 0; j < trial_size_nonzero; ++j)
            h_trial->Fill(h_isooverrho->GetRandom() * rho_average);

        double trial_quantile = get_cutoff_point(h_trial, cutoff);
        h_trial_cutoffs->Fill(trial_quantile);
    }

    return h_trial_cutoffs->GetRMS();
}

TH1D* find_90pc_cutoff(TH2D* h_iso_rho, TH1D* h_isooverrho, double cutoff,
                           int ieta, int nrhos, std::vector<float>& rhos,
                           int nisos, std::vector<float>& isos) {
    std::string index = std::to_string(ieta);
    TH1D* h_90pc_cutoff = new TH1D(("p_iso_rho_eta" + index).data(), "",
                                   nrhos, &rhos[0]);

    for (int i = 1; i <= nrhos; ++i) {
        auto slice = h_iso_rho->ProjectionY(Form("_slice%i", i), i, i);
        float rho_average = h_iso_rho->GetXaxis()->GetBinCenter(i);
        h_90pc_cutoff->SetBinContent(i, get_cutoff_point(slice, cutoff));
        h_90pc_cutoff->SetBinError(i, evaluate_cutoff_error(
            slice, h_isooverrho, ieta, i, rho_average, cutoff, nisos, isos));
    }

    return h_90pc_cutoff;
}

int area4(char const* config, char const* tag) {
    auto conf = new configurer(config);

    auto input = conf->get<std::string>("input");
    auto max_entries = conf->get<int64_t>("max_entries");

    auto etas = conf->get<std::vector<float>>("etas");
    auto rhos = conf->get<std::vector<float>>("rhos");
    auto isos = conf->get<std::vector<float>>("isos");

    auto netas = conf->get<int>("netas");
    auto nrhos = conf->get<int>("nrhos");
    auto nisos = conf->get<int>("nisos");

    auto min_pt = conf->get<float>("min_pt");
    auto min_rho_fit = conf->get<std::vector<float>>("min_rho_fit");
    auto max_rho_fit = conf->get<std::vector<float>>("max_rho_fit");

    auto use_90pc_eff = conf->get<bool>("use_90pc_eff");

    generate_bins_from(etas, netas, "eta");
    generate_bins_from(rhos, nrhos, "rho");
    generate_bins_from(isos, nisos, "iso");

    TFile* f = new TFile(input.data(), "read");
    TTree* t = (TTree*)f->Get("electrons");
    auto e = new electrontree(true, false, false, t);

    TFile* fout = new TFile(Form("eff-area-%s.root", tag), "recreate");

    TH1::SetDefaultSumw2();
    TH2D** h_iso_rho = new TH2D*[netas];
    TH1D** p_iso_rho = new TH1D*[netas];
    TH1D** h_isooverrho = new TH1D*[netas];
    TF1** f_iso_rho = new TF1*[netas];

    for (int i = 0; i < netas; ++i) {
        std::string index = std::to_string(i);
        h_iso_rho[i] = new TH2D(("h_iso_rho_eta" + index).data(),
                                "", nrhos, &rhos[0], nisos, &isos[0]);
        h_isooverrho[i] = new TH1D(("h_isooverrho_eta" + index).data(),
                                   "", nrhos, &rhos[0]);
        f_iso_rho[i] = new TF1(("f_iso_rho_eta" + index).data(),
                               "pol1", min_rho_fit[i], max_rho_fit[i]);
    }

    int64_t nentries = t->GetEntries();
    if (max_entries) nentries = std::min(nentries, max_entries);

    printf("entries: %li\n", nentries);
    for (int64_t i = 0; i < nentries; ++i) {
        t->GetEntry(i);

        if (i % 100000 == 0)
            printf("entry: %li\n", i);

        for (int j = 0; j < e->nEle; ++j) {
            if ((*e->elePt)[j] < min_pt)
                continue;
            if ((*e->eleGenMatchIndex)[j] < 0)
                continue;
            float eta = (*e->eleSCEta)[j];
            if (eta < etas[0] || eta > etas[netas])
                continue;

            /* i'm stupid and forgot to save rho */
            float area = area_for_abs_sceta(std::abs((*e->eleSCEta)[j]));
            float rho = (*e->eleEffAreaTimesRho)[j] / area;
            float sum_iso = (*e->elePFChIso04)[j] + (*e->elePFPhoIso04)[j]
                + (*e->elePFNeuIso04)[j];

            int ieta = 0;
            for (; ieta < netas && eta > etas[ieta + 1]; ++ieta);

            h_iso_rho[ieta]->Fill(rho, sum_iso);
            h_isooverrho[ieta]->Fill(sum_iso / rho);
        }
    }

    /* fit effective areas */
    float* eff_area = new float[netas];
    float* eff_area_err = new float[netas];

    for (int i = 0; i < netas; ++i) {
        std::string index = std::to_string(i);

        p_iso_rho[i] = !use_90pc_eff
            ? h_iso_rho[i]->ProfileX(("p_iso_rho_eta" + index).data())
            /* derive 90% efficiency cutoff */
            : find_90pc_cutoff(h_iso_rho[i], h_isooverrho[i], 0.9, i, nrhos,
                               rhos, nisos, isos);

        f_iso_rho[i]->SetParameters(1, 0.1);
        p_iso_rho[i]->Fit(("f_iso_rho_eta" + index).data(), "RME");

        eff_area[i] = f_iso_rho[i]->GetParameter(1);
        eff_area_err[i] = f_iso_rho[i]->GetParError(1);
    }

    for (int i = 0; i < netas; ++i) {
        printf("eta: [%.3f, %.3f] %.4f (%.4f)\n", etas[i], etas[i+1],
               eff_area[i], eff_area_err[i]);
    }

    /* painting */
    TCanvas* c1 = new TCanvas("c1", "", 400, 400);

    std::string type = use_90pc_eff ? "90pc" : "meaniso";

    {
        for (int i = 0; i < netas; ++i) {
            htitle(h_iso_rho[i], Form(
                "%.3f < |#eta| < %.3f;#rho;#Sigma#it{iso}", etas[i], etas[i+1]));
            h_iso_rho[i]->SetStats(0);
            h_iso_rho[i]->Draw("colz");

            hstyle(p_iso_rho[i], 21, 1, 0.6);
            p_iso_rho[i]->SetStats(0);
            p_iso_rho[i]->Draw("same pe");

            c1->SetLogz(1);
            c1->SaveAs(Form("a4_iso_rho_eta%i-%s-%s.pdf", i, type.data(), tag));
        }
    }

    {
        for (int i = 0; i < netas; ++i) {
            hstyle(p_iso_rho[i], 21, 1, 0.6);
            hformat(p_iso_rho[i], 0.f, rhos.back(), Form(
                "%.3f < |#eta| < %.3f;#rho;#LT#Sigma#it{iso}#GT", etas[i], etas[i+1]));
            p_iso_rho[i]->SetStats(0);
            p_iso_rho[i]->Draw("pe");

            c1->SetLogz(1);
            c1->SaveAs(Form("a4_p_iso_rho_eta%i-%s-%s.pdf", i, type.data(), tag));
        }
    }

    fout->Write("", TObject::kOverwrite);
    fout->Close();

    f->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return area4(argv[1], argv[2]);

    return 1;
}
