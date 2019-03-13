#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
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

void generate_bins_from(std::vector<float>& edges, int& count, std::string tag) {
    if (count && edges.size() != 2)
        printf("invalid bin configuration for %s\n", tag.data());

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

int area4(char const* config, char const* output) {
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

    auto min_rho_fit = conf->get<float>("min_rho_fit");
    auto max_rho_fit = conf->get<float>("max_rho_fit");

    generate_bins_from(etas, netas, "eta");
    generate_bins_from(rhos, nrhos, "rho");
    generate_bins_from(isos, nisos, "iso");

    TFile* f = new TFile(input.data(), "read");
    TTree* t = (TTree*)f->Get("electrons");
    auto e = new electrontree(true, false, false, t);

    TFile* fout = new TFile(output, "recreate");
    fout->cd();

    TH1::SetDefaultSumw2();

    TH2F** h_iso_rho = new TH2F*[netas];
    TProfile** p_iso_rho = new TProfile*[netas];
    TH1F** h_isooverrho = new TH1F*[netas];

    TF1** f_iso_rho = new TF1*[netas];

    for (int i = 0; i < netas; ++i) {
        std::string index = std::to_string(i);
        h_iso_rho[i] = new TH2F(("h_iso_rho_eta" + index).data(),
                                "", nrhos, &rhos[0], nisos, &isos[0]);
        p_iso_rho[i] = new TProfile(("p_iso_rho_eta" + index).data(),
                                    "", nrhos, &rhos[0]);
        h_isooverrho[i] = new TH1F(("h_isooverrho_eta" + index).data(),
                                   "", nrhos, &rhos[0]);
        f_iso_rho[i] = new TF1(("f_iso_rho_eta" + index).data(),
                               "pol1", min_rho_fit, max_rho_fit);
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
            p_iso_rho[ieta]->Fill(rho, sum_iso);
            h_isooverrho[ieta]->Fill(sum_iso / rho);
        }
    }

    /* fit effective areas */
    float* eff_area = new float[netas];
    float* eff_area_err = new float[netas];

    for (int i = 0; i < netas; ++i) {
        f_iso_rho[i]->SetParameters(1, 0.1);

        std::string index = std::to_string(i);
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

    {
        for (int i = 0; i < netas; ++i) {
            htitle(h_iso_rho[i], Form(
                "%.3f < |#eta| < %.3f;#rho;#Sigma#it{iso}", etas[i], etas[i+1]));
            h_iso_rho[i]->SetStats(0);
            h_iso_rho[i]->Draw("colz");

            c1->SaveAs(Form("a4_iso_rho_eta%i.png", i));
        }
    }

    {
        for (int i = 0; i < netas; ++i) {
            hstyle(p_iso_rho[i], 21, 1, 0.6);
            hformat(p_iso_rho[i], 0.f, 15.f, Form(
                "%.3f < |#eta| < %.3f;#rho;#LT#Sigma#it{iso}#GT", etas[i], etas[i+1]));
            p_iso_rho[i]->SetStats(0);
            p_iso_rho[i]->Draw("pe");

            c1->SaveAs(Form("a4_p_iso_rho_eta%i.png", i));
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
