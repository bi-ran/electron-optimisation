#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Reader.h"

#include "TMVA/deviations.h"
#include "TMVA/regression_averagedevs.h"

#include <map>
#include <string>
#include <vector>

#include "include/cosmetics.h"

int regress(char const* tag, char const* input, char const* output,
            int options) {
    TMVA::Tools::Instance();

    std::map<std::string, bool> use;
    use["BDT"] = true;
    use["BDTG"] = true;

    TFile* fout = TFile::Open(output, "recreate");

    TMVA::Factory* factory = new TMVA::Factory(
        "TMVARegression",
        fout,
        "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression");

    TMVA::DataLoader* loader = new TMVA::DataLoader(tag);

    loader->AddVariable("eleSCRawEn", "eleSCRawEn", "", 'F');
    loader->AddVariable("eleSCEtaWidth", "eleSCEtaWidth", "", 'F');
    loader->AddVariable("eleSCPhiWidth", "eleSCPhiWidth", "", 'F');
    loader->AddVariable("eleSeedEn/eleSCRawEn", "eleSeedEn", "", 'F');
    loader->AddVariable("eleSeedE5x5/eleSCRawEn", "eleSeedE5x5", "", 'F');
    loader->AddVariable("eleSeedE3x3/eleSCRawEn", "eleSeedE3x3", "", 'F');
    loader->AddVariable("eleHoverE", "eleHoverE", "", 'F');
    /* loader->AddVariable("rho", "rho", "", 'F'); */
    loader->AddVariable("eledEtaSCSeed", "eledEtaSCSeed", "", 'F');
    loader->AddVariable("eledPhiSCSeed", "eledPhiSCSeed", "", 'F');
    loader->AddVariable("eleSEE", "eleSEE", "", 'F');
    loader->AddVariable("eleSEP", "eleSEP", "", 'F');
    loader->AddVariable("eleSPP", "eleSPP", "", 'F');
    loader->AddVariable("eleSeedEMax/eleSeedE5x5", "eleSeedEMax", "", 'F');
    loader->AddVariable("eleSeedE2nd/eleSeedE5x5", "eleSeedE2nd", "", 'F');
    loader->AddVariable("eleSeedETop/eleSeedE5x5", "eleSeedETop", "", 'F');
    loader->AddVariable("eleSeedEBottom/eleSeedE5x5", "eleSeedEBottom", "", 'F');
    loader->AddVariable("eleSeedELeft/eleSeedE5x5", "eleSeedELeft", "", 'F');
    loader->AddVariable("eleSeedERight/eleSeedE5x5", "eleSeedERight", "", 'F');
    loader->AddVariable("eleSeedE2x5Max/eleSeedE5x5", "eleSeedE2x5Max", "", 'F');
    loader->AddVariable("eleSeedE2x5Top/eleSeedE5x5", "eleSeedE2x5Top", "", 'F');
    loader->AddVariable("eleSeedE2x5Bottom/eleSeedE5x5", "eleSeedE2x5Bottom", "", 'F');
    loader->AddVariable("eleSeedE2x5Left/eleSeedE5x5", "eleSeedE2x5Left", "", 'F');
    loader->AddVariable("eleSeedE2x5Right/eleSeedE5x5", "eleSeedE2x5Right", "", 'F');
    loader->AddVariable("NEcalClusters", "NEcalClusters", "", 'I');
    loader->AddVariable("eleSeedCryIeta", "eleSeedCryIeta", "", 'F');
    loader->AddVariable("eleSeedCryIphi", "eleSeedCryIphi", "", 'F');

    if (options)
        loader->AddVariable("eleESOverRaw", "eleESOverRaw", "", 'F');

    loader->AddSpectator("hiBin", "centrality", "", 'I');

    loader->AddTarget("eleRefE / (eleSCRawEn + eleESEn)");

    TFile* f = TFile::Open(input, "read");
    TTree* t = (TTree*)f->Get("electrons");

    loader->AddRegressionTree(t, 1.);

    /* loader->SetWeightExpression("", ""); */

    TCut match = "eleGenMatchIndex!=-1";
    TCut fiducial = "elePt>15 && eleTrkPt>5 && mcE>15";
    TCut barrel = "abs(eleSCEta)<1.442";
    TCut endcap = "abs(eleSCEta)>1.566 && abs(eleSCEta)<2.5";

    TCut total = match && fiducial && (options ? endcap : barrel);

    int nevents = t->Draw("elePt", total, "goff") * 0.7;
    loader->PrepareTrainingAndTestTree(
        total,
        Form("nTrain_Regression=%i:nTest_Regression=0:"
             "SplitMode=Random:NormMode=NumEvents:!V",
             nevents));

    if (use["BDT"])
        factory->BookMethod(
            loader,
            TMVA::Types::kBDT,
            "BDT",
            "!H:!V:NTrees=1000:MinNodeSize=1.0%:"
            "BoostType=AdaBoostR2:SeparationType=RegressionVariance:"
            "nCuts=20:PruneMethod=CostComplexity:PruneStrength=30");

    if (use["BDTG"])
        factory->BookMethod(
            loader,
            TMVA::Types::kBDT,
            "BDTG",
            "!H:!V:NTrees=2000::"
            "BoostType=Grad:Shrinkage=0.1:"
            "UseBaggedBoost:BaggedSampleFraction=0.5:"
            "nCuts=20:MaxDepth=4");

    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();

    fout->Close();

    delete factory;
    delete loader;

    TMVA::deviations(tag, output, TMVA::kMVAType, kTRUE);
    TMVA::deviations(tag, output, TMVA::kCompareType, kTRUE);
    TMVA::deviations(tag, output, TMVA::kMVAType, kFALSE);
    TMVA::deviations(tag, output, TMVA::kCompareType, kFALSE);
    TMVA::regression_averagedevs(tag, output);

    return 0;
}

int apply(char const* tag, char const* input, char const* output, int options,
          char const* method) {
    TMVA::Tools::Instance();

    TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

    float r_eleESOverRaw;
    float r_eleSeedCryIphi;
    float r_eleSeedCryIeta;
    float r_NEcalClusters;
    float r_eleSeedE2x5Right;
    float r_eleSeedE2x5Left;
    float r_eleSeedE2x5Bottom;
    float r_eleSeedE2x5Top;
    float r_eleSeedE2x5Max;
    float r_eleSeedERight;
    float r_eleSeedELeft;
    float r_eleSeedEBottom;
    float r_eleSeedETop;
    float r_eleSeedE2nd;
    float r_eleSeedEMax;
    float r_eleSPP;
    float r_eleSEP;
    float r_eleSEE;
    float r_eledPhiSCSeed;
    float r_eledEtaSCSeed;
    /* float r_rho; */
    float r_eleHoverE;
    float r_eleSeedE3x3;
    float r_eleSeedE5x5;
    float r_eleSeedEn;
    float r_eleSCPhiWidth;
    float r_eleSCEtaWidth;
    float r_eleSCRawEn;

    reader->AddVariable("eleSCRawEn", &r_eleSCRawEn);
    reader->AddVariable("eleSCEtaWidth", &r_eleSCEtaWidth);
    reader->AddVariable("eleSCPhiWidth", &r_eleSCPhiWidth);
    reader->AddVariable("eleSeedEn/eleSCRawEn", &r_eleSeedEn);
    reader->AddVariable("eleSeedE5x5/eleSCRawEn", &r_eleSeedE5x5);
    reader->AddVariable("eleSeedE3x3/eleSCRawEn", &r_eleSeedE3x3);
    reader->AddVariable("eleHoverE", &r_eleHoverE);
    /* reader->AddVariable("rho", &r_rho); */
    reader->AddVariable("eledEtaSCSeed", &r_eledEtaSCSeed);
    reader->AddVariable("eledPhiSCSeed", &r_eledPhiSCSeed);
    reader->AddVariable("eleSEE", &r_eleSEE);
    reader->AddVariable("eleSEP", &r_eleSEP);
    reader->AddVariable("eleSPP", &r_eleSPP);
    reader->AddVariable("eleSeedEMax/eleSeedE5x5", &r_eleSeedEMax);
    reader->AddVariable("eleSeedE2nd/eleSeedE5x5", &r_eleSeedE2nd);
    reader->AddVariable("eleSeedETop/eleSeedE5x5", &r_eleSeedETop);
    reader->AddVariable("eleSeedEBottom/eleSeedE5x5", &r_eleSeedEBottom);
    reader->AddVariable("eleSeedELeft/eleSeedE5x5", &r_eleSeedELeft);
    reader->AddVariable("eleSeedERight/eleSeedE5x5", &r_eleSeedERight);
    reader->AddVariable("eleSeedE2x5Max/eleSeedE5x5", &r_eleSeedE2x5Max);
    reader->AddVariable("eleSeedE2x5Top/eleSeedE5x5", &r_eleSeedE2x5Top);
    reader->AddVariable("eleSeedE2x5Bottom/eleSeedE5x5", &r_eleSeedE2x5Bottom);
    reader->AddVariable("eleSeedE2x5Left/eleSeedE5x5", &r_eleSeedE2x5Left);
    reader->AddVariable("eleSeedE2x5Right/eleSeedE5x5", &r_eleSeedE2x5Right);
    reader->AddVariable("NEcalClusters", &r_NEcalClusters);
    reader->AddVariable("eleSeedCryIeta", &r_eleSeedCryIeta);
    reader->AddVariable("eleSeedCryIphi", &r_eleSeedCryIphi);

    if (options)
        reader->AddVariable("eleESOverRaw", &r_eleESOverRaw);

    int hiBin;

    reader->AddSpectator("hiBin", &hiBin);

    reader->BookMVA(
        method, Form("%s/weights/TMVARegression_%s.weights.xml", tag, method));

    TFile* f = new TFile(input, "read");
    TTree* t = (TTree*)f->Get("electrons");

    int nEle;
    std::vector<float>* elePt = 0;
    std::vector<float>* eleTrkPt = 0;
    std::vector<float>* eleSCEta = 0;
    std::vector<float>* eleSCEn = 0;
    std::vector<float>* eleESEn = 0;
    std::vector<int>* eleGenMatchIndex = 0;
    std::vector<float>* mcE = 0;

    std::vector<float>* eleESOverRaw = 0;
    std::vector<float>* eleSeedCryIphi = 0;
    std::vector<float>* eleSeedCryIeta = 0;
    std::vector<int>* NEcalClusters = 0;
    std::vector<float>* eleSeedE2x5Right = 0;
    std::vector<float>* eleSeedE2x5Left = 0;
    std::vector<float>* eleSeedE2x5Bottom = 0;
    std::vector<float>* eleSeedE2x5Top = 0;
    std::vector<float>* eleSeedE2x5Max = 0;
    std::vector<float>* eleSeedERight = 0;
    std::vector<float>* eleSeedELeft = 0;
    std::vector<float>* eleSeedEBottom = 0;
    std::vector<float>* eleSeedETop = 0;
    std::vector<float>* eleSeedE2nd = 0;
    std::vector<float>* eleSeedEMax = 0;
    std::vector<float>* eleSPP = 0;
    std::vector<float>* eleSEP = 0;
    std::vector<float>* eleSEE = 0;
    std::vector<float>* eledPhiSCSeed = 0;
    std::vector<float>* eledEtaSCSeed = 0;
    /* float rho = -1; */
    std::vector<float>* eleHoverE = 0;
    std::vector<float>* eleSeedE3x3 = 0;
    std::vector<float>* eleSeedE5x5 = 0;
    std::vector<float>* eleSeedEn = 0;
    std::vector<float>* eleSCPhiWidth = 0;
    std::vector<float>* eleSCEtaWidth = 0;
    std::vector<float>* eleSCRawEn = 0;

    t->SetBranchAddress("nEle", &nEle);
    t->SetBranchAddress("elePt", &elePt);
    t->SetBranchAddress("eleTrkPt", &eleTrkPt);
    t->SetBranchAddress("eleSCEta", &eleSCEta);
    t->SetBranchAddress("eleSCEn", &eleSCEn);
    t->SetBranchAddress("eleESEn", &eleESEn);
    t->SetBranchAddress("eleGenMatchIndex", &eleGenMatchIndex);
    t->SetBranchAddress("mcE", &mcE);

    t->SetBranchAddress("eleSCRawEn", &eleSCRawEn);
    t->SetBranchAddress("eleSCEtaWidth", &eleSCEtaWidth);
    t->SetBranchAddress("eleSCPhiWidth", &eleSCPhiWidth);
    t->SetBranchAddress("eleSeedEn", &eleSeedEn);
    t->SetBranchAddress("eleSeedE5x5", &eleSeedE5x5);
    t->SetBranchAddress("eleSeedE3x3", &eleSeedE3x3);
    t->SetBranchAddress("eleHoverE", &eleHoverE);
    /* t->SetBranchAddress("rho", &rho); */
    t->SetBranchAddress("eledEtaSCSeed", &eledEtaSCSeed);
    t->SetBranchAddress("eledPhiSCSeed", &eledPhiSCSeed);
    t->SetBranchAddress("eleSEE", &eleSEE);
    t->SetBranchAddress("eleSEP", &eleSEP);
    t->SetBranchAddress("eleSPP", &eleSPP);
    t->SetBranchAddress("eleSeedEMax", &eleSeedEMax);
    t->SetBranchAddress("eleSeedE2nd", &eleSeedE2nd);
    t->SetBranchAddress("eleSeedETop", &eleSeedETop);
    t->SetBranchAddress("eleSeedEBottom", &eleSeedEBottom);
    t->SetBranchAddress("eleSeedELeft", &eleSeedELeft);
    t->SetBranchAddress("eleSeedERight", &eleSeedERight);
    t->SetBranchAddress("eleSeedE2x5Max", &eleSeedE2x5Max);
    t->SetBranchAddress("eleSeedE2x5Top", &eleSeedE2x5Top);
    t->SetBranchAddress("eleSeedE2x5Bottom", &eleSeedE2x5Bottom);
    t->SetBranchAddress("eleSeedE2x5Left", &eleSeedE2x5Left);
    t->SetBranchAddress("eleSeedE2x5Right", &eleSeedE2x5Right);
    t->SetBranchAddress("NEcalClusters", &NEcalClusters);
    t->SetBranchAddress("eleSeedCryIeta", &eleSeedCryIeta);
    t->SetBranchAddress("eleSeedCryIphi", &eleSeedCryIphi);
    t->SetBranchAddress("eleESOverRaw", &eleESOverRaw);

    t->SetBranchAddress("hiBin", &hiBin);

    TFile* fout = new TFile(output, "update");

    char const* title = ";#deltaE/E;";

    std::vector<int> clow  = {   0,  0, 20,  60, 100 };
    std::vector<int> chigh = { 200, 20, 60, 100, 200 };
    int ncent = clow.size();

    TH1F* before[ncent];
    TH1F* after[ncent];

    for (int c = 0; c < ncent; ++c) {
        before[c] = new TH1F(Form("before_%s_c%i", tag, c), title, 100, -0.4, 0.4);
        after[c] = new TH1F(Form("after_%s_c%i", tag, c), title, 100, -0.4, 0.4);

        before[c]->SetStats(0);
        after[c]->SetStats(0);
    }

    int nentries = t->GetEntries();
    for (int i = 0; i < nentries; ++i) {
        t->GetEntry(i);

        if (i % 10000 == 0) { printf("entry: %i/%i\n", i, nentries); }

        for (int j = 0; j < nEle; ++j) {
            if ((*eleGenMatchIndex)[j] < 0) { continue; }
            if ((*elePt)[j] < 15) { continue; }
            if ((*eleTrkPt)[j] < 5) { continue; }
            if (!options && !(std::abs((*eleSCEta)[j]) < 1.442)) { continue; }
            if (options && ((std::abs((*eleSCEta)[j]) < 1.566)
                || (std::abs((*eleSCEta)[j]) > 2.5))) { continue; }

            r_eleSCRawEn = (*eleSCRawEn)[j];
            r_eleSCEtaWidth = (*eleSCEtaWidth)[j];
            r_eleSCPhiWidth = (*eleSCPhiWidth)[j];
            r_eleSeedEn = (*eleSeedEn)[j] / (*eleSCRawEn)[j];
            r_eleSeedE5x5 = (*eleSeedE5x5)[j] / (*eleSCRawEn)[j];
            r_eleSeedE3x3 = (*eleSeedE3x3)[j] / (*eleSCRawEn)[j];
            r_eleHoverE = (*eleHoverE)[j];
            /* r_rho = rho; */
            r_eledEtaSCSeed = (*eledEtaSCSeed)[j];
            r_eledPhiSCSeed = (*eledPhiSCSeed)[j];
            r_eleSEE = (*eleSEE)[j];
            r_eleSEP = (*eleSEP)[j];
            r_eleSPP = (*eleSPP)[j];
            r_eleSeedEMax = (*eleSeedEMax)[j] / (*eleSeedE5x5)[j];
            r_eleSeedE2nd = (*eleSeedE2nd)[j] / (*eleSeedE5x5)[j];
            r_eleSeedETop = (*eleSeedETop)[j] / (*eleSeedE5x5)[j];
            r_eleSeedEBottom = (*eleSeedEBottom)[j] / (*eleSeedE5x5)[j];
            r_eleSeedELeft = (*eleSeedELeft)[j] / (*eleSeedE5x5)[j];
            r_eleSeedERight = (*eleSeedERight)[j] / (*eleSeedE5x5)[j];
            r_eleSeedE2x5Max = (*eleSeedE2x5Max)[j] / (*eleSeedE5x5)[j];
            r_eleSeedE2x5Top = (*eleSeedE2x5Top)[j] / (*eleSeedE5x5)[j];
            r_eleSeedE2x5Bottom = (*eleSeedE2x5Bottom)[j] / (*eleSeedE5x5)[j];
            r_eleSeedE2x5Left = (*eleSeedE2x5Left)[j] / (*eleSeedE5x5)[j];
            r_eleSeedE2x5Right = (*eleSeedE2x5Right)[j] / (*eleSeedE5x5)[j];
            r_NEcalClusters = (*NEcalClusters)[j];
            r_eleSeedCryIeta = (*eleSeedCryIeta)[j];
            r_eleSeedCryIphi = (*eleSeedCryIphi)[j];
            r_eleESOverRaw = (*eleESOverRaw)[j];

            float val = (reader->EvaluateRegression(Form("%s", method)))[0];
            float truth = (*mcE)[(*eleGenMatchIndex)[j]];

            for (int c = 0; c < ncent; ++c) {
                if (hiBin >= clow[c] && hiBin < chigh[c]) {
                    before[c]->Fill(((*eleSCEn)[j] - truth) / truth);
                    after[c]->Fill(
                        (val * ((*eleSCRawEn)[j] + (*eleESEn)[j]) - truth)
                        / truth);
                }
            }
        }
    }

    TCanvas* c1 = new TCanvas("c1", "", 400, 400);

    {
        before[0]->Scale(1. / before[0]->GetEntries());
        after[0]->Scale(1. / after[0]->GetEntries());

        hstyle(before[0], 21, 46, 0.6);
        hstyle(after[0], 21, 38, 0.6);

        before[0]->SetAxisRange(
            0,
            1.2 * std::max(before[0]->GetMaximum(), after[0]->GetMaximum()),
            "Y");

        before[0]->Draw("p same");
        after[0]->Draw("p same");

        TLegend* l1 = new TLegend(0.6, 0.72, 0.8, 0.84);
        lstyle(l1, 43, 12);
        l1->AddEntry(before[0], "default", "pl");
        l1->AddEntry(after[0], "applying regression", "pl");
        l1->Draw();

        c1->SaveAs(Form("resolution_%s_%s.png", options ? "endcap" : "barrel",
                        method));
    }

    {
        for (int c = 1; c < 5; ++c) {
            before[c]->Scale(1. / before[c]->GetEntries());
            after[c]->Scale(1. / after[c]->GetEntries());
        }

        hstyle(before[1], 20, 46, 0.6);
        hstyle(after[1], 24, 38, 0.6);
        hstyle(before[2], 21, 46, 0.6);
        hstyle(after[2], 25, 38, 0.6);
        hstyle(before[3], 22, 46, 0.6);
        hstyle(after[3], 26, 38, 0.6);
        hstyle(before[4], 23, 46, 0.6);
        hstyle(after[4], 32, 38, 0.6);

        before[1]->SetAxisRange(
            0,
            1.2 * std::max(before[1]->GetMaximum(), after[1]->GetMaximum()),
            "Y");

        before[1]->Draw("p same");
        after[1]->Draw("p same");
        before[2]->Draw("p same");
        after[2]->Draw("p same");
        before[3]->Draw("p same");
        after[3]->Draw("p same");
        before[4]->Draw("p same");
        after[4]->Draw("p same");

        TLegend* l1 = new TLegend(0.6, 0.72, 0.8, 0.84);
        lstyle(l1, 43, 12);
        l1->AddEntry(before[1], "default", "pl");
        l1->AddEntry(after[1], "applying regression", "pl");
        l1->Draw();

        c1->SaveAs(Form("resolution_centrality_%s_%s.png",
                        options ? "endcap" : "barrel",
                        method));
    }

    fout->Write("", TObject::kOverwrite);
    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 5)
        return regress(argv[1], argv[2], argv[3], atoi(argv[4]));

    if (argc == 6)
        return apply(argv[1], argv[2], argv[3], atoi(argv[4]), argv[5]);

    printf("usage: %s [tag] [input] [output] [options] [method]\n", argv[0]);
    printf("  options: [0] barrel | [1] endcap\n");

    return 1;
}
