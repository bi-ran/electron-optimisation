#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLorentzVector.h"

#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Reader.h"

#include "TMVA/deviations.h"
#include "TMVA/regression_averagedevs.h"

#include "include/electrontree.h"

#include <map>
#include <string>
#include <vector>

#include "include/cosmetics.h"

int passes_looseid_barrel(electrontree* t, int j) {
    return (fabs((*t->eleSCEta)[j]) < 1.4442
        && ((t->hiBin < 60
            && (*t->eleSigmaIEtaIEta_2012)[j] < 0.0161
            && fabs((*t->eledEtaAtVtx)[j]) < 0.0053
            && fabs((*t->eledPhiAtVtx)[j]) < 0.0288
            && (*t->eleHoverE)[j] < 0.1984
            && (*t->eleEoverPInv)[j] < 0.1129)
        || (t->hiBin > 59
            && (*t->eleSigmaIEtaIEta_2012)[j] < 0.0117
            && fabs((*t->eledEtaAtVtx)[j]) < 0.0071
            && fabs((*t->eledPhiAtVtx)[j]) < 0.0221
            && (*t->eleHoverE)[j] < 0.1892
            && (*t->eleEoverPInv)[j] < 0.0405))
        && fabs((*t->eleD0)[j]) < 0.01
        && fabs((*t->eleDz)[j]) < 0.04
        && (*t->eleMissHits)[j] < 2);
}

int passes_looseid_endcap(electrontree* t, int j) {
    return (fabs((*t->eleSCEta)[j]) > 1.566
        && fabs((*t->eleSCEta)[j]) < 2.5
        && ((t->hiBin < 60
            && (*t->eleSigmaIEtaIEta_2012)[j] < 0.0479
            && fabs((*t->eledEtaAtVtx)[j]) < 0.0145
            && fabs((*t->eledPhiAtVtx)[j]) < 0.0516
            && (*t->eleHoverE)[j] < 0.1910
            && (*t->eleEoverPInv)[j] < 0.0115)
        || (t->hiBin > 59
            && (*t->eleSigmaIEtaIEta_2012)[j] < 0.0447
            && fabs((*t->eledEtaAtVtx)[j]) < 0.0108
            && fabs((*t->eledPhiAtVtx)[j]) < 0.0301
            && (*t->eleHoverE)[j] < 0.1627
            && (*t->eleEoverPInv)[j] < 0.0281))
        && fabs((*t->eleD0)[j]) < 0.02
        && fabs((*t->eleDz)[j]) < 0.04
        && (*t->eleMissHits)[j] < 2);
}

#define VARIABLES(ACTION)                                               \
    ACTION(eleSCRawEn, 1, 0, 400)                                       \
    ACTION(eleSCEtaWidth, 1, 0, 0.04)                                   \
    ACTION(eleSCPhiWidth, 1, 0, 0.1)                                    \
    ACTION(eleSeedEn, eleSCRawEn, 0, 400)                               \
    ACTION(eleSeedE5x5, eleSCRawEn, 0, 400)                             \
    ACTION(eleSeedE3x3, eleSCRawEn, 0, 400)                             \
    ACTION(eleHoverE, 1, 0, 0.1)                                        \
    ACTION(eledEtaSCSeed, 1, 0, 0.02)                                   \
    ACTION(eledPhiSCSeed, 1, 0, 0.02)                                   \
    ACTION(eleSEE, 1, 0, 0.04)                                          \
    ACTION(eleSEP, 1, 0, 0.8)                                           \
    ACTION(eleSPP, 1, 0, 0.1)                                           \
    ACTION(eleSeedEMax, eleSeedE5x5, 0, 300)                            \
    ACTION(eleSeedE2nd, eleSeedE5x5, 0, 100)                            \
    ACTION(eleSeedETop, eleSeedE5x5, 0, 100)                            \
    ACTION(eleSeedEBottom, eleSeedE5x5, 0, 100)                         \
    ACTION(eleSeedELeft, eleSeedE5x5, 0, 100)                           \
    ACTION(eleSeedERight, eleSeedE5x5, 0, 100)                          \
    ACTION(eleSeedE2x5Max, eleSeedE5x5, 0, 400)                         \
    ACTION(eleSeedE2x5Top, eleSeedE5x5, 0, 100)                         \
    ACTION(eleSeedE2x5Bottom, eleSeedE5x5, 0, 100)                      \
    ACTION(eleSeedE2x5Left, eleSeedE5x5, 0, 100)                        \
    ACTION(eleSeedE2x5Right, eleSeedE5x5, 0, 100)                       \
    ACTION(NEcalClusters, 1, 0, 15)                                     \
    ACTION(eleSeedCryIeta, 1, -50, 50)                                  \
    ACTION(eleSeedCryIphi, 1, -50, 50)                                  \

#define BACKUP(ACTION)                                                  \
    ACTION(eleTrkQoverPMode, 1, -10, 10)                                \
    ACTION(eleTrkPtMode, 1, 0, 100)                                     \
    ACTION(eleTrkQoverPModeErr, 1, 0, 10)                               \
    ACTION(eleTrkPtModeErr, 1, 0, 10)                                   \
    ACTION(eleBrem, 1, 0, 1)                                            \
    ACTION(eleTrkLayers, 1, 0, 40)                                      \
    ACTION(eleTrkValidHits, 1, 0, 40)                                   \
    ACTION(eleTrkNormalizedChi2, 1, 0, 100)

int compare(char const* data, char const* sim) {
    TFile* fs = new TFile(sim, "read");
    TTree* ts = (TTree*)fs->Get("electrons");

#define BOOK(var, denom, min, max)                                      \
    TH1F* hdb##var = new TH1F("hdb" #var, "", 100, min, max);           \
    TH1F* hde##var = new TH1F("hde" #var, "", 100, min, max);           \
    TH1F* hsb##var = new TH1F("hsb" #var, "", 100, min, max);           \
    TH1F* hse##var = new TH1F("hse" #var, "", 100, min, max);

    VARIABLES(BOOK)

    TCut match = "eleGenMatchIndex!=-1";
    TCut fiducial = "elePt>15 && eleTrkPt>5 && mcE>15";
    TCut barrel = "abs(eleSCEta)<1.442";
    TCut endcap = "abs(eleSCEta)>1.566 && abs(eleSCEta)<2.5";

    TCut totalb = match && fiducial && barrel;
    TCut totale = match && fiducial && endcap;

#define DRAW(var, denom, ...)                                           \
    ts->Draw(#var "/" #denom ">>hsb" #var, totalb * "ncoll", "goff");   \
    ts->Draw(#var "/" #denom ">>hse" #var, totale * "ncoll", "goff");

    VARIABLES(DRAW)

    TFile* fd = new TFile(data, "read");
    TTree* td = (TTree*)fd->Get("electrons");
    auto ed = new electrontree(0, 0, 0, td);

    int64_t nentries = td->GetEntries();
    for (int64_t i = 0; i < nentries; ++i) {
        td->GetEntry(i);

        for (int j = 0; j < ed->nEle; ++j) {
            bool is_barrel = passes_looseid_barrel(ed, j);
            bool is_endcap = passes_looseid_endcap(ed, j);

            if (!is_barrel && !is_endcap) { continue; }

#define FILL(var, ...) {                                                \
    TH1F* h = is_barrel ? hdb##var : hde##var;                          \
    h->Fill((*ed->var)[j]); }

            VARIABLES(FILL);
        }
    }

    TCanvas* c1 = new TCanvas("c1", "", 400, 400);

#define PLOT(var, ...)                                                  \
    {                                                                   \
        hstyle(hdb##var, 20, 40, 0.6);                                  \
        hstyle(hde##var, 21, 46, 0.6);                                  \
        hstyle(hsb##var, 24, 40, 0.6);                                  \
        hstyle(hse##var, 25, 46, 0.6);                                  \
                                                                        \
        hdb##var->Scale(1. / hdb##var->Integral());                     \
        hde##var->Scale(1. / hde##var->Integral());                     \
        hsb##var->Scale(1. / hsb##var->Integral());                     \
        hse##var->Scale(1. / hse##var->Integral());                     \
                                                                        \
        hsb##var->Draw();                                               \
        hdb##var->Draw("same");                                         \
        hde##var->Draw("same");                                         \
        hsb##var->Draw("same");                                         \
        hse##var->Draw("same");                                         \
                                                                        \
        c1->SaveAs("comp_" #var ".pdf");                                \
    }

    VARIABLES(PLOT)

    return 0;
}

int regress(char const* tag, char const* input, char const* output,
            int options) {
    TMVA::Tools::Instance();

    std::map<std::string, bool> use;
    use["BDTG"] = true;

    TFile* fout = TFile::Open(output, "recreate");

    TMVA::Factory* factory = new TMVA::Factory(
        "TMVARegression",
        fout,
        "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression");

    TMVA::DataLoader* loader = new TMVA::DataLoader(tag);

    /* loader->AddVariable("eleSCRawEn", "eleSCRawEn", "", 'F'); */
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
    /* loader->AddVariable("eleTrkQoverPMode", "eleTrkQoverPMode", "", 'F'); */
    /* loader->AddVariable("eleTrkPtMode", "eleTrkPtMode", "", 'F'); */
    /* loader->AddVariable("eleTrkQoverPModeErr/eleTrkQoverPMode", "eleTrkQoverPModeErr",  "", 'F'); */
    /* loader->AddVariable("eleTrkPtModeErr/eleTrkPtMode", "eleTrkPtModeErr", "", 'F'); */
    /* loader->AddVariable("eleBrem", "eleBrem", "", 'F'); */
    /* loader->AddVariable("eleTrkLayers", "eleTrkLayers", "", 'I'); */
    /* loader->AddVariable("eleTrkValidHits", "eleTrkValidHits", "", 'I'); */
    /* loader->AddVariable("eleTrkNormalizedChi2", "eleTrkNormalizedChi2", "", 'F'); */

    if (options)
        loader->AddVariable("eleESOverRaw", "eleESOverRaw", "", 'F');

    loader->AddSpectator("hiBin", "centrality", "", 'I');

    loader->AddTarget("eleRefE / (eleSCRawEn + eleESEn)");

    TFile* f = TFile::Open(input, "read");
    TTree* t = (TTree*)f->Get("electrons");

    loader->AddRegressionTree(t, 1.);

    loader->SetWeightExpression("ncoll");

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

    /* if (use["BDT"]) */
    /*     factory->BookMethod( */
    /*         loader, */
    /*         TMVA::Types::kBDT, */
    /*         "BDT", */
    /*         "!H:!V:NTrees=1000:MinNodeSize=1.0%:" */
    /*         "BoostType=AdaBoostR2:SeparationType=RegressionVariance:" */
    /*         "nCuts=20:PruneMethod=CostComplexity:PruneStrength=30"); */

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

int validate(char const* tag, char const* input, char const* output,
             int options, char const* method) {
    TMVA::Tools::Instance();

    TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

    float r_eleESOverRaw;
    /* float r_eleTrkNormalizedChi2; */
    /* float r_eleTrkValidHits; */
    /* float r_eleTrkLayers; */
    /* float r_eleBrem; */
    /* float r_eleTrkPtModeErr; */
    /* float r_eleTrkQoverPModeErr; */
    /* float r_eleTrkPtMode; */
    /* float r_eleTrkQoverPMode; */
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
    /* float r_eleSCRawEn; */

    /* reader->AddVariable("eleSCRawEn", &r_eleSCRawEn); */
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
    /* reader->AddVariable("eleTrkQoverPMode", &r_eleTrkQoverPMode); */
    /* reader->AddVariable("eleTrkPtMode", &r_eleTrkPtMode); */
    /* reader->AddVariable("eleTrkQoverPModeErr/eleTrkQoverPMode", &r_eleTrkQoverPModeErr); */
    /* reader->AddVariable("eleTrkPtModeErr/eleTrkPtMode", &r_eleTrkPtModeErr); */
    /* reader->AddVariable("eleBrem", &r_eleBrem); */
    /* reader->AddVariable("eleTrkLayers", &r_eleTrkLayers); */
    /* reader->AddVariable("eleTrkValidHits", &r_eleTrkValidHits); */
    /* reader->AddVariable("eleTrkNormalizedChi2", &r_eleTrkNormalizedChi2); */

    if (options)
        reader->AddVariable("eleESOverRaw", &r_eleESOverRaw);

    int hiBin;

    reader->AddSpectator("hiBin", &hiBin);

    reader->BookMVA(
        method, Form("%s/weights/TMVARegression_%s.weights.xml", tag, method));

    TFile* f = new TFile(input, "read");
    TTree* t = (TTree*)f->Get("electrons");
    auto e = new electrontree(0, 0, 0, t);

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

        for (int j = 0; j < e->nEle; ++j) {
            if ((*e->eleGenMatchIndex)[j] < 0) { continue; }
            if ((*e->elePt)[j] < 15) { continue; }
            if ((*e->eleTrkPt)[j] < 5) { continue; }
            if (!options && !(std::abs((*e->eleSCEta)[j]) < 1.442)) { continue; }
            if (options && ((std::abs((*e->eleSCEta)[j]) < 1.566)
                || (std::abs((*e->eleSCEta)[j]) > 2.5))) { continue; }

            /* r_eleSCRawEn = (*e->eleSCRawEn)[j]; */
            r_eleSCEtaWidth = (*e->eleSCEtaWidth)[j];
            r_eleSCPhiWidth = (*e->eleSCPhiWidth)[j];
            r_eleSeedEn = (*e->eleSeedEn)[j] / (*e->eleSCRawEn)[j];
            r_eleSeedE5x5 = (*e->eleSeedE5x5)[j] / (*e->eleSCRawEn)[j];
            r_eleSeedE3x3 = (*e->eleSeedE3x3)[j] / (*e->eleSCRawEn)[j];
            r_eleHoverE = (*e->eleHoverE)[j];
            /* r_rho = e->rho; */
            r_eledEtaSCSeed = (*e->eledEtaSCSeed)[j];
            r_eledPhiSCSeed = (*e->eledPhiSCSeed)[j];
            r_eleSEE = (*e->eleSEE)[j];
            r_eleSEP = (*e->eleSEP)[j];
            r_eleSPP = (*e->eleSPP)[j];
            r_eleSeedEMax = (*e->eleSeedEMax)[j] / (*e->eleSeedE5x5)[j];
            r_eleSeedE2nd = (*e->eleSeedE2nd)[j] / (*e->eleSeedE5x5)[j];
            r_eleSeedETop = (*e->eleSeedETop)[j] / (*e->eleSeedE5x5)[j];
            r_eleSeedEBottom = (*e->eleSeedEBottom)[j] / (*e->eleSeedE5x5)[j];
            r_eleSeedELeft = (*e->eleSeedELeft)[j] / (*e->eleSeedE5x5)[j];
            r_eleSeedERight = (*e->eleSeedERight)[j] / (*e->eleSeedE5x5)[j];
            r_eleSeedE2x5Max = (*e->eleSeedE2x5Max)[j] / (*e->eleSeedE5x5)[j];
            r_eleSeedE2x5Top = (*e->eleSeedE2x5Top)[j] / (*e->eleSeedE5x5)[j];
            r_eleSeedE2x5Bottom = (*e->eleSeedE2x5Bottom)[j] / (*e->eleSeedE5x5)[j];
            r_eleSeedE2x5Left = (*e->eleSeedE2x5Left)[j] / (*e->eleSeedE5x5)[j];
            r_eleSeedE2x5Right = (*e->eleSeedE2x5Right)[j] / (*e->eleSeedE5x5)[j];
            r_NEcalClusters = (*e->NEcalClusters)[j];
            r_eleSeedCryIeta = (*e->eleSeedCryIeta)[j];
            r_eleSeedCryIphi = (*e->eleSeedCryIphi)[j];
            r_eleESOverRaw = (*e->eleESOverRaw)[j];
            /* r_eleTrkQoverPMode = (*e->eleTrkQoverPMode)[j]; */
            /* r_eleTrkPtMode = (*e->eleTrkPtMode)[j]; */
            /* r_eleTrkQoverPModeErr = (*e->eleTrkQoverPModeErr)[j] / (*e->eleTrkQoverPMode)[j]; */
            /* r_eleTrkPtModeErr = (*e->eleTrkPtModeErr)[j] / (*e->eleTrkPtMode)[j]; */
            /* r_eleBrem = (*e->eleBrem)[j]; */
            /* r_eleTrkLayers = (*e->eleTrkLayers)[j]; */
            /* r_eleTrkValidHits = (*e->eleTrkValidHits)[j]; */
            /* r_eleTrkNormalizedChi2 = (*e->eleTrkNormalizedChi2)[j]; */

            float val = (reader->EvaluateRegression(Form("%s", method)))[0];
            float truth = (*e->mcE)[(*e->eleGenMatchIndex)[j]];

            for (int c = 0; c < ncent; ++c) {
                if (hiBin >= clow[c] && hiBin < chigh[c]) {
                    before[c]->Fill(((*e->eleSCEn)[j] - truth) / truth);
                    after[c]->Fill(
                        (val * ((*e->eleSCRawEn)[j] + (*e->eleESEn)[j]) - truth)
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

        before[0]->Draw("p");
        after[0]->Draw("p same");

        TLegend* l1 = new TLegend(0.6, 0.72, 0.8, 0.84);
        lstyle(l1, 43, 12);
        l1->AddEntry(before[0], "default", "pl");
        l1->AddEntry(after[0], "applying regression", "pl");
        l1->Draw();

        c1->SaveAs(Form("resolution_%s_%s.pdf", options ? "endcap" : "barrel",
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

        before[1]->Draw("p");
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

        c1->SaveAs(Form("resolution_centrality_%s_%s.pdf",
                        options ? "endcap" : "barrel",
                        method));
    }

    fout->Write("", TObject::kOverwrite);
    fout->Close();

    return 0;
}

int apply(char const* tag, char const* input, char const* method) {
    char modtag[100];
    strcpy(modtag, tag);
    if (strstr(input, "HIHardProbes") != NULL)
        strcat(modtag, "data");

    TMVA::Tools::Instance();

    TMVA::Reader* breader = new TMVA::Reader("!Color:!Silent");
    TMVA::Reader* ereader = new TMVA::Reader("!Color:!Silent");

    float r_eleESOverRaw;
    /* float r_eleTrkNormalizedChi2; */
    /* float r_eleTrkValidHits; */
    /* float r_eleTrkLayers; */
    /* float r_eleBrem; */
    /* float r_eleTrkPtModeErr; */
    /* float r_eleTrkQoverPModeErr; */
    /* float r_eleTrkPtMode; */
    /* float r_eleTrkQoverPMode; */
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
    /* float r_eleSCRawEn; */

#define READER_ADD_VAR(reader)                                                          \
    reader->AddVariable("eleSCEtaWidth", &r_eleSCEtaWidth);                             \
    reader->AddVariable("eleSCPhiWidth", &r_eleSCPhiWidth);                             \
    reader->AddVariable("eleSeedEn/eleSCRawEn", &r_eleSeedEn);                          \
    reader->AddVariable("eleSeedE5x5/eleSCRawEn", &r_eleSeedE5x5);                      \
    reader->AddVariable("eleSeedE3x3/eleSCRawEn", &r_eleSeedE3x3);                      \
    reader->AddVariable("eleHoverE", &r_eleHoverE);                                     \
    /* reader->AddVariable("rho", &r_rho); */                                           \
    reader->AddVariable("eledEtaSCSeed", &r_eledEtaSCSeed);                             \
    reader->AddVariable("eledPhiSCSeed", &r_eledPhiSCSeed);                             \
    reader->AddVariable("eleSEE", &r_eleSEE);                                           \
    reader->AddVariable("eleSEP", &r_eleSEP);                                           \
    reader->AddVariable("eleSPP", &r_eleSPP);                                           \
    reader->AddVariable("eleSeedEMax/eleSeedE5x5", &r_eleSeedEMax);                     \
    reader->AddVariable("eleSeedE2nd/eleSeedE5x5", &r_eleSeedE2nd);                     \
    reader->AddVariable("eleSeedETop/eleSeedE5x5", &r_eleSeedETop);                     \
    reader->AddVariable("eleSeedEBottom/eleSeedE5x5", &r_eleSeedEBottom);               \
    reader->AddVariable("eleSeedELeft/eleSeedE5x5", &r_eleSeedELeft);                   \
    reader->AddVariable("eleSeedERight/eleSeedE5x5", &r_eleSeedERight);                 \
    reader->AddVariable("eleSeedE2x5Max/eleSeedE5x5", &r_eleSeedE2x5Max);               \
    reader->AddVariable("eleSeedE2x5Top/eleSeedE5x5", &r_eleSeedE2x5Top);               \
    reader->AddVariable("eleSeedE2x5Bottom/eleSeedE5x5", &r_eleSeedE2x5Bottom);         \
    reader->AddVariable("eleSeedE2x5Left/eleSeedE5x5", &r_eleSeedE2x5Left);             \
    reader->AddVariable("eleSeedE2x5Right/eleSeedE5x5", &r_eleSeedE2x5Right);           \
    reader->AddVariable("NEcalClusters", &r_NEcalClusters);                             \
    reader->AddVariable("eleSeedCryIeta", &r_eleSeedCryIeta);                           \
    reader->AddVariable("eleSeedCryIphi", &r_eleSeedCryIphi);                           \

#define TMP                                                                             \
    reader->AddVariable("eleSCRawEn", &r_eleSCRawEn);                                   \
    reader->AddVariable("eleTrkQoverPMode", &r_eleTrkQoverPMode);                       \
    reader->AddVariable("eleTrkPtMode", &r_eleTrkPtMode);                               \
    reader->AddVariable("eleTrkQoverPModeErr/eleTrkQoverPMode", &r_eleTrkQoverPModeErr);\
    reader->AddVariable("eleTrkPtModeErr/eleTrkPtMode", &r_eleTrkPtModeErr);            \
    reader->AddVariable("eleBrem", &r_eleBrem);                                         \
    reader->AddVariable("eleTrkLayers", &r_eleTrkLayers);                               \
    reader->AddVariable("eleTrkValidHits", &r_eleTrkValidHits);                         \
    reader->AddVariable("eleTrkNormalizedChi2", &r_eleTrkNormalizedChi2);

    READER_ADD_VAR(breader)
    READER_ADD_VAR(ereader)
    ereader->AddVariable("eleESOverRaw", &r_eleESOverRaw);

    int hiBin;

    breader->AddSpectator("hiBin", &hiBin);
    ereader->AddSpectator("hiBin", &hiBin);

    breader->BookMVA(
        method, Form("%s_barrel/weights/TMVARegression_%s.weights.xml", tag, method));
    ereader->BookMVA(
        method, Form("%s_endcap/weights/TMVARegression_%s.weights.xml", tag, method));

    TFile* f = new TFile(input, "read");
    TTree* t = (TTree*)f->Get("electrons");
    auto e = new electrontree(0, 0, 0, t);

    std::vector<int> clow  = {   0,  0, 20,  60, 100 };
    std::vector<int> chigh = { 200, 20, 60, 100, 200 };
    int ncent = clow.size();

    char const* label[2] = { "before", "after" };
    int marker[ncent][2] = {
        { 25, 21 }, { 24, 20 }, { 25, 21 }, { 26, 22 }, { 32, 23 } };

    TH1F* hminv[3][2][ncent][2];
    for (int c = 0; c < ncent; ++c) {
        for (int ab = 0; ab < 2; ++ab) {
            hminv[0][0][c][ab] = new TH1F(Form("hbb_os_c%i_%s", c, label[ab]), "", 30, 60, 120);
            hminv[0][1][c][ab] = new TH1F(Form("hbb_ss_c%i_%s", c, label[ab]), "", 30, 60, 120);
            hminv[1][0][c][ab] = new TH1F(Form("heb_os_c%i_%s", c, label[ab]), "", 30, 60, 120);
            hminv[1][1][c][ab] = new TH1F(Form("heb_ss_c%i_%s", c, label[ab]), "", 30, 60, 120);
            hminv[2][0][c][ab] = new TH1F(Form("hee_os_c%i_%s", c, label[ab]), "", 30, 60, 120);
            hminv[2][1][c][ab] = new TH1F(Form("hee_ss_c%i_%s", c, label[ab]), "", 30, 60, 120);

            hstyle(hminv[0][0][c][ab], marker[c][ab], 46, 0.6);
            hstyle(hminv[1][0][c][ab], marker[c][ab], 40, 0.6);
            hstyle(hminv[2][0][c][ab], marker[c][ab], 38, 0.6);
        }
    }

    int nentries = t->GetEntries();
    for (int i = 0; i < nentries; ++i) {
        t->GetEntry(i);

        if (i % 10000 == 0) { printf("entry: %i/%i\n", i, nentries); }

        std::vector<float> corrected(e->nEle, 0);

        for (int j = 0; j < e->nEle; ++j) {
            if ((*e->elePt)[j] < 15) { continue; }
            if ((*e->eleTrkPt)[j] < 5) { continue; }

            float eta = std::abs((*e->eleSCEta)[j]);
            if ((eta > 1.442 && eta < 1.566) || eta > 2.5) { continue; }

            /* r_eleSCRawEn = (*e->eleSCRawEn)[j]; */
            r_eleSCEtaWidth = (*e->eleSCEtaWidth)[j];
            r_eleSCPhiWidth = (*e->eleSCPhiWidth)[j];
            r_eleSeedEn = (*e->eleSeedEn)[j] / (*e->eleSCRawEn)[j];
            r_eleSeedE5x5 = (*e->eleSeedE5x5)[j] / (*e->eleSCRawEn)[j];
            r_eleSeedE3x3 = (*e->eleSeedE3x3)[j] / (*e->eleSCRawEn)[j];
            r_eleHoverE = (*e->eleHoverE)[j];
            /* r_rho = e->rho; */
            r_eledEtaSCSeed = (*e->eledEtaSCSeed)[j];
            r_eledPhiSCSeed = (*e->eledPhiSCSeed)[j];
            r_eleSEE = (*e->eleSEE)[j];
            r_eleSEP = (*e->eleSEP)[j];
            r_eleSPP = (*e->eleSPP)[j];
            r_eleSeedEMax = (*e->eleSeedEMax)[j] / (*e->eleSeedE5x5)[j];
            r_eleSeedE2nd = (*e->eleSeedE2nd)[j] / (*e->eleSeedE5x5)[j];
            r_eleSeedETop = (*e->eleSeedETop)[j] / (*e->eleSeedE5x5)[j];
            r_eleSeedEBottom = (*e->eleSeedEBottom)[j] / (*e->eleSeedE5x5)[j];
            r_eleSeedELeft = (*e->eleSeedELeft)[j] / (*e->eleSeedE5x5)[j];
            r_eleSeedERight = (*e->eleSeedERight)[j] / (*e->eleSeedE5x5)[j];
            r_eleSeedE2x5Max = (*e->eleSeedE2x5Max)[j] / (*e->eleSeedE5x5)[j];
            r_eleSeedE2x5Top = (*e->eleSeedE2x5Top)[j] / (*e->eleSeedE5x5)[j];
            r_eleSeedE2x5Bottom = (*e->eleSeedE2x5Bottom)[j] / (*e->eleSeedE5x5)[j];
            r_eleSeedE2x5Left = (*e->eleSeedE2x5Left)[j] / (*e->eleSeedE5x5)[j];
            r_eleSeedE2x5Right = (*e->eleSeedE2x5Right)[j] / (*e->eleSeedE5x5)[j];
            r_NEcalClusters = (*e->NEcalClusters)[j];
            r_eleSeedCryIeta = (*e->eleSeedCryIeta)[j];
            r_eleSeedCryIphi = (*e->eleSeedCryIphi)[j];
            r_eleESOverRaw = (*e->eleESOverRaw)[j];
            /* r_eleTrkQoverPMode = (*e->eleTrkQoverPMode)[j]; */
            /* r_eleTrkPtMode = (*e->eleTrkPtMode)[j]; */
            /* r_eleTrkQoverPModeErr = (*e->eleTrkQoverPModeErr)[j] / (*e->eleTrkQoverPMode)[j]; */
            /* r_eleTrkPtModeErr = (*e->eleTrkPtModeErr)[j] / (*e->eleTrkPtMode)[j]; */
            /* r_eleBrem = (*e->eleBrem)[j]; */
            /* r_eleTrkLayers = (*e->eleTrkLayers)[j]; */
            /* r_eleTrkValidHits = (*e->eleTrkValidHits)[j]; */
            /* r_eleTrkNormalizedChi2 = (*e->eleTrkNormalizedChi2)[j]; */

            TMVA::Reader* reader = eta < 1.442 ? breader : ereader;
            float val = (reader->EvaluateRegression(Form("%s", method)))[0];

            float energy = val * ((*e->eleSCRawEn)[j] + (*e->eleESEn)[j]);
            corrected[j] = energy / std::cosh(eta);
        }

        for (int j=0; j<e->nEle; ++j) {
            if ((*e->elePt)[j] < 15) { continue; }
            if ((*e->eleTrkPt)[j] < 5) { continue; }

            int is_1_barrel = passes_looseid_barrel(e, j);
            int is_1_endcap = passes_looseid_endcap(e, j);

            if (!is_1_barrel && !is_1_endcap)
                continue;

            /* double electron invariant mass */
            for (int k=j+1; k<e->nEle; ++k) {
                if ((*e->elePt)[k] < 15) { continue; }
                if ((*e->eleTrkPt)[k] < 5) { continue; }

                int is_2_barrel = passes_looseid_barrel(e, k);
                int is_2_endcap = passes_looseid_endcap(e, k);

                if (!is_2_barrel && !is_2_endcap)
                    continue;

                TLorentzVector e1; TLorentzVector e2;
                e1.SetPtEtaPhiM(corrected[j], (*e->eleEta)[j], (*e->elePhi)[j], 0.000511);
                e2.SetPtEtaPhiM(corrected[k], (*e->eleEta)[k], (*e->elePhi)[k], 0.000511);
                TLorentzVector zcand = e1 + e2;

                TLorentzVector e1b; TLorentzVector e2b;
                e1b.SetPtEtaPhiM((*e->elePt)[j], (*e->eleEta)[j], (*e->elePhi)[j], 0.000511);
                e2b.SetPtEtaPhiM((*e->elePt)[k], (*e->eleEta)[k], (*e->elePhi)[k], 0.000511);
                TLorentzVector zcandb = e1b + e2b;

                int det = is_1_endcap + is_2_endcap;
                int index = std::abs((*e->eleCharge)[j] + (*e->eleCharge)[k]) / 2;

                for (int c = 0; c < ncent; ++c) {
                    if (hiBin >= clow[c] && hiBin < chigh[c]) {
                        hminv[det][index][c][0]->Fill(zcand.M());
                        hminv[det][index][c][1]->Fill(zcandb.M());
                    }
                }
            }
        }
    }

    TCanvas* c1 = new TCanvas("c1", "", 400, 400);

#define DRAW_MINV(d, c)                                         \
    {                                                           \
        hminv[d][0][c][0]->SetStats(0);                         \
        hminv[d][0][c][1]->SetStats(0);                         \
                                                                \
        hminv[d][0][c][0]->SetAxisRange(                        \
            0, hminv[d][0][c][0]->GetMaximum() * 1.5, "Y");     \
        hminv[d][0][c][0]->Draw("pe");                          \
        hminv[d][0][c][1]->Draw("pe same");                     \
                                                                \
        TLegend* l1 = new TLegend(0.6, 0.72, 0.8, 0.84);        \
        lstyle(l1, 43, 12);                                     \
        l1->AddEntry(hminv[d][0][c][1], "before", "pl");        \
        l1->AddEntry(hminv[d][0][c][0], "after", "pl");         \
        l1->Draw();                                             \
                                                                \
        c1->SaveAs(Form("minv_%i_c%i_%s.pdf", d, c, modtag));   \
    }

    for (int d = 0; d < 3; ++d) {
        for (int c = 0; c < ncent; ++c) {
            DRAW_MINV(d, c)
        }
    }

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return compare(argv[1], argv[2]);

    if (argc == 4)
        return apply(argv[1], argv[2], argv[3]);

    if (argc == 5)
        return regress(argv[1], argv[2], argv[3], atoi(argv[4]));

    if (argc == 6)
        return validate(argv[1], argv[2], argv[3], atoi(argv[4]), argv[5]);

    printf("usage: %s [tag] [input] [output] [options] [method]\n", argv[0]);
    printf("  options: [0] barrel | [1] endcap\n");

    return 1;
}
