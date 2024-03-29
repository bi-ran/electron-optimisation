#ifndef COSMETICS_H
#define COSMETICS_H

#include "TLatex.h"

void watermark() {
   TLatex* lcms = new TLatex();
   lcms->SetTextFont(62);
   lcms->SetTextSize(0.052);
   lcms->SetTextAlign(13);
   lcms->DrawLatexNDC(0.135, 0.875, "CMS");

   TLatex* lprelim = new TLatex();
   lprelim->SetTextFont(52);
   lprelim->SetTextSize(0.032);
   lprelim->SetTextAlign(13);
   lprelim->DrawLatexNDC(0.135, 0.83, "Performance");

   TLatex* linfo = new TLatex();
   linfo->SetTextFont(42);
   linfo->SetTextSize(0.032);
   linfo->SetTextAlign(31);
   linfo->DrawLatexNDC(0.89, 0.92, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
}

#include "TH1.h"

template<typename T>
void hstyle_title_label_size(T* h, uint32_t title_size, uint32_t label_size) {
   float npixelspad = gPad->GetWh() * gPad->GetAbsHNDC();

   h->GetXaxis()->SetLabelSize(label_size / npixelspad);
   h->GetXaxis()->SetTitleSize(title_size / npixelspad);
   h->GetYaxis()->SetLabelSize(label_size / npixelspad);
   h->GetYaxis()->SetTitleSize(title_size / npixelspad);
   h->GetYaxis()->CenterTitle();
   h->GetYaxis()->SetTitleOffset(0.5);
}

template<typename T>
void hstyle(T* h, int style, int colour, float size) {
   hstyle(h, style, colour);
   h->SetMarkerSize(size);
}

template<typename T>
void hstyle(T* h, int style, int colour) {
   h->SetMarkerStyle(style);
   h->SetMarkerSize(1.2);
   h->SetMarkerColor(colour);
   h->SetLineColor(colour);
}

void htitle(TH1* h, const char* title) {
   h->SetTitle(title);
   h->GetXaxis()->CenterTitle();
   h->GetXaxis()->SetTitleOffset(1.2);
   h->GetYaxis()->CenterTitle();
   h->GetYaxis()->SetTitleOffset(1.2);
}

void haxes(TH1* h, float ymin, float ymax) {
   h->SetAxisRange(ymin, ymax, "Y");
}

void hformat(TH1* h, int style, int colour, const char* title) {
   h->SetStats(0);

   hstyle(h, style, colour);
   htitle(h, title);
}

void hformat(TH1* h, float ymin, float ymax, const char* title) {
   h->SetStats(0);

   htitle(h, title);
   haxes(h, ymin, ymax);
}

void hformat(TH1* h, int style, int colour, float ymin, float ymax,
             const char* title) {
   h->SetStats(0);

   hstyle(h, style, colour);
   htitle(h, title);
   haxes(h, ymin, ymax);
}

void gstyle(TH1* h, int style, int colour) {
   h->SetMarkerSize(0);
   h->SetMarkerColor(colour);
   h->SetLineStyle(style);
   h->SetLineColor(colour);
}

void gformat(TH1* h, int style, int colour) {
   h->SetStats(0);

   gstyle(h, style, colour);
}

#include "TLegend.h"

void lstyle(TLegend* l, int font, float size) {
   l->SetBorderSize(0);
   l->SetFillStyle(0);
   l->SetTextFont(font);
   l->SetTextSize(size);
}

#endif  /* COSMETICS_H */
