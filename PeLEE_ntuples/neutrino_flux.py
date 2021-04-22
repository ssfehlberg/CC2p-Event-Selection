import numpy as np
import os
import ROOT 

infile = ROOT.TFile("/pnfs/uboone/persistent/uboonebeam/bnb_gsimple/bnb_gsimple_fluxes_01.09.2019_463_hist/MCC9_FluxHist_volTPCActive.root","READ")

h_original = infile.Get("hEnumu_cv")
h_scaled = h_original.Clone("h_scaled")
h = h_original.Clone("h_scaled_w_fit")

outfile = ROOT.TFile("root_files/Run_all/neutrino_flux.root","RECREATE")
h_original.Write()
scale_factor= 1/(4997.*5e8)/(256.35*233);
h_scaled.Scale(scale_factor)
h_scaled.Write()

sample_name="#scale[0.6]{MicroBooNE In-Progress}";
t = ROOT.TLatex()
ROOT.gStyle.SetPaintTextFormat("4.2f")
ROOT.gStyle.SetHistMinimumZero(ROOT.kTRUE)
#ROOT.gStyle.SetHistMinimumZero(ROOT.kFALSE)
t.SetNDC();
t.SetTextAlign(22);

canv = ROOT.TCanvas("canv","canv",700,500)
scale_factor= 1/(4997.*5e8)/(256.35*233);
h.Scale(scale_factor)
h.Draw("HIST")
h.GetXaxis().SetRangeUser(0.0,4.0)
h.GetXaxis().SetTitle("E_{#nu} (GeV)")
h.GetYaxis().SetTitle("# #nu/4.183e+7 POT/GeV/cm^{2}")
g = ROOT.TF1("g","gaus")
fit = h.Fit(g,"Q")
parameters = g.GetParameters() 
ymax = 33.5E-12#h.GetMaximum()

mean = parameters[1]
mean_line = ROOT.TLine(mean,0,mean,ymax)
mean_line.Draw("SAME")
mean_line.SetLineColor(2)
mean_line.SetLineStyle(1)

sigma = parameters[2]
sigma_plus = mean+sigma
sigma_minus = mean-sigma

a_plus = ROOT.TLine(sigma_plus,0,sigma_plus,ymax)
a_plus.Draw("SAME")
a_plus.SetLineColor(4)
a_plus.SetLineStyle(9)

a_minus = ROOT.TLine(sigma_minus,0,sigma_minus,ymax)
a_minus.Draw("SAME")
a_minus.SetLineColor(4)
a_minus.SetLineStyle(9)

legend = ROOT.TLegend(0.6,0.6,0.9,0.9)
legend.AddEntry(h,"Simulated BNB #nu_{#mu} Flux","L")
legend.AddEntry(mean_line,"<E_{#nu}>= %f GeV"%mean,"L")
legend.AddEntry(a_plus,"+1#sigma Energy Range= %f GeV"%sigma_plus,"L")
legend.AddEntry(a_minus,"-1#sigma Energy Range= %f GeV"%sigma_minus,"L")
legend.Draw("SAME")
legend.SetBorderSize(0)
t.DrawLatex(0.8,0.93,"%s"%sample_name)

canv.Print("root_files/Run_all/neutrino_flux.pdf")
canv.Print("root_files/Run_all/neutrino_flux.png")

integral = h.Integral()
print("Integral of the Histogram: %E"%integral)

#write histogram and close output file
h.Write()
outfile.Close()
