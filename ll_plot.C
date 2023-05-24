void ll_plot() {
  TString filename = "run2_run3_ZGamma2JToGamma2L2J_EWK.root";
  TFile * file = new TFile(filename);
  TTree* tree_a = (TTree*)file->Get("tree_a");
  //tree_a->Draw("mass","(pid==1000023&&mass<150&&mass>50)*weight/abs(weight)", "hist");
  tree_a->Draw("mass","(pid==1000023&&mass<150&&mass>50)*weight/abs(weight)");
  ((TH1F*)gPad->GetPrimitive("htemp"))->SetName("htemp1");
  ((TH1F*)gPad->GetPrimitive("htemp1"))->SetLineColor(kRed);
  tree_a->Draw("mass","(pid==1000023&&mass<150&&mass>50&&MinIf$(pid,pid==23)==23)*weight/abs(weight)", "same hist");
  ((TH1F*)gPad->GetPrimitive("htemp"))->SetName("htemp2");
  ((TH1F*)gPad->GetPrimitive("htemp2"))->SetFillColor(kRed);
  ((TH1F*)gPad->GetPrimitive("htemp2"))->SetLineColor(kRed);
  gPad->Modified();

  TTree* tree_b = (TTree*)file->Get("tree_b");
  //tree_b->Draw("mass","(pid==1000023&&mass<150&&mass>50)*weight/abs(weight)", "hist same");
  tree_b->Draw("mass","(pid==1000023&&mass<150&&mass>50)*weight/abs(weight)", "same");
  ((TH1F*)gPad->GetPrimitive("htemp"))->SetName("htemp3");
  tree_b->Draw("mass","(pid==1000023&&mass<150&&mass>50&&MinIf$(pid,pid==23)==23)*weight/abs(weight)", "same hist");
  ((TH1F*)gPad->GetPrimitive("htemp"))->SetName("htemp4");
  ((TH1F*)gPad->GetPrimitive("htemp4"))->SetFillColor(kBlue);
  gPad->Modified();
}
