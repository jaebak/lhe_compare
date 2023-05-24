#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <tuple>
#include <regex>
#include <getopt.h>

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TH1D.h"
#include "THStack.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TGraph.h"
#include "TLine.h"

#include "HepMC3.1/LHEF.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::map;
using std::tuple;
using std::pair;
using std::min;
using std::max;
using std::make_tuple;
using std::get;
using std::to_string;

namespace{
  string lhe_a = "";
  string lhe_b = "";
  string lhe_a_tag = "a";
  string lhe_b_tag = "b";
  string output = "";
  long nEvents = -1;
  //long nEvents = 100;
}
void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"lhe_a", required_argument, 0, 'a'},
      {"lhe_b", required_argument, 0, 'b'},
      {"tag_a", required_argument, 0, 0},
      {"tag_b", required_argument, 0, 0},
      {"output", required_argument, 0, 'o'},
      {"nevents", optional_argument, 0, 'n'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "a:b:o:n:", long_options, &option_index);

    if( opt == -1) break;

    string optname;
    switch(opt){
    case 'a':
      lhe_a = optarg;
      break;
    case 'b':
      lhe_b = optarg;
      break;
    case 'o':
      output = optarg;
      break;
    case 'n':
      nEvents = atoi(optarg);
      break;
    case 0:
      optname = long_options[option_index].name;
      if (optname == "tag_a") {
        lhe_a_tag = optarg;
      } else if (optname == "tag_b") {
        lhe_b_tag = optarg;
      } else {
        printf("Bad option! Found option name %s\n", optname.c_str());
      }
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}

TCanvas * new_canvas(string const & name = "", int size = 500)
{
  TSeqCollection * canvases = gROOT->GetListOfCanvases();
  double iCanvas = canvases->GetEntries();
  string canvasName;
  if (name == "") canvasName = "c_g_" + to_string(iCanvas++);
  else canvasName = name;
  return new TCanvas(canvasName.c_str(), canvasName.c_str(), size, size);
}

void setMaximum(float max)
{
  TList * list = gPad->GetListOfPrimitives();
  TIter next(list);
  while (TObject * obj = next())
  {
    std::string className = obj->ClassName();
    if (className.find("TH1") != std::string::npos)
    {
      TH1 * th1 = static_cast<TH1*>(obj);
      th1->SetMaximum(max);
    }
    if (className.find("THStack") != std::string::npos)
    {
      THStack * thstack = static_cast<THStack*>(obj);
      thstack->SetMaximum(max);
    }
  }
  gPad->Modified();
  gPad->Update();

  //TH1 * obj = (TH1*)(gPad->GetListOfPrimitives()->First());
  //obj->SetMaximum(max);
  //gPad->Update();
}

double getMaximumTH1()
{
  TList * list = gPad->GetListOfPrimitives();
  TIter next(list);
  int index = 0;
  double max = 0;
  while (TObject * obj = next())
  {
    std::string className = obj->ClassName();
    if (className.find("TH1") != std::string::npos)
    {
      TH1 * th1 = static_cast<TH1*>(obj);
      double t_max = th1->GetMaximum();
      if (t_max>max || index==0) max = t_max;
    }
    index++;
  }
  return max;
}

double setMaximumTH1(double maxFraction = 1.05)
{
  double max = getMaximumTH1() * maxFraction;
  setMaximum(max);
  return 1;
}

void scale_histogram(TH1F* hist) {
  hist->Scale(1./hist->GetXaxis()->GetBinWidth(1)/hist->Integral());
}

string get_histname(long pid, string variable_name) {
  return "pid_"+to_string(pid)+"_"+variable_name;
}

class ParticleInfo {
  public:
  map<string, double> m_double;
  map<string, vector<int> > m_vint;
  map<string, vector<long> > m_vlong;
  map<string, vector<double> > m_vdouble;

  void set_vlong_branches(vector<string> const & var_names, TTree * tree) {
    for (auto var_name : var_names) {
      m_vlong[var_name];
      tree->Branch(var_name.c_str(), &m_vlong[var_name]);
    }
  }
  void set_vint_branches(vector<string> const & var_names, TTree * tree) {
    for (auto var_name : var_names) {
      m_vint[var_name];
      tree->Branch(var_name.c_str(), &m_vint[var_name]);
    }
  }
  void set_vdouble_branches(vector<string> const & var_names, TTree * tree) {
    for (auto var_name : var_names) {
      m_vdouble[var_name];
      tree->Branch(var_name.c_str(), &m_vdouble[var_name]);
    }
  }
  void set_double_branches(vector<string> const & var_names, TTree * tree) {
    for (auto var_name : var_names) {
      m_double[var_name];
      tree->Branch(var_name.c_str(), &m_double[var_name]);
    }
  }

  void resize (const unsigned nParticles) {
    for (auto it : m_vint) m_vint[it.first].resize(nParticles);
    for (auto it : m_vlong) m_vlong[it.first].resize(nParticles);
    for (auto it : m_vdouble) m_vdouble[it.first].resize(nParticles);
  }

};

void fill_hist_definitions(ParticleInfo const & particle_info, unsigned iParticle, map<string, tuple<string, long, double, double> > & hist_definitions) {

  vector<string> hist_names;
  for (auto var : particle_info.m_vlong) hist_names.push_back(get_histname(particle_info.m_vlong.at("pid").at(iParticle), var.first));
  for (auto var : particle_info.m_vint) hist_names.push_back(get_histname(particle_info.m_vlong.at("pid").at(iParticle), var.first));
  for (auto var : particle_info.m_vdouble) hist_names.push_back(get_histname(particle_info.m_vlong.at("pid").at(iParticle), var.first));

  // Create histogram definitions
  if (hist_definitions.find(hist_names[0]) == hist_definitions.end()) {
    for (auto var : particle_info.m_vlong) {
      long pid = particle_info.m_vlong.at("pid").at(iParticle);
      string histname_var = get_histname(pid, var.first);
      hist_definitions[histname_var] = {var.first, pid, var.second.at(iParticle), var.second.at(iParticle)};
    }
    for (auto var : particle_info.m_vint) {
      long pid = particle_info.m_vlong.at("pid").at(iParticle);
      string histname_var = get_histname(pid, var.first);
      hist_definitions[histname_var] = {var.first, pid, var.second.at(iParticle), var.second.at(iParticle)};
    }
    for (auto var : particle_info.m_vdouble) {
      long pid = particle_info.m_vlong.at("pid").at(iParticle);
      string histname_var = get_histname(pid, var.first);
      hist_definitions[histname_var] = {var.first, pid, var.second.at(iParticle), var.second.at(iParticle)};
    }
  } else {
    // Set min and max
    for (auto var : particle_info.m_vlong) {
      long pid = particle_info.m_vlong.at("pid").at(iParticle);
      string histname_var = get_histname(pid, var.first);
      if (var.second.at(iParticle) < get<2>(hist_definitions[histname_var])) get<2>(hist_definitions[histname_var]) = var.second.at(iParticle);
      if (var.second.at(iParticle) > get<3>(hist_definitions[histname_var])) get<3>(hist_definitions[histname_var]) = var.second.at(iParticle);
    }
    for (auto var : particle_info.m_vlong) {
      long pid = particle_info.m_vlong.at("pid").at(iParticle);
      string histname_var = get_histname(pid, var.first);
      if (var.second.at(iParticle) < get<2>(hist_definitions[histname_var])) get<2>(hist_definitions[histname_var]) = var.second.at(iParticle);
      if (var.second.at(iParticle) > get<3>(hist_definitions[histname_var])) get<3>(hist_definitions[histname_var]) = var.second.at(iParticle);
    }
    for (auto var : particle_info.m_vdouble) {
      long pid = particle_info.m_vlong.at("pid").at(iParticle);
      string histname_var = get_histname(pid, var.first);
      if (var.second.at(iParticle) < get<2>(hist_definitions[histname_var])) get<2>(hist_definitions[histname_var]) = var.second.at(iParticle);
      if (var.second.at(iParticle) > get<3>(hist_definitions[histname_var])) get<3>(hist_definitions[histname_var]) = var.second.at(iParticle);
    }
  }
}

void fill_tree(LHEF::Reader & reader, TTree * tree, ParticleInfo & particle_info, map<string, tuple<string, long, double, double> > & hist_definitions, bool verbose = false) {
  long iEvent = 0;

  particle_info.set_vlong_branches({"pid"}, tree);
  particle_info.set_vint_branches({"status_code", "mother_idx_first", "mother_idx_second", "colors_first", "colors_second"}, tree);
  particle_info.set_vdouble_branches({"px", "py", "pz", "pt", "mass", "rec_mass", "energy", "eta", "phi", "lifetime", "cos_spin_decay"}, tree);
  particle_info.set_double_branches({"weight", "has_z"}, tree);

  while (reader.readEvent()) {

    if (iEvent % 100 == 0) cout<<"Processing entry: "<<iEvent<<endl;

    if (verbose) cout<<"event: "<<iEvent<<endl;
    if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;
    if (verbose) cout<<"comments: "<<reader.eventComments<<endl;
    if (verbose) cout<<"hepeup nup: "<<reader.hepeup.NUP<<endl;

    TLorentzVector particle;

    // Reset storage
    unsigned nParticles = reader.hepeup.IDUP.size();
    particle_info.resize(nParticles);

    // Set event variables
    particle_info.m_double["weight"] = reader.hepeup.XWGTUP;

    // Loop over all particles to set particle variables
    for (unsigned iParticle = 0; iParticle < nParticles; ++iParticle) {
      particle_info.m_vlong["pid"][iParticle] = reader.hepeup.IDUP[iParticle];
      particle_info.m_vint["status_code"][iParticle] = reader.hepeup.ISTUP[iParticle];
      particle_info.m_vdouble["px"][iParticle] = reader.hepeup.PUP[iParticle][0];
      particle_info.m_vdouble["py"][iParticle] = reader.hepeup.PUP[iParticle][1];
      particle_info.m_vdouble["pz"][iParticle] = reader.hepeup.PUP[iParticle][2];
      particle_info.m_vdouble["energy"][iParticle] = reader.hepeup.PUP[iParticle][3];
      particle_info.m_vdouble["mass"][iParticle] = reader.hepeup.PUP[iParticle][4];
      particle_info.m_vdouble["lifetime"][iParticle] = reader.hepeup.VTIMUP[iParticle];
      particle_info.m_vdouble["cos_spin_decay"][iParticle] = reader.hepeup.SPINUP[iParticle];

      pair<int, int> mother_idx = reader.hepeup.MOTHUP[iParticle];
      // Synchronize index with iParticle
      mother_idx.first -= 1;
      mother_idx.second -= 1;
      particle_info.m_vint["mother_idx_first"][iParticle] = mother_idx.first;
      particle_info.m_vint["mother_idx_second"][iParticle] = mother_idx.second;

      pair<int, int> colors = reader.hepeup.ICOLUP[iParticle];
      particle_info.m_vint["colors_first"][iParticle] = colors.first;
      particle_info.m_vint["colors_second"][iParticle] = colors.second;

      //particle.SetPxPyPzE(particle_info.m_vdouble["px"][iParticle], particle_info.m_vdouble["py"][iParticle], particle_info.m_vdouble["pz"][iParticle], particle_info.m_vdouble["energy"][iParticle]);
      //particle_info.m_vdouble["pt"][iParticle] = particle.Pt();
      //particle_info.m_vdouble["eta"][iParticle] = particle.Eta();
      //particle_info.m_vdouble["phi"][iParticle] = particle.Phi();
      //particle_info.m_vdouble["rec_mass"][iParticle] = particle.M();

      if (verbose) cout<<"iParticle: "<<iParticle<<" PID: "<<particle_info.m_vlong["pid"][iParticle]<<" status code: "<<particle_info.m_vint["status_code"][iParticle]<<" mother1: "<<particle_info.m_vint["mother_idx_first"][iParticle]<<" mother2: "<<particle_info.m_vint["mother_idx_second"][iParticle]<<" Mass: "<<particle_info.m_vdouble["mass"][iParticle]<<" px, py, pz, e: "<<particle_info.m_vdouble["px"][iParticle]<<" "<<particle_info.m_vdouble["py"][iParticle]<<" "<<particle_info.m_vdouble["pz"][iParticle]<<" "<<particle_info.m_vdouble["energy"][iParticle]<<" color 1, 2: "<<particle_info.m_vint["colors_first"][iParticle]<<" "<<particle_info.m_vint["colors.second"][iParticle]<<" lifetime: "<<particle_info.m_vdouble["lifetime"][iParticle]<<" cos_spin_decay: "<<particle_info.m_vdouble["cos_spin_decay"][iParticle]<<endl;
      fill_hist_definitions(particle_info, iParticle, hist_definitions);
    }

    // Basic reconstruction

    // Reconstruct OSSF (Opposite sign same flavor) leptons
    // Find OSSF pair
    vector<long> lepton_pids = {11, 13};
    //cout<<"nParticles: "<<nParticles<<endl;
    vector<pair<unsigned, unsigned> > ossf_indices;
    for (unsigned iParticle = 0; iParticle < nParticles; ++iParticle) {
      // Find lepton
      long iParticle_pid = particle_info.m_vlong["pid"][iParticle];
      if (std::find(lepton_pids.begin(), lepton_pids.end(), abs(iParticle_pid)) == lepton_pids.end()) continue;
      for (unsigned jParticle = iParticle+1; jParticle < nParticles; ++jParticle) {
        // Find same sign, opposite flavor lepton
        long jParticle_pid = particle_info.m_vlong["pid"][jParticle];
        if (iParticle_pid != -jParticle_pid) continue;
        ossf_indices.push_back({iParticle, jParticle});

        // Reconstruct OSSF pair
        TLorentzVector lepton_a, lepton_b;
        lepton_a.SetPxPyPzE(particle_info.m_vdouble["px"][iParticle], particle_info.m_vdouble["py"][iParticle], particle_info.m_vdouble["pz"][iParticle], particle_info.m_vdouble["energy"][iParticle]);
        lepton_b.SetPxPyPzE(particle_info.m_vdouble["px"][jParticle], particle_info.m_vdouble["py"][jParticle], particle_info.m_vdouble["pz"][jParticle], particle_info.m_vdouble["energy"][jParticle]);
        TLorentzVector ossf_particle = lepton_a + lepton_b;
        // Save information
        particle_info.m_vlong["pid"].push_back(1000023);
        particle_info.m_vint["status_code"].push_back(0);
        particle_info.m_vdouble["px"].push_back(ossf_particle.Px());
        particle_info.m_vdouble["py"].push_back(ossf_particle.Py());
        particle_info.m_vdouble["pz"].push_back(ossf_particle.Pz());
        particle_info.m_vdouble["energy"].push_back(ossf_particle.E());
        particle_info.m_vdouble["mass"].push_back(ossf_particle.M());
        particle_info.m_vdouble["lifetime"].push_back(0);
        particle_info.m_vdouble["cos_spin_decay"].push_back(0);
        particle_info.m_vint["mother_idx_first"].push_back(0);
        particle_info.m_vint["mother_idx_second"].push_back(0);
        particle_info.m_vint["colors_first"].push_back(0);
        particle_info.m_vint["colors_second"].push_back(0);
        // Set hist definition
        vector<pair<string, string> > ossf_hists {{"ossf_px", "px"}, {"ossf_py", "py"}, {"ossf_pz", "pz"}, {"ossf_energy", "energy"}, {"ossf_mass", "mass"}};
        if (hist_definitions.find(ossf_hists[0].first) == hist_definitions.end()) {
          // Create histogram definitions
          for (auto ossf_hist : ossf_hists) {
            string & hist_name = ossf_hist.first;
            string & var_name = ossf_hist.second;
            double var_value = particle_info.m_vdouble[var_name].back();
            //hist_definitions[hist_name] = {var_name, 1000023, var_value, var_value};
            hist_definitions[hist_name] = {var_name, 1000023, 50, 130};
          }
        } 
        //else {
        //  // Set min and max
        //  for (auto ossf_hist : ossf_hists) {
        //    string & hist_name = ossf_hist.first;
        //    string & var_name = ossf_hist.second;
        //    double var_value = particle_info.m_vdouble[var_name].back();
        //    if (var_value < get<2>(hist_definitions[hist_name])) get<2>(hist_definitions[hist_name]) = var_value;
        //    if (var_value > get<3>(hist_definitions[hist_name])) get<3>(hist_definitions[hist_name]) = var_value;
        //  }
        //}

      }
    }
    //// Select best OSSF
    //double target_mass = 91.188;
    //pair<unsigned, unsigned> best_pair;
    //long best_mass = 0;
    //for (auto pair_indices : ossf_indices) {
    //  unsigned iParticle = pair_indices.first;
    //  unsigned jParticle = pair_indices.second;
    //  TLorentzVector lepton_a, lepton_b;
    //  lepton_a.SetPxPyPzE(particle_info.m_vdouble["px"][iParticle], particle_info.m_vdouble["py"][iParticle], particle_info.m_vdouble["pz"][iParticle], particle_info.m_vdouble["energy"][iParticle]);
    //  lepton_b.SetPxPyPzE(particle_info.m_vdouble["px"][jParticle], particle_info.m_vdouble["py"][jParticle], particle_info.m_vdouble["pz"][jParticle], particle_info.m_vdouble["energy"][jParticle]);
    //  TLorentzVector ossf_particle = lepton_a + lepton_b;
    //  if (fabs(ossf_particle.M()-target_mass) < fabs(best_mass-target_mass)) {
    //    if (best_mass != 0) cout<<"prev_best: "<<best_mass<<" curr_mass: "<<ossf_particle.M()<<endl;
    //    best_mass = ossf_particle.M();
    //    best_pair = {iParticle, jParticle};
    //  }
    //}
    //// Save best OSSF
    //if (best_mass != 0) {
    //  unsigned iParticle = best_pair.first;
    //  unsigned jParticle = best_pair.second;
    //  // Reconstruct OSSF pair
    //  TLorentzVector lepton_a, lepton_b;
    //  lepton_a.SetPxPyPzE(particle_info.m_vdouble["px"][iParticle], particle_info.m_vdouble["py"][iParticle], particle_info.m_vdouble["pz"][iParticle], particle_info.m_vdouble["energy"][iParticle]);
    //  lepton_b.SetPxPyPzE(particle_info.m_vdouble["px"][jParticle], particle_info.m_vdouble["py"][jParticle], particle_info.m_vdouble["pz"][jParticle], particle_info.m_vdouble["energy"][jParticle]);
    //  TLorentzVector ossf_particle = lepton_a + lepton_b;
    //  // Save information
    //  particle_info.m_vlong["pid"].push_back(1000023);
    //  particle_info.m_vint["status_code"].push_back(0);
    //  particle_info.m_vdouble["px"].push_back(ossf_particle.Px());
    //  particle_info.m_vdouble["py"].push_back(ossf_particle.Py());
    //  particle_info.m_vdouble["pz"].push_back(ossf_particle.Pz());
    //  particle_info.m_vdouble["energy"].push_back(ossf_particle.E());
    //  particle_info.m_vdouble["mass"].push_back(ossf_particle.M());
    //  particle_info.m_vdouble["lifetime"].push_back(0);
    //  particle_info.m_vdouble["cos_spin_decay"].push_back(0);
    //  particle_info.m_vint["mother_idx_first"].push_back(0);
    //  particle_info.m_vint["mother_idx_second"].push_back(0);
    //  particle_info.m_vint["colors_first"].push_back(0);
    //  particle_info.m_vint["colors_second"].push_back(0);
    //  // Set hist definition
    //  vector<pair<string, string> > ossf_hists {{"ossf_px", "px"}, {"ossf_py", "py"}, {"ossf_pz", "pz"}, {"ossf_energy", "energy"}, {"ossf_mass", "mass"}};
    //  if (hist_definitions.find(ossf_hists[0].first) == hist_definitions.end()) {
    //    // Create histogram definitions
    //    for (auto ossf_hist : ossf_hists) {
    //      string & hist_name = ossf_hist.first;
    //      string & var_name = ossf_hist.second;
    //      double var_value = particle_info.m_vdouble[var_name].back();
    //      hist_definitions[hist_name] = {var_name, 1000023, var_value, var_value};
    //    }
    //  } else {
    //    // Set min and max
    //    for (auto ossf_hist : ossf_hists) {
    //      string & hist_name = ossf_hist.first;
    //      string & var_name = ossf_hist.second;
    //      double var_value = particle_info.m_vdouble[var_name].back();
    //      if (var_value < get<2>(hist_definitions[hist_name])) get<2>(hist_definitions[hist_name]) = var_value;
    //      if (var_value > get<3>(hist_definitions[hist_name])) get<3>(hist_definitions[hist_name]) = var_value;
    //    }
    //  }
    //}

    // Make weight hist definition
    double weight = particle_info.m_double["weight"];
    if (hist_definitions.find("weight") == hist_definitions.end()) {
      hist_definitions["weight"] = {"weight", 0, weight, weight};
    } else {
      if (weight < get<2>(hist_definitions["weight"])) get<2>(hist_definitions["weight"]) = weight;
      if (weight > get<3>(hist_definitions["weight"])) get<3>(hist_definitions["weight"]) = weight;
    }
    //if (particle_info.m_double["weight"]<0) {
    //  cout<<particle_info.m_double["weight"]<<endl;
    //  cout<<"min: "<<get<2>(hist_definitions["weight"])<<" max: "<<get<3>(hist_definitions["weight"])<<endl;
    //}

    tree->Fill();
    iEvent++;
    if (nEvents != -1 && iEvent == nEvents) break;
  }

}

void fill_histograms(TTree * tree, map<string, tuple<string, long, double, double> > const & hist_definitions, map<string,TH1F*> & histograms, int nbins, string const & tag) {
  // Make histograms
  new_canvas();
  vector<string> ignore_variables = {"rec_mass", "colors_first", "colors_second", "mother_idx_first", "mother_idx_second", "rec_mass", "pt", "eta", "phi", "lifetime", "cos_spin_decay", "weight"};
  for (auto hist_definition : hist_definitions) {
    if (find(ignore_variables.begin(), ignore_variables.end(), get<0>(hist_definition.second)) != ignore_variables.end()) continue;
    histograms[hist_definition.first] = new TH1F((hist_definition.first+"_"+tag).c_str(), hist_definition.first.c_str(), nbins, get<2>(hist_definition.second), get<3>(hist_definition.second));
    //cout<<get<0>(hist_definition.second)<<" "<<hist_definition.first<<" pid=="+to_string(get<1>(hist_definition.second))<<endl;
    //tree->Scan((get<0>(hist_definition.second)+":pid:weight").c_str(), ("pid=="+to_string(get<1>(hist_definition.second))).c_str()); 
    tree->Draw((get<0>(hist_definition.second)+">>"+hist_definition.first+"_"+tag).c_str(), ("(pid=="+to_string(get<1>(hist_definition.second))+")*weight/abs(weight)").c_str());
    //cout<<hist_definition.first<<" "<<"pid=="+to_string(get<1>(hist_definition.second))<<endl;
  }
  //cout<<"min: "<<get<2>(hist_definitions.at("weight"))<<" max: "<<get<3>(hist_definitions.at("weight"))<<endl;
  //cout<<"tree min: "<<tree->GetMinimum("weight")<<" tree max: "<<tree->GetMaximum("weight")<<endl;
  histograms["weight"] = new TH1F(("weight_"+tag).c_str(), "weight", nbins, get<2>(hist_definitions.at("weight"))-1e-8, get<3>(hist_definitions.at("weight"))+1e-8);
  tree->Draw((get<0>(hist_definitions.at("weight"))+">>"+"weight_"+tag).c_str());
}

long get_pid(string const & name) {
  std::smatch pid_match;
  bool found = std::regex_search(name, pid_match, std::regex("pid_(-?[0-9]+)"));
  if (found) return stoi(pid_match[1]);
  else return -999;
}

string get_var(string const & name) {
  std::smatch var_match;
  std::regex_search(name, var_match, std::regex("[0-9]_(.*)"));
  return var_match[1];
}

bool compare_keys(string const & item1, string const & item2) {
  long pid_1 = get_pid(item1);
  long pid_2 = get_pid(item2);
  string var_1 = get_var(item1);
  string var_2 = get_var(item2);
  if (pid_1 == pid_2) return var_1 < var_2;
  return pid_1 < pid_2;
}

vector<string> get_sorted_keys(map<string, TH1F*> const & t_map) {
  vector<string> keys;
  keys.reserve(t_map.size());
  for (auto item : t_map) {
    keys.push_back(item.first);
  }
  sort(keys.begin(), keys.end(), compare_keys);
  return keys;
}

int main(int argc, char *argv[]){
  time_t begtime, endtime;
  time(&begtime);

  GetOptions(argc, argv);

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(1111111);
  //gStyle->SetOptStat(0);

  //string path_a = "lhe_files/ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8.lhe";
  //string path_b = "lhe_files/ZGTo2LG_2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8.lhe";
  string path_a = lhe_a;
  string path_b = lhe_b;

  LHEF::Reader * reader_a = new LHEF::Reader(path_a);
  LHEF::Reader * reader_b = 0;
  if (path_b != "") reader_b = new LHEF::Reader(path_b);

  //TFile * file = TFile::Open("compare_lhe.root","RECREATE");
  TFile * file = TFile::Open((output+".root").c_str(),"RECREATE");


  // Fill particle tree and hist definitions
  // hist_definitions[hist_name] = [var, pid, min, max]
  map<string, tuple<string, long, double, double> > hist_definitions;
  TTree * tree_a;
  tree_a = new TTree("tree_a", "tree_a");
  ParticleInfo particle_info_a;
  fill_tree(*reader_a, tree_a, particle_info_a, hist_definitions);
  tree_a->Write();

  TTree * tree_b = 0;
  ParticleInfo particle_info_b;
  if (path_b != "") {
    tree_b = new TTree("tree_b", "tree_b");
    fill_tree(*reader_b, tree_b, particle_info_b, hist_definitions);
    tree_b->Write();
  }

  //// Print hist definitions
  //for (auto hist_definition : hist_definitions) {
  //  cout<<hist_definition.first<<" "<<get<2>(hist_definition.second)<<" "<<get<3>(hist_definition.second)<<endl;
  //}

  // Fill histograms
  map<string,TH1F*> histograms_a;
  map<string,TH1F*> histograms_b;

  int nbins = 30;
  fill_histograms(tree_a, hist_definitions, histograms_a, nbins, /*tag*/ lhe_a_tag);
  if (path_b != "") fill_histograms(tree_b, hist_definitions, histograms_b, nbins, /*tag*/ lhe_b_tag);

  //for (auto & histogram_a : histograms_a) {
  //  cout<<histogram_a.first<<endl;
  //  //TCanvas * canvas = new_canvas();
  //  //histogram_a.second->SetMinimum(0);
  //  //histogram_a.second->Draw();
  //  //canvas->SaveAs(("plots/"+histogram_a.first+".pdf").c_str());
  //}

  vector<string> keys = get_sorted_keys(histograms_a);

  //string paper_name = "plots/compare.pdf";
  string paper_name = output+".pdf";
  TCanvas * paper = new_canvas();
  paper->Print((paper_name+"[").c_str());

  for (string const & key : keys) {
    cout<<key<<endl;
    TCanvas * canvas = new_canvas();
    TPad * up_pad = new TPad("up_pad", "up_pad", /*xlow*/0, /*ylow*/0.3, /*xhigh*/1., /*yhigh*/1);
    TPad * low_pad = new TPad("low_pad", "low_pad", /*xlow*/0, /*ylow*/0, /*xhigh*/1., /*yhigh*/0.3);

    up_pad->Draw();
    low_pad->Draw();

    up_pad->cd();
    histograms_a[key]->SetMinimum(0);
    histograms_a[key]->SetLineWidth(2);
    histograms_a[key]->SetLineColor(9);
    double x_min = histograms_a[key]->GetXaxis()->GetBinLowEdge(histograms_a[key]->GetXaxis()->GetFirst());
    double x_max = histograms_a[key]->GetXaxis()->GetBinUpEdge(histograms_a[key]->GetXaxis()->GetLast());
    char ytitle[64]; snprintf(ytitle, sizeof ytitle, "Events / %.3g", (x_max-x_min)/nbins);
    histograms_a[key]->GetYaxis()->SetTitle(ytitle);
    histograms_a[key]->Draw("HIST E");
    gPad->Update(); // To force draw for statbox
    TPaveStats *statbox_a = static_cast<TPaveStats*>(histograms_a[key]->FindObject("stats"));
    statbox_a->SetTextColor(9);
    statbox_a->SetY1NDC(0.7);
    statbox_a->SetY2NDC(0.9);

    if (path_b != "") {
      histograms_b[key]->SetLineColor(kRed);
      histograms_b[key]->Draw("HIST E sames");
      gPad->Update(); // To force draw for statbox
      TPaveStats *statbox_b = static_cast<TPaveStats*>(histograms_b[key]->FindObject("stats"));
      statbox_b->SetTextColor(kRed);
      statbox_b->SetY1NDC(0.5);
      statbox_b->SetY2NDC(0.7);
    }

    setMaximumTH1(1.3);

    //TLegend * legend = new TLegend(/*x1*/0.1, /*y1*/0.8, /*x2*/0.3, /*y2*/0.9);
    //legend->AddEntry(histograms_a[key], lhe_a_tag.c_str());
    //legend->AddEntry(histograms_b[key], lhe_b_tag.c_str());
    //legend->Draw();

    if (path_b != "") {
      low_pad->cd();
      if (fabs(x_max - x_min)>0.000001) { // Prevent residuals for one bin histograms
        vector<Double_t> res_y(static_cast<unsigned>(nbins));
        vector<Double_t> res_x(static_cast<unsigned>(nbins));
        for (unsigned index = 0; index < static_cast<unsigned>(nbins); ++index) res_x[index] = x_min + (x_max-x_min)/nbins * index + (x_max-x_min)/nbins/2 ;
        gErrorIgnoreLevel = kError;
        double pvalue = histograms_b[key]->Chi2Test(histograms_a[key], "WW", &res_y[0]);
        gErrorIgnoreLevel = kPrint;
        TGraph * residual = new TGraph(nbins, &res_x[0], &res_y[0]);
        residual->SetTitle(0);
        // Set range on computed graph. Set twice because TGraph's axis looks strange otherwise
        residual->GetXaxis()->SetLimits(x_min, x_max);
        residual->GetYaxis()->SetRangeUser(-3.5, 3.5);
        residual->GetYaxis()->SetTitle("Normalized residuals");
        residual->GetYaxis()->SetTitleSize(0.07);
        residual->GetYaxis()->SetTitleOffset(0.5);
        residual->GetYaxis()->CenterTitle();
        residual->SetMarkerStyle(21);
        residual->SetMarkerSize(0.3);
        residual->Draw("AP");
        TLine * line = new TLine(/*x1*/x_min,/*y1*/0,/*x2*/x_max,/*y2*/0);
        line->Draw();
        TPaveText * text_box = new TPaveText(/*x1*/0.12, /*y1*/0.8, /*x2*/0.5, /*y2*/0.9, "NDC NB");
        text_box->SetFillColorAlpha(0,0);
        char c_pvalue[64]; snprintf(c_pvalue, sizeof c_pvalue, "#chi^{2} test (%s vs %s) p-value: %.3g", lhe_a_tag.c_str(), lhe_b_tag.c_str(), pvalue);
        text_box->AddText(c_pvalue);
        text_box->Draw();
      }
    }

    canvas->Print(paper_name.c_str());
  }
  paper->Print((paper_name+"]").c_str());

  file->Close();

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
  return 0;
}

