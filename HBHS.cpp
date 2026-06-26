#include "HBHS.h"
#include "DecayTable.h"
#include "HiggsBounds.h"
#include "SM.h"
#include "THDM.h"
#include <array>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>

#if defined HiggsBounds


namespace {
double zeroIfNaN(double x) { return x == x ? x : 0.; }

constexpr bool debug = false;
} // namespace
namespace HP = Higgs::predictions;

HBHS::HBHS(int hs_pdf) {
  initialize_HiggsBounds(nHzero, nHplus, 3);

}

HBHS::~HBHS() {
  finish_HiggsBounds();
//  finish_HiggsSignals();
}

HBHSResult HBHS::check(THDM &model) const {

  neutral_input(model);
  charged_input(model);

  HBHSResult result{};

  const auto ResultsHB = bounds(pred);


//  cout << "\n\n" << ResultsHB << "\n\n\n";

  result.hb.result.at(0) = ResultsHB.allowed;

  vector<std::string> keys;
  for (auto a = ResultsHB.selectedLimits.begin(); a != ResultsHB.selectedLimits.end(); ++a){
    keys.push_back(a->first);
  }

  int most_sen = 0;
  for (int j=0; j<5; j++){
    result.hb.result.at(j) = -2;
    result.hb.obsratio.at(j) = -2;
    result.hb.channel.at(j) = "-2";
    result.hb.ncombined.at(j) = "-2";
  }
  for (auto key : keys){
    if (key == "h1"){
      result.hb.result.at(1) = (ResultsHB.selectedLimits.at("h1").obsRatio() < 1 ? 1 : 0);
      result.hb.obsratio.at(1) = ResultsHB.selectedLimits.at("h1").obsRatio();
      result.hb.channel.at(1) = ResultsHB.selectedLimits.at("h1").limit().get()->citeKey();
      result.hb.ncombined.at(1) = "";
      for (long unsigned int i=0; i<ResultsHB.selectedLimits.at("h1").contributingParticles().size(); i++){
        result.hb.ncombined.at(1).append(ResultsHB.selectedLimits.at("h1").contributingParticles().at(i));
      }

    }if (key == "h2"){
      result.hb.result.at(2) = (ResultsHB.selectedLimits.at("h2").obsRatio() < 1 ? 1 : 0);
      result.hb.obsratio.at(2) = ResultsHB.selectedLimits.at("h2").obsRatio();
      result.hb.channel.at(2) = ResultsHB.selectedLimits.at("h2").limit().get()->citeKey();
      result.hb.ncombined.at(2) = "";
      for (long unsigned int i=0; i<ResultsHB.selectedLimits.at("h2").contributingParticles().size(); i++){
        result.hb.ncombined.at(2).append(ResultsHB.selectedLimits.at("h2").contributingParticles().at(i));
      }

    }if (key == "a1"){
      result.hb.result.at(3) = (ResultsHB.selectedLimits.at("a1").obsRatio() < 1 ? 1 : 0);
      result.hb.obsratio.at(3) = ResultsHB.selectedLimits.at("a1").obsRatio();
      result.hb.channel.at(3) = ResultsHB.selectedLimits.at("a1").limit().get()->citeKey();
      result.hb.ncombined.at(3) = "";
      for (long unsigned int i=0; i<ResultsHB.selectedLimits.at("a1").contributingParticles().size(); i++){
        result.hb.ncombined.at(3).append(ResultsHB.selectedLimits.at("a1").contributingParticles().at(i));
      }
    }if (key == "hp"){
      result.hb.result.at(4) = (ResultsHB.selectedLimits.at("hp").obsRatio() < 1 ? 1 : 0);
      result.hb.obsratio.at(4) = ResultsHB.selectedLimits.at("hp").obsRatio();
      result.hb.channel.at(4) = ResultsHB.selectedLimits.at("hp").limit().get()->citeKey();
      result.hb.ncombined.at(4) = "";
      for (long unsigned int i=0; i<ResultsHB.selectedLimits.at("hp").contributingParticles().size(); i++){
        result.hb.ncombined.at(4).append(ResultsHB.selectedLimits.at("hp").contributingParticles().at(i));
      }
    }
  }
  for (int j=1; j<5; j++){
    if (result.hb.obsratio.at(0) < result.hb.obsratio.at(j)){
      result.hb.obsratio.at(0) = result.hb.obsratio.at(j);
      most_sen = j;
    }
  }
  result.hb.result.at(0) = result.hb.result.at(most_sen);
  result.hb.channel.at(0) = result.hb.channel.at(most_sen);
  result.hb.ncombined.at(0) = result.hb.ncombined.at(most_sen);



/*  if (debug) {
    for (size_t h = 1; h <= nHzero; ++h) {
      double BR_hjss, BR_hjcc, BR_hjbb, BR_hjtt, BR_hjmumu, BR_hjtautau,
          BR_hjWW, BR_hjZZ, BR_hjZga, BR_hjgaga, BR_hjgg;
      HiggsBounds_get_neutral_BR(h, &BR_hjss, &BR_hjcc, &BR_hjbb, &BR_hjtt,
                                 &BR_hjmumu, &BR_hjtautau, &BR_hjWW, &BR_hjZZ,
                                 &BR_hjZga, &BR_hjgaga, &BR_hjgg);
      printf(
          "HiggsBounds BRs:\n %2ld %16.8E %16.8E %16.8E %16.8E %16.8E %16.8E "
          "%16.8E %16.8E %16.8E %16.8E %16.8E\n",
          h, BR_hjss, BR_hjcc, BR_hjbb, BR_hjtt, BR_hjmumu, BR_hjtautau,
          BR_hjWW, BR_hjZZ, BR_hjZga, BR_hjgaga, BR_hjgg);
      double singleH, ggH, bbH, VBF, WH, ZH, ttH, tH_tchan, tH_schan, qqZH,
          ggZH;
      HiggsBounds_get_neutral_hadr_CS(h, LHC13, &singleH, &ggH, &bbH, &VBF, &WH,
                                      &ZH, &ttH, &tH_tchan, &tH_schan, &qqZH,
                                      &ggZH);
      printf(
          "HiggsBounds LHC13 cxns:\n %2ld %16.8E %16.8E %16.8E %16.8E %16.8E "
          "%16.8E %16.8E %16.8E %16.8E %16.8E %16.8E\n",
          h, singleH, ggH, bbH, VBF, WH, ZH, ttH, tH_tchan, tH_schan, qqZH,
          ggZH);
    }
  }*/


  result.hs.chisq = signals(pred);
  result.hs.nobs = signals.observableCount();

  return result;
}

HBHS::EffC HBHS::effective_couplings(THDM &model) const {
  EffC effC{};
  std::array<double, 3> Mh = {model.get_hmass(1), model.get_hmass(2),
                              model.get_hmass(3)};

  DecayTable table(model);
  THDM sm_like;
  for (size_t h = 1; h <= nHzero; ++h) {
    sm_like.set_param_sm(Mh[h - 1]);
    sm_like.set_yukawas_type(1);

    std::complex<double> cs, cp, cs_sm, cp_sm;
    // h -> ss
    model.get_coupling_hdd(h, 2, 2, cs, cp);
    sm_like.get_coupling_hdd(1, 2, 2, cs_sm, cp_sm);
    effC.ghjss_s[h - 1] = cs.imag() / cs_sm.imag();
    effC.ghjss_p[h - 1] = -cp.real() / cs_sm.imag();
    if (debug)
      printf("%2ld %5s %16.8E %16.8E\n", h, "ss", effC.ghjss_s[h - 1],
             effC.ghjss_p[h - 1]);

    model.get_coupling_hdd(h, 3, 3, cs, cp);
    sm_like.get_coupling_hdd(1, 3, 3, cs_sm, cp_sm);
    effC.ghjbb_s[h - 1] = cs.imag() / cs_sm.imag();
    effC.ghjbb_p[h - 1] = -cp.real() / cs_sm.imag();
    if (debug)
      printf("%2ld %5s %16.8E %16.8E\n", h, "bb", effC.ghjbb_s[h - 1],
             effC.ghjbb_p[h - 1]);

    model.get_coupling_huu(h, 2, 2, cs, cp);
    sm_like.get_coupling_huu(1, 2, 2, cs_sm, cp_sm);
    effC.ghjcc_s[h - 1] = cs.imag() / cs_sm.imag();
    effC.ghjcc_p[h - 1] = -cp.real() / cs_sm.imag();
    if (debug)
      printf("%2ld %5s %16.8E %16.8E\n", h, "cc", effC.ghjcc_s[h - 1],
             effC.ghjcc_p[h - 1]);

    model.get_coupling_huu(h, 3, 3, cs, cp);
    sm_like.get_coupling_huu(1, 3, 3, cs_sm, cp_sm);
    effC.ghjtt_s[h - 1] = cs.imag() / cs_sm.imag();
    effC.ghjtt_p[h - 1] = -cp.real() / cs_sm.imag();
    if (debug)
      printf("%2ld %5s %16.8E %16.8E\n", h, "tt", effC.ghjtt_s[h - 1],
             effC.ghjtt_p[h - 1]);

    model.get_coupling_hll(h, 2, 2, cs, cp);
    sm_like.get_coupling_hll(1, 2, 2, cs_sm, cp_sm);
    effC.ghjmumu_s[h - 1] = cs.imag() / cs_sm.imag();
    effC.ghjmumu_p[h - 1] = -cp.real() / cs_sm.imag();
    if (debug)
      printf("%2ld %5s %16.8E %16.8E\n", h, "mumu", effC.ghjmumu_s[h - 1],
             effC.ghjmumu_p[h - 1]);

    model.get_coupling_hll(h, 3, 3, cs, cp);
    sm_like.get_coupling_hll(1, 3, 3, cs_sm, cp_sm);
    effC.ghjtautau_s[h - 1] = cs.imag() / cs_sm.imag();
    effC.ghjtautau_p[h - 1] = -cp.real() / cs_sm.imag();
    if (debug)
      printf("%2ld %5s %16.8E %16.8E\n", h, "tata", effC.ghjtautau_s[h - 1],
             effC.ghjtautau_p[h - 1]);

    model.get_coupling_vvh(2, 2, h, cs);
    sm_like.get_coupling_vvh(2, 2, 1, cs_sm);
    effC.ghjZZ[h - 1] = cs.imag() / cs_sm.imag();
    if (debug)
      printf("%2ld %5s %16.8E\n", h, "ZZ", effC.ghjZZ[h - 1]);

    model.get_coupling_vvh(3, 3, h, cs);
    sm_like.get_coupling_vvh(3, 3, 1, cs_sm);
    effC.ghjWW[h - 1] = cs.imag() / cs_sm.imag();
    if (debug)
      printf("%2ld %5s %16.8E\n", h, "WW", effC.ghjWW[h - 1]);

    DecayTable sm_table(sm_like);
    double hgaga = table.get_gamma_hgaga(h);
    double hgaga_sm = sm_table.get_gamma_hgaga(1);
    effC.ghjgaga[h - 1] = sqrt(hgaga / hgaga_sm);
    if (debug)
      printf("%2ld %5s %16.8E\n", h, "gaga", effC.ghjgaga[h - 1]);

    double hZga = table.get_gamma_hZga(h);
    double hZga_sm = sm_table.get_gamma_hZga(1);
    effC.ghjZga[h - 1] = zeroIfNaN(sqrt(hZga / hZga_sm));
    if (debug)
      printf("%2ld %5s %16.8E\n", h, "Zga", effC.ghjZga[h - 1]);

    double hgg = table.get_gamma_hgg(h);
    double hgg_sm = sm_table.get_gamma_hgg(1);
    effC.ghjgg[h - 1] = sqrt(hgg / hgg_sm);
    if (debug)
      printf("%2ld %5s %16.8E\n", h, "gg", effC.ghjgg[h - 1]);

    for (size_t j = 1; j <= nHzero; ++j) {
      std::complex<double> c;
      model.get_coupling_vhh(2, h, j, c);
      effC.ghjhiZ[h - 1][j - 1] =
          c.real() / (model.get_SM().get_g() / 2. / model.get_SM().get_costw());
    }
    if (debug)
      printf("%2ld %5s %16.8E %16.8E %16.8E\n", h, "Zhj", effC.ghjhiZ[h - 1][0],
             effC.ghjhiZ[h - 1][1], effC.ghjhiZ[h - 1][2]);
  }
  return effC;
}

HBHS::NonSMBR HBHS::nonSM_branching_ratios(THDM &model) {
  DecayTable table(model);
  NonSMBR nonSMBR{};
  for (size_t h = 1; h <= nHzero; ++h) {
    nonSMBR.BR_hjinvisible[h - 1] = 0.0;
    nonSMBR.BR_hjHpW[h - 1][0] =
        zeroIfNaN(table.get_gamma_hvh(h, 3, 4) / table.get_gammatot_h(h));
    if (debug)
      printf("BR h%2ld -> H+- W-+ = %16.8E\n", h, nonSMBR.BR_hjHpW[h - 1][0]);
    nonSMBR.BR_hjemu[h - 1] =
        zeroIfNaN(table.get_gamma_hll(h, 1, 2) / table.get_gammatot_h(h));
    nonSMBR.BR_hjetau[h - 1] =
        zeroIfNaN(table.get_gamma_hll(h, 1, 3) / table.get_gammatot_h(h));
    nonSMBR.BR_hjmutau[h - 1] =
        zeroIfNaN(table.get_gamma_hll(h, 2, 3) / table.get_gammatot_h(h));
  }
  for (size_t j = 1; j <= nHzero; j++) {
    for (size_t i = 1; i <= nHzero; i++) {
      for (size_t k = 1; k <= nHzero; k++) {
        nonSMBR.BR_hkhjhi[k - 1][j - 1][i - 1] =
            zeroIfNaN(table.get_gamma_hhh(k, j, i) / table.get_gammatot_h(k));
        if (debug)
          printf("BR h%2ld -> h%2ld h%2ld = %16.8E\n", k, j, i,
                 nonSMBR.BR_hkhjhi[k - 1][j - 1][i - 1]);
      }
      nonSMBR.BR_hjhiZ[j - 1][i - 1] =
          zeroIfNaN(table.get_gamma_hvh(j, 2, i) / table.get_gammatot_h(j));
      if (debug)
        printf("BR h%2ld -> h%2ld Z = %16.8E\n", j, i,
               nonSMBR.BR_hjhiZ[j - 1][i - 1]);
    }
  }
  if (model.get_yukawas_type() == 0) {
    if (model.get_hmass(1) < model.get_hmass(3)) {
      nonSMBR.BR_hjinvisible[1] = nonSMBR.BR_hkhjhi[1][0][0];
      nonSMBR.BR_hkhjhi[1][0][0] = 0;
    } else if (model.get_hmass(1) > model.get_hmass(3)) {
      nonSMBR.BR_hjinvisible[1] = nonSMBR.BR_hkhjhi[1][2][2];
      nonSMBR.BR_hkhjhi[1][2][2] = 0;
    } else {
      nonSMBR.BR_hjinvisible[1] =
          nonSMBR.BR_hkhjhi[1][0][0] + nonSMBR.BR_hkhjhi[1][2][2];
      nonSMBR.BR_hkhjhi[1][0][0] = 0;
      nonSMBR.BR_hkhjhi[1][2][2] = 0;
    }
    if (debug)
      printf("BR_hSM->invisible = %16.8E\n", nonSMBR.BR_hjinvisible[1]);
  }
  return nonSMBR;
}

//void HBHS::neutral_input(THDM &model, Higgs::predictions::Predictions pred) const {
void HBHS::neutral_input(THDM &model) const {
  // properties



  h1.setMass(model.get_hmass(1));
  h2.setMass(model.get_hmass(2));
  a1.setMass(model.get_hmass(3));



  DecayTable table{model};
  vector<double> TotalWidths = {table.get_gammatot_h(1), table.get_gammatot_h(2), table.get_gammatot_h(3)};

  if (model.get_yukawas_type() == 0) {
    h1.setTotalWidth(TotalWidths.at(0));
    a1.setTotalWidth(TotalWidths.at(2));
  }

  // effective couplings
  const EffC effC = effective_couplings(model);

  auto clps_h1 = HP::NeutralEffectiveCouplings();
  auto clps_h2 = HP::NeutralEffectiveCouplings();
  auto clps_a1 = HP::NeutralEffectiveCouplings();


  clps_h1.ss     = effC.ghjss_s.at(0)     + effC.ghjss_p.at(0)*1i;
  clps_h1.cc     = effC.ghjcc_s.at(0)     + effC.ghjcc_p.at(0)*1i;
  clps_h1.bb     = effC.ghjbb_s.at(0)     + effC.ghjbb_p.at(0)*1i;
  clps_h1.tt     = effC.ghjtt_s.at(0)     + effC.ghjtt_p.at(0)*1i;
  clps_h1.mumu   = effC.ghjmumu_s.at(0)   + effC.ghjmumu_p.at(0)*1i;
  clps_h1.tautau = effC.ghjtautau_s.at(0) + effC.ghjtautau_p.at(0)*1i;
  clps_h1.WW     = effC.ghjWW.at(0);
  clps_h1.ZZ     = effC.ghjZZ.at(0);
  clps_h1.Zgam   = effC.ghjZga.at(0);
  clps_h1.gamgam = effC.ghjgaga.at(0);
  clps_h1.gg     = effC.ghjgg.at(0);

  clps_h2.ss     = effC.ghjss_s.at(1)     + effC.ghjss_p.at(1)*1i;
  clps_h2.cc     = effC.ghjcc_s.at(1)     + effC.ghjcc_p.at(1)*1i;
  clps_h2.bb     = effC.ghjbb_s.at(1)     + effC.ghjbb_p.at(1)*1i;
  clps_h2.tt     = effC.ghjtt_s.at(1)     + effC.ghjtt_p.at(1)*1i;
  clps_h2.mumu   = effC.ghjmumu_s.at(1)   + effC.ghjmumu_p.at(1)*1i;
  clps_h2.tautau = effC.ghjtautau_s.at(1) + effC.ghjtautau_p.at(1)*1i;
  clps_h2.WW     = effC.ghjWW.at(1);
  clps_h2.ZZ     = effC.ghjZZ.at(1);
  clps_h2.Zgam   = effC.ghjZga.at(1);
  clps_h2.gamgam = effC.ghjgaga.at(1);
  clps_h2.gg     = effC.ghjgg.at(1);

  clps_a1.ss     = effC.ghjss_s.at(2)     + effC.ghjss_p.at(2)*1i;
  clps_a1.cc     = effC.ghjcc_s.at(2)     + effC.ghjcc_p.at(2)*1i;
  clps_a1.bb     = effC.ghjbb_s.at(2)     + effC.ghjbb_p.at(2)*1i;
  clps_a1.tt     = effC.ghjtt_s.at(2)     + effC.ghjtt_p.at(2)*1i;
  clps_a1.mumu   = effC.ghjmumu_s.at(2)   + effC.ghjmumu_p.at(2)*1i;
  clps_a1.tautau = effC.ghjtautau_s.at(2) + effC.ghjtautau_p.at(2)*1i;
  clps_a1.WW     = effC.ghjWW.at(2);
  clps_a1.ZZ     = effC.ghjZZ.at(2);
  clps_a1.Zgam   = effC.ghjZga.at(2);
  clps_a1.gamgam = effC.ghjgaga.at(2);
  clps_a1.gg     = effC.ghjgg.at(2);




  HP::effectiveCouplingInput(h1, clps_h1, HP::ReferenceModel::SMHiggs);
  HP::effectiveCouplingInput(h2, clps_h2, HP::ReferenceModel::SMHiggs);
  HP::effectiveCouplingInput(a1, clps_a1, HP::ReferenceModel::SMHiggs);


  // non-SM decay branching ratios
  const auto nonSMBR = nonSM_branching_ratios(model);

  h1.setDecayWidth(HP::Decay::directInv,    nonSMBR.BR_hjinvisible.at(0)*TotalWidths.at(0));
  h1.setDecayWidth("h1","h1",               nonSMBR.BR_hkhjhi.at(0)[0][0]*TotalWidths.at(0));
  h1.setDecayWidth(HP::ChainDecay::Z, "h1", nonSMBR.BR_hjhiZ.at(0)[0]*TotalWidths.at(0));
  h1.setDecayWidth(HP::Decay::emu,          nonSMBR.BR_hjemu.at(0)*TotalWidths.at(0));
  h1.setDecayWidth(HP::Decay::etau,         nonSMBR.BR_hjetau.at(0)*TotalWidths.at(0));
  h1.setDecayWidth(HP::ChainDecay::W, "hp", nonSMBR.BR_hjHpW.at(0)[0]*TotalWidths.at(0));

  h2.setDecayWidth(HP::Decay::directInv,    nonSMBR.BR_hjinvisible.at(1)*TotalWidths.at(1));
  h2.setDecayWidth("h1","h1",               nonSMBR.BR_hkhjhi.at(1)[0][0]*TotalWidths.at(1));
  h2.setDecayWidth(HP::ChainDecay::Z, "h1", nonSMBR.BR_hjhiZ.at(1)[0]*TotalWidths.at(1));
  h2.setDecayWidth(HP::Decay::emu,          nonSMBR.BR_hjemu.at(1)*TotalWidths.at(1));
  h2.setDecayWidth(HP::Decay::etau,         nonSMBR.BR_hjetau.at(1)*TotalWidths.at(1));
  h2.setDecayWidth(HP::ChainDecay::W, "hp", nonSMBR.BR_hjHpW.at(1)[0]*TotalWidths.at(1));

  a1.setDecayWidth(HP::Decay::directInv,    nonSMBR.BR_hjinvisible.at(2)*TotalWidths.at(2));
  a1.setDecayWidth("h1","h1",               nonSMBR.BR_hkhjhi.at(2)[0][0]*TotalWidths.at(2));
  a1.setDecayWidth(HP::ChainDecay::Z, "h1", nonSMBR.BR_hjhiZ.at(2)[0]*TotalWidths.at(2));
  a1.setDecayWidth(HP::Decay::emu,          nonSMBR.BR_hjemu.at(2)*TotalWidths.at(2));
  a1.setDecayWidth(HP::Decay::etau,         nonSMBR.BR_hjetau.at(2)*TotalWidths.at(2));
  a1.setDecayWidth(HP::ChainDecay::W, "hp", nonSMBR.BR_hjHpW.at(2)[0]*TotalWidths.at(2));



}

//void HBHS::charged_input(THDM &model, Higgs::predictions::Predictions pred) const {
void HBHS::charged_input(THDM &model) const {


  DecayTable table{model};
  SM sm = model.get_SM();

  const std::array<double, nHplus> MHp{model.get_hmass(4)};
  const std::array<double, nHplus> HPlusGammaTot{table.get_gammatot_h(4)};
  const std::array<double, nHplus> CS_lep_HpjHmi_ratio{1.};
  const std::array<double, nHplus> BR_tHpjb{table.get_gamma_uhd(3, 4, 3) /
                                            table.get_gammatot_top()};
  const std::array<double, nHplus> BR_Hpjcs{table.get_gamma_hdu(4, 2, 2) /
                                            table.get_gammatot_h(4)};
  const std::array<double, nHplus> BR_Hpjcb{table.get_gamma_hdu(4, 3, 2) /
                                            table.get_gammatot_h(4)};
  const std::array<double, nHplus> BR_Hpjtaunu{table.get_gamma_hln(4, 3, 3) /
                                               table.get_gammatot_h(4)};
  const std::array<double, nHplus> BR_Hpjtb{table.get_gamma_hdu(4, 3, 3) /
                                            table.get_gammatot_h(4)};
  const std::array<double, nHplus> BR_Hpjud{table.get_gamma_hdu(4, 1, 1) /
                                            table.get_gammatot_h(4)};
  const std::array<double, nHplus> BR_Hpjus{table.get_gamma_hdu(4, 2, 1) /
                                            table.get_gammatot_h(4)};
  const std::array<double, nHplus> BR_Hpjub{table.get_gamma_hdu(4, 3, 1) /
                                            table.get_gammatot_h(4)};
  const std::array<double, nHplus> BR_Hpjcd{table.get_gamma_hdu(4, 1, 2) /
                                            table.get_gammatot_h(4)};
  const std::array<double, nHplus> BR_Hpjenu{table.get_gamma_hln(4, 1, 1) /
                                             table.get_gammatot_h(4)};
  const std::array<double, nHplus> BR_Hpjmunu{table.get_gamma_hln(4, 2, 2) /
                                              table.get_gammatot_h(4)};
  const std::array<double, nHplus> BR_HpWZ{0.};
  std::array<std::array<double, nHzero>, nHplus> BR_HphiW{};
  for (size_t i = 1; i <= nHzero; i++) {
    BR_HphiW[0][i - 1] = table.get_gamma_hvh(4, 3, i) / table.get_gammatot_h(4);
  }


  hp.setMass(MHp.at(0));
  hp.setTotalWidth(HPlusGammaTot.at(0));
  hp.setNormalizedCxn(HP::Collider::LEP, HP::Production::pair, CS_lep_HpjHmi_ratio.at(0), HP::ReferenceModel::SMHiggs);

  hp.setCxn(HP::Collider::LHC13,  HP::Production::brtHpb, BR_tHpjb.at(0));

  hp.setDecayWidth(HP::Decay::cs,    BR_Hpjcs.at(0)*HPlusGammaTot.at(0));
  hp.setDecayWidth(HP::Decay::cb,    BR_Hpjcb.at(0)*HPlusGammaTot.at(0));
  hp.setDecayWidth(HP::Decay::taunu, BR_Hpjtaunu.at(0)*HPlusGammaTot.at(0));
  hp.setDecayWidth(HP::Decay::tb,    BR_tHpjb.at(0)*HPlusGammaTot.at(0));
  hp.setDecayWidth(HP::Decay::WZ,    BR_HpWZ.at(0)*HPlusGammaTot.at(0));

  hp.setDecayWidth(HP::Decay::ud,    BR_Hpjud.at(0)*HPlusGammaTot.at(0));
  hp.setDecayWidth(HP::Decay::us,    BR_Hpjus.at(0)*HPlusGammaTot.at(0));
  hp.setDecayWidth(HP::Decay::ub,    BR_Hpjub.at(0)*HPlusGammaTot.at(0));
  hp.setDecayWidth(HP::Decay::cd,    BR_Hpjcd.at(0)*HPlusGammaTot.at(0));
  hp.setDecayWidth(HP::Decay::enu,   BR_Hpjenu.at(0)*HPlusGammaTot.at(0));
  hp.setDecayWidth(HP::Decay::munu,  BR_Hpjmunu.at(0)*HPlusGammaTot.at(0));





  // get charged Higgs couplings
  std::complex<double> cs, cp;
  model.get_coupling_hdu(4, 3, 3, cs, cp);
  const double mb = sm.get_dmass_MSbar(3);
  double mbrun = 0.;
  if (mb > 0) {
    mbrun = sm.run_qmass_MSbar(mb, mb, model.get_hmass(4), sm.get_qmass_pole(6),
                               sm.get_qmass_pole(5));
  }
  const double kappa_b =
      mb > 0 ? (cs + cp).imag() * sm.get_v() / (sqrt(2) * mbrun) : 0;

  const double mt = sm.get_umass_MSbar(3);
  const double mtrun = sm.run_qmass_MSbar(
      mt, mt, model.get_hmass(4), sm.get_qmass_pole(6), sm.get_qmass_pole(5));
  const double kappa_t = -(-cs + cp).imag() * sm.get_v() / (sqrt(2) * mtrun);
  if (debug){
    printf("Hccoups %16.8E %16.8E\n", kappa_t, kappa_b);
  }

  std::array<double, nHplus> CS_Hpmjtb{
      HCCS_tHc(model.get_hmass(4), kappa_t, kappa_b,
               table.get_gamma_uhd(3, 4, 3) / table.get_gammatot_top())};

  hp.setCxn(HP::Collider::LHC13, HP::Production::Hpmtb, CS_Hpmjtb.at(0));



}

void HBHS::set_mass_uncertainties(
    const std::array<double, HBHS::nHzero> &dMh,
    const std::array<double, HBHS::nHplus> &dMHp) {
  //HiggsBounds_set_mass_uncertainties(dMh.data(), dMHp.data());
  //HiggsSignals_neutral_input_MassUncertainty(dMh.data());
}

void HBResult::print() const {
  printf("\nHiggsBounds results:\n");
  printf("  Higgs  res  chan    \t\t\tratio   process\n");
  for (int i = 1; i <= 4; i++) {
    printf("%5d %5d %s  \t%16.8E %s   %s\n", i, result[i], channel[i].data(),
           obsratio[i], ncombined[i].data(), result[i] == 1 ? "Allowed" : "Excluded");
  }
  printf("------------------------------------------------------------\n");
  printf("  TOT %5d %s  \t%16.8E %s   %s\n", result[0], channel[0].data(), obsratio[0],
         ncombined[0].data(), result[0] == 1 ? "ALLOWED" : "EXCLUDED");
}

void HSResult::print() const {
  printf("\nHiggsSignals results:\n");
//  printf(" Chi^2 from rates: %16.8E\n", chisq_mu);
//  printf("  Chi^2 from mass: %16.8E\n", chisq_mh);
  printf("      Total chi^2: %16.8E\n", chisq);
  printf("    # observables: %16d\n\n", nobs);
}

#endif
