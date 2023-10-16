{
  const bool save = true;
  const TString saveDirectory = "/sbnd/app/users/lynnt/my_larsoft/v09_73_00/opt0_crumbs/nom_OpT0MeasuredPE";
  using namespace std;
  // gROOT->SetStyle("henrySBND");
  // gROOT->ForceStyle();

  //  TFile *file = new TFile("/sbnd/data/users/hlay/crumbs/training_data/crumbs_trees_combined.root","READ");
  // TFile *file = new TFile("/sbnd/data/users/hlay/crumbs/validation_data/crumbs_trees_combined.root","READ");
  TFile *file = new TFile("/sbnd/data/users/hlay/opt0_crumbs/production/combined_trees.root","READ");
  // TFile *file = new TFile("/sbnd/data/users/hlay/opt0_crumbs/production/intrnue_trees.root","READ");
  // TTree *all_slices = (TTree*) file->Get("crumbs/AllSliceTree");
  TTree *all_slices = (TTree*) file->Get("crumbsSCE/SliceTree");

  float nomReader_tpc_CRFracHitsInLongestTrack, nomReader_tpc_CRLongestTrackDeflection, nomReader_tpc_CRLongestTrackDirY, nomReader_tpc_CRNHitsMax,
    nomReader_tpc_NuEigenRatioInSphere, nomReader_tpc_NuNFinalStatePfos, nomReader_tpc_NuNHitsTotal, nomReader_tpc_NuNSpacePointsInSphere, nomReader_tpc_NuVertexY, nomReader_tpc_NuWeightedDirZ,
    nomReader_tpc_StoppingChi2Pol0, nomReader_tpc_StoppingChi2Exp, nomReader_tpc_StoppingChi2Ratio, nomReader_tpc_StoppingChi2CosmicPol0, nomReader_tpc_StoppingChi2CosmicExp, nomReader_tpc_StoppingChi2CosmicRatio,
    nomReader_pds_FMTotalScore, nomReader_pds_FMYScore, nomReader_pds_FMZScore, nomReader_pds_FMRRScore, nomReader_pds_FMRatioScore, nomReader_pds_FMPE, nomReader_pds_FMTime, nomReader_crt_TrackScore, nomReader_crt_HitScore,
    nomReader_crt_TrackTime, nomReader_crt_HitTime, nomReader_crt_nTrackMatches, nomReader_crt_nHitMatches;

  float optReader_tpc_CRFracHitsInLongestTrack, optReader_tpc_CRLongestTrackDeflection, optReader_tpc_CRLongestTrackDirY, optReader_tpc_CRNHitsMax,
    optReader_tpc_NuEigenRatioInSphere, optReader_tpc_NuNFinalStatePfos, optReader_tpc_NuNHitsTotal, optReader_tpc_NuNSpacePointsInSphere, optReader_tpc_NuVertexY, optReader_tpc_NuWeightedDirZ,
    optReader_tpc_StoppingChi2Pol0, optReader_tpc_StoppingChi2Exp, optReader_tpc_StoppingChi2Ratio, optReader_tpc_StoppingChi2CosmicPol0, optReader_tpc_StoppingChi2CosmicExp, optReader_tpc_StoppingChi2CosmicRatio,
    optReader_pds_OpT0Score, 
    // optReader_pds_OpT0Time,
    optReader_pds_OpT0MeasuredPE, 
    // optReader_pds_OpT0HypothesisPE,  
    optReader_crt_TrackScore, optReader_crt_HitScore,
    optReader_crt_TrackTime, optReader_crt_HitTime, optReader_crt_nTrackMatches, optReader_crt_nHitMatches;

  float tpc_NuScore, tpc_CRFracHitsInLongestTrack, tpc_CRLongestTrackDeflection, tpc_CRLongestTrackDirY, tpc_CRNHitsMax,
    tpc_NuEigenRatioInSphere, tpc_NuNFinalStatePfos, tpc_NuNHitsTotal, tpc_NuNSpacePointsInSphere, tpc_NuVertexY, tpc_NuWeightedDirZ,
    tpc_StoppingChi2Pol0, tpc_StoppingChi2Exp, tpc_StoppingChi2Ratio, tpc_StoppingChi2CosmicPol0, tpc_StoppingChi2CosmicExp, tpc_StoppingChi2CosmicRatio,
    pds_FMTotalScore, pds_FMYScore, pds_FMZScore, pds_FMRRScore, pds_FMRatioScore, pds_FMPE, pds_FMTime,
    pds_OpT0Score, pds_OpT0Time, pds_OpT0MeasuredPE, pds_OpT0HypothesisPE, 
    crt_TrackScore, crt_HitScore,
    crt_TrackTime, crt_HitTime;
  int crt_nTrackMatches, crt_nHitMatches;

  unsigned eventID, subRunID, runID, slicePDG, sliceIndex, matchedIndex;
  std::string *matchedType;
  double matchedPurity, matchedCompleteness;
 
  all_slices->SetBranchAddress("tpc_NuScore",&tpc_NuScore);
  all_slices->SetBranchAddress("tpc_CRFracHitsInLongestTrack",&tpc_CRFracHitsInLongestTrack);
  all_slices->SetBranchAddress("tpc_CRLongestTrackDeflection",&tpc_CRLongestTrackDeflection);
  all_slices->SetBranchAddress("tpc_CRLongestTrackDirY",&tpc_CRLongestTrackDirY);
  all_slices->SetBranchAddress("tpc_CRNHitsMax",&tpc_CRNHitsMax);
  all_slices->SetBranchAddress("tpc_NuEigenRatioInSphere",&tpc_NuEigenRatioInSphere);
  all_slices->SetBranchAddress("tpc_NuNFinalStatePfos",&tpc_NuNFinalStatePfos);
  all_slices->SetBranchAddress("tpc_NuNHitsTotal",&tpc_NuNHitsTotal);
  all_slices->SetBranchAddress("tpc_NuNSpacePointsInSphere",&tpc_NuNSpacePointsInSphere);
  all_slices->SetBranchAddress("tpc_NuVertexY",&tpc_NuVertexY);
  all_slices->SetBranchAddress("tpc_NuWeightedDirZ",&tpc_NuWeightedDirZ);
  all_slices->SetBranchAddress("tpc_StoppingChi2CosmicRatio",&tpc_StoppingChi2CosmicRatio);

  all_slices->SetBranchAddress("pds_FMTotalScore",&pds_FMTotalScore);
  all_slices->SetBranchAddress("pds_FMPE",&pds_FMPE);
  all_slices->SetBranchAddress("pds_FMTime",&pds_FMTime);

  all_slices->SetBranchAddress("pds_OpT0Score",&pds_OpT0Score);
  all_slices->SetBranchAddress("pds_OpT0Time",&pds_OpT0Time);
  all_slices->SetBranchAddress("pds_OpT0MeasuredPE",&pds_OpT0MeasuredPE);
  all_slices->SetBranchAddress("pds_OpT0HypothesisPE",&pds_OpT0HypothesisPE);

  all_slices->SetBranchAddress("crt_TrackScore",&crt_TrackScore);
  all_slices->SetBranchAddress("crt_HitScore",&crt_HitScore);
  all_slices->SetBranchAddress("crt_TrackTime",&crt_TrackTime);
  all_slices->SetBranchAddress("crt_HitTime",&crt_HitTime);

  all_slices->SetBranchAddress("eventID",&eventID);
  all_slices->SetBranchAddress("subRunID",&subRunID);
  all_slices->SetBranchAddress("runID",&runID);
  all_slices->SetBranchAddress("slicePDG",&slicePDG);
  all_slices->SetBranchAddress("matchedIndex",&matchedIndex);
  all_slices->SetBranchAddress("matchedType",&matchedType);
  all_slices->SetBranchAddress("matchedPurity",&matchedPurity);
  all_slices->SetBranchAddress("matchedCompleteness",&matchedCompleteness);

  TMVA::Reader *nomReader = new TMVA::Reader("!Color:!Silent");
  nomReader->AddVariable("tpc_CRFracHitsInLongestTrack",&nomReader_tpc_CRFracHitsInLongestTrack);
  nomReader->AddVariable("tpc_CRLongestTrackDeflection",&nomReader_tpc_CRLongestTrackDeflection);
  nomReader->AddVariable("tpc_CRLongestTrackDirY",&nomReader_tpc_CRLongestTrackDirY);
  nomReader->AddVariable("tpc_CRNHitsMax",&nomReader_tpc_CRNHitsMax);
  nomReader->AddVariable("tpc_NuEigenRatioInSphere",&nomReader_tpc_NuEigenRatioInSphere);
  nomReader->AddVariable("tpc_NuNFinalStatePfos",&nomReader_tpc_NuNFinalStatePfos);
  nomReader->AddVariable("tpc_NuNHitsTotal",&nomReader_tpc_NuNHitsTotal);
  nomReader->AddVariable("tpc_NuNSpacePointsInSphere",&nomReader_tpc_NuNSpacePointsInSphere);
  nomReader->AddVariable("tpc_NuVertexY",&nomReader_tpc_NuVertexY);
  nomReader->AddVariable("tpc_NuWeightedDirZ",&nomReader_tpc_NuWeightedDirZ);
  nomReader->AddVariable("tpc_StoppingChi2CosmicRatio",&nomReader_tpc_StoppingChi2CosmicRatio);
  
  nomReader->AddVariable("pds_FMTotalScore",&nomReader_pds_FMTotalScore);
  nomReader->AddVariable("pds_FMPE",&nomReader_pds_FMPE);
  nomReader->AddVariable("pds_FMTime",&nomReader_pds_FMTime);

  nomReader->AddVariable("crt_TrackScore",&nomReader_crt_TrackScore);
  nomReader->AddVariable("crt_HitScore",&nomReader_crt_HitScore);
  nomReader->AddVariable("crt_TrackTime",&nomReader_crt_TrackTime);
  nomReader->AddVariable("crt_HitTime",&nomReader_crt_HitTime);

  nomReader->BookMVA("BDT Method", "/sbnd/data/users/hlay/opt0_crumbs/training/CRUMBS_Nominal/weights/CrumbsTMVAClassification_BDT.weights.xml");

  TMVA::Reader *optReader = new TMVA::Reader("!Color:!Silent");
  optReader->AddVariable("tpc_CRFracHitsInLongestTrack",&optReader_tpc_CRFracHitsInLongestTrack);
  optReader->AddVariable("tpc_CRLongestTrackDeflection",&optReader_tpc_CRLongestTrackDeflection);
  optReader->AddVariable("tpc_CRLongestTrackDirY",&optReader_tpc_CRLongestTrackDirY);
  optReader->AddVariable("tpc_CRNHitsMax",&optReader_tpc_CRNHitsMax);
  optReader->AddVariable("tpc_NuEigenRatioInSphere",&optReader_tpc_NuEigenRatioInSphere);
  optReader->AddVariable("tpc_NuNFinalStatePfos",&optReader_tpc_NuNFinalStatePfos);
  optReader->AddVariable("tpc_NuNHitsTotal",&optReader_tpc_NuNHitsTotal);
  optReader->AddVariable("tpc_NuNSpacePointsInSphere",&optReader_tpc_NuNSpacePointsInSphere);
  optReader->AddVariable("tpc_NuVertexY",&optReader_tpc_NuVertexY);
  optReader->AddVariable("tpc_NuWeightedDirZ",&optReader_tpc_NuWeightedDirZ);
  optReader->AddVariable("tpc_StoppingChi2CosmicRatio",&optReader_tpc_StoppingChi2CosmicRatio);
  
  optReader->AddVariable("pds_OpT0Score",&optReader_pds_OpT0Score);
  // optReader->AddVariable("pds_OpT0Time",&optReader_pds_OpT0Time);
  optReader->AddVariable("pds_OpT0MeasuredPE",&optReader_pds_OpT0MeasuredPE);
  // optReader->AddVariable("pds_OpT0HypothesisPE",&optReader_pds_OpT0HypothesisPE);

  optReader->AddVariable("crt_TrackScore",&optReader_crt_TrackScore);
  optReader->AddVariable("crt_HitScore",&optReader_crt_HitScore);
  optReader->AddVariable("crt_TrackTime",&optReader_crt_TrackTime);
  optReader->AddVariable("crt_HitTime",&optReader_crt_HitTime);

  // optReader->BookMVA("BDT Method", "/sbnd/data/users/hlay/crumbs/TMVA/CrumbsDataset_41/weights/CrumbsTMVAClassification_BDT.weights.xml");
  // optReader->BookMVA("BDT Method", "/sbnd/data/users/hlay/opt0_crumbs/training/CRUMBS_AllOpT0/weights/CrumbsTMVAClassification_BDT.weights.xml");
  // optReader->BookMVA("BDT Method", "/sbnd/data/users/hlay/opt0_crumbs/training/CRUMBS_OpT0Score/weights/CrumbsTMVAClassification_BDT.weights.xml");
  optReader->BookMVA("BDT Method", "/sbnd/data/users/hlay/opt0_crumbs/training/CRUMBS_OpT0MeasuredPE/weights/CrumbsTMVAClassification_BDT.weights.xml");

  int N_all_slices = all_slices->GetEntries();

  TH1F *hNuSliceNuScore = new TH1F("hNuSliceNuScore",";Nu Score;Slices",50, 0, 1);
  TH1F *hNotNuSliceNuScore = new TH1F("hNotNuSliceNuScore",";Nu Score;Slices",50, 0, 1);
  TH1F *hNuSliceNomCrumbsScore = new TH1F("hNuSliceNomCrumbsScore",";CRUMBS Score;Slices",      50, -.7, .6);
  TH1F *hNotNuSliceNomCrumbsScore = new TH1F("hNotNuSliceNomCrumbsScore",";CRUMBS Score;Slices",50, -.7, .6);
  TH1F *hNuSliceOptCrumbsScore = new TH1F("hNuSliceOptCrumbsScore",";CRUMBS Score;Slices",      50, -.7, .6);
  TH1F *hNotNuSliceOptCrumbsScore = new TH1F("hNotNuSliceOptCrumbsScore",";CRUMBS Score;Slices",50, -.7, .6);
  
  for(int all_i = 0; all_i < N_all_slices; ++all_i){
    all_slices->GetEntry(all_i);
    nomReader_tpc_CRFracHitsInLongestTrack = tpc_CRFracHitsInLongestTrack;
    nomReader_tpc_CRLongestTrackDeflection = tpc_CRLongestTrackDeflection;
    nomReader_tpc_CRLongestTrackDirY = tpc_CRLongestTrackDirY;
    nomReader_tpc_CRNHitsMax = tpc_CRNHitsMax;
    nomReader_tpc_NuEigenRatioInSphere = tpc_NuEigenRatioInSphere;
    nomReader_tpc_NuNFinalStatePfos = tpc_NuNFinalStatePfos;
    nomReader_tpc_NuNHitsTotal = tpc_NuNHitsTotal;
    nomReader_tpc_NuNSpacePointsInSphere = tpc_NuNSpacePointsInSphere;
    nomReader_tpc_NuVertexY = tpc_NuVertexY; 
    nomReader_tpc_NuWeightedDirZ = tpc_NuWeightedDirZ;
    nomReader_tpc_StoppingChi2Pol0 = tpc_StoppingChi2Pol0;
    nomReader_tpc_StoppingChi2Exp = tpc_StoppingChi2Exp;
    nomReader_tpc_StoppingChi2Ratio = tpc_StoppingChi2Ratio;
    nomReader_tpc_StoppingChi2CosmicPol0 = tpc_StoppingChi2CosmicPol0;
    nomReader_tpc_StoppingChi2CosmicExp = tpc_StoppingChi2CosmicExp;
    nomReader_tpc_StoppingChi2CosmicRatio = tpc_StoppingChi2CosmicRatio;
    nomReader_pds_FMTotalScore = pds_FMTotalScore;
    nomReader_pds_FMYScore = pds_FMYScore;
    nomReader_pds_FMZScore = pds_FMZScore;
    nomReader_pds_FMRRScore = pds_FMRRScore;
    nomReader_pds_FMRatioScore = pds_FMRatioScore;
    nomReader_pds_FMPE = pds_FMPE;
    nomReader_pds_FMTime = pds_FMTime;
    nomReader_crt_TrackScore = crt_TrackScore;
    nomReader_crt_HitScore = crt_HitScore;
    nomReader_crt_TrackTime = crt_TrackTime;
    nomReader_crt_HitTime = crt_HitTime;
    nomReader_crt_nTrackMatches = crt_nTrackMatches;
    nomReader_crt_nHitMatches = crt_nHitMatches;

    float bdtscore = nomReader->EvaluateMVA("BDT Method");

    if(*matchedType == "Nu" && matchedPurity > 0.8 && matchedCompleteness > 0.8)
      hNuSliceNomCrumbsScore->Fill(bdtscore);
    else if(*matchedType != "Nu")
      hNotNuSliceNomCrumbsScore->Fill(bdtscore); 
  }

  for(int all_i = 0; all_i < N_all_slices; ++all_i){
    all_slices->GetEntry(all_i);
      if(*matchedType == "Nu" && matchedPurity > 0.8 && matchedCompleteness > 0.8)
      hNuSliceNuScore->Fill(tpc_NuScore);
    else if(*matchedType != "Nu")
      hNotNuSliceNuScore->Fill(tpc_NuScore);
    optReader_tpc_CRFracHitsInLongestTrack = tpc_CRFracHitsInLongestTrack;
    optReader_tpc_CRLongestTrackDeflection = tpc_CRLongestTrackDeflection;
    optReader_tpc_CRLongestTrackDirY = tpc_CRLongestTrackDirY;
    optReader_tpc_CRNHitsMax = tpc_CRNHitsMax;
    optReader_tpc_NuEigenRatioInSphere = tpc_NuEigenRatioInSphere;
    optReader_tpc_NuNFinalStatePfos = tpc_NuNFinalStatePfos;
    optReader_tpc_NuNHitsTotal = tpc_NuNHitsTotal;
    optReader_tpc_NuNSpacePointsInSphere = tpc_NuNSpacePointsInSphere;
    optReader_tpc_NuVertexY = tpc_NuVertexY; 
    optReader_tpc_NuWeightedDirZ = tpc_NuWeightedDirZ;
    optReader_tpc_StoppingChi2Pol0 = tpc_StoppingChi2Pol0;
    optReader_tpc_StoppingChi2Exp = tpc_StoppingChi2Exp;
    optReader_tpc_StoppingChi2Ratio = tpc_StoppingChi2Ratio;
    optReader_tpc_StoppingChi2CosmicPol0 = tpc_StoppingChi2CosmicPol0;
    optReader_tpc_StoppingChi2CosmicExp = tpc_StoppingChi2CosmicExp;
    optReader_tpc_StoppingChi2CosmicRatio = tpc_StoppingChi2CosmicRatio;
    optReader_pds_OpT0Score = pds_OpT0Score;
    // optReader_pds_OpT0Time = pds_OpT0Time;
    optReader_pds_OpT0MeasuredPE = pds_OpT0MeasuredPE;
    // optReader_pds_OpT0HypothesisPE = pds_OpT0HypothesisPE;
    optReader_crt_TrackScore = crt_TrackScore;
    optReader_crt_HitScore = crt_HitScore;
    optReader_crt_TrackTime = crt_TrackTime;
    optReader_crt_HitTime = crt_HitTime;
    optReader_crt_nTrackMatches = crt_nTrackMatches;
    optReader_crt_nHitMatches = crt_nHitMatches;

    float bdtscore = optReader->EvaluateMVA("BDT Method");

    if(*matchedType == "Nu" && matchedPurity > 0.8 && matchedCompleteness > 0.8)
      hNuSliceOptCrumbsScore->Fill(bdtscore);
    else if(*matchedType != "Nu")
      hNotNuSliceOptCrumbsScore->Fill(bdtscore); 
  }

  float nuEff[50], nuRej[50];
  int nuSignalSum = 0, nuBackSum = 0;
  const int nuSignalTotal = hNuSliceNuScore->GetEntries();
  const int nuBackTotal = hNotNuSliceNuScore->GetEntries();

  float nomcrumbsEff[50], nomcrumbsRej[50];
  int nomcrumbsSignalSum = 0, nomcrumbsBackSum = 0;
  const int nomcrumbsSignalTotal = hNuSliceNomCrumbsScore->GetEntries();
  const int nomcrumbsBackTotal = hNotNuSliceNomCrumbsScore->GetEntries();
 
  float optcrumbsEff[50], optcrumbsRej[50];
  int optcrumbsSignalSum = 0, optcrumbsBackSum = 0;
  const int optcrumbsSignalTotal = hNuSliceOptCrumbsScore->GetEntries();
  const int optcrumbsBackTotal = hNotNuSliceOptCrumbsScore->GetEntries();

  for(int i = 1; i < 51; ++i)
    {
      nuSignalSum += hNuSliceNuScore->GetBinContent(i);
      nuBackSum += hNotNuSliceNuScore->GetBinContent(i);

      nuEff[i-1] = (float) (nuSignalTotal - nuSignalSum) / (float) nuSignalTotal;
      nuRej[i-1] = (float) nuBackSum / (float) nuBackTotal;
    }      

  for(int i = 1; i < 51; ++i)
    {
      nomcrumbsSignalSum += hNuSliceNomCrumbsScore->GetBinContent(i);
      nomcrumbsBackSum += hNotNuSliceNomCrumbsScore->GetBinContent(i);

      nomcrumbsEff[i-1] = (float) (nomcrumbsSignalTotal - nomcrumbsSignalSum) / (float) nomcrumbsSignalTotal;
      nomcrumbsRej[i-1] = (float) nomcrumbsBackSum / (float) nomcrumbsBackTotal;
    }

  for(int i = 1; i < 51; ++i)
    {
      optcrumbsSignalSum += hNuSliceOptCrumbsScore->GetBinContent(i);
      optcrumbsBackSum += hNotNuSliceOptCrumbsScore->GetBinContent(i);

      optcrumbsEff[i-1] = (float) (optcrumbsSignalTotal - optcrumbsSignalSum) / (float) optcrumbsSignalTotal;
      optcrumbsRej[i-1] = (float) optcrumbsBackSum / (float) optcrumbsBackTotal;
    }

  TGraph *nuROC = new TGraph(50, nuRej, nuEff);
  TGraph *nomcrumbsROC = new TGraph(50, nomcrumbsRej, nomcrumbsEff);
  TGraph *optcrumbsROC = new TGraph(50, optcrumbsRej, optcrumbsEff);

  // TCanvas *cNuScore = new TCanvas("cNuScore","cNuScore");
  // cNuScore->cd();

  // hNuSliceNuScore->SetLineColor(kBlue+2);
  // hNotNuSliceNuScore->SetLineColor(kRed+2);
  // hNotNuSliceNuScore->GetYaxis()->SetTitleOffset(1.2);
  // hNotNuSliceNuScore->GetYaxis()->SetNdivisions(505);
  // hNotNuSliceNuScore->GetXaxis()->SetNdivisions(505);

  // hNotNuSliceNuScore->Draw("same");
  // hNuSliceNuScore->Draw("same");

  TLegend *leg = new TLegend(0.6,0.75,0.9,0.9);
  leg->AddEntry(hNuSliceOptCrumbsScore, "OpT0 (Score + Meas PE) #nu Slices", "l");
  leg->AddEntry(hNotNuSliceOptCrumbsScore, "OpT0 (Score + Meas PE) Other Slices", "l");
  leg->AddEntry(hNuSliceNomCrumbsScore, "Nominal #nu Slices", "l");
  leg->AddEntry(hNotNuSliceNomCrumbsScore, "Nominal Other Slices", "l");
  // leg->Draw();

  // TLine *line_nu = new TLine(0.4,0,0.4,1.04*hNotNuSliceNuScore->GetMaximum());
  // line_nu->SetLineColor(kGreen+3);
  // line_nu->SetLineWidth(5);
  // line_nu->SetLineStyle(9);
  //  line_nu->Draw();

  // if(save)
  //   {
  //     cNuScore->SaveAs(saveDirectory + "/nu_score_quality.png");
  //     cNuScore->SaveAs(saveDirectory + "/nu_score_quality.pdf");
  //   }

  TCanvas *cCrumbsScore = new TCanvas("cCrumbsScore","cCrumbsScore");
  cCrumbsScore->cd();

  hNuSliceOptCrumbsScore->SetLineColor(kTeal+2);
  hNuSliceOptCrumbsScore->SetLineWidth(3);
  hNotNuSliceOptCrumbsScore->SetLineColor(kMagenta+2);
  hNotNuSliceOptCrumbsScore->SetLineWidth(3);
  hNotNuSliceOptCrumbsScore->GetYaxis()->SetTitleOffset(1.2);
  hNotNuSliceOptCrumbsScore->GetYaxis()->SetNdivisions(505);
  hNotNuSliceOptCrumbsScore->GetXaxis()->SetNdivisions(505);
  hNotNuSliceOptCrumbsScore->Draw("same");
  hNuSliceOptCrumbsScore->Draw("same");

  hNuSliceNomCrumbsScore->SetLineColor(kBlue+2);
  hNuSliceNomCrumbsScore->SetLineWidth(3);
  hNuSliceNomCrumbsScore->SetLineStyle(2);
  hNotNuSliceNomCrumbsScore->SetLineColor(kRed+2);
  hNotNuSliceNomCrumbsScore->SetLineWidth(3);
  hNotNuSliceNomCrumbsScore->SetLineStyle(2);
  hNotNuSliceNomCrumbsScore->Draw("same");
  hNuSliceNomCrumbsScore->Draw("same");

  leg->Draw();

  TLine *line_crumbs = new TLine(-0.05,0,-0.05,1.04*hNotNuSliceNomCrumbsScore->GetMaximum());
  line_crumbs->SetLineColor(kGreen+3);
  line_crumbs->SetLineWidth(5);
  line_crumbs->SetLineStyle(9);
  //  line_crumbs->Draw();

  TLine *line_crumbs_0 = new TLine(-0.1,0,-0.1,1.04*hNotNuSliceNomCrumbsScore->GetMaximum());
  line_crumbs_0->SetLineColor(kGray+3);
  line_crumbs_0->SetLineWidth(3);
  line_crumbs_0->SetLineStyle(1);
  // line_crumbs_0->Draw();

  TLine *line_crumbs_1 = new TLine(0,0,0,1.04*hNotNuSliceNomCrumbsScore->GetMaximum());
  line_crumbs_1->SetLineColor(kGray+3);
  line_crumbs_1->SetLineWidth(3);
  line_crumbs_1->SetLineStyle(1);
  // line_crumbs_1->Draw();

 TLine *line_crumbs_2 = new TLine(0.05,0,0.05,1.04*hNotNuSliceNomCrumbsScore->GetMaximum());
  line_crumbs_2->SetLineColor(kGray+3);
  line_crumbs_2->SetLineWidth(3);
  line_crumbs_2->SetLineStyle(1);
  // line_crumbs_2->Draw();

  TLine *line_crumbs_3 = new TLine(0.1,0,0.1,1.04*hNotNuSliceNomCrumbsScore->GetMaximum());
  line_crumbs_3->SetLineColor(kGray+3);
  line_crumbs_3->SetLineWidth(3);
  line_crumbs_3->SetLineStyle(1);
  // line_crumbs_3->Draw();


  if(save)
    {
      cCrumbsScore->SaveAs(saveDirectory + "/crumbs_score_quality.png");
      cCrumbsScore->SaveAs(saveDirectory + "/crumbs_score_quality.pdf");
    }

  TCanvas *cROC = new TCanvas("cROC","cROC");
  cROC->SetGrid();
  cROC->cd();

  gStyle->SetLabelSize(0.03,"xy");

  TMultiGraph *multi = new TMultiGraph();
  nuROC->SetLineColor(kCyan-3);
  nuROC->SetLineWidth(3);
  nuROC->SetMarkerColor(kCyan-3);
  nuROC->SetMarkerStyle(3);
  multi->Add(nuROC);
  nomcrumbsROC->SetLineColor(kViolet-5);
  nomcrumbsROC->SetLineWidth(3);
  nomcrumbsROC->SetMarkerColor(kViolet-5);
  nomcrumbsROC->SetMarkerStyle(3);
  multi->Add(nomcrumbsROC);
  optcrumbsROC->SetLineColor(kAzure-5);
  optcrumbsROC->SetLineWidth(3);
  optcrumbsROC->SetMarkerColor(kAzure-5);
  optcrumbsROC->SetMarkerStyle(3);
  multi->Add(optcrumbsROC);
  multi->SetTitle(";Signal Efficiency;Background Rejection");
  multi->Draw("ACP");
  // multi->GetXaxis()->SetLimits(0.9,1.);
  // multi->GetYaxis()->SetLimits(0.9,1.);
  multi->GetXaxis()->SetNdivisions(20);
  multi->GetYaxis()->SetNdivisions(20);

  TLegend *rocLeg = new TLegend(.25,.25,.55,.45);
  rocLeg->AddEntry(nuROC,"Nu Score","lp");
  rocLeg->AddEntry(nomcrumbsROC,"Nom CRUMBS Score","lp");
  rocLeg->AddEntry(optcrumbsROC,"OpT0 (Score + Meas PE) CRUMBS Score","lp");
  rocLeg->Draw();

  if(save)
    {
      cROC->SaveAs(saveDirectory + "/roc_quality.png"); 
      cROC->SaveAs(saveDirectory + "/roc_quality.pdf");
    }

  TCanvas *cROC_zoomed = new TCanvas("cROC_zoomed","cROC_zoomed");
  cROC_zoomed->SetGrid();
  cROC_zoomed->cd();

  gStyle->SetLabelSize(0.03,"xy");

  TMultiGraph *multi_zoomed = new TMultiGraph();
  nuROC->SetLineColor(kCyan-3);
  nuROC->SetLineWidth(3);
  nuROC->SetMarkerColor(kCyan-3);
  nuROC->SetMarkerStyle(3);
  multi_zoomed->Add(nuROC);
  nomcrumbsROC->SetLineColor(kViolet-5);
  nomcrumbsROC->SetLineWidth(3);
  nomcrumbsROC->SetMarkerColor(kViolet-5);
  nomcrumbsROC->SetMarkerStyle(3);
  multi_zoomed->Add(nomcrumbsROC);
  optcrumbsROC->SetLineColor(kAzure-5);
  optcrumbsROC->SetLineWidth(3);
  optcrumbsROC->SetMarkerColor(kAzure-5);
  optcrumbsROC->SetMarkerStyle(3);
  multi_zoomed->Add(optcrumbsROC);
  multi_zoomed->SetTitle(";Signal Efficiency;Background Rejection");
  multi_zoomed->Draw("ACP");
  multi_zoomed->GetXaxis()->SetLimits(0.8,1.);
  multi_zoomed->SetMinimum(0.8);
  multi_zoomed->SetMaximum(1.);

  TLegend *rocZoomLeg = new TLegend(.25,.25,.55,.45);
  rocZoomLeg->AddEntry(nuROC,"Nu Score","lp");
  rocZoomLeg->AddEntry(nomcrumbsROC,"Nom CRUMBS Score","lp");
  rocZoomLeg->AddEntry(optcrumbsROC,"OpT0 (Score + Meas PE) CRUMBS Score","lp");
  rocZoomLeg->Draw();

  if(save)
    {
      cROC_zoomed->SaveAs(saveDirectory + "/roc_quality_zoomed.png");
      cROC_zoomed->SaveAs(saveDirectory + "/roc_quality_zoomed.pdf");
    }
}
