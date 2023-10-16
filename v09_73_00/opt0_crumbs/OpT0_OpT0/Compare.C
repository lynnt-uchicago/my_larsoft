{
  const bool save = true;
  const TString saveDirectory = "/sbnd/app/users/lynnt/my_larsoft/v09_73_00/opt0_crumbs/OpT0_OpT0";
  using namespace std;
  // gROOT->SetStyle("henrySBND");
  // gROOT->ForceStyle();

  TFile *file = new TFile("/sbnd/data/users/hlay/opt0_crumbs/production/combined_trees.root","READ");
  TTree *all_slices = (TTree*) file->Get("crumbsSCE/SliceTree");

  // all vars
  float optReader_tpc_CRFracHitsInLongestTrack, optReader_tpc_CRLongestTrackDeflection, optReader_tpc_CRLongestTrackDirY, optReader_tpc_CRNHitsMax,
    optReader_tpc_NuEigenRatioInSphere, optReader_tpc_NuNFinalStatePfos, optReader_tpc_NuNHitsTotal, optReader_tpc_NuNSpacePointsInSphere, optReader_tpc_NuVertexY, optReader_tpc_NuWeightedDirZ,
    optReader_tpc_StoppingChi2Pol0, optReader_tpc_StoppingChi2Exp, optReader_tpc_StoppingChi2Ratio, optReader_tpc_StoppingChi2CosmicPol0, optReader_tpc_StoppingChi2CosmicExp, optReader_tpc_StoppingChi2CosmicRatio,
    optReader_pds_OpT0Score, 
    optReader_pds_OpT0Time,
    optReader_pds_OpT0MeasuredPE, 
    optReader_pds_OpT0HypothesisPE,  
    optReader_crt_TrackScore, optReader_crt_HitScore,
    optReader_crt_TrackTime, optReader_crt_HitTime, optReader_crt_nTrackMatches, optReader_crt_nHitMatches;

  // score only
  float optScoreReader_tpc_CRFracHitsInLongestTrack, optScoreReader_tpc_CRLongestTrackDeflection, optScoreReader_tpc_CRLongestTrackDirY, optScoreReader_tpc_CRNHitsMax,
    optScoreReader_tpc_NuEigenRatioInSphere, optScoreReader_tpc_NuNFinalStatePfos, optScoreReader_tpc_NuNHitsTotal, optScoreReader_tpc_NuNSpacePointsInSphere, optScoreReader_tpc_NuVertexY, optScoreReader_tpc_NuWeightedDirZ,
    optScoreReader_tpc_StoppingChi2Pol0, optScoreReader_tpc_StoppingChi2Exp, optScoreReader_tpc_StoppingChi2Ratio, optScoreReader_tpc_StoppingChi2CosmicPol0, optScoreReader_tpc_StoppingChi2CosmicExp, optScoreReader_tpc_StoppingChi2CosmicRatio,
    optScoreReader_pds_OpT0Score, 
    optScoreReader_crt_TrackScore, optScoreReader_crt_HitScore,
    optScoreReader_crt_TrackTime, optScoreReader_crt_HitTime, optScoreReader_crt_nTrackMatches, optScoreReader_crt_nHitMatches;

  // score + meas PE 
  float optMeasReader_tpc_CRFracHitsInLongestTrack, optMeasReader_tpc_CRLongestTrackDeflection, optMeasReader_tpc_CRLongestTrackDirY, optMeasReader_tpc_CRNHitsMax,
    optMeasReader_tpc_NuEigenRatioInSphere, optMeasReader_tpc_NuNFinalStatePfos, optMeasReader_tpc_NuNHitsTotal, optMeasReader_tpc_NuNSpacePointsInSphere, optMeasReader_tpc_NuVertexY, optMeasReader_tpc_NuWeightedDirZ,
    optMeasReader_tpc_StoppingChi2Pol0, optMeasReader_tpc_StoppingChi2Exp, optMeasReader_tpc_StoppingChi2Ratio, optMeasReader_tpc_StoppingChi2CosmicPol0, optMeasReader_tpc_StoppingChi2CosmicExp, optMeasReader_tpc_StoppingChi2CosmicRatio,
    optMeasReader_pds_OpT0Score, 
    optMeasReader_pds_OpT0MeasuredPE, 
    optMeasReader_crt_TrackScore, optMeasReader_crt_HitScore,
    optMeasReader_crt_TrackTime, optMeasReader_crt_HitTime, optMeasReader_crt_nTrackMatches, optMeasReader_crt_nHitMatches;

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
  optReader->AddVariable("pds_OpT0Time",&optReader_pds_OpT0Time);
  optReader->AddVariable("pds_OpT0MeasuredPE",&optReader_pds_OpT0MeasuredPE);
  optReader->AddVariable("pds_OpT0HypothesisPE",&optReader_pds_OpT0HypothesisPE);
  optReader->AddVariable("crt_TrackScore",&optReader_crt_TrackScore);
  optReader->AddVariable("crt_HitScore",&optReader_crt_HitScore);
  optReader->AddVariable("crt_TrackTime",&optReader_crt_TrackTime);
  optReader->AddVariable("crt_HitTime",&optReader_crt_HitTime);
  optReader->BookMVA("BDT Method", "/sbnd/data/users/hlay/opt0_crumbs/training/CRUMBS_AllOpT0/weights/CrumbsTMVAClassification_BDT.weights.xml");


  TMVA::Reader *optScoreReader = new TMVA::Reader("!Color:!Silent");
  optScoreReader->AddVariable("tpc_CRFracHitsInLongestTrack",&optScoreReader_tpc_CRFracHitsInLongestTrack);
  optScoreReader->AddVariable("tpc_CRLongestTrackDeflection",&optScoreReader_tpc_CRLongestTrackDeflection);
  optScoreReader->AddVariable("tpc_CRLongestTrackDirY",&optScoreReader_tpc_CRLongestTrackDirY);
  optScoreReader->AddVariable("tpc_CRNHitsMax",&optScoreReader_tpc_CRNHitsMax);
  optScoreReader->AddVariable("tpc_NuEigenRatioInSphere",&optScoreReader_tpc_NuEigenRatioInSphere);
  optScoreReader->AddVariable("tpc_NuNFinalStatePfos",&optScoreReader_tpc_NuNFinalStatePfos);
  optScoreReader->AddVariable("tpc_NuNHitsTotal",&optScoreReader_tpc_NuNHitsTotal);
  optScoreReader->AddVariable("tpc_NuNSpacePointsInSphere",&optScoreReader_tpc_NuNSpacePointsInSphere);
  optScoreReader->AddVariable("tpc_NuVertexY",&optScoreReader_tpc_NuVertexY);
  optScoreReader->AddVariable("tpc_NuWeightedDirZ",&optScoreReader_tpc_NuWeightedDirZ);
  optScoreReader->AddVariable("tpc_StoppingChi2CosmicRatio",&optScoreReader_tpc_StoppingChi2CosmicRatio);
  optScoreReader->AddVariable("pds_OpT0Score",&optScoreReader_pds_OpT0Score);
  optScoreReader->AddVariable("crt_TrackScore",&optScoreReader_crt_TrackScore);
  optScoreReader->AddVariable("crt_HitScore",&optScoreReader_crt_HitScore);
  optScoreReader->AddVariable("crt_TrackTime",&optScoreReader_crt_TrackTime);
  optScoreReader->AddVariable("crt_HitTime",&optScoreReader_crt_HitTime);
  optScoreReader->BookMVA("BDT Method", "/sbnd/data/users/hlay/opt0_crumbs/training/CRUMBS_OpT0Score/weights/CrumbsTMVAClassification_BDT.weights.xml");


  TMVA::Reader *optMeasReader = new TMVA::Reader("!Color:!Silent");
  optMeasReader->AddVariable("tpc_CRFracHitsInLongestTrack",&optMeasReader_tpc_CRFracHitsInLongestTrack);
  optMeasReader->AddVariable("tpc_CRLongestTrackDeflection",&optMeasReader_tpc_CRLongestTrackDeflection);
  optMeasReader->AddVariable("tpc_CRLongestTrackDirY",&optMeasReader_tpc_CRLongestTrackDirY);
  optMeasReader->AddVariable("tpc_CRNHitsMax",&optMeasReader_tpc_CRNHitsMax);
  optMeasReader->AddVariable("tpc_NuEigenRatioInSphere",&optMeasReader_tpc_NuEigenRatioInSphere);
  optMeasReader->AddVariable("tpc_NuNFinalStatePfos",&optMeasReader_tpc_NuNFinalStatePfos);
  optMeasReader->AddVariable("tpc_NuNHitsTotal",&optMeasReader_tpc_NuNHitsTotal);
  optMeasReader->AddVariable("tpc_NuNSpacePointsInSphere",&optMeasReader_tpc_NuNSpacePointsInSphere);
  optMeasReader->AddVariable("tpc_NuVertexY",&optMeasReader_tpc_NuVertexY);
  optMeasReader->AddVariable("tpc_NuWeightedDirZ",&optMeasReader_tpc_NuWeightedDirZ);
  optMeasReader->AddVariable("tpc_StoppingChi2CosmicRatio",&optMeasReader_tpc_StoppingChi2CosmicRatio);
  optMeasReader->AddVariable("pds_OpT0Score",&optMeasReader_pds_OpT0Score);
  optMeasReader->AddVariable("pds_OpT0MeasuredPE",&optMeasReader_pds_OpT0MeasuredPE);
  optMeasReader->AddVariable("crt_TrackScore",&optMeasReader_crt_TrackScore);
  optMeasReader->AddVariable("crt_HitScore",&optMeasReader_crt_HitScore);
  optMeasReader->AddVariable("crt_TrackTime",&optMeasReader_crt_TrackTime);
  optMeasReader->AddVariable("crt_HitTime",&optMeasReader_crt_HitTime);
  optMeasReader->BookMVA("BDT Method", "/sbnd/data/users/hlay/opt0_crumbs/training/CRUMBS_OpT0MeasuredPE/weights/CrumbsTMVAClassification_BDT.weights.xml");

  int N_all_slices = all_slices->GetEntries();

  TH1F *hNuSliceOptCrumbsScore = new TH1F("hNuSliceOptCrumbsScore",";CRUMBS Score;Slices",      50, -.7, .6);
  TH1F *hNotNuSliceOptCrumbsScore = new TH1F("hNotNuSliceOptCrumbsScore",";CRUMBS Score;Slices",50, -.7, .6);

  TH1F *hNuSliceOptScoreCrumbsScore = new TH1F("hNuSliceOptScoreCrumbsScore",";CRUMBS Score;Slices",      50, -.7, .6);
  TH1F *hNotNuSliceOptScoreCrumbsScore = new TH1F("hNotNuSliceOptScoreCrumbsScore",";CRUMBS Score;Slices",50, -.7, .6);

  TH1F *hNuSliceOptMeasCrumbsScore = new TH1F("hNuSliceOptMeasCrumbsScore",";CRUMBS Score;Slices",      50, -.7, .6);
  TH1F *hNotNuSliceOptMeasCrumbsScore = new TH1F("hNotNuSliceOptMeasCrumbsScore",";CRUMBS Score;Slices",50, -.7, .6);

  for(int all_i = 0; all_i < N_all_slices; ++all_i){
    all_slices->GetEntry(all_i);
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
    optReader_pds_OpT0Time = pds_OpT0Time;
    optReader_pds_OpT0MeasuredPE = pds_OpT0MeasuredPE;
    optReader_pds_OpT0HypothesisPE = pds_OpT0HypothesisPE;
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

  for(int all_i = 0; all_i < N_all_slices; ++all_i){
    all_slices->GetEntry(all_i);
    optScoreReader_tpc_CRFracHitsInLongestTrack = tpc_CRFracHitsInLongestTrack;
    optScoreReader_tpc_CRLongestTrackDeflection = tpc_CRLongestTrackDeflection;
    optScoreReader_tpc_CRLongestTrackDirY = tpc_CRLongestTrackDirY;
    optScoreReader_tpc_CRNHitsMax = tpc_CRNHitsMax;
    optScoreReader_tpc_NuEigenRatioInSphere = tpc_NuEigenRatioInSphere;
    optScoreReader_tpc_NuNFinalStatePfos = tpc_NuNFinalStatePfos;
    optScoreReader_tpc_NuNHitsTotal = tpc_NuNHitsTotal;
    optScoreReader_tpc_NuNSpacePointsInSphere = tpc_NuNSpacePointsInSphere;
    optScoreReader_tpc_NuVertexY = tpc_NuVertexY; 
    optScoreReader_tpc_NuWeightedDirZ = tpc_NuWeightedDirZ;
    optScoreReader_tpc_StoppingChi2Pol0 = tpc_StoppingChi2Pol0;
    optScoreReader_tpc_StoppingChi2Exp = tpc_StoppingChi2Exp;
    optScoreReader_tpc_StoppingChi2Ratio = tpc_StoppingChi2Ratio;
    optScoreReader_tpc_StoppingChi2CosmicPol0 = tpc_StoppingChi2CosmicPol0;
    optScoreReader_tpc_StoppingChi2CosmicExp = tpc_StoppingChi2CosmicExp;
    optScoreReader_tpc_StoppingChi2CosmicRatio = tpc_StoppingChi2CosmicRatio;
    optScoreReader_pds_OpT0Score = pds_OpT0Score;
    optScoreReader_crt_TrackScore = crt_TrackScore;
    optScoreReader_crt_HitScore = crt_HitScore;
    optScoreReader_crt_TrackTime = crt_TrackTime;
    optScoreReader_crt_HitTime = crt_HitTime;
    optScoreReader_crt_nTrackMatches = crt_nTrackMatches;
    optScoreReader_crt_nHitMatches = crt_nHitMatches;

    float bdtscore = optScoreReader->EvaluateMVA("BDT Method");

    if(*matchedType == "Nu" && matchedPurity > 0.8 && matchedCompleteness > 0.8)
      hNuSliceOptScoreCrumbsScore->Fill(bdtscore);
    else if(*matchedType != "Nu")
      hNotNuSliceOptScoreCrumbsScore->Fill(bdtscore); 
  }

  for(int all_i = 0; all_i < N_all_slices; ++all_i){
    all_slices->GetEntry(all_i);
    optMeasReader_tpc_CRFracHitsInLongestTrack = tpc_CRFracHitsInLongestTrack;
    optMeasReader_tpc_CRLongestTrackDeflection = tpc_CRLongestTrackDeflection;
    optMeasReader_tpc_CRLongestTrackDirY = tpc_CRLongestTrackDirY;
    optMeasReader_tpc_CRNHitsMax = tpc_CRNHitsMax;
    optMeasReader_tpc_NuEigenRatioInSphere = tpc_NuEigenRatioInSphere;
    optMeasReader_tpc_NuNFinalStatePfos = tpc_NuNFinalStatePfos;
    optMeasReader_tpc_NuNHitsTotal = tpc_NuNHitsTotal;
    optMeasReader_tpc_NuNSpacePointsInSphere = tpc_NuNSpacePointsInSphere;
    optMeasReader_tpc_NuVertexY = tpc_NuVertexY; 
    optMeasReader_tpc_NuWeightedDirZ = tpc_NuWeightedDirZ;
    optMeasReader_tpc_StoppingChi2Pol0 = tpc_StoppingChi2Pol0;
    optMeasReader_tpc_StoppingChi2Exp = tpc_StoppingChi2Exp;
    optMeasReader_tpc_StoppingChi2Ratio = tpc_StoppingChi2Ratio;
    optMeasReader_tpc_StoppingChi2CosmicPol0 = tpc_StoppingChi2CosmicPol0;
    optMeasReader_tpc_StoppingChi2CosmicExp = tpc_StoppingChi2CosmicExp;
    optMeasReader_tpc_StoppingChi2CosmicRatio = tpc_StoppingChi2CosmicRatio;
    optMeasReader_pds_OpT0Score = pds_OpT0Score;
    optMeasReader_pds_OpT0MeasuredPE = pds_OpT0MeasuredPE;
    optMeasReader_crt_TrackScore = crt_TrackScore;
    optMeasReader_crt_HitScore = crt_HitScore;
    optMeasReader_crt_TrackTime = crt_TrackTime;
    optMeasReader_crt_HitTime = crt_HitTime;
    optMeasReader_crt_nTrackMatches = crt_nTrackMatches;
    optMeasReader_crt_nHitMatches = crt_nHitMatches;

    float bdtscore = optMeasReader->EvaluateMVA("BDT Method");

    if(*matchedType == "Nu" && matchedPurity > 0.8 && matchedCompleteness > 0.8)
      hNuSliceOptMeasCrumbsScore->Fill(bdtscore);
    else if(*matchedType != "Nu")
      hNotNuSliceOptMeasCrumbsScore->Fill(bdtscore); 
  }

  float optcrumbsEff[50], optcrumbsRej[50];
  int optcrumbsSignalSum = 0, optcrumbsBackSum = 0;
  const int optcrumbsSignalTotal = hNuSliceOptCrumbsScore->GetEntries();
  const int optcrumbsBackTotal = hNotNuSliceOptCrumbsScore->GetEntries();
  for(int i = 1; i < 51; ++i)
    {
      optcrumbsSignalSum += hNuSliceOptCrumbsScore->GetBinContent(i);
      optcrumbsBackSum += hNotNuSliceOptCrumbsScore->GetBinContent(i);

      optcrumbsEff[i-1] = (float) (optcrumbsSignalTotal - optcrumbsSignalSum) / (float) optcrumbsSignalTotal;
      optcrumbsRej[i-1] = (float) optcrumbsBackSum / (float) optcrumbsBackTotal;
    }


  float optScorecrumbsEff[50], optScorecrumbsRej[50];
  int optScorecrumbsSignalSum = 0, optScorecrumbsBackSum = 0;
  const int optScorecrumbsSignalTotal = hNuSliceOptScoreCrumbsScore->GetEntries();
  const int optScorecrumbsBackTotal = hNotNuSliceOptScoreCrumbsScore->GetEntries();
  for(int i = 1; i < 51; ++i)
    {
      optScorecrumbsSignalSum += hNuSliceOptScoreCrumbsScore->GetBinContent(i);
      optScorecrumbsBackSum += hNotNuSliceOptScoreCrumbsScore->GetBinContent(i);

      optScorecrumbsEff[i-1] = (float) (optScorecrumbsSignalTotal - optScorecrumbsSignalSum) / (float) optScorecrumbsSignalTotal;
      optScorecrumbsRej[i-1] = (float) optScorecrumbsBackSum / (float) optScorecrumbsBackTotal;
    }


  float optMeascrumbsEff[50], optMeascrumbsRej[50];
  int optMeascrumbsSignalSum = 0, optMeascrumbsBackSum = 0;
  const int optMeascrumbsSignalTotal = hNuSliceOptMeasCrumbsScore->GetEntries();
  const int optMeascrumbsBackTotal = hNotNuSliceOptMeasCrumbsScore->GetEntries();
  for(int i = 1; i < 51; ++i)
    {
      optMeascrumbsSignalSum += hNuSliceOptMeasCrumbsScore->GetBinContent(i);
      optMeascrumbsBackSum += hNotNuSliceOptMeasCrumbsScore->GetBinContent(i);

      optMeascrumbsEff[i-1] = (float) (optMeascrumbsSignalTotal - optMeascrumbsSignalSum) / (float) optMeascrumbsSignalTotal;
      optMeascrumbsRej[i-1] = (float) optMeascrumbsBackSum / (float) optMeascrumbsBackTotal;
    }


  TGraph *optcrumbsROC = new TGraph(50, optcrumbsRej, optcrumbsEff);
  TGraph *optScorecrumbsROC = new TGraph(50, optScorecrumbsRej, optScorecrumbsEff);
  TGraph *optMeascrumbsROC = new TGraph(50, optMeascrumbsRej, optMeascrumbsEff);


  TLegend *leg = new TLegend(0.6,0.75,0.9,0.9);
  leg->AddEntry(hNuSliceOptCrumbsScore, "OpT0 (All Vars) #nu Slices", "l");
  leg->AddEntry(hNotNuSliceOptCrumbsScore, "OpT0 (All Vars) Other Slices", "l");
  leg->AddEntry(hNuSliceOptScoreCrumbsScore, "OpT0 (Score Only) #nu Slices", "l");
  leg->AddEntry(hNotNuSliceOptScoreCrumbsScore, "OpT0 (Score Only) Other Slices", "l");
  leg->AddEntry(hNuSliceOptMeasCrumbsScore, "OpT0 (Score + Meas PE) #nu Slices", "l");
  leg->AddEntry(hNotNuSliceOptMeasCrumbsScore, "OpT0 (Score + Meas PE) Other Slices", "l");

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

  hNuSliceOptScoreCrumbsScore->SetLineColor(kCyan+2);
  hNuSliceOptScoreCrumbsScore->SetLineWidth(3);
  hNotNuSliceOptScoreCrumbsScore->SetLineColor(kPink-5);
  hNotNuSliceOptScoreCrumbsScore->SetLineWidth(3);
  hNotNuSliceOptScoreCrumbsScore->Draw("same");
  hNuSliceOptScoreCrumbsScore->Draw("same");

  hNuSliceOptMeasCrumbsScore->SetLineColor(kGreen+4);
  hNuSliceOptMeasCrumbsScore->SetLineWidth(3);
  hNotNuSliceOptMeasCrumbsScore->SetLineColor(kOrange+5);
  hNotNuSliceOptMeasCrumbsScore->SetLineWidth(3);
  hNotNuSliceOptMeasCrumbsScore->Draw("same");
  hNuSliceOptMeasCrumbsScore->Draw("same");

  leg->Draw();

  if(save)
    {
      cCrumbsScore->SaveAs(saveDirectory + "/crumbs_score_quality.png");
      cCrumbsScore->SaveAs(saveDirectory + "/crumbs_score_quality.pdf");
    }


  TCanvas *cROC_zoomed = new TCanvas("cROC_zoomed","cROC_zoomed");
  cROC_zoomed->SetGrid();
  cROC_zoomed->cd();

  gStyle->SetLabelSize(0.03,"xy");

  TMultiGraph *multi_zoomed = new TMultiGraph();

  optcrumbsROC->SetLineColor(kAzure-5);
  optcrumbsROC->SetLineWidth(3);
  optcrumbsROC->SetMarkerColor(kAzure-5);
  optcrumbsROC->SetMarkerStyle(3);
  multi_zoomed->Add(optcrumbsROC);


  optScorecrumbsROC->SetLineColor(kPink-5);
  optScorecrumbsROC->SetLineWidth(3);
  optScorecrumbsROC->SetMarkerColor(kPink-5);
  optScorecrumbsROC->SetMarkerStyle(34);
  multi_zoomed->Add(optScorecrumbsROC);

  optMeascrumbsROC->SetLineColor(kOrange+5);
  optMeascrumbsROC->SetLineWidth(3);
  optMeascrumbsROC->SetMarkerColor(kOrange+5);
  optMeascrumbsROC->SetMarkerStyle(21);
  multi_zoomed->Add(optMeascrumbsROC);

  multi_zoomed->SetTitle(";Signal Efficiency;Background Rejection");
  multi_zoomed->Draw("ACP");
  multi_zoomed->GetXaxis()->SetLimits(0.9,1.);
  multi_zoomed->SetMinimum(0.9);
  multi_zoomed->SetMaximum(1.);

  TLegend *rocZoomLeg = new TLegend(.10,.25,.45,.45);
  rocZoomLeg->AddEntry(optcrumbsROC,"OpT0 (All Vars) CRUMBS Score","lp");
  rocZoomLeg->AddEntry(optScorecrumbsROC,"OpT0 (Score Only) CRUMBS Score","lp");
  rocZoomLeg->AddEntry(optMeascrumbsROC,"OpT0 (Score + Meas PE) CRUMBS Score","lp");
  rocZoomLeg->Draw();

  if(save)
    {
      cROC_zoomed->SaveAs(saveDirectory + "/roc_quality_zoomed.png");
      cROC_zoomed->SaveAs(saveDirectory + "/roc_quality_zoomed.pdf");
    }
}
