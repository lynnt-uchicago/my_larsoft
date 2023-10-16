// edited MARCH 20 2023 

int anodeId(int channelId){
  return (int)(channelId / 5632);
}

double median(TH1F* h1){
  if(h1->GetEntries()==0) return -1;
  double mean = h1->GetSum()/h1->GetNbinsX();
  double rms = 0;
  int nbin = h1->GetNbinsX();
  for (int j=0;j!=nbin;j++){
    rms += pow(h1->GetBinContent(j+1)-mean,2);
  }
  rms = sqrt(rms/h1->GetNbinsX());

  TH1F h2("h2_tmp","h2_tmp",100,mean-10*rms,mean+10*rms);
  for (int j=0;j!=nbin;j++){
    if (fabs(h1->GetBinContent(j+1)-mean)<6*rms)
      h2.Fill(h1->GetBinContent(j+1));
  }

  if(h2.GetEntries()==0) return -1;

  double par[3];
  double xq = 0.5;
  h2.GetQuantiles(1,&par[1],&xq);
  xq = 0.5 + 0.34;
  h2.GetQuantiles(1,&par[0],&xq);
  xq = 0.5 - 0.34;
  h2.GetQuantiles(1,&par[2],&xq);
  mean = par[1];
  // rms = std::sqrt((pow(par[0]-par[1],2)+pow(par[2]-par[1],2))/2.);
  return mean;
}

// double median(TH1I* h1){
//   if(h1->GetEntries()==0) return -1;
//   double mean = h1->GetSum()/h1->GetNbinsX();
//   double rms = 0;
//   int nbin = h1->GetNbinsX();
//   for (int j=0;j!=nbin;j++){
//     rms += pow(h1->GetBinContent(j+1)-mean,2);
//   }
//   rms = sqrt(rms/h1->GetNbinsX());

//   TH1F h2("h2_tmp","h2_tmp",100,mean-10*rms,mean+10*rms);
//   for (int j=0;j!=nbin;j++){
//     if (fabs(h1->GetBinContent(j+1)-mean)<6*rms)
//       h2.Fill(h1->GetBinContent(j+1));
//   }

//   if(h2.GetEntries()==0) return -1;

//   double par[3];
//   double xq = 0.5;
//   h2.GetQuantiles(1,&par[1],&xq);
//   xq = 0.5 + 0.34;
//   h2.GetQuantiles(1,&par[0],&xq);
//   xq = 0.5 - 0.34;
//   h2.GetQuantiles(1,&par[2],&xq);
//   mean = par[1];
//   // rms = std::sqrt((pow(par[0]-par[1],2)+pow(par[2]-par[1],2))/2.);
//   return mean;
// }


void dump_waveform_all(int chmin=0, int chmax=11264, bool trim_baseline=true, int nticks=3400, std::string filename="celltree.root"){ // [chmin, chmax]
  std::cout << "start function" << std::endl;
  auto in = TFile::Open(filename.c_str());
  auto t = dynamic_cast<TTree*>(in->Get("Event/Sim"));

  auto ofile = new TFile("waveform.root","recreate");
  auto otree = new TTree("wvfm_tree","wvfm_tree");


  TDirectory *cd_wvfms = ofile->mkdir("wvfm");
  cd_wvfms->cd();
  TDirectory *cd_raw   = cd_wvfms->mkdir("raw_wvfm");
  TDirectory *cd_decon = cd_wvfms->mkdir("decon_wvfm");
  TDirectory *cd_sim   = cd_wvfms->mkdir("sim_wvfm");

  int run;
  int sub;
  int evt;
  t->SetBranchAddress("runNo",&run);
  t->SetBranchAddress("subRunNo",&sub);
  t->SetBranchAddress("eventNo",&evt);

  int raw_nChannel;
  std::vector<int>* raw_channelId = 0;
  TClonesArray* raw_wf = 0;
  t->SetBranchAddress("raw_nChannel", &raw_nChannel);
  t->SetBranchAddress("raw_channelId", &raw_channelId);
  t->SetBranchAddress("raw_wf", &raw_wf);
 
  int calib_nChannel;
  std::vector<int>* calib_channelId = 0;
  TClonesArray* calib_wf = 0;
  t->SetBranchAddress("calib_nChannel", &calib_nChannel);
  t->SetBranchAddress("calib_channelId", &calib_channelId);
  t->SetBranchAddress("calib_wf", &calib_wf);

  int simide_size;
  std::vector<int>* simide_channelIdY = 0;
  std::vector<unsigned short>* simide_tdc = 0;
  std::vector<float>* simide_numElectrons = 0;
  t->SetBranchAddress("simide_size", &simide_size);
  t->SetBranchAddress("simide_channelIdY", &simide_channelIdY);
  t->SetBranchAddress("simide_tdc", &simide_tdc);
  t->SetBranchAddress("simide_numElectrons", &simide_numElectrons);

  std::vector<float> decon_sum;
  std::vector<float> sim_sum;

  otree->Branch("run",&run,"run/I");
  otree->Branch("sub",&sub,"sub/I");
  otree->Branch("evt",&evt,"evt/I");
  otree->Branch("decon_sum",&decon_sum);
  otree->Branch("sim_sum",&sim_sum);


  int nentries = t->GetEntries(); 
  std::cout << "Loading " << nentries << " entries..." << std::endl;

  for (int entry = 0; entry < nentries; entry++){
    t->GetEntry(entry);
    std::cout << "run: " << run << ", subrun: " << sub << ", event: " << evt << std::endl;
    if (entry%50==0) std::cout << "Entry " << entry << " out of " << nentries << std::endl;
    decon_sum = std::vector<float>(chmax,0.0);
    sim_sum = std::vector<float>(chmax,0.0);

    cd_raw->cd();

    std::stringstream h1_name; 
    h1_name.str(std::string());
    h1_name << "run_" << run
            << "_sub_" << sub 
            << "_evt_" << evt
            << "_raw"; 
           
    auto h1 = new TH2F(h1_name.str().c_str(),"raw", chmax-chmin+1, chmin-0.5, chmax+0.5, nticks, 0,nticks);
    h1->GetXaxis()->SetTitle("chID");
    h1->GetYaxis()->SetTitle("Ticks");

    // raw_wf
    for (int ich=0; ich<raw_nChannel; ich++){
      int channelId = raw_channelId->at(ich);
      // std::cout << "channelId: " << channelId << " anode: " << anodeId(channelId) << std::endl;
      if (channelId>=chmin and channelId<=chmax) {
        auto hraw  = (TH1F*)raw_wf->At(ich);
              if (nticks!= hraw->GetNbinsX()) std::cout << "Warning: nticks != " << nticks << std::endl;
              if (channelId==0) {
                  cout << "channel: " << channelId << " , RawDigit size: " << hraw->GetNbinsX() << endl;
                  cout << "bin 0: " << hraw->GetBinContent(0) << endl;
              }
              int baseline = 0;
              if (trim_baseline) baseline = median(hraw);
        // std::cout << "baseline: " << baseline << std::endl;
        for (int it=0; it<nticks; it++){
        h1->SetBinContent(channelId-chmin+1, it+1, hraw->GetBinContent(it+1) - int(baseline)); 
        }
      }
    }
    h1->Write();
    cd_wvfms->cd();

    cd_decon->cd();

    std::stringstream h2_name; 
    h2_name.str(std::string());
    h2_name << "run_"  << run
            << "_sub_" << sub 
            << "_evt_" << evt
            << "_decon"; 

    auto h2 = new TH2F(h2_name.str().c_str(),"decon", chmax-chmin+1, chmin-0.5, chmax+0.5, nticks, 0,nticks);
    h1->GetXaxis()->SetTitle("chID");
    h1->GetYaxis()->SetTitle("Ticks");

        // calib_wf
    for (int ich=0; ich<calib_nChannel; ich++){
      int channelId = calib_channelId->at(ich);
      // std::cout << "channelId: " << channelId << " anode: " << anodeId(channelId) << std::endl;
      if (channelId>=chmin and channelId<=chmax) {
        auto hcalib  = (TH1F*)calib_wf->At(ich);
        if (nticks!= hcalib->GetNbinsX()) std::cout << "Warning: nticks != " << nticks << std::endl;
        if (channelId==0) {
                    cout << "channel: " << channelId << " , recob::Wire size: " << hcalib->GetNbinsX() << endl;
                    cout << "bin 0: " << hcalib->GetBinContent(0) << endl;
                }
        for (int it=0; it<nticks; it++){
          h2->SetBinContent(channelId-chmin+1, it+1, hcalib->GetBinContent(it+1)); 
          decon_sum.at(channelId) += hcalib->GetBinContent(it+1);
        }
      }
    }
    h2->Write();
    cd_wvfms->cd();

    cd_sim->cd();

    std::stringstream h3_name; 
    h3_name.str(std::string());
    h3_name << "run_" << run
            << "_sub_" << sub 
            << "_evt_" << evt
            << "_sim"; 
    auto h3 = new TH2F(h3_name.str().c_str(),"sim", chmax-chmin+1, chmin-0.5, chmax+0.5, nticks, 0,nticks);
    // simide
    for (int idx=0; idx<simide_size; idx++){
      int channelId = simide_channelIdY->at(idx);
      unsigned int tdc = simide_tdc->at(idx);
            tdc -= 3400 - 400; // -1.7 ms g4 ref time; -0.2 ms TPC readout t0
      float nele = simide_numElectrons->at(idx);
      if (channelId>=chmin and channelId<=chmax) {
                  //cout << "channel: " << channelId << " tdc: " << tdc << " nele: " << nele << endl;
        int tbin = (int)(tdc);
        h3->SetBinContent(channelId-chmin+1, tbin, nele + h3->GetBinContent(channelId-chmin+1, tbin));
        sim_sum.at(channelId) += nele;
      }
    }
    h3->Write();
    cd_wvfms->cd();

    otree->Fill();
    h1->Reset();
    h2->Reset();
    h3->Reset();

    delete h1;
    delete h2;
    delete h3;
  }
  std::cout << "end function" << std::endl;
  ofile->Write();
  ofile->Close();
  in->Close();

}
