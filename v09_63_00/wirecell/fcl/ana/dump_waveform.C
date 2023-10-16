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

void dump_waveform(int chmin=0, int chmax=11264, bool trim_baseline=true, int nticks=3400, std::string filename="celltree.root"){ // [chmin, chmax]
  auto in = TFile::Open(filename.c_str());
  auto t = dynamic_cast<TTree*>(in->Get("Event/Sim"));

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


  t->GetEntry(0);
  std::cout << calib_nChannel << std::endl;
  // std::cout << simide_size << std::endl; 
  // std::cout << simide_channelIdY->size() << std::endl;

  auto h1 = new TH2F("h1","raw", chmax-chmin+1, chmin-0.5, chmax+0.5, nticks, 0,nticks);
  h1->GetXaxis()->SetTitle("chID");
  h1->GetYaxis()->SetTitle("Ticks");

  auto h2 = new TH2F("h2","decon", chmax-chmin+1, chmin-0.5, chmax+0.5, nticks, 0,nticks);
  h2->GetXaxis()->SetTitle("chID");
  h2->GetYaxis()->SetTitle("Ticks");



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

                double baseline = 0;
                if (trim_baseline) baseline = median(hraw);
        	std::cout << "baseline: " << baseline << std::endl;
  		for (int it=0; it<nticks; it++){
  			h1->SetBinContent(channelId-chmin+1, it+1, hraw->GetBinContent(it+1) - baseline); 
  		}

        }
  }


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
  		}

        }
  }

   auto h3 = new TH2F("h3","simChannel", chmax-chmin+1, chmin-0.5, chmax+0.5, nticks, 0,nticks);

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
   	}
   }

  auto ofile = new TFile("waveform.root","recreate");
  h1->Write();
   h2->Write();
   h3->Write();
  ofile->Close();

}
