#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ValidatedPluginMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "vdt/vdtMath.h"
#include "TRandom3.h"

#include "RecoVertex/PrimaryVertexProducer/interface/VertexTimeAlgorithmFromTracksPID.h"

#ifdef PVTX_DEBUG
#define LOG edm::LogPrint("VertexTimeAlgorithmFromTracksPID")
#else
#define LOG LogDebug("VertexTimeAlgorithmFromTracksPID")
#endif

using namespace std;

VertexTimeAlgorithmFromTracksPID::VertexTimeAlgorithmFromTracksPID(edm::ParameterSet const& iConfig,
                                                                   edm::ConsumesCollector& iCC)
    : VertexTimeAlgorithmBase(iConfig, iCC),
      trackMTDTimeToken_(iCC.consumes(iConfig.getParameter<edm::InputTag>("trackMTDTimeVMapTag"))),
      trackMTDTimeErrorToken_(iCC.consumes(iConfig.getParameter<edm::InputTag>("trackMTDTimeErrorVMapTag"))),
      trackMTDTimeQualityToken_(iCC.consumes(iConfig.getParameter<edm::InputTag>("trackMTDTimeQualityVMapTag"))),
      trackMTDTofPiToken_(iCC.consumes(iConfig.getParameter<edm::InputTag>("trackMTDTofPiVMapTag"))),
      trackMTDTofKToken_(iCC.consumes(iConfig.getParameter<edm::InputTag>("trackMTDTofKVMapTag"))),
      trackMTDTofPToken_(iCC.consumes(iConfig.getParameter<edm::InputTag>("trackMTDTofPVMapTag"))),
      trackMTDSigmaTofPiToken_(iCC.consumes(iConfig.getParameter<edm::InputTag>("trackMTDSigmaTofPiVMapTag"))),
      trackMTDSigmaTofKToken_(iCC.consumes(iConfig.getParameter<edm::InputTag>("trackMTDSigmaTofKVMapTag"))),
      trackMTDSigmaTofPToken_(iCC.consumes(iConfig.getParameter<edm::InputTag>("trackMTDSigmaTofPVMapTag"))),
      minTrackVtxWeight_(iConfig.getParameter<double>("minTrackVtxWeight")),
      minTrackTimeQuality_(iConfig.getParameter<double>("minTrackTimeQuality")),
      probPion_(iConfig.getParameter<double>("probPion")),
      probKaon_(iConfig.getParameter<double>("probKaon")),
      probProton_(iConfig.getParameter<double>("probProton")),
      Tstart_(iConfig.getParameter<double>("Tstart")),
      coolingFactor_(iConfig.getParameter<double>("coolingFactor")),
      populationSize_(45),
      Nm_(3){}
      

void VertexTimeAlgorithmFromTracksPID::fillPSetDescription(edm::ParameterSetDescription& iDesc) {
  VertexTimeAlgorithmBase::fillPSetDescription(iDesc);

  iDesc.add<edm::InputTag>("trackMTDTimeVMapTag", edm::InputTag("trackExtenderWithMTD:generalTracktmtd"))
      ->setComment("Input ValueMap for track time at MTD");
  iDesc.add<edm::InputTag>("trackMTDTimeErrorVMapTag", edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd"))
      ->setComment("Input ValueMap for track time uncertainty at MTD");
  iDesc.add<edm::InputTag>("trackMTDTimeQualityVMapTag", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"))
      ->setComment("Input ValueMap for track MVA quality value");
  iDesc.add<edm::InputTag>("trackMTDTofPiVMapTag", edm::InputTag("trackExtenderWithMTD:generalTrackTofPi"))
      ->setComment("Input ValueMap for track tof as pion");
  iDesc.add<edm::InputTag>("trackMTDTofKVMapTag", edm::InputTag("trackExtenderWithMTD:generalTrackTofK"))
      ->setComment("Input ValueMap for track tof as kaon");
  iDesc.add<edm::InputTag>("trackMTDTofPVMapTag", edm::InputTag("trackExtenderWithMTD:generalTrackTofP"))
      ->setComment("Input ValueMap for track tof as proton");
  iDesc.add<edm::InputTag>("trackMTDSigmaTofPiVMapTag", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofPi"))
      ->setComment("Input ValueMap for track tof uncertainty as pion");
  iDesc.add<edm::InputTag>("trackMTDSigmaTofKVMapTag", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofK"))
      ->setComment("Input ValueMap for track tof uncertainty as kaon");
  iDesc.add<edm::InputTag>("trackMTDSigmaTofPVMapTag", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofP"))
      ->setComment("Input ValueMap for track tof uncertainty as proton");

  iDesc.add<double>("minTrackVtxWeight", 0.5)->setComment("Minimum track weight");
  iDesc.add<double>("minTrackTimeQuality", 0.8)->setComment("Minimum MVA Quality selection on tracks");

  iDesc.add<double>("probPion", 0.7)->setComment("A priori probability pions");
  iDesc.add<double>("probKaon", 0.2)->setComment("A priori probability kaons");
  iDesc.add<double>("probProton", 0.1)->setComment("A priori probability protons");

  iDesc.add<double>("Tstart", 256.)->setComment("DA initial temperature T");
  iDesc.add<double>("coolingFactor", 0.5)->setComment("DA cooling factor");
}

void VertexTimeAlgorithmFromTracksPID::setEvent(edm::Event& iEvent, edm::EventSetup const&) {
  // additional collections required for vertex-time calculation
  trackMTDTimes_ = iEvent.get(trackMTDTimeToken_);
  trackMTDTimeErrors_ = iEvent.get(trackMTDTimeErrorToken_);
  trackMTDTimeQualities_ = iEvent.get(trackMTDTimeQualityToken_);
  trackMTDTofPi_ = iEvent.get(trackMTDTofPiToken_);
  trackMTDTofK_ = iEvent.get(trackMTDTofKToken_);
  trackMTDTofP_ = iEvent.get(trackMTDTofPToken_);
  trackMTDSigmaTofPi_ = iEvent.get(trackMTDSigmaTofPiToken_);
  trackMTDSigmaTofK_ = iEvent.get(trackMTDSigmaTofKToken_);
  trackMTDSigmaTofP_ = iEvent.get(trackMTDSigmaTofPToken_);
  mRan=new TRandom3();
}

bool VertexTimeAlgorithmFromTracksPID::vertexTime(float& vtxTime,
                                                  float& vtxTimeError,
                                                  const TransientVertex& vtx) const {
  if (vtx.originalTracks().empty()) {
    return false;
  }

  auto const vtxTime_init = vtxTime;
  auto const vtxTimeError_init = vtxTimeError;
  const int max_steps=1500;//todo:maybe add protection max_step<(3^Ntrk-populationSize_)

  //TRandom3 mRan(123);//fix the random seed
  
  vector<TrackInfo> v_trackInfo;
  v_trackInfo.reserve(vtx.originalTracks().size());
  double tmp_ndof=0.0;

  //Loop over the tracks associated with the vertex
  for (const auto& trk : vtx.originalTracks()) {
    auto const trkWeight = vtx.trackWeight(trk);
    if (trkWeight > minTrackVtxWeight_) {
      auto const trkTimeQuality = trackMTDTimeQualities_[trk.trackBaseRef()];

      if (trkTimeQuality >= minTrackTimeQuality_) {
        tmp_ndof+=trkWeight;
        auto const trkTime = trackMTDTimes_[trk.trackBaseRef()];
        auto const trkTimeError = trackMTDTimeErrors_[trk.trackBaseRef()];

        v_trackInfo.emplace_back();//insert a new element at the end of the vector, 
        //effectively increases the container size by one.      
        auto& trkInfo = v_trackInfo.back();//Returns a reference to the last element in the vector

        //trkInfo.trkWeight = trkWeight;
        trkInfo.trkTimeHyp[0] = trkTime - trackMTDTofPi_[trk.trackBaseRef()];
        trkInfo.trkTimeHyp[1] = trkTime - trackMTDTofK_[trk.trackBaseRef()];
        trkInfo.trkTimeHyp[2] = trkTime - trackMTDTofP_[trk.trackBaseRef()];
        
        trkInfo.trkTimeErrorHyp[0] =
            std::sqrt(trkTimeError * trkTimeError +
                      trackMTDSigmaTofPi_[trk.trackBaseRef()] * trackMTDSigmaTofPi_[trk.trackBaseRef()]);
        trkInfo.trkTimeErrorHyp[1] =
            std::sqrt(trkTimeError * trkTimeError +
                      trackMTDSigmaTofK_[trk.trackBaseRef()] * trackMTDSigmaTofK_[trk.trackBaseRef()]);
        trkInfo.trkTimeErrorHyp[2] =
            std::sqrt(trkTimeError * trkTimeError +
                      trackMTDSigmaTofP_[trk.trackBaseRef()] * trackMTDSigmaTofP_[trk.trackBaseRef()]);
        /*
        trkInfo.trkTimeHyp[0] = trkTime - trackMTDTofPi_[trk.trackBaseRef()];
        trkInfo.trkTimeHyp[1] = trkTime - trackMTDTofK_[trk.trackBaseRef()];
        trkInfo.trkTimeHyp[2] = trkTime - trackMTDTofP_[trk.trackBaseRef()];
        
        trkInfo.trkTimeErrorHyp[0] =
            std::sqrt(trkTimeError * trkTimeError +
                      trackMTDSigmaTofPi_[trk.trackBaseRef()] * trackMTDSigmaTofPi_[trk.trackBaseRef()]);
        trkInfo.trkTimeErrorHyp[1] =
            std::sqrt(trkTimeError * trkTimeError +
                      trackMTDSigmaTofK_[trk.trackBaseRef()] * trackMTDSigmaTofK_[trk.trackBaseRef()]);
        trkInfo.trkTimeErrorHyp[2] =
            std::sqrt(trkTimeError * trkTimeError +
                      trackMTDSigmaTofP_[trk.trackBaseRef()] * trackMTDSigmaTofP_[trk.trackBaseRef()]);
        */
        trkInfo.weight[0]=trkWeight / (trkInfo.trkTimeErrorHyp[0] * trkInfo.trkTimeErrorHyp[0]);
        trkInfo.weight[1]=trkWeight / (trkInfo.trkTimeErrorHyp[1] * trkInfo.trkTimeErrorHyp[1]);
        trkInfo.weight[2]=trkWeight / (trkInfo.trkTimeErrorHyp[2] * trkInfo.trkTimeErrorHyp[2]);
        /*
        LOG << "vertexTimeFromTracks:     track"
            << " pt=" << trk.track().pt() << " eta=" << trk.track().eta() << " phi=" << trk.track().phi()
            << " vtxWeight=" << trkWeight << " time=" << trkTime << " timeError=" << trkTimeError
            << " timeQuality=" << trkTimeQuality << " timeHyp[pion]=" << trkInfo.trkTimeHyp[0] << " +/- "
            << trkInfo.trkTimeErrorHyp[0] << " timeHyp[kaon]=" << trkInfo.trkTimeHyp[1] << " +/- "
            << trkInfo.trkTimeErrorHyp[1] << " timeHyp[proton]=" << trkInfo.trkTimeHyp[2] << " +/- "
            << trkInfo.trkTimeErrorHyp[2];
        */
      }
    }
  }//End looping over tracks
  
  //----------Initiate the population of solution vectors----------
  int Ntrks=v_trackInfo.size();
  
  //Add protection for population<=Nm^Ntrks, where Nm is the possible number of particle spieces.
  if(populationSize_>=pow(3,Ntrks)){
    //LOG <<"Population size larger than Nm_^(Ntrks)";
    vtxTime = vtxTime_init;
    vtxTimeError = vtxTimeError_init;
    //cout<<"vz="<<vtx.position().z()<<" Ntrk="<<Ntrks<<" ndof="<<tmp_ndof<<" vtxTime="<<vtxTime<<" vtxTimeError="<<vtxTimeError<<endl;
    return false;
  }

  vector<vector<int>> populationVector;
  for(int i=0;i<populationSize_;i++){
    populationVector.emplace_back();
  }
  vector<int> massTracker;
  int count=0;
  while(count<populationSize_){
    populationVector[count].clear();
    for(int j=0;j<Ntrks;j++){
      int tmp=mRan->Integer(3);//random integer in [0,2]
      populationVector[count].push_back(tmp);
      if(count==0){
        massTracker.push_back(0);
      }
      else if(count!=0&&massTracker[j]==0){
        if(populationVector[count][j]!=populationVector[count-1][j]){
          massTracker[j]=1;
        }
      }
    }
    //make sure all population members are unique.
    if(count!=0){
      for(int k=count-1;k>=0;k--){
        if(populationVector[count]==populationVector[k]){
          count--;
          break;
        }
      }
    }
    //make sure at least two different masses per track;
    //only need to check and update the last vector in the population.
    if(count==populationSize_-1){
      int tmp_sum=0;
      for_each(massTracker.begin(), massTracker.end(), [&](int i) {
		    tmp_sum += i;
	    });
      if(tmp_sum!=Ntrks) count--;
    }
    count++; 
  }
   
  //calculate and record the vertex time and chi2 of the population
  vector<double> population_chi2;
  vector<double> population_tv;
  vector<double> population_sigmatv;

  for(int i=0;i<populationSize_;i++){
    double sum_parent=.0;
    double sumw2sigma2_parent=.0;
    double wgt_parent=.0;
    for(int j=0;j<Ntrks;j++){
      int tmp_pid=populationVector[i][j];
      sum_parent+=v_trackInfo[j].weight[tmp_pid]*v_trackInfo[j].trkTimeHyp[tmp_pid];
      sumw2sigma2_parent+=v_trackInfo[j].weight[tmp_pid]*v_trackInfo[j].weight[tmp_pid]
      *v_trackInfo[j].trkTimeErrorHyp[tmp_pid]*v_trackInfo[j].trkTimeErrorHyp[tmp_pid];
      wgt_parent+=v_trackInfo[j].weight[tmp_pid];
    }
    double tv_parent=sum_parent/wgt_parent;
    double sigmatv_parent=sqrt(sumw2sigma2_parent)/wgt_parent;
    population_tv.push_back(tv_parent);
    population_sigmatv.push_back(sigmatv_parent);

    double chi2_parent=.0;
    for(int j=0;j<Ntrks;j++){
      int tmp_pid=populationVector[i][j];
      chi2_parent+=v_trackInfo[j].weight[tmp_pid]*(v_trackInfo[j].trkTimeHyp[tmp_pid]-tv_parent)*(v_trackInfo[j].trkTimeHyp[tmp_pid]-tv_parent);
    }
    population_chi2.push_back(chi2_parent);
  }
  
  //make a map of population index and chi2
  //from now on, if solution i in the population gets updated, 
  //populationVector[i], population_tv[i], population_chi2[i], population_sigmatv[i] all get updated.
  //Also, population index i get mapped to a new chi2 
  vector<pair<int, double>> popindex_chi2_pairs;
  for(int i=0;i<populationSize_;i++){
    popindex_chi2_pairs.push_back(pair<int, double>{i,population_chi2[i]});
    /*
    for(int j=0;j<(int)populationVector[i].size();j++){
      cout<<populationVector[i][j]<<" ";
    }
    cout<<" chi2="<<popindex_chi2_pairs[i].second<<endl;
    */
  }

  sort(popindex_chi2_pairs.begin(),popindex_chi2_pairs.end(),[](auto&a, auto &b){
      return a.second<b.second;
  });

  //For now, set the derminating criterion to the number of iterations.
  int nstep = 0;
  //===Iteration starts===
  vector<int> v_mutant; 
  int ids[3]={-1,-1,-1};
  //cout<<(popindex_chi2_pairs[7].second-popindex_chi2_pairs[0].second)/(double)Ntrks<<endl;
  //while((popindex_chi2_pairs[7].second-popindex_chi2_pairs[0].second)/(double)Ntrks>1e-2){
    /*
    if(nstep>max_steps){
      vtxTime = vtxTime_init;
      int tmp_index=popindex_chi2_pairs[0].first;
      vtxTime = population_tv[tmp_index];
      vtxTimeError = vtxTimeError_init;
      int popindex_chi2min=(*popindex_chi2_pairs.begin()).first;
      for(int i=0;i<populationSize_;i++){
        cout<<popindex_chi2_pairs[i].first<<" "<<popindex_chi2_pairs[i].second<<endl;
      }
      LOG <<"ADEGA fails to converge, Ntrks="<<Ntrks<<" chi2_min="<<population_chi2[popindex_chi2min]<<" Nsteps="<<nstep
      <<" vertex_time="<<vtxTime<<" vertex_time_error="<<vtxTimeError;
      return false;
    }
    */
    while((nstep++)<max_steps){
      //----------Offspring generation----------
      int munique=0;
      while(!munique){
        //randomly generate three distinct ids
        gen3Id(ids);
        v_mutant.clear();
        //calculate the mutant vector and get it in range
        //for(auto &element:populationVector[0]){
        for(int element=0;element<(int)populationVector[0].size();element++){
          double tmp=populationVector[ids[0]][element]+populationVector[ids[1]][element]-populationVector[ids[2]][element];
          v_mutant.push_back(tmp);
        }
      /*
      for(int j=0;j<(int)v_mutant.size();j++){
        cout<<v_mutant[j]<<" ";
      }
      cout<<endl;
      */
        for(int i=0;i<Ntrks;i++){
          v_mutant[i]=getInRange(v_mutant[i]);
        }
      /*
      for(int j=0;j<(int)v_mutant.size();j++){
        cout<<v_mutant[j]<<" ";
      }
      cout<<endl;
      */
        munique=1;
        for(int i=0;i<populationSize_;i++){
          if(v_mutant==populationVector[i]){
            munique=0;
            break;
          }
        }
      }
      /*
      for(int j=0;j<(int)v_mutant.size();j++){
        cout<<v_mutant[j]<<" ";
      }
      cout<<endl;
      */
      //----------Natural selection----------
      //calcualte vtx time and chi2 of v_mutant
      double sum_mut=.0;
      double sumw2sigma2_mut=.0;
      double wgt_mut=.0;
      for(int i=0;i<Ntrks;i++){
        int tmp_pid=v_mutant[i];
        sum_mut+=v_trackInfo[i].weight[tmp_pid]*v_trackInfo[i].trkTimeHyp[tmp_pid];
        sumw2sigma2_mut+=v_trackInfo[i].weight[tmp_pid]*v_trackInfo[i].weight[tmp_pid]*v_trackInfo[i].trkTimeErrorHyp[tmp_pid]*v_trackInfo[i].trkTimeErrorHyp[tmp_pid];
        wgt_mut+=v_trackInfo[i].weight[tmp_pid];
      }
      double tv_mut=sum_mut/wgt_mut;
      double sigmatv_mut=sqrt(sumw2sigma2_mut)/wgt_mut;
      
      double chi2_mut=.0;
      for(int i=0;i<Ntrks;i++){
        int tmp_pid=v_mutant[i];
        chi2_mut+=v_trackInfo[i].weight[tmp_pid]*(v_trackInfo[i].trkTimeHyp[tmp_pid]-tv_mut)*(v_trackInfo[i].trkTimeHyp[tmp_pid]-tv_mut);
      }
      int tmp_parent_index=ids[0];
      /* 
      for(int j=0;j<(int)v_mutant.size();j++){
        cout<<v_mutant[j]<<" ";
      }
      cout<<endl;
      cout<<"chi2_mut="<<chi2_mut<<" chi2_parent="<<population_chi2[tmp_parent_index]<<endl;
      */
      //compare chi2_mut with chi2_parent.
      if(chi2_mut<population_chi2[tmp_parent_index]){
        //cout<<"interation "<<nstep<<" population "<<tmp_parent_index<<" chi2="<<population_chi2[tmp_parent_index]<<" -> chi2="<<chi2_mut<<endl;
        populationVector[tmp_parent_index]=v_mutant;
        population_tv[tmp_parent_index]=tv_mut;
        population_chi2[tmp_parent_index]=chi2_mut;
        population_sigmatv[tmp_parent_index]=sigmatv_mut;
        vector<pair<int,double>>::iterator it;
        it=find_if(popindex_chi2_pairs.begin(), popindex_chi2_pairs.end(),[tmp_parent_index](const auto& p){return p.first==tmp_parent_index;});
        (*it).second=population_chi2[tmp_parent_index];
      }
      //else{
      //  cout<<"population "<<tmp_parent_index<<" chi2="<<population_chi2[tmp_parent_index]<<endl;
      //}

      //LOG << "vertexTimeFromTracks:   iteration=" << nstep << " chi2_mut="<<chi2_mut<<" chi2_parent="<<population_chi2[tmp_parent_index];
      /* 
      for(int i=0;i<populationSize_;i++){

        for(int j=0;j<(int)populationVector[i].size();j++){
          cout<<populationVector[i][j]<<" ";
        }
        cout<<" chi2="<<popindex_chi2_pairs[i].second<<endl;
      }
      */
      std::sort(popindex_chi2_pairs.begin(),popindex_chi2_pairs.end(),[](auto&a, auto &b){
          return a.second<b.second;
      });
      /*
      for(int i=0;i<populationSize_;i++){
        int tmp_index_ordered=popindex_chi2_pairs[i].first;
        for(int j=0;j<(int)populationVector[tmp_index_ordered].size();j++){
          cout<<populationVector[tmp_index_ordered][j]<<" ";
        }
        cout<<" chi2="<<popindex_chi2_pairs[i].second<<endl;
      }
      */
      LOG << "vertexTimeFromTracks:   iteration=" << nstep <<" chi2_min/Ntrk="<<popindex_chi2_pairs[0].second/(double)Ntrks;
      //nstep++;
    }//=== Iteration ends===
     
    int popindex_chi2min=(*popindex_chi2_pairs.begin()).first;
    vtxTime=population_tv[popindex_chi2min];
    vtxTimeError=population_sigmatv[popindex_chi2min];
    /*
    for(int i=0;i<populationSize_;i++){
      cout<<popindex_chi2_pairs[i].first<<" "<<popindex_chi2_pairs[i].second<<endl;
    }
    */
    /*
    LOG <<"Ntrks="<<Ntrks<<" chi2_min="<<population_chi2[popindex_chi2min]<<" Nsteps="<<nstep
    <<" vertex_time="<<vtxTime<<" vertex_time_error="<<vtxTimeError;
    */
    //cout<<Ntrks<<" "<<population_chi2[popindex_chi2min]<<" "<<vtx.position().z()<<endl;
    //cout<<tmp_ndof<<" "<<population_chi2[popindex_chi2min]<<" "<<vtx.position().z()<<endl;
    //cout<<"vz="<<vtx.position().z()<<" Ntrk="<<Ntrks<<" ndof="<<tmp_ndof<<" chi2="<<population_chi2[popindex_chi2min]<<endl;
    //cout<<vtx.position().z()<<" "<<Ntrks<<" "<<tmp_ndof<<" "<<population_chi2[popindex_chi2min]<<endl;
    //ref e.g:https://github.com/cms-sw/cmssw/blob/master/Configuration/Skimming/src/LeptonSkimming.cc

    return true;  
  //}

}

void VertexTimeAlgorithmFromTracksPID::gen3Id(int (&ids)[3]) const {
  //TRandom3 mRan_tmp(123);//cannot fix random seed!!!
  //randomly choose a parent vector
  int parent_id=mRan->Integer(populationSize_);
  ids[0]=parent_id;
  //randomly choose two distinct vectors that are different from the parent vector.
  int delta_v_id1=mRan->Integer(populationSize_);
  while(delta_v_id1==parent_id){
    delta_v_id1=mRan->Integer(populationSize_);
  }
  ids[1]=delta_v_id1;
  int delta_v_id2=mRan->Integer(populationSize_);
  while(delta_v_id2==delta_v_id1||delta_v_id2==parent_id){
    delta_v_id2=mRan->Integer(populationSize_);
  }
  ids[2]=delta_v_id2;
}

int VertexTimeAlgorithmFromTracksPID::getInRange(int pid) const {
  if(pid>Nm_-1) return Nm_-1;
  else if(pid<0) return 0;
  else return pid;
}
