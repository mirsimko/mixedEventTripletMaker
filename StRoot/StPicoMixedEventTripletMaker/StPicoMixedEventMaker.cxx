#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TList.h"
#include "TString.h" // needed for the Form(...)

#include "StarClassLibrary/StThreeVectorF.hh"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StPicoMixedEventMaker.h"
#include "StPicoEventMixer.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StPicoHFMaker/StHFCuts.h"

#include <vector>
#include <string>

ClassImp(StPicoMixedEventMaker)


// _________________________________________________________
StPicoMixedEventMaker::StPicoMixedEventMaker(char const* name, StPicoDstMaker* picoMaker, StRefMultCorr* grefmultCorrUtil, StHFCuts* hfCuts,
        char const* outputBaseFileName,  char const* inputHFListHFtree = "") :
    StMaker(name), 
    mPicoDst(NULL), 
    mPicoDstMaker(picoMaker),  
    mPicoEvent(NULL),
    mGRefMultCorrUtil(grefmultCorrUtil),
    mHFCuts(hfCuts),
    mOuputFileBaseName(outputBaseFileName), 
    mInputFileName(inputHFListHFtree),
    mEventCounter(0), 
    mBufferSize(defaultBufferSize), 
    mSETuple(NULL), 
    mMETuple(NULL), 
    mOutputFileTree(NULL),
    mSingePartHists(NULL)
{

  TH1::AddDirectory(false);
  // -- create OutputTree
  mOutputFileTree = new TFile(Form("%s.picoMEtree.root", mOuputFileBaseName.Data()), "RECREATE");
  mOutputFileTree->SetCompressionLevel(1);
  mOutputFileTree->cd();

  const string varList = "p1pt:p2pt:p3pt:"
			 "charges:"
			 "m:pt:eta:phi:"
			 "cosPntAngle:dLength:"
			 "DCAtoPV:"
			 "p1Dca:p2Dca:p3Dca:"
			 "dcaDaughters12:dcaDaughters23:dcaDaughters31:"
			 "KNSigma:pNSigma:piNSigma:"
			 "KTOFinvBetaDiff:pTOFinvBetaDiff:piTOFinvBetaDiff:"
			 "KEta:pEta:piEta:"
			 "KPhi:pPhi:piPhi:"
			 "maxVertexDist:"
			 "centrality:centralityCorrection";

  mSETuple = new TNtuple("sameEvent", "SameEvent", varList.data() );
  mMETuple = new TNtuple("mixedEvent", "MixedEvent", varList.data() );

  mSingePartHists = new TList();

  mSingePartHists->SetOwner(true);
  mSingePartHists->SetName("HFSinglePartHists");

  // create single particle hists
  const std::string evtNames[2] = {"SE", "ME"};
  const std::string partNames[3] = {"pi", "p", "K"};

  for (int i = 0; i < 2; ++i)
  {
    // mSingePartHists->Add(new TH1D(Form("centrality%s", evtNames[i].data()),Form("centrality%s", evtNames[i].data()), 10, -1.5, 8.5));
    // mSingePartHists->Add(new TH1D(Form("centralityCorrection%s", evtNames[i].data()),Form("centrality corrected %s", evtNames[i].data()), 10, -1.5, 8.5));
    // mSingePartHists->Add(new TH1D(Form("refMult%s", evtNames[i].data()), Form("corrected refferernce multiplicity %s", evtNames[i].data()), 100, 0, 800));

    for (int iPart = 0; iPart < 3; ++iPart)
    {
      // eta phi
      mSingePartHists->Add(new TH2D(Form("%sEtaPhi%s",partNames[ iPart ].data(), evtNames[i].data()),
				    Form("%s Eta phi distribution %s",partNames[ iPart ].data(), evtNames[i].data()), 
				    100, -TMath::Pi(), TMath::Pi(), 100, -1.1, 1.1));

      // phi vs pT
      mSingePartHists->Add(new TH2D(Form("%sPhiPt%s",partNames[ iPart ].data(), evtNames[i].data()), 
				    Form("%s phi vs pT %s",partNames[ iPart ].data(), evtNames[i].data()), 
				    100, 0, 15, 100, -TMath::Pi(), TMath::Pi()));

      // DCA
      mSingePartHists->Add(new TH1D(Form("%sDCA%s",partNames[ iPart ].data(), evtNames[i].data()),
				    Form("%s DCA %s",partNames[ iPart ].data(), evtNames[i].data()), 
				    200, 0, 0.02));
      // nTracks
      mSingePartHists->Add(new TH1D(Form("%stracks%s",partNames[ iPart ].data(), evtNames[i].data()),
				    Form("Number of %s tracks %s",partNames[ iPart ].data(), evtNames[i].data()), 
				    100, -0.5, 99.5));
    }
  }

  // loop over all histograms to set Sumw2
  TH1* hist = static_cast<TH1*>(mSingePartHists->First());
  hist->Sumw2();
  while(hist != static_cast<TH1*>(mSingePartHists->Last()))
  {
    hist = static_cast<TH1*>(mSingePartHists->After(hist));
    hist->Sumw2();
  }

}

// _________________________________________________________
StPicoMixedEventMaker::~StPicoMixedEventMaker() {

  for(int iVz =0 ; iVz < 10 ; ++iVz){
    for(int iCentrality = 0 ; iCentrality < 9 ; ++iCentrality){
      delete mPicoEventMixer[iVz][iCentrality];
    }
  }
  mOutputFileTree->Close();
}
// _________________________________________________________
bool StPicoMixedEventMaker::loadEventPlaneCorr(Int_t const run) {
    //needs to implement, will currently break maker
    return false;
}
// _________________________________________________________
Int_t StPicoMixedEventMaker::Init() {
    mOutputFileTree->cd();
    for(int iVz =0 ; iVz < 10 ; ++iVz){
      for(int iCentrality = 0 ; iCentrality < 9 ; ++iCentrality){
	mPicoEventMixer[iVz][iCentrality] = new StPicoEventMixer(Form("Cent_%i_Vz_%i",iCentrality,iVz));
	mPicoEventMixer[iVz][iCentrality]->setEventBuffer(mBufferSize);
	mPicoEventMixer[iVz][iCentrality]->setHFCuts(mHFCuts);
	mPicoEventMixer[iVz][iCentrality]->setSameEvtNtuple(mSETuple);
	mPicoEventMixer[iVz][iCentrality]->setMixedEvtNtuple(mMETuple);
	mPicoEventMixer[iVz][iCentrality]->setSinglePartHistsList(mSingePartHists);
	mPicoEventMixer[iVz][iCentrality]->setFillSinglePartHists(fillSingleTrackHistos);
      }
    }
    // if(!LoadEventPlaneCorr(mRunId)){
    // LOG_WARN << "Event plane calculations unavalable! Skipping"<<endm;
    // return kStOk;
    // }

    // -- reset event to be in a defined state
    //resetEvent();

    return kStOK;
}

// _________________________________________________________
Int_t StPicoMixedEventMaker::Finish() {
    mOutputFileTree->cd();
    for(int iVz =0 ; iVz < 10 ; ++iVz){
      for(int iCentrality = 0 ; iCentrality < 9 ; ++iCentrality){
	mPicoEventMixer[iVz][iCentrality]->finish();
	//delete mPicoEventMixer[iVz][iCentrality];
      }
    }
    mSingePartHists->Write();
    return kStOK;
}
// _________________________________________________________
void StPicoMixedEventMaker::Clear(Option_t* opt) {
}
// _________________________________________________________
Int_t StPicoMixedEventMaker::Make() {

    if(!mPicoDstMaker) {
        LOG_WARN << "No PicoDstMaker! Skipping! "<<endm;
        return kStWarn;
    }

    StPicoDst const* picoDst = mPicoDstMaker->picoDst();
    if (!picoDst) {
        LOG_WARN << "No picoDst ! Skipping! "<<endm;
        return kStWarn;
    }
    // - GRefMultiplicty
    if(!mGRefMultCorrUtil) {
        LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
        return kStWarn;
    }
    if (!mHFCuts->isGoodEvent(picoDst))
      return kStOk;
    StThreeVectorF const pVtx = picoDst->event()->primaryVertex();
    if( fabs(pVtx.z()) >=6.0 )
      return kStOk;
    mGRefMultCorrUtil->init(picoDst->event()->runId());
    mGRefMultCorrUtil->initEvent(picoDst->event()->grefMult(),pVtx.z(),picoDst->event()->ZDCx()) ;
    int const centrality  = mGRefMultCorrUtil->getCentralityBin9();
    if(centrality < 0 || centrality >8 ) return kStOk;
    int const vz_bin = (int)((6 +pVtx.z())/1.2) ;
    //     Bin       Centrality (16)   Centrality (9)
    //     -1           80-100%           80-100% // this one should be rejected in your centrality related analysis
    //     0            75-80%            70-80%
    //     1            70-75%            60-70%
    //     2            65-70%            50-60%
    //     3            60-65%            40-50%
    //     4            55-60%            30-40%
    //     5            50-55%            20-30%
    //     6            45-50%            10-20%
    //     7            40-45%             5-10%
    //     8            35-40%             0- 5%

    if( mPicoEventMixer[vz_bin][centrality] -> addPicoEvent(picoDst, mGRefMultCorrUtil->getWeight()) ==  true )
      mPicoEventMixer[vz_bin][centrality]->mixEvents();

    return kStOk;
}
// _________________________________________________________
Int_t StPicoMixedEventMaker::SetCategories() {
    return kStOk;
}
// _________________________________________________________
int StPicoMixedEventMaker::categorize(StPicoDst const * picoDst ) {
    StThreeVectorF pVertex = (picoDst->event())->primaryVertex();
    if( fabs(pVertex.z())>6.0 ) return -99;
    int bin = -6.0 + (pVertex.z()+6.0)/1.2;
    return bin;
}
