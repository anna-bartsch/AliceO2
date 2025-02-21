// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   Controller.h
/// @author ruben.shahoyan@cern.ch, michael.lettrich@cern.ch
/// @since  2021-02-01
/// @brief  Steering class for the global alignment

#include "Align/Controller.h"
#include "Align/AlignConfig.h"
#include "Framework/Logger.h"
#include "Align/utils.h"
#include "Align/AlignmentPoint.h"
#include "Align/AlignableDetector.h"
#include "Align/AlignableVolume.h"
#include "Align/AlignableDetectorITS.h"
//#include "Align/AlignableDetectorTPC.h"
//#include "Align/AlignableDetectorTRD.h"
//#include "Align/AlignableDetectorTOF.h"
#include "Align/EventVertex.h"
#include "Align/ResidualsControllerFast.h"
#include "Align/GeometricalConstraint.h"
#include "Align/DOFStatistics.h"
//#include "AliTrackerBase.h"
//#include "AliESDCosmicTrack.h"
//#include "AliESDtrack.h"
//#include "AliESDEvent.h"
//#include "AliESDVertex.h"
//#include "AliRecoParam.h"
//#include "AliCDBRunRange.h"
//#include "AliCDBManager.h"
//#include "AliCDBEntry.h"
#include "MathUtils/Utils.h"

#include <TMath.h>
#include <TString.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TRandom.h>
#include <TH1F.h>
#include <TList.h>
#include <cstdio>
#include <TGeoGlobalMagField.h>
#include "CommonUtils/NameConf.h"
#include "DataFormatsParameters/GRPObject.h"

using namespace TMath;
using namespace o2::align::utils;
using namespace o2::dataformats;
using namespace o2::globaltracking;

using std::ifstream;

namespace o2
{
namespace align
{

void Controller::ProcStat::print() const
{
  // TODO RS
  //  const Char_t* Controller::sStatClName[Controller::kNStatCl] = {"Inp: ", "Acc: "};
  //  const Char_t* Controller::sStatName[Controller::kMaxStat] =
  //  {"runs", "Ev.Coll", "Ev.Cosm", "Trc.Coll", "Trc.Cosm"};
}

const Char_t* Controller::sMPDataExt = ".mille";

const Char_t* Controller::sDetectorName[Controller::kNDetectors] = {"ITS", "TPC", "TRD", "TOF", "HMPID"}; //RSREM

//const int Controller::mgkSkipLayers[Controller::kNLrSkip] = {AliGeomManager::kPHOS1, AliGeomManager::kPHOS2,
//                                                                 AliGeomManager::kMUON, AliGeomManager::kEMCAL}; TODO(milettri, shahoian): needs detector IDs previously stored in AliGeomManager
const int Controller::sSkipLayers[Controller::kNLrSkip] = {0, 0, 0, 0}; // TODO(milettri, shahoian): needs AliGeomManager - remove this line after fix.


const Char_t* Controller::sHStatName[Controller::kNHVars] = {
  "Runs", "Ev.Inp", "Ev.VtxOK", "Tr.Inp", "Tr.2Fit", "Tr.2FitVC", "Tr.2PrMat", "Tr.2ResDer", "Tr.Stored", "Tr.Acc", "Tr.ContRes"};

//________________________________________________________________
Controller::Controller(DetID::mask_t detmask)
  : mDetMask(detmask)
{
  // def c-tor

  //  SetOutCDBRunRange();   FIXME(milettri): needs OCDB
  init();

  // run config macro if provided
  if (!getInitDOFsDone()) {
    initDOFs();
  }
  if (!getNDOFs()) {
    LOG(fatal) << "No DOFs found, initialization failed";
  }
}

//________________________________________________________________
Controller::~Controller()
{
  // d-tor
  if (mMPRecFile) {
    closeMPRecOutput();
  }
  if (mMille) {
    closeMilleOutput();
  }
  if (mResidFile) {
    closeResidOutput();
  }
  //
  for (int i = 0; i < DetID::nDetectors; i++) {
    delete mDetectors[i];
  }
  delete mHistoStat;
  //
}

//________________________________________________________________
void Controller::init()
{
  if (mDetMask[DetID::ITS]) {
    addDetector(new AlignableDetectorITS(this));
  }
  mVtxSens = std::make_unique<EventVertex>(this);
}

//________________________________________________________________
void Controller::process(const RecoContainer& recodata)
{
  /*
  auto creator = [](auto& _tr, GTrackID _origID, float t0, float terr) {
    if (!_origID.includesDet(DetID::ITS)) {
      return true; // just in case this selection was not done on RecoContainer filling level
    }
    if constexpr (isITSTrack<decltype(_tr)>()) {
      t0 += halfROFITS;  // ITS time is supplied in \mus as beginning of ROF
      terr *= hw2ErrITS; // error is supplied as a half-ROF duration, convert to \mus
    }
    // for all other tracks the time is in \mus with gaussian error
    if constexpr (std::is_base_of_v<o2::track::TrackParCov, std::decay_t<decltype(_tr)>>) {
      if (terr < maxTrackTimeError) {
        tracks.emplace_back(TrackWithTimeStamp{_tr, {t0, terr}});
        gids.emplace_back(_origID);
      }
    }

    return true;
  };
  recoData.createTracksVariadic(creator); // create track sample considered for vertexing
*/
}

//________________________________________________________________
void Controller::initDetectors()
{
  // init all detectors geometry
  //
  if (getInitGeomDone()) {
    return;
  }
  //

  //
  mAlgTrack = std::make_unique<AlignmentTrack>();
  mRefPoint = std::make_unique<AlignmentPoint>();
  //
  int dofCnt = 0;
  // special fake sensor for vertex constraint point
  // it has special T2L matrix adjusted for each track, no need to init it here
  mVtxSens = std::make_unique<EventVertex>(this);
  mVtxSens->setInternalID(1);
  mVtxSens->prepareMatrixL2G();
  mVtxSens->prepareMatrixL2GIdeal();
  dofCnt += mVtxSens->getNDOFs();
  //
  for (auto id = DetID::First; id <= DetID::Last; id++) {
    auto* det = getDetector(id);
    if (det) {
      dofCnt += det->initGeom();
    }
  }
  if (!dofCnt) {
    LOG(fatal) << "No DOFs found";
  }
  //
  //
  for (auto id = DetID::First; id <= DetID::Last; id++) {
    auto* det = getDetector(id);
    if (!det || det->isDisabled()) {
      continue;
    }
    det->cacheReferenceOCDB();
  }
  //
  assignDOFs();
  LOG(info) << "Booked " << dofCnt << " global parameters";
  //
  setInitGeomDone();
  //
}

//________________________________________________________________
void Controller::initDOFs()
{
  // scan all free global parameters, link detectors to array of params
  //
  if (getInitDOFsDone()) {
    LOG(info) << "initDOFs was already done, just reassigning " << getNDOFs() << "DOFs arrays/labels";
    assignDOFs();
    return;
  }
  const auto& conf = AlignConfig::Instance();

  mNDOFs = 0;
  int ndfAct = 0;
  assignDOFs();
  int nact = 0;
  mVtxSens->initDOFs();
  for (auto id = DetID::First; id <= DetID::Last; id++) {
    AlignableDetector* det = getDetector(id);
    if (det && !det->isDisabled()) {
      det->initDOFs();
      nact++;
      ndfAct += det->getNDOFsTot();
    }
  }
  for (int i = 0; i < NTrackTypes; i++) {
    if (nact < conf.minDetAcc[i]) {
      LOG(fatal) << nact << " detectors are active, while " << conf.minDetAcc[i] << " in track are asked";
    }
  }
  LOG(info) << mNDOFs << " global parameters " << mNDet << " detectors, " << ndfAct << " in " << nact << " active detectors";
  addAutoConstraints();
  setInitDOFsDone();
}

//________________________________________________________________
void Controller::assignDOFs()
{
  // add parameters/labels arrays to volumes. If the Controller is read from the file, this method need
  // to be called (of initDOFs should be called)
  //
  int ndfOld = -1;
  if (mNDOFs > 0) {
    ndfOld = mNDOFs;
  }
  mNDOFs = 0;
  //
  // reserve
  int ndofTOT = mVtxSens->getNDOFs();
  for (auto id = DetID::First; id <= DetID::Last; id++) {
    AlignableDetector* det = getDetector(id);
    if (!det) {
      continue;
    }
    ndofTOT += det->getNDOFsTot();
  }
  mGloParVal.clear();
  mGloParErr.clear();
  mGloParLab.clear();
  mOrderedLbl.clear();
  mLbl2ID.clear();
  mGloParVal.reserve(ndofTOT);
  mGloParErr.reserve(ndofTOT);
  mGloParLab.reserve(ndofTOT);
  mOrderedLbl.reserve(ndofTOT);
  mLbl2ID.reserve(ndofTOT);

  mVtxSens->assignDOFs();

  for (auto id = DetID::First; id <= DetID::Last; id++) {
    AlignableDetector* det = getDetector(id);
    if (!det) {
      continue;
    }
    mNDOFs += det->assignDOFs();
  }
  LOG(info) << "Assigned parameters/labels arrays for " << mNDOFs << " DOFs";
  if (ndfOld > -1 && ndfOld != mNDOFs) {
    LOG(error) << "Recalculated NDOFs=" << mNDOFs << " not equal to saved NDOFs=" << ndfOld;
  }
  //
  // build Lbl <-> parID table
  /* FIXME RS TODO
  Sort(mNDOFs, mGloParLab, mLbl2ID, false); // sort in increasing order
  for (int i = mNDOFs; i--;) {
    mOrderedLbl[i] = mGloParLab[mLbl2ID[i]];
  }
  */
  //
}

//_________________________________________________________
void Controller::addDetector(AlignableDetector* det)
{
  // add detector constructed externally to alignment framework
  mDetectors[det->getDetID()] = det;
  mNDet++;
}

//_________________________________________________________
bool Controller::checkDetectorPattern(DetID::mask_t patt) const
{
  //validate detector pattern
  return ((patt & mObligatoryDetPattern[mTracksType]) == mObligatoryDetPattern[mTracksType]) &&
         patt.count() >= AlignConfig::Instance().minDetAcc[mTracksType];
}

//_________________________________________________________
bool Controller::checkDetectorPoints(const int* npsel) const
{
  //validate detectors pattern according to number of selected points
  int ndOK = 0;
  for (auto id = DetID::First; id <= DetID::Last; id++) {
    AlignableDetector* det = getDetector(id);
    if (!det || det->isDisabled(mTracksType)) {
      continue;
    }
    if (npsel[id] < det->getNPointsSel(mTracksType)) {
      if (det->isObligatory(mTracksType)) {
        return false;
      }
      continue;
    }
    ndOK++;
  }
  return ndOK >= AlignConfig::Instance().minDetAcc[mTracksType];
}

//FIXME(milettri): needs AliESDtrack
////_________________________________________________________
//uint32_t Controller::AcceptTrack(const AliESDtrack* esdTr, bool strict) const
//{
//  // decide if the track should be processed
//  AlignableDetector* det = 0;
//  uint32_t detAcc = 0;
//  if (mFieldOn && esdTr->Pt() < mPtMin[mTracksType]){
//    return 0;}
//  if (Abs(esdTr->Eta()) > mEtaMax[mTracksType]){
//    return 0;}
//  //
//  for (auto id=DetID::First; id<=DetID::Last; id++) {
//
//    if (!(det = getDetector(id)) || det->isDisabled(mTracksType)){
//      continue;}
//    if (!det->AcceptTrack(esdTr, mTracksType)) {
//      if (strict && det->isObligatory(mTracksType)){
//        return 0;}
//      else
//        continue;
//    }
//    //
//    detAcc |= 0x1 << idet;
//  }
//  if (numberOfBitsSet(detAcc) < mMinDetAcc[mTracksType]){
//    return 0;}
//  return detAcc;
//  //
//}

//FIXME(milettri): needs AliESDtrack
////_________________________________________________________
//uint32_t Controller::AcceptTrackCosmic(const AliESDtrack* esdPairCosm[kNCosmLegs]) const
//{
//  // decide if the pair of tracks making cosmic track should be processed
//  uint32_t detAcc = 0, detAccLeg;
//  for (int i = kNCosmLegs; i--;) {
//    detAccLeg = AcceptTrack(esdPairCosm[i], mCosmicSelStrict); // missing obligatory detectors in one leg might be allowed
//    if (!detAccLeg){
//      return 0;}
//    detAcc |= detAccLeg;
//  }
//  if (mCosmicSelStrict){
//    return detAcc;}
//  //
//  // for non-stric selection check convolution of detector presence
//  if (!checkDetectorPattern(detAcc)){
//    return 0;}
//  return detAcc;
//  //
//}

//FIXME(milettri): needs AliESDEvent
////_________________________________________________________
//void Controller::SetESDEvent(const AliESDEvent* ev)
//{
//  // attach event to analyse
//  fESDEvent = ev;
//  // setup magnetic field if needed
//  if (fESDEvent &&
//      (!TGeoGlobalMagField::Instance()->GetField() ||
//       !smallerAbs(fESDEvent->GetMagneticField() - AliTrackerBase::GetBz(), 5e-4))) {
//    fESDEvent->InitMagneticField();
//  }
//}

//FIXME(milettri): needs AliESDEvent
////_________________________________________________________
//bool Controller::ProcessEvent(const AliESDEvent* esdEv)
//{
//  // process event
//  const int kProcStatFreq = 100;
//  static int evCount = 0;
//  if (!(evCount % kProcStatFreq)) {
//    ProcInfo_t procInfo;
//    gSystem->GetProcInfo(&procInfo);
//    LOG(info) << "ProcStat: CPUusr:" << int(procInfo.fCpuUser) << " CPUsys:" << int(procInfo.fCpuSys) << " RMem:" << int(procInfo.fMemResident / 1024) << " VMem:" << int(procInfo.fMemVirtual / 1024);
//  }
//  evCount++;
//  //
//  SetESDEvent(esdEv);
//  //
//  if (esdEv->getRunNumber() != getRunNumber()){
//    SetRunNumber(esdEv->getRunNumber());
//  }
//  //
//  setCosmic(esdEv->GetEventSpecie() == AliRecoParam::kCosmic ||
//            (esdEv->GetNumberOfCosmicTracks() > 0 && !esdEv->GetPrimaryVertexTracks()->GetStatus()));
//  //
//  fillStatHisto(kEvInp);
//  //
//#if DEBUG > 2
//  LOG << "Processing event " << esdEv->GetEventNumberInFile() << " of ev.specie " << esdEv->GetEventSpecie() << " -> Ntr: " << esdEv->GetNumberOfTracks() << " NtrCosm: " << esdEv->GetNumberOfCosmicTracks();
//#endif
//  //
//  setFieldOn(Abs(esdEv->GetMagneticField()) > kAlmost0Field);
//  if (!isCosmic() && !CheckSetVertex(esdEv->GetPrimaryVertexTracks())){
//    return false;}
//  fillStatHisto(kEvVtx);
//  //
//  int ntr = 0, accTr = 0;
//  if (isCosmic()) {
//    mStat[kInpStat][kEventCosm]++;
//    ntr = esdEv->GetNumberOfCosmicTracks();
//    fillStatHisto(kTrackInp, ntr);
//    for (int itr = 0; itr < ntr; itr++) {
//      accTr += ProcessTrack(esdEv->GetCosmicTrack(itr));
//    }
//    if (accTr){
//      mStat[kAccStat][kEventCosm]++;}
//  } else {
//    mStat[kInpStat][kEventColl]++;
//    ntr = esdEv->GetNumberOfTracks();
//    fillStatHisto(kTrackInp, ntr);
//    for (int itr = 0; itr < ntr; itr++) {
//      //      int accTrOld = accTr;
//      accTr += ProcessTrack(esdEv->GetTrack(itr));
//      /*
//      if (accTr>accTrOld && mCResid) {
// int ndf = mCResid.getNPoints()*2-5;
// if (mCResid.getChi2()/ndf>20 || !mCResid.getKalmanDone()
//     || mCResid.getChi2K()/ndf>20) {
//   printf("BAD FIT for %d\n",itr);
// }
//        mCResid.Print("er");
//      }
//      */
//    }
//    if (accTr){
//      mStat[kAccStat][kEventColl]++;}
//  }
//  //
//  fillStatHisto(kTrackAcc, accTr);
//  //
//  if (accTr) {
//    LOG(info) << "Processed event " << esdEv->GetEventNumberInFile() << " of ev.specie " << esdEv->GetEventSpecie() << " -> Accepted: " << accTr << " of " << ntr << " tracks";
//  }
//  return true;
//}

//FIXME(milettri): needs AliESDtrack
//_________________________________________________________
//bool Controller::ProcessTrack(const AliESDtrack* esdTr)
//{
//  // process single track
//  //
//  mStat[kInpStat][kTrackColl]++;
//  fESDTrack[0] = esdTr;
//  fESDTrack[1] = 0;
//  //
//  int nPnt = 0;
//  const AliESDfriendTrack* trF = esdTr->GetFriendTrack();
//  if (!trF){
//    return false;}
//  const AliTrackPointArray* trPoints = trF->GetTrackPointArray();
//  if (!trPoints || (nPnt = trPoints->GetNPoints()) < 1){
//    return false;}
//  //
//  uint32_t detAcc = AcceptTrack(esdTr);
//  if (!detAcc){
//    return false;}
//  //
//  resetDetectors();
//  mAlgTrack->Clear();
//  //
//  // process the track points for each detector,
//  AlignableDetector* det = 0;
//  for (auto id=DetID::First; id<=DetID::Last; id++) {
//    if (!(detAcc & (0x1 << idet))){  // RS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//      continue;}
//    det = getDetector(id);
//    if (det->ProcessPoints(esdTr, mAlgTrack) < det->getNPointsSel(kColl)) {
//      detAcc &= ~(0x1 << idet); // did not survive, suppress detector in the track
//      if (det->isObligatory(kColl)){
//        return false;}
//    }
//    if (numberOfBitsSet(detAcc) < mMinDetAcc[kColl]){
//      return false;} // abandon track
//  }
//  //
//  if (mAlgTrack->getNPoints() < getMinPoints()){
//    return false;}
//  // fill needed points (tracking frame) in the mAlgTrack
//  mRefPoint->setContainsMeasurement(false);
//  mRefPoint->setContainsMaterial(false);
//  mAlgTrack->addPoint(mRefPoint); // reference point which the track will refer to
//  //
//  mAlgTrack->copyFrom(esdTr);
//  if (!getFieldOn()){
//    mAlgTrack->imposePtBOff(mDefPtBOff[utils::kColl]);}
//  mAlgTrack->setFieldON(getFieldOn());
//  mAlgTrack->sortPoints();
//  //
//  // at this stage the points are sorted from maxX to minX, the latter corresponding to
//  // reference point (e.g. vertex) with X~0. The mAlgTrack->getInnerPointID() points on it,
//  // hence mAlgTrack->getInnerPointID() is the 1st really measured point. We will set the
//  // alpha of the reference point to alpha of the barrel sector corresponding to this
//  // 1st measured point
//  int pntMeas = mAlgTrack->getInnerPointID() - 1;
//  if (pntMeas < 0) { // this should not happen
//    mAlgTrack->Print("p meas");
//    LOG(fatal) << "AlignmentTrack->getInnerPointID() cannot be 0";
//  }
//  // do we want to add the vertex as a measured point ?
//  if (!addVertexConstraint()) { // no constrain, just reference point w/o measurement
//    mRefPoint->setXYZTracking(0, 0, 0);
//    mRefPoint->setAlphaSens(sector2Alpha(mAlgTrack->getPoint(pntMeas)->getAliceSector()));
//  } else
//    fillStatHisto(kTrackFitInpVC);
//  //
//  fillStatHisto(kTrackFitInp);
//  if (!mAlgTrack->iniFit()){
//    return false;}
//  fillStatHisto(kTrackProcMatInp);
//  if (!mAlgTrack->processMaterials()){
//    return false;}
//  mAlgTrack->defineDOFs();
//  //
//  fillStatHisto(kTrackResDerInp);
//  if (!mAlgTrack->calcResidDeriv()){
//    return false;}
//  //
//  if (!storeProcessedTrack(mMPOutType & ~kContR)){
//    return false;} // store derivatives for MP
//  //
//  if (getProduceControlRes() &&                                   // need control residuals, ignore selection fraction if this is the
//      (mMPOutType == kContR || gRandom->Rndm() < mControlFrac)) { // output requested
//    if (!testLocalSolution() || !storeProcessedTrack(kContR)){
//      return false;}
//  }
//  //
//  fillStatHisto(kTrackStore);
//  //
//  mStat[kAccStat][kTrackColl]++;
//  //
//  return true;
//}

//FIXME(milettri): needs AliESDVertex
////_________________________________________________________
//bool Controller::CheckSetVertex(const AliESDVertex* vtx)
//{
//  // vertex selection/constraint check
//  if (!vtx) {
//    fVertex = 0;
//    return true;
//  }
//  int ncont = vtx->GetNContributors();
//  if (mVtxMinCont > 0 && mVtxMinCont > ncont) {
//#if DEBUG > 2
//    LOG(info) << "Rejecting event with " << % d << " vertex contributors (min " << % d << " asked)", ncont, mVtxMinCont);
//#endif
//    return false;
//  }
//  if (mVtxMaxCont > 0 && ncont > mVtxMaxCont) {
//#if DEBUG > 2
//    LOG(info) << "Rejecting event with " << % d << " vertex contributors (max " << % d << " asked)", ncont, mVtxMaxCont);
//#endif
//    return false;
//  }
//  fVertex = (ncont >= mVtxMinContVC) ? vtx : 0; // use vertex as a constraint
//  return true;
//}

//FIXME(milettri): needs AliESDCosmicTrack
////_________________________________________________________
//bool Controller::ProcessTrack(const AliESDCosmicTrack* cosmTr)
//{
//  // process single cosmic track
//  //
//  mStat[kInpStat][kTrackCosm]++;
//  int nPnt = 0;
//  fESDTrack[0] = 0;
//  fESDTrack[1] = 0;
//  //
//  for (int leg = kNCosmLegs; leg--;) {
//    const AliESDtrack* esdTr =
//      fESDEvent->GetTrack(leg == kCosmLow ? cosmTr->GetESDLowerTrackIndex() : cosmTr->GetESDUpperTrackIndex());
//    const AliESDfriendTrack* trF = esdTr->GetFriendTrack();
//    if (!trF){
//      return false;}
//    const AliTrackPointArray* trPoints = trF->GetTrackPointArray();
//    if (!trPoints || (nPnt += trPoints->GetNPoints()) < 1){
//      return false;}
//    //
//    fESDTrack[leg] = esdTr;
//  }
//  //
//  uint32_t detAcc = AcceptTrackCosmic(fESDTrack);
//  if (!detAcc){
//    return false;}
//  //
//  resetDetectors();
//  mAlgTrack->Clear();
//  mAlgTrack->setCosmic(true);
//  //
//  // process the track points for each detector,
//  // fill needed points (tracking frame) in the mAlgTrack
//  mRefPoint->setContainsMeasurement(false);
//  mRefPoint->setContainsMaterial(false);
//  mAlgTrack->addPoint(mRefPoint); // reference point which the track will refer to
//  //
//  AlignableDetector* det = 0;
//  int npsel[kNDetectors] = {0};
//  for (int nPleg = 0, leg = kNCosmLegs; leg--;) {
//    for (auto id=DetID::First; id<=DetID::Last; id++) {
//      if (!(detAcc & (0x1 << idet))){
//        continue;}
//      det = getDetector(id);
//      //
//      // upper leg points marked as the track going in inverse direction
//      int np = det->ProcessPoints(fESDTrack[leg], mAlgTrack, leg == kCosmUp);
//      if (np < det->getNPointsSel(Cosm) && mCosmicSelStrict &&
//          det->isObligatory(Cosm))
//        return false;
//      npsel[id] += np;
//      nPleg += np;
//    }
//    if (nPleg < getMinPoints()){
//      return false;}
//  }
//  // last check on legs-combined patter
//  if (!checkDetectorPoints(npsel)){
//    return false;}
//  //
//  mAlgTrack->copyFrom(cosmTr);
//  if (!getFieldOn()){
//    mAlgTrack->imposePtBOff(mDefPtBOff[utils::Cosm]);}
//  mAlgTrack->setFieldON(getFieldOn());
//  mAlgTrack->sortPoints();
//  //
//  // at this stage the points are sorted from maxX to minX, the latter corresponding to
//  // reference point (e.g. vertex) with X~0. The mAlgTrack->getInnerPointID() points on it,
//  // hence mAlgTrack->getInnerPointID() is the 1st really measured point. We will set the
//  // alpha of the reference point to alpha of the barrel sector corresponding to this
//  // 1st measured point
//  int pntMeas = mAlgTrack->getInnerPointID() - 1;
//  if (pntMeas < 0) { // this should not happen
//    mAlgTrack->Print("p meas");
//    LOG(fatal) << "AlignmentTrack->getInnerPointID() cannot be 0";
//  }
//  mRefPoint->setAlphaSens(sector2Alpha(mAlgTrack->getPoint(pntMeas)->getAliceSector()));
//  //
//  fillStatHisto(kTrackFitInp);
//  if (!mAlgTrack->iniFit()){
//    return false;}
//  //
//  fillStatHisto(kTrackProcMatInp);
//  if (!mAlgTrack->processMaterials()){
//    return false;}
//  mAlgTrack->defineDOFs();
//  //
//  fillStatHisto(kTrackResDerInp);
//  if (!mAlgTrack->calcResidDeriv()){
//    return false;}
//  //
//  if (!storeProcessedTrack(mMPOutType & ~kContR)){
//    return false;} // store derivatives for MP
//  //
//  if (getProduceControlRes() &&                                   // need control residuals, ignore selection fraction if this is the
//      (mMPOutType == kContR || gRandom->Rndm() < mControlFrac)) { // output requested
//    if (!testLocalSolution() || !storeProcessedTrack(kContR)){
//      return false;}
//  }
//  //
//  fillStatHisto(kTrackStore);
//  mStat[kAccStat][kTrackCosm]++;
//  return true;
//}

//_________________________________________________________
bool Controller::storeProcessedTrack(int what)
{
  // write alignment track
  bool res = true;
  if ((what & kMille)) {
    res &= fillMilleData();
  }
  if ((what & kMPRec)) {
    res &= fillMPRecData();
  }
  if ((what & kContR)) {
    res &= fillControlData();
  }
  //
  return res;
}

//_________________________________________________________
bool Controller::fillMilleData()
{
  // store MP2 data in Mille format
  if (!mMille) {
    mMille = std::make_unique<Mille>(fmt::format("{}{}", mMPDatFileName, sMPDataExt).c_str());
  }
  //
  if (!mAlgTrack->getDerivDone()) {
    LOG(error) << "Track derivatives are not yet evaluated";
    return false;
  }
  int np(mAlgTrack->getNPoints()), nDGloTot(0); // total number global derivatives stored
  int nParETP(mAlgTrack->getNLocExtPar());      // numnber of local parameters for reference track param
  int nVarLoc(mAlgTrack->getNLocPar());         // number of local degrees of freedom in the track
  float *buffDL(nullptr), *buffDG(nullptr);     // faster acces arrays
  int* buffI(nullptr);
  //
  const int* gloParID(mAlgTrack->getGloParID()); // IDs of global DOFs this track depends on
  for (int ip = 0; ip < np; ip++) {
    AlignmentPoint* pnt = mAlgTrack->getPoint(ip);
    if (pnt->containsMeasurement()) {
      int gloOffs = pnt->getDGloOffs(); // 1st entry of global derivatives for this point
      int nDGlo = pnt->getNGloDOFs();   // number of global derivatives (number of DOFs it depends on)
      if (!pnt->isStatOK()) {
        pnt->incrementStat();
      }
      // check buffer sizes
      {
        if (mMilleDBuffer.GetSize() < nVarLoc + nDGlo) {
          mMilleDBuffer.Set(100 + nVarLoc + nDGlo);
        }
        if (mMilleIBuffer.GetSize() < nDGlo) {
          mMilleIBuffer.Set(100 + nDGlo);
        }
        buffDL = mMilleDBuffer.GetArray(); // faster acces
        buffDG = buffDL + nVarLoc;         // faster acces
        buffI = mMilleIBuffer.GetArray();  // faster acces
      }
      // local der. array cannot be 0-suppressed by Mille construction, need to reset all to 0
      //
      for (int idim = 0; idim < 2; idim++) { // 2 dimensional orthogonal measurement
        memset(buffDL, 0, nVarLoc * sizeof(float));
        const double* deriv = mAlgTrack->getDResDLoc(idim, ip); // array of Dresidual/Dparams_loc
        // derivatives over reference track parameters
        for (int j = 0; j < nParETP; j++) {
          buffDL[j] = (isZeroAbs(deriv[j])) ? 0 : deriv[j];
        }
        //
        // point may depend on material variables within these limits
        int lp0 = pnt->getMinLocVarID(), lp1 = pnt->getMaxLocVarID();
        for (int j = lp0; j < lp1; j++) {
          buffDL[j] = (isZeroAbs(deriv[j])) ? 0 : deriv[j];
        }
        //
        // derivatives over global params: this array can be 0-suppressed, no need to reset
        int nGlo(0);
        deriv = mAlgTrack->getDResDGlo(idim, gloOffs);
        const int* gloIDP(gloParID + gloOffs);
        for (int j = 0; j < nDGlo; j++) {
          if (!isZeroAbs(deriv[j])) {
            buffDG[nGlo] = deriv[j];                 // value of derivative
            buffI[nGlo++] = getGloParLab(gloIDP[j]); // global DOF ID + 1 (Millepede needs positive labels)
          }
        }
        mMille->mille(nVarLoc, buffDL, nGlo, buffDG, buffI,
                      mAlgTrack->getResidual(idim, ip), Sqrt(pnt->getErrDiag(idim)));
        nDGloTot += nGlo;
        //
      }
    }
    if (pnt->containsMaterial()) { // material point can add 4 or 5 otrhogonal pseudo-measurements
      memset(buffDL, 0, nVarLoc * sizeof(float));
      int nmatpar = pnt->getNMatPar(); // residuals (correction expectation value)
      //      const float* expMatCorr = pnt->getMatCorrExp(); // expected corrections (diagonalized)
      const float* expMatCov = pnt->getMatCorrCov(); // their diagonalized error matrix
      int offs = pnt->getMaxLocVarID() - nmatpar;    // start of material variables
      // here all derivatives are 1 = dx/dx
      for (int j = 0; j < nmatpar; j++) { // mat. "measurements" don't depend on global params
        int j1 = j + offs;
        buffDL[j1] = 1.0; // only 1 non-0 derivative
        //mMille->mille(nVarLoc,buffDL,0,buffDG,buffI,expMatCorr[j],Sqrt(expMatCov[j]));
        // expectation for MS effect is 0
        mMille->mille(nVarLoc, buffDL, 0, buffDG, buffI, 0, Sqrt(expMatCov[j]));
        buffDL[j1] = 0.0; // reset buffer
      }
    } // material "measurement"
  }   // loop over points
  //
  if (!nDGloTot) {
    LOG(info) << "Track does not depend on free global parameters, discard";
    mMille->kill();
    return false;
  }
  mMille->end(); // store the record
  return true;
}

//_________________________________________________________
bool Controller::fillMPRecData()
{
  LOG(fatal) << __PRETTY_FUNCTION__ << " is disabled";
  //FIXME(milettri): needs AliESDEvent
  //  // store MP2 in MPRecord format
  //  if (!mMPRecord){
  //    initMPRecOutput();}
  //  //
  //  mMPRecord->Clear();
  //  if (!mMPRecord->fillTrack(mAlgTrack, mGloParLab)){
  //    return false;}
  //  mMPRecord->SetRun(mRunNumber);
  //  mMPRecord->setTimeStamp(fESDEvent->GetTimeStamp());
  //  uint32_t tID = 0xffff & uint(fESDTrack[0]->GetID());
  //  if (isCosmic()){
  //    tID |= (0xffff & uint(fESDTrack[1]->GetID())) << 16;}
  //  mMPRecord->setTrackID(tID);
  //  mMPRecTree->Fill();
  return true;
}

//_________________________________________________________
bool Controller::fillControlData()
{
  LOG(fatal) << __PRETTY_FUNCTION__ << " is disabled";
  //FIXME(milettri): needs AliESDEvent
  //  // store control residuals
  //  if (!mCResid){
  //    initResidOutput();}
  //  //
  //  int nps, np = mAlgTrack->getNPoints();
  //  nps = (!mRefPoint->containsMeasurement()) ? np - 1 : np; // ref point is dummy?
  //  if (nps < 0){
  //    return true;}
  //  //
  //  mCResid.Clear();
  //  if (!mCResid.fillTrack(mAlgTrack, mDoKalmanResid)){
  //    return false;}
  //  mCResid.setRun(mRunNumber);
  //  mCResid.setTimeStamp(fESDEvent->GetTimeStamp());
  //  mCResid.setBz(fESDEvent->GetMagneticField());
  //  uint32_t tID = 0xffff & uint(fESDTrack[0]->GetID());
  //  if (isCosmic()){
  //    tID |= (0xffff & uint(fESDTrack[1]->GetID())) << 16;}
  //  mCResid.setTrackID(tID);
  //  //
  //  mResidTree->Fill();
  //  fillStatHisto(kTrackControl);
  //  //
  return true;
}

//_________________________________________________________
void Controller::setRunNumber(int run)
{
  if (run == mRunNumber) {
    return;
  } // nothing to do
  //
  acknowledgeNewRun(run);
}

//_________________________________________________________
void Controller::acknowledgeNewRun(int run)
{
  LOG(warning) << __PRETTY_FUNCTION__ << " yet incomplete";

  o2::base::GeometryManager::loadGeometry();
  o2::base::PropagatorImpl<double>::initFieldFromGRP();
  std::unique_ptr<o2::parameters::GRPObject> grp{o2::parameters::GRPObject::loadFrom()};

  //FIXME(milettri): needs AliESDEvent
  //  // load needed info for new run
  //  if (run == mRunNumber){
  //    return;} // nothing to do
  //  if (run > 0) {
  //    mStat[kAccStat][kRun]++;
  //  }
  //  if (mRunNumber > 0){
  //    fillStatHisto(kRunDone);}
  //  mRunNumber = run;
  //  LOG(info) << "Processing new run " << mRunNumber;
  //  //
  //  // setup magnetic field
  //  if (fESDEvent &&
  //      (!TGeoGlobalMagField::Instance()->GetField() ||
  //       !smallerAbs(fESDEvent->GetMagneticField() - AliTrackerBase::GetBz(), 5e-4))) {
  //    fESDEvent->InitMagneticField();
  //  }
  //  //
  //  if (!mUseRecoOCDB) {
  //    LOG(warning) << "Reco-time OCDB will NOT be preloaded";
  //    return;
  //  }
  //  LoadRecoTimeOCDB();
  //  //
  //  for (auto id=DetID::First; id<=DetID::Last; id++) {
  //    AlignableDetector* det = getDetector(id);
  //    if (!det->isDisabled()){
  //      det->acknowledgeNewRun(run);}
  //  }
  //  //
  //  // bring to virgin state
  //  // CleanOCDB();
  //  //
  //  // LoadRefOCDB(); //??? we need to get back reference OCDB ???
  //  //
  //  mStat[kInpStat][kRun]++;
  //  //
}

// FIXME(milettri): needs OCDB
////_________________________________________________________
//bool Controller::LoadRecoTimeOCDB()
//{
//  // Load OCDB paths used for the reconstruction of data being processed
//  // In order to avoid unnecessary uploads, the objects are not actually
//  // loaded/cached but just added as specific paths with version
//  LOG(info) << "Preloading Reco-Time OCDB for run " << mRunNumber << " from ESD UserInfo list";
//  //
//  CleanOCDB();
//  //
//  if (!mRecoOCDBConf.IsNull() && !gSystem->AccessPathName(mRecoOCDBConf.c_str(), kFileExists)) {
//    LOG(info) << "Executing reco-time OCDB setup macro " << mRecoOCDBConf.c_str();
//    gROOT->ProcessLine(Form(".x %s(%d)", mRecoOCDBConf.c_str(), mRunNumber));
//    if (AliCDBManager::Instance()->IsDefaultStorageSet()){
//      return true;}
//    LOG(fatal) << "macro " << mRecoOCDBConf.c_str() << " failed to configure reco-time OCDB";
//  } else
//    LOG(warning) << "No reco-time OCDB config macro" << mRecoOCDBConf.c_str() << "  is found, will use ESD:UserInfo";
//  //
//  if (!mESDTree){
//    LOG(fatal) << "Cannot preload Reco-Time OCDB since the ESD tree is not set";}
//  const TTree* tr = mESDTree; // go the the real ESD tree
//  while (tr->GetTree() && tr->GetTree() != tr)
//    tr = tr->GetTree();
//  //
//  const TList* userInfo = const_cast<TTree*>(tr)->GetUserInfo();
//  TMap* cdbMap = (TMap*)userInfo->FindObject("cdbMap");
//  TList* cdbList = (TList*)userInfo->FindObject("cdbList");
//  //
//  if (!cdbMap || !cdbList) {
//    userInfo->Print();
//    LOG(fatal) << "Failed to extract cdbMap and cdbList from UserInfo list";
//  }
//  //
//  return PreloadOCDB(mRunNumber, cdbMap, cdbList);
//}

//____________________________________________
void Controller::Print(const Option_t* opt) const
{
  // print info
  TString opts = opt;
  opts.ToLower();
  printf("%5d DOFs in %d detectors\n", mNDOFs, mNDet);
  if (getMPAlignDone()) {
    printf("ALIGNMENT FROM MILLEPEDE SOLUTION IS APPLIED\n");
  }
  //
  for (auto id = DetID::First; id <= DetID::Last; id++) {
    AlignableDetector* det = getDetector(id);
    if (!det) {
      continue;
    }
    det->Print(opt);
  }
  if (!opts.IsNull()) {
    printf("\nSpecial sensor for Vertex Constraint\n");
    mVtxSens->Print(opt);
  }
  //
  if (mRefRunNumber >= 0) {
    printf("(%d)", mRefRunNumber);
  }
  //
  printf("%-40s:\t%s\n", "Output OCDB path", mOutCDBPath.c_str());
  printf("%-40s:\t%s/%s\n", "Output OCDB comment/responsible",
         mOutCDBComment.c_str(), mOutCDBResponsible.c_str());
  printf("%-40s:\t%6d:%6d\n", "Output OCDB run range", mOutCDBRunRange[0], mOutCDBRunRange[1]);
  //
  printf("%-40s:\t%s\n", "Filename for MillePede steering", mMPSteerFileName.c_str());
  printf("%-40s:\t%s\n", "Filename for MillePede parameters", mMPParFileName.c_str());
  printf("%-40s:\t%s\n", "Filename for MillePede constraints", mMPConFileName.c_str());
  printf("%-40s:\t%s\n", "Filename for control residuals:", mResidFileName.c_str());
  printf("%-40s:\t%.3f\n", "Fraction of control tracks", mControlFrac);
  printf("MPData output :\t");
  if (getProduceMPData()) {
    printf("%s%s ", mMPDatFileName.c_str(), sMPDataExt);
  }
  if (getProduceMPRecord()) {
    printf("%s%s ", mMPDatFileName.c_str(), ".root");
  }
  printf("\n");
  //
  if (opts.Contains("stat")) {
    printStatistics();
  }

  if (opts.Contains("conf")) {
    AlignConfig::Instance().printKeyValues(true);
  }
}

//________________________________________________________
void Controller::printStatistics() const
{
  // print processing stat
  mStat.print();
}

//________________________________________________________
void Controller::resetDetectors()
{
  // reset detectors for next track
  mRefPoint->Clear();
  for (auto id = DetID::First; id <= DetID::Last; id++) {
    AlignableDetector* det = getDetector(id);
    if (det) {
      det->resetPool(); // reset used alignment points
    }
  }
}

//____________________________________________
bool Controller::testLocalSolution()
{
  LOG(fatal) << __PRETTY_FUNCTION__ << " is disabled";
  //FIXME(milettri): needs AliSymMatrix
  //  // test track local solution
  //  TVectorD rhs;
  //  AliSymMatrix* mat = BuildMatrix(rhs);
  //  if (!mat){
  //    return false;}
  //  //  mat->Print("long data");
  //  //  rhs.Print();
  //  TVectorD vsl(rhs);
  //  if (!mat->SolveChol(rhs, vsl, true)) {
  //    delete mat;
  //    return false;
  //  }
  //  //
  //  /*
  //  // print solution vector
  //  int nlocpar = mAlgTrack->getNLocPar();
  //  int nlocparETP = mAlgTrack->getNLocExtPar(); // parameters of external track param
  //  printf("ETP Update: ");
  //  for (int i=0;i<nlocparETP;i++) printf("%+.2e(%+.2e) ",vsl[i],Sqrt((*mat)(i,i))); printf("\n");
  //  //
  //  if (nlocpar>nlocparETP) printf("Mat.Corr. update:\n");
  //  for (int ip=mAlgTrack->getNPoints();ip--;) {
  //    AlignmentPoint* pnt = mAlgTrack->getPoint(ip);
  //    int npm = pnt->getNMatPar();
  //    const float* expMatCov  = pnt->getMatCorrCov(); // its error
  //    int offs  = pnt->getMaxLocVarID() - npm;
  //    for (int ipar=0;ipar<npm;ipar++) {
  //      int parI = offs + ipar;
  //      double err = Sqrt(expMatCov[ipar]);
  //      printf("Pnt:%3d MatVar:%d DOF %3d | %+.3e(%+.3e) -> sig:%+.3e -> pull: %+.2e\n",
  //            ip,ipar,parI,vsl[parI],Sqrt((*mat)(parI,parI)), err,vsl[parI]/err);
  //    }
  //  }
  //  */
  //  //
  //  // increment current params by new solution
  //  rhs.SetElements(mAlgTrack->getLocPars());
  //  vsl += rhs;
  //  mAlgTrack->setLocPars(vsl.GetMatrixArray());
  //  mAlgTrack->calcResiduals();
  //  delete mat;
  //  //
  return true;
}

//FIXME(milettri): needs AliSymMatrix
////____________________________________________
//AliSymMatrix* Controller::BuildMatrix(TVectorD& vec)
//{
//  // build matrix/vector for local track
//  int npnt = mAlgTrack->getNPoints();
//  int nlocpar = mAlgTrack->getNLocPar();
//  //
//  vec.ResizeTo(nlocpar);
//  memset(vec.GetMatrixArray(), 0, nlocpar * sizeof(double));
//  AliSymMatrix* matp = new AliSymMatrix(nlocpar);
//  AliSymMatrix& mat = *matp;
//  //
//  for (int ip = npnt; ip--;) {
//    AlignmentPoint* pnt = mAlgTrack->getPoint(ip);
//    //
//    if (pnt->containsMeasurement()) {
//      //      pnt->Print("meas");
//      for (int idim = 2; idim--;) {                       // each point has 2 position residuals
//        double sigma2 = pnt->getErrDiag(idim);            // residual error
//        double resid = mAlgTrack->getResidual(idim, ip);  // residual
//        double* deriv = mAlgTrack->getDResDLoc(idim, ip); // array of Dresidual/Dparams
//        //
//        double sg2inv = 1. / sigma2;
//        for (int parI = nlocpar; parI--;) {
//          vec[parI] -= deriv[parI] * resid * sg2inv;
//          //  printf("%d %d %d %+e %+e %+e -> %+e\n",ip,idim,parI,sg2inv,deriv[parI],resid,vec[parI]);
//          //   for (int parJ=nlocpar;parJ--;) {
//          for (int parJ = parI + 1; parJ--;) {
//            mat(parI, parJ) += deriv[parI] * deriv[parJ] * sg2inv;
//          }
//        }
//      } // loop over 2 orthogonal measurements at the point
//    }   // derivarives at measured points
//    //
//    // if the point contains material, consider its expected kinks, eloss
//    // as measurements
//    if (pnt->containsMaterial()) {
//      // at least 4 parameters: 2 spatial + 2 angular kinks with 0 expectaction
//      int npm = pnt->getNMatPar();
//      // const float* expMatCorr = pnt->getMatCorrExp(); // expected correction (diagonalized)
//      const float* expMatCov = pnt->getMatCorrCov(); // its error
//      int offs = pnt->getMaxLocVarID() - npm;
//      for (int ipar = 0; ipar < npm; ipar++) {
//        int parI = offs + ipar;
//        // expected
//        // vec[parI] -= expMatCorr[ipar]/expMatCov[ipar]; // consider expectation as measurement
//        mat(parI, parI) += 1. / expMatCov[ipar]; // this measurement is orthogonal to all others
//                                                 //printf("Pnt:%3d MatVar:%d DOF %3d | ExpVal: %+e Cov: %+e\n",ip,ipar,parI, expMatCorr[ipar], expMatCov[ipar]);
//      }
//    } // material effect descripotion params
//    //
//  } // loop over track points
//  //
//  return matp;
//}

//____________________________________________
void Controller::initMPRecOutput()
{
  // prepare MP record output
  mMPRecFile.reset(TFile::Open(fmt::format("{}{}", mMPDatFileName, ".root").c_str(), "recreate"));
  mMPRecTree = std::make_unique<TTree>("mpTree", "MPrecord Tree");
  mMPRecTree->Branch("mprec", "Millepede2Record", &mMPRecordPtr);
  //
}

//____________________________________________
void Controller::initResidOutput()
{
  // prepare residual output
  mResidFile.reset(TFile::Open(mResidFileName.c_str(), "recreate"));
  mResidTree = std::make_unique<TTree>("res", "Control Residuals");
  mResidTree->Branch("t", "ResidualsController", &mCResidPtr);
  //
}

//____________________________________________
void Controller::closeMPRecOutput()
{
  // close output
  if (!mMPRecFile) {
    return;
  }
  LOG(info) << "Closing " << mMPRecFile->GetName();
  mMPRecFile->cd();
  mMPRecTree->Write();
  mMPRecTree.reset();
  mMPRecFile->Close();
  mMPRecFile.reset();
}

//____________________________________________
void Controller::closeResidOutput()
{
  // close output
  if (!mResidFile) {
    return;
  }
  LOG(info) << "Closing " << mResidFile->GetName();
  mResidFile->cd();
  mResidTree->Write();
  mResidTree.reset();
  mResidFile->Close();
  mResidFile.reset();
  mCResid.Clear();
}

//____________________________________________
void Controller::closeMilleOutput()
{
  // close output
  if (mMille) {
    LOG(info) << "Closing " << mMPDatFileName.c_str() << sMPDataExt;
  }
  mMille.reset();
}

//____________________________________________
void Controller::setMPDatFileName(const char* name)
{
  // set output file name
  mMPDatFileName = name;
  if (mMPDatFileName.empty()) {
    mMPDatFileName = "mpData";
  }
  //
}

//____________________________________________
void Controller::setMPParFileName(const char* name)
{
  // set MP params output file name
  mMPParFileName = name;
  if (mMPParFileName.empty()) {
    mMPParFileName = "mpParams.txt";
  }
  //
}

//____________________________________________
void Controller::setMPConFileName(const char* name)
{
  // set MP constraints output file name
  mMPConFileName = name;
  if (mMPConFileName.empty()) {
    mMPConFileName = "mpConstraints.txt";
  }
  //
}

//____________________________________________
void Controller::setMPSteerFileName(const char* name)
{
  // set MP constraints output file name
  mMPSteerFileName = name;
  if (mMPSteerFileName.empty()) {
    mMPSteerFileName = "mpConstraints.txt";
  }
  //
}

//____________________________________________
void Controller::setResidFileName(const char* name)
{
  // set output file name
  mResidFileName = name;
  if (mResidFileName.empty()) {
    mResidFileName = "mpControlRes.root";
  }
  //
}

//____________________________________________
void Controller::setOutCDBPath(const char* name)
{
  // set output storage name
  mOutCDBPath = name;
  if (mOutCDBPath.empty()) {
    mOutCDBPath = "local://outOCDB";
  }
  //
}

//____________________________________________
void Controller::setObligatoryDetector(DetID detID, int trtype, bool v)
{
  // mark detector presence obligatory in the track of given type
  AlignableDetector* det = getDetector(detID);
  if (!det) {
    LOG(error) << "Detector " << detID << " is not defined";
  }
  if (v) {
    mObligatoryDetPattern[trtype] |= detID.getMask();
  } else {
    mObligatoryDetPattern[trtype] &= ~detID.getMask();
  }
  if (det->isObligatory(trtype) != v) {
    det->setObligatory(trtype, v);
  }
  //
}

//____________________________________________
bool Controller::addVertexConstraint()
{
  LOG(fatal) << __PRETTY_FUNCTION__ << " is disabled";
  //FIXME(milettri): needs AliESDtrack
  //  // if vertex is set and if particle is primary, add vertex as a meared point
  //  //
  //  const AliESDtrack* esdTr = fESDTrack[0];
  //  if (!fVertex || !esdTr){
  //    return false;}
  //  //
  //  if (esdTr->GetNcls(0) < mMinITSClforVC){
  //    return false;} // not enough ITS clusters
  //  if (!AlignableDetectorITS::CheckHitPattern(esdTr, mITSPattforVC)){
  //    return false;}
  //  //
  //  AliExternalTrackParam trc = *esdTr;
  //  double dz[2], dzCov[3];
  //  if (!trc.PropagateToDCA(fVertex, AliTrackerBase::GetBz(), 2 * mMaxDCAforVC[0], dz, dzCov)){
  //    return false;}
  //  //
  //  // check if primary candidate
  //  if (Abs(dz[0]) > mMaxDCAforVC[0] || Abs(dz[1]) > mMaxDCAforVC[1]){
  //    return false;}
  //  double covar[6];
  //  fVertex->GetCovMatrix(covar);
  //  double p[2] = {trc.GetParameter()[0] - dz[0], trc.GetParameter()[1] - dz[1]};
  //  double c[3] = {0.5 * (covar[0] + covar[2]), 0., covar[5]};
  //  double chi2 = trc.GetPredictedChi2(p, c);
  //  if (chi2 > mMaxChi2forVC){
  //    return false;}
  //  //
  //  // assing measured vertex rotated to VtxSens frame as reference point
  //  double xyz[3], xyzT[3];
  //  fVertex->GetXYZ(xyz);
  //  mVtxSens->setAlpha(trc.GetAlpha());
  //  // usually translation from GLO to TRA frame should go via matrix T2G
  //  // but for the VertexSensor Local and Global are the same frames
  //  mVtxSens->applyCorrection(xyz);
  //  mVtxSens->getMatrixT2L().MasterToLocal(xyz, xyzT);
  //  mRefPoint->setSensor(mVtxSens);
  //  mRefPoint->setAlphaSens(mVtxSens->getAlpTracking());
  //  mRefPoint->setXYZTracking(xyzT);
  //  mRefPoint->setYZErrTracking(c);
  //  mRefPoint->setContainsMeasurement(true);
  //  mRefPoint->init();
  //  //
  return true;
}

//FIXME(milettri): needs OCDB
////______________________________________________________
//void Controller::writeCalibrationResults() const
//{
//  // writes output calibration
//  CleanOCDB();
//  AliCDBManager::Instance()->SetDefaultStorage(mOutCDBPath.c_str());
//  //
//  AlignableDetector* det;
//  for (auto id=DetID::First; id<=DetID::Last; id++) {
//    if (!(det = getDetector(id)) || det->isDisabled()){
//      continue;
//    }
//    det->writeCalibrationResults();
//  }
//  //
//}

//FIXME(milettri): needs OCDB
////______________________________________________________
//void Controller::SetOutCDBRunRange(int rmin, int rmax)
//{
//  // set output run range
//  mOutCDBRunRange[0] = rmin >= 0 ? rmin : 0;
//  mOutCDBRunRange[1] = rmax > mOutCDBRunRange[0] ? rmax : AliCDBRunRange::Infinity();
//}

//FIXME(milettri): needs OCDB
////______________________________________________________
//bool Controller::LoadRefOCDB()
//{
//  // setup OCDB whose objects will be used as a reference with respect to which the
//  // alignment/calibration will prodice its corrections.
//  // Detectors which need some reference calibration data must use this one
//  //
//  //
//  LOG(info) << "Loading reference OCDB");
//  CleanOCDB();
//  AliCDBManager* man = AliCDBManager::Instance();
//  //
//  if (!mRefOCDBConf.IsNull() && !gSystem->AccessPathName(mRefOCDBConf.c_str(), kFileExists)) {
//    LOG(info) << "Executing reference OCDB setup macro %s", mRefOCDBConf.c_str());
//    if (mRefRunNumber > 0){
//      gROOT->ProcessLine(Form(".x %s(%d)", mRefOCDBConf.c_str(), mRefRunNumber));}
//    else
//      gROOT->ProcessLine(Form(".x %s", mRefOCDBConf.c_str()));
//  } else {
//    LOG(warning) << "No reference OCDB config macro "<<mRefOCDBConf.c_str()<<" is found, assume raw:// with run " << AliCDBRunRange::Infinity();
//    man->SetRaw(true);
//    man->SetRun(AliCDBRunRange::Infinity());
//  }
//  //
//  if (AliGeomManager::GetGeometry()) {
//    LOG(info) << "Destroying current geometry before loading reference one");
//    AliGeomManager::Destroy();
//  }
//  AliGeomManager::LoadGeometry("geometry.root");
//  if (!AliGeomManager::GetGeometry()){
//    LOG(fatal) << "Failed to load geometry, cannot run");}
//  //
//  TString detList = "";
//  for (int i = 0; i < kNDetectors; i++) {
//    detList += getDetNameByDetID(i);
//    detList += " ";
//  }
//  AliGeomManager::ApplyAlignObjsFromCDB(detList.c_str());
//  //
//  mRefOCDBLoaded++;
//  //
//  return true;
//}

//________________________________________________________
AlignableDetector* Controller::getDetOfDOFID(int id) const
{
  // return detector owning DOF with this ID
  for (auto id = DetID::First; id <= DetID::Last; id++) {
    AlignableDetector* det = getDetector(id);
    if (det && det->ownsDOFID(id)) {
      return det;
    }
  }
  return nullptr;
}

//________________________________________________________
AlignableVolume* Controller::getVolOfDOFID(int id) const
{
  // return volume owning DOF with this ID
  for (auto id = DetID::First; id <= DetID::Last; id++) {
    AlignableDetector* det = getDetector(id);
    if (det && det->ownsDOFID(id)) {
      return det->getVolOfDOFID(id);
    }
  }
  if (mVtxSens && mVtxSens->ownsDOFID(id)) {
    return mVtxSens.get();
  }
  return nullptr;
}

//________________________________________________________
void Controller::terminate(bool doStat)
{
  // finalize processing
  if (mRunNumber > 0) {
    fillStatHisto(kRunDone);
  }
  if (doStat) {
    if (mVtxSens) {
      mVtxSens->fillDOFStat(mDOFStat);
    }
  }
  //
  for (auto id = DetID::First; id <= DetID::Last; id++) {
    if (getDetector(id)) {
      getDetector(id)->terminate();
    }
  }
  closeMPRecOutput();
  closeMilleOutput();
  closeResidOutput();
  Print("stat");
  //
}

//________________________________________________________
Char_t* Controller::getDOFLabelTxt(int idf) const
{
  // get DOF full label
  AlignableVolume* vol = getVolOfDOFID(idf);
  if (vol) {
    return Form("%d_%s_%s", getGloParLab(idf), vol->getSymName(),
                vol->getDOFName(idf - vol->getFirstParGloID()));
  }
  //
  // this might be detector-specific calibration dof
  AlignableDetector* det = getDetOfDOFID(idf);
  if (det) {
    return Form("%d_%s_%s", getGloParLab(idf), det->GetName(),
                det->getCalibDOFName(idf - det->getFirstParGloID()));
  }
  return nullptr;
}

//********************* interaction with PEDE **********************

//______________________________________________________
void Controller::genPedeSteerFile(const Option_t* opt) const
{
  // produce steering file template for PEDE + params and constraints
  //
  enum { kOff,
         kOn,
         kOnOn };
  const char* cmt[3] = {"  ", "! ", "!!"};
  const char* kSolMeth[] = {"inversion", "diagonalization", "fullGMRES", "sparseGMRES", "cholesky", "HIP"};
  const int kDefNIter = 3;     // default number of iterations to ask
  const float kDefDelta = 0.1; // def. delta to exit
  TString opts = opt;
  opts.ToLower();
  LOG(info) << "Generating MP2 templates:\n "
            << "Steering   :\t" << mMPSteerFileName << "\n"
            << "Parameters :\t" << mMPParFileName << "\n"
            << "Constraints:\t" << mMPConFileName << "\n";
  //
  FILE* parFl = fopen(mMPParFileName.c_str(), "w+");
  FILE* strFl = fopen(mMPSteerFileName.c_str(), "w+");
  //
  // --- template of steering file
  fprintf(strFl, "%-20s%s %s\n", mMPParFileName.c_str(), cmt[kOnOn], "parameters template");
  fprintf(strFl, "%-20s%s %s\n", mMPConFileName.c_str(), cmt[kOnOn], "constraints template");
  //
  fprintf(strFl, "\n\n%s %s\n", cmt[kOnOn], "MUST uncomment 1 solving methods and tune it");
  //
  int nm = sizeof(kSolMeth) / sizeof(char*);
  for (int i = 0; i < nm; i++) {
    fprintf(strFl, "%s%s %-20s %2d %.2f %s\n", cmt[kOn], "method", kSolMeth[i], kDefNIter, kDefDelta, cmt[kOnOn]);
  }
  //
  const float kDefChi2F0 = 20., kDefChi2F = 3.; // chi2 factors for 1st and following iterations
  const float kDefDWFrac = 0.2;                 // cut outliers with downweighting above this factor
  const int kDefOutlierDW = 4;                  // start Cauchy function downweighting from iteration
  const int kDefEntries = 25;                   // min entries per DOF to allow its variation
  //
  fprintf(strFl, "\n\n%s %s\n", cmt[kOnOn], "Optional settings");
  fprintf(strFl, "\n%s%-20s %.2f %.2f %s %s\n", cmt[kOn], "chisqcut", kDefChi2F0, kDefChi2F,
          cmt[kOnOn], "chi2 cut factors for 1st and next iterations");
  fprintf(strFl, "%s%-20s %2d %s %s\n", cmt[kOn], "outlierdownweighting", kDefOutlierDW,
          cmt[kOnOn], "iteration for outliers downweighting with Cauchi factor");
  fprintf(strFl, "%s%-20s %.3f %s %s\n", cmt[kOn], "dwfractioncut", kDefDWFrac,
          cmt[kOnOn], "cut outliers with downweighting above this factor");
  fprintf(strFl, "%s%-20s %2d %s %s\n", cmt[kOn], "entries", kDefEntries,
          cmt[kOnOn], "min entries per DOF to allow its variation");
  //
  fprintf(strFl, "\n\n\n%s%-20s %s %s\n\n\n", cmt[kOff], "CFiles", cmt[kOnOn], "put below *.mille files list");
  //
  if (mVtxSens) {
    mVtxSens->writePedeInfo(parFl, opt);
  }
  //
  for (auto id = DetID::First; id <= DetID::Last; id++) {
    AlignableDetector* det = getDetector(id);
    if (!det || det->isDisabled()) {
      continue;
    }
    det->writePedeInfo(parFl, opt);
    //
  }
  //
  writePedeConstraints();
  //
  fclose(strFl);
  fclose(parFl);
  //
}

//___________________________________________________________
bool Controller::readParameters(const char* parfile, bool useErrors)
{
  // read parameters file (millepede output)
  if (mNDOFs < 1 || mGloParVal.size() || mGloParErr.size()) {
    LOG(error) << "Something is wrong in init: mNDOFs=" << mNDOFs << " N GloParVal=" << mGloParVal.size() << " N GloParErr=" << mGloParErr.size();
  }
  ifstream inpf(parfile);
  if (!inpf.good()) {
    printf("Failed on input filename %s\n", parfile);
    return false;
  }
  mGloParVal.resize(mNDOFs);
  if (useErrors) {
    mGloParErr.resize(mNDOFs);
  }
  int cnt = 0;
  TString fline;
  fline.ReadLine(inpf);
  fline = fline.Strip(TString::kBoth, ' ');
  fline.ToLower();
  if (!fline.BeginsWith("parameter")) {
    LOG(error) << "First line is not parameter keyword: " << fline.Data();
    return false;
  }
  double v0, v1, v2;
  int lab, asg = 0, asg0 = 0;
  while (fline.ReadLine(inpf)) {
    cnt++;
    fline = fline.Strip(TString::kBoth, ' ');
    if (fline.BeginsWith("!") || fline.BeginsWith("*")) {
      continue;
    } // ignore comment
    int nr = sscanf(fline.Data(), "%d%lf%lf%lf", &lab, &v0, &v1, &v2);
    if (nr < 3) {
      LOG(error) << "Expected to read at least 3 numbers, got " << nr << ", this is NOT milleped output";
      LOG(error) << "line (" << cnt << ") was:\n " << fline.Data();
      return false;
    }
    if (nr == 3) {
      asg0++;
    }
    int parID = label2ParID(lab);
    if (parID < 0 || parID >= mNDOFs) {
      LOG(error) << "Invalid label " << lab << " at line " << cnt << " -> ParID=" << parID;
      return false;
    }
    mGloParVal[parID] = -v0;
    if (useErrors) {
      mGloParErr[parID] = v1;
    }
    asg++;
    //
  };
  LOG(info) << "Read " << cnt << " lines, assigned " << asg << " values, " << asg0 << " dummy";
  //
  return true;
}

//______________________________________________________
void Controller::checkConstraints(const char* params)
{
  // check how the constraints are satisfied with already uploaded or provided params
  //
  if (params && !readParameters(params)) {
    LOG(error) << "Failed to load parameters from " << params;
    return;
  }
  //
  int ncon = getNConstraints();
  for (int icon = 0; icon < ncon; icon++) {
    const GeometricalConstraint* con = getConstraint(icon);
    con->checkConstraint();
  }
  //
}

//___________________________________________________________
void Controller::mPRec2Mille(const char* mprecfile, const char* millefile, bool bindata)
{
  // converts MPRecord tree to millepede binary format
  TFile* flmpr = TFile::Open(mprecfile);
  if (!flmpr) {
    LOG(error) << "Failed to open MPRecord file " << mprecfile;
    return;
  }
  TTree* mprTree = (TTree*)flmpr->Get("mpTree");
  if (!mprTree) {
    LOG(error) << "No mpTree in xMPRecord file " << mprecfile;
    return;
  }
  mPRec2Mille(mprTree, millefile, bindata);
  delete mprTree;
  flmpr->Close();
  delete flmpr;
}

//___________________________________________________________
void Controller::mPRec2Mille(TTree* mprTree, const char* millefile, bool bindata)
{
  // converts MPRecord tree to millepede binary format
  //
  TBranch* br = mprTree->GetBranch("mprec");
  if (!br) {
    LOG(error) << "provided tree does not contain branch mprec";
    return;
  }
  Millepede2Record* rec = new Millepede2Record();
  br->SetAddress(&rec);
  int nent = mprTree->GetEntries();
  TString mlname = millefile;
  if (mlname.IsNull()) {
    mlname = "mpRec2mpData";
  }
  if (!mlname.EndsWith(sMPDataExt)) {
    mlname += sMPDataExt;
  }
  Mille* mille = new Mille(mlname, bindata);
  TArrayF buffDLoc;
  for (int i = 0; i < nent; i++) {
    br->GetEntry(i);
    int nr = rec->getNResid(); // number of residual records
    int nloc = rec->getNVarLoc();
    if (buffDLoc.GetSize() < nloc) {
      buffDLoc.Set(nloc + 100);
    }
    float* buffLocV = buffDLoc.GetArray();
    const float* recDGlo = rec->getArrGlo();
    const float* recDLoc = rec->getArrLoc();
    const short* recLabLoc = rec->getArrLabLoc();
    const int* recLabGlo = rec->getArrLabGlo();
    //
    for (int ir = 0; ir < nr; ir++) {
      memset(buffLocV, 0, nloc * sizeof(float));
      int ndglo = rec->getNDGlo(ir);
      int ndloc = rec->getNDLoc(ir);
      // fill 0-suppressed array from MPRecord to non-0-suppressed array of Mille
      for (int l = ndloc; l--;) {
        buffLocV[recLabLoc[l]] = recDLoc[l];
      }
      //
      mille->mille(nloc, buffLocV, ndglo, recDGlo, recLabGlo, rec->getResid(ir), rec->getResErr(ir));
      //
      recLabGlo += ndglo; // next record
      recDGlo += ndglo;
      recLabLoc += ndloc;
      recDLoc += ndloc;
    }
    mille->end();
  }
  delete mille;
  br->SetAddress(nullptr);
  delete rec;
}

//____________________________________________________________
void Controller::fillStatHisto(int type, float w)
{
  if (!mHistoStat) {
    createStatHisto();
  }
  mHistoStat->Fill((isCosmic() ? kNHVars : 0) + type, w);
}

//____________________________________________________________
void Controller::createStatHisto()
{
  mHistoStat = new TH1F("stat", "stat", 2 * kNHVars, -0.5, 2 * kNHVars - 0.5);
  mHistoStat->SetDirectory(nullptr);
  TAxis* xax = mHistoStat->GetXaxis();
  for (int j = 0; j < 2; j++) {
    for (int i = 0; i < kNHVars; i++) {
      xax->SetBinLabel(j * kNHVars + i + 1, Form("%s.%s", j ? "CSM" : "COL", sHStatName[i]));
    }
  }
}

//____________________________________________________________
void Controller::printLabels() const
{
  // print global IDs and Labels
  for (int i = 0; i < mNDOFs; i++) {
    printf("%5d %s\n", i, getDOFLabelTxt(i));
  }
}

//____________________________________________________________
int Controller::label2ParID(int lab) const
{
  // convert Mille label to ParID (slow)
  int ind = 0; // FIXME RS TODO // findKeyIndex(lab, mOrderedLbl, mNDOFs);
  if (ind < 0) {
    return -1;
  }
  return mLbl2ID[ind];
}

//____________________________________________________________
void Controller::addAutoConstraints()
{
  // add default constraints on children cumulative corrections within the volumes
  for (auto id = DetID::First; id <= DetID::Last; id++) {
    AlignableDetector* det = getDetector(id);
    if (!det || det->isDisabled()) {
      continue;
    }
    det->addAutoConstraints();
  }
  LOG(info) << "Added " << getNConstraints() << " automatic constraints";
}

//____________________________________________________________
void Controller::writePedeConstraints() const
{
  // write constraints file
  FILE* conFl = fopen(mMPConFileName.c_str(), "w+");
  //
  int nconstr = getNConstraints();
  for (int icon = 0; icon < nconstr; icon++) {
    getConstraint(icon)->writeChildrenConstraints(conFl);
  }
  //
  fclose(conFl);
}

//____________________________________________________________
void Controller::fixLowStatFromDOFStat(int thresh)
{
  // fix DOFs having stat below threshold
  //
  if (mNDOFs != mDOFStat.getNDOFs()) {
    LOG(error) << "Discrepancy between NDOFs=" << mNDOFs << " of and statistics object: " << mDOFStat.getNDOFs();
    return;
  }
  for (int parID = 0; parID < mNDOFs; parID++) {
    if (mDOFStat.getStat(parID) >= thresh) {
      continue;
    }
    mGloParErr[parID] = -999.;
  }
  //
}

//____________________________________________________________
void Controller::loadStat(const char* flname)
{
  // load statistics histos from external file produced by alignment task
  TFile* fl = TFile::Open(flname);
  //
  TObject *hdfO = nullptr, *hstO = nullptr;
  TList* lst = (TList*)fl->Get("clist");
  if (lst) {
    hdfO = lst->FindObject("DOFstat");
    if (hdfO) {
      lst->Remove(hdfO);
    }
    hstO = lst->FindObject("stat");
    if (hstO) {
      lst->Remove(hstO);
    }
    delete lst;
  } else {
    hdfO = fl->Get("DOFstat");
    hstO = fl->Get("stat");
  }
  TH1F* hst = nullptr;
  if (hstO && (hst = dynamic_cast<TH1F*>(hstO))) {
    hst->SetDirectory(nullptr);
  } else {
    LOG(warning) << "did not find stat histo";
  }
  //
  DOFStatistics* dofSt = nullptr;
  if (!hdfO || !(dofSt = dynamic_cast<DOFStatistics*>(hdfO))) {
    LOG(warning) << "did not find DOFstat object";
  }
  //
  setHistoStat(hst);
  setDOFStat(*dofSt); // FIXME RS TODO
  //
  fl->Close();
  delete fl;
}

//______________________________________________
void Controller::checkSol(TTree* mpRecTree, bool store,
                          bool verbose, bool loc, const char* outName)
{
  // do fast check of pede solution with MPRecord tree
  ResidualsControllerFast* rLG = store ? new ResidualsControllerFast() : nullptr;
  ResidualsControllerFast* rL = store && loc ? new ResidualsControllerFast() : nullptr;
  TTree *trLG = nullptr, *trL = nullptr;
  TFile* outFile = nullptr;
  if (store) {
    TString outNS = outName;
    if (outNS.IsNull()) {
      outNS = "resFast";
    }
    if (!outNS.EndsWith(".root")) {
      outNS += ".root";
    }
    outFile = TFile::Open(outNS.Data(), "recreate");
    trLG = new TTree("resFLG", "Fast residuals with LG correction");
    trLG->Branch("rLG", "ResidualsControllerFast", &rLG);
    //
    if (rL) {
      trL = new TTree("resFL", "Fast residuals with L correction");
      trL->Branch("rL", "ResidualsControllerFast", &rL);
    }
  }
  //
  Millepede2Record* rec = new Millepede2Record();
  mpRecTree->SetBranchAddress("mprec", &rec);
  int nrec = mpRecTree->GetEntriesFast();
  for (int irec = 0; irec < nrec; irec++) {
    mpRecTree->GetEntry(irec);
    checkSol(rec, rLG, rL, verbose, loc);
    // store even in case of failure, to have the trees aligned with controlRes
    if (trLG) {
      trLG->Fill();
    }
    if (trL) {
      trL->Fill();
    }
  }
  //
  // save
  if (trLG) {
    outFile->cd();
    trLG->Write();
    delete trLG;
    if (trL) {
      trL->Write();
      delete trL;
    }
    outFile->Close();
    delete outFile;
  }
  //
}

//______________________________________________
bool Controller::checkSol(Millepede2Record* rec,
                          ResidualsControllerFast* rLG, ResidualsControllerFast* rL,
                          bool verbose, bool loc)
{
  LOG(fatal) << __PRETTY_FUNCTION__ << " is disabled";
  //FIXME(milettri): needs AliSymMatrix
  //  // Check pede solution using derivates, rather than updated geometry
  //  // If loc==true, also produces residuals for current geometry,
  //  // neglecting global corrections
  //  //
  //  if (rL){
  //    loc = true;} // if local sol. tree asked, always evaluate it
  //  //
  //  int nres = rec->getNResid();
  //  //
  //  const float* recDGlo = rec->getArrGlo();
  //  const float* recDLoc = rec->getArrLoc();
  //  const short* recLabLoc = rec->getArrLabLoc();
  //  const int* recLabGlo = rec->getArrLabGlo();
  //  int nvloc = rec->getNVarLoc();
  //  //
  //  // count number of real measurement duplets and material correction fake 4-plets
  //  int nPoints = 0;
  //  int nMatCorr = 0;
  //  for (int irs = 0; irs < nres; irs++) {
  //    if (rec->getNDGlo(irs) > 0) {
  //      if (irs == nres - 1 || rec->getNDGlo(irs + 1) == 0){
  //        LOG(fatal) << ("Real coordinate measurements must come in pairs");}
  //      nPoints++;
  //      irs++; // skip 2nd
  //      continue;
  //    } else if (rec->getResid(irs) == 0 && rec->getVolID(irs) == -1) { // material corrections have 0 residual
  //      nMatCorr++;
  //    } else { // might be fixed parameter, global derivs are skept
  //      nPoints++;
  //      irs++; // skip 2nd
  //      continue;
  //    }
  //  }
  //  //
  //  if (nMatCorr % 4){
  //    LOG(warning) << "Error? NMatCorr=" << nMatCorr << " is not multiple of 4";}
  //  //
  //  if (rLG) {
  //    rLG->Clear();
  //    rLG->setNPoints(nPoints);
  //    rLG->setNMatSol(nMatCorr);
  //    rLG->setCosmic(rec->isCosmic());
  //  }
  //  if (rL) {
  //    rL->Clear();
  //    rL->setNPoints(nPoints);
  //    rL->setNMatSol(nMatCorr);
  //    rL->setCosmic(rec->isCosmic());
  //  }
  //  //
  //  AliSymMatrix* matpG = new AliSymMatrix(nvloc);
  //  TVectorD *vecp = 0, *vecpG = new TVectorD(nvloc);
  //  //
  //  if (loc){
  //    vecp = new TVectorD(nvloc);}
  //  //
  //  float chi2Ini = 0, chi2L = 0, chi2LG = 0;
  //  //
  //  // residuals, accounting for global solution
  //  double* resid = new double[nres];
  //  int* volID = new int[nres];
  //  for (int irs = 0; irs < nres; irs++) {
  //    double resOr = rec->getResid(irs);
  //    resid[irs] = resOr;
  //    //
  //    int ndglo = rec->getNDGlo(irs);
  //    int ndloc = rec->getNDLoc(irs);
  //    volID[irs] = 0;
  //    for (int ig = 0; ig < ndglo; ig++) {
  //      int lbI = recLabGlo[ig];
  //      int idP = label2ParID(lbI);
  //      if (idP < 0){
  //        LOG(fatal) << "Did not find parameted for label " << lbI;}
  //      double parVal = getGloParVal()[idP];
  //      //      resid[irs] -= parVal*recDGlo[ig];
  //      resid[irs] += parVal * recDGlo[ig];
  //      if (!ig) {
  //        AlignableVolume* vol = getVolOfDOFID(idP);
  //        if (vol){
  //          volID[irs] = vol->getVolID();}
  //        else
  //          volID[irs] = -2; // calibration DOF !!! TODO
  //      }
  //    }
  //    //
  //    double sg2inv = rec->getResErr(irs);
  //    sg2inv = 1. / (sg2inv * sg2inv);
  //    //
  //    chi2Ini += resid[irs] * resid[irs] * sg2inv; // chi accounting for global solution only
  //    //
  //    // Build matrix to solve local parameters
  //    for (int il = 0; il < ndloc; il++) {
  //      int lbLI = recLabLoc[il]; // id of local variable
  //      (*vecpG)[lbLI] -= recDLoc[il] * resid[irs] * sg2inv;
  //      if (loc){
  //        (*vecp)[lbLI] -= recDLoc[il] * resOr * sg2inv;}
  //      for (int jl = il + 1; jl--;) {
  //        int lbLJ = recLabLoc[jl]; // id of local variable
  //        (*matpG)(lbLI, lbLJ) += recDLoc[il] * recDLoc[jl] * sg2inv;
  //      }
  //    }
  //    //
  //    recLabGlo += ndglo; // prepare for next record
  //    recDGlo += ndglo;
  //    recLabLoc += ndloc;
  //    recDLoc += ndloc;
  //    //
  //  }
  //  //
  //  if (rL){
  //    rL->setChi2Ini(chi2Ini);}
  //  if (rLG){
  //    rLG->setChi2Ini(chi2Ini);}
  //  //
  //  TVectorD vecSol(nvloc);
  //  TVectorD vecSolG(nvloc);
  //  //
  //  if (!matpG->SolveChol(*vecpG, vecSolG, false)) {
  //    LOG(info) << "Failed to solve track corrected for globals";
  //    delete matpG;
  //    matpG = 0;
  //  } else if (loc) { // solution with local correction only
  //    if (!matpG->SolveChol(*vecp, vecSol, false)) {
  //      LOG(info) << "Failed to solve track corrected for globals";
  //      delete matpG;
  //      matpG = 0;
  //    }
  //  }
  //  delete vecpG;
  //  delete vecp;
  //  if (!matpG) { // failed
  //    delete[] resid;
  //    delete[] volID;
  //    if (rLG){
  //      rLG->Clear();}
  //    if (rL){
  //      rL->Clear();}
  //    return false;
  //  }
  //  // check
  //  recDGlo = rec->getArrGlo();
  //  recDLoc = rec->getArrLoc();
  //  recLabLoc = rec->getArrLabLoc();
  //  recLabGlo = rec->getArrLabGlo();
  //  //
  //  if (verbose) {
  //    printf(loc ? "Sol L/LG:\n" : "Sol LG:\n");
  //    int nExtP = (nvloc % 4) ? 5 : 4;
  //    for (int i = 0; i < nExtP; i++){
  //      loc ? printf("%+.3e/%+.3e ", vecSol[i], vecSolG[i]) : printf("%+.3e ", vecSolG[i]);}
  //    printf("\n");
  //    bool nln = true;
  //    int cntL = 0;
  //    for (int i = nExtP; i < nvloc; i++) {
  //      nln = true;
  //      loc ? printf("%+.3e/%+.3e ", vecSol[i], vecSolG[i]) : printf("%+.3e ", vecSolG[i]);
  //      if (((++cntL) % 4) == 0) {
  //        printf("\n");
  //        nln = false;
  //      }
  //    }
  //    if (!nln){
  //      printf("\n");}
  //    if (loc){
  //      printf("%3s (%9s) %6s | [ %7s:%7s ] [ %7s:%7s ]\n", "Pnt", "Label",
  //             "Sigma", "resid", "pull/L ", "resid", "pull/LG");}
  //    else{
  //      printf("%3s (%9s) %6s | [ %7s:%7s ]\n", "Pnt", "Label",
  //             "Sigma", "resid", "pull/LG");}
  //  }
  //  int idMeas = -1, pntID = -1, matID = -1;
  //  for (int irs = 0; irs < nres; irs++) {
  //    double resOr = rec->getResid(irs);
  //    double resL = resOr;
  //    double resLG = resid[irs];
  //    double sg = rec->getResErr(irs);
  //    double sg2Inv = 1 / (sg * sg);
  //    //
  //    int ndglo = rec->getNDGlo(irs);
  //    int ndloc = rec->getNDLoc(irs);
  //    //
  //    for (int il = 0; il < ndloc; il++) {
  //      int lbLI = recLabLoc[il]; // id of local variable
  //      resL += recDLoc[il] * vecSol[lbLI];
  //      resLG += recDLoc[il] * vecSolG[lbLI];
  //    }
  //    //
  //    chi2L += resL * resL * sg2Inv;    // chi accounting for global solution only
  //    chi2LG += resLG * resLG * sg2Inv; // chi accounting for global solution only
  //    //
  //    if (ndglo || resOr != 0) { // real measurement
  //      idMeas++;
  //      if (idMeas > 1){
  //        idMeas = 0;}
  //      if (idMeas == 0){
  //        pntID++;} // measurements come in pairs
  //      int lbl = rec->getVolID(irs);
  //      lbl = ndglo ? recLabGlo[0] : 0; // TMP, until VolID is filled // RS!!!!
  //      if (rLG) {
  //        rLG->setResSigMeas(pntID, idMeas, resLG, sg);
  //        if (idMeas == 0){
  //          rLG->setLabel(pntID, lbl, volID[irs]);}
  //      }
  //      if (rL) {
  //        rL->setResSigMeas(pntID, idMeas, resL, sg);
  //        if (idMeas == 0){
  //          rL->setLabel(pntID, lbl, volID[irs]);}
  //      }
  //    } else {
  //      matID++; // mat.correcitons come in 4-plets, but we fill each separately
  //      //
  //      if (rLG){
  //        rLG->setMatCorr(matID, resLG, sg);}
  //      if (rL){
  //        rL->setMatCorr(matID, resL, sg);}
  //    }
  //    //
  //    if (verbose) {
  //      int lbl = rec->getVolID(irs);
  //      lbl = ndglo ? recLabGlo[0] : (resOr == 0 ? -1 : 0); // TMP, until VolID is filled // RS!!!!
  //      if (loc){
  //        printf("%3d (%9d) %6.4f | [%+.2e:%+7.2f] [%+.2e:%+7.2f]\n",
  //               irs, lbl, sg, resL, resL / sg, resLG, resLG / sg);}
  //      else
  //        printf("%3d (%9d) %6.4f | [%+.2e:%+7.2f]\n",
  //               irs, lbl, sg, resLG, resLG / sg);
  //    }
  //    //
  //    recLabGlo += ndglo; // prepare for next record
  //    recDGlo += ndglo;
  //    recLabLoc += ndloc;
  //    recDLoc += ndloc;
  //  }
  //  if (rL){
  //    rL->setChi2(chi2L);}
  //  if (rLG){
  //    rLG->setChi2(chi2LG);}
  //  //
  //  if (verbose) {
  //    printf("Chi: G = %e | LG = %e", chi2Ini, chi2LG);
  //    if (loc){
  //      printf(" | L = %e", chi2L);}
  //    printf("\n");
  //  }
  //  // store track corrections
  //  int nTrCor = nvloc - matID - 1;
  //  for (int i = 0; i < nTrCor; i++) {
  //    if (rLG){
  //      rLG->getTrCor()[i] = vecSolG[i];}
  //    if (rL){
  //      rL->getTrCor()[i] = vecSol[i];}
  //  }
  //  //
  //  delete[] resid;
  //  delete[] volID;
  return true;
}

//______________________________________________
void Controller::applyAlignmentFromMPSol()
{
  // apply alignment from millepede solution array to reference alignment level
  LOG(info) << "Applying alignment from Millepede solution";
  for (auto id = DetID::First; id <= DetID::Last; id++) {
    AlignableDetector* det = getDetector(id);
    if (!det || det->isDisabled()) {
      continue;
    }
    det->applyAlignmentFromMPSol();
  }
  setMPAlignDone();
  //
}

//______________________________________________
void Controller::expandGlobalsBy(int n)
{
  // expand global param contaiers by n
  int snew = n + mGloParVal.size();
  mGloParVal.resize(snew);
  mGloParErr.resize(snew);
  mGloParLab.resize(snew);
  mOrderedLbl.resize(snew);
  mLbl2ID.resize(snew);
  mNDOFs += n;
}

} // namespace align
} // namespace o2
