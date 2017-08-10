#ifndef StMixerClosePair_hh
#define StMixerClosePair_hh

/* **************************************************
 *  ClosePairClass for storing pairs for a future calculation of a triplet.
 *  Some topological cuts can be performed on the pair to save computaion time.
 *
 *  Make sure that this does not happen (no virtual destructiors in the base class):
 *  StMixerClosePair *ClosePair = new StMixerClosePair(...);
 *  StHFClosePair *pair = ClosePair;
 *  delete pair;
 *
 *
 * **************************************************
 *
 *  Initial Authors: 
 *          **Miroslav Simko  (msimko@bnl.gov)
 *
 *  ** Code Maintainer 
 *
 * **************************************************
 */

#include "TObject.h"
#include "StarClassLibrary/StLorentzVectorF.hh"
#include "StarClassLibrary/StThreeVectorF.hh"
#include "StMixerTrack.h"
#include "StPicoHFMaker/StHFClosePair.h"

class StMixerTrack;
class StHFClosePair;

class StMixerClosePair : public StHFClosePair
{
 public:
  StMixerClosePair(StMixerTrack const &  particle1, StMixerTrack const & particle2,
		   float p1MassHypo, float p2MassHypo,
		   const StThreeVectorF & vtx1, const StThreeVectorF & vtx2,
		   float bField);

  ~StMixerClosePair() {;}
  
 private:
  StMixerClosePair(StMixerClosePair const &);
  StMixerClosePair& operator=(StMixerClosePair const &);
  ClassDef(StMixerClosePair,1)
};
#endif
