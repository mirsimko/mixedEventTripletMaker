#ifndef StMixerTriplet_hh
#define StMixerTriplet_hh

/* **************************************************
 *  Generic class calculating and storing triplets in Event Mixing wrapped arround the StHFTriplet class
 *  Allows to combine:
 *  - three particles, using
 *      StMixerTriplet(StPicoTrack const * particle1, StPicoTrack const * particle2, ...
 *
 *  Make sure that this does not happen (no virtual destructiors in the base class):
 *  StMixerTriplet *triplet = new StMixerTriplet(...);
 *  StHFTriplet *tr = triplet;
 *  delete tr;
 *
 *
 * **************************************************
 *
 *  Initial Authors: 
 *            Michael Lomnitz (mrlomnitz@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
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
#include "StPicoHFMaker/StHFTriplet.h"
#include "StMixerClosePair.h"

class StMixerTrack;
class StHFTriplet;
class StMixerClosePair;

class StMixerTriplet : public StHFTriplet
{
 public:
  StMixerTriplet(StMixerClosePair & pair, StMixerTrack const & particle3,
		 float p3MassHypo,
		 StThreeVectorF const & vtx1, StThreeVectorF const & vtx3,
		 float bField);

  ~StMixerTriplet() {;}

 private:
  StMixerTriplet(const StMixerTriplet &);
  StMixerTriplet& operator=(StMixerTriplet const &);
  ClassDef(StMixerTriplet,1)
};
#endif
