#include <limits>
#include <cmath>

#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "StarClassLibrary/SystemOfUnits.h"

#include "StMixerTriplet.h"
#include "StPicoHFMaker/StHFTriplet.h"
#include "StMixerClosePair.h"
#include "StMixerTrack.h"

ClassImp(StMixerTriplet)

// _________________________________________________________
StMixerTriplet::StMixerTriplet(StMixerClosePair &  pair, const StMixerTrack & particle3,
			       float p3MassHypo,
			       const StThreeVectorF & vtx1, const StThreeVectorF & vtx3,
			       float bField) :
  StHFTriplet()
{
  // -- Create pair out of 2 tracks
  //     prefixes code:
  //      p1 means particle 1
  //      p2 means particle 2
  //      pair means particle1-particle2  pair

  StThreeVectorF dVtx13 = vtx1 - vtx3;

  StPhysicalHelixD p3Helix(particle3.gMom(), particle3.origin() + dVtx13, bField*kilogauss, particle3.charge());

  calculateTopology(&pair, p3Helix, p3MassHypo, particle3.charge(), 0, vtx1, bField);
}
