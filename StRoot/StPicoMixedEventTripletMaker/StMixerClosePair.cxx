#include <limits>
#include <cmath>

#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "StarClassLibrary/SystemOfUnits.h"

#include "StMixerClosePair.h"
#include "StPicoHFMaker/StHFClosePair.h"
#include "StMixerTrack.h"

ClassImp(StMixerClosePair)

// _________________________________________________________
StMixerClosePair::StMixerClosePair(const StMixerTrack &  particle1, const StMixerTrack & particle2,
			       float p1MassHypo, float p2MassHypo,
			       const StThreeVectorF & vtx1, const StThreeVectorF & vtx2,
			       float bField) :
  StHFClosePair()
  // mParticle1Mom(StThreeVectorF()),
  // mParticle2Mom(StThreeVectorF())
{
  // -- Create pair out of 2 tracks
  //     prefixes code:
  //      p1 means particle 1
  //      p2 means particle 2
  //      pair means particle1-particle2  pair

  if (vtx1 == vtx2) // for the same event
  {
    // see if the particles are the same
    if( particle1.gMom() == particle2.gMom()) { 
      // if they have exactly the same momentum, it means that they are the same tracks
      return;
    }
  }

  StThreeVectorF dVtx12 = vtx1 - vtx2;
  StPhysicalHelixD *p1Helix = new StPhysicalHelixD( particle1.gMom(), particle1.origin(),bField*kilogauss, particle1.charge());
  StPhysicalHelixD *p2Helix = new StPhysicalHelixD( particle2.gMom(), particle2.origin() + dVtx12, bField*kilogauss,  particle2.charge());
  if (!p1Helix || !p2Helix)
  {
    cerr << "StMixerClosePair::StMixerClosePair(...): Helices not initiated" << endl;
    return;
  }

  calculateTopology(p1Helix, p2Helix, p1MassHypo, p2MassHypo, particle1.charge(), particle2.charge(), 0,0, vtx1, bField, true);
}
