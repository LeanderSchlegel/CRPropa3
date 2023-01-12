#ifndef SEARCHMAGMIRROR_H
#define SEARCHMAGMIRROR_H

#include "crpropa/Module.h"
#include "crpropa/Units.h"

namespace crpropa {

/**
 @class SearchMagMirror
 @brief Runtime analysis/reduction searching for magnetic mirroring events.

 This module implements a runtime analysis for magnetic mirroring events.
 The parallel momentum of a candidate is checked for a flip in sign
 relative to a reference position and eventually reports a mirrorposition in the output.
 */
 
class SearchMagMirror: public Module {
private:
	double tolerance; 
	int timescale;
	bool isMaxTimescale;
	ref_ptr<Module> TxtOut; 

public:
	SearchMagMirror(double tolerance = (0.0), int timescale = (0), bool isMaxTimescale = true);
	void process(Candidate *c) const;
	void addOutput(ref_ptr<Module> TxtOut);
};

} // namespace crpropa

#endif // SEARCHMAGMIRROR_H

