#include "crpropa/module/SearchMagMirror.h"
#include "crpropa/Vector3.h"

#include <unistd.h>

namespace crpropa{

SearchMagMirror::SearchMagMirror(double tolerance, int timescale, bool isMaxTimescale) : tolerance(tolerance), timescale(timescale), isMaxTimescale(isMaxTimescale) {
}

void SearchMagMirror::addOutput(ref_ptr<Module> TxtOut) {
	this->TxtOut = TxtOut;
}

void SearchMagMirror::process(Candidate *c) const {

	// check if candidate has already enough positions in the vector
	double zcurrent = (c->current.getPosition()).getZ();
	std::vector<double> positionList = c -> lastPositions;
	if(positionList.size() < timescale) {
		if(isMaxTimescale) {	// only one module in the simulation is allowed do add to the last positions
			//std::cout << "add to vector \n";
			//sleep(.1);
			c-> lastPositions.push_back(zcurrent); // update list with new position
		}
		return;
	}

	//get current par momentum, reference momentum and counter values.
	Vector3d Pcurrent = c->current.getDirection();
	Vector3d Bcurrent = c->getBField();
	double pparcurrent = Pcurrent.dot(Bcurrent) / Bcurrent.getR();
	double Pref = c -> lastParallelMomentum;
	//std::cout << "Pref\n";
	//sleep(0.1);
	double zref = positionList[0];
	
	// check for sign change in parallel momentum
	if(Pref * pparcurrent < 0) {
		// check for significant change in z position
		if (fabs(zcurrent - zref) > tolerance) {
			//std::cout << "output! \n";
			//sleep(.1);
			TxtOut -> process(c);
		}
	} 

	// update reference momentum and position list	(only for one module in the simulation)
	if(isMaxTimescale) {
		//std::cout << "start with erasing \n";
		//sleep(.1);
		positionList.push_back(zcurrent);
		//std::cout << "push_back done \n";
		//sleep(.1);
		positionList.erase(positionList.begin());
		//std::cout << "erase done \n";
		//sleep(.1);
		c -> lastParallelMomentum = pparcurrent;
		//std::cout << "ppar done \n";
		//sleep(.1);
	}
}


} // namespace crpropa

