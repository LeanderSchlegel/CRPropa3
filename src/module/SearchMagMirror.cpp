#include "crpropa/module/SearchMagMirror.h"
#include "crpropa/Vector3.h"


namespace crpropa{

SearchMagMirror::SearchMagMirror(double tolerance, int timescale) : tolerance(tolerance), timescale(timescale){
}

void SearchMagMirror::addOutput(Module *TxtOut) {
	this->TxtOut = TxtOut;
}

void SearchMagMirror::process(Candidate *c) const {

	// check if candidate has already enough positions in the vector
	std::vector<double> lastField = c -> lastMagneticField;
	std::vector<double> muList = c -> lastPitchAngle;

	Vector3d Bcurrent = c->getBField();
	double B = Bcurrent.getR();
	Vector3d Pcurrent = c->current.getDirection();
	double mu = Pcurrent.dot(Bcurrent) / Pcurrent.getR() / B;

	if(lastField.size() < timescale) {
		c -> lastMagneticField.push_back(B); // update list with new position
		c -> lastPitchAngle.push_back(mu);
		return;
	}

	// check reference values
	double B_ref = lastField[0];
	double mu_ref = muList[0];

	// check for sign change in z-movement
	if(Pcurrent.getZ() * c -> previous.getDirection().getZ() < 0) {
		// check criterium for mirroring 
		double criteria = (1 - mu_ref * mu_ref) / (B_ref / B) - 1;
		if(fabs(criteria) < tolerance)
			TxtOut -> process(c);
	} 

	// update reference momentum and position list	
	lastField.erase(lastField.begin());
	lastField.push_back(B);
	c -> lastMagneticField = lastField;

	muList.erase(muList.begin());
	muList.push_back(mu);
	c -> lastPitchAngle = muList;
}


} // namespace crpropa

