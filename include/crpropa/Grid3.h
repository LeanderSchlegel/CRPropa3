#ifndef CRPROPA_GRID_H
#define CRPROPA_GRID_H

#include "crpropa/Referenced.h"
#include "crpropa/Vector3.h"
#include <vector>
#include <typeinfo>
//#include <fstream>

namespace crpropa {

/** Lower and upper neighbor in a periodically continued unit grid */
inline void periodicClamp(double x, int n, int &lo, int &hi) {
	lo = ((int(floor(x)) % n) + n) % n;
	hi = (lo + 1) % n;
}

/** Lower and upper neighbor in a reflectively repeated unit grid */
inline void reflectiveClamp(double x, int n, int &lo, int &hi) {
	while ((x < 0) or (x > n))
		x = 2 * n * (x > n) - x;
	lo = floor(x);
	hi = lo + (lo < n-1);
}

/** Symmetrical round */
inline double round(double r) {
    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}

/**
 @class Grid
 @brief Template class for fields on a periodic grid with trilinear interpolation

 The grid spacing is constant and equal along all three axes.
 Values are calculated by trilinear interpolation of the surrounding 8 grid points.
 The grid is periodically (default) or reflectively extended.
 The grid sample positions are at 1/2 * size/N, 3/2 * size/N ... (2N-1)/2 * size/N.
 */
template<typename T>
class Grid: public Referenced {
	std::vector<T> grid;
	size_t Nx, Ny, Nz; /**< Number of grid points */
	Vector3d origin; /**< Origin of the volume that is represented by the grid. */
	Vector3d gridOrigin; /**< Grid origin */
	double spacing; /**< Distance between grid points, determines the extension of the grid */
	bool reflective; /**< If set to true, the grid is repeated reflectively instead of periodically */

public:
	/** Constructor for cubic grid
	 @param	origin	Position of the lower left front corner of the volume
	 @param	N		Number of grid points in one direction
	 @param spacing	Spacing between grid points
	 */
	Grid(Vector3d origin, size_t N, double spacing) {
		setOrigin(origin);
		setGridSize(N, N, N);
		setSpacing(spacing);
		setReflective(false);
	}

	/** Constructor for non-cubic grid
	 @param	origin	Position of the lower left front corner of the volume
	 @param	Nx		Number of grid points in x-direction
	 @param	Ny		Number of grid points in y-direction
	 @param	Nz		Number of grid points in z-direction
	 @param spacing	Spacing between grid points
	 */
	Grid(Vector3d origin, size_t Nx, size_t Ny, size_t Nz, double spacing) {
		setOrigin(origin);
		setGridSize(Nx, Ny, Nz);
		setSpacing(spacing);
		setReflective(false);
	}

	void setOrigin(Vector3d origin) {
		this->origin = origin;
		this->gridOrigin = origin + Vector3d(spacing/2);
	}

	/** Resize grid, also enlarges the volume as the spacing stays constant */
	void setGridSize(size_t Nx, size_t Ny, size_t Nz) {
		this->Nx = Nx;
		this->Ny = Ny;
		this->Nz = Nz;
		grid.resize(Nx * Ny * Nz);
		setOrigin(origin);
	}

	void setSpacing(double spacing) {
		this->spacing = spacing;
		setOrigin(origin);
	}

	void setReflective(bool b) {
		reflective = b;
	}

	Vector3d getOrigin() const {
		return origin;
	}
	size_t getNx() const {
		return Nx;
	}

	size_t getNy() const {
		return Ny;
	}

	size_t getNz() const {
		return Nz;
	}

	double getSpacing() const {
		return spacing;
	}

	bool isReflective() const {
		return reflective;
	}

	/** Inspector & Mutator */
	T &get(size_t ix, size_t iy, size_t iz) {
		return grid[ix * Ny * Nz + iy * Nz + iz];
	}
	
	//Set Grid Value
	//~ T &setval(double x, double y,double z) {
		//~ grid.get(0, 0, 0) = Vector3f(x, y, z);
	//~ }

	/** Inspector */
	const T &get(size_t ix, size_t iy, size_t iz) const {
		return grid[ix * Ny * Nz + iy * Nz + iz];
	}

	T getValue(size_t ix, size_t iy, size_t iz) {
		return grid[ix * Ny * Nz + iy * Nz + iz];
	}

	/** Return a reference to the grid values */
	std::vector<T> &getGrid() {
		return grid;
	}

	/** Position of the grid point of a given index */
	Vector3d positionFromIndex(int index) const {
		int ix = index / (Ny * Nz);
		int iy = (index / Nz) % Ny;
		int iz = index % Nz;
		return Vector3d(ix, iy, iz) * spacing + gridOrigin;
	}

	/** Value of a grid point that is closest to a given position */
	T closestValue(const Vector3d &position) const {
		Vector3d r = (position - gridOrigin) / spacing;
		int ix = round(r.x);
		int iy = round(r.y);
		int iz = round(r.z);
		if (reflective) {
			while ((ix < 0) or (ix > Nx))
				ix = 2 * Nx * (ix > Nx) - ix;
			while ((iy < 0) or (iy > Ny))
				iy = 2 * Ny * (iy > Ny) - iy;
			while ((iz < 0) or (iz > Nz))
				iz = 2 * Nz * (iz > Nz) - iz;
		} else {
			ix = ((ix % Nx) + Nx) % Nx;
			iy = ((iy % Ny) + Ny) % Ny;
			iz = ((iz % Nz) + Nz) % Nz;
		}
		return get(ix, iy, iz);
	}


/** Interpolate the grid at a given position */ //OLD FUNCTION FOR COMPARISON
	T interpolate_old(const Vector3d &position) const {
		// position on a unit grid
		Vector3d r = (position - gridOrigin) / spacing;

		// indices of lower and upper neighbors
		int ix, iX, iy, iY, iz, iZ;
		if (reflective) {
			reflectiveClamp(r.x, Nx, ix, iX);
			reflectiveClamp(r.y, Ny, iy, iY);
			reflectiveClamp(r.z, Nz, iz, iZ);
		} else {
			periodicClamp(r.x, Nx, ix, iX);
			periodicClamp(r.y, Ny, iy, iY);
			periodicClamp(r.z, Nz, iz, iZ);
		}

		// linear fraction to lower and upper neighbors
		double fx = r.x - floor(r.x);
		double fX = 1 - fx;
		double fy = r.y - floor(r.y);
		double fY = 1 - fy;
		double fz = r.z - floor(r.z);
		double fZ = 1 - fz;

		// trilinear interpolation (see http://paulbourke.net/miscellaneous/interpolation)
		T b(0.);
		//V000 (1 - x) (1 - y) (1 - z) +
		b += get(ix, iy, iz) * fX * fY * fZ;
		//V100 x (1 - y) (1 - z) +
		b += get(iX, iy, iz) * fx * fY * fZ;
		//V010 (1 - x) y (1 - z) +
		b += get(ix, iY, iz) * fX * fy * fZ;
		//V001 (1 - x) (1 - y) z +
		b += get(ix, iy, iZ) * fX * fY * fz;
		//V101 x (1 - y) z +
		b += get(iX, iy, iZ) * fx * fY * fz;
		//V011 (1 - x) y z +
		b += get(ix, iY, iZ) * fX * fy * fz;
		//V110 x y (1 - z) +
		b += get(iX, iY, iz) * fx * fy * fZ;
		//V111 x y z
		b += get(iX, iY, iZ) * fx * fy * fz;

		return b;
	}

T interpolate(const Vector3d &position) {
return interpolate(T(),position);
}

private:

Vector3d interpolate(Vector3d, const Vector3d &position) {
		//position on a unit grid
		Vector3d r = (position - gridOrigin) / spacing;
		
		//std::cout << "Vector3d interpolate(Vector3d,const Vector3d &position) { //VEKTOR "<< std::endl;
		//std::cout << "V"<< std::endl;
		
		//FILE * myfile("testlog.txt","a");
		//myfile << "Vector3d interpolate(Vector3d,const Vector3d &position) { //VEKTOR "<< std::endl;
		
		//~ std::cout << flush();
// indices of lower and upper neighbors
		int ix, iX, iy, iY, iz, iZ;
		if (reflective) {
			reflectiveClamp(r.x, Nx, ix, iX);
			reflectiveClamp(r.y, Ny, iy, iY);
			reflectiveClamp(r.z, Nz, iz, iZ);
		} else {
			periodicClamp(r.x, Nx, ix, iX);
			periodicClamp(r.y, Ny, iy, iY);
			periodicClamp(r.z, Nz, iz, iZ);
		}

		// linear fraction to lower and upper neighbors
		double fx = r.x - floor(r.x);
		double fX = 1 - fx;
		double fy = r.y - floor(r.y);
		double fY = 1 - fy;
		double fz = r.z - floor(r.z);
		double fZ = 1 - fz;
		
		//Vektor wie vorher bestimmen
		Vector3d v = Vector3d(0,0,0);
		// trilinear interpolation (see http://paulbourke.net/miscellaneous/interpolation)
		//T v(0.);
		//V000 (1 - x) (1 - y) (1 - z) +
		v += get(ix, iy, iz) * fX * fY * fZ;
		//V100 x (1 - y) (1 - z) +
		v += get(iX, iy, iz) * fx * fY * fZ;
		//V010 (1 - x) y (1 - z) +
		v += get(ix, iY, iz) * fX * fy * fZ;
		//V001 (1 - x) (1 - y) z +
		v += get(ix, iy, iZ) * fX * fY * fz;
		//V101 x (1 - y) z +
		v += get(iX, iy, iZ) * fx * fY * fz;
		//V011 (1 - x) y z +
		v += get(ix, iY, iZ) * fX * fy * fz;
		//V110 x y (1 - z) +
		v += get(iX, iY, iz) * fx * fy * fZ;
		//V111 x y z
		v += get(iX, iY, iZ) * fx * fy * fz;

		//return b;
		
		//Betrag anpassen durch neuen Algorithmus
		
		//Trilinear interpolation - improved algorithm
		
			
		
			Vector3d v000 = get(ix,iy,iz);
			Vector3d v100 = get(iX,iy,iz);
			Vector3d v010 = get(ix,iY,iz);
			Vector3d v001 = get(ix,iy,iZ);
			Vector3d v101 = get(iX,iy,iZ);
			Vector3d v011 = get(ix,iY,iZ);
			Vector3d v110 = get(iX,iY,iz);
			Vector3d v111 = get(iX,iY,iZ);
			
			double laenge = (1/pow((fX*fY*fZ),11/3)+pow(fx*fY*fZ,11/3)+pow(fX*fy*fZ,11/3)+pow(fX*fY*fz,11/3)+pow(fx*fY*fz,11/3) + pow(fX*fy*fz,11/3) + pow(fx*fy*fZ,11/3) + pow(fx*fy*fz,11/3))    *v000.getR() * pow((fX*fY*fZ),11/3) + v100.getR() * pow(fx*fY*fZ,11/3) + v010.getR() * pow(fX*fy*fZ,11/3) + v001.getR() * pow(fX*fY*fz,11/3) +v101.getR() *pow(fx*fY*fz,11/3) + v011.getR() *pow(fX*fy*fz,11/3) + v110.getR() *pow(fx*fy*fZ,11/3) + v111.getR() *pow(fx*fy*fz,11/3);
			
			//Normieren auf neuen Betrag, falls nicht Nullvektor
			if (v.getR() != 0)
			{
			v = v * laenge/(v.getR());
			}
			//double theta =  v000.getTheta() * fX*fY*fZ + v100.getTheta() * fx*fY*fZ + v010.getTheta() * fX*fy*fZ + v001.getTheta() * fX*fY*fz +v101.getTheta() *fx*fY*fz + v011.getTheta() *fX*fy*fz + v110.getTheta() *fx*fy*fZ + v111.getTheta() *fx*fy*fz;
			
			//double phi =  v000.getPhi() * fX*fY*fZ + v100.getPhi() * fx*fY*fZ + v010.getPhi() * fX*fy*fZ + v001.getPhi() * fX*fY*fz +v101.getPhi() *fx*fY*fz + v011.getPhi() *fX*fy*fz + v110.getPhi() *fx*fy*fZ + v111.getPhi() *fx*fy*fz;
			
			//double z = v000.getZ() * fX*fY*fZ + v100.getZ() * fx*fY*fZ + v010.getZ() * fX*fy*fZ + v001.getZ() * fX*fY*fz +v101.getZ() *fx*fY*fz + v011.getZ() *fX*fy*fz + v110.getZ() *fx*fy*fZ + v111.getZ() *fx*fy*fz;
			//double theta = acos(z/laenge);
			
			//double x = v000.getX() * fX*fY*fZ + v100.getX() * fx*fY*fZ + v010.getX() * fX*fy*fZ + v001.getX() * fX*fY*fz +v101.getX() *fx*fY*fz + v011.getX() *fX*fy*fz + v110.getX() *fx*fy*fZ + v111.getX() *fx*fy*fz;
			//double y = v000.getR() * fX*fY*fZ + v100.getY() * fx*fY*fZ + v010.getY() * fX*fy*fZ + v001.getY() * fX*fY*fz +v101.getY() *fx*fY*fz + v011.getY() *fX*fy*fz + v110.getY() *fx*fy*fZ + v111.getY() *fx*fy*fz;
			//double phi = atan2(y,x);
			
			//~ std::cout << laenge;
			//~ std::cout << theta;
			//~ std::cout << phi;
			
			//v.setRThetaPhi(laenge,theta,phi);

		return v;

   }

double interpolate(double,const Vector3d &position) { //SKALAR
//std::cout << "T interpolate(T,const Vector3d &position) { //SKALAR" << std::endl;
//std::cout << "S"<< std::endl;
//~ std::cout << flush();
// position on a unit grid
		Vector3d r = (position - gridOrigin) / spacing;

		// indices of lower and upper neighbors
		int ix, iX, iy, iY, iz, iZ;
		if (reflective) {
			reflectiveClamp(r.x, Nx, ix, iX);
			reflectiveClamp(r.y, Ny, iy, iY);
			reflectiveClamp(r.z, Nz, iz, iZ);
		} else {
			periodicClamp(r.x, Nx, ix, iX);
			periodicClamp(r.y, Ny, iy, iY);
			periodicClamp(r.z, Nz, iz, iZ);
		}

		// linear fraction to lower and upper neighbors
		double fx = r.x - floor(r.x);
		double fX = 1 - fx;
		double fy = r.y - floor(r.y);
		double fY = 1 - fy;
		double fz = r.z - floor(r.z);
		double fZ = 1 - fz;

		// trilinear interpolation (see http://paulbourke.net/miscellaneous/interpolation)
		T b(0.);
		//V000 (1 - x) (1 - y) (1 - z) +
		b += get(ix, iy, iz) * fX * fY * fZ;
		//V100 x (1 - y) (1 - z) +
		b += get(iX, iy, iz) * fx * fY * fZ;
		//V010 (1 - x) y (1 - z) +
		b += get(ix, iY, iz) * fX * fy * fZ;
		//V001 (1 - x) (1 - y) z +
		b += get(ix, iy, iZ) * fX * fY * fz;
		//V101 x (1 - y) z +
		b += get(iX, iy, iZ) * fx * fY * fz;
		//V011 (1 - x) y z +
		b += get(ix, iY, iZ) * fX * fy * fz;
		//V110 x y (1 - z) +
		b += get(iX, iY, iz) * fx * fy * fZ;
		//V111 x y z
		b += get(iX, iY, iZ) * fx * fy * fz;

		return b;
	}

};

typedef Grid<Vector3f> VectorGrid;
typedef Grid<float> ScalarGrid;

} // namespace crpropa

#endif // CRPROPA_GRID_H
