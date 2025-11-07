/************************************************
 * AhdcDetector.cpp
 *
 * @author Felix Touchte Codjo
 * @date November 28, 2024
 * **********************************************/

#include "AhdcDetector.h"
#include <cmath>
#include <iostream>

double futils::toRadian(double degree) { return M_PI*degree/180;}

/*******************************
 * Ahdcwire
 * ****************************/

/** Constructor */
AhdcWire::AhdcWire(int _id, Point3D _top, Point3D _bot) : id(_id), top(_top), bot(_bot) {
	set_z(-150.0);
}

/** (Default) Constructor */
/*AhdcWire::AhdcWire(int _id, double x1, double y1, double z1, double x2, double y2, double z2) : id(_id) { 
	top = Point3D(x1,y1,z1); 
	bot = Point3D(x2,y2,z2);
	set_z(0.0);
}*/

AhdcWire & AhdcWire::operator=(const AhdcWire & obj) {
	if (this != &obj) {
		id = obj.id;
		top = obj.top;
		bot = obj.bot;
		x = obj.x;
		y = obj.y;
		z = obj.z;
		occ = obj.occ;
	}
	return *this;
}

void AhdcWire::set_z(double _z) {
	double ux = bot.x - top.x;
	double uy = bot.y - top.y;
	double uz = bot.z - top.z;
	double rho = sqrt(ux*ux + uy*uy + uz*uz);
	ux = ux/rho;
	uy = uy/rho;
	uz = uz/rho;
	z = _z;
	x = top.x + (z+150)*ux;
	y = top.y + (z+150)*uy;
	//printf("%lf %lf %lf\n", x, y, z);
}

/*******************************
 * AhdcLayer
 * ****************************/

/** (Default) Constructor */
AhdcLayer::AhdcLayer(int _id, int _nwires, double _rlayer, double _stereoangle) : id(_id), nwires(_nwires), rlayer(_rlayer), stereoangle(_stereoangle) {
	if (nwires > 0) {
		ptrWires = new AhdcWire[nwires];
		for (int i = 0; i < nwires; i++){
			// All distance is in mm !
			// From coatjava
			// start at phi=0
                        // in each layer the first wire is the first at phi>=0, i.e.
                        // 0.5 0 0.5 0 0.5 0.5 0 0.5 for layer 1 to 8
			// geometry from April 22, 2025
                        double wirePhiIndex = i + 0.5*(nwires%2) + 0.5*(id-1)*(1-2*(nwires%2)); // here layer id starts at 1
			//printf("layer id : %d\n", id);
			Point3D top, bot;
			double alpha = futils::toRadian(360.0/nwires);
			top.x = rlayer*cos(alpha*wirePhiIndex); 
			top.y = rlayer*sin(alpha*wirePhiIndex);
			top.z = -150; // sometimes set to -150
			bot.x = rlayer*cos(alpha*wirePhiIndex + futils::toRadian(stereoangle));
			bot.y = rlayer*sin(alpha*wirePhiIndex + futils::toRadian(stereoangle));
			bot.z = 150; // sometimes set to +150
			ptrWires[i] = AhdcWire(i+1,top,bot);
		}
	}
}

/** Copy constructor */
AhdcLayer::AhdcLayer(const AhdcLayer & obj){
	id = obj.id;
	nwires = obj.nwires;
	rlayer = obj.rlayer;
	stereoangle = obj.stereoangle;
	ptrWires = new AhdcWire[nwires];
	for (int i = 0; i < nwires; i++){
		ptrWires[i] = obj.ptrWires[i];
	}
}


/** Affectation operator */
AhdcLayer & AhdcLayer::operator=(const AhdcLayer & obj){
	if (this != &obj) {
		id = obj.id;
		nwires = obj.nwires;
		rlayer = obj.rlayer;
		stereoangle = obj.stereoangle;
		if (ptrWires != nullptr)
			delete[] ptrWires;
		ptrWires = new AhdcWire[nwires];
		for (int i = 0; i < nwires; i++){
			ptrWires[i] = obj.ptrWires[i];
		}
	}
	return *this;
}

/** Destructor */
AhdcLayer::~AhdcLayer(){
	if (ptrWires != nullptr)
		delete[] ptrWires;
}

/** Return a pointer to the wire number i 
 * @class assuming the numerotation starts at 0 !
 */
AhdcWire* AhdcLayer::GetWire(int i) {return ptrWires + i;} // arithmétiques des pointeurs !

/** Return the number of wires */
int AhdcLayer::GetNumberOfWires() {return nwires;}


/*******************************
 * AhdcSuperLayer
 * ****************************/

/** (Default) Constructor */
AhdcSuperLayer::AhdcSuperLayer(int _id, int _nlayers, int _nwires, double _rsuperlayer, int _orientation) : id(_id), nlayers(_nlayers), nwires(_nwires), rsuperlayer(_rsuperlayer), orientation(_orientation) {
	if (nlayers > 0) {
		ptrLayers = new AhdcLayer[nlayers];
		for (int i = 0; i < nlayers; i++){
			std::cout << "BUILD LAYER " << i + 1 << std::endl;
			ptrLayers[i] = AhdcLayer(i+1,nwires,rsuperlayer + 4.0*i,orientation*(-20)); // 4.0 is the distance between two layer of a given superlayer in AHDC
		}
	}
}

/** Copy constructor */
AhdcSuperLayer::AhdcSuperLayer(const AhdcSuperLayer & obj){
	id = obj.id;
	nlayers = obj.nlayers;
	nwires = obj.nwires;
	rsuperlayer = obj.rsuperlayer;
	orientation = obj.orientation;
	ptrLayers = new AhdcLayer[nlayers];
	for (int i = 0; i < nlayers; i++){
		ptrLayers[i] = obj.ptrLayers[i];
	}
}

/** Affectation operator */
AhdcSuperLayer & AhdcSuperLayer::operator=(const AhdcSuperLayer & obj){
	if (this != &obj) {
		id = obj.id;
		nlayers = obj.nlayers;
		nwires = obj.nwires;
		rsuperlayer = obj.rsuperlayer;
		orientation = obj.orientation;
		if (ptrLayers != nullptr) 
			delete[] ptrLayers;
		ptrLayers = new AhdcLayer[nlayers];
		for (int i = 0; i < nlayers; i++){
			ptrLayers[i] = obj.ptrLayers[i];
		}
	}
	return *this;
}
/** Destructor */
AhdcSuperLayer::~AhdcSuperLayer(){
	if (ptrLayers != nullptr)
		delete[] ptrLayers;
}

/** Return a pointer to the layer number i 
 * @class assuming the numerotation starts at 0 !
 */
AhdcLayer* AhdcSuperLayer::GetLayer(int i) {return ptrLayers + i;}

/** Return the number of layer */
int AhdcSuperLayer::GetNumberOfLayers(){ return nlayers;}



/*******************************
 * AhdcSector
 * ****************************/

AhdcSector::AhdcSector(int _nsuperlayers) : nsuperlayers(_nsuperlayers) {
	if (nsuperlayers > 0) {
		ptrSuperLayers = new AhdcSuperLayer[nsuperlayers];
		for (int i = 0; i < nsuperlayers; i++) {
			std::cout << "BUILD SUPERLAYER " << i + 1 << std::endl;
			int nlayers; 
			int nwires;
			double rsuperlayer;
			int orientation = pow(-1.0,i);
			if (i == 0) {
				nlayers     = 1;
				nwires      = 47;
				rsuperlayer = 32.0;
			} else if (i == 1) {
				nlayers     = 2;
				nwires      = 56;
				rsuperlayer = 38.0;
			} else if (i == 2) {
				nlayers     = 2;
				nwires      = 72;
				rsuperlayer = 48;
			} else if (i == 3) {
				nlayers     = 2;
				nwires      = 87;
				rsuperlayer = 58.0;
			} else { // (i == 4)
				nlayers     = 1;
				nwires      = 99;
				rsuperlayer = 68;
			}
			ptrSuperLayers[i] = AhdcSuperLayer(i+1,nlayers,nwires,rsuperlayer,orientation);
		}
	}
}

/** Copy constructor */
AhdcSector::AhdcSector(const AhdcSector & obj){
	nsuperlayers = obj.nsuperlayers;
	ptrSuperLayers = new AhdcSuperLayer[nsuperlayers];
	for (int i = 0; i < nsuperlayers; i++){
		ptrSuperLayers[i] = obj.ptrSuperLayers[i];
	}
}

/** Affectation operator */
AhdcSector & AhdcSector::operator=(const AhdcSector & obj){
	if (this != &obj) {
		nsuperlayers = obj.nsuperlayers;
		if (ptrSuperLayers != nullptr)
			delete[] ptrSuperLayers;
		ptrSuperLayers = new AhdcSuperLayer[nsuperlayers];
		for (int i = 0; i < nsuperlayers; i++){
			ptrSuperLayers[i] = obj.ptrSuperLayers[i];
		}
	}
	return *this;
}

/** Destructor */
AhdcSector::~AhdcSector(){
	if (ptrSuperLayers != nullptr)
		delete[] ptrSuperLayers;
}

/** Return a pointer to the super layer number i 
 * @class assuming the numerotation starts at 0 !
 */
AhdcSuperLayer* AhdcSector::GetSuperLayer(int i){return ptrSuperLayers + i;}

/** Return the number of super layer */
int AhdcSector::GetNumberOfSuperLayers() {return nsuperlayers;}



/*******************************
 * AhdcDetector
 * ****************************/

/** (Default) constructor */
AhdcDetector::AhdcDetector() : nsectors(1) {
	if (nsectors > 0) {
		std::cout << "BUILD DETECTOR : AHDC" << std::endl;
		ptrSectors = new AhdcSector[nsectors];
		for (int i = 0; i < nsectors; i++) {
			std::cout << "BUILD SECTOR " << i + 1 << std::endl;
			ptrSectors[i] = AhdcSector(5);
		}
	}
}

/** Copy constructor */
AhdcDetector::AhdcDetector(const AhdcDetector & obj){
	nsectors = obj.nsectors;
	ptrSectors = new AhdcSector[nsectors];
	for (int i = 0; i < nsectors; i++) {
		ptrSectors[i] = obj.ptrSectors[i];
	}
}

/** Affectation operator */
AhdcDetector & AhdcDetector::operator=(const AhdcDetector & obj){
	if (this != &obj) {
		nsectors = obj.nsectors;
		if (ptrSectors != nullptr)
			delete[] ptrSectors;
		ptrSectors = new AhdcSector[nsectors];
		for (int i = 0; i < nsectors; i++) {
			ptrSectors[i] = obj.ptrSectors[i];
		}
	}
	return *this;
}

/** Destructor */
AhdcDetector::~AhdcDetector() {
	if (ptrSectors != nullptr)
		delete[] ptrSectors;
}

/** Return a pointer to the super layer number i 
 * @class assuming the numerotation starts at 0 !
 */
AhdcSector * AhdcDetector::GetSector(int i) {return ptrSectors + i;}

/** Return the number of super layer */
int AhdcDetector::GetNumberOfSectors() {return nsectors;}



