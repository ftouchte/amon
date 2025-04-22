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
AhdcWire::AhdcWire(Point3D _top, Point3D _bot) : top(_top), bot(_bot) {
	set_z(0.0);
}

/** (Default) Constructor */
AhdcWire::AhdcWire(double x1, double y1, double z1, double x2, double y2, double z2) { 
	top = Point3D(x1,y1,z1); 
	bot = Point3D(x2,y2,z2);
	set_z(0.0);
}

AhdcWire & AhdcWire::operator=(const AhdcWire & obj) {
	if (this != &obj) {
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
	x = top.x + z*ux;
	y = top.y + z*uy;
	//printf("%lf %lf %lf\n", x, y, z);
}

/*******************************
 * AhdcLayer
 * ****************************/

/** (Default) Constructor */
AhdcLayer::AhdcLayer(int _nwires, double _rlayer, double _stereoangle) : nwires(_nwires), rlayer(_rlayer), stereoangle(_stereoangle) {
	if (nwires > 0) {
		ptrWires = new AhdcWire[nwires];
		for (int id = 0; id < nwires; id++){
			// All distance is in mm !
			Point3D top, bot;
			double alpha = futils::toRadian(360.0/nwires);
			top.x = -rlayer*sin(alpha*id); 
			top.y = -rlayer*cos(alpha*id);
			top.z = 0; // sometimes set to -150
			bot.x = -rlayer*sin(alpha*id + futils::toRadian(stereoangle));
			bot.y = -rlayer*cos(alpha*id + futils::toRadian(stereoangle));
			bot.z = 300; // sometimes set to +150
			ptrWires[id] = AhdcWire(top,bot);
		}
	}
}

/** Copy constructor */
AhdcLayer::AhdcLayer(const AhdcLayer & obj){
	nwires = obj.nwires;
	rlayer = obj.rlayer;
	stereoangle = obj.stereoangle;
	ptrWires = new AhdcWire[nwires];
	for (int id = 0; id < nwires; id++){
		ptrWires[id] = obj.ptrWires[id];
	}
}


/** Affectation operator */
AhdcLayer & AhdcLayer::operator=(const AhdcLayer & obj){
	if (this != &obj) {
		nwires = obj.nwires;
		rlayer = obj.rlayer;
		stereoangle = obj.stereoangle;
		if (ptrWires != nullptr)
			delete[] ptrWires;
		ptrWires = new AhdcWire[nwires];
		for (int id = 0; id < nwires; id++){
			ptrWires[id] = obj.ptrWires[id];
		}
	}
	return *this;
}

/** Destructor */
AhdcLayer::~AhdcLayer(){
	if (ptrWires != nullptr)
		delete[] ptrWires;
}

/** Return a pointer to the wire number id 
 * @class assuming the numerotation starts at 0 !
 */
AhdcWire* AhdcLayer::GetWire(int id) {return ptrWires + id;} // arithmétiques des pointeurs !

/** Return the number of wires */
int AhdcLayer::GetNumberOfWires() {return nwires;}


/*******************************
 * AhdcSuperLayer
 * ****************************/

/** (Default) Constructor */
AhdcSuperLayer::AhdcSuperLayer(int _nlayers, int _nwires, double _rsuperlayer, int _orientation) : nlayers(_nlayers), nwires(_nwires), rsuperlayer(_rsuperlayer), orientation(_orientation) {
	if (nlayers > 0) {
		ptrLayers = new AhdcLayer[nlayers];
		for (int id = 0; id < nlayers; id++){
			std::cout << "BUILD LAYER " << id + 1 << std::endl;
			ptrLayers[id] = AhdcLayer(nwires,rsuperlayer + 4.0*id,orientation*20); // 4.0 is the distance between two layer of a superlayer in AHDC
		}
	}
}

/** Copy constructor */
AhdcSuperLayer::AhdcSuperLayer(const AhdcSuperLayer & obj){
	nlayers = obj.nlayers;
	nwires = obj.nwires;
	rsuperlayer = obj.rsuperlayer;
	orientation = obj.orientation;
	ptrLayers = new AhdcLayer[nlayers];
	for (int id = 0; id < nlayers; id++){
		ptrLayers[id] = obj.ptrLayers[id];
	}
}

/** Affectation operator */
AhdcSuperLayer & AhdcSuperLayer::operator=(const AhdcSuperLayer & obj){
	if (this != &obj) {
		nlayers = obj.nlayers;
		nwires = obj.nwires;
		rsuperlayer = obj.rsuperlayer;
		orientation = obj.orientation;
		if (ptrLayers != nullptr) 
			delete[] ptrLayers;
		ptrLayers = new AhdcLayer[nlayers];
		for (int id = 0; id < nlayers; id++){
			ptrLayers[id] = obj.ptrLayers[id];
		}
	}
	return *this;
}
/** Destructor */
AhdcSuperLayer::~AhdcSuperLayer(){
	if (ptrLayers != nullptr)
		delete[] ptrLayers;
}

/** Return a pointer to the layer number id 
 * @class assuming the numerotation starts at 0 !
 */
AhdcLayer* AhdcSuperLayer::GetLayer(int id) {return ptrLayers + id;}

/** Return the number of layer */
int AhdcSuperLayer::GetNumberOfLayers(){ return nlayers;}



/*******************************
 * AhdcSector
 * ****************************/

AhdcSector::AhdcSector(int _nsuperlayers) : nsuperlayers(_nsuperlayers) {
	if (nsuperlayers > 0) {
		ptrSuperLayers = new AhdcSuperLayer[nsuperlayers];
		for (int id = 0; id < nsuperlayers; id++) {
			std::cout << "BUILD SUPERLAYER " << id + 1 << std::endl;
			int nlayers; 
			int nwires;
			double rsuperlayer;
			int orientation = pow(-1.0,id);
			if (id == 0) {
				nlayers     = 1;
				nwires      = 47;
				rsuperlayer = 32.0;
			} else if (id == 1) {
				nlayers     = 2;
				nwires      = 56;
				rsuperlayer = 38.0;
			} else if (id == 2) {
				nlayers     = 2;
				nwires      = 72;
				rsuperlayer = 48;
			} else if (id == 3) {
				nlayers     = 2;
				nwires      = 87;
				rsuperlayer = 58.0;
			} else { // (id == 4)
				nlayers     = 1;
				nwires      = 99;
				rsuperlayer = 68;
			}
			ptrSuperLayers[id] = AhdcSuperLayer(nlayers,nwires,rsuperlayer,orientation);
		}
	}
}

/** Copy constructor */
AhdcSector::AhdcSector(const AhdcSector & obj){
	nsuperlayers = obj.nsuperlayers;
	ptrSuperLayers = new AhdcSuperLayer[nsuperlayers];
	for (int id = 0; id < nsuperlayers; id++){
		ptrSuperLayers[id] = obj.ptrSuperLayers[id];
	}
}

/** Affectation operator */
AhdcSector & AhdcSector::operator=(const AhdcSector & obj){
	if (this != &obj) {
		nsuperlayers = obj.nsuperlayers;
		if (ptrSuperLayers != nullptr)
			delete[] ptrSuperLayers;
		ptrSuperLayers = new AhdcSuperLayer[nsuperlayers];
		for (int id = 0; id < nsuperlayers; id++){
			ptrSuperLayers[id] = obj.ptrSuperLayers[id];
		}
	}
	return *this;
}

/** Destructor */
AhdcSector::~AhdcSector(){
	if (ptrSuperLayers != nullptr)
		delete[] ptrSuperLayers;
}

/** Return a pointer to the super layer number id 
 * @class assuming the numerotation starts at 0 !
 */
AhdcSuperLayer* AhdcSector::GetSuperLayer(int id){return ptrSuperLayers + id;}

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
		for (int id = 0; id < nsectors; id++) {
			std::cout << "BUILD SECTOR " << id + 1 << std::endl;
			ptrSectors[id] = AhdcSector(5);
		}
	}
}

/** Copy constructor */
AhdcDetector::AhdcDetector(const AhdcDetector & obj){
	nsectors = obj.nsectors;
	ptrSectors = new AhdcSector[nsectors];
	for (int id = 0; id < nsectors; id++) {
		ptrSectors[id] = obj.ptrSectors[id];
	}
}

/** Affectation operator */
AhdcDetector & AhdcDetector::operator=(const AhdcDetector & obj){
	if (this != &obj) {
		nsectors = obj.nsectors;
		if (ptrSectors != nullptr)
			delete[] ptrSectors;
		ptrSectors = new AhdcSector[nsectors];
		for (int id = 0; id < nsectors; id++) {
			ptrSectors[id] = obj.ptrSectors[id];
		}
	}
	return *this;
}

/** Destructor */
AhdcDetector::~AhdcDetector() {
	if (ptrSectors != nullptr)
		delete[] ptrSectors;
}

/** Return a pointer to the super layer number id 
 * @class assuming the numerotation starts at 0 !
 */
AhdcSector * AhdcDetector::GetSector(int id) {return ptrSectors + id;}

/** Return the number of super layer */
int AhdcDetector::GetNumberOfSectors() {return nsectors;}



