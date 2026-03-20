/************************************
 * ATOF Detector
 * 
 * @date February 26, 2026
 * @author Felix Touchte Codjo
 ************************************/

#ifndef ATOF_DETECTOR_H
#define ATOF_DETECTOR_H

#include "AbstractDetector.h"
#include <string>
#include <filesystem>

class AtofComponent : public AbstractComponent {
public:
    AtofComponent(int _sector, int _superlayer, int _layer, int _component) : AbstractComponent(_sector,_superlayer,_layer,_component) {}
    AtofComponent() = delete;
};

class AtofLayer : public AbstractLayer {
public: 
    AtofLayer(int _sector, int _superlayer, int _layer) : AbstractLayer(_sector,_superlayer,_layer) {}
    AtofLayer() = delete;
};

class AtofSuperLayer : public AbstractSuperLayer  {
    AtofSuperLayer(int _sector, int _superlayer) : AbstractSuperLayer(_sector,_superlayer) {}
    AtofSuperLayer() = delete;
};

class AtofSector : public AbstractSector {
    AtofSector(int _sector) : AbstractSector(_sector) {}
    AtofSector() = delete;
};

class AtofDetector : public AbstractDetector {
public:
    AtofDetector() : AbstractDetector("ATOF") {}
    AtofDetector() = delete;
};

class AtofFactory : public DetectorFactory {
    std::string path2geometry = "../ressources/ATOF_geometry.txt";
public:
    std::string GetGeometryDefinition() {return path2geometry;}
    const AtofDetector* create() override;
};


#endif