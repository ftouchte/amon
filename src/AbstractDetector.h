/************************************
 * Abstract Detector
 * 
 * @date February 22, 2026
 * @author Felix Touchte Codjo
 ************************************/

#ifndef ABSTRACT_DETECTOR_H
#define ABSTRACT_DETECTOR_H

#include <cmath>
#include <vector>
#include <memory>
#include <map>

class Point3D {
    double dim[3];
public:
    Point3D(double x, double y, double z) {
        dim[0] = x;
        dim[1] = y;
        dim[2] = z;
    }

    double x() {return dim[0];}
    double y() {return dim[1];}
    double z() {return dim[2];}
    double r() {return sqrt(pow(x(),2) + pow(y(),2) + pow(z(),2));}
    double phi() {return atan2(x(), y());}
    double theta() {return acos(z()/r());}
};

class AbstractComponent {
    std::vector<Point3D> listOfPoints;
    int sector, superlayer, layer, component;
public:
    AbstractComponent(int _sector, int _superlayer, int _layer, int _component) {
        sector = _sector;
        superlayer = _superlayer;
        layer = _layer;
        component = _component;
    }
    AbstractComponent() = delete;
    int GetSectorId() const {return sector;}
    int GetSuperLayerId() const {return superlayer;}
    int GetLayerId() const {return layer;}
    int GetComponentId() const {return component;}
    void addPoint(const Point3D& pt) {listOfPoints.push_back(pt);}
    Point3D GetPoint(int i) const {return listOfPoints.at(i);}
    const std::vector<Point3D>& GetPoints() const {return listOfPoints;}
    int GetNumberOfPoints() const {return listOfPoints.size();}
    
    virtual ~AbstractComponent() = default;
};

class AbstractLayer {
    std::map<int, std::unique_ptr<AbstractComponent>> map_components;
    int sector, superlayer, layer;
public:
    AbstractLayer(int _sector, int _superlayer, int _layer) {
        sector = _sector;
        superlayer = _superlayer;
        layer = _layer;
    }
    AbstractLayer() = delete;
    int GetSectorId() const {return sector;}
    int GetSuperLayerId() const {return superlayer;}
    int GetLayerId() const {return layer;}
    void addComponent(int i, std::unique_ptr<AbstractComponent> comp) {
        map_components[i] = std::move(comp);
    }
    const AbstractComponent* GetComponent(int i) {
        auto it = map_components.find(i);
        return (it != map_components.end()) ? it->second.get() : nullptr;
    }
    int GetNumberOfComponents() const {return map_components.size();}

    virtual ~AbstractLayer() = default;
};

class AbstractSuperLayer {
    std::map<int, std::unique_ptr<AbstractLayer>> map_layers;
    int sector, superlayer;
public:
    AbstractSuperLayer(int _sector, int _superlayer) {
        sector = _sector;
        superlayer = _superlayer;
    }
    AbstractSuperLayer() = delete;
    int GetSectorId() const {return sector;}
    int GetSuperLayerId() const {return superlayer;}
    void addLayer(int i, std::unique_ptr<AbstractLayer> comp) {
        map_layers[i] = std::move(comp);
    }
    const AbstractLayer* GetLayer(int i) {
        auto it = map_layers.find(i);
        return (it != map_layers.end()) ? it->second.get() : nullptr;
    }
    int GetNumberOfLayers() const {return map_layers.size();}

    virtual ~AbstractSuperLayer() = default;
};


class AbstractSector {
    std::map<int, std::unique_ptr<AbstractSuperLayer>> map_superlayers;
    int sector;
public:
    AbstractSector(int _sector) { sector = _sector;}
    AbstractSector() = delete;
    int GetSectorId() const {return sector;}
    void addSuperLayer(int i, std::unique_ptr<AbstractSuperLayer> comp) {
        map_superlayers[i] = std::move(comp);
    }
    const AbstractSuperLayer* GetSuperLayer(int i) {
        auto it = map_superlayers.find(i);
        return (it != map_superlayers.end()) ? it->second.get() : nullptr;
    }
    int GetNumberOfSuperLayers() const {return map_superlayers.size();}

    virtual ~AbstractSector() = default;
};


class AbstractDetector {
    std::map<int, std::unique_ptr<AbstractSector>> map_sectors;
    std::string name;
public:
    AbstractDetector(std::string _name) : name(_name) {}
    std::string GetName() {return name;}
    void addSector(int i, std::unique_ptr<AbstractSector> comp) {
        map_sectors[i] = std::move(comp);
    }
    const AbstractSector* GetSector(int i) {
        auto it = map_sectors.find(i);
        return (it != map_sectors.end()) ? it->second.get() : nullptr;
    }
    int GetNumberOfSectors() const {return map_sectors.size();}

    virtual ~AbstractDetector() = default;
};

class DetectorFactory {
public:
    virtual ~DetectorFactory() = default;
    virtual const AbstractDetector* create() = 0;
};



#endif
