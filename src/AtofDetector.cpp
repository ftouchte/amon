/************************************
 * ATOF Detector
 * 
 * @date February 27, 2026
 * @author Felix Touchte Codjo
 ************************************/

#include "AtofDetector.h"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>


const AtofDetector * AtofFactory::create() {
    std::ifstream file(path2geometry.c_str());
    if (!file) {
        throw std::runtime_error("Cannot openn the file : " + path2geometry);
    }

    // Intialise the ATOF detector
    auto ATOFdet = std::unique_ptr<AtofDetector>();

    // list of components
    std::vector<std::unique_ptr<AtofComponent>> all_components;

    std::string line;
    int row = 0;
    int sector, superlayer, layer, component;
    double z_top,z_bottom,x0,y0,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7;

    // loop over lines and fill the list of components
    while (std::getline(file, line)) {
        // ignore the first line
        if (row == 0) continue;
        row++;
        // use this line as a flux
        std::istringstream iss(line);
        // load values
        iss >> z_top >> z_bottom >> x0 >> y0 >> x1 >> y1 >> x2 >> y2 >> x3 >> y3 >> x4 >> y4 >> x5 >> y5 >> x6 >> y6 >> x7 >> y7;
        Point3D p0(x0, y0, z_bottom);
        Point3D p1(x1, y1, z_bottom);
        Point3D p2(x2, y2, z_bottom);
        Point3D p3(x3, y3, z_bottom);
        Point3D p4(x4, y4, z_top);
        Point3D p5(x5, y5, z_top);
        Point3D p6(x6, y6, z_top);
        Point3D p7(x7, y7, z_top);
        auto atof_comp = std::make_unique<AtofComponent>(sector, superlayer, layer, component);
        atof_comp->addPoint(p0);
        atof_comp->addPoint(p1);
        atof_comp->addPoint(p2);
        atof_comp->addPoint(p3);
        atof_comp->addPoint(p4);
        atof_comp->addPoint(p5);
        atof_comp->addPoint(p6);
        atof_comp->addPoint(p7);
        all_components.push_back(std::move(atof_comp));
    }

    //////////////////////////////////////////
    /// Astuce: faire des maps de map<sectorId, uniqueOrderedId> et stocker les infos dans des std::vector
    /// ce sera beaucoup plus lisible
    //////////////////////////////////////////
    


    ////////////////////
    // sort by sector
    ////////////////////
    std::map<int, std::vector<std::unique_ptr<AtofComponent>>> map_components_by_sectors;
    for (size_t i = 0; i < all_components.size(); i++) {
        int sectorId = all_components[i]->GetSectorId();
        auto it = map_components_by_sectors.find(sectorId);
        if (it != map_components_by_sectors.end()) {
            it->second.push_back(std::move(all_components[i]));
        } else {
            // if it is a new sector, initialise the vector
            map_components_by_sectors[i] = std::vector<std::unique_ptr<AtofComponent>>({std::move(all_components[i])});
        }
    }
    // here all th entries should be nullptr
    all_components.clear();

    ////////////////////
    // sort by sector
    ////////////////////

    // créer une fonction: sort_by_id (sector, superlayer, layer) à partir d'un vecteur de component
    for (auto& [sectorId, vec] : map_components_by_sectors) {
        auto atof_sector = std::make_unique<AtofSector>(sectorId);
        //std::cout << "Sector " << sectorId << " has " << vec.size() << " components\n";

        // sort by superlayer
        std::map<int, std::vector<std::unique_ptr<AtofComponent>>> map_superlayers;
        for (auto& comp : vec) {

            //std::cout << "  Component ID: " << comp->GetId() << "\n";
        }
    }

    

}
