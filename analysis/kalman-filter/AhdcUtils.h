/***************************************************
 * AhdcUtils
 * 
 * @author ftouchte
 * @date April 22, 2026
 ***************************************************/

#ifndef AHDC_UTILS_H
#define AHDC_UTILS_H

#include <vector>


namespace AhdcUtils {

    /**
     * Convert wire number (number from 0 to 575) to (sector,layer,component) ids
     * 
     * This is the invert operation of  @link slc2wire(int, int, int) @endlink 
     * 
     * @param wire wire number between 0 and  576 (excluded)
     * @return a triplet (sector, layer, component) in std::vector<int>
     */
    std::vector<int> wire2slc(int wire);

    /**
     * Convert (sector, layer, component) to a unique wire id (number betwwen 0 and 575)
     * 
     * @param sector (not used)
     * @param layer 
     * @param component 
     * @return unique wire id, return -1 if the slc is not recognized
     */
    int slc2wire(int sector, int layer, int component);

    /**
     * Convert the digit-layer (11,21,...,51) to layer number between 1 and 8
     * 
     * @param digit 
     * @return layer number
     */
    int layer2number(int digit);

    /**
     * Convert layer number (from 1 to 8) to the superlayer-layer id (11,21,...,51)
     * 
     * @param digit between 1 and 8
     * @return layer number between (11,21,...,51), return 0 otherwise
     */
    int number2layer(int num);

    /**
     * 
     * @param _layer (number 11, 21, 22, ..., 51)
     * @return the radius of the _layer
     */
    double layer2Radius(int _layer);

    /**
     * 
     * @param _layer_num between 1 and 8
     * @return the radius of the _layer. See {@link #layer2Radius(int)}
     */
    double layerNum2Radius(int _layer_num);

}

#endif