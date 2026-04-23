/***************************************************
 * AhdcUtils
 * 
 * @author ftouchte
 * @date April 22, 2026
 ***************************************************/

#include "AhdcUtils.h"

std::vector<int> AhdcUtils::wire2slc(int wire) {
        //int sector = 1;
        int layer = -1;
        int component = -1;
        if (wire < 47) {
            layer = 11;
            component = wire + 1;
        }
        else if ((47 <= wire) && (wire < 47 + 56)) {
            layer = 21;
            component = wire - 47 + 1;
        }
        else if ((47 + 56 <= wire) && (wire < 47 + 56 + 56)) {
            layer = 22;
            component = wire - 47 - 56 + 1;
        }
        else if ((47 + 56 + 56 <= wire) && (wire < 47 + 56 + 56 + 72)) {
            layer = 31;
            component = wire - 47 - 56 - 56 + 1;
        }
        else if ((47 + 56 + 56 + 72 <= wire) && (wire < 47 + 56 + 56 + 72 + 72)) {
            layer = 32;
            component = wire - 47 - 56 - 56 - 72 + 1;
        }
        else if ((47 + 56 + 56 + 72 + 72 <= wire) && (wire < 47 + 56 + 56 + 72 + 72 + 87)) {
            layer = 41;
            component = wire - 47 - 56 - 56 - 72 - 72 + 1;
        }
        else if ((47 + 56 + 56 + 72 + 72 + 87 <= wire) && (wire < 47 + 56 + 56 + 72 + 72 + 87 + 87)) {
            layer = 42;
            component = wire - 47 - 56 - 56 - 72 - 72 - 87 + 1;
        }
        else { // ((47 + 56 + 56 + 72 + 72 + 87 + 87 <= wire) && (wire < 47 + 56 + 56 + 72 + 72 + 87 + 87 + 99)) {
            layer = 51;
            component = wire - 47 - 56 - 56 - 72 - 72 - 87 - 87 + 1;
        }
        return {1, layer, component};
    }

int AhdcUtils::slc2wire(int sector, int layer, int component) {
    if (layer == 11) {
        return component - 1;
    } 
    else if (layer == 21) {
        return 47 + component - 1;
    } 
    else if (layer == 22) {
        return 47 + 56 + component - 1;
    } 
    else if (layer == 31) {
        return 47 + 56 + 56 + component - 1;
    } 
    else if (layer == 32) {
        return 47 + 56 + 56 + 72 + component - 1;
    } 
    else if (layer == 41) {
        return 47 + 56 + 56 + 72 + 72 + component - 1;
    } 
    else if (layer == 42) {
        return 47 + 56 + 56 + 72 + 72 + 87 + component - 1;
    } 
    else if (layer == 51) {
        return 47 + 56 + 56 + 72 + 72 + 87 + 87 + component - 1;
    } else {
        return -1; // not a ahdc wire
    }
}

int AhdcUtils::layer2number(int digit) {
    if      (digit == 11) {
        return 1;
    } 
    else if (digit == 21) {
        return 2;
    } 
    else if (digit == 22) {
        return 3;
    } 
    else if (digit == 31) {
        return 4;
    } 
    else if (digit == 32) {
        return 5;
    } 
    else if (digit == 41) {
        return 6;
    } 
    else if (digit == 42) {
        return 7;
    } 
    else if (digit == 51) {
        return 8;
    } else {
        return 0; // not a layer, can encode all layers
    }
}

int AhdcUtils::number2layer(int num) {
    if      (num == 1) {
        return 11;
    } 
    else if (num == 2) {
        return 21;
    } 
    else if (num == 3) {
        return 22;
    } 
    else if (num == 4) {
        return 31;
    } 
    else if (num == 5) {
        return 32;
    } 
    else if (num == 6) {
        return 41;
    } 
    else if (num == 7) {
        return 42;
    } 
    else if (num == 8) {
        return 51;
    } else {
        return 0; // not a layer, can encode all layers
    }
}

double AhdcUtils::layer2Radius(int _layer) {
    if (_layer == 11) {
        return 32.0;
    }
    else if (_layer == 21) {
        return 38.0;
    }
    else if (_layer == 22) {
        return 42.0;
    }
    else if (_layer == 31) {
        return 48.0;
    }
    else if (_layer == 32) {
        return 52.0;
    }
    else if (_layer == 41) {
        return 58.0;
    }
    else if (_layer == 42) {
        return 62.0;
    }
    else if (_layer == 51) {
        return 68.0;
    } else {
        return 0.0;
    }
}

double  AhdcUtils::layerNum2Radius(int _layer_num) {
    return  layer2Radius(number2layer(_layer_num));
}

int AhdcUtils::layerNbWires(int _layer) {
    if (_layer == 11) {
        return 47;
    }
    else if (_layer == 21) {
        return 56;
    }
    else if (_layer == 22) {
        return 56;
    }
    else if (_layer == 31) {
        return 72;
    }
    else if (_layer == 32) {
        return 72;
    }
    else if (_layer == 41) {
        return 87;
    }
    else if (_layer == 42) {
        return 87;
    }
    else if (_layer == 51) {
        return 99;
    } else {
        return 0; // not a lyaer
    }
}
