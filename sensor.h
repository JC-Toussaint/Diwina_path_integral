#ifndef sensor_h
#define sensor_h

#include <vector>
#include "config.h"

/** \struct SENSOR
 */

struct SENSOR
    {
    /** \f$ zoomFactor \f$ for image sensor */
    double zoomFactor;

    /** \f$ meshSize \f$ for 2d regular triangular mesh   */
    double meshSize;
    };

#endif
