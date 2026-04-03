#ifndef sensor_h
#define sensor_h

#include <vector>
#include "config.h"

/** \struct SENSOR
 */

struct SENSOR
    {
    /** \f$ relative_size \f$ for image sensor */
    double relative_size;

    /** \f$ pixel_size \f$ for 2d regular triangular mesh   */
    double pixel_size;
    };

#endif
