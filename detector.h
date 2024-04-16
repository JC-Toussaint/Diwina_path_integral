#ifndef detector_h
#define detector_h

#include <vector>
#include "config.h"

/** \struct DETECTOR
 */

struct DETECTOR
    {
    /** \f$ zoomFactor \f$ for detector image  */
    double zoomFactor;

    /** \f$ meshSize \f$ for 2d regular triangular mesh   */
    double meshSize;
    };

#endif
