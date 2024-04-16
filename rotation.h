#ifndef rotation_h
#define rotation_h

#include <vector>
#include "config.h"

/** \struct ROT
 */

struct ROT
    {
    /** \f$ angles \f$ of rotation  */
    double angle1, angle2;

    /** \f$ axes \f$ of rotation   */
    std::vector<double> axe1, axe2;

    /** printing function */
    inline void infos(void)
        {
        std::cout << "Axe1   : " << axe1[0] << ';' << axe1[1] << ';' << axe1[2] << std::endl;
	    std::cout << "Angle1 : " << angle1 << "deg\n";
        std::cout << "Axe2   : " << axe2[0] << ';' << axe2[1] << ';' << axe2[2] << std::endl;
	    std::cout << "Angle2 : " << angle2 << "deg\n";
        }
    };

#endif
