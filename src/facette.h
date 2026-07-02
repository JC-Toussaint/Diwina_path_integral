#ifndef facette_h
#define facette_h

/** \file facette.h
  \brief contains namespace Facette
  header containing Fac class, and some constants and a less_than operator to redo orientation of
  triangular faces
 */

#include "element.h"

/** \namespace Facette
 to grab altogether some constants and calculation functions for class Fac
 */
namespace Facette
    {
const int N = 3;   /**< number of sommits */

/** \class prm
region number and material constants
*/
struct prm
    {
    std::string regName;   /**< region name */

    /** print the struct parameters */
    inline void infos()
        {
        std::cout << "surface region name = " << regName << std::endl;
        };
    };

/** \class Fac
Face is a class containing the index references to nodes, its surface and its normal unit vector, it has a triangular shape and should not
be degenerated, orientation must be defined in adequation with the mesh
*/
class Fac : public element<N>
    {
public:
    /** constructor used by readMesh */
    inline Fac(const std::vector<Nodes::Node> &_p_node /**< [in] vector of nodes */,
               const int _NOD /**< [in] nb nodes */,
               const int _idx /**< [in] region index in region vector */,
               std::initializer_list<int> _i /**< [in] node index */)
        : element<N>(_p_node,_idx,_i)
        {
        if (_NOD > 0)
            {
            zeroBasing();
            surf = calc_surf();
            n = calc_norm();
            }
        else
            {
            surf = 0.0;
            n = Eigen::Vector3d(0,0,0);
            }  // no index shift here if NOD == 0 : usefull while reordering face indices
        }

    /** surface of the element */
    double surf;

    /** normal vector (unit vector) */
    Eigen::Vector3d n;

    /** lexicographic order on indices */
    inline bool operator<(const Fac &f) const
        {
        return (this->ind[0] < f.ind[0])
               || ((this->ind[0] == f.ind[0])
                   && ((this->ind[1] < f.ind[1])
                       || ((this->ind[1] == f.ind[1]) && (this->ind[2] < f.ind[2]))));
        }

    /** computes the norm to the face, returns a unit vector */
    inline Eigen::Vector3d calc_norm(void) const
        {
        Eigen::Vector3d _n = normal_vect();
        _n.normalize();
        return _n;
        }

private:
    void orientate(void) {}// orientation is done in mesh::indexReorder

    /** computes surface of the face */
    inline double calc_surf(void) const { return 0.5 * normal_vect().norm(); }

    /** return normal to the triangular face, not normalized */
    inline Eigen::Vector3d normal_vect() const
        {
        Eigen::Vector3d p0p1 = getNode(1).p - getNode(0).p;
        Eigen::Vector3d p0p2 = getNode(2).p - getNode(0).p;

        return p0p1.cross(p0p2);
        }
    };  // end class Fac

    }  // namespace Facette

#endif /* facette_h */
