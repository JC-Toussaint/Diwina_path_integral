/**
  Elementary matrix Calculation for a tetrahedron element
 */

#include <set>

#include "tetra.h"
#include "facette.h"

using namespace Tetra;
using namespace Nodes;

double Tet::Jacobian(Eigen::Ref<Eigen::Matrix3d> J)
    {
    Eigen::Vector3d p0p1 = getNode(1).p - getNode(0).p;
    Eigen::Vector3d p0p2 = getNode(2).p - getNode(0).p;
    Eigen::Vector3d p0p3 = getNode(3).p - getNode(0).p;
    J(0,0) = p0p1.x();
    J(0,1) = p0p2.x();
    J(0,2) = p0p3.x();
    J(1,0) = p0p1.y();
    J(1,1) = p0p2.y();
    J(1,2) = p0p3.y();
    J(2,0) = p0p1.z();
    J(2,1) = p0p2.z();
    J(2,2) = p0p3.z();
    return J.determinant();
    }

double Tet::calc_vol(void) const
    {
    Eigen::Vector3d p0p1 = getNode(1).p - getNode(0).p;
    Eigen::Vector3d p0p2 = getNode(2).p - getNode(0).p;
    Eigen::Vector3d p0p3 = getNode(3).p - getNode(0).p;

    return p0p1.dot(p0p2.cross(p0p3))/6.0;
    }

