#ifndef node_h
#define node_h

/** \file node.h
\brief header to define struct Node
*/

//#include <memory>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include "config.h"

/**
 \namespace Nodes
 to grab altogether some dedicated functions and data of the nodes
 */
namespace Nodes
    {
     /** space dimension */
    const int DIM = 3;
    
    /** size of the array of dataNode */
    const int NB_DATANODE = 2;

    /** convenient enum to avoid direct index manipulation in the array of struct dataNode */
    enum step
        {
        CURRENT = 0,
        NEXT = 1
        };

    /** convenient enum to specify what coordinate in calculations */
    enum index
        {
        IDX_UNDEF = -1,
        IDX_X = 0,
        IDX_Y = 1,
        IDX_Z = 2
        };

/** \return \f$ x^2 \f$ */
inline double sq(const double x) { return x * x; }

/** \struct dataNode
contains the vector field u
*/

struct dataNode
    {
    Eigen::Vector3d u;  /**< magnetization */
    };

/** \struct Node
Node is containing physical point of coordinates \f$ p = (x,y,z) \f$, magnetization value at \f$
m(p,t) \f$.
*/
struct Node
    {
    Eigen::Vector3d p;  /**< Physical position p=(x,y,z)  of the node */
    
    /** datas associated to position p
    step CURRENT (0) : start of the time step
    step NEXT (1) : after the current time step 
    */ 
    dataNode d[NB_DATANODE];

    /** getter for u at step k */
    inline const Eigen::Vector3d get_u(step k /**< [in] */) const
    { return d[k].u; }

    };  // end struct node

/** getter for magnetizations */
template <step K>
Eigen::Vector3d get_u(Node const &n /**< [in] */) { return n.d[K].u; }

    }  // end namespace Nodes

#endif /* node_h */
