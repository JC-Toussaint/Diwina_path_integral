#ifndef tetra_h
#define tetra_h

/** \file tetra.h
  \brief namespace Tetra
  header containing Tet class, some constants, and integrales
 */

#include <set>

#include "facette.h"
#include "node.h"
#include "element.h"

/** \namespace Tetra
 to grab altogether some constants for struct Tet
 */
namespace Tetra
    {
const int N = 4;   /**< number of sommits */

const double epsilon =
        EPSILON; /**< this constant is defined from a macro in config.h.in, it is used to check the
                    validity of the tetrahedreon, a degeneracy test */

/** \class prm
region number and material constants
*/
struct prm
    {
    std::string regName; /**< region name */
    double J;            /**< \f$ M_s = \nu_0 J \f$ */

    double lbp, lbm; /**< XMCD absorption parameters */

    /** print the struct parameters */
    inline void infos()
        {
        std::cout << "volume region name = " << regName << std::endl;
        std::cout << "J = " << J << std::endl;
        };
    };

/** \class Tet
Tet is a tetrahedron, containing the index references to nodes, must not be flat <br>
indices convention is<br>
```
                        v
                      .
                    ,/
                   /
                2(ic)                                 2
              ,/|`\                                 ,/|`\
            ,/  |  `\                             ,/  |  `\
          ,/    '.   `\                         ,6    '.   `5
        ,/       |     `\                     ,/       8     `\
      ,/         |       `\                 ,/         |       `\
     0(ia)-------'.--------1(ib) --> u     0--------4--'.--------1
      `\.         |      ,/                 `\.         |      ,/
         `\.      |    ,/                      `\.      |    ,9
            `\.   '. ,/                           `7.   '. ,/
               `\. |/                                `\. |/
                  `3(id)                                `3
                     `\.
                        ` w
```
*/
class Tet : public element<N>
    {
public:
    /** constructor. It initializes dad(x|y|z) if \f$ | detJ | < \epsilon \f$ jacobian
    is considered degenerated unit tests : Tet_constructor; Tet_inner_tables
    */
    inline Tet(const std::vector<Nodes::Node> &_p_node /**< vector of nodes */,
               const int _idx /**< [in] region index in region vector */,
               std::initializer_list<int> _i /**< [in] node index */)
        : element<N>(_p_node,_idx,_i), idx(0)
        {
        zeroBasing();
        da.setZero();

        if (existNodes())
            {
            orientate();

            Eigen::Matrix3d J;
            // we have to rebuild the jacobian in case of ill oriented tetrahedron
            double detJ = Jacobian(J);

            if (fabs(detJ) < Tetra::epsilon)
                {
                std::cerr << "Singular jacobian in tetrahedron" << std::endl;
                element::infos();
                SYSTEM_ERROR;
                }
            Eigen::Matrix<double,N,Nodes::DIM> dadu;
            dadu << -1., -1., -1., 1., 0., 0., 0., 1., 0., 0., 0., 1.;
            da = dadu * J.inverse();
            }
        }

    /** \return \f$ |J| \f$ build Jacobian \f$ J \f$ */
    double Jacobian(Eigen::Ref<Eigen::Matrix3d> J);

    /** computes volume of the tetrahedron ; unit test Tet_calc_vol */
    double calc_vol(void) const;

    /** idx is the index of the tetrahedron in the vector of tetrahedron */
    int idx;

private:
    /** local hat functions matrix, initialized by constructor: da = dadu * inverse(Jacobian) */
    Eigen::Matrix<double,N,Nodes::DIM> da;

    void orientate(void)
        {
        if (calc_vol() < 0.0)
                std::swap(ind[2], ind[3]);
        }

    };  // end class Tetra

    }   // end namespace Tetra

#endif /* tetra_h */
