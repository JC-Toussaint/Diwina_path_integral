#ifndef element_h
#define element_h

#include <execution>

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>

#include "node.h"

/** \class element
\brief Template abstract class, mother class for tetraedrons and facettes.

template parameter is N number of sommits. It contains a list of indices to the N nodes of the element, a reference to the full nodes vector, and index refering to the associated material parameters. All indices are zero based, derived class constructor should call zerobasing() if needed.
orientate() is a pure virtual function, it should manipulate indices to orientate positively the element.
*/

template <int N>
class element
    {
    /** constructor */
    public:
    explicit element(const std::vector<Nodes::Node> &_p_node /**< vector of nodes */,
            const int _idx /**< index to params */,
            std::initializer_list<int> & _i /**< indices to the nodes */
            ) : idxPrm(_idx), refNode(_p_node)
        {
        if(_i.size() == N)
            { ind.assign(_i); }
        else
            {
            std::cout<<"Warning: element constructor is given an init list with size() != N\n";
            }
        }

    /** indices to the nodes */
    std::vector<int> ind;
    
    /** index of the material parameters of the element */
    int idxPrm;

    /** getter for N */
    inline constexpr int getN(void) const { return N; }

    /** info: print node indices of the element and the vector index of the associated param */
    void infos() const
        {
        std::cout << "idxPrm: " << idxPrm << " ind: (";
        for(unsigned int i = 0; i < N-1; i++)
            { std::cout << ind[i] << ", "; }
        std::cout << ind[N-1] << ")\n";
        };

    protected:
        /** returns reference to node at ind[i] from mesh node vector */
        inline const Nodes::Node & getNode(const int i) const { return refNode[ind[i]]; }

        /** returns true if mesh node vector is not empty */
        inline bool existNodes(void)
        { return (refNode.size() > 0); }

        /** zeroBasing: index convention Matlab/msh (one based) -> C++ (zero based) */
        inline void zeroBasing(void)
            { std::for_each(ind.begin(),ind.end(),[](int & _i){ _i--; } ); }

    private:
        /** vector of nodes */
        const std::vector<Nodes::Node> & refNode;

        /** a method to orientate the element */
        virtual void orientate() = 0;
    };
    
#endif
