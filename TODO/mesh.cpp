#include "mesh.h"
#include <utility>

using namespace Mesh;

void mesh::infos(void) const
    {
    std::cout << "mesh:\n";
    std::cout << "  nodes:              " << getNbNodes() << '\n';
    std::cout << "  faces:              " << fac.size() << '\n';
    std::cout << "  tetraedrons:        " << tet.size() << '\n';
    std::cout << "  total volume:       " << vol << '\n';
    }

/**
 * @brief Computes the weighted average of a quantity over a region of the mesh.
 *
 * This function loops over all tetrahedral elements (`Tetra::Tet`) and computes
 * a weighted average of a scalar field. The field values are provided by the `getter` function,
 * which takes a node and a coordinate index, and returns the corresponding scalar.
 *
 * The contribution of each tetrahedron is computed using interpolation at its integration points,
 * and the result is weighted using the dot product with predefined weights.
 *
 * @param getter   Function that extracts a scalar value from a node and a coordinate index.
 * @param d        Index of the data/component to extract (e.g., 0 for x, 1 for y...).
 * @param region   Index of the region to consider, or -1 to include all regions.
 *
 * @return         The weighted average of the scalar field over the specified region.
 *
 * @note Example usage — average of the X coordinate over the whole mesh:
 * @code
 * double avgX = mesh.avg([](const Nodes::Node& n, Nodes::index idx) {
 *     return n.p(idx); // Return coordinate at index 'idx'
 * }, 0, -1); // 0 for x, -1 for all regions
 * @endcode
 */
double mesh::avg(std::function<double(Nodes::Node, Nodes::index)> getter,
                 Nodes::index d,
                 int region) const
{
    double sum = std::transform_reduce(
        EXEC_POL, tet.begin(), tet.end(), 0.0, std::plus<>(),
        [getter, &d, region](const Tetra::Tet& te) {
            if (te.idxPrm != region && region != -1)
                return 0.0;

            Eigen::Matrix<double, Tetra::NPI, 1> val;
            te.interpolation(getter, d, val);  // Interpolates the scalar field at integration points
            return te.weight.dot(val);         // Weighted sum contribution
        });

    double volume = (region == -1) ? vol : paramTetra[region].volume;
    return sum / volume;
}


/**
 * @brief Iterates over all mesh nodes and applies a custom predicate to compute a result.
 *
 * This function loops through each node in the mesh, extracts the coordinate at the given index `coord`,
 * and applies the user-defined predicate `whatToDo` to determine whether the current node's value should
 * update the accumulator `result`.
 *
 * It is useful for operations like finding the minimum or maximum value across nodes, or for applying
 * any custom comparison logic.
 *
 * @param init_val   Initial value for the accumulator `result`.
 * @param coord      Coordinate index to extract from each node (e.g., 0 = x, 1 = y, etc.).
 * @param whatToDo   A boolean function taking (val, result) and returning true if `val` should replace `result`.
 *
 * @return           The final value of `result` after applying the predicate to all nodes.
 *
 * @note Example usage — finding the minimum X coordinate:
 * @code
 * double minX = mesh.doOnNodes(DBL_MAX, 0, [](double val, double result) {
 *     return val < result;
 * });
 * @endcode
 */
double mesh::doOnNodes(const double init_val, const Nodes::index coord,
                       std::function<bool(double, double)> whatToDo) const
{
    double result = init_val;

    std::for_each(node.begin(), node.end(),
                  [&result, coord, &whatToDo](const Nodes::Node& n) {
                      double val = n.p(coord);
                      if (whatToDo(val, result)) {
                          result = val;
                      }
                  });

    return result;
}


/**
 * @brief Reorders and updates facet orientation and values based on tetrahedral connectivity.
 *
 * This function performs a two-step process:
 * 
 * 1. **Facet Set Construction**:
 *    It loops through all tetrahedra and generates the set of their faces (`Facette::Fac`),
 *    inserting them into a `std::set` to exploit the custom `operator<` defined in `Fac`.
 * 
 * 2. **Face Matching and Orientation**:
 *    Then, it loops through the existing mesh facets and attempts to match them with those
 *    in the set. If a match is found (after testing cyclic permutations and a swap of vertex order),
 *    it computes the dot product between vectors on the matched triangle and the normal of the facet
 *    to determine the orientation sign. The value `fa.dMs` is then updated using `copysign(...)`
 *    based on the Jacobian and magnetic permeability.
 *
 * @details This is crucial for correctly orienting boundary faces and assigning physical values
 * like surface integrals, normal vectors, or boundary conditions.
 */
void mesh::indexReorder()
{
    std::set<Facette::Fac> sf;  // Set of unique faces from tetrahedra

    // Step 1: Build the set of tetrahedron faces
    std::for_each(tet.begin(), tet.end(), [this, &sf](const Tetra::Tet& te) {
        const int ia = te.ind[0];
        const int ib = te.ind[1];
        const int ic = te.ind[2];
        const int id = te.ind[3];

        // Insert 4 faces of the tetrahedron into the set (with consistent winding)
        sf.insert(Facette::Fac(node, 0, te.idxPrm, {ia, ic, ib}));
        sf.insert(Facette::Fac(node, 0, te.idxPrm, {ib, ic, id}));
        sf.insert(Facette::Fac(node, 0, te.idxPrm, {ia, id, ic}));
        sf.insert(Facette::Fac(node, 0, te.idxPrm, {ia, ib, id}));
    });

    // Step 2: Match and reorient existing mesh facets
    std::for_each(fac.begin(), fac.end(), [this, &sf](Facette::Fac& fa) {
        int i0 = fa.ind[0], i1 = fa.ind[1], i2 = fa.ind[2];
        auto it = sf.end();

        for (int perm = 0; perm < 2; ++perm) {
            for (int nrot = 0; nrot < 3; ++nrot) {
                Facette::Fac fc(node, 0, 0, {0, 0, 0});
                fc.ind[(0 + nrot) % 3] = i0;
                fc.ind[(1 + nrot) % 3] = i1;
                fc.ind[(2 + nrot) % 3] = i2;

                it = sf.find(fc);
                if (it != sf.end()) break;
            }

            if (it != sf.end()) {
                // Found a matching face in the set
                Eigen::Vector3d p0p1 = node[it->ind[1]].p - node[it->ind[0]].p;
                Eigen::Vector3d p0p2 = node[it->ind[2]].p - node[it->ind[0]].p;

                // Compute signed area-related value using the normal direction
                double sign = p0p1.dot(p0p2.cross(fa.calc_norm()));
                fa.dMs = std::copysign(paramTetra[it->idxPrm].J / mu0, sign);
            }

            // Swap i1 and i2 to try the opposite winding
            std::swap(i1, i2);

            // Update the normal vector after possible index reordering
            fa.n = fa.calc_norm();
        }
    });
}

/**
 * @brief Sorts mesh nodes along a specified axis and updates all related connectivity.
 *
 * This function performs a stable sort of the mesh's nodes along the given coordinate axis (`long_axis`)
 * by using a permutation vector. After sorting, it updates:
 * - The `node` array to reflect the new order
 * - The `node_index` lookup vector to map old indices to new ones
 * - All tetrahedron and facet indices that reference nodes, ensuring mesh integrity
 *
 * This kind of operation can be useful for:
 * - Optimizing data layout for cache efficiency
 * - Ensuring reproducibility
 * - Preprocessing for spatial partitioning or hierarchical algorithms
 *
 * @param long_axis The coordinate axis (e.g. 0 = X, 1 = Y, 2 = Z) used for sorting.
 */
void mesh::sortNodes(Nodes::index long_axis)
{
    // Step 1: Create a permutation vector [0, 1, ..., N-1]
    std::vector<int> permutation(node.size());
    std::iota(permutation.begin(), permutation.end(), 0);

    // Step 2: Sort the permutation based on the node coordinates along the given axis
    std::sort(permutation.begin(), permutation.end(),
              [this, long_axis](int a, int b) {
                  return node[a].p(long_axis) < node[b].p(long_axis);
              });

    // Step 3: Build a reverse mapping from old index to new index
    node_index.resize(node.size());
    for (size_t i = 0; i < node.size(); ++i) {
        node_index[permutation[i]] = i;
    }

    // Step 4: Reorder the actual node array according to the sorted permutation
    std::vector<Nodes::Node> node_copy(node);
    for (size_t i = 0; i < node.size(); ++i) {
        node[i] = node_copy[permutation[i]];
    }

    // Step 5: Update node indices in tetrahedra
    std::for_each(tet.begin(), tet.end(),
                  [this](Tetra::Tet& tetrahedron) {
                      tetrahedron.ind[0] = node_index[tetrahedron.ind[0]];
                      tetrahedron.ind[1] = node_index[tetrahedron.ind[1]];
                      tetrahedron.ind[2] = node_index[tetrahedron.ind[2]];
                      tetrahedron.ind[3] = node_index[tetrahedron.ind[3]];
                  });

    // Step 6: Update node indices in facets
    std::for_each(fac.begin(), fac.end(),
                  [this](Facette::Fac& facette) {
                      facette.ind[0] = node_index[facette.ind[0]];
                      facette.ind[1] = node_index[facette.ind[1]];
                      facette.ind[2] = node_index[facette.ind[2]];
                  });
}

