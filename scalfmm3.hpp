#pragma once

#include <array>


//
#include "scalfmm/container/particle.hpp"
#include "scalfmm/container/particle_container.hpp"
#include "scalfmm/container/point.hpp"
//
// Tree
#include "scalfmm/tree/box.hpp"
#include "scalfmm/tree/cell.hpp"
#include "scalfmm/tree/for_each.hpp"
#include "scalfmm/tree/group_tree_view.hpp"
#include "scalfmm/tree/leaf_view.hpp"
#include "scalfmm/tree/utils.hpp"
//
#include "scalfmm/meta/utils.hpp"
#include "scalfmm/tools/fma_loader.hpp"
#include "scalfmm/utils/source_target.hpp"
//  Algorithms 
#include "scalfmm/algorithms/fmm.hpp"
#ifdef _OPENMP
#include "scalfmm/algorithms/omp/utils.hpp"
#endif
#include "scalfmm/lists/lists.hpp"

#include "scalfmm/algorithms/sequential/utils.hpp"
#include "scalfmm/matrix_kernels/mk_common.hpp"

#include "scalfmm/interpolation/interpolation.hpp"
#include "scalfmm/operators/fmm_operators.hpp"
//
#include "scalfmm/utils/parameters.hpp"
// 
#include <cpp_tools/cl_parser/help_descriptor.hpp>
#include <cpp_tools/cl_parser/tcli.hpp>
#include <cpp_tools/colors/colorized.hpp>
#include <cpp_tools/timers/simple_timer.hpp>
    //
    ////////////////////////////////////////////////////////
    //
static const int dimension = 2 ;
using value_type = double;

///////
/// \brief The holographic struct corresponds to the  \f$ (d/dr 1/r)^T \f$
/// kernel
///
///   The kernel \f$k(x,y): R^{2} -> R^{1}\f$  with\f$ kn =1 ;  km = 2\f$
///           \f$k(x,y) = ( (x-y)/| x - y |^2 )\f$
///
/// The kernel is homogeneous \f$k( a x,a y) = 1/a k(  x, y) \f$ with coefficient 1/a
///
struct holographic {
  static constexpr auto homogeneity_tag{
      scalfmm::matrix_kernels::homogeneity::homogenous};
  static constexpr auto symmetry_tag{
      scalfmm::matrix_kernels::symmetry::non_symmetric};
  static constexpr std::size_t km{2};
  static constexpr std::size_t kn{1};
  template <typename ValueType>
  using matrix_type = std::array<ValueType, kn * km>;
  template <typename ValueType> using vector_type = std::array<ValueType, kn>;

  const std::string name() const { return std::string("holographic"); }

  // Pas utile
  template <typename ValueType>
  [[nodiscard]] inline constexpr auto mutual_coefficient() const {
            using decayed_type = typename std::decay_t<ValueType>;
            return vector_type<decayed_type>{decayed_type(1.)};
  }
  template <typename ValueType1, typename ValueType2>
  [[nodiscard]] inline auto
  evaluate(scalfmm::container::point<ValueType1, 2> const &x,
           scalfmm::container::point<ValueType2, 2> const &y) const noexcept
      -> std::enable_if_t<
          std::is_same_v<std::decay_t<ValueType1>, std::decay_t<ValueType2>>,
          matrix_type<std::decay_t<ValueType1>>> {
    using decayed_type = typename std::decay_t<ValueType1>;

    auto diff = y - x ;
    decayed_type tmp =
        decayed_type(1.0) / (diff.at(0) * diff.at(0) + diff.at(1) * diff.at(1));
    return matrix_type<decayed_type>{tmp * diff.at(0), tmp * diff.at(1)};
  }

  template <typename ValueType>
  [[nodiscard]] inline auto scale_factor(ValueType cell_width) const noexcept {
    return vector_type<ValueType>{ValueType(1. / cell_width)};
  }

  static constexpr int separation_criterion{1};
};
using options =
      scalfmm::options::chebyshev_<scalfmm::options::low_rank_>;
//    scalfmm::options::uniform_<scalfmm::options::fft_>;
//  near field
using near_matrix_kernel_type = holographic;
using near_field_type =
    scalfmm::operators::near_field_operator<near_matrix_kernel_type>;
// far field
//  Uniform interpolation and fft optimization
using far_matrix_kernel_type = holographic;
using interpolation_type =
    scalfmm::interpolation::interpolator<value_type, dimension,
                                         far_matrix_kernel_type, options>;
using far_field_type =
    scalfmm::operators::far_field_operator<interpolation_type>;
//
using fmm_operators_type =
    scalfmm::operators::fmm_operators<near_field_type, far_field_type>;
//
// 
static constexpr int nb_inputs_source{2};
static constexpr int nb_output_source{1};        // should be zero
static constexpr int nb_input_target{1};         // should be zero
static constexpr int nb_output_target{1};       

using point_type = scalfmm::container::point<value_type, dimension>;
// 
//  source : seulement l'input est important
// point (2d), input (2 valeurs), output (0 valeurs ), variable (indice)
using particle_source_type =
    scalfmm::container::particle<value_type, dimension, value_type,
                                 nb_inputs_source, value_type, nb_output_source,
                                 std::size_t>;
//  target : seulement l'output est important
// point (3d), input (1 valeur), output (1 valeur ), variable (indice)
using particle_target_type =
    scalfmm::container::particle<value_type, dimension, value_type,
                                 nb_input_target, value_type, nb_output_target,
                                 std::size_t>;
// Definiton des feuilles
using leaf_source_type = scalfmm::component::leaf_view<particle_source_type>;
using leaf_target_type = scalfmm::component::leaf_view<particle_target_type>;
// this struct is crucial for having the good type for iterator and group types
// (source) in the target leaves and groups
namespace scalfmm::meta {
template <>
struct inject<scalfmm::component::group_of_particles<
    leaf_target_type, particle_target_type>> {
  using type =
      std::tuple<typename scalfmm::component::group_of_particles<
                     leaf_source_type, particle_source_type>::iterator_type,
                 scalfmm::component::group_of_particles<
                     leaf_source_type, particle_source_type>>;
};
    using box_type = scalfmm::component::box<point_type>;
    //
    using cell_type =
        scalfmm::component::cell<typename interpolation_type::storage_type>;
    using tree_source_type = scalfmm::component::group_tree_view<cell_type, leaf_source_type, box_type>;
    using tree_target_type = scalfmm::component::group_tree_view<cell_type, leaf_target_type, box_type>;
  
} // namespace scalfmm::meta
using box_type = scalfmm::component::box<point_type>;
//
using cell_type =
    scalfmm::component::cell<typename interpolation_type::storage_type>;
using tree_source_type =
    scalfmm::component::group_tree_view<cell_type, leaf_source_type, box_type>;
using tree_target_type =
    scalfmm::component::group_tree_view<cell_type, leaf_target_type, box_type>;


