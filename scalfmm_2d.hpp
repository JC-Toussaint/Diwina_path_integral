#pragma once

#include <vector>
#include <iostream>

#include "scalfmm/tree/for_each.hpp"
#include "scalfmm/utils/accurater.hpp"
#include <cpp_tools/colors/colorized.hpp>
#include <cpp_tools/timers/simple_timer.hpp>
#include "scalfmm3.hpp"
#include "scalfmm/algorithms/full_direct.hpp"
#include <cpp_tools/colors/colorized.hpp>
#include <cpp_tools/timers/simple_timer.hpp>

template<typename Tree, typename Container>
auto check_output(Container const& part, Tree const& tree)
{
    scalfmm::utils::accurater<
      typename scalfmm::container::particle_traits<typename Container::value_type>::outputs_value_type>
      error;
    static constexpr std::size_t nb_out = Container::value_type::outputs_size;
    int n{0};
    scalfmm::component::for_each_leaf(std::cbegin(tree), std::cend(tree),
                                      [&part, &error,&n](auto& leaf)
                                      {
                                          for(auto const p_tuple_ref: leaf)
                                          {
                                              const auto& p = typename Tree::leaf_type::const_proxy_type(p_tuple_ref);
					      //                                              const auto& idx = std::get<0>(p.variables());
                                              auto& output = p.outputs();
                                              const auto output_ref = part.at(n).outputs();
                                              for(int i{0}; i < nb_out; ++i)
                                              {
						//						std::cout <<n << " dir " <<  part.at(n) << " tree  " << p <<std::endl;
                                                  error.add(output_ref.at(i), output.at(i));
                                              }
					      ++n;
                                          }
                                      });

    return error;
}

template <typename T>
void scalfmm_init_tree(int &ierr, int &iprec, std::vector<T> const &source,
                     std::vector<T> &dipstr, std::vector<T> &dipvec,
                     std::vector<T> &target, std::vector<T> &pottarg) {


  std::cout << "ToDo \n";
};

template <typename T>
void scalfmm_execute(int &ierr, int &iprec, tree_source_type &tree_source,
                     tree_target_type &tree_target, std::vector<T> &pottarg) {

  far_matrix_kernel_type mk_far{};
  auto box_width = tree_target.box().width(0);
  auto tree_height = tree_target.height();
  auto order = tree_target.order();
  value_type extend_cell{0.0};

  interpolation_type interpolator(mk_far, order, tree_height, box_width,extend_cell);

   typename fmm_operators_type::far_field_type far_field(interpolator);
  // Near field
  near_matrix_kernel_type mk_near{};
  typename fmm_operators_type::near_field_type near_field(  mk_near, false /* no mutual for source!= target*/);
  // run

  std::cout << cpp_tools::colors::blue << "Fmm with kernels: " << std::endl
            << "       near      " << mk_near.name() << std::endl
            << "       far       " << mk_far.name() << std::endl
            // << "          option " << op.value() << std::endl
            << cpp_tools::colors::reset;
  ;

   fmm_operators_type fmm_operator(near_field, far_field);
   auto neighbour_separation = fmm_operator.near_field().separation_criterion();
   //  auto neighbour_separation = near_field.separation_criterion();
  //
   tree_target.build_interaction_lists(tree_source, neighbour_separation, false);

   auto operator_to_proceed = scalfmm::algorithms::all;
   std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
   std::chrono::duration<double,std::milli> milli;

    start = std::chrono::high_resolution_clock::now();
#ifdef SCAL_USE_SEQ
    scalfmm::algorithms::fmm[scalfmm::options::_s(scalfmm::options::seq_timit)](tree_source, tree_target, fmm_operator,
                                                                                operator_to_proceed);
#else

    scalfmm::algorithms::fmm[scalfmm::options::_s(scalfmm::options::omp)](tree_source, tree_target, fmm_operator,
                                                                                operator_to_proceed);
#endif
    end = std::chrono::high_resolution_clock::now();
    milli = end-start;
    std::cout << std::endl << boost::format("%5t fmm sum %50T. ");
    std::cout << milli.count() << " [ms]\n";

   // set pot

    start = std::chrono::high_resolution_clock::now();
    scalfmm::component::for_each_leaf(
       std::cbegin(tree_target), std::cend(tree_target),
       [&pottarg](auto const &leaf) {
         for (auto const p_tuple_ref : leaf) {
           const auto p =
               typename tree_target_type::leaf_type::const_proxy_type(
                   p_tuple_ref);

           auto const &out = p.outputs();
           const auto &idx = std::get<0>(p.variables());
           pottarg[idx] = out[0];
         }
       });
    end = std::chrono::high_resolution_clock::now();
    milli = end-start;
    std::cout << std::endl << boost::format("%5t get pot %50T. ");
    std::cout << milli.count() << " [ms]\n";


}

template <typename TS, typename TT>
void scalfmm_compare(std::vector<TS> const &source, std::vector<TT> &target,
                     tree_target_type &tree_target) 

       {
    cpp_tools::timers::timer<std::chrono::minutes> time{};
	 near_matrix_kernel_type mk_near{};

	 std::cout << cpp_tools::colors::green << "full interaction computation  with kernel: " << mk_near.name()
                  << std::endl
                  << cpp_tools::colors::reset;

        time.tic();

        scalfmm::algorithms::full_direct(source, target, mk_near);
        time.tac();
        std::cout << cpp_tools::colors::green << "... Done.\n" << cpp_tools::colors::reset;
        std::cout << cpp_tools::colors::yellow << "Computation done in " << time.elapsed() << " min\n"
                  << cpp_tools::colors::reset;
        // check the two containers
        // std::cout << "Final target container\n";
        // std::cout << container_target << std::endl;
        // Compare with the FMM computation
        auto error{check_output(target, tree_target).get_relative_l2_norm()};
        std::cout << cpp_tools::colors::magenta << "relative L2 error: " << error << '\n' << cpp_tools::colors::reset;
  
    }

