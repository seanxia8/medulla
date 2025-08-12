/**
 * @file variables_nue.h    
 * @brief Variables for the nue analysis.
 * @details This file contains definitions of analysis variables which can be used
 * to extract information from interactions specific to the nue analysis. Each
 * variable is implemented as a function which takes an interaction object as an
 * argument and returns a double. These are the building blocks for producing
 * high-level plots of the selected interactions.
*/
#ifndef VARS_NUE_H
#define VARS_NUE_H

#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRInteractionDLP.h"
#include "sbnanaobj/StandardRecord/SRInteractionTruthDLP.h"

#include "include/utilities.h"
#include "include/cuts.h"
// #include "include/bnb_nue/cuts_nue.h"


/**
 * @namespace vars::nue
 * @brief Namespace for organizing variables specific to the nue analysis.
 * @details This namespace is intended to be used for organizing variables which
 * act on interactions specific to the nue analysis. Each variable is
 * implemented as a function which takes an interaction object as an argument
 * and returns a double. The function should be templated on the type of
 * interaction object if the variable is intended to be used on both true and
 * reconstructed interactions.
 * @note The namespace is intended to be used in conjunction with the vars
 * namespace, which is used for organizing generic variables which act on
 * interactions.
 */
namespace vars::nue
{
    template<class T>
        double leading_shower_ke(const T & obj)
        {
            double leading_ke(0);
            auto & ke(obj.particles[utilities::leading_primary_shower_index(obj)].ke);
            return ke;
        }

    template<class T>
        double leading_shower_start_dedx(const T & obj)
        {
            double leading_dedx(0);
            auto & dedx(obj.particles[utilities::leading_primary_shower_index(obj)].start_dedx);
            return dedx;
        }
    template<class T>
        double leading_shower_vertex_distance(const T & obj)
        {
            double leading_distance(0);
            auto & distance(obj.particles[utilities::leading_primary_shower_index(obj)].vertex_distance);
            return distance;
        }
    template<class T>
        double leading_shower_directional_spread(const T & obj)
        {
            double leading_spread(0);
            auto & spread(obj.particles[utilities::leading_primary_shower_index(obj)].directional_spread);
            return spread;
        }
    template<class T>
        double leading_shower_axial_spread(const T & obj)
        {
            double leading_spread(0);
            auto & spread(obj.particles[utilities::leading_primary_shower_index(obj)].axial_spread);
            return spread;
        }
    template<class T>
        double leading_shower_start_straightness(const T & obj)
        {
            double leading_straightness(0);
            auto & straightness(obj.particles[utilities::leading_primary_shower_index(obj)].start_straightness);
            return straightness;
        }
    template<class T>
        double opening_angle(const T & obj)
        {
            auto & m(obj.particles[utilities::leading_primary_particle_index(obj, 1)]);
            auto & p(obj.particles[utilities::leading_primary_particle_index(obj, 4)]);
            return std::acos(m.start_dir[0] * p.start_dir[0] + m.start_dir[1] * p.start_dir[1] + m.start_dir[2] * p.start_dir[2]);
        }
}
#endif // VARS_NUE_H