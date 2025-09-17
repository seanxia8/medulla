/**
 * @file variables_nue.h    
 * @brief Variables for the nue analysis.
 * @details This file contains definitions of analysis variables which can be used
 * to extract information from interactions specific to the nue analysis. Each
 * variable is implemented as a function which takes an interaction object as an
 * argument and returns a double. These are the building blocks for producing
 * high-level plots of the selected interactions.
*/
#ifndef VARIABLES_NUE_H
#define VARIABLES_NUE_H

#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRInteractionDLP.h"
#include "sbnanaobj/StandardRecord/SRInteractionTruthDLP.h"
#include "sbnanaobj/StandardRecord/SRParticleDLP.h"

#include "include/utilities.h"
#include "include/cuts.h"
#include "include/nue2025/cuts_nue2025.h"


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
namespace vars::nue2025
{
    template<class T>
        double leading_shower_ke(const T & p)
        {
            return p.ke;
        }
        REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, leading_shower_ke, leading_shower_ke);

    template<class T>
        double leading_shower_start_dedx(const T & p)
        {
            return p.start_dedx;
        }
        REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, leading_shower_start_dedx, leading_shower_start_dedx);

    template<class T>
        double leading_shower_vertex_distance(const T & p)
        {
            return p.vertex_distance;
        }
        REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, leading_shower_vertex_distance, leading_shower_vertex_distance);

    template<class T>
        double leading_shower_directional_spread(const T & p)
        {
            return p.directional_spread;
        }
        REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, leading_shower_directional_spread, leading_shower_directional_spread);

    template<class T>
        double leading_shower_axial_spread(const T & p)
        {
            return p.axial_spread;
        }
        REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, leading_shower_axial_spread, leading_shower_axial_spread);

    template<class T>
        double leading_shower_start_straightness(const T & p)
        {
            return p.start_straightness;
        }
        REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, leading_shower_start_straightness, leading_shower_start_straightness);
 }
#endif // VARS_NUE_H