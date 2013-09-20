// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012-2013 by Andreas Lauser                               *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \copydoc Opm::FluidSystems::Spe5
 */
#ifndef OPM_SPE5_FLUID_SYSTEM_HH
#define OPM_SPE5_FLUID_SYSTEM_HH

#include "basefluidsystem.hh"
#include "spe5parametercache.hh"

#include <opm/core/utility/Spline.hpp>
#include <opm/material/constants.hh>
#include <opm/material/eos/pengrobinsonmixture.hh>

namespace Opm {
namespace FluidSystems {
/*!
 * \ingroup Fluidsystems
 * \brief The fluid system for the oil, gas and water phases of the
 *        SPE5 problem.
 *
 * This problem comprises \f$H_2O\f$, \f$C_1\f$, \f$C_3\f$, \f$C_6\f$,
 * \f$C_10\f$, \f$C_15\f$ and \f$C_20\f$ as components.
 *
 * See:
 *
 * J.E. Killough, et al.: Fifth Comparative Solution Project:
 * Evaluation of Miscible Flood Simulators, Ninth SPE Symposium on
 * Reservoir Simulation, 1987
 */
template <class Scalar>
class Spe5
    : public BaseFluidSystem<Scalar, Spe5<Scalar> >
{
    typedef Opm::FluidSystems::Spe5<Scalar> ThisType;

    typedef typename Opm::PengRobinsonMixture<Scalar, ThisType> PengRobinsonMixture;
    typedef typename Opm::PengRobinson<Scalar> PengRobinson;

    static const Scalar R;

public:
    //! \copydoc BaseFluidSystem::ParameterCache
    typedef Opm::Spe5ParameterCache<Scalar, ThisType> ParameterCache;

    /****************************************
     * Fluid phase parameters
     ****************************************/

    //! \copydoc BaseFluidSystem::numPhases
    static const int numPhases = 3;

    //! Index of the gas phase
    static const int gPhaseIdx = 0;
    //! Index of the water phase
    static const int wPhaseIdx = 1;
    //! Index of the oil phase
    static const int oPhaseIdx = 2;

    //! The component for pure water to be used
    typedef Opm::H2O<Scalar> H2O;

    //! \copydoc BaseFluidSystem::phaseName
    static const char *phaseName(int phaseIdx)
    {
        static const char *name[] = {
            "g",
            "w",
            "o",
        };

        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return name[phaseIdx];
    }

    //! \copydoc BaseFluidSystem::isLiquid
    static bool isLiquid(int phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        return phaseIdx != gPhaseIdx;
    }

    /*!
     * \copydoc BaseFluidSystem::isCompressible
     *
     * In the SPE-5 problems all fluids are compressible...
     */
    static bool isCompressible(int phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        return true;
    }

    //! \copydoc BaseFluidSystem::isIdealGas
    static bool isIdealGas(int phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        return false; // gas is not ideal here!
    }

    //! \copydoc BaseFluidSystem::isIdealMixture
    static bool isIdealMixture(int phaseIdx)
    {
        // always use the reference oil for the fugacity coefficents,
        // so they cannot be dependent on composition and they the
        // phases thus always an ideal mixture
        return phaseIdx == wPhaseIdx;
    }

    /****************************************
     * Component related parameters
     ****************************************/

    //! \copydoc BaseFluidSystem::numComponents
    static const int numComponents = 7;

    static const int H2OIdx = 0; //!< Index of the water component
    static const int C1Idx = 1; //!< Index of the C1 component
    static const int C3Idx = 2; //!< Index of the C3 component
    static const int C6Idx = 3; //!< Index of the C6 component
    static const int C10Idx = 4; //!< Index of the C10 component
    static const int C15Idx = 5; //!< Index of the C15 component
    static const int C20Idx = 6; //!< Index of the C20 component

    //! \copydoc BaseFluidSystem::componentName
    static const char *componentName(int compIdx)
    {
        static const char *name[] = {
            H2O::name(),
            "C1",
            "C3",
            "C6",
            "C10",
            "C15",
            "C20"
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return name[compIdx];
    }

    //! \copydoc BaseFluidSystem::molarMass
    static Scalar molarMass(int compIdx)
    {
        return
            (compIdx == H2OIdx)
            ? H2O::molarMass()
            : (compIdx == C1Idx)
            ? 16.04e-3
            : (compIdx == C3Idx)
            ? 44.10e-3
            : (compIdx == C6Idx)
            ? 86.18e-3
            : (compIdx == C10Idx)
            ? 142.29e-3
            : (compIdx == C15Idx)
            ? 206.00e-3
            : (compIdx == C20Idx)
            ? 282.00e-3
            : 1e100;
    }

    /*!
     * \brief Critical temperature of a component [K].
     */
    static Scalar criticalTemperature(int compIdx)
    {
        return
            (compIdx == H2OIdx)
            ? H2O::criticalTemperature()
            : (compIdx == C1Idx)
            ? 343.0*5/9
            : (compIdx == C3Idx)
            ? 665.7*5/9
            : (compIdx == C6Idx)
            ? 913.4*5/9
            : (compIdx == C10Idx)
            ? 1111.8*5/9
            : (compIdx == C15Idx)
            ? 1270.0*5/9
            : (compIdx == C20Idx)
            ? 1380.0*5/9
            : 1e100;
    }

    /*!
     * \brief Critical pressure of a component [Pa].
     */
    static Scalar criticalPressure(int compIdx)
    {
        return
            (compIdx == H2OIdx)
            ? H2O::criticalPressure()
            : (compIdx == C1Idx)
            ? 667.8 * 6894.7573
            : (compIdx == C3Idx)
            ? 616.3 * 6894.7573
            : (compIdx == C6Idx)
            ? 436.9 * 6894.7573
            : (compIdx == C10Idx)
            ? 304.0 * 6894.7573
            : (compIdx == C15Idx)
            ? 200.0 * 6894.7573
            : (compIdx == C20Idx)
            ? 162.0 * 6894.7573
            : 1e100;
    }

    /*!
     * \brief Molar volume of a component at the critical point [m^3/mol].
     */
    static Scalar criticalMolarVolume(int compIdx)
    {
        return
            (compIdx == H2OIdx)
            ? H2O::criticalMolarVolume()
            : (compIdx == C1Idx)
            ? 0.290*R*criticalTemperature(C1Idx)/criticalPressure(C1Idx)
            : (compIdx == C3Idx)
            ? 0.277*R*criticalTemperature(C3Idx)/criticalPressure(C3Idx)
            : (compIdx == C6Idx)
            ? 0.264*R*criticalTemperature(C6Idx)/criticalPressure(C6Idx)
            : (compIdx == C10Idx)
            ? 0.257*R*criticalTemperature(C10Idx)/criticalPressure(C10Idx)
            : (compIdx == C15Idx)
            ? 0.245*R*criticalTemperature(C15Idx)/criticalPressure(C15Idx)
            : (compIdx == C20Idx)
            ? 0.235*R*criticalTemperature(C20Idx)/criticalPressure(C20Idx)
            : 1e100;
    }

    /*!
     * \brief The acentric factor of a component [].
     */
    static Scalar acentricFactor(int compIdx)
    {
        return
            (compIdx == H2OIdx)
            ? H2O::acentricFactor()
            : (compIdx == C1Idx)
            ? 0.0130
            : (compIdx == C3Idx)
            ? 0.1524
            : (compIdx == C6Idx)
            ? 0.3007
            : (compIdx == C10Idx)
            ? 0.4885
            : (compIdx == C15Idx)
            ? 0.6500
            : (compIdx == C20Idx)
            ? 0.8500
            : 1e100;
    }

    /*!
     * \brief Returns the interaction coefficient for two components.
     *
     * The values are given by the SPE5 paper.
     */
    static Scalar interactionCoefficient(int comp1Idx, int comp2Idx)
    {
        int i = std::min(comp1Idx, comp2Idx);
        int j = std::max(comp1Idx, comp2Idx);
        if (i == C1Idx && (j == C15Idx || j == C20Idx))
            return 0.05;
        else if (i == C3Idx && (j == C15Idx || j == C20Idx))
            return 0.005;
        return 0;
    }

    /****************************************
     * Methods which compute stuff
     ****************************************/

    /*!
     * \brief \copydoc BaseFluidSystem::init
     *
     * \param minT The minimum temperature possibly encountered during the simulation
     * \param maxT The maximum temperature possibly encountered during the simulation
     * \param minP The minimum pressure possibly encountered during the simulation
     * \param maxP The maximum pressure possibly encountered during the simulation
     */
    static void init(Scalar minT = 273.15,
                     Scalar maxT = 373.15,
                     Scalar minP = 1e4,
                     Scalar maxP = 100e6)
    {
        Opm::PengRobinsonParamsMixture<Scalar, ThisType, gPhaseIdx, /*useSpe5=*/true> prParams;

        // find envelopes of the 'a' and 'b' parameters for the range
        // minT <= T <= maxT and minP <= p <= maxP. For
        // this we take advantage of the fact that 'a' and 'b' for
        // mixtures is just a convex combination of the attractive and
        // repulsive parameters of the pure components

        Scalar minA = 1e100, maxA = -1e100;
        Scalar minB = 1e100, maxB = -1e100;

        prParams.updatePure(minT, minP);
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            minA = std::min(prParams.pureParams(compIdx).a(), minA);
            maxA = std::max(prParams.pureParams(compIdx).a(), maxA);
            minB = std::min(prParams.pureParams(compIdx).b(), minB);
            maxB = std::max(prParams.pureParams(compIdx).b(), maxB);
        };

        prParams.updatePure(maxT, minP);
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            minA = std::min(prParams.pureParams(compIdx).a(), minA);
            maxA = std::max(prParams.pureParams(compIdx).a(), maxA);
            minB = std::min(prParams.pureParams(compIdx).b(), minB);
            maxB = std::max(prParams.pureParams(compIdx).b(), maxB);
        };

        prParams.updatePure(minT, maxP);
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            minA = std::min(prParams.pureParams(compIdx).a(), minA);
            maxA = std::max(prParams.pureParams(compIdx).a(), maxA);
            minB = std::min(prParams.pureParams(compIdx).b(), minB);
            maxB = std::max(prParams.pureParams(compIdx).b(), maxB);
        };

        prParams.updatePure(maxT, maxP);
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            minA = std::min(prParams.pureParams(compIdx).a(), minA);
            maxA = std::max(prParams.pureParams(compIdx).a(), maxA);
            minB = std::min(prParams.pureParams(compIdx).b(), minB);
            maxB = std::max(prParams.pureParams(compIdx).b(), maxB);
        };

        PengRobinson::init(/*aMin=*/minA, /*aMax=*/maxA, /*na=*/100,
                           /*bMin=*/minB, /*bMax=*/maxB, /*nb=*/200);
    }

    //! \copydoc BaseFluidSystem::density
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          const ParameterCache &paramCache,
                          int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        return fluidState.averageMolarMass(phaseIdx)/paramCache.molarVolume(phaseIdx);
    }

    //! \copydoc BaseFluidSystem::viscosity
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            const ParameterCache &paramCache,
                            int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx <= numPhases);

        if (phaseIdx == gPhaseIdx) {
            // given by SPE-5 in table on page 64. we use a constant
            // viscosity, though...
            return 0.0170e-2 * 0.1;
        }
        else if (phaseIdx == wPhaseIdx)
            // given by SPE-5: 0.7 centi-Poise  = 0.0007 Pa s
            return 0.7e-2 * 0.1;
        else {
            assert(phaseIdx == oPhaseIdx);
            // given by SPE-5 in table on page 64. we use a constant
            // viscosity, though...
            return 0.208e-2 * 0.1;
        }
    }

    //! \copydoc BaseFluidSystem::fugacityCoefficient
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      const ParameterCache &paramCache,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx <= numPhases);
        assert(0 <= compIdx  && compIdx <= numComponents);

        if (phaseIdx == oPhaseIdx || phaseIdx == gPhaseIdx)
            return PengRobinsonMixture::computeFugacityCoefficient(fluidState,
                                                                   paramCache,
                                                                   phaseIdx,
                                                                   compIdx);
        else {
            assert(phaseIdx == wPhaseIdx);
            return
                henryCoeffWater_(compIdx, fluidState.temperature(wPhaseIdx))
                / fluidState.pressure(wPhaseIdx);
        }
    }

protected:
    static Scalar henryCoeffWater_(int compIdx, Scalar temperature)
    {
        // use henry's law for the solutes and the vapor pressure for
        // the solvent.
        switch (compIdx) {
        case H2OIdx: return H2O::vaporPressure(temperature);

            // the values of the Henry constant for the solutes have
            // are faked so far...
        case C1Idx: return 5.57601e+09;
        case C3Idx: return 1e10;
        case C6Idx: return 1e10;
        case C10Idx: return 1e10;
        case C15Idx: return 1e10;
        case C20Idx: return 1e10;
        default: OPM_THROW(std::logic_error, "Unknown component index " << compIdx);
        }
    }
};

template <class Scalar>
const Scalar Spe5<Scalar>::R = Opm::Constants<Scalar>::R;

} // end namepace
} // end namepace

#endif
