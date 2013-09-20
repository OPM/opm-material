// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
 *   Copyright (C) 2011-2012 by Bernd Flemisch                               *
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
 * \copydoc Opm::FluidSystems::TwoPImmiscible
 */
#ifndef OPM_2P_IMMISCIBLE_FLUID_SYSTEM_HH
#define OPM_2P_IMMISCIBLE_FLUID_SYSTEM_HH

#include <limits>
#include <cassert>

#include <opm/material/fluidsystems/liquidphase.hh>
#include <opm/material/fluidsystems/gasphase.hh>
#include <opm/material/fluidstates/immisciblefluidstate.hh>

#include "basefluidsystem.hh"
#include "nullparametercache.hh"

namespace Opm {
namespace FluidSystems {

/*!
 * \ingroup Fluidsystems
 *
 * \brief A fluid system for two-phase models assuming immiscibility and
 *        thermodynamic equilibrium
 *
 * The wetting and the non-wetting phase can be defined individually
 * via <tt>Opm::LiquidPhase<Component></tt> and
 * <tt>Opm::GasPhase<Component></tt>. These phases consist of one pure
 * component. With the help of this adapter class, the phase
 * properties can be accessed. This is suitable for pure two-phase
 * systems without compositional effects.
 */
template <class Scalar, class WettingPhase, class NonwettingPhase>
class TwoPImmiscible
: public BaseFluidSystem<Scalar, TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase> >
{
    // do not try to instanciate this class, it has only static members!
    TwoPImmiscible()
    {}

    typedef TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase> ThisType;
    typedef BaseFluidSystem<Scalar, ThisType> Base;
public:
    //! \copydoc BaseFluidSystem::ParameterCache
    typedef NullParameterCache ParameterCache;

    /****************************************
     * Fluid phase related static parameters
     ****************************************/

    //! \copydoc BaseFluidSystem::numPhases
    static const int numPhases = 2;

    //! Index of the wetting phase
    static const int wPhaseIdx = 0;
    //! Index of the non-wetting phase
    static const int nPhaseIdx = 1;

    //! \copydoc BaseFluidSystem::phaseName
    static const char *phaseName(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        static const char *name[] = {
            "w",
            "n"
        };
        return name[phaseIdx];
    }

    //! \copydoc BaseFluidSystem::isLiquid
    static bool isLiquid(int phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        return
            (phaseIdx == wPhaseIdx)
            ? WettingPhase::isLiquid()
            : NonwettingPhase::isLiquid();
    }

    //! \copydoc BaseFluidSystem::isCompressible
    static bool isCompressible(int phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);

        return
            (phaseIdx == wPhaseIdx)
            ? WettingPhase::isCompressible()
            : NonwettingPhase::isCompressible();
    }

    //! \copydoc BaseFluidSystem::isIdealGas
    static bool isIdealGas(int phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);

        // let the fluids decide
        return
            (phaseIdx == wPhaseIdx)
            ? WettingPhase::isIdealGas()
            : NonwettingPhase::isIdealGas();
    }

    //! \copydoc BaseFluidSystem::isIdealMixture
    static bool isIdealMixture(int phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);

        // we assume immisibility
        return true;
    }

    /****************************************
     * Component related static parameters
     ****************************************/

    //! \copydoc BaseFluidSystem::numComponents
    static const int numComponents = 2;

    //! Index of the wetting phase's component
    static const int wCompIdx = 0;
    //! Index of the non-wetting phase's component
    static const int nCompIdx = 1;

    //! \copydoc BaseFluidSystem::componentName
    static const char *componentName(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        if (compIdx == wCompIdx)
            return WettingPhase::name();
        return NonwettingPhase::name();
    }

    //! \copydoc BaseFluidSystem::molarMass
    static Scalar molarMass(int compIdx)
    {
        //assert(0 <= compIdx && compIdx < numComponents);

        // let the fluids decide
        return
            (compIdx == wCompIdx)
            ? WettingPhase::molarMass()
            : NonwettingPhase::molarMass();
    }

    /*!
     * \brief Critical temperature of a component [K].
     */
    static Scalar criticalTemperature(int compIdx)
    {
        //assert(0 <= compIdx && compIdx < numComponents);
        // let the fluids decide
        return
            (compIdx == wCompIdx)
            ? WettingPhase::criticalTemperature()
            : NonwettingPhase::criticalTemperature();
    }

    /*!
     * \brief Critical pressure of a component [Pa].
     */
    static Scalar criticalPressure(int compIdx)
    {
        //assert(0 <= compIdx && compIdx < numComponents);
        // let the fluids decide
        return
            (compIdx == wCompIdx)
            ? WettingPhase::criticalPressure()
            : NonwettingPhase::criticalPressure();
    }

    /*!
     * \brief The acentric factor of a component [].
     */
    static Scalar acentricFactor(int compIdx)
    {
        //assert(0 <= compIdx && compIdx < numComponents);
        // let the fluids decide
        return
            (compIdx == wCompIdx)
            ? WettingPhase::acentricFactor()
            : NonwettingPhase::acentricFactor();
    }

    /****************************************
     * thermodynamic relations
     ****************************************/

    //! \copydoc BaseFluidSystem::init
    static void init()
    {
        // two gaseous phases at once do not make sense physically!
        // (But two liquids are fine)
        assert(WettingPhase::isLiquid() || NonwettingPhase::isLiquid());
    }

    //! \copydoc BaseFluidSystem::density
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          const ParameterCache &paramCache,
                          int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == wPhaseIdx)
            return WettingPhase::density(temperature, pressure);
        return NonwettingPhase::density(temperature, pressure);
    }

    //! \copydoc BaseFluidSystem::viscosity
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            const ParameterCache &paramCache,
                            int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == wPhaseIdx)
            return WettingPhase::viscosity(temperature, pressure);
        return NonwettingPhase::viscosity(temperature, pressure);
    }

    //! \copydoc BaseFluidSystem::fugacityCoefficient
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      const ParameterCache &paramCache,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        assert(0 <= compIdx  && compIdx < numComponents);

        if (phaseIdx == compIdx)
            // TODO (?): calculate the real fugacity coefficient of
            // the component in the fluid. Probably that's not worth
            // the effort, since the fugacity coefficient of the other
            // component is infinite anyway...
            return 1.0;
        return std::numeric_limits<Scalar>::infinity();
    }

    //! \copydoc BaseFluidSystem::enthalpy
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                           const ParameterCache &paramCache,
                           int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == wPhaseIdx)
            return WettingPhase::enthalpy(temperature, pressure);
        return NonwettingPhase::enthalpy(temperature, pressure);
    }

    //! \copydoc BaseFluidSystem::thermalConductivity
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      const ParameterCache &paramCache,
                                      int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == wPhaseIdx)
            return WettingPhase::thermalConductivity(temperature, pressure);
        return NonwettingPhase::thermalConductivity(temperature, pressure);
    }

    //! \copydoc BaseFluidSystem::heatCapacity
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               const ParameterCache &paramCache,
                               int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == wPhaseIdx)
            return WettingPhase::heatCapacity(temperature, pressure);
        return NonwettingPhase::heatCapacity(temperature, pressure);
    }
};

} // end namepace FluidSystems

} // end namepace

#endif
