// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
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
 * \copydoc Opm::ParameterCacheBase
 */
#ifndef OPM_PARAMETER_CACHE_BASE_HH
#define OPM_PARAMETER_CACHE_BASE_HH

namespace Opm {

/*!
 * \ingroup Fluidsystems
 * \brief The base class of the parameter caches of fluid systems
 */
template <class Implementation>
class ParameterCacheBase
{
public:
    /*!
     * \brief Constants for ORing the quantities of the fluid state that have not changed since the last update.
     */
    enum ExceptQuantities {
        None = 0, //!< All quantities have been (potentially) modified.
        Temperature = 1, //< The temperature has not been modified
        Pressure = 2, //< The pressures have not been modified
        Composition = 4  //< The compositions have not been modified
    };

    ParameterCacheBase()
    {}

    /*!
     * \brief Update the quantities of the parameter cache for all phases
     *
     * \param fluidState The representation of the thermodynamic system of interest.
     * \param exceptQuantities The quantities of the fluid state that have not changed since the last update.
     */
    template <class FluidState>
    void updateAll(const FluidState &fluidState, int exceptQuantities = None)
    {
        for (int phaseIdx = 0; phaseIdx < FluidState::numPhases; ++phaseIdx)
            asImp_().updatePhase(fluidState, phaseIdx);
    }

    /*!
     * \brief Update pressure dependent quantities of the parameter
     *        cache for all phases
     *
     * This method should be called if _only_ the phase pressures
     * changed since the last call to an update() method.
     *
     * \param fluidState The representation of the thermodynamic system of interest.
     */
    template <class FluidState>
    void updateAllPressures(const FluidState &fluidState)
    { asImp_().updateAll(fluidState, Temperature | Composition); }

    /*!
     * \brief Update temperature dependent quantities of the parameter
     *        cache for all phases
     *
     * This method should be called if _only_ the phase temperatures
     * changed since the last call to an update() method.
     *
     * \param fluidState The representation of the thermodynamic system of interest.
     */
    template <class FluidState>
    void updateAllTemperatures(const FluidState &fluidState)
    {
        for (int phaseIdx = 0; phaseIdx < FluidState::numPhases; ++phaseIdx)
            asImp_().updatePhase(fluidState, phaseIdx);
    }


    /*!
     * \brief Update all cached parameters of a specific fluid phase.
     *
     * \param fluidState The representation of the thermodynamic system of interest.
     * \param phaseIdx The index of the fluid phase of interest.
     * \param exceptQuantities The quantities of the fluid state that have not changed since the last update.
     */
    template <class FluidState>
    void updatePhase(const FluidState &fluidState, int phaseIdx, int exceptQuantities = None)
    {}

    /*!
     * \brief Update all cached parameters of a specific fluid phase
     *        which depend on temperature
     *
     * *Only* use this method if only the temperature of a phase
     * changed between two update*() calls. If more changed, call
     * updatePhase()!
     *
     * \param fluidState The representation of the thermodynamic system of interest.
     * \param phaseIdx The index of the fluid phase of interest.
     */
    template <class FluidState>
    void updateTemperature(const FluidState &fluidState, int phaseIdx)
    {
        asImp_().updatePhase(fluidState, phaseIdx);
    }

    /*!
     * \brief Update all cached parameters of a specific fluid phase
     *        which depend on pressure
     *
     * *Only* use this method if only the pressure of a phase changed
     * between two update*() calls. If more changed, call
     * updatePhase()!
     *
     * \param fluidState The representation of the thermodynamic system of interest.
     * \param phaseIdx The index of the fluid phase of interest.
     */
    template <class FluidState>
    void updatePressure(const FluidState &fluidState, int phaseIdx)
    {
        asImp_().updatePhase(fluidState, phaseIdx);
    }

    /*!
     * \brief Update all cached parameters of a specific fluid phase
     *        which depend on composition
     *
     * *Only* use this method if neither the pressure nor the
     * temperature of the phase changed between two update*()
     * calls. If more changed, call updatePhase()!
     *
     * \param fluidState The representation of the thermodynamic system of interest.
     * \param phaseIdx The index of the fluid phase of interest.
     */
    template <class FluidState>
    void updateComposition(const FluidState &fluidState, int phaseIdx)
    {
        asImp_().updatePhase(fluidState, phaseIdx, /*except=*/Temperature | Pressure);
    }

    /*!
     * \brief Update all cached parameters of a specific fluid phase
     *        which depend on the mole fraction of a single component
     *
     * *Only* use this method if just a single component's
     * concentration changed between two update*() calls. If more than
     * one concentration changed, call updatePhaseComposition() of
     * updatePhase()!
     *
     * \param fluidState The representation of the thermodynamic system of interest.
     * \param phaseIdx The index of the fluid phase of interest.
     * \param compIdx The component index of the component for which the mole fraction was modified in the fluid phase of interest.
     */
    template <class FluidState>
    void updateSingleMoleFraction(const FluidState &fluidState,
                                  int phaseIdx,
                                  int compIdx)
    {
        asImp_().updateComposition(fluidState, phaseIdx);
    }

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
};

} // end namepace

#endif
