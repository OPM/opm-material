// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
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
 * \copydoc Opm::ModularFluidState
 */
#ifndef OPM_MODULAR_FLUID_STATE_HH
#define OPM_MODULAR_FLUID_STATE_HH

#include "FluidStatePressureModules.hpp"
#include "FluidStateTemperatureModules.hpp"
#include "FluidStateCompositionModules.hpp"
#include "FluidStateFugacityModules.hpp"
#include "FluidStateSaturationModules.hpp"
#include "FluidStateDensityModules.hpp"
#include "FluidStateViscosityModules.hpp"
#include "FluidStateEnthalpyModules.hpp"

#include <opm/material/Valgrind.hpp>
#include <algorithm>

namespace Opm {

/*!
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase, multi-component fluid system assuming
 *        thermodynamic equilibrium.
 *
 * This class uses a modular approach which results in storing only a
 * set of requested thermodynamic quantities.
 */
template <class ScalarT,
          class FluidSystem,
          class PressureModule,
          class TemperatureModule,
          class CompositionModule,
          class FugacityModule,
          class SaturationModule,
          class DensityModule,
          class ViscosityModule,
          class EnthalpyModule>
class ModularFluidState
    : public PressureModule
    , public TemperatureModule
    , public CompositionModule
    , public FugacityModule
    , public SaturationModule
    , public DensityModule
    , public ViscosityModule
    , public EnthalpyModule
{
public:
    typedef ScalarT Scalar;
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };

    /*!
     * \brief Make sure that all attributes are defined.
     *
     * This method does not do anything if the program is not run
     * under valgrind. If it is, then valgrind will print an error
     * message if some attributes of the object have not been properly
     * defined.
     */
    void checkDefined() const
    {
        PressureModule::checkDefined();
        TemperatureModule::checkDefined();
        CompositionModule::checkDefined();
        SaturationModule::checkDefined();
        DensityModule::checkDefined();
        ViscosityModule::checkDefined();
        EnthalpyModule::checkDefined();
    }

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     */
    template <class FluidState>
    void assign(const FluidState& fs)
    {
        PressureModule::assign(fs);
        TemperatureModule::assign(fs);
        CompositionModule::assign(fs);
        SaturationModule::assign(fs);
        DensityModule::assign(fs);
        ViscosityModule::assign(fs);
        EnthalpyModule::assign(fs);
    }
};

} // end namespace Opm

#endif
