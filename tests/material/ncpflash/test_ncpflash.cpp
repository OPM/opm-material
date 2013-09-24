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
 *
 * \brief This is a program to test the flash calculation which uses
 *        non-linear complementarity problems (NCP)
 *
 * A flash calculation determines the pressures, saturations and
 * composition of all phases given the total mass (or, as in this case
 * the total number of moles) in a given amount of pore space.
 */
#include "config.h"

#include <opm/material/constraintsolvers/MiscibleMultiPhaseComposition.hpp>
#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>
#include <opm/material/constraintsolvers/NcpFlash.hpp>

#include <opm/material/fluidstates/CompositionalFluidState.hpp>

#include <opm/material/fluidsystems/H2ON2FluidSystem.hpp>

#include <opm/material/fluidmatrixinteractions/mp/MpLinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/mp/2pAdapter.hpp>
#include <opm/material/fluidmatrixinteractions/2p/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/2p/RegularizedLinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/2p/RegularizedBrooksCorey.hpp>
#include <opm/material/fluidmatrixinteractions/2p/EffToAbsLaw.hpp>

template <class Scalar, class FluidState>
void checkSame(const FluidState &fsRef, const FluidState &fsFlash)
{
    enum { numPhases = FluidState::numPhases };
    enum { numComponents = FluidState::numComponents };

    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
        Scalar error;

        // check the pressures
        error = 1 - fsRef.pressure(phaseIdx)/fsFlash.pressure(phaseIdx);
        if (std::abs(error) > 1e-6) {
            std::cout << "pressure error phase " << phaseIdx << ": "
                      << fsFlash.pressure(phaseIdx)  << " flash vs "
                      << fsRef.pressure(phaseIdx) << " reference"
                      << " error=" << error << "\n";
        }

        // check the saturations
        error = fsRef.saturation(phaseIdx) - fsFlash.saturation(phaseIdx);
        if (std::abs(error) > 1e-6)
            std::cout << "saturation error phase " << phaseIdx << ": "
                      << fsFlash.saturation(phaseIdx) << " flash vs "
                      << fsRef.saturation(phaseIdx) << " reference"
                      << " error=" << error << "\n";

        // check the compositions
        for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
            error = fsRef.moleFraction(phaseIdx, compIdx) - fsFlash.moleFraction(phaseIdx, compIdx);
            if (std::abs(error) > 1e-6)
                std::cout << "composition error phase " << phaseIdx << ", component " << compIdx << ": "
                          << fsFlash.moleFraction(phaseIdx, compIdx) << " flash vs "
                          << fsRef.moleFraction(phaseIdx, compIdx) << " reference"
                          << " error=" << error << "\n";
        }
    }
}

template <class Scalar, class FluidSystem, class MaterialLaw, class FluidState>
void checkNcpFlash(const FluidState &fsRef,
                   typename MaterialLaw::Params &matParams)
{
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;

    // calculate the total amount of stuff in the reference fluid
    // phase
    ComponentVector globalMolarities(0.0);
    for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            globalMolarities[compIdx] +=
                fsRef.saturation(phaseIdx)*fsRef.molarity(phaseIdx, compIdx);
        }
    }

    // initialize the fluid state for the flash calculation
    typedef Opm::NcpFlash<Scalar, FluidSystem> NcpFlash;
    FluidState fsFlash;

    fsFlash.setTemperature(fsRef.temperature(/*phaseIdx=*/0));

    // run the flash calculation
    typename FluidSystem::ParameterCache paramCache;
    NcpFlash::guessInitial(fsFlash, paramCache, globalMolarities);
    NcpFlash::template solve<MaterialLaw>(fsFlash, paramCache, matParams, globalMolarities);

    // compare the "flashed" fluid state with the reference one
    checkSame<Scalar>(fsRef, fsFlash);
}


template <class Scalar, class FluidSystem, class MaterialLaw, class FluidState>
void completeReferenceFluidState(FluidState &fs,
                                 typename MaterialLaw::Params &matParams,
                                 int refPhaseIdx)
{
    enum { numPhases = FluidSystem::numPhases };

    typedef Opm::ComputeFromReferencePhase<Scalar, FluidSystem> ComputeFromReferencePhase;
    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;

    int otherPhaseIdx = 1 - refPhaseIdx;

    // calculate the other saturation
    fs.setSaturation(otherPhaseIdx, 1.0 - fs.saturation(refPhaseIdx));

    // calulate the capillary pressure
    PhaseVector pC;
    MaterialLaw::capillaryPressures(pC, matParams, fs);
    fs.setPressure(otherPhaseIdx,
                   fs.pressure(refPhaseIdx)
                   + (pC[otherPhaseIdx] - pC[refPhaseIdx]));

    // make the fluid state consistent with local thermodynamic
    // equilibrium
    typename FluidSystem::ParameterCache paramCache;
    ComputeFromReferencePhase::solve(fs,
                                     paramCache,
                                     refPhaseIdx,
                                     /*setViscosity=*/false,
                                     /*setEnthalpy=*/false);
}


int main()
{
    typedef double Scalar;
    typedef Opm::FluidSystems::H2ON2<Scalar, false> FluidSystem;
    typedef Opm::CompositionalFluidState<Scalar, FluidSystem> CompositionalFluidState;

    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    enum { lPhaseIdx = FluidSystem::lPhaseIdx };
    enum { gPhaseIdx = FluidSystem::gPhaseIdx };

    enum { H2OIdx = FluidSystem::H2OIdx };
    enum { N2Idx = FluidSystem::N2Idx };

    typedef Opm::RegularizedBrooksCorey<Scalar> EffMaterialLaw;
    typedef Opm::EffToAbsLaw<EffMaterialLaw> TwoPMaterialLaw;
    typedef Opm::TwoPAdapter<lPhaseIdx, TwoPMaterialLaw> MaterialLaw;
    typedef MaterialLaw::Params MaterialLawParams;

    Scalar T = 273.15 + 25;

    // initialize the tables of the fluid system
    Scalar Tmin = T - 1.0;
    Scalar Tmax = T + 1.0;
    int nT = 3;

    Scalar pmin = 0.0;
    Scalar pmax = 1.25 * 2e6;
    int np = 100;

    FluidSystem::init(Tmin, Tmax, nT, pmin, pmax, np);

    // set the parameters for the capillary pressure law
    MaterialLawParams matParams;
    matParams.setSwr(0.0);
    matParams.setSnr(0.0);
    matParams.setPe(0);
    matParams.setLambda(2.0);

    CompositionalFluidState fsRef;

    // create an fluid state which is consistent

    // set the fluid temperatures
    fsRef.setTemperature(T);

    ////////////////
    // only liquid
    ////////////////
    std::cout << "testing single-phase liquid\n";

    // set liquid saturation
    fsRef.setSaturation(lPhaseIdx, 1.0);

    // set pressure of the liquid phase
    fsRef.setPressure(lPhaseIdx, 2e5);

    // set the liquid composition to pure water
    fsRef.setMoleFraction(lPhaseIdx, N2Idx, 0.0);
    fsRef.setMoleFraction(lPhaseIdx, H2OIdx, 1.0 - fsRef.moleFraction(lPhaseIdx, N2Idx));

    // "complete" the fluid state
    completeReferenceFluidState<Scalar, FluidSystem, MaterialLaw>(fsRef, matParams, lPhaseIdx);

    // check the flash calculation
    checkNcpFlash<Scalar, FluidSystem, MaterialLaw>(fsRef, matParams);

    ////////////////
    // only gas
    ////////////////
    std::cout << "testing single-phase gas\n";
    // set gas saturation
    fsRef.setSaturation(gPhaseIdx, 1.0);

    // set pressure of the gas phase
    fsRef.setPressure(gPhaseIdx, 1e6);

    // set the gas composition to 99.9% nitrogen and 0.1% water
    fsRef.setMoleFraction(gPhaseIdx, N2Idx, 0.999);
    fsRef.setMoleFraction(gPhaseIdx, H2OIdx, 0.001);

    // "complete" the fluid state
    completeReferenceFluidState<Scalar, FluidSystem, MaterialLaw>(fsRef, matParams, gPhaseIdx);

    // check the flash calculation
    checkNcpFlash<Scalar, FluidSystem, MaterialLaw>(fsRef, matParams);

    ////////////////
    // both phases
    ////////////////
    std::cout << "testing two-phase\n";

    // set saturations
    fsRef.setSaturation(lPhaseIdx, 0.5);
    fsRef.setSaturation(gPhaseIdx, 0.5);

    // set pressures
    fsRef.setPressure(lPhaseIdx, 1e6);
    fsRef.setPressure(gPhaseIdx, 1e6);

    FluidSystem::ParameterCache paramCache;
    typedef Opm::MiscibleMultiPhaseComposition<Scalar, FluidSystem> MiscibleMultiPhaseComposition;
    MiscibleMultiPhaseComposition::solve(fsRef, paramCache,
                                         /*setViscosity=*/false,
                                         /*setEnthalpy=*/false);

    // check the flash calculation
    checkNcpFlash<Scalar, FluidSystem, MaterialLaw>(fsRef, matParams);

    ////////////////
    // with capillary pressure
    ////////////////

    MaterialLawParams matParams2;
    matParams2.setSwr(0.0);
    matParams2.setSnr(0.0);
    matParams2.setPe(1e3);
    matParams2.setLambda(2.0);

    // set gas saturation
    fsRef.setSaturation(gPhaseIdx, 0.5);
    fsRef.setSaturation(lPhaseIdx, 0.5);

    // set pressure of the liquid phase
    fsRef.setPressure(lPhaseIdx, 1e6);

    // calulate the capillary pressure
    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;
    PhaseVector pC;
    MaterialLaw::capillaryPressures(pC, matParams2, fsRef);
    fsRef.setPressure(gPhaseIdx,
                      fsRef.pressure(lPhaseIdx)
                      + (pC[gPhaseIdx] - pC[lPhaseIdx]));

    typedef Opm::MiscibleMultiPhaseComposition<Scalar, FluidSystem> MiscibleMultiPhaseComposition;
    MiscibleMultiPhaseComposition::solve(fsRef, paramCache,
                                         /*setViscosity=*/false,
                                         /*setEnthalpy=*/false);


    // check the flash calculation
    checkNcpFlash<Scalar, FluidSystem, MaterialLaw>(fsRef, matParams2);

    return 0;
}
