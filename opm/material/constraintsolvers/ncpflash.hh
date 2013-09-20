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
 * \copydoc Opm::NcpFlash
 */
#ifndef OPM_NCP_FLASH_HH
#define OPM_NCP_FLASH_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/Exceptions.hpp>
#include <opm/core/utility/Average.hpp>
#include <opm/material/valgrind.hh>

#include <limits>
#include <iostream>

namespace Opm {

/*!
 * \brief Determines the phase compositions, pressures and saturations
 *        given the total mass of all components.
 *
 * In a M-phase, N-component context, we have the following
 * unknowns:
 *
 * - M pressures
 * - M saturations
 * - M*N mole fractions
 *
 * This sums up to M*(N + 2). On the equations side of things,
 * we have:
 *
 * - (M - 1)*N equation stemming from the fact that the
 *   fugacity of any component is the same in all phases
 * - 1 equation from the closure condition of all saturations
 *   (they sum up to 1)
 * - M - 1 constraints from the capillary pressures
 *   \f$(-> p_\beta = p_\alpha + p_c\alpha,\beta)\f$
 * - N constraints from the fact that the total mass of each
 *   component is given \f$(-> sum_\alpha rhoMolar_\alpha *
 *   x_\alpha^\kappa = const)\f$
 * - M model constraints. Here we use the NCP constraints
 *   (-> 0 = min \f$ {S_\alpha, 1 - \sum_\kappa x_\alpha^\kappa}\f$)
 *
 * this also sums up to M*(N + 2).
 *
 * We use the following catches: Capillary pressures are taken
 * into accout expicitly, so that only the pressure of the first
 * phase is solved implicitly, also the closure condition for the
 * saturations is taken into account explicitly, which means, that
 * we don't need to implicitly solve for the last
 * saturation. These two measures reduce the number of unknowns to
 * M*(N + 1), namely:
 *
 * - 1 pressure
 * - M - 1 saturations
 * - M*N mole fractions
 */
template <class Scalar, class FluidSystem>
class NcpFlash
{
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };

    typedef typename FluidSystem::ParameterCache ParameterCache;

    static const int numEq = numPhases*(numComponents + 1);

    typedef Dune::FieldMatrix<Scalar, numEq, numEq> Matrix;
    typedef Dune::FieldVector<Scalar, numEq> Vector;

public:
    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;

    /*!
     * \brief Guess initial values for all quantities.
     */
    template <class FluidState>
    static void guessInitial(FluidState &fluidState,
                             ParameterCache &paramCache,
                             const ComponentVector &globalMolarities)
    {
        // the sum of all molarities
        Scalar sumMoles = 0;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            sumMoles += globalMolarities[compIdx];

        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            // composition
            for (int compIdx = 0; compIdx < numComponents; ++ compIdx)
                fluidState.setMoleFraction(phaseIdx,
                                           compIdx,
                                           globalMolarities[compIdx]/sumMoles);

            // pressure. use atmospheric pressure as initial guess
            fluidState.setPressure(phaseIdx, 1.0135e5);

            // saturation. assume all fluids to be equally distributed
            fluidState.setSaturation(phaseIdx, 1.0/numPhases);
        }

        // set the fugacity coefficients of all components in all phases
        paramCache.updateAll(fluidState);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
                Scalar phi = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx);
                fluidState.setFugacityCoefficient(phaseIdx, compIdx, phi);
            }
        }
    }

    /*!
     * \brief Calculates the chemical equilibrium from the component
     *        fugacities in a phase.
     *
     * The phase's fugacities must already be set.
     */
    template <class MaterialLaw, class FluidState>
    static void solve(FluidState &fluidState,
                      ParameterCache &paramCache,
                      const typename MaterialLaw::Params &matParams,
                      const ComponentVector &globalMolarities,
                      Scalar tolerance = 0.0)
    {
        Dune::FMatrixPrecision<Scalar>::set_singular_limit(1e-35);

        if (tolerance <= 0.0) {
            tolerance = std::min(1e-10,
                                 Opm::utils::geometricAverage(Scalar(1.0),
                                                              std::numeric_limits<Scalar>::epsilon()));
        }

        /////////////////////////
        // Newton method
        /////////////////////////

        // Jacobian matrix
        Matrix J;
        // solution, i.e. phase composition
        Vector deltaX;
        // right hand side
        Vector b;

        Valgrind::SetUndefined(J);
        Valgrind::SetUndefined(deltaX);
        Valgrind::SetUndefined(b);

        // make the fluid state consistent with the fluid system.
        completeFluidState_<MaterialLaw>(fluidState,
                                         paramCache,
                                         matParams);

        /*
        std::cout << "--------------------\n";
        std::cout << "globalMolarities: ";
        for (int compIdx = 0; compIdx < numComponents; ++ compIdx)
            std::cout << globalMolarities[compIdx] << " ";
        std::cout << "\n";
        */
        const int nMax = 50; // <- maximum number of newton iterations
        for (int nIdx = 0; nIdx < nMax; ++nIdx) {
            // calculate Jacobian matrix and right hand side
            linearize_<MaterialLaw>(J,
                                    b,
                                    fluidState,
                                    paramCache,
                                    matParams,
                                    globalMolarities);
            Valgrind::CheckDefined(J);
            Valgrind::CheckDefined(b);

            // Solve J*x = b
            deltaX = 0;

            try { J.solve(deltaX, b); }
            catch (Dune::FMatrixError e)
            {
                /*
                printFluidState_(fluidState);
                std::cout << "error: " << e << "\n";
                std::cout << "b: " << b << "\n";
                std::cout << "J: " << J << "\n";
                */

                throw Opm::NumericalProblem(e.what());
            }
            Valgrind::CheckDefined(deltaX);

            /*
            printFluidState_(fluidState);
            std::cout << "J:\n";
            for (int i = 0; i < numEq; ++i) {
                for (int j = 0; j < numEq; ++j) {
                    std::ostringstream os;
                    os << J[i][j];

                    std::string s(os.str());
                    do {
                        s += " ";
                    } while (s.size() < 20);
                    std::cout << s;
                }
                std::cout << "\n";
            }

            std::cout << "deltaX: " << deltaX << "\n";
            std::cout << "---------------\n";
            */

            // update the fluid quantities.
            //update_<MaterialLaw>(fluidState, paramCache, matParams, deltaX);
            Scalar relError = update_<MaterialLaw>(fluidState, paramCache, matParams, deltaX);

            if (relError < 1e-9)
                return;
        }

        /*
        printFluidState_(fluidState);
        std::cout << "globalMolarities: ";
        for (int compIdx = 0; compIdx < numComponents; ++ compIdx)
            std::cout << globalMolarities[compIdx] << " ";
        std::cout << "\n";
        */

        OPM_THROW(NumericalProblem,
                  "Flash calculation failed."
                  " {c_alpha^kappa} = {" << globalMolarities << "}, T = "
                  << fluidState.temperature(/*phaseIdx=*/0));
    }


protected:
    template <class FluidState>
    static void printFluidState_(const FluidState &fluidState)
    {
        std::cout << "saturations: ";
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            std::cout << fluidState.saturation(phaseIdx) << " ";
        std::cout << "\n";

        std::cout << "pressures: ";
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            std::cout << fluidState.pressure(phaseIdx) << " ";
        std::cout << "\n";

        std::cout << "densities: ";
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            std::cout << fluidState.density(phaseIdx) << " ";
        std::cout << "\n";

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            std::cout << "composition " << FluidSystem::phaseName(phaseIdx) << "Phase: ";
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                std::cout << fluidState.moleFraction(phaseIdx, compIdx) << " ";
            }
            std::cout << "\n";
        }

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            std::cout << "fugacities " << FluidSystem::phaseName(phaseIdx) << "Phase: ";
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                std::cout << fluidState.fugacity(phaseIdx, compIdx) << " ";
            }
            std::cout << "\n";
        }

        std::cout << "global component molarities: ";
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            Scalar sum = 0;
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                sum += fluidState.saturation(phaseIdx)*fluidState.molarity(phaseIdx, compIdx);
            }
            std::cout << sum << " ";
        }
        std::cout << "\n";
    }

    template <class MaterialLaw, class FluidState>
    static Scalar linearize_(Matrix &J,
                             Vector &b,
                             FluidState &fluidState,
                             ParameterCache &paramCache,
                             const typename MaterialLaw::Params &matParams,
                             const ComponentVector &globalMolarities)
    {
        FluidState origFluidState(fluidState);
        ParameterCache origParamCache(paramCache);

        Vector tmp;

        // reset jacobian
        J = 0;

        Valgrind::SetUndefined(b);
        calculateDefect_(b, fluidState, fluidState, globalMolarities);
        Valgrind::CheckDefined(b);

        ///////
        // calculate the absolute error
        ///////
        Scalar absError = 0;

#if 0
        // fugacities are equal
        int eqIdx = 0;
        for (int compIdx = 0; compIdx < numComponents; ++ compIdx)
            absError = std::max(std::abs(b[eqIdx + compIdx]), absError);
        eqIdx += numComponents;

        // sum of concentrations are given
        for (int compIdx = 0; compIdx < numComponents; ++ compIdx)
            absError = std::max(std::abs(b[eqIdx + compIdx]), absError);
        eqIdx += numComponents;

        // NCP model assumptions
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            absError = std::max(std::abs(b[eqIdx + phaseIdx]), absError);
        eqIdx += numPhases;
#endif
        ///////
        // assemble jacobian matrix
        ///////
        for (int pvIdx = 0; pvIdx < numEq; ++ pvIdx) {
            ////////
            // approximately calculate partial derivatives of the
            // fugacity defect of all components in regard to the mole
            // fraction of the i-th component. This is done via
            // forward differences

            // deviate the mole fraction of the i-th component
            Scalar x_i = getQuantity_(fluidState, pvIdx);
            const Scalar eps = std::numeric_limits<Scalar>::epsilon()*1e7/(quantityWeight_(fluidState, pvIdx));

            setQuantity_<MaterialLaw>(fluidState, paramCache, matParams, pvIdx, x_i + eps);
            assert(std::abs(getQuantity_(fluidState, pvIdx) - (x_i + eps))
                   <= std::max(1.0, std::abs(x_i))*std::numeric_limits<Scalar>::epsilon()*100);

            // compute derivative of the defect
            calculateDefect_(tmp, origFluidState, fluidState, globalMolarities);
            tmp -= b;
            tmp /= eps;

            // store derivative in jacobian matrix
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                J[eqIdx][pvIdx] = tmp[eqIdx];

            // fluid state and parameter cache to their original values
            fluidState = origFluidState;
            paramCache = origParamCache;

            // end forward differences
            ////////
        }

        return absError;
    }

    template <class FluidState>
    static void calculateDefect_(Vector &b,
                                 const FluidState &fluidStateEval,
                                 const FluidState &fluidState,
                                 const ComponentVector &globalMolarities)
    {
        int eqIdx = 0;

        // fugacity of any component must be equal in all phases
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            for (int phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx) {
                b[eqIdx] =
                    fluidState.fugacity(/*phaseIdx=*/0, compIdx)
                    - fluidState.fugacity(phaseIdx, compIdx);
                ++eqIdx;
            }
        }

        assert(eqIdx == numComponents*(numPhases - 1));

        // the fact saturations must sum up to 1 is included explicitly!

        // capillary pressures are explicitly included!

        // global molarities are given
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            b[eqIdx] = 0.0;
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                b[eqIdx] +=
                    fluidState.saturation(phaseIdx)
                    * fluidState.molarity(phaseIdx, compIdx);
            }

            b[eqIdx] -= globalMolarities[compIdx];
            ++eqIdx;
        }

        // model assumptions (-> non-linear complementarity functions)
        // must be adhered
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar sumMoleFracEval = 0.0;
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                sumMoleFracEval += fluidStateEval.moleFraction(phaseIdx, compIdx);

            if (1.0 - sumMoleFracEval > fluidStateEval.saturation(phaseIdx)) {
                b[eqIdx] = fluidState.saturation(phaseIdx);
            }
            else {
                Scalar sumMoleFrac = 0.0;
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                    sumMoleFrac += fluidState.moleFraction(phaseIdx, compIdx);
                b[eqIdx] = 1.0 - sumMoleFrac;
            }

            ++eqIdx;
        }
    }

    template <class MaterialLaw, class FluidState>
    static Scalar update_(FluidState &fluidState,
                          ParameterCache &paramCache,
                          const typename MaterialLaw::Params &matParams,
                          const Vector &deltaX)
    {
        // make sure we don't swallow non-finite update vectors
        assert(std::isfinite(deltaX.two_norm()));

        Scalar relError = 0;
        for (int pvIdx = 0; pvIdx < numEq; ++ pvIdx) {
            Scalar tmp = getQuantity_(fluidState, pvIdx);
            Scalar delta = deltaX[pvIdx];

            relError = std::max(relError, std::abs(delta)*quantityWeight_(fluidState, pvIdx));

            if (isSaturationIdx_(pvIdx)) {
                // dampen to at most 25% change in saturation per
                // iteration
                delta = std::min(0.25, std::max(-0.25, delta));
            }
            else if (isMoleFracIdx_(pvIdx)) {
                // dampen to at most 20% change in mole fraction per
                // iteration
                delta = std::min(0.20, std::max(-0.20, delta));
            }
            else if (isPressureIdx_(pvIdx)) {
                // dampen to at most 50% change in pressure per
                // iteration
                delta = std::min(0.5*fluidState.pressure(0), std::max(-0.5*fluidState.pressure(0), delta));
            }

            setQuantityRaw_(fluidState, pvIdx, tmp - delta);
        }

        /*
        // make sure all saturations, pressures and mole fractions are non-negative
        Scalar sumSat = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar value = fluidState.saturation(phaseIdx);
            if (value < -0.05) {
                value = -0.05;
                fluidState.setSaturation(phaseIdx, value);
            }
            sumSat += value;

            value = fluidState.pressure(phaseIdx);
            if (value < 0)
                fluidState.setPressure(phaseIdx, 0.0);

            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                value = fluidState.moleFraction(phaseIdx, compIdx);
                if (value < 0)
                    fluidState.setMoleFraction(phaseIdx, compIdx, 0.0);
            }
        }

        // last saturation
        if (sumSat > 1.05) {
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                Scalar value = fluidState.saturation(phaseIdx)/(0.95*sumSat);
                fluidState.setSaturation(phaseIdx, value);
            }
        }
        */

        completeFluidState_<MaterialLaw>(fluidState, paramCache, matParams);

        return relError;
    }

    template <class MaterialLaw, class FluidState>
    static void completeFluidState_(FluidState &fluidState,
                                    ParameterCache &paramCache,
                                    const typename MaterialLaw::Params &matParams)
    {
        // calculate the saturation of the last phase as a function of
        // the other saturations
        Scalar sumSat = 0.0;
        for (int phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx)
            sumSat += fluidState.saturation(phaseIdx);
        fluidState.setSaturation(/*phaseIdx=*/numPhases - 1, 1.0 - sumSat);

        // update the pressures using the material law (saturations
        // and first pressure are already set because it is implicitly
        // solved for.)
        ComponentVector pC;
        MaterialLaw::capillaryPressures(pC, matParams, fluidState);
        for (int phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx)
            fluidState.setPressure(phaseIdx,
                                   fluidState.pressure(0)
                                   + (pC[phaseIdx] - pC[0]));

        // update the parameter cache
        paramCache.updateAll(fluidState, /*except=*/ParameterCache::Temperature);

        // update all densities and fugacity coefficients
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            fluidState.setDensity(phaseIdx, rho);

            for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
                Scalar phi = FluidSystem::fugacityCoefficient( fluidState, paramCache, phaseIdx, compIdx);
                fluidState.setFugacityCoefficient(phaseIdx, compIdx, phi);
            }
        }
    }

    static bool isPressureIdx_(int pvIdx)
    { return pvIdx == 0; }

    static bool isSaturationIdx_(int pvIdx)
    { return 1 <= pvIdx && pvIdx < numPhases; }

    static bool isMoleFracIdx_(int pvIdx)
    { return numPhases <= pvIdx && pvIdx < numPhases + numPhases*numComponents; }

    // retrieves a quantity from the fluid state
    template <class FluidState>
    static Scalar getQuantity_(const FluidState &fluidState, int pvIdx)
    {
        assert(pvIdx < numEq);

        // first pressure
        if (pvIdx < 1) {
            int phaseIdx = 0;
            return fluidState.pressure(phaseIdx);
        }
        // first M - 1 saturations
        else if (pvIdx < numPhases) {
            int phaseIdx = pvIdx - 1;
            return fluidState.saturation(phaseIdx);
        }
        // mole fractions
        else // if (pvIdx < numPhases + numPhases*numComponents)
        {
            int phaseIdx = (pvIdx - numPhases)/numComponents;
            int compIdx = (pvIdx - numPhases)%numComponents;
            return fluidState.moleFraction(phaseIdx, compIdx);
        }
    }

    // set a quantity in the fluid state
    template <class MaterialLaw, class FluidState>
    static void setQuantity_(FluidState &fluidState,
                             ParameterCache &paramCache,
                             const typename MaterialLaw::Params &matParams,
                             int pvIdx,
                             Scalar value)
    {
        assert(0 <= pvIdx && pvIdx < numEq);

        if (pvIdx < 1) { // <- first pressure
            Scalar delta = value - fluidState.pressure(0);

            // set all pressures. here we assume that the capillary
            // pressure does not depend on absolute pressure.
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                fluidState.setPressure(phaseIdx, fluidState.pressure(phaseIdx) + delta);
            paramCache.updateAllPressures(fluidState);

            // update all densities and fugacity coefficients
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
                fluidState.setDensity(phaseIdx, rho);

                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    Scalar phi = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx);
                    fluidState.setFugacityCoefficient(phaseIdx, compIdx, phi);
                }
            }
        }
        else if (pvIdx < numPhases) { // <- first M - 1 saturations
            Scalar delta = value - fluidState.saturation(/*phaseIdx=*/pvIdx - 1);
            fluidState.setSaturation(/*phaseIdx=*/pvIdx - 1, value);

            // set last saturation (-> minus the change of the
            // satuation of the other phase)
            fluidState.setSaturation(/*phaseIdx=*/numPhases - 1,
                                     fluidState.saturation(numPhases - 1) - delta);

            // update all fluid pressures using the capillary pressure
            // law
            ComponentVector pC;
            MaterialLaw::capillaryPressures(pC, matParams, fluidState);
            for (int phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx)
                fluidState.setPressure(phaseIdx,
                               fluidState.pressure(0)
                               + (pC[phaseIdx] - pC[0]));
            paramCache.updateAllPressures(fluidState);

            // update all densities and fugacity coefficients
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
                fluidState.setDensity(phaseIdx, rho);

                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    Scalar phi = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx);
                    fluidState.setFugacityCoefficient(phaseIdx, compIdx, phi);
                }
            }
        }
        else if (pvIdx < numPhases + numPhases*numComponents) // <- mole fractions
        {
            int phaseIdx = (pvIdx - numPhases)/numComponents;
            int compIdx = (pvIdx - numPhases)%numComponents;

            fluidState.setMoleFraction(phaseIdx, compIdx, value);
            paramCache.updateSingleMoleFraction(fluidState, phaseIdx, compIdx);

            // update the density of the phase
            Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            fluidState.setDensity(phaseIdx, rho);

            // if the phase's fugacity coefficients are composition
            // dependent, update them as well.
            if (!FluidSystem::isIdealMixture(phaseIdx)) {
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    Scalar phi = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx);
                    fluidState.setFugacityCoefficient(phaseIdx, compIdx, phi);
                }
            }
        }
        else {
            assert(false);
        }
    }

    // set a quantity in the fluid state
    template <class FluidState>
    static void setQuantityRaw_(FluidState &fluidState, int pvIdx, Scalar value)
    {
        assert(pvIdx < numEq);

        // first pressure
        if (pvIdx < 1) {
            int phaseIdx = 0;
            fluidState.setPressure(phaseIdx, value);
        }
        // first M - 1 saturations
        else if (pvIdx < numPhases) {
            int phaseIdx = pvIdx - 1;
            fluidState.setSaturation(phaseIdx, value);
        }
        // mole fractions
        else // if (pvIdx < numPhases + numPhases*numComponents)
        {
            int phaseIdx = (pvIdx - numPhases)/numComponents;
            int compIdx = (pvIdx - numPhases)%numComponents;
            fluidState.setMoleFraction(phaseIdx, compIdx, value);
        }
    }

    template <class FluidState>
    static Scalar quantityWeight_(const FluidState &fluidState, int pvIdx)
    {
        // first pressure
        if (pvIdx < 1)
            return 1/1e5;
        // first M - 1 saturations
        else if (pvIdx < numPhases)
            return 1.0;
        // mole fractions
        else // if (pvIdx < numPhases + numPhases*numComponents)
            return 1.0;
    }
};

} // end namespace Opm

#endif
