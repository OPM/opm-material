// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::NcpFlash
 */
#ifndef OPM_NCP_FLASH_HPP
#define OPM_NCP_FLASH_HPP

#include <opm/material/fluidmatrixinteractions/NullMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/common/MathToolbox.hpp>
#include <opm/common/Valgrind.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

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
 * We use the following catches: Capillary pressures are taken into
 * account explicitly, so that only the pressure of the first phase is
 * solved implicitly, also the closure condition for the saturations
 * is taken into account explicitly, which means that we don't need to
 * implicitly solve for the last saturation. These two measures reduce
 * the number of unknowns to M*(N + 1), namely:
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

    enum {
        p0PvIdx = 0,
        S0PvIdx = 1,
        x00PvIdx = S0PvIdx + numPhases - 1
    };

    static const int numEq = numPhases*(numComponents + 1);

public:
    /*!
     * \brief Guess initial values for all quantities.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static void guessInitial(FluidState& fluidState,
                             const Dune::FieldVector<Evaluation, numComponents>& globalMolarities)
    {
        // the sum of all molarities
        Evaluation sumMoles = 0;
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
            sumMoles += globalMolarities[compIdx];

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            // composition
            for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx)
                fluidState.setMoleFraction(phaseIdx,
                                           compIdx,
                                           globalMolarities[compIdx]/sumMoles);

            // pressure. use atmospheric pressure as initial guess
            fluidState.setPressure(phaseIdx, 1.0135e5);

            // saturation. assume all fluids to be equally distributed
            fluidState.setSaturation(phaseIdx, 1.0/numPhases);
        }

        // set the fugacity coefficients of all components in all phases
        typename FluidSystem::template ParameterCache<Evaluation> paramCache;
        paramCache.updateAll(fluidState);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx) {
                const typename FluidState::Scalar phi =
                    FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx);
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
    static void solve(FluidState& fluidState,
                      const typename MaterialLaw::Params& matParams,
                      typename FluidSystem::template ParameterCache<typename FluidState::Scalar>& paramCache,
                      const Dune::FieldVector<typename FluidState::Scalar, numComponents>& globalMolarities,
                      Scalar tolerance = -1.0)
    {
        typedef typename FluidState::Scalar InputEval;

        typedef Dune::FieldMatrix<InputEval, numEq, numEq> Matrix;
        typedef Dune::FieldVector<InputEval, numEq> Vector;

        typedef Opm::DenseAd::Evaluation</*Scalar=*/InputEval,
                                         /*numDerivs=*/numEq> FlashEval;

        typedef Dune::FieldVector<FlashEval, numEq> FlashDefectVector;
        typedef Opm::CompositionalFluidState<FlashEval, FluidSystem, /*energy=*/false> FlashFluidState;

        Dune::FMatrixPrecision<InputEval>::set_singular_limit(1e-35);

        if (tolerance <= 0)
            tolerance = std::min<Scalar>(1e-3,
                                         1e8*std::numeric_limits<Scalar>::epsilon());

        typename FluidSystem::template ParameterCache<FlashEval> flashParamCache;
        flashParamCache.assignPersistentData(paramCache);

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

        FlashFluidState flashFluidState;
        assignFlashFluidState_<MaterialLaw>(fluidState, flashFluidState, matParams, flashParamCache);

        // copy the global molarities to a vector of evaluations. Remember that the
        // global molarities are constants. (but we need to copy them to a vector of
        // FlashEvals anyway in order to avoid getting into hell's kitchen.)
        Dune::FieldVector<FlashEval, numComponents> flashGlobalMolarities;
        for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx)
            flashGlobalMolarities[compIdx] = globalMolarities[compIdx];

        FlashDefectVector defect;
        const unsigned nMax = 50; // <- maximum number of newton iterations
        for (unsigned nIdx = 0; nIdx < nMax; ++nIdx) {
            // calculate the defect of the flash equations and their derivatives
            evalDefect_(defect, flashFluidState, flashGlobalMolarities);
            Valgrind::CheckDefined(defect);

            // create field matrices and vectors out of the evaluation vector to solve
            // the linear system of equations.
            for (unsigned eqIdx = 0; eqIdx < numEq; ++ eqIdx) {
                for (unsigned pvIdx = 0; pvIdx < numEq; ++ pvIdx)
                    J[eqIdx][pvIdx] = defect[eqIdx].derivative(pvIdx);

                b[eqIdx] = defect[eqIdx].value();
            }
            Valgrind::CheckDefined(J);
            Valgrind::CheckDefined(b);

            // Solve J*x = b
            deltaX = 0.0;
            try { J.solve(deltaX, b); }
            catch (const Dune::FMatrixError& e) {
                throw Opm::NumericalProblem(e.what());
            }
            Valgrind::CheckDefined(deltaX);

            // update the fluid quantities.
            Scalar relError = update_<MaterialLaw>(flashFluidState, matParams, flashParamCache, deltaX);

            if (relError < tolerance) {
                assignOutputFluidState_(flashFluidState, fluidState);
                return;
            }
        }

        OPM_THROW(NumericalProblem,
                  "NcpFlash solver failed: "
                  "{c_alpha^kappa} = {" << globalMolarities << "}, "
                  << "T = " << fluidState.temperature(/*phaseIdx=*/0));
    }

    /*!
     * \brief Calculates the chemical equilibrium from the component
     *        fugacities in a phase.
     *
     * This is a convenience method which assumes that the capillary pressure is
     * zero...
     */
    template <class FluidState, class ComponentVector>
    static void solve(FluidState& fluidState,
                      const ComponentVector& globalMolarities,
                      Scalar tolerance = 0.0)
    {
        typedef NullMaterialTraits<Scalar, numPhases> MaterialTraits;
        typedef NullMaterial<MaterialTraits> MaterialLaw;
        typedef typename MaterialLaw::Params MaterialLawParams;

        MaterialLawParams matParams;
        solve<MaterialLaw>(fluidState, matParams, globalMolarities, tolerance);
    }


protected:
    template <class FluidState>
    static void printFluidState_(const FluidState& fluidState)
    {
        typedef typename FluidState::Scalar FsScalar;

        std::cout << "saturations: ";
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            std::cout << fluidState.saturation(phaseIdx) << " ";
        std::cout << "\n";

        std::cout << "pressures: ";
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            std::cout << fluidState.pressure(phaseIdx) << " ";
        std::cout << "\n";

        std::cout << "densities: ";
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            std::cout << fluidState.density(phaseIdx) << " ";
        std::cout << "\n";

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            std::cout << "composition " << FluidSystem::phaseName(phaseIdx) << "Phase: ";
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                std::cout << fluidState.moleFraction(phaseIdx, compIdx) << " ";
            }
            std::cout << "\n";
        }

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            std::cout << "fugacities " << FluidSystem::phaseName(phaseIdx) << "Phase: ";
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                std::cout << fluidState.fugacity(phaseIdx, compIdx) << " ";
            }
            std::cout << "\n";
        }

        std::cout << "global component molarities: ";
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            FsScalar sum = 0;
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                sum += fluidState.saturation(phaseIdx)*fluidState.molarity(phaseIdx, compIdx);
            }
            std::cout << sum << " ";
        }
        std::cout << "\n";
    }

    template <class MaterialLaw, class InputFluidState, class FlashFluidState>
    static void assignFlashFluidState_(const InputFluidState& inputFluidState,
                                       FlashFluidState& flashFluidState,
                                       const typename MaterialLaw::Params& matParams,
                                       typename FluidSystem::template ParameterCache<typename FlashFluidState::Scalar>& flashParamCache)
    {
        typedef typename FlashFluidState::Scalar FlashEval;

        // copy the temperature: even though the model which uses the flash solver might
        // be non-isothermal, the flash solver does not consider energy. (it could be
        // modified to do so relatively easily, but it would come at increased
        // computational cost and normally temperature instead of "total internal energy
        // of the fluids" is specified.)
        flashFluidState.setTemperature(inputFluidState.temperature(/*phaseIdx=*/0));

        // copy the saturations: the first N-1 phases are primary variables, the last one
        // is one minus the sum of the former.
        FlashEval Slast = 1.0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx) {
            FlashEval S = inputFluidState.saturation(phaseIdx);
            S.setDerivative(S0PvIdx + phaseIdx, 1.0);

            Slast -= S;

            flashFluidState.setSaturation(phaseIdx, S);
        }
        flashFluidState.setSaturation(numPhases - 1, Slast);

        // copy the pressures: the first pressure is the first primary variable, the
        // remaining ones are given as p_beta = p_alpha + p_calpha,beta
        FlashEval p0 = inputFluidState.pressure(0);
        p0.setDerivative(p0PvIdx, 1.0);

        std::array<FlashEval, numPhases> pc;
        MaterialLaw::capillaryPressures(pc, matParams, flashFluidState);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            flashFluidState.setPressure(phaseIdx, p0 + (pc[phaseIdx] - pc[0]));

        // copy the mole fractions: all of them are primary variables
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                FlashEval x = inputFluidState.moleFraction(phaseIdx, compIdx);
                x.setDerivative(x00PvIdx + phaseIdx*numComponents + compIdx, 1.0);
                flashFluidState.setMoleFraction(phaseIdx, compIdx, x);
            }
        }

        flashParamCache.updateAll(flashFluidState);

        // compute the density of each phase and the fugacity coefficient of each
        // component in each phase.
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            const FlashEval& rho = FluidSystem::density(flashFluidState, flashParamCache, phaseIdx);
            flashFluidState.setDensity(phaseIdx, rho);

            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                const FlashEval& fugCoeff = FluidSystem::fugacityCoefficient(flashFluidState, flashParamCache, phaseIdx, compIdx);
                flashFluidState.setFugacityCoefficient(phaseIdx, compIdx, fugCoeff);
            }
        }
    }

    template <class FlashFluidState, class OutputFluidState>
    static void assignOutputFluidState_(const FlashFluidState& flashFluidState,
                                        OutputFluidState& outputFluidState)
    {
        outputFluidState.setTemperature(flashFluidState.temperature(/*phaseIdx=*/0).value());

        // copy the saturations, pressures and densities
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            const auto& S = flashFluidState.saturation(phaseIdx).value();
            outputFluidState.setSaturation(phaseIdx, S);

            const auto& p = flashFluidState.pressure(phaseIdx).value();
            outputFluidState.setPressure(phaseIdx, p);

            const auto& rho = flashFluidState.density(phaseIdx).value();
            outputFluidState.setDensity(phaseIdx, rho);
        }

        // copy the mole fractions and fugacity coefficients
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                const auto& moleFrac =
                    flashFluidState.moleFraction(phaseIdx, compIdx).value();
                outputFluidState.setMoleFraction(phaseIdx, compIdx, moleFrac);

                const auto& fugCoeff =
                    flashFluidState.fugacityCoefficient(phaseIdx, compIdx).value();
                outputFluidState.setFugacityCoefficient(phaseIdx, compIdx, fugCoeff);
            }
        }
    }

    template <class FlashFluidState, class FlashDefectVector, class FlashComponentVector>
    static void evalDefect_(FlashDefectVector& b,
                            const FlashFluidState& fluidState,
                            const FlashComponentVector& globalMolarities)
    {
        typedef typename FlashFluidState::Scalar FlashEval;

        unsigned eqIdx = 0;

        // fugacity of any component must be equal in all phases
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            for (unsigned phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx) {
                b[eqIdx] =
                    fluidState.fugacity(/*phaseIdx=*/0, compIdx)
                    - fluidState.fugacity(phaseIdx, compIdx);
                ++eqIdx;
            }
        }
        assert(eqIdx == numComponents*(numPhases - 1));

        // the fact saturations must sum up to 1 is included implicitly and also,
        // capillary pressures are treated implicitly!

        // global molarities are given
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            b[eqIdx] = 0.0;
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                b[eqIdx] +=
                    fluidState.saturation(phaseIdx)
                    * fluidState.molarity(phaseIdx, compIdx);
            }

            b[eqIdx] -= globalMolarities[compIdx];
            ++eqIdx;
        }

        // model assumptions (-> non-linear complementarity functions) must be adhered
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            FlashEval oneMinusSumMoleFrac = 1.0;
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
                oneMinusSumMoleFrac -= fluidState.moleFraction(phaseIdx, compIdx);

            if (oneMinusSumMoleFrac > fluidState.saturation(phaseIdx))
                b[eqIdx] = fluidState.saturation(phaseIdx);
            else
                b[eqIdx] = oneMinusSumMoleFrac;

            ++eqIdx;
        }
    }

    template <class MaterialLaw, class FlashFluidState, class EvalVector>
    static Scalar update_(FlashFluidState& fluidState,
                          const typename MaterialLaw::Params& matParams,
                          typename FluidSystem::template ParameterCache<typename FlashFluidState::Scalar>& paramCache,
                          const EvalVector& deltaX)
    {
        // note that it is possible that FlashEval::Scalar is an Evaluation itself
        typedef typename FlashFluidState::Scalar FlashEval;
        typedef typename FlashEval::ValueType InnerEval;

#ifndef NDEBUG
        // make sure we don't swallow non-finite update vectors
        assert(deltaX.dimension == numEq);
        for (unsigned i = 0; i < numEq; ++i)
            assert(std::isfinite(Opm::scalarValue(deltaX[i])));
#endif

        Scalar relError = 0;
        for (unsigned pvIdx = 0; pvIdx < numEq; ++ pvIdx) {
            FlashEval tmp = getQuantity_(fluidState, pvIdx);
            InnerEval delta = deltaX[pvIdx];

            relError = std::max(relError,
                                std::abs(Opm::scalarValue(delta))
                                * quantityWeight_(fluidState, pvIdx));

            if (isSaturationIdx_(pvIdx)) {
                // dampen to at most 25% change in saturation per iteration
                delta = Opm::min(0.25, Opm::max(-0.25, delta));
            }
            else if (isMoleFracIdx_(pvIdx)) {
                // dampen to at most 20% change in mole fraction per iteration
                delta = Opm::min(0.20, Opm::max(-0.20, delta));
            }
            else if (isPressureIdx_(pvIdx)) {
                // dampen to at most 50% change in pressure per iteration
                delta = Opm::min(0.5*fluidState.pressure(0).value(),
                                 Opm::max(-0.5*fluidState.pressure(0).value(),
                                          delta));
            }

            tmp -= delta;
            setQuantity_(fluidState, pvIdx, tmp);
        }

        completeFluidState_<MaterialLaw>(fluidState, paramCache, matParams);

        return relError;
    }

    template <class MaterialLaw, class FlashFluidState>
    static void completeFluidState_(FlashFluidState& flashFluidState,
                                    typename FluidSystem::template ParameterCache<typename FlashFluidState::Scalar>& paramCache,
                                    const typename MaterialLaw::Params& matParams)
    {
        typedef typename FluidSystem::template ParameterCache<typename FlashFluidState::Scalar> ParamCache;

        typedef typename FlashFluidState::Scalar FlashEval;

        // calculate the saturation of the last phase as a function of
        // the other saturations
        FlashEval sumSat = 0.0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx)
            sumSat += flashFluidState.saturation(phaseIdx);
        flashFluidState.setSaturation(/*phaseIdx=*/numPhases - 1, 1.0 - sumSat);

        // update the pressures using the material law (saturations
        // and first pressure are already set because it is implicitly
        // solved for.)
        Dune::FieldVector<FlashEval, numPhases> pC;
        MaterialLaw::capillaryPressures(pC, matParams, flashFluidState);
        for (unsigned phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx)
            flashFluidState.setPressure(phaseIdx,
                                        flashFluidState.pressure(0)
                                        + (pC[phaseIdx] - pC[0]));

        // update the parameter cache
        paramCache.updateAll(flashFluidState, /*except=*/ParamCache::Temperature);

        // update all densities and fugacity coefficients
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            const FlashEval& rho = FluidSystem::density(flashFluidState, paramCache, phaseIdx);
            flashFluidState.setDensity(phaseIdx, rho);

            for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx) {
                const FlashEval& phi = FluidSystem::fugacityCoefficient(flashFluidState, paramCache, phaseIdx, compIdx);
                flashFluidState.setFugacityCoefficient(phaseIdx, compIdx, phi);
            }
        }
    }

    static bool isPressureIdx_(unsigned pvIdx)
    { return pvIdx == 0; }

    static bool isSaturationIdx_(unsigned pvIdx)
    { return 1 <= pvIdx && pvIdx < numPhases; }

    static bool isMoleFracIdx_(unsigned pvIdx)
    { return numPhases <= pvIdx && pvIdx < numPhases + numPhases*numComponents; }

    // retrieves a quantity from the fluid state
    template <class FluidState>
    static const typename FluidState::Scalar& getQuantity_(const FluidState& fluidState, unsigned pvIdx)
    {
        assert(pvIdx < numEq);

        // first pressure
        if (pvIdx < 1) {
            unsigned phaseIdx = 0;
            return fluidState.pressure(phaseIdx);
        }
        // first M - 1 saturations
        else if (pvIdx < numPhases) {
            unsigned phaseIdx = pvIdx - 1;
            return fluidState.saturation(phaseIdx);
        }
        // mole fractions
        else // if (pvIdx < numPhases + numPhases*numComponents)
        {
            unsigned phaseIdx = (pvIdx - numPhases)/numComponents;
            unsigned compIdx = (pvIdx - numPhases)%numComponents;
            return fluidState.moleFraction(phaseIdx, compIdx);
        }
    }

    // set a quantity in the fluid state
    template <class FluidState>
    static void setQuantity_(FluidState& fluidState,
                             unsigned pvIdx,
                             const typename FluidState::Scalar& value)
    {
        assert(pvIdx < numEq);

        Valgrind::CheckDefined(value);
        // first pressure
        if (pvIdx < 1) {
            unsigned phaseIdx = 0;
            fluidState.setPressure(phaseIdx, value);
        }
        // first M - 1 saturations
        else if (pvIdx < numPhases) {
            unsigned phaseIdx = pvIdx - 1;
            fluidState.setSaturation(phaseIdx, value);
        }
        // mole fractions
        else {
            assert(pvIdx < numPhases + numPhases*numComponents);
            unsigned phaseIdx = (pvIdx - numPhases)/numComponents;
            unsigned compIdx = (pvIdx - numPhases)%numComponents;
            fluidState.setMoleFraction(phaseIdx, compIdx, value);
        }
    }

    template <class FluidState>
    static Scalar quantityWeight_(const FluidState& /*fluidState*/, unsigned pvIdx)
    {
        // first pressure
        if (pvIdx < 1)
            return 1e-6;
        // first M - 1 saturations
        else if (pvIdx < numPhases)
            return 1.0;
        // mole fractions
        else // if (pvIdx < numPhases + numPhases*numComponents)
            return 1.0;
    }
};

} // namespace Opm

#endif
