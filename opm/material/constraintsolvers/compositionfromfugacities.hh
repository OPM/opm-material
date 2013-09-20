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
 * \copydoc Opm::CompositionFromFugacities
 */
#ifndef OPM_COMPOSITION_FROM_FUGACITIES_HH
#define OPM_COMPOSITION_FROM_FUGACITIES_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/Exceptions.hpp>
#include <opm/material/valgrind.hh>

#include <limits>

namespace Opm {

/*!
 * \brief Calculates the chemical equilibrium from the component
 *        fugacities in a phase.
 */
template <class Scalar, class FluidSystem>
class CompositionFromFugacities
{
    enum { numComponents = FluidSystem::numComponents };

    typedef typename FluidSystem::ParameterCache ParameterCache;

public:
    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;

    /*!
     * \brief Guess an initial value for the composition of the phase.
     */
    template <class FluidState>
    static void guessInitial(FluidState &fluidState,
                             ParameterCache &paramCache,
                             int phaseIdx,
                             const ComponentVector &fugVec)
    {
        if (FluidSystem::isIdealMixture(phaseIdx))
            return;

        // Pure component fugacities
        for (int i = 0; i < numComponents; ++ i) {
            //std::cout << f << " -> " << mutParams.fugacity(phaseIdx, i)/f << "\n";
            fluidState.setMoleFraction(phaseIdx,
                                   i,
                                   1.0/numComponents);
        }
    }

    /*!
     * \brief Calculates the chemical equilibrium from the component
     *        fugacities in a phase.
     *
     * The phase's fugacities must already be set.
     */
    template <class FluidState>
    static void solve(FluidState &fluidState,
                      ParameterCache &paramCache,
                      int phaseIdx,
                      const ComponentVector &targetFug)
    {
        // use a much more efficient method in case the phase is an
        // ideal mixture
        if (FluidSystem::isIdealMixture(phaseIdx)) {
            solveIdealMix_(fluidState, paramCache, phaseIdx, targetFug);
            return;
        }

        Dune::FMatrixPrecision<Scalar>::set_singular_limit(1e-25);

        // save initial composition in case something goes wrong
        Dune::FieldVector<Scalar, numComponents> xInit;
        for (int i = 0; i < numComponents; ++i) {
            xInit[i] = fluidState.moleFraction(phaseIdx, i);
        }

        /////////////////////////
        // Newton method
        /////////////////////////

        // Jacobian matrix
        Dune::FieldMatrix<Scalar, numComponents, numComponents> J;
        // solution, i.e. phase composition
        Dune::FieldVector<Scalar, numComponents> x;
        // right hand side
        Dune::FieldVector<Scalar, numComponents> b;

        paramCache.updatePhase(fluidState, phaseIdx);

        // maximum number of iterations
        const int nMax = 25;
        for (int nIdx = 0; nIdx < nMax; ++nIdx) {
            // calculate Jacobian matrix and right hand side
            linearize_(J, b, fluidState, paramCache, phaseIdx, targetFug);
            Valgrind::CheckDefined(J);
            Valgrind::CheckDefined(b);

            /*
            std::cout << FluidSystem::phaseName(phaseIdx) << "Phase composition: ";
            for (int i = 0; i < FluidSystem::numComponents; ++i)
                std::cout << fluidState.moleFraction(phaseIdx, i) << " ";
            std::cout << "\n";
            std::cout << FluidSystem::phaseName(phaseIdx) << "Phase phi: ";
            for (int i = 0; i < FluidSystem::numComponents; ++i)
                std::cout << fluidState.fugacityCoefficient(phaseIdx, i) << " ";
            std::cout << "\n";
            */

            // Solve J*x = b
            x = 0;
            try { J.solve(x, b); }
            catch (Dune::FMatrixError e)
            { throw Opm::NumericalProblem(e.what()); }

            //std::cout << "original delta: " << x << "\n";

            Valgrind::CheckDefined(x);

            /*
            std::cout << FluidSystem::phaseName(phaseIdx) << "Phase composition: ";
            for (int i = 0; i < FluidSystem::numComponents; ++i)
                std::cout << fluidState.moleFraction(phaseIdx, i) << " ";
            std::cout << "\n";
            std::cout << "J: " << J << "\n";
            std::cout << "rho: " << fluidState.density(phaseIdx) << "\n";
            std::cout << "delta: " << x << "\n";
            std::cout << "defect: " << b << "\n";

            std::cout << "J: " << J << "\n";

            std::cout << "---------------------------\n";
            */

            // update the fluid composition. b is also used to store
            // the defect for the next iteration.
            Scalar relError = update_(fluidState, paramCache, x, b, phaseIdx, targetFug);

            if (relError < 1e-9) {
                Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
                fluidState.setDensity(phaseIdx, rho);

                //std::cout << "num iterations: " << nIdx << "\n";
                return;
            }
        }

        OPM_THROW(Opm::NumericalProblem,
                  "Calculating the " << FluidSystem::phaseName(phaseIdx)
                  << "Phase composition failed. Initial {x} = {"
                  << xInit
                  << "}, {fug_t} = {" << targetFug << "}, p = " << fluidState.pressure(phaseIdx)
                  << ", T = " << fluidState.temperature(phaseIdx));
    }


protected:
    // update the phase composition in case the phase is an ideal
    // mixture, i.e. the component's fugacity coefficients are
    // independent of the phase's composition.
    template <class FluidState>
    static void solveIdealMix_(FluidState &fluidState,
                               ParameterCache &paramCache,
                               int phaseIdx,
                               const ComponentVector &fugacities)
    {
        for (int i = 0; i < numComponents; ++ i) {
            Scalar phi = FluidSystem::fugacityCoefficient(fluidState,
                                                          paramCache,
                                                          phaseIdx,
                                                          i);
            Scalar gamma = phi * fluidState.pressure(phaseIdx);
            Valgrind::CheckDefined(phi);
            Valgrind::CheckDefined(gamma);
            Valgrind::CheckDefined(fugacities[i]);
            fluidState.setFugacityCoefficient(phaseIdx, i, phi);
            fluidState.setMoleFraction(phaseIdx, i, fugacities[i]/gamma);
        };

        paramCache.updatePhase(fluidState, phaseIdx);

        Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
        fluidState.setDensity(phaseIdx, rho);
        return;
    }

    template <class FluidState>
    static Scalar linearize_(Dune::FieldMatrix<Scalar, numComponents, numComponents> &J,
                             Dune::FieldVector<Scalar, numComponents> &defect,
                             FluidState &fluidState,
                             ParameterCache &paramCache,
                             int phaseIdx,
                             const ComponentVector &targetFug)
    {
        // reset jacobian
        J = 0;

        Scalar absError = 0;
        // calculate the defect (deviation of the current fugacities
        // from the target fugacities)
        for (int i = 0; i < numComponents; ++ i) {
            Scalar phi = FluidSystem::fugacityCoefficient(fluidState,
                                                          paramCache,
                                                          phaseIdx,
                                                          i);
            Scalar f = phi*fluidState.pressure(phaseIdx)*fluidState.moleFraction(phaseIdx, i);
            fluidState.setFugacityCoefficient(phaseIdx, i, phi);

            defect[i] = targetFug[i] - f;
            absError = std::max(absError, std::abs(defect[i]));
        }

        // assemble jacobian matrix of the constraints for the composition
        static const Scalar eps = std::numeric_limits<Scalar>::epsilon()*1e6;
        for (int i = 0; i < numComponents; ++ i) {
            ////////
            // approximately calculate partial derivatives of the
            // fugacity defect of all components in regard to the mole
            // fraction of the i-th component. This is done via
            // forward differences

            // deviate the mole fraction of the i-th component
            Scalar xI = fluidState.moleFraction(phaseIdx, i);
            fluidState.setMoleFraction(phaseIdx, i, xI + eps);
            paramCache.updateSingleMoleFraction(fluidState, phaseIdx, i);

            // compute new defect and derivative for all component
            // fugacities
            for (int j = 0; j < numComponents; ++j) {
                // compute the j-th component's fugacity coefficient ...
                Scalar phi = FluidSystem::fugacityCoefficient(fluidState,
                                                              paramCache,
                                                              phaseIdx,
                                                              j);
                // ... and its fugacity ...
                Scalar f =
                    phi *
                    fluidState.pressure(phaseIdx) *
                    fluidState.moleFraction(phaseIdx, j);
                // as well as the defect for this fugacity
                Scalar defJPlusEps = targetFug[j] - f;

                // use forward differences to calculate the defect's
                // derivative
                J[j][i] = (defJPlusEps - defect[j])/eps;
            }

            // reset composition to original value
            fluidState.setMoleFraction(phaseIdx, i, xI);
            paramCache.updateSingleMoleFraction(fluidState, phaseIdx, i);

            // end forward differences
            ////////
        }

        return absError;
    }

    template <class FluidState>
    static Scalar update_(FluidState &fluidState,
                          ParameterCache &paramCache,
                          Dune::FieldVector<Scalar, numComponents> &x,
                          Dune::FieldVector<Scalar, numComponents> &b,
                          int phaseIdx,
                          const Dune::FieldVector<Scalar, numComponents> &targetFug)
    {
        // store original composition and calculate relative error
        Dune::FieldVector<Scalar, numComponents> origComp;
        Scalar relError = 0;
        Scalar sumDelta = 0;
        Scalar sumx = 0;
        for (int i = 0; i < numComponents; ++i) {
            origComp[i] = fluidState.moleFraction(phaseIdx, i);
            relError = std::max(relError, std::abs(x[i]));

            sumx += std::abs(fluidState.moleFraction(phaseIdx, i));
            sumDelta += std::abs(x[i]);
        }

#if 1
        // chop update to at most 20% change in composition
        const Scalar maxDelta = 0.2;
        if (sumDelta > maxDelta)
            x /= (sumDelta/maxDelta);
#endif

        //Scalar curDefect = calculateDefect_(fluidState, phaseIdx, targetFug);
        //Scalar nextDefect;
        //Scalar sumMoleFrac = 0.0;
        //ComponentVector newB(1e100);
        //for (int numTries = 0; numTries < 1; ++numTries) {
            // change composition
            for (int i = 0; i < numComponents; ++i) {
                Scalar newx = origComp[i] - x[i];
#if 1
                // only allow negative mole fractions if the target fugacity is negative
                if (targetFug[i] > 0)
                    newx = std::max(0.0, newx);
                // only allow positive mole fractions if the target fugacity is positive
                else if (targetFug[i] < 0)
                    newx = std::min(0.0, newx);
                // if the target fugacity is zero, the mole fraction must also be zero
                else
                    newx = 0;
#endif
                fluidState.setMoleFraction(phaseIdx, i, newx);
                //sumMoleFrac += std::abs(newx);
            }

            paramCache.updateComposition(fluidState, phaseIdx);

            /*
            // if the sum of the mole fractions gets 0, we take the
            // original composition divided by 100
            if (sumMoleFrac < 1e-10) {
                for (int i = 0; i < numComponents; ++i) {
                    fluidState.setMoleFraction(phaseIdx, i, origComp[i]/100);
                }
                return relError;
            }
            */

            /*
            // calculate new residual
            for (int i = 0; i < numComponents; ++i) {
                Scalar phi = FluidSystem::computeFugacityCoeff(fluidState,
                                                               phaseIdx,
                                                               i);
                fluidState.setFugacityCoeff(phaseIdx, i, phi);
            }

            nextDefect = calculateDefect_(fluidState, phaseIdx, targetFug);
            //std::cout << "try delta: " << x << "\n";
            //std::cout << "defect: old=" << curDefect << " new=" << nextDefect << "\n";
            if (nextDefect <= curDefect)
                break;

             // divide delta vector
             x /= 2;
        }
            */

        return relError;
    }

    template <class FluidState>
    static Scalar calculateDefect_(const FluidState &params,
                                   int phaseIdx,
                                   const ComponentVector &targetFug)
    {
        Scalar result = 0.0;
        for (int i = 0; i < numComponents; ++i) {
            // sum of the fugacity defect weighted by the inverse
            // fugacity coefficient
            result += std::abs(
                (targetFug[i] - params.fugacity(phaseIdx, i))
                /
                params.fugacityCoefficient(phaseIdx, i) );
        };
        return result;
    }
};

} // end namespace Opm

#endif
