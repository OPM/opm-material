// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2015 by Andreas Lauser

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
*/
/*!
 * \file
 * \copydoc Opm::WetGasPvt
 */
#ifndef OPM_WET_GAS_PVT_HPP
#define OPM_WET_GAS_PVT_HPP

#include <opm/material/Constants.hpp>

#include <opm/material/common/OpmFinal.hpp>
#include <opm/material/common/UniformXTabulated2DFunction.hpp>
#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/common/Spline.hpp>

#if HAVE_OPM_PARSER
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#endif

namespace Opm {

/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the gas phas
 *        with vaporized oil.
 */
template <class Scalar>
class WetGasPvt
{
    typedef Opm::UniformXTabulated2DFunction<Scalar> TabulatedTwoDFunction;
    typedef Opm::Tabulated1DFunction<Scalar> TabulatedOneDFunction;
    typedef Opm::Spline<Scalar> Spline;
    typedef std::vector<std::pair<Scalar, Scalar> > SamplingPoints;

public:
#if HAVE_OPM_PARSER
    /*!
     * \brief Initialize the parameters for wet gas using an ECL deck.
     *
     * This method assumes that the deck features valid DENSITY and PVTG keywords.
     */
    void initFromDeck(DeckConstPtr deck, EclipseStateConstPtr eclState)
    {
        const auto& pvtgTables = eclState->getTableManager()->getPvtgTables();
        DeckKeywordConstPtr densityKeyword = deck->getKeyword("DENSITY");

        assert(pvtgTables.size() == densityKeyword->size());

        size_t numRegions = pvtgTables.size();
        setNumRegions(numRegions);

        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            Scalar rhoRefO = densityKeyword->getRecord(regionIdx)->getItem("OIL")->getSIDouble(0);
            Scalar rhoRefG = densityKeyword->getRecord(regionIdx)->getItem("GAS")->getSIDouble(0);
            Scalar rhoRefW = densityKeyword->getRecord(regionIdx)->getItem("WATER")->getSIDouble(0);

            setReferenceDensities(regionIdx, rhoRefO, rhoRefG, rhoRefW);
        }

        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            const auto& pvtgTable = pvtgTables[regionIdx];

            const auto saturatedTable = pvtgTable.getSaturatedTable();
            assert(saturatedTable.numRows() > 1);

            auto& gasMu = gasMu_[regionIdx];
            auto& invGasB = inverseGasB_[regionIdx];
            auto& invSatGasB = inverseSaturatedGasB_[regionIdx];
            auto& invSatGasBMu = inverseSaturatedGasBMu_[regionIdx];
            auto& oilVaporizationFac = saturatedOilVaporizationFactorTable_[regionIdx];

            {
                std::vector<double> pressure = saturatedTable.getColumn("PG").vectorCopy( );
                std::vector<double> rv = saturatedTable.getColumn("RV").vectorCopy( );
                oilVaporizationFac.setXYArrays(saturatedTable.numRows(),
                                               pressure , rv );
            }

            std::vector<Scalar> invSatGasBArray;
            std::vector<Scalar> invSatGasBMuArray;

            // extract the table for the gas dissolution and the oil formation volume factors
            for (unsigned outerIdx = 0; outerIdx < saturatedTable.numRows(); ++ outerIdx) {
                Scalar pg = saturatedTable.get("PG" , outerIdx);
                Scalar B = saturatedTable.get("BG" , outerIdx);
                Scalar mu = saturatedTable.get("MUG" , outerIdx);

                invGasB.appendXPos(pg);
                gasMu.appendXPos(pg);

                invSatGasBArray.push_back(1.0/B);
                invSatGasBMuArray.push_back(1.0/(mu*B));

                assert(invGasB.numX() == outerIdx + 1);
                assert(gasMu.numX() == outerIdx + 1);

                const auto underSaturatedTable = pvtgTable.getUnderSaturatedTable(outerIdx);
                size_t numRows = underSaturatedTable.numRows();
                for (size_t innerIdx = 0; innerIdx < numRows; ++ innerIdx) {
                    Scalar Rv = underSaturatedTable.get("RV" , innerIdx);
                    Scalar Bg = underSaturatedTable.get("BG" , innerIdx);
                    Scalar mug = underSaturatedTable.get("MUG" , innerIdx);

                    invGasB.appendSamplePoint(outerIdx, Rv, 1.0/Bg);
                    gasMu.appendSamplePoint(outerIdx, Rv, mug);
                }
            }

            {
                std::vector<double> tmpPressure =  saturatedTable.getColumn("PG").vectorCopy( );

                invSatGasB.setXYContainers(tmpPressure, invSatGasBArray);
                invSatGasBMu.setXYContainers(tmpPressure, invSatGasBMuArray);
            }

            // make sure to have at least two sample points per gas pressure value
            for (unsigned xIdx = 0; xIdx < invGasB.numX(); ++xIdx) {
               // a single sample point is definitely needed
                assert(invGasB.numY(xIdx) > 0);

                // everything is fine if the current table has two or more sampling points
                // for a given mole fraction
                if (invGasB.numY(xIdx) > 1)
                    continue;

                // find the master table which will be used as a template to extend the
                // current line. We define master table as the first table which has values
                // for undersaturated gas...
                size_t masterTableIdx = xIdx + 1;
                for (; masterTableIdx < saturatedTable.numRows(); ++masterTableIdx)
                {
                    if (pvtgTable.getUnderSaturatedTable(masterTableIdx).numRows() > 1)
                        break;
                }

                if (masterTableIdx >= saturatedTable.numRows())
                    OPM_THROW(std::runtime_error,
                              "PVTG tables are invalid: The last table must exhibit at least one "
                              "entry for undersaturated gas!");


                // extend the current table using the master table.
                extendPvtgTable_(regionIdx,
                                 xIdx,
                                 pvtgTable.getUnderSaturatedTable(xIdx),
                                 pvtgTable.getUnderSaturatedTable(masterTableIdx));
            }
        }
    }

private:
    void extendPvtgTable_(unsigned regionIdx,
                          unsigned xIdx,
                          const SimpleTable& curTable,
                          const SimpleTable& masterTable)
    {
        std::vector<double> RvArray = curTable.getColumn("RV").vectorCopy();
        std::vector<double> gasBArray = curTable.getColumn("BG").vectorCopy();
        std::vector<double> gasMuArray = curTable.getColumn("MUG").vectorCopy();

        auto& invGasB = inverseGasB_[regionIdx];
        auto& gasMu = gasMu_[regionIdx];

        for (size_t newRowIdx = 1; newRowIdx < masterTable.numRows(); ++ newRowIdx) {
            const auto& RVColumn = masterTable.getColumn("RV");
            const auto& BGColumn = masterTable.getColumn("BG");
            const auto& viscosityColumn = masterTable.getColumn("MUG");

            // compute the gas pressure for the new entry
            Scalar diffRv = RVColumn[newRowIdx] - RVColumn[newRowIdx - 1];
            Scalar newRv = RvArray.back() + diffRv;

            // calculate the compressibility of the master table
            Scalar B1 = BGColumn[newRowIdx];
            Scalar B2 = BGColumn[newRowIdx - 1];
            Scalar x = (B1 - B2)/( (B1 + B2)/2.0 );

            // calculate the gas formation volume factor which exhibits the same
            // "compressibility" for the new value of Rv
            Scalar newBg = gasBArray.back()*(1.0 + x/2.0)/(1.0 - x/2.0);

            // calculate the "viscosibility" of the master table
            Scalar mu1 = viscosityColumn[newRowIdx];
            Scalar mu2 = viscosityColumn[newRowIdx - 1];
            Scalar xMu = (mu1 - mu2)/( (mu1 + mu2)/2.0 );

            // calculate the gas formation volume factor which exhibits the same
            // compressibility for the new pressure
            Scalar newMug = gasMuArray.back()*(1.0 + xMu/2)/(1.0 - xMu/2.0);

            // append the new values to the arrays which we use to compute the additional
            // values ...
            RvArray.push_back(newRv);
            gasBArray.push_back(newBg);
            gasMuArray.push_back(newMug);

            // ... and register them with the internal table objects
            invGasB.appendSamplePoint(xIdx, newRv, 1.0/newBg);
            gasMu.appendSamplePoint(xIdx, newRv, newMug);
        }
    }

public:
#endif // HAVE_OPM_PARSER

    void setNumRegions(size_t numRegions)
    {
        oilReferenceDensity_.resize(numRegions);
        gasReferenceDensity_.resize(numRegions);
        inverseGasB_.resize(numRegions);
        inverseGasBMu_.resize(numRegions);
        inverseSaturatedGasB_.resize(numRegions);
        inverseSaturatedGasBMu_.resize(numRegions);
        gasMu_.resize(numRegions);
        saturatedOilVaporizationFactorTable_.resize(numRegions);
        saturationPressureSpline_.resize(numRegions);
    }

    /*!
     * \brief Initialize the reference densities of all fluids for a given PVT region
     */
    void setReferenceDensities(unsigned regionIdx,
                               Scalar rhoRefOil,
                               Scalar rhoRefGas,
                               Scalar /*rhoRefWater*/)
    {
        oilReferenceDensity_[regionIdx] = rhoRefOil;
        gasReferenceDensity_[regionIdx] = rhoRefGas;
    }

    /*!
     * \brief Initialize the function for the oil vaporization factor \f$R_v\f$
     *
     * \param samplePoints A container of (x,y) values.
     */
    void setSaturatedGasOilVaporizationFactor(unsigned regionIdx, const SamplingPoints &samplePoints)
    { saturatedOilVaporizationFactorTable_[regionIdx].setContainerOfTuples(samplePoints); }

    /*!
     * \brief Initialize the function for the gas formation volume factor
     *
     * The gas formation volume factor \f$B_g\f$ is a function of \f$(p_g, X_g^O)\f$ and
     * represents the partial density of the oil component in the gas phase at a given
     * pressure. This method only requires the volume factor of oil-saturated gas (which
     * only depends on pressure) while the dependence on the oil mass fraction is
     * guesstimated...
     */
    void setSaturatedGasFormationVolumeFactor(unsigned regionIdx, const SamplingPoints &samplePoints)
    {
        auto& invGasB = inverseGasB_[regionIdx];

        auto &RvTable = saturatedOilVaporizationFactorTable_[regionIdx];

        Scalar T = 273.15 + 15.56; // [K]

        Scalar RvMin = 0.0;
        Scalar RvMax = RvTable.eval(saturatedOilVaporizationFactorTable_[regionIdx].xMax(), /*extrapolate=*/true);

        Scalar poMin = samplePoints.front().first;
        Scalar poMax = samplePoints.back().first;

        size_t nRv = 20;
        size_t nP = samplePoints.size()*2;

        Scalar rhogRef = gasReferenceDensity_[regionIdx];
        Scalar rhooRef = oilReferenceDensity_[regionIdx];

        Spline gasFormationVolumeFactorSpline;
        gasFormationVolumeFactorSpline.setContainerOfTuples(samplePoints, /*type=*/Spline::Monotonic);

        updateSaturationPressureSpline_(regionIdx);

        // calculate a table of estimated densities depending on pressure and gas mass
        // fraction. note that this assumes oil of constant compressibility. (having said
        // that, if only the saturated gas densities are available, there's not much
        // choice.)
        for (size_t RvIdx = 0; RvIdx < nRv; ++RvIdx) {
            Scalar Rv = RvMin + (RvMax - RvMin)*RvIdx/nRv;

            invGasB.appendXPos(Rv);

            for (size_t pIdx = 0; pIdx < nP; ++pIdx) {
                Scalar pg = poMin + (poMax - poMin)*pIdx/nP;

                Scalar poSat = saturationPressure(regionIdx, T, Rv);
                Scalar BgSat = gasFormationVolumeFactorSpline.eval(poSat, /*extrapolate=*/true);
                Scalar drhoo_dp = (1.1200 - 1.1189)/((5000 - 4000)*6894.76);
                Scalar rhoo = rhooRef/BgSat*(1 + drhoo_dp*(pg - poSat));

                Scalar Bg = rhooRef/rhoo;

                invGasB.appendSamplePoint(RvIdx, pg, 1.0/Bg);
            }
        }
    }

    /*!
     * \brief Initialize the function for the gas formation volume factor
     *
     * The gas formation volume factor \f$B_g\f$ is a function of \f$(p_g, X_g^O)\f$ and
     * represents the partial density of the oil component in the gas phase at a given
     * pressure.
     *
     * This method sets \f$1/B_g(R_v, p_g)\f$. Note that instead of the mass fraction of
     * the oil component in the gas phase, this function depends on the gas dissolution
     * factor. Also note, that the order of the arguments needs to be \f$(R_s, p_o)\f$
     * and not the other way around.
     */
    void setInverseGasFormationVolumeFactor(unsigned regionIdx, const TabulatedTwoDFunction& invBg)
    { inverseGasB_[regionIdx] = invBg; }

    /*!
     * \brief Initialize the viscosity of the gas phase.
     *
     * This is a function of \f$(R_s, p_o)\f$...
     */
    void setGasViscosity(unsigned regionIdx, const TabulatedTwoDFunction& mug)
    { gasMu_[regionIdx] = mug; }

    /*!
     * \brief Initialize the phase viscosity for oil saturated gas
     *
     * The gas viscosity is a function of \f$(p_g, X_g^O)\f$, but this method only
     * requires the viscosity of oil-saturated gas (which only depends on pressure) while
     * there is assumed to be no dependence on the gas mass fraction...
     */
    void setSaturatedGasViscosity(unsigned regionIdx, const SamplingPoints &samplePoints  )
    {
        auto& oilVaporizationFac = saturatedOilVaporizationFactorTable_[regionIdx];

        Scalar RvMin = 0.0;
        Scalar RvMax = oilVaporizationFac.eval(saturatedOilVaporizationFactorTable_[regionIdx].xMax(), /*extrapolate=*/true);

        Scalar poMin = samplePoints.front().first;
        Scalar poMax = samplePoints.back().first;

        size_t nRv = 20;
        size_t nP = samplePoints.size()*2;

        Spline mugSpline;
        mugSpline.setContainerOfTuples(samplePoints, /*type=*/Spline::Monotonic);

        // calculate a table of estimated densities depending on pressure and gas mass
        // fraction
        for (size_t RvIdx = 0; RvIdx < nRv; ++RvIdx) {
            Scalar Rv = RvMin + (RvMax - RvMin)*RvIdx/nRv;

            gasMu_[regionIdx].appendXPos(Rv);

            for (size_t pIdx = 0; pIdx < nP; ++pIdx) {
                Scalar pg = poMin + (poMax - poMin)*pIdx/nP;
                Scalar mug = mugSpline.eval(pg, /*extrapolate=*/true);

                gasMu_[regionIdx].appendSamplePoint(RvIdx, pg, mug);
            }
        }
    }

    /*!
     * \brief Finish initializing the gas phase PVT properties.
     */
    void initEnd()
    {
        // calculate the final 2D functions which are used for interpolation.
        size_t numRegions = gasMu_.size();
        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            // calculate the table which stores the inverse of the product of the gas
            // formation volume factor and the gas viscosity
            const auto& gasMu = gasMu_[regionIdx];
            const auto& invGasB = inverseGasB_[regionIdx];
            assert(gasMu.numX() == invGasB.numX());

            auto& invGasBMu = inverseGasBMu_[regionIdx];
            auto& invSatGasB = inverseSaturatedGasB_[regionIdx];
            auto& invSatGasBMu = inverseSaturatedGasBMu_[regionIdx];

            std::vector<Scalar> satPressuresArray;
            std::vector<Scalar> invSatGasBArray;
            std::vector<Scalar> invSatGasBMuArray;
            for (size_t pIdx = 0; pIdx < gasMu.numX(); ++pIdx) {
                invGasBMu.appendXPos(gasMu.xAt(pIdx));

                assert(gasMu.numY(pIdx) == invGasB.numY(pIdx));

                size_t numRv = gasMu.numY(pIdx);
                for (size_t rvIdx = 0; rvIdx < numRv; ++rvIdx)
                    invGasBMu.appendSamplePoint(pIdx,
                                                gasMu.yAt(pIdx, rvIdx),
                                                invGasB.valueAt(pIdx, rvIdx)
                                                / gasMu.valueAt(pIdx, rvIdx));

                // the sampling points in UniformXTabulated2DFunction are always sorted
                // in ascending order. Thus, the value for saturated gas is the last one
                // (i.e., the one with the largest Rv value)
                satPressuresArray.push_back(gasMu.xAt(pIdx));
                invSatGasBArray.push_back(invGasB.valueAt(pIdx, numRv - 1));
                invSatGasBMuArray.push_back(invGasBMu.valueAt(pIdx, numRv - 1));
            }

            invSatGasB.setXYContainers(satPressuresArray, invSatGasBArray);
            invSatGasBMu.setXYContainers(satPressuresArray, invSatGasBMuArray);

            updateSaturationPressureSpline_(regionIdx);
        }
    }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& /*temperature*/,
                         const Evaluation& pressure,
                         const Evaluation& Rv) const
    {
        const Evaluation& invBg = inverseGasB_[regionIdx].eval(pressure, Rv, /*extrapolate=*/true);
        const Evaluation& invMugBg = inverseGasBMu_[regionIdx].eval(pressure, Rv, /*extrapolate=*/true);

        return invBg/invMugBg;
    }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of oil saturated gas at a given pressure.
     */
    template <class Evaluation>
    Evaluation saturatedViscosity(unsigned regionIdx,
                                  const Evaluation& /*temperature*/,
                                  const Evaluation& pressure) const
    {
        const Evaluation& invBg = inverseSaturatedGasB_[regionIdx].eval(pressure, /*extrapolate=*/true);
        const Evaluation& invMugBg = inverseSaturatedGasBMu_[regionIdx].eval(pressure, /*extrapolate=*/true);

        return invBg/invMugBg;
    }

    /*!
     * \brief Returns the density [kg/m^3] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation density(unsigned regionIdx,
                       const Evaluation& temperature,
                       const Evaluation& pressure,
                       const Evaluation& Rv) const
    {
        const Evaluation& Bg = formationVolumeFactor(regionIdx, temperature, pressure, Rv);

        Scalar rhooRef = oilReferenceDensity_[regionIdx];
        Scalar rhogRef = gasReferenceDensity_[regionIdx];
        Evaluation rhog = rhogRef/Bg;

        // the oil formation volume factor just represents the partial density of the gas
        // component in the gas phase. to get the total density of the phase, we have to
        // add the partial density of the oil component.
        rhog += (rhooRef*Rv)/Bg;

        return rhog;
    }

    /*!
     * \brief Returns the density [kg/m^3] of oil saturated gas at a given pressure.
     */
    template <class Evaluation>
    Evaluation saturatedDensity(unsigned regionIdx,
                                const Evaluation& temperature,
                                const Evaluation& pressure) const
    {
        Scalar rhooRef = oilReferenceDensity_[regionIdx];
        Scalar rhogRef = gasReferenceDensity_[regionIdx];
        Valgrind::CheckDefined(rhooRef);
        Valgrind::CheckDefined(rhogRef);

        const Evaluation& Bg = saturatedFormationVolumeFactor(regionIdx, temperature, pressure);

        Evaluation rhog = rhogRef/Bg;
        const Evaluation& RvSat = saturatedOilVaporizationFactor(regionIdx, temperature, pressure);
        rhog += (rhooRef*RvSat)/Bg;

        return rhog;
    }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation formationVolumeFactor(unsigned regionIdx,
                                     const Evaluation& /*temperature*/,
                                     const Evaluation& pressure,
                                     const Evaluation& Rv) const
    { return 1.0 / inverseGasB_[regionIdx].eval(pressure, Rv, /*extrapolate=*/true); }

    /*!
     * \brief Returns the formation volume factor [-] of oil saturated gas at a given pressure.
     */
    template <class Evaluation>
    Evaluation saturatedFormationVolumeFactor(unsigned regionIdx,
                                              const Evaluation& /*temperature*/,
                                              const Evaluation& pressure) const
    { return 1.0 / inverseSaturatedGasB_[regionIdx].eval(pressure, /*extrapolate=*/true); }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of the oil phase.
     */
    template <class Evaluation>
    Evaluation saturatedOilVaporizationFactor(unsigned regionIdx,
                                              const Evaluation& /*temperature*/,
                                              const Evaluation& pressure) const
    { return saturatedOilVaporizationFactorTable_[regionIdx].eval(pressure, /*extrapolate=*/true); }

    /*!
     * \brief Returns the saturation pressure of the gas phase [Pa]
     *        depending on its mass fraction of the oil component
     *
     * \param Rv The surface volume of oil component dissolved in what will yield one cubic meter of gas at the surface [-]
     */
    template <class Evaluation>
    Evaluation saturationPressure(unsigned regionIdx,
                                  const Evaluation& temperature,
                                  const Evaluation& Rv) const
    {
        typedef Opm::MathToolbox<Evaluation> Toolbox;

        // use the saturation pressure spline to get a pretty good initial value
        Evaluation pSat = saturationPressureSpline_[regionIdx].eval(Rv, /*extrapolate=*/true);
        const Evaluation& eps = pSat*1e-11;

        // Newton method to do the remaining work. If the initial
        // value is good, this should only take two to three
        // iterations...
        for (unsigned i = 0; i < 20; ++i) {
            const Evaluation& f = saturatedOilVaporizationFactor(regionIdx, temperature, pSat) - Rv;
            const Evaluation& fPrime = ((saturatedOilVaporizationFactor(regionIdx, temperature, pSat + eps) - Rv) - f)/eps;

            const Evaluation& delta = f/fPrime;
            pSat -= delta;

            if (std::abs(Toolbox::value(delta)) < std::abs(Toolbox::value(pSat)) * 1e-10)
                return pSat;
        }

        OPM_THROW(NumericalProblem, "Could find the gas saturation pressure for X_g^O = " << Rv);
    }

private:
    void updateSaturationPressureSpline_(unsigned regionIdx)
    {
        auto& oilVaporizationFac = saturatedOilVaporizationFactorTable_[regionIdx];

        // create the spline representing saturation pressure
        // depending of the mass fraction in gas
        size_t n = oilVaporizationFac.numSamples()*5;
        Scalar delta = (oilVaporizationFac.xMax() - oilVaporizationFac.xMin())/(n + 1);

        SamplingPoints pSatSamplePoints;
        Scalar Rv = 0;
        for (size_t i = 0; i <= n; ++ i) {
            Scalar pSat = oilVaporizationFac.xMin() + i*delta;
            Rv = saturatedOilVaporizationFactor(regionIdx, /*temperature=*/Scalar(1e100), pSat);

            std::pair<Scalar, Scalar> val(Rv, pSat);
            pSatSamplePoints.push_back(val);
        }
        saturationPressureSpline_[regionIdx].setContainerOfTuples(pSatSamplePoints,
                                                                  /*type=*/Spline::Monotonic);
    }

    std::vector<Scalar> gasReferenceDensity_;
    std::vector<Scalar> oilReferenceDensity_;
    std::vector<TabulatedTwoDFunction> inverseGasB_;
    std::vector<TabulatedOneDFunction> inverseSaturatedGasB_;
    std::vector<TabulatedTwoDFunction> gasMu_;
    std::vector<TabulatedTwoDFunction> inverseGasBMu_;
    std::vector<TabulatedOneDFunction> inverseSaturatedGasBMu_;
    std::vector<TabulatedOneDFunction> saturatedOilVaporizationFactorTable_;
    std::vector<Spline> saturationPressureSpline_;
};

} // namespace Opm

#endif
