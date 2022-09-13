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

#include <opm/material/common/MathToolbox.hpp>

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#endif

#if HAVE_OPM_COMMON
#include <opm/common/OpmLog/OpmLog.hpp>
#endif

namespace Opm {

template<class Scalar>
WetHumidGasPvt<Scalar>::
WetHumidGasPvt(const std::vector<Scalar>& gasReferenceDensity,
               const std::vector<Scalar>& oilReferenceDensity,
               const std::vector<Scalar>& waterReferenceDensity,
               const std::vector<TabulatedTwoDFunction>& inverseGasBRvwSat,
               const std::vector<TabulatedTwoDFunction>& inverseGasBRvSat,
               const std::vector<TabulatedOneDFunction>& inverseSaturatedGasB,
               const std::vector<TabulatedTwoDFunction>& gasMuRvwSat,
               const std::vector<TabulatedTwoDFunction>& gasMuRvSat,
               const std::vector<TabulatedTwoDFunction>& inverseGasBMuRvwSat,
               const std::vector<TabulatedTwoDFunction>& inverseGasBMuRvSat,
               const std::vector<TabulatedOneDFunction>& inverseSaturatedGasBMu,
               const std::vector<TabulatedOneDFunction>& saturatedWaterVaporizationFactorTable,
               const std::vector<TabulatedOneDFunction>& saturatedOilVaporizationFactorTable,
               const std::vector<TabulatedOneDFunction>& saturationPressure,
               Scalar vapPar1)
    : gasReferenceDensity_(gasReferenceDensity)
    , oilReferenceDensity_(oilReferenceDensity)
    , waterReferenceDensity_(waterReferenceDensity)
    , inverseGasBRvwSat_(inverseGasBRvwSat) // inverse of Bg evaluated at saturated water-gas ratio (Rvw) values; pvtg
    , inverseGasBRvSat_(inverseGasBRvSat) // inverse of Bg evaluated at saturated oil-gas ratio (Rv) values; pvtgw
    , inverseSaturatedGasB_(inverseSaturatedGasB) // evaluated at saturated water-gas ratio (Rvw) and oil-gas ratio (Rv) values; pvtgw
    , gasMuRvwSat_(gasMuRvwSat) // Mug evaluated at saturated water-gas ratio (Rvw) values; pvtg
    , gasMuRvSat_(gasMuRvSat) // Mug evaluated at saturated oil-gas ratio (Rv) values; pvtgw
    , inverseGasBMuRvwSat_(inverseGasBMuRvwSat) // Bg^-1*Mug evaluated at saturated water-gas ratio (Rvw) values; pvtg
    , inverseGasBMuRvSat_(inverseGasBMuRvSat) // Bg^-1*Mug evaluated at saturated oil-gas ratio (Rv) values; pvtgw
    , inverseSaturatedGasBMu_(inverseSaturatedGasBMu) //pvtgw
    , saturatedWaterVaporizationFactorTable_(saturatedWaterVaporizationFactorTable) //pvtgw
    , saturatedOilVaporizationFactorTable_(saturatedOilVaporizationFactorTable) //pvtg
    , saturationPressure_(saturationPressure)
    , vapPar1_(vapPar1)
{
}

#if HAVE_ECL_INPUT
template<class Scalar>
void WetHumidGasPvt<Scalar>::
initFromState(const EclipseState& eclState, const Schedule& schedule)
{
    const auto& pvtgwTables = eclState.getTableManager().getPvtgwTables();
    const auto& pvtgTables = eclState.getTableManager().getPvtgTables();
    const auto& densityTable = eclState.getTableManager().getDensityTable();

    assert(pvtgwTables.size() == densityTable.size());
    assert(pvtgTables.size() == densityTable.size());

    size_t numRegions = pvtgwTables.size();
    setNumRegions(numRegions);

    for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
        Scalar rhoRefO = densityTable[regionIdx].oil;
        Scalar rhoRefG = densityTable[regionIdx].gas;
        Scalar rhoRefW = densityTable[regionIdx].water;

        setReferenceDensities(regionIdx, rhoRefO, rhoRefG, rhoRefW);
    }
    // Table PVTGW
    for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
        const auto& pvtgwTable = pvtgwTables[regionIdx];

        const auto& saturatedTable = pvtgwTable.getSaturatedTable();
        assert(saturatedTable.numRows() > 1);

        // PVTGW table contains values at saturated Rv
        auto& gasMuRvSat = gasMuRvSat_[regionIdx];
        auto& invGasBRvSat = inverseGasBRvSat_[regionIdx];
        auto& invSatGasB = inverseSaturatedGasB_[regionIdx];
        auto& invSatGasBMu = inverseSaturatedGasBMu_[regionIdx];
        auto& waterVaporizationFac = saturatedWaterVaporizationFactorTable_[regionIdx];

        waterVaporizationFac.setXYArrays(saturatedTable.numRows(),
                                       saturatedTable.getColumn("PG"),
                                       saturatedTable.getColumn("RW"));

        std::vector<Scalar> invSatGasBArray;
        std::vector<Scalar> invSatGasBMuArray;

        // extract the table for the gas viscosity and formation volume factors
        for (unsigned outerIdx = 0; outerIdx < saturatedTable.numRows(); ++ outerIdx) {
            Scalar pg = saturatedTable.get("PG" , outerIdx);
            Scalar B = saturatedTable.get("BG" , outerIdx);
            Scalar mu = saturatedTable.get("MUG" , outerIdx);

            invGasBRvSat.appendXPos(pg);
            gasMuRvSat.appendXPos(pg);

            invSatGasBArray.push_back(1.0/B);
            invSatGasBMuArray.push_back(1.0/(mu*B));

            assert(invGasBRvSat.numX() == outerIdx + 1);
            assert(gasMuRvSat.numX() == outerIdx + 1);

            const auto& underSaturatedTable = pvtgwTable.getUnderSaturatedTable(outerIdx);
            size_t numRows = underSaturatedTable.numRows();
            for (size_t innerIdx = 0; innerIdx < numRows; ++ innerIdx) {
                Scalar Rw = underSaturatedTable.get("RW" , innerIdx);
                Scalar Bg = underSaturatedTable.get("BG" , innerIdx);
                Scalar mug = underSaturatedTable.get("MUG" , innerIdx);

                invGasBRvSat.appendSamplePoint(outerIdx, Rw, 1.0/Bg);
                gasMuRvSat.appendSamplePoint(outerIdx, Rw, mug);
            }
        }

        {
            std::vector<double> tmpPressure =  saturatedTable.getColumn("PG").vectorCopy( );

            invSatGasB.setXYContainers(tmpPressure, invSatGasBArray);
            invSatGasBMu.setXYContainers(tmpPressure, invSatGasBMuArray);
        }

        // make sure to have at least two sample points per gas pressure value
        for (unsigned xIdx = 0; xIdx < invGasBRvSat.numX(); ++xIdx) {
           // a single sample point is definitely needed
            assert(invGasBRvSat.numY(xIdx) > 0);

            // everything is fine if the current table has two or more sampling points
            // for a given mole fraction
            if (invGasBRvSat.numY(xIdx) > 1)
                continue;

            // find the master table which will be used as a template to extend the
            // current line. We define master table as the first table which has values
            // for undersaturated gas...
            size_t masterTableIdx = xIdx + 1;
            for (; masterTableIdx < saturatedTable.numRows(); ++masterTableIdx)
            {
                if (pvtgwTable.getUnderSaturatedTable(masterTableIdx).numRows() > 1)
                    break;
            }

            if (masterTableIdx >= saturatedTable.numRows())
                throw std::runtime_error("PVTGW tables are invalid: The last table must exhibit at least one "
                          "entry for undersaturated gas!");


            // extend the current table using the master table.
            extendPvtgwTable_(regionIdx,
                             xIdx,
                             pvtgwTable.getUnderSaturatedTable(xIdx),
                             pvtgwTable.getUnderSaturatedTable(masterTableIdx));
        }
    }

    // Table PVTG
    for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
        const auto& pvtgTable = pvtgTables[regionIdx];

        const auto& saturatedTable = pvtgTable.getSaturatedTable();
        assert(saturatedTable.numRows() > 1);
        // PVTG table contains values at saturated Rvw
        auto& gasMuRvwSat = gasMuRvwSat_[regionIdx];
        auto& invGasBRvwSat = inverseGasBRvwSat_[regionIdx];
        auto& invSatGasB = inverseSaturatedGasB_[regionIdx];
        auto& invSatGasBMu = inverseSaturatedGasBMu_[regionIdx];
        auto& oilVaporizationFac = saturatedOilVaporizationFactorTable_[regionIdx];

        oilVaporizationFac.setXYArrays(saturatedTable.numRows(),
                                       saturatedTable.getColumn("PG"),
                                       saturatedTable.getColumn("RV"));

        std::vector<Scalar> invSatGasBArray;
        std::vector<Scalar> invSatGasBMuArray;

        //// extract the table for the gas viscosity and formation volume factors
        for (unsigned outerIdx = 0; outerIdx < saturatedTable.numRows(); ++ outerIdx) {
            Scalar pg = saturatedTable.get("PG" , outerIdx);
            Scalar B = saturatedTable.get("BG" , outerIdx);
            Scalar mu = saturatedTable.get("MUG" , outerIdx);

            invGasBRvwSat.appendXPos(pg);
            gasMuRvwSat.appendXPos(pg);

            invSatGasBArray.push_back(1.0/B);
            invSatGasBMuArray.push_back(1.0/(mu*B));

            assert(invGasBRvwSat.numX() == outerIdx + 1);
            assert(gasMuRvwSat.numX() == outerIdx + 1);

            const auto& underSaturatedTable = pvtgTable.getUnderSaturatedTable(outerIdx);
            size_t numRows = underSaturatedTable.numRows();
            for (size_t innerIdx = 0; innerIdx < numRows; ++ innerIdx) {
                Scalar Rv = underSaturatedTable.get("RV" , innerIdx);
                Scalar Bg = underSaturatedTable.get("BG" , innerIdx);
                Scalar mug = underSaturatedTable.get("MUG" , innerIdx);

                invGasBRvwSat.appendSamplePoint(outerIdx, Rv, 1.0/Bg);
                gasMuRvwSat.appendSamplePoint(outerIdx, Rv, mug);
            }
        }

        {
            std::vector<double> tmpPressure =  saturatedTable.getColumn("PG").vectorCopy( );

            invSatGasB.setXYContainers(tmpPressure, invSatGasBArray);
            invSatGasBMu.setXYContainers(tmpPressure, invSatGasBMuArray);
        }

        // make sure to have at least two sample points per gas pressure value
        for (unsigned xIdx = 0; xIdx < invGasBRvwSat.numX(); ++xIdx) {
           // a single sample point is definitely needed
            assert(invGasBRvwSat.numY(xIdx) > 0);

            // everything is fine if the current table has two or more sampling points
            // for a given mole fraction
            if (invGasBRvwSat.numY(xIdx) > 1)
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
                throw std::runtime_error("PVTG tables are invalid: The last table must exhibit at least one "
                          "entry for undersaturated gas!");


            // extend the current table using the master table.
            extendPvtgTable_(regionIdx,
                             xIdx,
                             pvtgTable.getUnderSaturatedTable(xIdx),
                             pvtgTable.getUnderSaturatedTable(masterTableIdx));
        }
    } //end PVTGW
    vapPar1_ = 0.0;
    const auto& oilVap = schedule[0].oilvap();
    if (oilVap.getType() == OilVaporizationProperties::OilVaporization::VAPPARS) {
        vapPar1_ = oilVap.vap1();
    }

    initEnd();
}

template<class Scalar>
void WetHumidGasPvt<Scalar>::
extendPvtgwTable_(unsigned regionIdx,
                  unsigned xIdx,
                  const SimpleTable& curTable,
                  const SimpleTable& masterTable)
{
    std::vector<double> RvArray = curTable.getColumn("RW").vectorCopy();
    std::vector<double> gasBArray = curTable.getColumn("BG").vectorCopy();
    std::vector<double> gasMuArray = curTable.getColumn("MUG").vectorCopy();

    auto& invGasBRvSat = inverseGasBRvSat_[regionIdx];
    auto& gasMuRvSat = gasMuRvSat_[regionIdx];

    for (size_t newRowIdx = 1; newRowIdx < masterTable.numRows(); ++ newRowIdx) {
        const auto& RVColumn = masterTable.getColumn("RW");
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
        invGasBRvSat.appendSamplePoint(xIdx, newRv, 1.0/newBg);
        gasMuRvSat.appendSamplePoint(xIdx, newRv, newMug);
    }
}

template<class Scalar>
void WetHumidGasPvt<Scalar>::
extendPvtgTable_(unsigned regionIdx,
                 unsigned xIdx,
                 const SimpleTable& curTable,
                 const SimpleTable& masterTable)
{
    std::vector<double> RvArray = curTable.getColumn("RV").vectorCopy();
    std::vector<double> gasBArray = curTable.getColumn("BG").vectorCopy();
    std::vector<double> gasMuArray = curTable.getColumn("MUG").vectorCopy();

    auto& invGasBRvwSat= inverseGasBRvwSat_[regionIdx];
    auto& gasMuRvwSat = gasMuRvwSat_[regionIdx];

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
        invGasBRvwSat.appendSamplePoint(xIdx, newRv, 1.0/newBg);
        gasMuRvwSat.appendSamplePoint(xIdx, newRv, newMug);
    }
}
#endif

template<class Scalar>
void WetHumidGasPvt<Scalar>::setNumRegions(size_t numRegions)
{
    waterReferenceDensity_.resize(numRegions);
    oilReferenceDensity_.resize(numRegions);
    gasReferenceDensity_.resize(numRegions);
    inverseGasBRvwSat_.resize(numRegions, TabulatedTwoDFunction{TabulatedTwoDFunction::InterpolationPolicy::RightExtreme});
    inverseGasBRvSat_.resize(numRegions, TabulatedTwoDFunction{TabulatedTwoDFunction::InterpolationPolicy::RightExtreme});
    inverseGasBMuRvwSat_.resize(numRegions, TabulatedTwoDFunction{TabulatedTwoDFunction::InterpolationPolicy::RightExtreme});
    inverseGasBMuRvSat_.resize(numRegions, TabulatedTwoDFunction{TabulatedTwoDFunction::InterpolationPolicy::RightExtreme});
    inverseSaturatedGasB_.resize(numRegions);
    inverseSaturatedGasBMu_.resize(numRegions);
    gasMuRvwSat_.resize(numRegions, TabulatedTwoDFunction{TabulatedTwoDFunction::InterpolationPolicy::RightExtreme});
    gasMuRvSat_.resize(numRegions, TabulatedTwoDFunction{TabulatedTwoDFunction::InterpolationPolicy::RightExtreme});
    saturatedWaterVaporizationFactorTable_.resize(numRegions);
    saturatedOilVaporizationFactorTable_.resize(numRegions);
    saturationPressure_.resize(numRegions);
}

template<class Scalar>
void WetHumidGasPvt<Scalar>::setReferenceDensities(unsigned regionIdx,
                                                   Scalar rhoRefOil,
                                                   Scalar rhoRefGas,
                                                   Scalar rhoRefWater)
{
    waterReferenceDensity_[regionIdx] = rhoRefWater;
    oilReferenceDensity_[regionIdx] = rhoRefOil;
    gasReferenceDensity_[regionIdx] = rhoRefGas;
}

template<class Scalar>
void WetHumidGasPvt<Scalar>::
setSaturatedGasWaterVaporizationFactor(unsigned regionIdx,
                                       const SamplingPoints& samplePoints)
{
    saturatedWaterVaporizationFactorTable_[regionIdx].setContainerOfTuples(samplePoints);
}

template<class Scalar>
void WetHumidGasPvt<Scalar>::
setSaturatedGasOilVaporizationFactor(unsigned regionIdx,
                                     const SamplingPoints& samplePoints)
{
    saturatedOilVaporizationFactorTable_[regionIdx].setContainerOfTuples(samplePoints);
}

template<class Scalar>
void WetHumidGasPvt<Scalar>::initEnd()
{

    //PVTGW
    // calculate the final 2D functions which are used for interpolation.
    size_t numRegions = gasMuRvSat_.size();
    for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
        // calculate the table which stores the inverse of the product of the gas
        // formation volume factor and the gas viscosity
        const auto& gasMuRvSat = gasMuRvSat_[regionIdx];
        const auto& invGasBRvSat = inverseGasBRvSat_[regionIdx];
        assert(gasMuRvSat.numX() == invGasBRvSat.numX());

        auto& invGasBMuRvSat = inverseGasBMuRvSat_[regionIdx];
        auto& invSatGasB = inverseSaturatedGasB_[regionIdx];
        auto& invSatGasBMu = inverseSaturatedGasBMu_[regionIdx];

        std::vector<Scalar> satPressuresArray;
        std::vector<Scalar> invSatGasBArray;
        std::vector<Scalar> invSatGasBMuArray;
        for (size_t pIdx = 0; pIdx < gasMuRvSat.numX(); ++pIdx) {
            invGasBMuRvSat.appendXPos(gasMuRvSat.xAt(pIdx));

            assert(gasMuRvSat.numY(pIdx) == invGasBRvSat.numY(pIdx));

            size_t numRw = gasMuRvSat.numY(pIdx);
            for (size_t RwIdx = 0; RwIdx < numRw; ++RwIdx)
                invGasBMuRvSat.appendSamplePoint(pIdx,
                                            gasMuRvSat.yAt(pIdx, RwIdx),
                                            invGasBRvSat.valueAt(pIdx, RwIdx)
                                            / gasMuRvSat.valueAt(pIdx, RwIdx));

            // the sampling points in UniformXTabulated2DFunction are always sorted
            // in ascending order. Thus, the value for saturated gas is the last one
            // (i.e., the one with the largest Rw value)
            satPressuresArray.push_back(gasMuRvSat.xAt(pIdx));
            invSatGasBArray.push_back(invGasBRvSat.valueAt(pIdx, numRw - 1));
            invSatGasBMuArray.push_back(invGasBMuRvSat.valueAt(pIdx, numRw - 1));
        }

        invSatGasB.setXYContainers(satPressuresArray, invSatGasBArray);
        invSatGasBMu.setXYContainers(satPressuresArray, invSatGasBMuArray);
    }

    //PVTG
    // calculate the final 2D functions which are used for interpolation.
    //size_t numRegions = gasMuRvwSat_.size();
    for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
        // calculate the table which stores the inverse of the product of the gas
        // formation volume factor and the gas viscosity
        const auto& gasMuRvwSat = gasMuRvwSat_[regionIdx];
        const auto& invGasBRvwSat = inverseGasBRvwSat_[regionIdx];
        assert(gasMuRvwSat.numX() == invGasBRvwSat.numX());

        auto& invGasBMuRvwSat = inverseGasBMuRvwSat_[regionIdx];
        auto& invSatGasB = inverseSaturatedGasB_[regionIdx];
        auto& invSatGasBMu = inverseSaturatedGasBMu_[regionIdx];

        std::vector<Scalar> satPressuresArray;
        std::vector<Scalar> invSatGasBArray;
        std::vector<Scalar> invSatGasBMuArray;
        for (size_t pIdx = 0; pIdx < gasMuRvwSat.numX(); ++pIdx) {
            invGasBMuRvwSat.appendXPos(gasMuRvwSat.xAt(pIdx));

            assert(gasMuRvwSat.numY(pIdx) == invGasBRvwSat.numY(pIdx));

            size_t numRw = gasMuRvwSat.numY(pIdx);
            for (size_t RwIdx = 0; RwIdx < numRw; ++RwIdx)
                invGasBMuRvwSat.appendSamplePoint(pIdx,
                                            gasMuRvwSat.yAt(pIdx, RwIdx),
                                            invGasBRvwSat.valueAt(pIdx, RwIdx)
                                            / gasMuRvwSat.valueAt(pIdx, RwIdx));

            // the sampling points in UniformXTabulated2DFunction are always sorted
            // in ascending order. Thus, the value for saturated gas is the last one
            // (i.e., the one with the largest Rw value)
            satPressuresArray.push_back(gasMuRvwSat.xAt(pIdx));
            invSatGasBArray.push_back(invGasBRvwSat.valueAt(pIdx, numRw - 1));
            invSatGasBMuArray.push_back(invGasBMuRvwSat.valueAt(pIdx, numRw - 1));
        }

        invSatGasB.setXYContainers(satPressuresArray, invSatGasBArray);
        invSatGasBMu.setXYContainers(satPressuresArray, invSatGasBMuArray);

        updateSaturationPressure_(regionIdx);
    }
}

template<class Scalar>
bool WetHumidGasPvt<Scalar>::operator==(const WetHumidGasPvt<Scalar>& data) const
{
    return this->gasReferenceDensity_ == data.gasReferenceDensity_ &&
           this->oilReferenceDensity_ == data.oilReferenceDensity_ &&
           this->waterReferenceDensity_ == data.waterReferenceDensity_ &&
           this->inverseGasB() == data.inverseGasB() &&
           this->inverseSaturatedGasB() == data.inverseSaturatedGasB() &&
           this->gasMu() == data.gasMu() &&
           this->inverseGasBMu() == data.inverseGasBMu() &&
           this->inverseSaturatedGasBMu() == data.inverseSaturatedGasBMu() &&
           this->saturatedWaterVaporizationFactorTable() == data.saturatedWaterVaporizationFactorTable() &&
           this->saturationPressure() == data.saturationPressure() &&
           this->vapPar1() == data.vapPar1();
}

template<class Scalar>
void WetHumidGasPvt<Scalar>::updateSaturationPressure_(unsigned regionIdx)
{
    const auto& oilVaporizationFac = saturatedOilVaporizationFactorTable_[regionIdx];

    // create the taublated function representing saturation pressure depending of
    // Rv
    size_t n = oilVaporizationFac.numSamples();
    Scalar delta = (oilVaporizationFac.xMax() - oilVaporizationFac.xMin())/Scalar(n + 1);

    SamplingPoints pSatSamplePoints;
    Scalar Rv = 0;
    for (size_t i = 0; i <= n; ++ i) {
        Scalar pSat = oilVaporizationFac.xMin() + Scalar(i)*delta;
        Rv = saturatedOilVaporizationFactor(regionIdx, /*temperature=*/Scalar(1e30), pSat);

        pSatSamplePoints.emplace_back(Rv, pSat);
    }

    //Prune duplicate Rv values (can occur, and will cause problems in further interpolation)
    auto x_coord_comparator = [](const auto& a, const auto& b) { return a.first == b.first; };
    auto last = std::unique(pSatSamplePoints.begin(), pSatSamplePoints.end(), x_coord_comparator);
    pSatSamplePoints.erase(last, pSatSamplePoints.end());

    saturationPressure_[regionIdx].setContainerOfTuples(pSatSamplePoints);
}

template<class Scalar>
template<class Evaluation>
Evaluation WetHumidGasPvt<Scalar>::internalEnergy(unsigned,
                                                  const Evaluation&,
                                                  const Evaluation&,
                                                  const Evaluation&) const
{
    throw std::runtime_error("Requested the enthalpy of gas but the thermal option is not enabled");
}

template<class Scalar>
template<class Evaluation>
Evaluation WetHumidGasPvt<Scalar>::viscosity(unsigned regionIdx,
                                             const Evaluation& /*temperature*/,
                                             const Evaluation& pressure,
                                             const Evaluation& Rv,
                                             const Evaluation& Rvw) const
{
    const Evaluation& temperature = 1E30;

    if (Rv >= (1.0 - 1e-10)*saturatedOilVaporizationFactor(regionIdx, temperature, pressure)) {
        const Evaluation& invBg = inverseGasBRvSat_[regionIdx].eval(pressure, Rvw, /*extrapolate=*/true);
        const Evaluation& invMugBg = inverseGasBMuRvSat_[regionIdx].eval(pressure, Rvw, /*extrapolate=*/true);
        return invBg/invMugBg;
    }
    else {
        // for Rv undersaturated viscosity is evaluated at saturated Rvw values
        const Evaluation& invBg = inverseGasBRvwSat_[regionIdx].eval(pressure, Rv, /*extrapolate=*/true);
        const Evaluation& invMugBg = inverseGasBMuRvwSat_[regionIdx].eval(pressure, Rv, /*extrapolate=*/true);
        return invBg/invMugBg;
    }
}

template<class Scalar>
template<class Evaluation>
Evaluation WetHumidGasPvt<Scalar>::
saturatedViscosity(unsigned regionIdx,
                   const Evaluation& /*temperature*/,
                   const Evaluation& pressure) const
{
    const Evaluation& invBg = inverseSaturatedGasB_[regionIdx].eval(pressure, /*extrapolate=*/true);
    const Evaluation& invMugBg = inverseSaturatedGasBMu_[regionIdx].eval(pressure, /*extrapolate=*/true);

    return invBg/invMugBg;
}

template<class Scalar>
template<class Evaluation>
Evaluation WetHumidGasPvt<Scalar>::
inverseFormationVolumeFactor(unsigned regionIdx,
                             const Evaluation& /*temperature*/,
                             const Evaluation& pressure,
                             const Evaluation& Rv,
                             const Evaluation& Rvw) const
{
    const Evaluation& temperature = 1E30;

    if (Rv >= (1.0 - 1e-10)*saturatedOilVaporizationFactor(regionIdx, temperature, pressure)) {
        return inverseGasBRvSat_[regionIdx].eval(pressure, Rvw, /*extrapolate=*/true);
    }
    else {
        // for Rv undersaturated Bg^-1 is evaluated at saturated Rvw values
        return inverseGasBRvwSat_[regionIdx].eval(pressure, Rv, /*extrapolate=*/true);
    }

}

template<class Scalar>
template<class Evaluation>
Evaluation WetHumidGasPvt<Scalar>::
saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                      const Evaluation& /*temperature*/,
                                      const Evaluation& pressure) const
{
    return inverseSaturatedGasB_[regionIdx].eval(pressure, /*extrapolate=*/true);
}

template<class Scalar>
template<class Evaluation>
Evaluation WetHumidGasPvt<Scalar>::
saturatedWaterVaporizationFactor(unsigned regionIdx,
                                 const Evaluation& /*temperature*/,
                                 const Evaluation& pressure) const
{
    return saturatedWaterVaporizationFactorTable_[regionIdx].eval(pressure, /*extrapolate=*/true);
}

template<class Scalar>
template<class Evaluation>
Evaluation WetHumidGasPvt<Scalar>::
saturatedOilVaporizationFactor(unsigned regionIdx,
                               const Evaluation& /*temperature*/,
                               const Evaluation& pressure) const
{
    return saturatedOilVaporizationFactorTable_[regionIdx].eval(pressure, /*extrapolate=*/true);
}

template<class Scalar>
template<class Evaluation>
Evaluation WetHumidGasPvt<Scalar>::
saturatedOilVaporizationFactor(unsigned regionIdx,
                               const Evaluation& /*temperature*/,
                               const Evaluation& pressure,
                               const Evaluation& oilSaturation,
                               Evaluation maxOilSaturation) const
{
    Evaluation tmp =
        saturatedOilVaporizationFactorTable_[regionIdx].eval(pressure, /*extrapolate=*/true);

    // apply the vaporization parameters for the gas phase (cf. the Eclipse VAPPARS
    // keyword)
    maxOilSaturation = min(maxOilSaturation, Scalar(1.0));
    if (vapPar1_ > 0.0 && maxOilSaturation > 0.01 && oilSaturation < maxOilSaturation) {
        constexpr const Scalar eps = 0.001;
        const Evaluation& So = max(oilSaturation, eps);
        tmp *= max(1e-3, pow(So/maxOilSaturation, vapPar1_));
    }

    return tmp;
}

template<class Scalar>
template<class Evaluation>
Evaluation WetHumidGasPvt<Scalar>::
saturationPressure(unsigned regionIdx,
                   const Evaluation&,
                   const Evaluation& Rw) const
{
    using Toolbox = MathToolbox<Evaluation>;

    const auto& RwTable = saturatedWaterVaporizationFactorTable_[regionIdx];
    constexpr const Scalar eps = std::numeric_limits<typename Toolbox::Scalar>::epsilon()*1e6;

    // use the tabulated saturation pressure function to get a pretty good initial value
    Evaluation pSat = saturationPressure_[regionIdx].eval(Rw, /*extrapolate=*/true);

    // Newton method to do the remaining work. If the initial
    // value is good, this should only take two to three
    // iterations...
    bool onProbation = false;
    for (unsigned i = 0; i < 20; ++i) {
        const Evaluation& f = RwTable.eval(pSat, /*extrapolate=*/true) - Rw;
        const Evaluation& fPrime = RwTable.evalDerivative(pSat, /*extrapolate=*/true);

        // If the derivative is "zero" Newton will not converge,
        // so simply return our initial guess.
        if (std::abs(scalarValue(fPrime)) < 1.0e-30) {
            return pSat;
        }

        const Evaluation& delta = f/fPrime;

        pSat -= delta;

        if (pSat < 0.0) {
            // if the pressure is lower than 0 Pascals, we set it back to 0. if this
            // happens twice, we give up and just return 0 Pa...
            if (onProbation)
                return 0.0;

            onProbation = true;
            pSat = 0.0;
        }

        if (std::abs(scalarValue(delta)) < std::abs(scalarValue(pSat))*eps)
            return pSat;
    }

    std::stringstream errlog;
    errlog << "Finding saturation pressure did not converge:"
           << " pSat = " << pSat
           << ", Rw = " << Rw;
#if HAVE_OPM_COMMON
    OpmLog::debug("Wet gas saturation pressure", errlog.str());
#endif
    throw NumericalIssue(errlog.str());
}

template<class Scalar>
template<class Evaluation>
Evaluation WetHumidGasPvt<Scalar>::
diffusionCoefficient(const Evaluation& /*temperature*/,
                     const Evaluation& /*pressure*/,
                     unsigned /*compIdx*/) const
{
    throw std::runtime_error("Not implemented: The PVT model does not provide a diffusionCoefficient()");
}

} // namespace Opm
