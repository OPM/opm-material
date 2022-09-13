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
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/input/eclipse/Schedule/OilVaporizationProperties.hpp>
#endif

#if HAVE_OPM_COMMON
#include <opm/common/OpmLog/OpmLog.hpp>
#endif

namespace Opm {

template<class Scalar>
WetGasPvt<Scalar>::WetGasPvt(const std::vector<Scalar>& gasReferenceDensity,
                             const std::vector<Scalar>& oilReferenceDensity,
                             const std::vector<TabulatedTwoDFunction>& inverseGasB,
                             const std::vector<TabulatedOneDFunction>& inverseSaturatedGasB,
                             const std::vector<TabulatedTwoDFunction>& gasMu,
                             const std::vector<TabulatedTwoDFunction>& inverseGasBMu,
                             const std::vector<TabulatedOneDFunction>& inverseSaturatedGasBMu,
                             const std::vector<TabulatedOneDFunction>& saturatedOilVaporizationFactorTable,
                             const std::vector<TabulatedOneDFunction>& saturationPressure,
                             Scalar vapPar1)
    : gasReferenceDensity_(gasReferenceDensity)
    , oilReferenceDensity_(oilReferenceDensity)
    , inverseGasB_(inverseGasB)
    , inverseSaturatedGasB_(inverseSaturatedGasB)
    , gasMu_(gasMu)
    , inverseGasBMu_(inverseGasBMu)
    , inverseSaturatedGasBMu_(inverseSaturatedGasBMu)
    , saturatedOilVaporizationFactorTable_(saturatedOilVaporizationFactorTable)
    , saturationPressure_(saturationPressure)
    , vapPar1_(vapPar1)
{
}

#if HAVE_ECL_INPUT
template<class Scalar>
void WetGasPvt<Scalar>::initFromState(const EclipseState& eclState,
                                      const Schedule& schedule)
{
    const auto& pvtgTables = eclState.getTableManager().getPvtgTables();
    const auto& densityTable = eclState.getTableManager().getDensityTable();

    assert(pvtgTables.size() == densityTable.size());

    size_t numRegions = pvtgTables.size();
    setNumRegions(numRegions);

    for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
        Scalar rhoRefO = densityTable[regionIdx].oil;
        Scalar rhoRefG = densityTable[regionIdx].gas;
        Scalar rhoRefW = densityTable[regionIdx].water;

        setReferenceDensities(regionIdx, rhoRefO, rhoRefG, rhoRefW);
    }

    for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
        const auto& pvtgTable = pvtgTables[regionIdx];

        const auto& saturatedTable = pvtgTable.getSaturatedTable();
        assert(saturatedTable.numRows() > 1);

        auto& gasMu = gasMu_[regionIdx];
        auto& invGasB = inverseGasB_[regionIdx];
        auto& invSatGasB = inverseSaturatedGasB_[regionIdx];
        auto& invSatGasBMu = inverseSaturatedGasBMu_[regionIdx];
        auto& oilVaporizationFac = saturatedOilVaporizationFactorTable_[regionIdx];

        oilVaporizationFac.setXYArrays(saturatedTable.numRows(),
                                       saturatedTable.getColumn("PG"),
                                       saturatedTable.getColumn("RV"));

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

            const auto& underSaturatedTable = pvtgTable.getUnderSaturatedTable(outerIdx);
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
                throw std::runtime_error("PVTG tables are invalid: The last table must exhibit at least one "
                          "entry for undersaturated gas!");


            // extend the current table using the master table.
            extendPvtgTable_(regionIdx,
                             xIdx,
                             pvtgTable.getUnderSaturatedTable(xIdx),
                             pvtgTable.getUnderSaturatedTable(masterTableIdx));
        }
    }

    vapPar1_ = 0.0;
    const auto& oilVap = schedule[0].oilvap();
    if (oilVap.getType() == OilVaporizationProperties::OilVaporization::VAPPARS) {
        vapPar1_ = oilVap.vap1();
    }

    initEnd();
}

template<class Scalar>
void WetGasPvt<Scalar>::extendPvtgTable_(unsigned regionIdx,
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
#endif // HAVE_ECL_INPUT

template<class Scalar>
void WetGasPvt<Scalar>::setNumRegions(size_t numRegions)
{
    oilReferenceDensity_.resize(numRegions);
    gasReferenceDensity_.resize(numRegions);
    inverseGasB_.resize(numRegions, TabulatedTwoDFunction{TabulatedTwoDFunction::InterpolationPolicy::RightExtreme});
    inverseGasBMu_.resize(numRegions, TabulatedTwoDFunction{TabulatedTwoDFunction::InterpolationPolicy::RightExtreme});
    inverseSaturatedGasB_.resize(numRegions);
    inverseSaturatedGasBMu_.resize(numRegions);
    gasMu_.resize(numRegions, TabulatedTwoDFunction{TabulatedTwoDFunction::InterpolationPolicy::RightExtreme});
    saturatedOilVaporizationFactorTable_.resize(numRegions);
    saturationPressure_.resize(numRegions);
}

template<class Scalar>
void WetGasPvt<Scalar>::setReferenceDensities(unsigned regionIdx,
                                              Scalar rhoRefOil,
                                              Scalar rhoRefGas,
                                              Scalar /*rhoRefWater*/)
{
    oilReferenceDensity_[regionIdx] = rhoRefOil;
    gasReferenceDensity_[regionIdx] = rhoRefGas;
}

template<class Scalar>
void WetGasPvt<Scalar>::
setSaturatedGasOilVaporizationFactor(unsigned regionIdx,
                                     const SamplingPoints& samplePoints)
{
    saturatedOilVaporizationFactorTable_[regionIdx].setContainerOfTuples(samplePoints);
}

template<class Scalar>
void WetGasPvt<Scalar>::
setSaturatedGasFormationVolumeFactor(unsigned regionIdx,
                                     const SamplingPoints& samplePoints)
{
    auto& invGasB = inverseGasB_[regionIdx];

    const auto& RvTable = saturatedOilVaporizationFactorTable_[regionIdx];

    constexpr const Scalar T = 273.15 + 15.56; // [K]

    constexpr const Scalar RvMin = 0.0;
    Scalar RvMax = RvTable.eval(saturatedOilVaporizationFactorTable_[regionIdx].xMax(), /*extrapolate=*/true);

    Scalar poMin = samplePoints.front().first;
    Scalar poMax = samplePoints.back().first;

    constexpr const size_t nRv = 20;
    size_t nP = samplePoints.size()*2;

    Scalar rhooRef = oilReferenceDensity_[regionIdx];

    TabulatedOneDFunction gasFormationVolumeFactor;
    gasFormationVolumeFactor.setContainerOfTuples(samplePoints);

    updateSaturationPressure_(regionIdx);

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
            Scalar BgSat = gasFormationVolumeFactor.eval(poSat, /*extrapolate=*/true);
            Scalar drhoo_dp = (1.1200 - 1.1189)/((5000 - 4000)*6894.76);
            Scalar rhoo = rhooRef/BgSat*(1 + drhoo_dp*(pg - poSat));

            Scalar Bg = rhooRef/rhoo;

            invGasB.appendSamplePoint(RvIdx, pg, 1.0/Bg);
        }
    }
}

template<class Scalar>
void WetGasPvt<Scalar>::
setInverseGasFormationVolumeFactor(unsigned regionIdx,
                                   const TabulatedTwoDFunction& invBg)
{
    inverseGasB_[regionIdx] = invBg;
}

template<class Scalar>
void WetGasPvt<Scalar>::
setGasViscosity(unsigned regionIdx,
                const TabulatedTwoDFunction& mug)
{
    gasMu_[regionIdx] = mug;
}

template<class Scalar>
void WetGasPvt<Scalar>::
setSaturatedGasViscosity(unsigned regionIdx,
                         const SamplingPoints& samplePoints  )
{
    auto& oilVaporizationFac = saturatedOilVaporizationFactorTable_[regionIdx];

    constexpr const Scalar RvMin = 0.0;
    Scalar RvMax = oilVaporizationFac.eval(saturatedOilVaporizationFactorTable_[regionIdx].xMax(), /*extrapolate=*/true);

    Scalar poMin = samplePoints.front().first;
    Scalar poMax = samplePoints.back().first;

    constexpr const size_t nRv = 20;
    size_t nP = samplePoints.size()*2;

    TabulatedOneDFunction mugTable;
    mugTable.setContainerOfTuples(samplePoints);

    // calculate a table of estimated densities depending on pressure and gas mass
    // fraction
    for (size_t RvIdx = 0; RvIdx < nRv; ++RvIdx) {
        Scalar Rv = RvMin + (RvMax - RvMin)*RvIdx/nRv;

        gasMu_[regionIdx].appendXPos(Rv);

        for (size_t pIdx = 0; pIdx < nP; ++pIdx) {
            Scalar pg = poMin + (poMax - poMin)*pIdx/nP;
            Scalar mug = mugTable.eval(pg, /*extrapolate=*/true);

            gasMu_[regionIdx].appendSamplePoint(RvIdx, pg, mug);
        }
    }
}

template<class Scalar>
void WetGasPvt<Scalar>::initEnd()
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

        updateSaturationPressure_(regionIdx);
    }
}

template<class Scalar>
bool WetGasPvt<Scalar>::operator==(const WetGasPvt<Scalar>& data) const
{
    return this->gasReferenceDensity_ == data.gasReferenceDensity_ &&
           this->oilReferenceDensity_ == data.oilReferenceDensity_ &&
           this->inverseGasB() == data.inverseGasB() &&
           this->inverseSaturatedGasB() == data.inverseSaturatedGasB() &&
           this->gasMu() == data.gasMu() &&
           this->inverseGasBMu() == data.inverseGasBMu() &&
           this->inverseSaturatedGasBMu() == data.inverseSaturatedGasBMu() &&
           this->saturatedOilVaporizationFactorTable() == data.saturatedOilVaporizationFactorTable() &&
           this->saturationPressure() == data.saturationPressure() &&
           this->vapPar1() == data.vapPar1();
}

template<class Scalar>
void WetGasPvt<Scalar>::updateSaturationPressure_(unsigned regionIdx)
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
    if (std::distance(pSatSamplePoints.begin(), last) > 1) // only remove them if there are more than two points
        pSatSamplePoints.erase(last, pSatSamplePoints.end());

    saturationPressure_[regionIdx].setContainerOfTuples(pSatSamplePoints);
}

template<class Scalar>
template<class Evaluation>
Evaluation WetGasPvt<Scalar>::internalEnergy(unsigned,
                                             const Evaluation&,
                                             const Evaluation&,
                                             const Evaluation&) const
{
    throw std::runtime_error("Requested the enthalpy of gas but the thermal option is not enabled");
}

template<class Scalar>
template <class Evaluation>
Evaluation WetGasPvt<Scalar>::viscosity(unsigned regionIdx,
                                        const Evaluation& /*temperature*/,
                                        const Evaluation& pressure,
                                        const Evaluation& Rv,
                                        const Evaluation& /*Rvw*/) const
{
    const Evaluation& invBg = inverseGasB_[regionIdx].eval(pressure, Rv, /*extrapolate=*/true);
    const Evaluation& invMugBg = inverseGasBMu_[regionIdx].eval(pressure, Rv, /*extrapolate=*/true);

    return invBg/invMugBg;
}

template<class Scalar>
template<class Evaluation>
Evaluation WetGasPvt<Scalar>::saturatedViscosity(unsigned regionIdx,
                                                 const Evaluation& /*temperature*/,
                                                 const Evaluation& pressure) const
{
    const Evaluation& invBg = inverseSaturatedGasB_[regionIdx].eval(pressure, /*extrapolate=*/true);
    const Evaluation& invMugBg = inverseSaturatedGasBMu_[regionIdx].eval(pressure, /*extrapolate=*/true);

    return invBg/invMugBg;
}

template<class Scalar>
template<class Evaluation>
Evaluation WetGasPvt<Scalar>::
inverseFormationVolumeFactor(unsigned regionIdx,
                             const Evaluation& /*temperature*/,
                             const Evaluation& pressure,
                             const Evaluation& Rv,
                             const Evaluation& /*Rvw*/) const
{
    return inverseGasB_[regionIdx].eval(pressure, Rv, /*extrapolate=*/true);
}

template<class Scalar>
template<class Evaluation>
Evaluation WetGasPvt<Scalar>::
saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                      const Evaluation& /*temperature*/,
                                      const Evaluation& pressure) const
{
    return inverseSaturatedGasB_[regionIdx].eval(pressure, /*extrapolate=*/true);
}

template<class Scalar>
template<class Evaluation>
Evaluation WetGasPvt<Scalar>::
saturatedWaterVaporizationFactor(unsigned /*regionIdx*/,
                                 const Evaluation& /*temperature*/,
                                 const Evaluation& /*pressure*/) const
{
    return 0.0; /* this is non-humid gas! */
}

template<class Scalar>
template <class Evaluation>
Evaluation WetGasPvt<Scalar>::
saturatedOilVaporizationFactor(unsigned regionIdx,
                               const Evaluation& /*temperature*/,
                               const Evaluation& pressure) const
{
    return saturatedOilVaporizationFactorTable_[regionIdx].eval(pressure, /*extrapolate=*/true);
}

template<class Scalar>
template<class Evaluation>
Evaluation WetGasPvt<Scalar>::
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
Evaluation WetGasPvt<Scalar>::
saturationPressure(unsigned regionIdx,
                   const Evaluation&,
                   const Evaluation& Rv) const
{
    using Toolbox = MathToolbox<Evaluation>;

    const auto& RvTable = saturatedOilVaporizationFactorTable_[regionIdx];
    constexpr const Scalar eps = std::numeric_limits<typename Toolbox::Scalar>::epsilon()*1e6;

    // use the tabulated saturation pressure function to get a pretty good initial value
    Evaluation pSat = saturationPressure_[regionIdx].eval(Rv, /*extrapolate=*/true);

    // Newton method to do the remaining work. If the initial
    // value is good, this should only take two to three
    // iterations...
    bool onProbation = false;
    for (unsigned i = 0; i < 20; ++i) {
        const Evaluation& f = RvTable.eval(pSat, /*extrapolate=*/true) - Rv;
        const Evaluation& fPrime = RvTable.evalDerivative(pSat, /*extrapolate=*/true);

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
           << ", Rv = " << Rv;
#if HAVE_OPM_COMMON
    OpmLog::debug("Wet gas saturation pressure", errlog.str());
#endif
    throw NumericalIssue(errlog.str());
}

template<class Scalar>
template<class Evaluation>
Evaluation WetGasPvt<Scalar>::
diffusionCoefficient(const Evaluation& /*temperature*/,
                     const Evaluation& /*pressure*/,
                     unsigned /*compIdx*/) const
{
    throw std::runtime_error("Not implemented: The PVT model does not provide a diffusionCoefficient()");
}

} // namespace Opm
