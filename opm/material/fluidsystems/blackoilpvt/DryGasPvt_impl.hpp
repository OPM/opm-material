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

#include <opm/material/Constants.hpp>

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PvdgTable.hpp>
#endif

namespace Opm {

template<class Scalar>
DryGasPvt<Scalar>::DryGasPvt(const std::vector<Scalar>& gasReferenceDensity,
                             const std::vector<TabulatedOneDFunction>& inverseGasB,
                             const std::vector<TabulatedOneDFunction>& gasMu,
                             const std::vector<TabulatedOneDFunction>& inverseGasBMu)
  : gasReferenceDensity_(gasReferenceDensity)
  , inverseGasB_(inverseGasB)
  , gasMu_(gasMu)
  , inverseGasBMu_(inverseGasBMu)
{
}

#if HAVE_ECL_INPUT
template<class Scalar>
void DryGasPvt<Scalar>::initFromState(const EclipseState& eclState, const Schedule&)
{
    const auto& pvdgTables = eclState.getTableManager().getPvdgTables();
    const auto& densityTable = eclState.getTableManager().getDensityTable();

    assert(pvdgTables.size() == densityTable.size());

    size_t numRegions = pvdgTables.size();
    setNumRegions(numRegions);

    for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
        Scalar rhoRefO = densityTable[regionIdx].oil;
        Scalar rhoRefG = densityTable[regionIdx].gas;
        Scalar rhoRefW = densityTable[regionIdx].water;

        setReferenceDensities(regionIdx, rhoRefO, rhoRefG, rhoRefW);

        // determine the molar masses of the components
        constexpr Scalar p = 1.01325e5; // surface pressure, [Pa]
        constexpr Scalar T = 273.15 + 15.56; // surface temperature, [K]
        constexpr Scalar MO = 175e-3; // [kg/mol]
        Scalar MG = Constants<Scalar>::R*T*rhoRefG / p; // [kg/mol], consequence of the ideal gas law
        constexpr Scalar MW = 18.0e-3; // [kg/mol]
        // TODO (?): the molar mass of the components can possibly specified
        // explicitly in the deck.
        setMolarMasses(regionIdx, MO, MG, MW);

        const auto& pvdgTable = pvdgTables.getTable<PvdgTable>(regionIdx);

        // say 99.97% of all time: "premature optimization is the root of all
        // evil". Eclipse does this "optimization" for apparently no good reason!
        std::vector<Scalar> invB(pvdgTable.numRows());
        const auto& Bg = pvdgTable.getFormationFactorColumn();
        for (unsigned i = 0; i < Bg.size(); ++ i) {
            invB[i] = 1.0/Bg[i];
        }

        size_t numSamples = invB.size();
        inverseGasB_[regionIdx].setXYArrays(numSamples, pvdgTable.getPressureColumn(), invB);
        gasMu_[regionIdx].setXYArrays(numSamples, pvdgTable.getPressureColumn(), pvdgTable.getViscosityColumn());
    }

    initEnd();
}
#endif

template<class Scalar>
void DryGasPvt<Scalar>::setNumRegions(size_t numRegions)
{
    gasReferenceDensity_.resize(numRegions);
    inverseGasB_.resize(numRegions);
    inverseGasBMu_.resize(numRegions);
    gasMu_.resize(numRegions);
}

template<class Scalar>
void DryGasPvt<Scalar>::setReferenceDensities(unsigned regionIdx,
                                              Scalar /*rhoRefOil*/,
                                              Scalar rhoRefGas,
                                              Scalar /*rhoRefWater*/)
{
    gasReferenceDensity_[regionIdx] = rhoRefGas;
}

template<class Scalar>
void DryGasPvt<Scalar>::setGasViscosity(unsigned regionIdx,
                                        const TabulatedOneDFunction& mug)
{
    gasMu_[regionIdx] = mug;
}

template<class Scalar>
void DryGasPvt<Scalar>::setGasFormationVolumeFactor(unsigned regionIdx,
                                                    const SamplingPoints& samplePoints)
{
    SamplingPoints tmp(samplePoints);
    auto it = tmp.begin();
    const auto& endIt = tmp.end();
    for (; it != endIt; ++ it)
        std::get<1>(*it) = 1.0/std::get<1>(*it);

    inverseGasB_[regionIdx].setContainerOfTuples(tmp);
    assert(inverseGasB_[regionIdx].monotonic());
}

template<class Scalar>
void DryGasPvt<Scalar>::initEnd()
{
    // calculate the final 2D functions which are used for interpolation.
    size_t numRegions = gasMu_.size();
    for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
        // calculate the table which stores the inverse of the product of the gas
        // formation volume factor and the gas viscosity
        const auto& gasMu = gasMu_[regionIdx];
        const auto& invGasB = inverseGasB_[regionIdx];
        assert(gasMu.numSamples() == invGasB.numSamples());

        std::vector<Scalar> pressureValues(gasMu.numSamples());
        std::vector<Scalar> invGasBMuValues(gasMu.numSamples());
        for (unsigned pIdx = 0; pIdx < gasMu.numSamples(); ++pIdx) {
            pressureValues[pIdx] = invGasB.xAt(pIdx);
            invGasBMuValues[pIdx] = invGasB.valueAt(pIdx) * (1.0/gasMu.valueAt(pIdx));
        }

        inverseGasBMu_[regionIdx].setXYContainers(pressureValues, invGasBMuValues);
    }
}

template<class Scalar>
bool DryGasPvt<Scalar>::operator==(const DryGasPvt<Scalar>& data) const
{
  return gasReferenceDensity_ == data.gasReferenceDensity_ &&
      inverseGasB_ == data.inverseGasB_ &&
      gasMu_ == data.gasMu_ &&
      inverseGasBMu_ == data.inverseGasBMu_;
}

template<class Scalar>
template<class Evaluation>
Evaluation DryGasPvt<Scalar>::
internalEnergy(unsigned,
               const Evaluation&,
               const Evaluation&,
               const Evaluation&) const
{
    throw std::runtime_error("Requested the enthalpy of gas but the thermal option is not enabled");
}

template<class Scalar>
template<class Evaluation>
Evaluation DryGasPvt<Scalar>::
viscosity(unsigned regionIdx,
          const Evaluation& temperature,
          const Evaluation& pressure,
          const Evaluation& /*Rv*/,
          const Evaluation& /*Rvw*/) const
{
    return saturatedViscosity(regionIdx, temperature, pressure);
}

template<class Scalar>
template<class Evaluation>
Evaluation DryGasPvt<Scalar>::
saturatedViscosity(unsigned regionIdx,
                   const Evaluation& /*temperature*/,
                   const Evaluation& pressure) const
{
    const Evaluation& invBg = inverseGasB_[regionIdx].eval(pressure, /*extrapolate=*/true);
    const Evaluation& invMugBg = inverseGasBMu_[regionIdx].eval(pressure, /*extrapolate=*/true);

    return invBg/invMugBg;
}

template<class Scalar>
template<class Evaluation>
Evaluation DryGasPvt<Scalar>::
inverseFormationVolumeFactor(unsigned regionIdx,
                             const Evaluation& temperature,
                             const Evaluation& pressure,
                             const Evaluation& /*Rv*/,
                             const Evaluation& /*Rvw*/) const
{
    return saturatedInverseFormationVolumeFactor(regionIdx, temperature, pressure);
}

template<class Scalar>
template<class Evaluation>
Evaluation DryGasPvt<Scalar>::
saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                      const Evaluation& /*temperature*/,
                                      const Evaluation& pressure) const
{
    return inverseGasB_[regionIdx].eval(pressure, /*extrapolate=*/true);
}

template<class Scalar>
template<class Evaluation>
Evaluation DryGasPvt<Scalar>::
saturationPressure(unsigned /*regionIdx*/,
                   const Evaluation& /*temperature*/,
                   const Evaluation& /*Rv*/) const
{
    return 0.0; /* this is dry gas! */
}

template<class Scalar>
template<class Evaluation>
Evaluation DryGasPvt<Scalar>::
saturatedWaterVaporizationFactor(unsigned /*regionIdx*/,
                                 const Evaluation& /*temperature*/,
                                 const Evaluation& /*pressure*/) const
{
    return 0.0; /* this is non-humid gas! */
}

template<class Scalar>
template<class Evaluation>
Evaluation DryGasPvt<Scalar>::
saturatedOilVaporizationFactor(unsigned /*regionIdx*/,
                               const Evaluation& /*temperature*/,
                               const Evaluation& /*pressure*/,
                               const Evaluation& /*oilSaturation*/,
                               const Evaluation& /*maxOilSaturation*/) const
{
    return 0.0; /* this is dry gas! */
}

template<class Scalar>
template<class Evaluation>
Evaluation DryGasPvt<Scalar>::
saturatedOilVaporizationFactor(unsigned /*regionIdx*/,
                               const Evaluation& /*temperature*/,
                               const Evaluation& /*pressure*/) const
{
    return 0.0; /* this is dry gas! */
}

template<class Scalar>
template<class Evaluation>
Evaluation DryGasPvt<Scalar>::
diffusionCoefficient(const Evaluation& /*temperature*/,
                     const Evaluation& /*pressure*/,
                     unsigned /*compIdx*/) const
{
    throw std::runtime_error("Not implemented: The PVT model does not provide a diffusionCoefficient()");
}

} // namespace Opm
