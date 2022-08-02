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

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PvdoTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#endif

namespace Opm {

template<class Scalar>
DeadOilPvt<Scalar>::DeadOilPvt(const std::vector<Scalar>& oilReferenceDensity,
                               const std::vector<TabulatedOneDFunction>& inverseOilB,
                               const std::vector<TabulatedOneDFunction>& oilMu,
                               const std::vector<TabulatedOneDFunction>& inverseOilBMu)
    : oilReferenceDensity_(oilReferenceDensity)
    , inverseOilB_(inverseOilB)
    , oilMu_(oilMu)
    , inverseOilBMu_(inverseOilBMu)
{
}

#if HAVE_ECL_INPUT
template<class Scalar>
void DeadOilPvt<Scalar>::
initFromState(const EclipseState& eclState, const Schedule&)
{
    const auto& pvdoTables = eclState.getTableManager().getPvdoTables();
    const auto& densityTable = eclState.getTableManager().getDensityTable();

    assert(pvdoTables.size() == densityTable.size());

    size_t numRegions = pvdoTables.size();
    setNumRegions(numRegions);

    for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
        Scalar rhoRefO = densityTable[regionIdx].oil;
        Scalar rhoRefG = densityTable[regionIdx].gas;
        Scalar rhoRefW = densityTable[regionIdx].water;

        setReferenceDensities(regionIdx, rhoRefO, rhoRefG, rhoRefW);

        const auto& pvdoTable = pvdoTables.getTable<PvdoTable>(regionIdx);

        const auto& BColumn(pvdoTable.getFormationFactorColumn());
        std::vector<Scalar> invBColumn(BColumn.size());
        for (unsigned i = 0; i < invBColumn.size(); ++i)
            invBColumn[i] = 1/BColumn[i];

        inverseOilB_[regionIdx].setXYArrays(pvdoTable.numRows(),
                                            pvdoTable.getPressureColumn(),
                                            invBColumn);
        oilMu_[regionIdx].setXYArrays(pvdoTable.numRows(),
                                      pvdoTable.getPressureColumn(),
                                      pvdoTable.getViscosityColumn());
    }

    initEnd();
}
#endif // HAVE_ECL_INPUT

template<class Scalar>
void DeadOilPvt<Scalar>::setNumRegions(size_t numRegions)
{
    oilReferenceDensity_.resize(numRegions);
    inverseOilB_.resize(numRegions);
    inverseOilBMu_.resize(numRegions);
    oilMu_.resize(numRegions);
}

template<class Scalar>
void DeadOilPvt<Scalar>::setReferenceDensities(unsigned regionIdx,
                                               Scalar rhoRefOil,
                                               Scalar /*rhoRefGas*/,
                                               Scalar /*rhoRefWater*/)
{
    oilReferenceDensity_[regionIdx] = rhoRefOil;
}

template<class Scalar>
void DeadOilPvt<Scalar>::
setInverseOilFormationVolumeFactor(unsigned regionIdx,
                                   const TabulatedOneDFunction& invBo)
{
    inverseOilB_[regionIdx] = invBo;
}

template<class Scalar>
void DeadOilPvt<Scalar>::
setOilViscosity(unsigned regionIdx,
                const TabulatedOneDFunction& muo)
{
    oilMu_[regionIdx] = muo;
}

template<class Scalar>
void DeadOilPvt<Scalar>::initEnd()
{
    // calculate the final 2D functions which are used for interpolation.
    size_t numRegions = oilMu_.size();
    for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
        // calculate the table which stores the inverse of the product of the oil
        // formation volume factor and the oil viscosity
        const auto& oilMu = oilMu_[regionIdx];
        const auto& invOilB = inverseOilB_[regionIdx];
        assert(oilMu.numSamples() == invOilB.numSamples());

        std::vector<Scalar> invBMuColumn;
        std::vector<Scalar> pressureColumn;
        invBMuColumn.resize(oilMu.numSamples());
        pressureColumn.resize(oilMu.numSamples());

        for (unsigned pIdx = 0; pIdx < oilMu.numSamples(); ++pIdx) {
            pressureColumn[pIdx] = invOilB.xAt(pIdx);
            invBMuColumn[pIdx] = invOilB.valueAt(pIdx)*1/oilMu.valueAt(pIdx);
        }

        inverseOilBMu_[regionIdx].setXYArrays(pressureColumn.size(),
                                              pressureColumn,
                                              invBMuColumn);
    }
}

template<class Scalar>
bool DeadOilPvt<Scalar>::operator==(const DeadOilPvt<Scalar>& data) const
{
    return this->oilReferenceDensity_ == data.oilReferenceDensity_ &&
           this->inverseOilB() == data.inverseOilB() &&
           this->oilMu() == data.oilMu() &&
           this->inverseOilBMu() == data.inverseOilBMu();
}

template<class Scalar>
template<class Evaluation>
Evaluation DeadOilPvt<Scalar>::
internalEnergy(unsigned,
               const Evaluation&,
               const Evaluation&,
               const Evaluation&) const
{
    throw std::runtime_error("Requested the enthalpy of oil but the thermal option is not enabled");
}

template<class Scalar>
template<class Evaluation>
Evaluation DeadOilPvt<Scalar>::
viscosity(unsigned regionIdx,
          const Evaluation& temperature,
          const Evaluation& pressure,
          const Evaluation& /*Rs*/) const
{
    return saturatedViscosity(regionIdx, temperature, pressure);
}

template<class Scalar>
template<class Evaluation>
Evaluation DeadOilPvt<Scalar>::
saturatedViscosity(unsigned regionIdx,
                   const Evaluation& /*temperature*/,
                   const Evaluation& pressure) const
{
    const Evaluation& invBo = inverseOilB_[regionIdx].eval(pressure, /*extrapolate=*/true);
    const Evaluation& invMuoBo = inverseOilBMu_[regionIdx].eval(pressure, /*extrapolate=*/true);

    return invBo/invMuoBo;
}

template<class Scalar>
template<class Evaluation>
Evaluation DeadOilPvt<Scalar>::
inverseFormationVolumeFactor(unsigned regionIdx,
                             const Evaluation& /*temperature*/,
                             const Evaluation& pressure,
                             const Evaluation& /*Rs*/) const
{
    return inverseOilB_[regionIdx].eval(pressure, /*extrapolate=*/true);
}

template<class Scalar>
template<class Evaluation>
Evaluation DeadOilPvt<Scalar>::
saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                      const Evaluation& /*temperature*/,
                                      const Evaluation& pressure) const
{
    return inverseOilB_[regionIdx].eval(pressure, /*extrapolate=*/true);
}

template<class Scalar>
template<class Evaluation>
Evaluation DeadOilPvt<Scalar>::
saturatedGasDissolutionFactor(unsigned /*regionIdx*/,
                              const Evaluation& /*temperature*/,
                              const Evaluation& /*pressure*/) const
{
    return 0.0; /* this is dead oil! */
}

template<class Scalar>
template<class Evaluation>
Evaluation DeadOilPvt<Scalar>::
saturatedGasDissolutionFactor(unsigned /*regionIdx*/,
                              const Evaluation& /*temperature*/,
                              const Evaluation& /*pressure*/,
                              const Evaluation& /*oilSaturation*/,
                              const Evaluation& /*maxOilSaturation*/) const
{
    return 0.0; /* this is dead oil! */
}

template<class Scalar>
template<class Evaluation>
Evaluation DeadOilPvt<Scalar>::
saturationPressure(unsigned /*regionIdx*/,
                   const Evaluation& /*temperature*/,
                   const Evaluation& /*Rs*/) const
{
    return 0.0; /* this is dead oil, so there isn't any meaningful saturation pressure! */
}

template<class Scalar>
template<class Evaluation>
Evaluation DeadOilPvt<Scalar>::
diffusionCoefficient(const Evaluation& /*temperature*/,
                     const Evaluation& /*pressure*/,
                     unsigned /*compIdx*/) const
{
    throw std::runtime_error("Not implemented: The PVT model does not provide a diffusionCoefficient()");
}

} // namespace Opm
