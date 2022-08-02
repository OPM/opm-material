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
#endif

namespace Opm {

template<class Scalar>
ConstantCompressibilityWaterPvt<Scalar>::
ConstantCompressibilityWaterPvt(const std::vector<Scalar>& waterReferenceDensity,
                                const std::vector<Scalar>& waterReferencePressure,
                                const std::vector<Scalar>& waterReferenceFormationVolumeFactor,
                                const std::vector<Scalar>& waterCompressibility,
                                const std::vector<Scalar>& waterViscosity,
                                const std::vector<Scalar>& waterViscosibility)
    : waterReferenceDensity_(waterReferenceDensity)
    , waterReferencePressure_(waterReferencePressure)
    , waterReferenceFormationVolumeFactor_(waterReferenceFormationVolumeFactor)
    , waterCompressibility_(waterCompressibility)
    , waterViscosity_(waterViscosity)
    , waterViscosibility_(waterViscosibility)
{
}

#if HAVE_ECL_INPUT
template<class Scalar>
void ConstantCompressibilityWaterPvt<Scalar>::
initFromState(const EclipseState& eclState, const Schedule&)
{
    const auto& pvtwTable = eclState.getTableManager().getPvtwTable();
    const auto& densityTable = eclState.getTableManager().getDensityTable();

    assert(pvtwTable.size() == densityTable.size());

    size_t numRegions = pvtwTable.size();
    setNumRegions(numRegions);

    for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
        waterReferenceDensity_[regionIdx] = densityTable[regionIdx].water;

        waterReferencePressure_[regionIdx] = pvtwTable[regionIdx].reference_pressure;
        waterReferenceFormationVolumeFactor_[regionIdx] = pvtwTable[regionIdx].volume_factor;
        waterCompressibility_[regionIdx] = pvtwTable[regionIdx].compressibility;
        waterViscosity_[regionIdx] = pvtwTable[regionIdx].viscosity;
        waterViscosibility_[regionIdx] = pvtwTable[regionIdx].viscosibility;
    }

    initEnd();
}
#endif

template<class Scalar>
void ConstantCompressibilityWaterPvt<Scalar>::setNumRegions(size_t numRegions)
{
    waterReferenceDensity_.resize(numRegions);
    waterReferencePressure_.resize(numRegions);
    waterReferenceFormationVolumeFactor_.resize(numRegions);
    waterCompressibility_.resize(numRegions);
    waterViscosity_.resize(numRegions);
    waterViscosibility_.resize(numRegions);

    for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
        setReferenceDensities(regionIdx, 650.0, 1.0, 1000.0);
        setReferenceFormationVolumeFactor(regionIdx, 1.0);
        setReferencePressure(regionIdx, 1e5);
    }
}

template<class Scalar>
void ConstantCompressibilityWaterPvt<Scalar>::
setReferenceDensities(unsigned regionIdx,
                      Scalar /*rhoRefOil*/,
                      Scalar /*rhoRefGas*/,
                      Scalar rhoRefWater)
{
    waterReferenceDensity_[regionIdx] = rhoRefWater;
}

template<class Scalar>
void ConstantCompressibilityWaterPvt<Scalar>::
setReferencePressure(unsigned regionIdx, Scalar p)
{
    waterReferencePressure_[regionIdx] = p;
}

template<class Scalar>
void ConstantCompressibilityWaterPvt<Scalar>::
setViscosity(unsigned regionIdx, Scalar muw, Scalar waterViscosibility)
{
    waterViscosity_[regionIdx] = muw;
    waterViscosibility_[regionIdx] = waterViscosibility;
}

template<class Scalar>
void ConstantCompressibilityWaterPvt<Scalar>::
setCompressibility(unsigned regionIdx, Scalar waterCompressibility)
{
    waterCompressibility_[regionIdx] = waterCompressibility;
}

template<class Scalar>
void ConstantCompressibilityWaterPvt<Scalar>::
setReferenceFormationVolumeFactor(unsigned regionIdx, Scalar BwRef)
{
    waterReferenceFormationVolumeFactor_[regionIdx] = BwRef;
}

template<class Scalar>
void ConstantCompressibilityWaterPvt<Scalar>::
setViscosibility(unsigned regionIdx, Scalar muComp)
{
    waterViscosibility_[regionIdx] = muComp;
}

template<class Scalar>
bool ConstantCompressibilityWaterPvt<Scalar>::
operator==(const ConstantCompressibilityWaterPvt<Scalar>& data) const
{
    return this->waterReferenceDensity_ == data.waterReferenceDensity_ &&
           this->waterReferencePressure() == data.waterReferencePressure() &&
           this->waterReferenceFormationVolumeFactor() == data.waterReferenceFormationVolumeFactor() &&
           this->waterCompressibility() == data.waterCompressibility() &&
           this->waterViscosity() == data.waterViscosity() &&
           this->waterViscosibility() == data.waterViscosibility();
}

template<class Scalar>
template<class Evaluation>
Evaluation ConstantCompressibilityWaterPvt<Scalar>::
internalEnergy(unsigned,
               const Evaluation&,
               const Evaluation&,
               const Evaluation&) const
{
  throw std::runtime_error("Requested the enthalpy of water but the thermal option is not enabled");
}

template<class Scalar>
template<class Evaluation>
Evaluation ConstantCompressibilityWaterPvt<Scalar>::
viscosity(unsigned regionIdx,
          const Evaluation& temperature,
          const Evaluation& pressure,
          const Evaluation& saltconcentration) const
{
    Scalar BwMuwRef = waterViscosity_[regionIdx]*waterReferenceFormationVolumeFactor_[regionIdx];
    const Evaluation& bw = inverseFormationVolumeFactor(regionIdx, temperature, pressure, saltconcentration);

    Scalar pRef = waterReferencePressure_[regionIdx];
    const Evaluation& Y =
        (waterCompressibility_[regionIdx] - waterViscosibility_[regionIdx])
        * (pressure - pRef);
    return BwMuwRef*bw/(1 + Y*(1 + Y/2));
}

template<class Scalar>
template<class Evaluation>
Evaluation ConstantCompressibilityWaterPvt<Scalar>::
inverseFormationVolumeFactor(unsigned regionIdx,
                             const Evaluation& /*temperature*/,
                             const Evaluation& pressure,
                             const Evaluation& /*saltconcentration*/) const
{
    // cf. ECLiPSE 2011 technical description, p. 116
    Scalar pRef = waterReferencePressure_[regionIdx];
    const Evaluation& X = waterCompressibility_[regionIdx]*(pressure - pRef);

    Scalar BwRef = waterReferenceFormationVolumeFactor_[regionIdx];

    // TODO (?): consider the salt concentration of the brine
    return (1.0 + X*(1.0 + X/2.0))/BwRef;
}

} // namespace Opm
