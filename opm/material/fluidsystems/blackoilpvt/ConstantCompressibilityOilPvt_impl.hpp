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
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#endif

namespace Opm {

template<class Scalar>
ConstantCompressibilityOilPvt<Scalar>::
ConstantCompressibilityOilPvt(const std::vector<Scalar>& oilReferenceDensity,
                              const std::vector<Scalar>& oilReferencePressure,
                              const std::vector<Scalar>& oilReferenceFormationVolumeFactor,
                              const std::vector<Scalar>& oilCompressibility,
                              const std::vector<Scalar>& oilViscosity,
                              const std::vector<Scalar>& oilViscosibility)
    : oilReferenceDensity_(oilReferenceDensity)
    , oilReferencePressure_(oilReferencePressure)
    , oilReferenceFormationVolumeFactor_(oilReferenceFormationVolumeFactor)
    , oilCompressibility_(oilCompressibility)
    , oilViscosity_(oilViscosity)
    , oilViscosibility_(oilViscosibility)
{
}

#if HAVE_ECL_INPUT
template<class Scalar>
void ConstantCompressibilityOilPvt<Scalar>::
initFromState(const EclipseState& eclState, const Schedule&)
{
    const auto& pvcdoTable = eclState.getTableManager().getPvcdoTable();
    const auto& densityTable = eclState.getTableManager().getDensityTable();

    assert(pvcdoTable.size() == densityTable.size());

    size_t numRegions = pvcdoTable.size();
    setNumRegions(numRegions);

    for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
        Scalar rhoRefO = densityTable[regionIdx].oil;
        Scalar rhoRefG = densityTable[regionIdx].gas;
        Scalar rhoRefW = densityTable[regionIdx].water;

        setReferenceDensities(regionIdx, rhoRefO, rhoRefG, rhoRefW);

        oilReferencePressure_[regionIdx] = pvcdoTable[regionIdx].reference_pressure;
        oilReferenceFormationVolumeFactor_[regionIdx] = pvcdoTable[regionIdx].volume_factor;
        oilCompressibility_[regionIdx] = pvcdoTable[regionIdx].compressibility;
        oilViscosity_[regionIdx] = pvcdoTable[regionIdx].viscosity;
        oilViscosibility_[regionIdx] = pvcdoTable[regionIdx].viscosibility;
    }

    initEnd();
}
#endif

template<class Scalar>
void ConstantCompressibilityOilPvt<Scalar>::setNumRegions(size_t numRegions)
{
    oilReferenceDensity_.resize(numRegions);
    oilReferencePressure_.resize(numRegions);
    oilReferenceFormationVolumeFactor_.resize(numRegions);
    oilCompressibility_.resize(numRegions);
    oilViscosity_.resize(numRegions);
    oilViscosibility_.resize(numRegions);

    for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
        setReferenceFormationVolumeFactor(regionIdx, 1.0);
        setReferencePressure(regionIdx, 1.03125);
    }
}

template<class Scalar>
void ConstantCompressibilityOilPvt<Scalar>::
setReferenceDensities(unsigned regionIdx,
                      Scalar rhoRefOil,
                      Scalar /*rhoRefGas*/,
                      Scalar /*rhoRefWater*/)
{
    oilReferenceDensity_[regionIdx] = rhoRefOil;
}

template<class Scalar>
void ConstantCompressibilityOilPvt<Scalar>::
setViscosity(unsigned regionIdx, Scalar muo, Scalar oilViscosibility)
{
    oilViscosity_[regionIdx] = muo;
    oilViscosibility_[regionIdx] = oilViscosibility;
}

template<class Scalar>
void ConstantCompressibilityOilPvt<Scalar>::
setCompressibility(unsigned regionIdx, Scalar oilCompressibility)
{
    oilCompressibility_[regionIdx] = oilCompressibility;
}

template<class Scalar>
void ConstantCompressibilityOilPvt<Scalar>::
setReferencePressure(unsigned regionIdx, Scalar p)
{
    oilReferencePressure_[regionIdx] = p;
}

template<class Scalar>
void ConstantCompressibilityOilPvt<Scalar>::
setReferenceFormationVolumeFactor(unsigned regionIdx, Scalar BoRef)
{
    oilReferenceFormationVolumeFactor_[regionIdx] = BoRef;
}

template<class Scalar>
void ConstantCompressibilityOilPvt<Scalar>::
setViscosibility(unsigned regionIdx, Scalar muComp)
{
    oilViscosibility_[regionIdx] = muComp;
}

template<class Scalar>
bool ConstantCompressibilityOilPvt<Scalar>::
operator==(const ConstantCompressibilityOilPvt<Scalar>& data) const
{
    return this->oilReferenceDensity_ == data.oilReferenceDensity_ &&
           this->oilReferencePressure_ == data.oilReferencePressure_ &&
           this->oilReferenceFormationVolumeFactor() == data.oilReferenceFormationVolumeFactor() &&
           this->oilCompressibility() == data.oilCompressibility() &&
           this->oilViscosity() == data.oilViscosity() &&
           this->oilViscosibility() == data.oilViscosibility();
}

template<class Scalar>
template<class Evaluation>
Evaluation ConstantCompressibilityOilPvt<Scalar>::
internalEnergy(unsigned,
               const Evaluation&,
               const Evaluation&,
               const Evaluation&) const
{
    throw std::runtime_error("Requested the enthalpy of oil but the thermal option is not enabled");
}

template<class Scalar>
template<class Evaluation>
Evaluation ConstantCompressibilityOilPvt<Scalar>::
viscosity(unsigned regionIdx,
          const Evaluation& temperature,
          const Evaluation& pressure,
          const Evaluation& /*Rs*/) const
{
    return saturatedViscosity(regionIdx, temperature, pressure);
}

template<class Scalar>
template<class Evaluation>
Evaluation ConstantCompressibilityOilPvt<Scalar>::
saturatedViscosity(unsigned regionIdx,
                   const Evaluation& temperature,
                   const Evaluation& pressure) const
{
    Scalar BoMuoRef = oilViscosity_[regionIdx]*oilReferenceFormationVolumeFactor_[regionIdx];
    const Evaluation& bo = saturatedInverseFormationVolumeFactor(regionIdx, temperature, pressure);

    Scalar pRef = oilReferencePressure_[regionIdx];
    const Evaluation& Y =
        (oilCompressibility_[regionIdx] - oilViscosibility_[regionIdx])
        * (pressure - pRef);
    return BoMuoRef*bo/(1.0 + Y*(1.0 + Y/2.0));
}

template<class Scalar>
template<class Evaluation>
Evaluation ConstantCompressibilityOilPvt<Scalar>::
inverseFormationVolumeFactor(unsigned regionIdx,
                             const Evaluation& temperature,
                             const Evaluation& pressure,
                             const Evaluation& /*Rs*/) const
{
    return saturatedInverseFormationVolumeFactor(regionIdx, temperature, pressure);
}

template<class Scalar>
template<class Evaluation>
Evaluation ConstantCompressibilityOilPvt<Scalar>::
saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                      const Evaluation& /*temperature*/,
                                      const Evaluation& pressure) const
{
    // cf. ECLiPSE 2011 technical description, p. 116
    Scalar pRef = oilReferencePressure_[regionIdx];
    const Evaluation& X = oilCompressibility_[regionIdx]*(pressure - pRef);

    Scalar BoRef = oilReferenceFormationVolumeFactor_[regionIdx];
    return (1 + X*(1 + X/2))/BoRef;
}

template<class Scalar>
template<class Evaluation>
Evaluation ConstantCompressibilityOilPvt<Scalar>::
saturatedGasDissolutionFactor(unsigned /*regionIdx*/,
                              const Evaluation& /*temperature*/,
                              const Evaluation& /*pressure*/) const
{
    return 0.0; /* this is dead oil! */
}

template<class Scalar>
template<class Evaluation>
Evaluation ConstantCompressibilityOilPvt<Scalar>::
saturationPressure(unsigned /*regionIdx*/,
                   const Evaluation& /*temperature*/,
                   const Evaluation& /*Rs*/) const
{
    return 0.0; /* this is dead oil, so there isn't any meaningful saturation pressure! */
}

template<class Scalar>
template<class Evaluation>
Evaluation ConstantCompressibilityOilPvt<Scalar>::
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
Evaluation ConstantCompressibilityOilPvt<Scalar>::
diffusionCoefficient(const Evaluation& /*temperature*/,
                     const Evaluation& /*pressure*/,
                     unsigned /*compIdx*/) const
{
    throw std::runtime_error("Not implemented: The PVT model does not provide a diffusionCoefficient()");
}

} // namespace Opm
