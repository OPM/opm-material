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

#include <opm/material/binarycoefficients/Brine_CO2.hpp>
#include <opm/material/components/CO2.hpp>
#include <opm/material/components/SimpleHuDuanH2O.hpp>
#include <opm/material/common/UniformTabulated2DFunction.hpp>
#include <opm/material/components/co2tables.inc>

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#endif

namespace Opm {

template<class Scalar> using Co2 = ::Opm::CO2<Scalar,CO2Tables>;

template<class Scalar>
Co2GasPvt<Scalar>::Co2GasPvt(const std::vector<Scalar>& gasReferenceDensity)
    : gasReferenceDensity_(gasReferenceDensity)
{
}

template<class Scalar>
Co2GasPvt<Scalar>::Co2GasPvt(size_t numRegions,
                             Scalar T_ref, Scalar P_ref)
{
    setNumRegions(numRegions);
    for (size_t i = 0; i < numRegions; ++i) {
        gasReferenceDensity_[i] = Co2<Scalar>::gasDensity(T_ref, P_ref, extrapolate);
    }
}

#if HAVE_ECL_INPUT
template<class Scalar>
void Co2GasPvt<Scalar>::initFromState(const EclipseState& eclState, const Schedule&)
{
    if( !eclState.getTableManager().getDensityTable().empty()) {
        std::cerr << "WARNING: CO2STOR is enabled but DENSITY is in the deck. \n" <<
                     "The surface density is computed based on CO2-BRINE PVT at standard conditions (STCOND) and DENSITY is ignored " << std::endl;
    }

    if( eclState.getTableManager().hasTables("PVDG") || !eclState.getTableManager().getPvtgTables().empty()) {
        std::cerr << "WARNING: CO2STOR is enabled but PVDG or PVTG is in the deck. \n" <<
                     "CO2 PVT properties are computed based on the Span-Wagner pvt model and PVDG/PVTG input is ignored. " << std::endl;
    }

    // We only supported single pvt region for the co2-brine module
    size_t numRegions = 1;
    setNumRegions(numRegions);
    size_t regionIdx = 0;
    Scalar T_ref = eclState.getTableManager().stCond().temperature;
    Scalar P_ref = eclState.getTableManager().stCond().pressure;
    gasReferenceDensity_[regionIdx] = Co2<Scalar>::gasDensity(T_ref, P_ref, extrapolate);
    initEnd();
}
#endif

template<class Scalar>
void Co2GasPvt<Scalar>::setNumRegions(size_t numRegions)
{
    gasReferenceDensity_.resize(numRegions);
}

template<class Scalar>
void Co2GasPvt<Scalar>::setReferenceDensities(unsigned regionIdx,
                                              Scalar /*rhoRefOil*/,
                                              Scalar rhoRefGas,
                                              Scalar /*rhoRefWater*/)
{
    gasReferenceDensity_[regionIdx] = rhoRefGas;
}

template<class Scalar>
bool Co2GasPvt<Scalar>::operator==(const Co2GasPvt<Scalar>& data) const
{
    return gasReferenceDensity_ == data.gasReferenceDensity_;
}

template<class Scalar>
template <class Evaluation>
Evaluation Co2GasPvt<Scalar>::
internalEnergy(unsigned,
               const Evaluation& temperature,
               const Evaluation& pressure,
               const Evaluation&) const
{
    return Co2<Scalar>::gasInternalEnergy(temperature, pressure, extrapolate);
}

template<class Scalar>
template<class Evaluation>
Evaluation Co2GasPvt<Scalar>::
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
Evaluation Co2GasPvt<Scalar>::
saturatedViscosity(unsigned /*regionIdx*/,
                   const Evaluation& temperature,
                   const Evaluation& pressure) const
{
    return Co2<Scalar>::gasViscosity(temperature, pressure, extrapolate);
}

template<class Scalar>
template<class Evaluation>
Evaluation Co2GasPvt<Scalar>::
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
Evaluation Co2GasPvt<Scalar>::
saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                      const Evaluation& temperature,
                                      const Evaluation& pressure) const
{
    return Co2<Scalar>::gasDensity(temperature, pressure, extrapolate)/gasReferenceDensity_[regionIdx];
}

template<class Scalar>
template<class Evaluation>
Evaluation Co2GasPvt<Scalar>::
saturationPressure(unsigned /*regionIdx*/,
                   const Evaluation& /*temperature*/,
                   const Evaluation& /*Rv*/) const
{
    return 0.0; /* this is dry gas! */
}

template<class Scalar>
template<class Evaluation>
Evaluation Co2GasPvt<Scalar>::
saturatedWaterVaporizationFactor(unsigned /*regionIdx*/,
                                 const Evaluation& /*temperature*/,
                                 const Evaluation& /*pressure*/) const
{
    return 0.0; /* this is non-humid gas! */
}

template<class Scalar>
template<class Evaluation>
Evaluation Co2GasPvt<Scalar>::
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
Evaluation Co2GasPvt<Scalar>::
saturatedOilVaporizationFactor(unsigned /*regionIdx*/,
                               const Evaluation& /*temperature*/,
                               const Evaluation& /*pressure*/) const
{
    return 0.0; /* this is dry gas! */
}

template<class Scalar>
template<class Evaluation>
Evaluation Co2GasPvt<Scalar>::
diffusionCoefficient(const Evaluation& temperature,
                     const Evaluation& pressure,
                     unsigned /*compIdx*/) const
{
    using H2O = SimpleHuDuanH2O<Scalar>;
    using BinaryCoeffBrineCO2 = BinaryCoeff::Brine_CO2<Scalar,H2O,Co2<Scalar>>;
    return BinaryCoeffBrineCO2::gasDiffCoeff(temperature, pressure, extrapolate);
}

} // namespace Opm
