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

#define OPM_GAS_PVT_MULTIPLEXER_CALL(codeToCall)                          \
    switch (gasPvtApproach_) {                                            \
    case GasPvtApproach::DryGasPvt: {                                     \
        auto& pvtImpl = getRealPvt<GasPvtApproach::DryGasPvt>();          \
        codeToCall;                                                       \
        break;                                                            \
    }                                                                     \
    case GasPvtApproach::DryHumidGasPvt: {                                \
        auto& pvtImpl = getRealPvt<GasPvtApproach::DryHumidGasPvt>();     \
        codeToCall;                                                       \
        break;                                                            \
    }                                                                     \
    case GasPvtApproach::WetHumidGasPvt: {                                \
        auto& pvtImpl = getRealPvt<GasPvtApproach::WetHumidGasPvt>();     \
        codeToCall;                                                       \
        break;                                                            \
    }                                                                     \
    case GasPvtApproach::WetGasPvt: {                                     \
        auto& pvtImpl = getRealPvt<GasPvtApproach::WetGasPvt>();          \
        codeToCall;                                                       \
        break;                                                            \
    }                                                                     \
    case GasPvtApproach::ThermalGasPvt: {                                 \
        auto& pvtImpl = getRealPvt<GasPvtApproach::ThermalGasPvt>();      \
        codeToCall;                                                       \
        break;                                                            \
    }                                                                     \
    case GasPvtApproach::Co2GasPvt: {                                     \
        auto& pvtImpl = getRealPvt<GasPvtApproach::Co2GasPvt>();          \
        codeToCall;                                                       \
        break;                                                            \
    }                                                                     \
    case GasPvtApproach::NoGasPvt:                                        \
        throw std::logic_error("Not implemented: Gas PVT of this deck!"); \
    } \

namespace Opm {

template<class Scalar, bool enableThermal>
GasPvtMultiplexer<Scalar,enableThermal>::GasPvtMultiplexer()
{
    gasPvtApproach_ = GasPvtApproach::NoGasPvt;
    realGasPvt_ = nullptr;
}

template<class Scalar, bool enableThermal>
GasPvtMultiplexer<Scalar,enableThermal>::
GasPvtMultiplexer(GasPvtApproach approach, void* realGasPvt)
    : gasPvtApproach_(approach)
    , realGasPvt_(realGasPvt)
{
}

template<class Scalar, bool enableThermal>
GasPvtMultiplexer<Scalar,enableThermal>::
GasPvtMultiplexer(const GasPvtMultiplexer<Scalar,enableThermal>& data)
{
    *this = data;
}

template<class Scalar, bool enableThermal>
GasPvtMultiplexer<Scalar,enableThermal>::~GasPvtMultiplexer()
{
    switch (gasPvtApproach_) {
    case GasPvtApproach::DryGasPvt: {
        delete &getRealPvt<GasPvtApproach::DryGasPvt>();
        break;
    }
    case GasPvtApproach::DryHumidGasPvt: {
        delete &getRealPvt<GasPvtApproach::DryHumidGasPvt>();
        break;
    }
    case GasPvtApproach::WetHumidGasPvt: {
        delete &getRealPvt<GasPvtApproach::WetHumidGasPvt>();
        break;
    }
    case GasPvtApproach::WetGasPvt: {
        delete &getRealPvt<GasPvtApproach::WetGasPvt>();
        break;
    }
    case GasPvtApproach::ThermalGasPvt: {
        delete &getRealPvt<GasPvtApproach::ThermalGasPvt>();
        break;
    }
    case GasPvtApproach::Co2GasPvt: {
        delete &getRealPvt<GasPvtApproach::Co2GasPvt>();
        break;
    }
    case GasPvtApproach::NoGasPvt:
        break;
    }
}

#if HAVE_ECL_INPUT
template<class Scalar, bool enableThermal>
void GasPvtMultiplexer<Scalar,enableThermal>::
initFromState(const EclipseState& eclState, const Schedule& schedule)
{
    if (!eclState.runspec().phases().active(Phase::GAS))
        return;
    if (eclState.runspec().co2Storage())
        setApproach(GasPvtApproach::Co2GasPvt);
    else if (enableThermal && eclState.getSimulationConfig().isThermal())
        setApproach(GasPvtApproach::ThermalGasPvt);
    else if (!eclState.getTableManager().getPvtgwTables().empty() && !eclState.getTableManager().getPvtgTables().empty())
        setApproach(GasPvtApproach::WetHumidGasPvt);
    else if (!eclState.getTableManager().getPvtgTables().empty())
        setApproach(GasPvtApproach::WetGasPvt);
    else if (eclState.getTableManager().hasTables("PVDG"))
        setApproach(GasPvtApproach::DryGasPvt);
    else if (!eclState.getTableManager().getPvtgwTables().empty())
        setApproach(GasPvtApproach::DryHumidGasPvt);


    OPM_GAS_PVT_MULTIPLEXER_CALL(pvtImpl.initFromState(eclState, schedule));
}
#endif // HAVE_ECL_INPUT

template<class Scalar, bool enableThermal>
void GasPvtMultiplexer<Scalar,enableThermal>::setApproach(GasPvtApproach gasPvtAppr)
{
    switch (gasPvtAppr) {
    case GasPvtApproach::DryGasPvt:
        realGasPvt_ = new DryGasPvt<Scalar>;
        break;

    case GasPvtApproach::DryHumidGasPvt:
        realGasPvt_ = new DryHumidGasPvt<Scalar>;
        break;

    case GasPvtApproach::WetHumidGasPvt:
        realGasPvt_ = new WetHumidGasPvt<Scalar>;
        break;

    case GasPvtApproach::WetGasPvt:
        realGasPvt_ = new WetGasPvt<Scalar>;
        break;

    case GasPvtApproach::ThermalGasPvt:
        realGasPvt_ = new GasPvtThermal<Scalar>;
        break;

    case GasPvtApproach::Co2GasPvt:
        realGasPvt_ = new Co2GasPvt<Scalar>;
        break;

    case GasPvtApproach::NoGasPvt:
        throw std::logic_error("Not implemented: Gas PVT of this deck!");
    }

    gasPvtApproach_ = gasPvtAppr;
}

template<class Scalar, bool enableThermal>
void GasPvtMultiplexer<Scalar,enableThermal>::initEnd()
{
    OPM_GAS_PVT_MULTIPLEXER_CALL(pvtImpl.initEnd());
}

template<class Scalar, bool enableThermal>
unsigned GasPvtMultiplexer<Scalar,enableThermal>::numRegions() const
{
    OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.numRegions());
    return 1;
}

template<class Scalar, bool enableThermal>
const Scalar GasPvtMultiplexer<Scalar,enableThermal>::
gasReferenceDensity(unsigned regionIdx)
{
    OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.gasReferenceDensity(regionIdx));
    return 2.;
}

template<class Scalar, bool enableThermal>
bool GasPvtMultiplexer<Scalar,enableThermal>::
operator==(const GasPvtMultiplexer<Scalar,enableThermal>& data) const
{
    if (this->gasPvtApproach() != data.gasPvtApproach())
        return false;

    switch (gasPvtApproach_) {
    case GasPvtApproach::DryGasPvt:
        return *static_cast<const DryGasPvt<Scalar>*>(realGasPvt_) ==
               *static_cast<const DryGasPvt<Scalar>*>(data.realGasPvt_);
    case GasPvtApproach::DryHumidGasPvt:
        return *static_cast<const DryHumidGasPvt<Scalar>*>(realGasPvt_) ==
               *static_cast<const DryHumidGasPvt<Scalar>*>(data.realGasPvt_);
    case GasPvtApproach::WetHumidGasPvt:
        return *static_cast<const WetHumidGasPvt<Scalar>*>(realGasPvt_) ==
               *static_cast<const WetHumidGasPvt<Scalar>*>(data.realGasPvt_);
    case GasPvtApproach::WetGasPvt:
        return *static_cast<const WetGasPvt<Scalar>*>(realGasPvt_) ==
               *static_cast<const WetGasPvt<Scalar>*>(data.realGasPvt_);
    case GasPvtApproach::ThermalGasPvt:
        return *static_cast<const GasPvtThermal<Scalar>*>(realGasPvt_) ==
               *static_cast<const GasPvtThermal<Scalar>*>(data.realGasPvt_);
    case GasPvtApproach::Co2GasPvt:
        return *static_cast<const Co2GasPvt<Scalar>*>(realGasPvt_) ==
                *static_cast<const Co2GasPvt<Scalar>*>(data.realGasPvt_);
    default:
        return true;
    }
}

template<class Scalar, bool enableThermal>
GasPvtMultiplexer<Scalar,enableThermal>&
GasPvtMultiplexer<Scalar,enableThermal>::
operator=(const GasPvtMultiplexer<Scalar,enableThermal>& data)
{
    gasPvtApproach_ = data.gasPvtApproach_;
    switch (gasPvtApproach_) {
    case GasPvtApproach::DryGasPvt:
        realGasPvt_ = new DryGasPvt<Scalar>(*static_cast<const DryGasPvt<Scalar>*>(data.realGasPvt_));
        break;
    case GasPvtApproach::DryHumidGasPvt:
        realGasPvt_ = new DryHumidGasPvt<Scalar>(*static_cast<const DryHumidGasPvt<Scalar>*>(data.realGasPvt_));
        break;
    case GasPvtApproach::WetHumidGasPvt:
        realGasPvt_ = new WetHumidGasPvt<Scalar>(*static_cast<const WetHumidGasPvt<Scalar>*>(data.realGasPvt_));
        break;
    case GasPvtApproach::WetGasPvt:
        realGasPvt_ = new WetGasPvt<Scalar>(*static_cast<const WetGasPvt<Scalar>*>(data.realGasPvt_));
        break;
    case GasPvtApproach::ThermalGasPvt:
        realGasPvt_ = new GasPvtThermal<Scalar>(*static_cast<const GasPvtThermal<Scalar>*>(data.realGasPvt_));
        break;
    case GasPvtApproach::Co2GasPvt:
        realGasPvt_ = new Co2GasPvt<Scalar>(*static_cast<const Co2GasPvt<Scalar>*>(data.realGasPvt_));
        break;
    default:
        break;
    }

    return *this;
}

template<class Scalar, bool enableThermal>
template<class Evaluation>
Evaluation GasPvtMultiplexer<Scalar,enableThermal>::
internalEnergy(unsigned regionIdx,
               const Evaluation& temperature,
               const Evaluation& pressure,
               const Evaluation& Rv) const
{
    OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.internalEnergy(regionIdx, temperature, pressure, Rv));
    return 0;
}

template<class Scalar, bool enableThermal>
template<class Evaluation>
Evaluation GasPvtMultiplexer<Scalar,enableThermal>::
viscosity(unsigned regionIdx,
          const Evaluation& temperature,
          const Evaluation& pressure,
          const Evaluation& Rv,
          const Evaluation& Rvw ) const
{
    OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.viscosity(regionIdx, temperature, pressure, Rv, Rvw));
    return 0;
}

template<class Scalar, bool enableThermal>
template<class Evaluation>
Evaluation GasPvtMultiplexer<Scalar,enableThermal>::
saturatedViscosity(unsigned regionIdx,
                   const Evaluation& temperature,
                   const Evaluation& pressure) const
{
    OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedViscosity(regionIdx, temperature, pressure));
    return 0;
}

template<class Scalar, bool enableThermal>
template<class Evaluation>
Evaluation GasPvtMultiplexer<Scalar,enableThermal>::
inverseFormationVolumeFactor(unsigned regionIdx,
                             const Evaluation& temperature,
                             const Evaluation& pressure,
                             const Evaluation& Rv,
                             const Evaluation& Rvw) const
{
    OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.inverseFormationVolumeFactor(regionIdx, temperature, pressure, Rv, Rvw));
    return 0;
}

template<class Scalar, bool enableThermal>
template<class Evaluation>
Evaluation GasPvtMultiplexer<Scalar,enableThermal>::
saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                      const Evaluation& temperature,
                                      const Evaluation& pressure) const
{
    OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedInverseFormationVolumeFactor(regionIdx, temperature, pressure));
    return 0;
}

template<class Scalar, bool enableThermal>
template<class Evaluation>
Evaluation GasPvtMultiplexer<Scalar,enableThermal>::
saturatedOilVaporizationFactor(unsigned regionIdx,
                               const Evaluation& temperature,
                               const Evaluation& pressure) const
{
    OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedOilVaporizationFactor(regionIdx, temperature, pressure));
    return 0;
}

template<class Scalar, bool enableThermal>
template<class Evaluation>
Evaluation GasPvtMultiplexer<Scalar,enableThermal>::
saturatedOilVaporizationFactor(unsigned regionIdx,
                               const Evaluation& temperature,
                               const Evaluation& pressure,
                               const Evaluation& oilSaturation,
                               const Evaluation& maxOilSaturation) const
{
    OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedOilVaporizationFactor(regionIdx, temperature, pressure, oilSaturation, maxOilSaturation));
    return 0;
}

template<class Scalar, bool enableThermal>
template<class Evaluation>
Evaluation GasPvtMultiplexer<Scalar,enableThermal>::
saturatedWaterVaporizationFactor(unsigned regionIdx,
                                 const Evaluation& temperature,
                                 const Evaluation& pressure) const
{
    OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedWaterVaporizationFactor(regionIdx, temperature, pressure));
    return 0;
}

template<class Scalar, bool enableThermal>
template<class Evaluation>
Evaluation GasPvtMultiplexer<Scalar,enableThermal>::
saturationPressure(unsigned regionIdx,
                   const Evaluation& temperature,
                   const Evaluation& Rv) const
{
    OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.saturationPressure(regionIdx, temperature, Rv));
    return 0;
}

template<class Scalar, bool enableThermal>
template<class Evaluation>
Evaluation GasPvtMultiplexer<Scalar,enableThermal>::
diffusionCoefficient(const Evaluation& temperature,
                     const Evaluation& pressure,
                     unsigned compIdx) const
{
  OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.diffusionCoefficient(temperature, pressure, compIdx));
  return 0;
}

} // namespace Opm

#undef OPM_GAS_PVT_MULTIPLEXER_CALL
