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
#include <opm/input/eclipse/EclipseState/Runspec.hpp>
#endif

#define OPM_OIL_PVT_MULTIPLEXER_CALL(codeToCall)                                     \
    switch (approach_) {                                                             \
    case OilPvtApproach::ConstantCompressibilityOilPvt: {                            \
        auto& pvtImpl = getRealPvt<OilPvtApproach::ConstantCompressibilityOilPvt>(); \
        codeToCall;                                                                  \
        break;                                                                       \
    }                                                                                \
    case OilPvtApproach::DeadOilPvt: {                                               \
        auto& pvtImpl = getRealPvt<OilPvtApproach::DeadOilPvt>();                    \
        codeToCall;                                                                  \
        break;                                                                       \
    }                                                                                \
    case OilPvtApproach::LiveOilPvt: {                                               \
        auto& pvtImpl = getRealPvt<OilPvtApproach::LiveOilPvt>();                    \
        codeToCall;                                                                  \
        break;                                                                       \
    }                                                                                \
    case OilPvtApproach::ThermalOilPvt: {                                            \
        auto& pvtImpl = getRealPvt<OilPvtApproach::ThermalOilPvt>();                 \
        codeToCall;                                                                  \
        break;                                                                       \
    }                                                                                \
    case OilPvtApproach::BrineCo2Pvt: {                                              \
        auto& pvtImpl = getRealPvt<OilPvtApproach::BrineCo2Pvt>();                   \
        codeToCall;                                                                  \
        break;                                                                       \
    }                                                                                \
    case OilPvtApproach::NoOilPvt:                                                   \
        throw std::logic_error("Not implemented: Oil PVT of this deck!");            \
    }                                                                                \

namespace Opm {

template<class Scalar, bool enableThermal>
OilPvtMultiplexer<Scalar,enableThermal>::OilPvtMultiplexer()
{
    approach_ = OilPvtApproach::NoOilPvt;
    realOilPvt_ = nullptr;
}

template<class Scalar, bool enableThermal>
OilPvtMultiplexer<Scalar,enableThermal>::
OilPvtMultiplexer(OilPvtApproach approach, void* realOilPvt)
    : approach_(approach)
    , realOilPvt_(realOilPvt)
{
}

template<class Scalar, bool enableThermal>
OilPvtMultiplexer<Scalar,enableThermal>::
OilPvtMultiplexer(const OilPvtMultiplexer<Scalar,enableThermal>& data)
{
    *this = data;
}

template<class Scalar, bool enableThermal>
OilPvtMultiplexer<Scalar,enableThermal>::~OilPvtMultiplexer()
{
    switch (approach_) {
    case OilPvtApproach::LiveOilPvt: {
        delete &getRealPvt<OilPvtApproach::LiveOilPvt>();
        break;
    }
    case OilPvtApproach::DeadOilPvt: {
        delete &getRealPvt<OilPvtApproach::DeadOilPvt>();
        break;
    }
    case OilPvtApproach::ConstantCompressibilityOilPvt: {
        delete &getRealPvt<OilPvtApproach::ConstantCompressibilityOilPvt>();
        break;
    }
    case OilPvtApproach::ThermalOilPvt: {
        delete &getRealPvt<OilPvtApproach::ThermalOilPvt>();
        break;
    }
    case OilPvtApproach::BrineCo2Pvt: {
        delete &getRealPvt<OilPvtApproach::BrineCo2Pvt>();
        break;
    }

    case OilPvtApproach::NoOilPvt:
        break;
    }
}

#if HAVE_ECL_INPUT
template<class Scalar, bool enableThermal>
void OilPvtMultiplexer<Scalar,enableThermal>::
initFromState(const EclipseState& eclState, const Schedule& schedule)
{
    if (!eclState.runspec().phases().active(Phase::OIL))
        return;
    // TODO move the BrineCo2 approach to the waterPvtMultiplexer
    // when a proper gas-water simulator is supported
    if (eclState.runspec().co2Storage())
        setApproach(OilPvtApproach::BrineCo2Pvt);
    else if (enableThermal && eclState.getSimulationConfig().isThermal())
        setApproach(OilPvtApproach::ThermalOilPvt);
    else if (!eclState.getTableManager().getPvcdoTable().empty())
        setApproach(OilPvtApproach::ConstantCompressibilityOilPvt);
    else if (eclState.getTableManager().hasTables("PVDO"))
        setApproach(OilPvtApproach::DeadOilPvt);
    else if (!eclState.getTableManager().getPvtoTables().empty())
        setApproach(OilPvtApproach::LiveOilPvt);

    OPM_OIL_PVT_MULTIPLEXER_CALL(pvtImpl.initFromState(eclState, schedule));
}
#endif // HAVE_ECL_INPUT

template<class Scalar, bool enableThermal>
void OilPvtMultiplexer<Scalar,enableThermal>::initEnd()
{
    OPM_OIL_PVT_MULTIPLEXER_CALL(pvtImpl.initEnd());
}


template<class Scalar, bool enableThermal>
unsigned OilPvtMultiplexer<Scalar,enableThermal>::numRegions() const
{
    OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.numRegions());
    return 1;
}

template<class Scalar, bool enableThermal>
const Scalar OilPvtMultiplexer<Scalar,enableThermal>::
oilReferenceDensity(unsigned regionIdx) const
{
    OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.oilReferenceDensity(regionIdx));
    return 700.;
}

template<class Scalar, bool enableThermal>
void OilPvtMultiplexer<Scalar,enableThermal>::setApproach(OilPvtApproach appr)
{
    switch (appr) {
    case OilPvtApproach::LiveOilPvt:
        realOilPvt_ = new LiveOilPvt<Scalar>;
        break;

    case OilPvtApproach::DeadOilPvt:
        realOilPvt_ = new DeadOilPvt<Scalar>;
        break;

    case OilPvtApproach::ConstantCompressibilityOilPvt:
        realOilPvt_ = new ConstantCompressibilityOilPvt<Scalar>;
        break;

    case OilPvtApproach::ThermalOilPvt:
        realOilPvt_ = new OilPvtThermal<Scalar>;
        break;

    case OilPvtApproach::BrineCo2Pvt:
        realOilPvt_ = new BrineCo2Pvt<Scalar>;
        break;

    case OilPvtApproach::NoOilPvt:
        throw std::logic_error("Not implemented: Oil PVT of this deck!");
    }

    approach_ = appr;
}

template<class Scalar, bool enableThermal>
bool OilPvtMultiplexer<Scalar,enableThermal>::
operator==(const OilPvtMultiplexer<Scalar,enableThermal>& data) const
{
    if (this->approach() != data.approach())
        return false;

    switch (approach_) {
    case OilPvtApproach::ConstantCompressibilityOilPvt:
        return *static_cast<const ConstantCompressibilityOilPvt<Scalar>*>(realOilPvt_) ==
               *static_cast<const ConstantCompressibilityOilPvt<Scalar>*>(data.realOilPvt_);
    case OilPvtApproach::DeadOilPvt:
        return *static_cast<const DeadOilPvt<Scalar>*>(realOilPvt_) ==
               *static_cast<const DeadOilPvt<Scalar>*>(data.realOilPvt_);
    case OilPvtApproach::LiveOilPvt:
        return *static_cast<const LiveOilPvt<Scalar>*>(realOilPvt_) ==
               *static_cast<const LiveOilPvt<Scalar>*>(data.realOilPvt_);
    case OilPvtApproach::ThermalOilPvt:
        return *static_cast<const OilPvtThermal<Scalar>*>(realOilPvt_) ==
               *static_cast<const OilPvtThermal<Scalar>*>(data.realOilPvt_);
    case OilPvtApproach::BrineCo2Pvt:
        return *static_cast<const BrineCo2Pvt<Scalar>*>(realOilPvt_) ==
                *static_cast<const BrineCo2Pvt<Scalar>*>(data.realOilPvt_);
    default:
        return true;
    }
}

template<class Scalar, bool enableThermal>
OilPvtMultiplexer<Scalar,enableThermal>&
OilPvtMultiplexer<Scalar,enableThermal>::
operator=(const OilPvtMultiplexer<Scalar,enableThermal>& data)
{
    approach_ = data.approach_;
    switch (approach_) {
    case OilPvtApproach::ConstantCompressibilityOilPvt:
        realOilPvt_ = new ConstantCompressibilityOilPvt<Scalar>(*static_cast<const ConstantCompressibilityOilPvt<Scalar>*>(data.realOilPvt_));
        break;
    case OilPvtApproach::DeadOilPvt:
        realOilPvt_ = new DeadOilPvt<Scalar>(*static_cast<const DeadOilPvt<Scalar>*>(data.realOilPvt_));
        break;
    case OilPvtApproach::LiveOilPvt:
        realOilPvt_ = new LiveOilPvt<Scalar>(*static_cast<const LiveOilPvt<Scalar>*>(data.realOilPvt_));
        break;
    case OilPvtApproach::ThermalOilPvt:
        realOilPvt_ = new OilPvtThermal<Scalar>(*static_cast<const OilPvtThermal<Scalar>*>(data.realOilPvt_));
        break;
    case OilPvtApproach::BrineCo2Pvt:
        realOilPvt_ = new BrineCo2Pvt<Scalar>(*static_cast<const BrineCo2Pvt<Scalar>*>(data.realOilPvt_));
        break;
    default:
        break;
    }

    return *this;
}

template<class Scalar, bool enableThermal>
template<class Evaluation>
Evaluation OilPvtMultiplexer<Scalar,enableThermal>::
internalEnergy(unsigned regionIdx,
               const Evaluation& temperature,
               const Evaluation& pressure,
               const Evaluation& Rs) const
{
    OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.internalEnergy(regionIdx, temperature, pressure, Rs));
    return 0;
}

template<class Scalar, bool enableThermal>
template<class Evaluation>
Evaluation OilPvtMultiplexer<Scalar,enableThermal>::
viscosity(unsigned regionIdx,
          const Evaluation& temperature,
          const Evaluation& pressure,
          const Evaluation& Rs) const
{
    OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.viscosity(regionIdx, temperature, pressure, Rs));
    return 0;
}

template<class Scalar, bool enableThermal>
template<class Evaluation>
Evaluation OilPvtMultiplexer<Scalar,enableThermal>::
saturatedViscosity(unsigned regionIdx,
                   const Evaluation& temperature,
                   const Evaluation& pressure) const
{
    OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedViscosity(regionIdx, temperature, pressure));
    return 0;
}

template<class Scalar, bool enableThermal>
template<class Evaluation>
Evaluation OilPvtMultiplexer<Scalar,enableThermal>::
inverseFormationVolumeFactor(unsigned regionIdx,
                             const Evaluation& temperature,
                             const Evaluation& pressure,
                             const Evaluation& Rs) const
{
    OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.inverseFormationVolumeFactor(regionIdx, temperature, pressure, Rs));
    return 0;
}

template<class Scalar, bool enableThermal>
template<class Evaluation>
Evaluation OilPvtMultiplexer<Scalar,enableThermal>::
saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                      const Evaluation& temperature,
                                      const Evaluation& pressure) const
{
    OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedInverseFormationVolumeFactor(regionIdx, temperature, pressure));
    return 0;
}

template<class Scalar, bool enableThermal>
template<class Evaluation>
Evaluation OilPvtMultiplexer<Scalar,enableThermal>::
saturatedGasDissolutionFactor(unsigned regionIdx,
                              const Evaluation& temperature,
                              const Evaluation& pressure) const
{
    OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedGasDissolutionFactor(regionIdx, temperature, pressure));
    return 0;
}

template<class Scalar, bool enableThermal>
template<class Evaluation>
Evaluation OilPvtMultiplexer<Scalar,enableThermal>::
saturatedGasDissolutionFactor(unsigned regionIdx,
                              const Evaluation& temperature,
                              const Evaluation& pressure,
                              const Evaluation& oilSaturation,
                              const Evaluation& maxOilSaturation) const
{
    OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedGasDissolutionFactor(regionIdx, temperature, pressure, oilSaturation, maxOilSaturation));
    return 0;
}

template<class Scalar, bool enableThermal>
template<class Evaluation>
Evaluation OilPvtMultiplexer<Scalar,enableThermal>::
saturationPressure(unsigned regionIdx,
                   const Evaluation& temperature,
                   const Evaluation& Rs) const
{
    OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.saturationPressure(regionIdx, temperature, Rs));
    return 0;
}

template<class Scalar, bool enableThermal>
template<class Evaluation>
Evaluation OilPvtMultiplexer<Scalar,enableThermal>::
diffusionCoefficient(const Evaluation& temperature,
                     const Evaluation& pressure,
                     unsigned compIdx) const
{
  OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.diffusionCoefficient(temperature, pressure, compIdx));
  return 0;
}

} // namespace Opm

#undef OPM_OIL_PVT_MULTIPLEXER_CALL
