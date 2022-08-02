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

#define OPM_WATER_PVT_MULTIPLEXER_CALL(codeToCall)                      \
    switch (approach_) {                                                \
    case WaterPvtApproach::ConstantCompressibilityWaterPvt: {           \
        auto& pvtImpl = getRealPvt<WaterPvtApproach::ConstantCompressibilityWaterPvt>();  \
        codeToCall;                                                     \
        break;                                                          \
    }                                                                   \
    case WaterPvtApproach::ConstantCompressibilityBrinePvt: {           \
        auto& pvtImpl = getRealPvt<WaterPvtApproach::ConstantCompressibilityBrinePvt>();  \
        codeToCall;                                                     \
        break;                                                          \
    }                                                                   \
    case WaterPvtApproach::ThermalWaterPvt: {                           \
        auto& pvtImpl = getRealPvt<WaterPvtApproach::ThermalWaterPvt>();                  \
        codeToCall;                                                     \
        break;                                                          \
    }                                                                   \
    case WaterPvtApproach::NoWaterPvt:                                  \
        throw std::logic_error("Not implemented: Water PVT of this deck!"); \
    }

namespace Opm {

template<class Scalar, bool enableThermal, bool enableBrine>
WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>::WaterPvtMultiplexer()
{
    approach_ = WaterPvtApproach::NoWaterPvt;
    realWaterPvt_ = nullptr;
}

template<class Scalar, bool enableThermal, bool enableBrine>
WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>::
WaterPvtMultiplexer(WaterPvtApproach approach, void* realWaterPvt)
    : approach_(approach)
    , realWaterPvt_(realWaterPvt)
{
}

template<class Scalar, bool enableThermal, bool enableBrine>
 WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>::
 WaterPvtMultiplexer(const WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>& data)
{
    *this = data;
}

template<class Scalar, bool enableThermal, bool enableBrine>
WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>::~WaterPvtMultiplexer()
{
    switch (approach_) {
    case WaterPvtApproach::ConstantCompressibilityWaterPvt: {
        delete &getRealPvt<WaterPvtApproach::ConstantCompressibilityWaterPvt>();
        break;
    }
    case WaterPvtApproach::ConstantCompressibilityBrinePvt: {
        delete &getRealPvt<WaterPvtApproach::ConstantCompressibilityBrinePvt>();
        break;
    }
    case WaterPvtApproach::ThermalWaterPvt: {
        delete &getRealPvt<WaterPvtApproach::ThermalWaterPvt>();
        break;
    }
    case WaterPvtApproach::NoWaterPvt:
        break;
    }
}

#if HAVE_ECL_INPUT
template<class Scalar, bool enableThermal, bool enableBrine>
void WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>::
initFromState(const EclipseState& eclState, const Schedule& schedule)
{
    if (!eclState.runspec().phases().active(Phase::WATER))
        return;

    if (enableThermal && eclState.getSimulationConfig().isThermal())
        setApproach(WaterPvtApproach::ThermalWaterPvt);
    else if (!eclState.getTableManager().getPvtwTable().empty())
        setApproach(WaterPvtApproach::ConstantCompressibilityWaterPvt);
    else if (enableBrine && !eclState.getTableManager().getPvtwSaltTables().empty())
        setApproach(WaterPvtApproach::ConstantCompressibilityBrinePvt);

    OPM_WATER_PVT_MULTIPLEXER_CALL(pvtImpl.initFromState(eclState, schedule));
}
#endif // HAVE_ECL_INPUT

template<class Scalar, bool enableThermal, bool enableBrine>
void WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>::initEnd()
{
    OPM_WATER_PVT_MULTIPLEXER_CALL(pvtImpl.initEnd());
}

template<class Scalar, bool enableThermal, bool enableBrine>
unsigned WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>::numRegions() const
{
    OPM_WATER_PVT_MULTIPLEXER_CALL(return pvtImpl.numRegions());
    return 1;
}

template<class Scalar, bool enableThermal, bool enableBrine>
const Scalar WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>::
waterReferenceDensity(unsigned regionIdx)
{
     OPM_WATER_PVT_MULTIPLEXER_CALL(return pvtImpl.waterReferenceDensity(regionIdx));
     return 1000.;
}

template<class Scalar, bool enableThermal, bool enableBrine>
void WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>::
setApproach(WaterPvtApproach appr)
{
    switch (appr) {
    case WaterPvtApproach::ConstantCompressibilityWaterPvt:
        realWaterPvt_ = new ConstantCompressibilityWaterPvt<Scalar>;
        break;

    case WaterPvtApproach::ConstantCompressibilityBrinePvt:
        realWaterPvt_ = new ConstantCompressibilityBrinePvt<Scalar>;
        break;

    case WaterPvtApproach::ThermalWaterPvt:
        realWaterPvt_ = new WaterPvtThermal<Scalar, enableBrine>;
        break;

    case WaterPvtApproach::NoWaterPvt:
        throw std::logic_error("Not implemented: Water PVT of this deck!");
    }

    approach_ = appr;
}

template<class Scalar, bool enableThermal, bool enableBrine>
bool WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>::
operator==(const WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>& data) const
{
    if (this->approach() != data.approach())
        return false;

    switch (approach_) {
    case WaterPvtApproach::ConstantCompressibilityWaterPvt:
        return *static_cast<const ConstantCompressibilityWaterPvt<Scalar>*>(realWaterPvt_) ==
               *static_cast<const ConstantCompressibilityWaterPvt<Scalar>*>(data.realWaterPvt_);
    case WaterPvtApproach::ConstantCompressibilityBrinePvt:
        return *static_cast<const ConstantCompressibilityBrinePvt<Scalar>*>(realWaterPvt_) ==
               *static_cast<const ConstantCompressibilityBrinePvt<Scalar>*>(data.realWaterPvt_);
    case WaterPvtApproach::ThermalWaterPvt:
        return *static_cast<const WaterPvtThermal<Scalar, enableBrine>*>(realWaterPvt_) ==
               *static_cast<const WaterPvtThermal<Scalar, enableBrine>*>(data.realWaterPvt_);
    default:
        return true;
    }
}

template<class Scalar, bool enableThermal, bool enableBrine>
WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>&
WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>::
operator=(const WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>& data)
{
    approach_ = data.approach_;
    switch (approach_) {
    case WaterPvtApproach::ConstantCompressibilityWaterPvt:
        realWaterPvt_ = new ConstantCompressibilityWaterPvt<Scalar>(*static_cast<const ConstantCompressibilityWaterPvt<Scalar>*>(data.realWaterPvt_));
        break;
    case WaterPvtApproach::ConstantCompressibilityBrinePvt:
        realWaterPvt_ = new ConstantCompressibilityBrinePvt<Scalar>(*static_cast<const ConstantCompressibilityBrinePvt<Scalar>*>(data.realWaterPvt_));
        break;
    case WaterPvtApproach::ThermalWaterPvt:
        realWaterPvt_ = new WaterPvtThermal<Scalar, enableBrine>(*static_cast<const WaterPvtThermal<Scalar, enableBrine>*>(data.realWaterPvt_));
        break;
    default:
        break;
    }

    return *this;
}

template<class Scalar, bool enableThermal, bool enableBrine>
template<class Evaluation>
Evaluation WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>::
internalEnergy(unsigned regionIdx,
               const Evaluation& temperature,
               const Evaluation& pressure,
               const Evaluation& saltconcentration) const
{
    OPM_WATER_PVT_MULTIPLEXER_CALL(return pvtImpl.internalEnergy(regionIdx, temperature, pressure, saltconcentration));
    return 0;
}

template<class Scalar, bool enableThermal, bool enableBrine>
template<class Evaluation>
Evaluation WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>::
viscosity(unsigned regionIdx,
          const Evaluation& temperature,
          const Evaluation& pressure,
          const Evaluation& saltconcentration) const
{
    OPM_WATER_PVT_MULTIPLEXER_CALL(return pvtImpl.viscosity(regionIdx, temperature, pressure, saltconcentration));
    return 0;
}

template<class Scalar, bool enableThermal, bool enableBrine>
template<class Evaluation>
Evaluation WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>::
inverseFormationVolumeFactor(unsigned regionIdx,
                             const Evaluation& temperature,
                             const Evaluation& pressure,
                             const Evaluation& saltconcentration) const
{
    OPM_WATER_PVT_MULTIPLEXER_CALL(return pvtImpl.inverseFormationVolumeFactor(regionIdx, temperature, pressure, saltconcentration));
    return 0;
}

} // namespace Opm

#undef OPM_WATER_PVT_MULTIPLEXER_CALL
