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

#include <opm/material/fluidsystems/blackoilpvt/GasPvtMultiplexer.hpp>

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SimpleTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#endif

namespace Opm {

template<class Scalar>
GasPvtThermal<Scalar>::
GasPvtThermal(IsothermalPvt* isothermalPvt,
              const std::vector<TabulatedOneDFunction>& gasvisctCurves,
              const std::vector<Scalar>& gasdentRefTemp,
              const std::vector<Scalar>& gasdentCT1,
              const std::vector<Scalar>& gasdentCT2,
              const std::vector<Scalar>& gasJTRefPres,
              const std::vector<Scalar>& gasJTC,
              const std::vector<TabulatedOneDFunction>& internalEnergyCurves,
              bool enableThermalDensity,
              bool enableJouleThomson,
              bool enableThermalViscosity,
              bool enableInternalEnergy)
    : isothermalPvt_(isothermalPvt)
    , gasvisctCurves_(gasvisctCurves)
    , gasdentRefTemp_(gasdentRefTemp)
    , gasdentCT1_(gasdentCT1)
    , gasdentCT2_(gasdentCT2)
    , gasJTRefPres_(gasJTRefPres)
    , gasJTC_(gasJTC)
    , internalEnergyCurves_(internalEnergyCurves)
    , enableThermalDensity_(enableThermalDensity)
    , enableJouleThomson_(enableJouleThomson)
    , enableThermalViscosity_(enableThermalViscosity)
    , enableInternalEnergy_(enableInternalEnergy)
{
}

template<class Scalar>
GasPvtThermal<Scalar>::GasPvtThermal(const GasPvtThermal& data)
{
    *this = data;
}

template<class Scalar>
GasPvtThermal<Scalar>::~GasPvtThermal()
{
    delete isothermalPvt_;
}

#if HAVE_ECL_INPUT
template<class Scalar>
void GasPvtThermal<Scalar>::
initFromState(const EclipseState& eclState, const Schedule& schedule)
{
    //////
    // initialize the isothermal part
    //////
    isothermalPvt_ = new IsothermalPvt;
    isothermalPvt_->initFromState(eclState, schedule);

    //////
    // initialize the thermal part
    //////
    const auto& tables = eclState.getTableManager();

    enableThermalDensity_ = tables.GasDenT().size() > 0;
    enableJouleThomson_ = tables.GasJT().size() > 0;
    enableThermalViscosity_ = tables.hasTables("GASVISCT");
    enableInternalEnergy_ = tables.hasTables("SPECHEAT");

    unsigned numRegions = isothermalPvt_->numRegions();
    setNumRegions(numRegions);

    // viscosity
    if (enableThermalViscosity_) {
        const auto& gasvisctTables = tables.getGasvisctTables();
        auto gasCompIdx = tables.gas_comp_index();
        std::string gasvisctColumnName = "Viscosity" + std::to_string(static_cast<long long>(gasCompIdx));

        for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            const auto& T = gasvisctTables[regionIdx].getColumn("Temperature").vectorCopy();
            const auto& mu = gasvisctTables[regionIdx].getColumn(gasvisctColumnName).vectorCopy();
            gasvisctCurves_[regionIdx].setXYContainers(T, mu);
        }
    }

    // temperature dependence of gas density
    if (enableThermalDensity_) {
        const auto& gasDenT = tables.GasDenT();

        assert(gasDenT.size() == numRegions);
        for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            const auto& record = gasDenT[regionIdx];

            gasdentRefTemp_[regionIdx] = record.T0;
            gasdentCT1_[regionIdx] = record.C1;
            gasdentCT2_[regionIdx] = record.C2;
        }
    }

    // Joule Thomson
    if (enableJouleThomson_) {
        const auto& gasJT = tables.GasJT();

        assert(gasJT.size() == numRegions);
        for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            const auto& record = gasJT[regionIdx];

            gasJTRefPres_[regionIdx] =  record.P0;
            gasJTC_[regionIdx] = record.C1;
        }

        const auto& densityTable = eclState.getTableManager().getDensityTable();

        assert(densityTable.size() == numRegions);
        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
             rhoRefO_[regionIdx] = densityTable[regionIdx].oil;
        }
    }

    if (enableInternalEnergy_) {
        // the specific internal energy of gas. be aware that ecl only specifies the heat capacity
        // (via the SPECHEAT keyword) and we need to integrate it ourselfs to get the
        // internal energy
        for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            const auto& specHeatTable = tables.getSpecheatTables()[regionIdx];
            const auto& temperatureColumn = specHeatTable.getColumn("TEMPERATURE");
            const auto& cvGasColumn = specHeatTable.getColumn("CV_GAS");

            std::vector<double> uSamples(temperatureColumn.size());

            // the specific enthalpy of vaporization. since ECL does not seem to
            // feature a proper way to specify this quantity, we use the value for
            // methane. A proper model would also need to consider the enthalpy
            // change due to dissolution, i.e. the enthalpies of the gas and oil
            // phases should depend on the phase composition
            constexpr const Scalar hVap = 480.6e3; // [J / kg]

            Scalar u = temperatureColumn[0]*cvGasColumn[0] + hVap;
            for (size_t i = 0;; ++i) {
                uSamples[i] = u;

                if (i >= temperatureColumn.size() - 1)
                    break;

                // integrate to the heat capacity from the current sampling point to the next
                // one. this leads to a quadratic polynomial.
                Scalar c_v0 = cvGasColumn[i];
                Scalar c_v1 = cvGasColumn[i + 1];
                Scalar T0 = temperatureColumn[i];
                Scalar T1 = temperatureColumn[i + 1];
                u += 0.5*(c_v0 + c_v1)*(T1 - T0);
            }

            internalEnergyCurves_[regionIdx].setXYContainers(temperatureColumn.vectorCopy(), uSamples);
        }
    }
}
#endif // HAVE_ECL_INPUT

template<class Scalar>
void GasPvtThermal<Scalar>::setNumRegions(size_t numRegions)
{
    gasvisctCurves_.resize(numRegions);
    internalEnergyCurves_.resize(numRegions);
    gasdentRefTemp_.resize(numRegions);
    gasdentCT1_.resize(numRegions);
    gasdentCT2_.resize(numRegions);
    gasJTRefPres_.resize(numRegions);
    gasJTC_.resize(numRegions);
    rhoRefO_.resize(numRegions);
}

template<class Scalar>
bool GasPvtThermal<Scalar>::operator==(const GasPvtThermal<Scalar>& data) const
{
    if (isothermalPvt_ && !data.isothermalPvt_)
        return false;
    if (!isothermalPvt_ && data.isothermalPvt_)
        return false;

    return (!this->isoThermalPvt() ||
            (*this->isoThermalPvt() == *data.isoThermalPvt())) &&
            this->gasvisctCurves() == data.gasvisctCurves() &&
            this->gasdentRefTemp() == data.gasdentRefTemp() &&
            this->gasdentCT1() == data.gasdentCT1() &&
            this->gasdentCT2() == data.gasdentCT2() &&
            this->gasJTRefPres() == data.gasJTRefPres() &&
            this->gasJTC() == data.gasJTC() &&
            this->internalEnergyCurves() == data.internalEnergyCurves() &&
            this->enableThermalDensity() == data.enableThermalDensity() &&
            this->enableJouleThomson() == data.enableJouleThomson() &&
            this->enableThermalViscosity() == data.enableThermalViscosity() &&
            this->enableInternalEnergy() == data.enableInternalEnergy();
}

template<class Scalar>
GasPvtThermal<Scalar>&
GasPvtThermal<Scalar>::operator=(const GasPvtThermal<Scalar>& data)
{
    if (data.isothermalPvt_)
        isothermalPvt_ = new IsothermalPvt(*data.isothermalPvt_);
    else
        isothermalPvt_ = nullptr;
    gasvisctCurves_ = data.gasvisctCurves_;
    gasdentRefTemp_ = data.gasdentRefTemp_;
    gasdentCT1_ = data.gasdentCT1_;
    gasdentCT2_ = data.gasdentCT2_;
    gasJTRefPres_ =  data.gasJTRefPres_;
    gasJTC_ =  data.gasJTC_;
    internalEnergyCurves_ = data.internalEnergyCurves_;
    enableThermalDensity_ = data.enableThermalDensity_;
    enableJouleThomson_ = data.enableJouleThomson_;
    enableThermalViscosity_ = data.enableThermalViscosity_;
    enableInternalEnergy_ = data.enableInternalEnergy_;

    return *this;
}

template<class Scalar>
template<class Evaluation>
Evaluation GasPvtThermal<Scalar>::
internalEnergy(unsigned regionIdx,
               const Evaluation& temperature,
               const Evaluation& pressure,
               const Evaluation& Rv) const
{
    if (!enableInternalEnergy_)
         throw std::runtime_error("Requested the internal energy of gas but it is disabled");

    if (!enableJouleThomson_) {
        // compute the specific internal energy for the specified tempature. We use linear
        // interpolation here despite the fact that the underlying heat capacities are
        // piecewise linear (which leads to a quadratic function)
        return internalEnergyCurves_[regionIdx].eval(temperature, /*extrapolate=*/true);
    }
    else {
        Evaluation Tref = gasdentRefTemp_[regionIdx];
        Evaluation Pref = gasJTRefPres_[regionIdx];
        Scalar JTC = gasJTC_[regionIdx]; // if JTC is default then JTC is calculated
        Evaluation Rvw = 0.0;

        Evaluation invB = inverseFormationVolumeFactor(regionIdx, temperature, pressure, Rv, Rvw);
        const Scalar hVap = 480.6e3; // [J / kg]
        Evaluation Cp = (internalEnergyCurves_[regionIdx].eval(temperature, /*extrapolate=*/true) - hVap)/temperature;
        Evaluation density = invB * (gasReferenceDensity(regionIdx) + Rv * rhoRefO_[regionIdx]);

        Evaluation enthalpyPres;
        if  (JTC != 0) {
            enthalpyPres = -Cp * JTC * (pressure -Pref);
        }
        else if(enableThermalDensity_) {
            Scalar c1T = gasdentCT1_[regionIdx];
            Scalar c2T = gasdentCT2_[regionIdx];

            Evaluation alpha = (c1T + 2 * c2T * (temperature - Tref)) /
                (1 + c1T  *(temperature - Tref) + c2T * (temperature - Tref) * (temperature - Tref));

            constexpr const int N = 100; // value is experimental
            Evaluation deltaP = (pressure - Pref)/N;
            Evaluation enthalpyPresPrev = 0;
            for (size_t i = 0; i < N; ++i) {
                Evaluation Pnew = Pref + i * deltaP;
                Evaluation rho = inverseFormationVolumeFactor(regionIdx, temperature, Pnew, Rv, Rvw) *
                                 (gasReferenceDensity(regionIdx) + Rv * rhoRefO_[regionIdx]);
                // see e.g.https://en.wikipedia.org/wiki/Joule-Thomson_effect for a derivation of the Joule-Thomson coeff.
                Evaluation jouleThomsonCoefficient = -(1.0/Cp) * (1.0 - alpha * temperature)/rho;
                Evaluation deltaEnthalpyPres = -Cp * jouleThomsonCoefficient * deltaP;
                enthalpyPres = enthalpyPresPrev + deltaEnthalpyPres;
                enthalpyPresPrev = enthalpyPres;
            }
        }
        else {
              throw std::runtime_error("Requested Joule-thomson calculation but thermal gas density (GASDENT) is not provided");
        }

        Evaluation enthalpy = Cp * (temperature - Tref) + enthalpyPres;

        return enthalpy - pressure/density;
    }
}

template<class Scalar>
template<class Evaluation>
Evaluation GasPvtThermal<Scalar>::
viscosity(unsigned regionIdx,
          const Evaluation& temperature,
          const Evaluation& pressure,
          const Evaluation& Rv,
          const Evaluation& Rvw) const
{
    if (!enableThermalViscosity())
        return isothermalPvt_->viscosity(regionIdx, temperature, pressure, Rv, Rvw);

    // compute the viscosity deviation due to temperature
    const auto& muGasvisct = gasvisctCurves_[regionIdx].eval(temperature);
    return muGasvisct;
}

template<class Scalar>
template<class Evaluation>
Evaluation GasPvtThermal<Scalar>::
saturatedViscosity(unsigned regionIdx,
                   const Evaluation& temperature,
                   const Evaluation& pressure) const
{
    if (!enableThermalViscosity())
        return isothermalPvt_->saturatedViscosity(regionIdx, temperature, pressure);

    // compute the viscosity deviation due to temperature
    const auto& muGasvisct = gasvisctCurves_[regionIdx].eval(temperature, true);
    return muGasvisct;
}

template<class Scalar>
template<class Evaluation>
Evaluation GasPvtThermal<Scalar>::
inverseFormationVolumeFactor(unsigned regionIdx,
                             const Evaluation& temperature,
                             const Evaluation& pressure,
                             const Evaluation& Rv,
                             const Evaluation& /*Rvw*/) const
{
    const Evaluation& Rvw = 0.0;
    const auto& b =
        isothermalPvt_->inverseFormationVolumeFactor(regionIdx, temperature, pressure, Rv, Rvw);

    if (!enableThermalDensity())
        return b;

    // we use the same approach as for the for water here, but with the OPM-specific
    // GASDENT keyword.
    //
    // TODO: Since gas is quite a bit more compressible than water, it might be
    //       necessary to make GASDENT to a table keyword. If the current temperature
    //       is relatively close to the reference temperature, the current approach
    //       should be good enough, though.
    Scalar TRef = gasdentRefTemp_[regionIdx];
    Scalar cT1 = gasdentCT1_[regionIdx];
    Scalar cT2 = gasdentCT2_[regionIdx];
    const Evaluation& Y = temperature - TRef;

    return b/(1 + (cT1 + cT2*Y)*Y);
}

template<class Scalar>
template<class Evaluation>
Evaluation GasPvtThermal<Scalar>::
saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                      const Evaluation& temperature,
                                      const Evaluation& pressure) const
{
    const auto& b =
        isothermalPvt_->saturatedInverseFormationVolumeFactor(regionIdx, temperature, pressure);

    if (!enableThermalDensity())
        return b;

    // we use the same approach as for the for water here, but with the OPM-specific
    // GASDENT keyword.
    //
    // TODO: Since gas is quite a bit more compressible than water, it might be
    //       necessary to make GASDENT to a table keyword. If the current temperature
    //       is relatively close to the reference temperature, the current approach
    //       should be good enough, though.
    Scalar TRef = gasdentRefTemp_[regionIdx];
    Scalar cT1 = gasdentCT1_[regionIdx];
    Scalar cT2 = gasdentCT2_[regionIdx];
    const Evaluation& Y = temperature - TRef;

    return b/(1 + (cT1 + cT2*Y)*Y);
}

template<class Scalar>
template<class Evaluation>
Evaluation GasPvtThermal<Scalar>::
saturatedWaterVaporizationFactor(unsigned /*regionIdx*/,
                                 const Evaluation& /*temperature*/,
                                 const Evaluation& /*pressure*/) const
{
    return 0.0;
}

template<class Scalar>
template<class Evaluation>
Evaluation GasPvtThermal<Scalar>::
saturatedOilVaporizationFactor(unsigned regionIdx,
                               const Evaluation& temperature,
                               const Evaluation& pressure) const
{
    return isothermalPvt_->saturatedOilVaporizationFactor(regionIdx, temperature, pressure);
}

template<class Scalar>
template<class Evaluation>
Evaluation GasPvtThermal<Scalar>::
saturatedOilVaporizationFactor(unsigned regionIdx,
                               const Evaluation& temperature,
                               const Evaluation& pressure,
                               const Evaluation& oilSaturation,
                               const Evaluation& maxOilSaturation) const
{
    return isothermalPvt_->saturatedOilVaporizationFactor(regionIdx, temperature, pressure, oilSaturation, maxOilSaturation);
}

template<class Scalar>
template<class Evaluation>
Evaluation GasPvtThermal<Scalar>::
saturationPressure(unsigned regionIdx,
                   const Evaluation& temperature,
                   const Evaluation& pressure) const
{
    return isothermalPvt_->saturationPressure(regionIdx, temperature, pressure);
}

template<class Scalar>
template<class Evaluation>
Evaluation GasPvtThermal<Scalar>::
diffusionCoefficient(const Evaluation& temperature,
                     const Evaluation& pressure,
                     unsigned compIdx) const
{
    return isothermalPvt_->diffusionCoefficient(temperature, pressure, compIdx);
}

} // namespace Opm
