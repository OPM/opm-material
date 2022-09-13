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

#include <opm/material/fluidsystems/blackoilpvt/OilPvtMultiplexer.hpp>

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SimpleTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#endif

namespace Opm {

template<class Scalar>
OilPvtThermal<Scalar>::
OilPvtThermal(IsothermalPvt* isothermalPvt,
              const std::vector<TabulatedOneDFunction>& oilvisctCurves,
              const std::vector<Scalar>& viscrefPress,
              const std::vector<Scalar>& viscrefRs,
              const std::vector<Scalar>& viscRef,
              const std::vector<Scalar>& oildentRefTemp,
              const std::vector<Scalar>& oildentCT1,
              const std::vector<Scalar>& oildentCT2,
              const std::vector<Scalar>& oilJTRefPres,
              const std::vector<Scalar>& oilJTC,
              const std::vector<TabulatedOneDFunction>& internalEnergyCurves,
              bool enableThermalDensity,
              bool enableJouleThomson,
              bool enableThermalViscosity,
              bool enableInternalEnergy)
    : isothermalPvt_(isothermalPvt)
    , oilvisctCurves_(oilvisctCurves)
    , viscrefPress_(viscrefPress)
    , viscrefRs_(viscrefRs)
    , viscRef_(viscRef)
    , oildentRefTemp_(oildentRefTemp)
    , oildentCT1_(oildentCT1)
    , oildentCT2_(oildentCT2)
    , oilJTRefPres_(oilJTRefPres)
    , oilJTC_(oilJTC)
    , internalEnergyCurves_(internalEnergyCurves)
    , enableThermalDensity_(enableThermalDensity)
    , enableJouleThomson_(enableJouleThomson)
    , enableThermalViscosity_(enableThermalViscosity)
    , enableInternalEnergy_(enableInternalEnergy)
{
}

template<class Scalar>
 OilPvtThermal<Scalar>::OilPvtThermal(const OilPvtThermal& data)
{
     *this = data;
 }

template<class Scalar>
OilPvtThermal<Scalar>::~OilPvtThermal()
{
    delete isothermalPvt_;
}

#if HAVE_ECL_INPUT
template<class Scalar>
void OilPvtThermal<Scalar>::
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

    enableThermalDensity_ = tables.OilDenT().size() > 0;
    enableThermalViscosity_ = tables.hasTables("OILVISCT");
    enableInternalEnergy_ = tables.hasTables("SPECHEAT");

    unsigned numRegions = isothermalPvt_->numRegions();
    setNumRegions(numRegions);

    // viscosity
    if (enableThermalViscosity_) {
        if (tables.getViscrefTable().empty())
            throw std::runtime_error("VISCREF is required when OILVISCT is present");

        const auto& oilvisctTables = tables.getOilvisctTables();
        const auto& viscrefTable = tables.getViscrefTable();

        assert(oilvisctTables.size() == numRegions);
        assert(viscrefTable.size() == numRegions);

        for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            const auto& TCol = oilvisctTables[regionIdx].getColumn("Temperature").vectorCopy();
            const auto& muCol = oilvisctTables[regionIdx].getColumn("Viscosity").vectorCopy();
            oilvisctCurves_[regionIdx].setXYContainers(TCol, muCol);

            viscrefPress_[regionIdx] = viscrefTable[regionIdx].reference_pressure;
            viscrefRs_[regionIdx] = viscrefTable[regionIdx].reference_rs;

            // temperature used to calculate the reference viscosity [K]. the
            // value does not really matter if the underlying PVT object really
            // is isothermal...
            constexpr const Scalar Tref = 273.15 + 20;

            // compute the reference viscosity using the isothermal PVT object.
            viscRef_[regionIdx] =
                isothermalPvt_->viscosity(regionIdx,
                                          Tref,
                                          viscrefPress_[regionIdx],
                                          viscrefRs_[regionIdx]);
        }
    }

    // temperature dependence of oil density
    const auto& oilDenT = tables.OilDenT();
    if (oilDenT.size() > 0) {
        assert(oilDenT.size() == numRegions);
        for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            const auto& record = oilDenT[regionIdx];

            oildentRefTemp_[regionIdx] = record.T0;
            oildentCT1_[regionIdx] = record.C1;
            oildentCT2_[regionIdx] = record.C2;
        }
    }

    // Joule Thomson
    if (enableJouleThomson_) {
        const auto& oilJT = tables.OilJT();

        assert(oilJT.size() == numRegions);
        for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            const auto& record = oilJT[regionIdx];

            oilJTRefPres_[regionIdx] =  record.P0;
            oilJTC_[regionIdx] = record.C1;
        }

        const auto& densityTable = eclState.getTableManager().getDensityTable();

        assert(densityTable.size() == numRegions);
        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
             rhoRefG_[regionIdx] = densityTable[regionIdx].gas;
        }
    }

    if (enableInternalEnergy_) {
        // the specific internal energy of liquid oil. be aware that ecl only specifies the
        // heat capacity (via the SPECHEAT keyword) and we need to integrate it
        // ourselfs to get the internal energy
        for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            const auto& specheatTable = tables.getSpecheatTables()[regionIdx];
            const auto& temperatureColumn = specheatTable.getColumn("TEMPERATURE");
            const auto& cvOilColumn = specheatTable.getColumn("CV_OIL");

            std::vector<double> uSamples(temperatureColumn.size());

            Scalar u = temperatureColumn[0]*cvOilColumn[0];
            for (size_t i = 0;; ++i) {
                uSamples[i] = u;

                if (i >= temperatureColumn.size() - 1)
                    break;

                // integrate to the heat capacity from the current sampling point to the next
                // one. this leads to a quadratic polynomial.
                Scalar c_v0 = cvOilColumn[i];
                Scalar c_v1 = cvOilColumn[i + 1];
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
void OilPvtThermal<Scalar>::setNumRegions(size_t numRegions)
{
    oilvisctCurves_.resize(numRegions);
    viscrefPress_.resize(numRegions);
    viscrefRs_.resize(numRegions);
    viscRef_.resize(numRegions);
    internalEnergyCurves_.resize(numRegions);
    oildentRefTemp_.resize(numRegions);
    oildentCT1_.resize(numRegions);
    oildentCT2_.resize(numRegions);
    oilJTRefPres_.resize(numRegions);
    oilJTC_.resize(numRegions);
    rhoRefG_.resize(numRegions);
}

template<class Scalar>
bool OilPvtThermal<Scalar>::operator==(const OilPvtThermal<Scalar>& data) const
{
    if (isothermalPvt_ && !data.isothermalPvt_)
        return false;
    if (!isothermalPvt_ && data.isothermalPvt_)
        return false;

    return (!this->isoThermalPvt() ||
            (*this->isoThermalPvt() == *data.isoThermalPvt())) &&
            this->oilvisctCurves() == data.oilvisctCurves() &&
            this->viscrefPress() == data.viscrefPress() &&
            this->viscrefRs() == data.viscrefRs() &&
            this->viscRef() == data.viscRef() &&
            this->oildentRefTemp() == data.oildentRefTemp() &&
            this->oildentCT1() == data.oildentCT1() &&
            this->oildentCT2() == data.oildentCT2() &&
            this->oilJTRefPres() == data.oilJTRefPres() &&
            this->oilJTC() == data.oilJTC() &&
            this->internalEnergyCurves() == data.internalEnergyCurves() &&
            this->enableThermalDensity() == data.enableThermalDensity() &&
            this->enableJouleThomson() == data.enableJouleThomson() &&
            this->enableThermalViscosity() == data.enableThermalViscosity() &&
            this->enableInternalEnergy() == data.enableInternalEnergy();
}

template<class Scalar>
OilPvtThermal<Scalar>& OilPvtThermal<Scalar>::
operator=(const OilPvtThermal<Scalar>& data)
{
    if (data.isothermalPvt_)
        isothermalPvt_ = new IsothermalPvt(*data.isothermalPvt_);
    else
        isothermalPvt_ = nullptr;
    oilvisctCurves_ = data.oilvisctCurves_;
    viscrefPress_ = data.viscrefPress_;
    viscrefRs_ = data.viscrefRs_;
    viscRef_ = data.viscRef_;
    oildentRefTemp_ = data.oildentRefTemp_;
    oildentCT1_ = data.oildentCT1_;
    oildentCT2_ = data.oildentCT2_;
    oilJTRefPres_ =  data.oilJTRefPres_;
    oilJTC_ =  data.oilJTC_;
    internalEnergyCurves_ = data.internalEnergyCurves_;
    enableThermalDensity_ = data.enableThermalDensity_;
    enableJouleThomson_ = data.enableJouleThomson_;
    enableThermalViscosity_ = data.enableThermalViscosity_;
    enableInternalEnergy_ = data.enableInternalEnergy_;

    return *this;
}

template<class Scalar>
template<class Evaluation>
Evaluation OilPvtThermal<Scalar>::
internalEnergy(unsigned regionIdx,
               const Evaluation& temperature,
               const Evaluation& pressure,
               const Evaluation& Rs) const
{
    if (!enableInternalEnergy_)
         throw std::runtime_error("Requested the internal energy of oil but it is disabled");

    if (!enableJouleThomson_) {
        // compute the specific internal energy for the specified tempature. We use linear
        // interpolation here despite the fact that the underlying heat capacities are
        // piecewise linear (which leads to a quadratic function)
        return internalEnergyCurves_[regionIdx].eval(temperature, /*extrapolate=*/true);
    }
    else {
        Evaluation Tref = oildentRefTemp_[regionIdx];
        Evaluation Pref = oilJTRefPres_[regionIdx];
        Scalar JTC = oilJTC_[regionIdx]; // if JTC is default then JTC is calculated

        Evaluation invB = inverseFormationVolumeFactor(regionIdx, temperature, pressure, Rs);
        Evaluation Cp = internalEnergyCurves_[regionIdx].eval(temperature, /*extrapolate=*/true)/temperature;
        Evaluation density = invB * (oilReferenceDensity(regionIdx) + Rs * rhoRefG_[regionIdx]);

        Evaluation enthalpyPres;
        if  (JTC != 0) {
            enthalpyPres = -Cp * JTC * (pressure -Pref);
        }
        else if(enableThermalDensity_) {
            Scalar c1T = oildentCT1_[regionIdx];
            Scalar c2T = oildentCT2_[regionIdx];

            Evaluation alpha = (c1T + 2 * c2T * (temperature - Tref)) /
                (1 + c1T  *(temperature - Tref) + c2T * (temperature - Tref) * (temperature - Tref));

            const int N = 100; // value is experimental
            Evaluation deltaP = (pressure - Pref)/N;
            Evaluation enthalpyPresPrev = 0;
            for (size_t i = 0; i < N; ++i) {
                Evaluation Pnew = Pref + i * deltaP;
                Evaluation rho = inverseFormationVolumeFactor(regionIdx, temperature, Pnew, Rs) *
                                 (oilReferenceDensity(regionIdx) + Rs * rhoRefG_[regionIdx]) ;
                // see e.g.https://en.wikipedia.org/wiki/Joule-Thomson_effect for a derivation of the Joule-Thomson coeff.
                Evaluation jouleThomsonCoefficient = -(1.0/Cp) * (1.0 - alpha * temperature)/rho;
                Evaluation deltaEnthalpyPres = -Cp * jouleThomsonCoefficient * deltaP;
                enthalpyPres = enthalpyPresPrev + deltaEnthalpyPres;
                enthalpyPresPrev = enthalpyPres;
            }
        }
        else {
              throw std::runtime_error("Requested Joule-thomson calculation but thermal oil density (OILDENT) is not provided");
        }

        Evaluation enthalpy = Cp * (temperature - Tref) + enthalpyPres;

        return enthalpy - pressure/density;
    }
}

template<class Scalar>
template<class Evaluation>
Evaluation OilPvtThermal<Scalar>::
viscosity(unsigned regionIdx,
          const Evaluation& temperature,
          const Evaluation& pressure,
          const Evaluation& Rs) const
{
    const auto& isothermalMu = isothermalPvt_->viscosity(regionIdx, temperature, pressure, Rs);
    if (!enableThermalViscosity())
        return isothermalMu;

    // compute the viscosity deviation due to temperature
    const auto& muOilvisct = oilvisctCurves_[regionIdx].eval(temperature, /*extrapolate=*/true);
    return muOilvisct/viscRef_[regionIdx]*isothermalMu;
}

template<class Scalar>
template<class Evaluation>
Evaluation OilPvtThermal<Scalar>::
saturatedViscosity(unsigned regionIdx,
                   const Evaluation& temperature,
                   const Evaluation& pressure) const
{
    const auto& isothermalMu = isothermalPvt_->saturatedViscosity(regionIdx, temperature, pressure);
    if (!enableThermalViscosity())
        return isothermalMu;

    // compute the viscosity deviation due to temperature
    const auto& muOilvisct = oilvisctCurves_[regionIdx].eval(temperature, /*extrapolate=*/true);
    return muOilvisct/viscRef_[regionIdx]*isothermalMu;
}

template<class Scalar>
template<class Evaluation>
Evaluation OilPvtThermal<Scalar>::
inverseFormationVolumeFactor(unsigned regionIdx,
                             const Evaluation& temperature,
                             const Evaluation& pressure,
                             const Evaluation& Rs) const
{
    const auto& b =
        isothermalPvt_->inverseFormationVolumeFactor(regionIdx, temperature, pressure, Rs);

    if (!enableThermalDensity())
        return b;

    // we use the same approach as for the for water here, but with the OPM-specific
    // OILDENT keyword.
    Scalar TRef = oildentRefTemp_[regionIdx];
    Scalar cT1 = oildentCT1_[regionIdx];
    Scalar cT2 = oildentCT2_[regionIdx];
    const Evaluation& Y = temperature - TRef;

    return b/(1 + (cT1 + cT2*Y)*Y);
}

template<class Scalar>
template<class Evaluation>
Evaluation OilPvtThermal<Scalar>::
saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                      const Evaluation& temperature,
                                      const Evaluation& pressure) const
{
    const auto& b =
        isothermalPvt_->saturatedInverseFormationVolumeFactor(regionIdx, temperature, pressure);

    if (!enableThermalDensity())
        return b;

    // we use the same approach as for the for water here, but with the OPM-specific
    // OILDENT keyword.
    Scalar TRef = oildentRefTemp_[regionIdx];
    Scalar cT1 = oildentCT1_[regionIdx];
    Scalar cT2 = oildentCT2_[regionIdx];
    const Evaluation& Y = temperature - TRef;

    return b/(1 + (cT1 + cT2*Y)*Y);
}

template<class Scalar>
template<class Evaluation>
Evaluation OilPvtThermal<Scalar>::
saturatedGasDissolutionFactor(unsigned regionIdx,
                              const Evaluation& temperature,
                              const Evaluation& pressure) const
{
    return isothermalPvt_->saturatedGasDissolutionFactor(regionIdx, temperature, pressure);
}

template<class Scalar>
template<class Evaluation>
Evaluation OilPvtThermal<Scalar>::
saturatedGasDissolutionFactor(unsigned regionIdx,
                              const Evaluation& temperature,
                              const Evaluation& pressure,
                              const Evaluation& oilSaturation,
                              const Evaluation& maxOilSaturation) const
{
    return isothermalPvt_->saturatedGasDissolutionFactor(regionIdx, temperature, pressure, oilSaturation, maxOilSaturation);
}

template<class Scalar>
template<class Evaluation>
Evaluation OilPvtThermal<Scalar>::
saturationPressure(unsigned regionIdx,
                   const Evaluation& temperature,
                   const Evaluation& pressure) const
{
    return isothermalPvt_->saturationPressure(regionIdx, temperature, pressure);
}

template<class Scalar>
template<class Evaluation>
Evaluation OilPvtThermal<Scalar>::
diffusionCoefficient(const Evaluation& temperature,
                     const Evaluation& pressure,
                     unsigned compIdx) const
{
    return isothermalPvt_->diffusionCoefficient(temperature, pressure, compIdx);
}

} // namespace Opm
