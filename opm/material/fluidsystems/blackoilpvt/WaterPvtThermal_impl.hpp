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

#include <opm/material/fluidsystems/blackoilpvt/WaterPvtMultiplexer.hpp>

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#endif

namespace Opm {

template <class Scalar, bool enableBrine>
WaterPvtThermal<Scalar,enableBrine>::
WaterPvtThermal(IsothermalPvt* isothermalPvt,
                const std::vector<Scalar>& viscrefPress,
                const std::vector<Scalar>& watdentRefTemp,
                const std::vector<Scalar>& watdentCT1,
                const std::vector<Scalar>& watdentCT2,
                const std::vector<Scalar>& watJTRefPres,
                const std::vector<Scalar>& watJTC,
                const std::vector<Scalar>& pvtwRefPress,
                const std::vector<Scalar>& pvtwRefB,
                const std::vector<Scalar>& pvtwCompressibility,
                const std::vector<Scalar>& pvtwViscosity,
                const std::vector<Scalar>& pvtwViscosibility,
                const std::vector<TabulatedOneDFunction>& watvisctCurves,
                const std::vector<TabulatedOneDFunction>& internalEnergyCurves,
                bool enableThermalDensity,
                bool enableJouleThomson,
                bool enableThermalViscosity,
                bool enableInternalEnergy)
    : isothermalPvt_(isothermalPvt)
    , viscrefPress_(viscrefPress)
    , watdentRefTemp_(watdentRefTemp)
    , watdentCT1_(watdentCT1)
    , watdentCT2_(watdentCT2)
    , watJTRefPres_(watJTRefPres)
    , watJTC_(watJTC)
    , pvtwRefPress_(pvtwRefPress)
    , pvtwRefB_(pvtwRefB)
    , pvtwCompressibility_(pvtwCompressibility)
    , pvtwViscosity_(pvtwViscosity)
    , pvtwViscosibility_(pvtwViscosibility)
    , watvisctCurves_(watvisctCurves)
    , internalEnergyCurves_(internalEnergyCurves)
    , enableThermalDensity_(enableThermalDensity)
    , enableJouleThomson_(enableJouleThomson)
    , enableThermalViscosity_(enableThermalViscosity)
    , enableInternalEnergy_(enableInternalEnergy)
{
}

template<class Scalar, bool enableBrine>
WaterPvtThermal<Scalar,enableBrine>::WaterPvtThermal(const WaterPvtThermal& data)
{
    *this = data;
}

template<class Scalar, bool enableBrine>
WaterPvtThermal<Scalar,enableBrine>::~WaterPvtThermal()
{
    delete isothermalPvt_;
}

#if HAVE_ECL_INPUT
template<class Scalar, bool enableBrine>
void WaterPvtThermal<Scalar,enableBrine>::
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

    enableThermalDensity_ = tables.WatDenT().size() > 0;
    enableJouleThomson_ = tables.WatJT().size() > 0;
    enableThermalViscosity_ = tables.hasTables("WATVISCT");
    enableInternalEnergy_ = tables.hasTables("SPECHEAT");

    unsigned numRegions = isothermalPvt_->numRegions();
    setNumRegions(numRegions);

    if (enableThermalDensity_) {
        const auto& watDenT = tables.WatDenT();

        assert(watDenT.size() == numRegions);
        for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            const auto& record = watDenT[regionIdx];

            watdentRefTemp_[regionIdx] = record.T0;
            watdentCT1_[regionIdx] = record.C1;
            watdentCT2_[regionIdx] = record.C2;
        }

        const auto& pvtwTables = tables.getPvtwTable();

        assert(pvtwTables.size() == numRegions);

        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            pvtwRefPress_[regionIdx] = pvtwTables[regionIdx].reference_pressure;
            pvtwRefB_[regionIdx] = pvtwTables[regionIdx].volume_factor;
        }
    }

    // Joule Thomson
    if (enableJouleThomson_) {
         const auto& watJT = tables.WatJT();

        assert(watJT.size() == numRegions);
        for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            const auto& record = watJT[regionIdx];

            watJTRefPres_[regionIdx] =  record.P0;
            watJTC_[regionIdx] = record.C1;
        }
    }

    if (enableThermalViscosity_) {
        if (tables.getViscrefTable().empty())
            throw std::runtime_error("VISCREF is required when WATVISCT is present");

        const auto& watvisctTables = tables.getWatvisctTables();
        const auto& viscrefTables = tables.getViscrefTable();

        const auto& pvtwTables = tables.getPvtwTable();

        assert(pvtwTables.size() == numRegions);
        assert(watvisctTables.size() == numRegions);
        assert(viscrefTables.size() == numRegions);

        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            const auto& T = watvisctTables[regionIdx].getColumn("Temperature").vectorCopy();
            const auto& mu = watvisctTables[regionIdx].getColumn("Viscosity").vectorCopy();
            watvisctCurves_[regionIdx].setXYContainers(T, mu);

            viscrefPress_[regionIdx] = viscrefTables[regionIdx].reference_pressure;
        }

        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            pvtwViscosity_[regionIdx] = pvtwTables[regionIdx].viscosity;
            pvtwViscosibility_[regionIdx] = pvtwTables[regionIdx].viscosibility;
        }
    }

    if (enableInternalEnergy_) {
        // the specific internal energy of liquid water. be aware that ecl only specifies the heat capacity
        // (via the SPECHEAT keyword) and we need to integrate it ourselfs to get the
        // internal energy
        for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            const auto& specHeatTable = tables.getSpecheatTables()[regionIdx];
            const auto& temperatureColumn = specHeatTable.getColumn("TEMPERATURE");
            const auto& cvWaterColumn = specHeatTable.getColumn("CV_WATER");

            std::vector<double> uSamples(temperatureColumn.size());

            Scalar u = temperatureColumn[0]*cvWaterColumn[0];
            for (size_t i = 0;; ++i) {
                uSamples[i] = u;

                if (i >= temperatureColumn.size() - 1)
                    break;

                // integrate to the heat capacity from the current sampling point to the next
                // one. this leads to a quadratic polynomial.
                Scalar c_v0 = cvWaterColumn[i];
                Scalar c_v1 = cvWaterColumn[i + 1];
                Scalar T0 = temperatureColumn[i];
                Scalar T1 = temperatureColumn[i + 1];
                u += 0.5*(c_v0 + c_v1)*(T1 - T0);
            }

            internalEnergyCurves_[regionIdx].setXYContainers(temperatureColumn.vectorCopy(), uSamples);
        }
    }
}
#endif // HAVE_ECL_INPUT

template<class Scalar, bool enableBrine>
void WaterPvtThermal<Scalar,enableBrine>::setNumRegions(size_t numRegions)
{
    pvtwRefPress_.resize(numRegions);
    pvtwRefB_.resize(numRegions);
    pvtwCompressibility_.resize(numRegions);
    pvtwViscosity_.resize(numRegions);
    pvtwViscosibility_.resize(numRegions);
    viscrefPress_.resize(numRegions);
    watvisctCurves_.resize(numRegions);
    watdentRefTemp_.resize(numRegions);
    watdentCT1_.resize(numRegions);
    watdentCT2_.resize(numRegions);
    watJTRefPres_.resize(numRegions);
    watJTC_.resize(numRegions);
    internalEnergyCurves_.resize(numRegions);
}

template<class Scalar, bool enableBrine>
bool WaterPvtThermal<Scalar,enableBrine>::
operator==(const WaterPvtThermal<Scalar,enableBrine>& data) const
{
    if (isothermalPvt_ && !data.isothermalPvt_)
        return false;
    if (!isothermalPvt_ && data.isothermalPvt_)
        return false;

    return (!this->isoThermalPvt() ||
           (*this->isoThermalPvt() == *data.isoThermalPvt())) &&
           this->viscrefPress() == data.viscrefPress() &&
           this->watdentRefTemp() == data.watdentRefTemp() &&
           this->watdentCT1() == data.watdentCT1() &&
           this->watdentCT2() == data.watdentCT2() &&
           this->watdentCT2() == data.watdentCT2() &&
           this->watJTRefPres() == data.watJTRefPres() &&
           this->watJTC() == data.watJTC() &&
           this->pvtwRefPress() == data.pvtwRefPress() &&
           this->pvtwRefB() == data.pvtwRefB() &&
           this->pvtwCompressibility() == data.pvtwCompressibility() &&
           this->pvtwViscosity() == data.pvtwViscosity() &&
           this->pvtwViscosibility() == data.pvtwViscosibility() &&
           this->watvisctCurves() == data.watvisctCurves() &&
           this->internalEnergyCurves() == data.internalEnergyCurves() &&
           this->enableThermalDensity() == data.enableThermalDensity() &&
           this->enableJouleThomson() == data.enableJouleThomson() &&
           this->enableThermalViscosity() == data.enableThermalViscosity() &&
           this->enableInternalEnergy() == data.enableInternalEnergy();
}

template<class Scalar, bool enableBrine>
WaterPvtThermal<Scalar,enableBrine>&
WaterPvtThermal<Scalar,enableBrine>::
operator=(const WaterPvtThermal<Scalar,enableBrine>& data)
{
    if (data.isothermalPvt_)
        isothermalPvt_ = new IsothermalPvt(*data.isothermalPvt_);
    else
        isothermalPvt_ = nullptr;
    viscrefPress_ = data.viscrefPress_;
    watdentRefTemp_ = data.watdentRefTemp_;
    watdentCT1_ = data.watdentCT1_;
    watdentCT2_ = data.watdentCT2_;
    watJTRefPres_ =  data.watJTRefPres_;
    watJTC_ =  data.watJTC_;
    pvtwRefPress_ = data.pvtwRefPress_;
    pvtwRefB_ = data.pvtwRefB_;
    pvtwCompressibility_ = data.pvtwCompressibility_;
    pvtwViscosity_ = data.pvtwViscosity_;
    pvtwViscosibility_ = data.pvtwViscosibility_;
    watvisctCurves_ = data.watvisctCurves_;
    internalEnergyCurves_ = data.internalEnergyCurves_;
    enableThermalDensity_ = data.enableThermalDensity_;
    enableJouleThomson_ = data.enableJouleThomson_;
    enableThermalViscosity_ = data.enableThermalViscosity_;
    enableInternalEnergy_ = data.enableInternalEnergy_;

    return *this;
}

template<class Scalar, bool enableBrine>
template<class Evaluation>
Evaluation WaterPvtThermal<Scalar,enableBrine>::
internalEnergy(unsigned regionIdx,
               const Evaluation& temperature,
               const Evaluation& pressure,
               const Evaluation& saltconcentration) const
{
    if (!enableInternalEnergy_)
        throw std::runtime_error("Requested the internal energy of water but it is disabled");

    if (!enableJouleThomson_) {
        // compute the specific internal energy for the specified tempature. We use linear
        // interpolation here despite the fact that the underlying heat capacities are
        // piecewise linear (which leads to a quadratic function)
        return internalEnergyCurves_[regionIdx].eval(temperature, /*extrapolate=*/true);
    }
    else {
        Evaluation Tref = watdentRefTemp_[regionIdx];
        Evaluation Pref = watJTRefPres_[regionIdx];
        Scalar JTC =watJTC_[regionIdx]; // if JTC is default then JTC is calculated

        Evaluation invB = inverseFormationVolumeFactor(regionIdx, temperature, pressure, saltconcentration);
        Evaluation Cp = internalEnergyCurves_[regionIdx].eval(temperature, /*extrapolate=*/true)/temperature;
        Evaluation density = invB * waterReferenceDensity(regionIdx);

        Evaluation enthalpyPres;
        if  (JTC != 0) {
            enthalpyPres = -Cp * JTC * (pressure -Pref);
        }
        else if(enableThermalDensity_) {
            Scalar c1T = watdentCT1_[regionIdx];
            Scalar c2T = watdentCT2_[regionIdx];

            Evaluation alpha = (c1T + 2 * c2T * (temperature - Tref)) /
                (1 + c1T  *(temperature - Tref) + c2T * (temperature - Tref) * (temperature - Tref));

            constexpr const int N = 100; // value is experimental
            Evaluation deltaP = (pressure - Pref)/N;
            Evaluation enthalpyPresPrev = 0;
            for (size_t i = 0; i < N; ++i) {
                Evaluation Pnew = Pref + i * deltaP;
                Evaluation rho = inverseFormationVolumeFactor(regionIdx, temperature, Pnew, saltconcentration) * waterReferenceDensity(regionIdx);
                Evaluation jouleThomsonCoefficient = -(1.0/Cp) * (1.0 - alpha * temperature)/rho;
                Evaluation deltaEnthalpyPres = -Cp * jouleThomsonCoefficient * deltaP;
                enthalpyPres = enthalpyPresPrev + deltaEnthalpyPres;
                enthalpyPresPrev = enthalpyPres;
            }
        }
        else {
              throw std::runtime_error("Requested Joule-thomson calculation but thermal water density (WATDENT) is not provided");
        }

        Evaluation enthalpy = Cp * (temperature - Tref) + enthalpyPres;

        return enthalpy - pressure/density;
    }
}

template<class Scalar, bool enableBrine>
template<class Evaluation>
Evaluation WaterPvtThermal<Scalar,enableBrine>::
viscosity(unsigned regionIdx,
          const Evaluation& temperature,
          const Evaluation& pressure,
          const Evaluation& saltconcentration) const
{
    const auto& isothermalMu = isothermalPvt_->viscosity(regionIdx, temperature, pressure, saltconcentration);
    if (!enableThermalViscosity())
        return isothermalMu;

    Scalar x = -pvtwViscosibility_[regionIdx]*(viscrefPress_[regionIdx] - pvtwRefPress_[regionIdx]);
    Scalar muRef = pvtwViscosity_[regionIdx]/(1.0 + x + 0.5*x*x);

    // compute the viscosity deviation due to temperature
    const auto& muWatvisct = watvisctCurves_[regionIdx].eval(temperature, true);
    return isothermalMu * muWatvisct/muRef;
}

template<class Scalar, bool enableBrine>
template<class Evaluation>
Evaluation WaterPvtThermal<Scalar,enableBrine>::
inverseFormationVolumeFactor(unsigned regionIdx,
                             const Evaluation& temperature,
                             const Evaluation& pressure,
                             const Evaluation& saltconcentration) const
{
    if (!enableThermalDensity())
        return isothermalPvt_->inverseFormationVolumeFactor(regionIdx, temperature, pressure, saltconcentration);

    Scalar BwRef = pvtwRefB_[regionIdx];
    Scalar TRef = watdentRefTemp_[regionIdx];
    const Evaluation& X = pvtwCompressibility_[regionIdx]*(pressure - pvtwRefPress_[regionIdx]);
    Scalar cT1 = watdentCT1_[regionIdx];
    Scalar cT2 = watdentCT2_[regionIdx];
    const Evaluation& Y = temperature - TRef;

    // this is inconsistent with the density calculation of water in the isothermal
    // case (it misses the quadratic pressure term), but it is the equation given in
    // the documentation.
    return 1.0/(((1 - X)*(1 + cT1*Y + cT2*Y*Y))*BwRef);
}

} // namespace Opm
