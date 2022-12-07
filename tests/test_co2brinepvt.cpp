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
/*!
 * \file
 *
 * \brief This is the unit test for the co2 brine PVT model
 *
 * This test requires the presence of opm-common.
 */
#include "config.h"

#if !HAVE_ECL_INPUT
#error "The test for the co2 brine PVT classes requires eclipse input support in opm-common"
#endif

//#include <opm/material/fluidsystems/blackoilpvt/Co2GasPvt.hpp>
//#include <opm/material/fluidsystems/blackoilpvt/BrineCo2Pvt.hpp>

#include <opm/material/fluidsystems/blackoilpvt/GasPvtMultiplexer.hpp>
#include <opm/material/fluidsystems/blackoilpvt/OilPvtMultiplexer.hpp>
#include <opm/material/fluidsystems/blackoilpvt/WaterPvtMultiplexer.hpp>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>

#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Python/Python.hpp>

#include <dune/common/parallel/mpihelper.hh>

// values of strings based on the first SPE1 test case of opm-data.  note that in the
// real world it does not make much sense to specify a fluid phase using more than a
// single keyword, but for a unit test, this saves a lot of boiler-plate code.
static const char* deckString1 =
    "RUNSPEC\n"
    "\n"
    "DIMENS\n"
    "   10 10 3 /\n"
    "\n"
    "TABDIMS\n"
    " * 1 /\n"
    "\n"
    "OIL\n"
    "GAS\n"
    "CO2STOR\n"
    "\n"
    "DISGAS\n"
    "\n"
    "METRIC\n"
    "\n"
    "GRID\n"
    "\n"
    "DX\n"
    "   	300*1000 /\n"
    "DY\n"
    "	300*1000 /\n"
    "DZ\n"
    "	100*20 100*30 100*50 /\n"
    "\n"
    "TOPS\n"
    "	100*1234 /\n"
    "\n"
    "PORO\n"
    "  300*0.15 /\n"
    "PROPS\n"
    "\n";


static const char* deckString2 =
    "RUNSPEC\n"
    "\n"
    "DIMENS\n"
    "   10 10 3 /\n"
    "\n"
    "TABDIMS\n"
    " * 1 /\n"
    "\n"
    "WATER\n"
    "GAS\n"
    "CO2STORE\n"
    "\n"
    "DISGAS\n"
    "\n"
    "METRIC\n"
    "\n"
    "GRID\n"
    "\n"
    "DX\n"
    "   	300*1000 /\n"
    "DY\n"
    "	300*1000 /\n"
    "DZ\n"
    "	100*20 100*30 100*50 /\n"
    "\n"
    "TOPS\n"
    "	100*1234 /\n"
    "\n"
    "PORO\n"
    "  300*0.15 /\n"
    "PROPS\n"
    "\n";

template <class Evaluation, class BrinePvt>
void ensurePvtApiBrine(const BrinePvt& brinePvt)
{
    // we don't want to run this, we just want to make sure that it compiles
    while (0) {
        Evaluation temperature = 273.15 + 20.0;
        Evaluation pressure = 1e5;
        Evaluation saltconcentration = 0.0;
        Evaluation rs = 0.0;
        Evaluation tmp;

        ////
        // Water PVT API
        /////
        tmp = brinePvt.viscosity(/*regionIdx=*/0,
                                 temperature,
                                 pressure,
                                 rs,
                                 saltconcentration);
        tmp = brinePvt.inverseFormationVolumeFactor(/*regionIdx=*/0,
                                                    temperature,
                                                    pressure,
                                                    rs,
                                                    saltconcentration);

        // prevent GCC from producing a "variable assigned but unused" warning
        tmp = 2.0*tmp;
    }
}

template <class Evaluation, class Co2Pvt>
void ensurePvtApiGas(const Co2Pvt& co2Pvt)
{
    // we don't want to run this, we just want to make sure that it compiles
    while (0) {
        Evaluation temperature = 273.15 + 20.0;
        Evaluation pressure = 1e5;
        Evaluation Rv = 0.0;
        Evaluation Rvw = 0.0;
        Evaluation So = 0.5;
        Evaluation maxSo = 1.0;
        Evaluation tmp;

        /////
        // co2 PVT API
        /////
        tmp = co2Pvt.viscosity(/*regionIdx=*/0,
                               temperature,
                               pressure,
                               Rv,
                               Rvw);
        tmp = co2Pvt.inverseFormationVolumeFactor(/*regionIdx=*/0,
                                                  temperature,
                                                  pressure,
                                                  Rv,
                                                  Rvw);
        tmp = co2Pvt.saturatedViscosity(/*regionIdx=*/0,
                                        temperature,
                                        pressure);
        tmp = co2Pvt.saturatedInverseFormationVolumeFactor(/*regionIdx=*/0,
                                                           temperature,
                                                           pressure);
        tmp = co2Pvt.saturationPressure(/*regionIdx=*/0,
                                        temperature,
                                        Rv);
        tmp = co2Pvt.saturatedOilVaporizationFactor(/*regionIdx=*/0,
                                                    temperature,
                                                    pressure);
        tmp = co2Pvt.saturatedOilVaporizationFactor(/*regionIdx=*/0,
                                                    temperature,
                                                    pressure,
                                                    So,
                                                    maxSo);

        // prevent GCC from producing a "variable assigned but unused" warning
        tmp = 2.0*tmp;
    }
}

template <class Evaluation, class BrinePvt>
void ensurePvtApiBrineOil(const BrinePvt& brinePvt)
{
    // we don't want to run this, we just want to make sure that it compiles
    while (0) {
        Evaluation temperature = 273.15 + 20.0;
        Evaluation pressure = 1e5;
        Evaluation Rs = 0.0;
        Evaluation So = 0.5;
        Evaluation maxSo = 1.0;
        Evaluation tmp;

        /////
        // brine PVT API
        /////
        tmp = brinePvt.viscosity(/*regionIdx=*/0,
                               temperature,
                               pressure,
                               Rs);
        tmp = brinePvt.inverseFormationVolumeFactor(/*regionIdx=*/0,
                                                  temperature,
                                                  pressure,
                                                  Rs);
        tmp = brinePvt.saturatedViscosity(/*regionIdx=*/0,
                                        temperature,
                                        pressure);
        tmp = brinePvt.saturatedInverseFormationVolumeFactor(/*regionIdx=*/0,
                                                           temperature,
                                                           pressure);
        tmp = brinePvt.saturationPressure(/*regionIdx=*/0,
                                        temperature,
                                        Rs);
        tmp = brinePvt.saturatedGasDissolutionFactor(/*regionIdx=*/0,
                                                   temperature,
                                                   pressure);
        tmp = brinePvt.saturatedGasDissolutionFactor(/*regionIdx=*/0,
                                                   temperature,
                                                   pressure,
                                                   So,
                                                   maxSo);
    }
}

template <class Scalar>
inline void testAll()
{
    Opm::Parser parser;
    auto python = std::make_shared<Opm::Python>();

    auto deck1 = parser.parseString(deckString1);
    Opm::EclipseState eclState1(deck1);
    Opm::Schedule schedule1(deck1, eclState1, python);

    Opm::GasPvtMultiplexer<Scalar> co2Pvt_oil;
    Opm::OilPvtMultiplexer<Scalar> brinePvt_oil;

    co2Pvt_oil.initFromState(eclState1, schedule1);
    brinePvt_oil.initFromState(eclState1, schedule1);

    typedef Opm::DenseAd::Evaluation<Scalar, 1> FooEval;
    ensurePvtApiGas<Scalar>(co2Pvt_oil);
    ensurePvtApiBrineOil<FooEval>(brinePvt_oil);

    auto deck2 = parser.parseString(deckString2);
    Opm::EclipseState eclState2(deck2);
    Opm::Schedule schedule2(deck2, eclState2, python);

    Opm::GasPvtMultiplexer<Scalar> co2Pvt;
    Opm::WaterPvtMultiplexer<Scalar> brinePvt;

    co2Pvt.initFromState(eclState2, schedule2);
    brinePvt.initFromState(eclState2, schedule2);

    typedef Opm::DenseAd::Evaluation<Scalar, 1> FooEval;
    ensurePvtApiGas<Scalar>(co2Pvt);
    ensurePvtApiBrine<FooEval>(brinePvt);
}

int main(int argc, char **argv)
{
    Dune::MPIHelper::instance(argc, argv);

    testAll<double>();
    testAll<float>();
    return 0;
}
