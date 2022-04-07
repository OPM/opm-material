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
 */
#include "config.h"
#include <fstream>
#include <opm/output/eclipse/EclipseIO.hpp>

#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Parser/ParseContext.hpp>
#include <opm/input/eclipse/Parser/ErrorGuard.hpp>

#include <opm/input/eclipse/Python/Python.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>
#include <opm/material/fluidstates/SimpleModularFluidState.hpp>

int main(int argc, char** argv) {
    std::string deck_file = argv[1];
    int cellIdx = 0;
    if (argc > 2)
        cellIdx = atol(argv[2]);

    std::string outputname = deck_file.substr(0, deck_file.find(".DATA"));
    Opm::Parser parser;
    Opm::ParseContext parse_context;
    Opm::ErrorGuard error_guard;
    auto python = std::make_shared<Opm::Python>();

    Opm::OpmLog::setupSimpleDefaultLogging();
    Opm::Deck deck = parser.parseFile(deck_file, parse_context, error_guard);
    Opm::EclipseState state(deck);
    const Opm::EclipseState eclState(deck);
    const auto& eclGrid = eclState.getInputGrid();
    size_t n = eclGrid.getCartesianSize();

    enum { numPhases = 3 };
    enum { waterPhaseIdx = 0 };
    enum { oilPhaseIdx = 1 };
    enum { gasPhaseIdx = 2 };
    typedef Opm::ThreePhaseMaterialTraits<double,
                                          /*wettingPhaseIdx=*/waterPhaseIdx,
                                          /*nonWettingPhaseIdx=*/oilPhaseIdx,
                                          /*gasPhaseIdx=*/gasPhaseIdx> MaterialTraits;

    typedef Opm::SimpleModularFluidState<double,
                                         /*numPhases=*/3,
                                         /*numComponents=*/3,
                                         void,
                                         /*storePressure=*/false,
                                         /*storeTemperature=*/false,
                                         /*storeComposition=*/false,
                                         /*storeFugacity=*/false,
                                         /*storeSaturation=*/true,
                                         /*storeDensity=*/false,
                                         /*storeViscosity=*/false,
                                         /*storeEnthalpy=*/false> FluidState;

    typedef Opm::EclMaterialLawManager<MaterialTraits> MaterialLawManager;
    typedef typename MaterialLawManager::MaterialLaw MaterialLaw;

    MaterialLawManager materialLawManager;
    materialLawManager.initFromState(eclState);
    materialLawManager.initParamsForElements(eclState, n);
    double Swco = 0.1;
    int numPoints = 100;
    double dSo = 1.0 / numPoints;
    double So = 0.0;

    std::ofstream kr_og(outputname + "_kr_og.dat");
    std::ofstream kr_ow(outputname + "_kr_ow.dat");
    std::ofstream pc_og(outputname + "_pc_og.dat");
    std::ofstream pc_ow(outputname + "_pc_ow.dat");
    for (int i = 0; i < numPoints; ++i)
    {
        So += dSo;

        // Oil in gas and conate water
        double Sw = Swco;
        double Sg = 1 - Sw - So;
        FluidState fs;
        fs.setSaturation(waterPhaseIdx, Sw);
        fs.setSaturation(oilPhaseIdx, So);
        fs.setSaturation(gasPhaseIdx, Sg);

        double pc[numPhases] = {0.0, 0.0, 0.0};
        MaterialLaw::capillaryPressures(pc,
                                        materialLawManager.materialLawParams(cellIdx),
                                        fs);

        pc_og << So << " " << pc[gasPhaseIdx] << std::endl;

        double kr[numPhases] = {0.0, 0.0, 0.0};
        MaterialLaw::relativePermeabilities(kr,
                                            materialLawManager.materialLawParams(cellIdx),
                                            fs);
        kr_og << So << " " << kr[oilPhaseIdx] << " " << kr[gasPhaseIdx] << std::endl;

        // Oil in water
        Sg = 0.0;
        Sw = 1.0 - So - Sg;
        fs.setSaturation(waterPhaseIdx, Sw);
        fs.setSaturation(oilPhaseIdx, So);
        fs.setSaturation(gasPhaseIdx, Sg);

        MaterialLaw::capillaryPressures(pc,
                                        materialLawManager.materialLawParams(cellIdx),
                                        fs);
        pc_ow << So << " " << pc[waterPhaseIdx] << std::endl;
        MaterialLaw::relativePermeabilities(kr,
                                            materialLawManager.materialLawParams(cellIdx),
                                            fs);
        kr_ow << So << " " << kr[oilPhaseIdx] << " " << kr[waterPhaseIdx] << std::endl;
    }
}
