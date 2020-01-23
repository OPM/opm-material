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
 * \copydoc Opm::EclThermalLawManager
 */
#if ! HAVE_ECL_INPUT
#error "Eclipse input support in opm-common is required to use the ECL thermal law manager!"
#endif

#ifndef OPM_ECL_THERMAL_LAW_MANAGER_HPP
#define OPM_ECL_THERMAL_LAW_MANAGER_HPP

#include "EclSolidEnergyLawMultiplexer.hpp"
#include "EclSolidEnergyLawMultiplexerParams.hpp"

#include "EclThermalConductionLawMultiplexer.hpp"
#include "EclThermalConductionLawMultiplexerParams.hpp"

#include <opm/material/common/Exceptions.hpp>

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

namespace Opm {

/*!
 * \ingroup fluidmatrixinteractions
 *
 * \brief Provides an simple way to create and manage the thermal law objects
 *        for a complete ECL deck.
 */
template <class Scalar, class FluidSystem>
class EclThermalLawManager
{
public:
    typedef EclSolidEnergyLawMultiplexer<Scalar, FluidSystem> SolidEnergyLaw;
    typedef typename SolidEnergyLaw::Params SolidEnergyLawParams;
    typedef typename SolidEnergyLawParams::HeatcrLawParams HeatcrLawParams;
    typedef typename SolidEnergyLawParams::SpecrockLawParams SpecrockLawParams;

    typedef EclThermalConductionLawMultiplexer<Scalar, FluidSystem> ThermalConductionLaw;
    typedef typename ThermalConductionLaw::Params ThermalConductionLawParams;

    EclThermalLawManager()
    {
        solidEnergyApproach_ = SolidEnergyLawParams::undefinedApproach;
        thermalConductivityApproach_ = ThermalConductionLawParams::undefinedApproach;
    }

    void initFromDeck(const Opm::Deck& deck,
                      const Opm::EclipseState& eclState,
                      const std::vector<int>& compressedToCartesianElemIdx)
    {
        const auto& fp = eclState.fieldProps();
        bool has_heatcr = fp.has<double>("HEATCR");
        bool has_thconr = fp.has<double>("THCONR");
        bool has_thc = fp.has<double>("THCROCK") || fp.has<double>("THCOIL") || fp.has<double>("THCGAS") || fp.has<double>("THCWATER");

        if (has_heatcr)
            initHeatcr_(deck, eclState, compressedToCartesianElemIdx);
        else if (deck.hasKeyword("SPECROCK"))
            initSpecrock_(deck, eclState, compressedToCartesianElemIdx);
        else
            initNullRockEnergy_(deck, eclState, compressedToCartesianElemIdx);


        if (has_thconr)
            initThconr_(deck, eclState, compressedToCartesianElemIdx);
        else if (has_thc)
            initThc_(deck, eclState, compressedToCartesianElemIdx);
        else
            initNullCond_(deck, eclState, compressedToCartesianElemIdx);
    }

    const SolidEnergyLawParams& solidEnergyLawParams(unsigned elemIdx) const
    {
        switch (solidEnergyApproach_) {
        case SolidEnergyLawParams::heatcrApproach:
            assert(0 <= elemIdx && elemIdx <  solidEnergyLawParams_.size());
            return solidEnergyLawParams_[elemIdx];

        case SolidEnergyLawParams::specrockApproach:
        {
            assert(0 <= elemIdx && elemIdx <  elemToSatnumIdx_.size());
            unsigned satnumIdx = elemToSatnumIdx_[elemIdx];
            assert(0 <= satnumIdx && satnumIdx <  solidEnergyLawParams_.size());
            return solidEnergyLawParams_[satnumIdx];
        }

        case SolidEnergyLawParams::nullApproach:
            return solidEnergyLawParams_[0];

        default:
            throw std::runtime_error("Attempting to retrieve solid energy storage parameters "
                                     "without a known approach being defined by the deck.");
        }
    }

    const ThermalConductionLawParams& thermalConductionLawParams(unsigned elemIdx) const
    {
        switch (thermalConductivityApproach_) {
        case ThermalConductionLawParams::thconrApproach:
        case ThermalConductionLawParams::thcApproach:
            assert(0 <= elemIdx && elemIdx <  thermalConductionLawParams_.size());
            return thermalConductionLawParams_[elemIdx];

        case ThermalConductionLawParams::nullApproach:
            return thermalConductionLawParams_[0];

        default:
            throw std::runtime_error("Attempting to retrieve thermal conduction parameters without "
                                     "a known approach being defined by the deck.");
        }
    }

private:
    /*!
     * \brief Initialize the parameters for the solid energy law using using HEATCR and friends.
     */
    void initHeatcr_(const Opm::Deck& deck OPM_UNUSED,
                     const Opm::EclipseState& eclState,
                     const std::vector<int>& compressedToCartesianElemIdx)
    {
        solidEnergyApproach_ = SolidEnergyLawParams::heatcrApproach;
        // actually the value of the reference temperature does not matter for energy
        // conservation. We set it anyway to faciliate comparisons with ECL
        HeatcrLawParams::setReferenceTemperature(FluidSystem::surfaceTemperature);

        const auto& fp = eclState.fieldProps();
        const auto& indexmap = fp.indexmap();
        const std::vector<double>& heatcrData  = fp.get<double>("HEATCR");
        const std::vector<double>& heatcrtData = fp.get<double>("HEATCRT");
        unsigned numElems = compressedToCartesianElemIdx.size();
        solidEnergyLawParams_.resize(numElems);
        for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            auto& elemParam = solidEnergyLawParams_[elemIdx];
            elemParam.setSolidEnergyApproach(SolidEnergyLawParams::heatcrApproach);
            auto& heatcrElemParams = elemParam.template getRealParams<SolidEnergyLawParams::heatcrApproach>();
            const auto& global_index = compressedToCartesianElemIdx[elemIdx];
            const auto& input_index = indexmap[global_index];

            heatcrElemParams.setReferenceRockHeatCapacity(heatcrData[input_index]);
            heatcrElemParams.setDRockHeatCapacity_dT(heatcrtData[input_index]);
            heatcrElemParams.finalize();
            elemParam.finalize();
        }
    }

    /*!
     * \brief Initialize the parameters for the solid energy law using using SPECROCK and friends.
     */
    void initSpecrock_(const Opm::Deck& deck OPM_UNUSED,
                       const Opm::EclipseState& eclState,
                       const std::vector<int>& compressedToCartesianElemIdx)
    {
        solidEnergyApproach_ = SolidEnergyLawParams::specrockApproach;

        // initialize the element index -> SATNUM index mapping
        const auto& fp = eclState.fieldProps();
        const auto& indexmap = fp.indexmap();
        const std::vector<int>& satnumData = fp.get<int>("SATNUM");
        elemToSatnumIdx_.resize(compressedToCartesianElemIdx.size());
        for (unsigned elemIdx = 0; elemIdx < compressedToCartesianElemIdx.size(); ++ elemIdx) {
            unsigned cartesianElemIdx = compressedToCartesianElemIdx[elemIdx];

            // satnumData contains Fortran-style indices, i.e., they start with 1 instead
            // of 0!
            elemToSatnumIdx_[elemIdx] = satnumData[indexmap[cartesianElemIdx]] - 1;
        }
        // internalize the SPECROCK table
        unsigned numSatRegions = eclState.runspec().tabdims().getNumSatTables();
        const auto& tableManager = eclState.getTableManager();
        solidEnergyLawParams_.resize(numSatRegions);
        for (unsigned satnumIdx = 0; satnumIdx < numSatRegions; ++satnumIdx) {
            const auto& specrockTable = tableManager.getSpecrockTables()[satnumIdx];

            auto& multiplexerParams = solidEnergyLawParams_[satnumIdx];

            multiplexerParams.setSolidEnergyApproach(SolidEnergyLawParams::specrockApproach);

            auto& specrockParams = multiplexerParams.template getRealParams<SolidEnergyLawParams::specrockApproach>();
            const auto& temperatureColumn = specrockTable.getColumn("TEMPERATURE");
            const auto& cvRockColumn = specrockTable.getColumn("CV_ROCK");
            specrockParams.setHeatCapacities(temperatureColumn, cvRockColumn);
            specrockParams.finalize();

            multiplexerParams.finalize();
        }
    }

    /*!
     * \brief Specify the solid energy law by setting heat capacity of rock to 0
     */
    void initNullRockEnergy_(const Opm::Deck& deck OPM_UNUSED,
                             const Opm::EclipseState& eclState OPM_UNUSED,
                             const std::vector<int>& compressedToCartesianElemIdx OPM_UNUSED)
    {
        solidEnergyApproach_ = SolidEnergyLawParams::nullApproach;

        solidEnergyLawParams_.resize(1);
        solidEnergyLawParams_[0].finalize();
    }

    /*!
     * \brief Initialize the parameters for the thermal conduction law using THCONR and friends.
     */
    void initThconr_(const Opm::Deck& deck OPM_UNUSED,
                     const Opm::EclipseState& eclState,
                     const std::vector<int>& compressedToCartesianElemIdx)
    {
        thermalConductivityApproach_ = ThermalConductionLawParams::thconrApproach;

        const auto& fp = eclState.fieldProps();
        const auto& indexmap = fp.indexmap();
        std::vector<double> thconrData(fp.active_size(), 0);
        std::vector<double> thconsfData(fp.active_size(), 0);
        if (fp.has<double>("THCONR"))
            thconrData  = fp.get<double>("THCONR");

        if (fp.has<double>("THCONSF"))
            thconsfData = fp.get<double>("THCONSF");

        unsigned numElems = compressedToCartesianElemIdx.size();
        thermalConductionLawParams_.resize(numElems);
        for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            auto& elemParams = thermalConductionLawParams_[elemIdx];
            elemParams.setThermalConductionApproach(ThermalConductionLawParams::thconrApproach);
            auto& thconrElemParams = elemParams.template getRealParams<ThermalConductionLawParams::thconrApproach>();

            int cartElemIdx = compressedToCartesianElemIdx[elemIdx];
            const auto& input_index = indexmap[cartElemIdx];
            thconrElemParams.setReferenceTotalThermalConductivity(thconrData[input_index]);
            thconrElemParams.setDTotalThermalConductivity_dSg(thconsfData[input_index]);

            thconrElemParams.finalize();
            elemParams.finalize();
        }
    }

    /*!
     * \brief Initialize the parameters for the thermal conduction law using THCROCK and friends.
     */
    void initThc_(const Opm::Deck& deck OPM_UNUSED,
                  const Opm::EclipseState& eclState,
                  const std::vector<int>& compressedToCartesianElemIdx)
    {
        thermalConductivityApproach_ = ThermalConductionLawParams::thcApproach;

        const auto& fp = eclState.fieldProps();
        const auto& indexmap = fp.indexmap();
        std::vector<double> thcrockData(fp.active_size(),0);
        std::vector<double> thcoilData(fp.active_size(),0);
        std::vector<double> thcgasData(fp.active_size(),0);
        std::vector<double> thcwaterData = fp.get<double>("THCWATER");

        if (fp.has<double>("THCROCK"))
            thcrockData = fp.get<double>("THCROCK");

        if (fp.has<double>("THCOIL"))
            thcoilData = fp.get<double>("THCOIL");

        if (fp.has<double>("THCGAS"))
            thcgasData = fp.get<double>("THCGAS");

        if (fp.has<double>("THCWATER"))
            thcwaterData = fp.get<double>("THCWATER");

        const std::vector<double>& poroData = fp.get<double>("PORO");

        unsigned numElems = compressedToCartesianElemIdx.size();
        thermalConductionLawParams_.resize(numElems);
        for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            auto& elemParams = thermalConductionLawParams_[elemIdx];
            elemParams.setThermalConductionApproach(ThermalConductionLawParams::thcApproach);
            auto& thcElemParams = elemParams.template getRealParams<ThermalConductionLawParams::thcApproach>();

            int cartElemIdx = compressedToCartesianElemIdx[elemIdx];
            const auto& input_index = indexmap[cartElemIdx];
            thcElemParams.setPorosity(poroData[input_index]);
            thcElemParams.setThcrock(thcrockData[input_index]);
            thcElemParams.setThcoil(thcoilData[input_index]);
            thcElemParams.setThcgas(thcgasData[input_index]);
            thcElemParams.setThcwater(thcwaterData[input_index]);

            thcElemParams.finalize();
            elemParams.finalize();
        }
    }

    /*!
     * \brief Disable thermal conductivity
     */
    void initNullCond_(const Opm::Deck& deck OPM_UNUSED,
                       const Opm::EclipseState& eclState OPM_UNUSED,
                       const std::vector<int>& compressedToCartesianElemIdx OPM_UNUSED)
    {
        thermalConductivityApproach_ = ThermalConductionLawParams::nullApproach;

        thermalConductionLawParams_.resize(1);
        thermalConductionLawParams_[0].finalize();
    }

private:
    typename ThermalConductionLawParams::ThermalConductionApproach thermalConductivityApproach_;
    typename SolidEnergyLawParams::SolidEnergyApproach solidEnergyApproach_;

    std::vector<unsigned> elemToSatnumIdx_;

    std::vector<SolidEnergyLawParams> solidEnergyLawParams_;
    std::vector<ThermalConductionLawParams> thermalConductionLawParams_;
};
} // namespace Opm

#endif
