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
 * \copydoc Opm::EclEpsTwoPhaseLawPoints
 */
#ifndef OPM_ECL_EPS_GRID_PROPERTIES_HPP
#define OPM_ECL_EPS_GRID_PROPERTIES_HPP

#include "EclEpsConfig.hpp"

#if HAVE_ECL_INPUT
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Deck/DeckRecord.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SgfnTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SgofTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SlgofTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/Sof2Table.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/Sof3Table.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SwfnTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SwofTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>
#endif

#include <opm/material/common/Means.hpp>

#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <memory>
namespace Opm {
/*!
 * \brief Collects all grid properties which are relevant for end point scaling.
 *
 * This class is used for both, the drainage and the imbibition variants of the ECL
 * keywords.
 */

namespace {
template <typename T>
std::vector<T> compressed_copy(const std::vector<T>& input_vector, const std::vector<int>& indexmap, const std::vector<int>& compressedToCartesianElemIdx) {
    std::vector<T> compressed = std::vector<T>(compressedToCartesianElemIdx.size());

    for (std::size_t elem_index = 0; elem_index < compressedToCartesianElemIdx.size(); elem_index++) {
        auto global_index = compressedToCartesianElemIdx[elem_index];
        compressed[elem_index] = input_vector[indexmap[global_index]];
    }

    return compressed;
}



std::vector<double> try_get(const FieldPropsManager& fp, const std::string& keyword, const std::vector<int>& compressedToCartesianElemIdx) {
    if (fp.has<double>(keyword))
        return compressed_copy(fp.get<double>(keyword), fp.indexmap(), compressedToCartesianElemIdx);

    return {};
}

}



class EclEpsGridProperties
{

public:
#if HAVE_ECL_INPUT

    EclEpsGridProperties(const Opm::EclipseState& eclState,
                         bool useImbibition,
                         const std::vector<int>& compressedToCartesianElemIdx)
    {
        std::string kwPrefix = useImbibition?"I":"";

        const auto& fp = eclState.fieldProps();
        const auto& indexmap = fp.indexmap();

        if (useImbibition)
            compressed_satnum = compressed_copy(fp.get<int>("IMBNUM"), indexmap, compressedToCartesianElemIdx);
        else
            compressed_satnum = compressed_copy(fp.get<int>("SATNUM"), indexmap, compressedToCartesianElemIdx);

        this->compressed_swl = try_get( fp, kwPrefix+"SWL", compressedToCartesianElemIdx);
        this->compressed_sgl = try_get( fp, kwPrefix+"SGL", compressedToCartesianElemIdx);
        this->compressed_swcr = try_get( fp, kwPrefix+"SWCR", compressedToCartesianElemIdx);
        this->compressed_sgcr = try_get( fp, kwPrefix+"SGCR", compressedToCartesianElemIdx);
        this->compressed_sowcr = try_get( fp, kwPrefix+"SOWCR", compressedToCartesianElemIdx);
        this->compressed_sogcr = try_get( fp, kwPrefix+"SOGCR", compressedToCartesianElemIdx);
        this->compressed_swu = try_get( fp, kwPrefix+"SWU", compressedToCartesianElemIdx);
        this->compressed_sgu = try_get( fp, kwPrefix+"SGU", compressedToCartesianElemIdx);
        this->compressed_pcw = try_get( fp, kwPrefix+"PCW", compressedToCartesianElemIdx);
        this->compressed_pcg = try_get( fp, kwPrefix+"PCG", compressedToCartesianElemIdx);
        this->compressed_krw = try_get( fp, kwPrefix+"KRW", compressedToCartesianElemIdx);
        this->compressed_kro = try_get( fp, kwPrefix+"KRO", compressedToCartesianElemIdx);
        this->compressed_krg = try_get( fp, kwPrefix+"KRG", compressedToCartesianElemIdx);

        // _may_ be needed to calculate the Leverett capillary pressure scaling factor
        if (fp.has<double>("PORO"))
            this->compressed_poro = compressed_copy(fp.get<double>("PORO"), indexmap, compressedToCartesianElemIdx);

        if (fp.has<double>("PERMX"))
            this->compressed_permx = compressed_copy(fp.get<double>("PERMX"), indexmap, compressedToCartesianElemIdx);
        else
            this->compressed_permx = std::vector<double>(this->compressed_satnum.size());

        if (fp.has<double>("PERMY"))
            this->compressed_permy= compressed_copy(fp.get<double>("PERMY"), indexmap, compressedToCartesianElemIdx);
        else
            this->compressed_permy= this->compressed_permx;

        if (fp.has<double>("PERMZ"))
            this->compressed_permz= compressed_copy(fp.get<double>("PERMZ"), indexmap, compressedToCartesianElemIdx);
        else
            this->compressed_permz= this->compressed_permx;
    }

#endif



    unsigned satRegion(std::size_t active_index) const {
        return this->compressed_satnum[active_index] - 1;
    }

    double permx(std::size_t active_index) const {
        return this->compressed_permx[active_index];
    }

    double permy(std::size_t active_index) const {
        return this->compressed_permy[active_index];
    }

    double permz(std::size_t active_index) const {
        return this->compressed_permz[active_index];
    }

    double poro(std::size_t active_index) const {
        return this->compressed_poro[active_index];
    }

    const double * swl(std::size_t active_index) const {
        return this->satfunc(this->compressed_swl, active_index);
    }

    const double * sgl(std::size_t active_index) const {
        return this->satfunc(this->compressed_sgl, active_index);
    }

    const double * swcr(std::size_t active_index) const {
        return this->satfunc(this->compressed_swcr, active_index);
    }

    const double * sgcr(std::size_t active_index) const {
        return this->satfunc(this->compressed_sgcr, active_index);
    }

    const double * sowcr(std::size_t active_index) const {
        return this->satfunc(this->compressed_sowcr, active_index);
    }

    const double * sogcr(std::size_t active_index) const {
        return this->satfunc(this->compressed_sogcr, active_index);
    }

    const double * swu(std::size_t active_index) const {
        return this->satfunc(this->compressed_swu, active_index);
    }

    const double * sgu(std::size_t active_index) const {
        return this->satfunc(this->compressed_sgu, active_index);
    }

    const double * pcw(std::size_t active_index) const {
        return this->satfunc(this->compressed_pcw, active_index);
    }

    const double * pcg(std::size_t active_index) const {
        return this->satfunc(this->compressed_pcg, active_index);
    }

    const double * krw(std::size_t active_index) const {
        return this->satfunc(this->compressed_krw, active_index);
    }

    const double * krg(std::size_t active_index) const {
        return this->satfunc(this->compressed_krg, active_index);
    }

    const double * kro(std::size_t active_index) const {
        return this->satfunc(this->compressed_kro, active_index);
    }


private:

    const double * satfunc(const std::vector<double>& data, std::size_t active_index) const {
        if (data.empty())
            return nullptr;
        return &(data[active_index]);
    }


    std::vector<int> compressed_satnum;
    std::vector<double> compressed_swl;
    std::vector<double> compressed_sgl;
    std::vector<double> compressed_swcr;
    std::vector<double> compressed_sgcr;
    std::vector<double> compressed_sowcr;
    std::vector<double> compressed_sogcr;
    std::vector<double> compressed_swu;
    std::vector<double> compressed_sgu;
    std::vector<double> compressed_pcw;
    std::vector<double> compressed_pcg;
    std::vector<double> compressed_krw;
    std::vector<double> compressed_kro;
    std::vector<double> compressed_krg;

    std::vector<double> compressed_permx;
    std::vector<double> compressed_permy;
    std::vector<double> compressed_permz;
    std::vector<double> compressed_poro;
};
}
#endif

