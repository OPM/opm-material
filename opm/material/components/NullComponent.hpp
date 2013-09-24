// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \copydoc Opm::NullComponent
 */
#ifndef OPM_NULL_COMPONENT_HH
#define OPM_NULL_COMPONENT_HH

#include "Component.hpp"

namespace Opm
{
/*!
 * \ingroup Components
 *
 * \brief A component that only throws exceptions.
 *
 * Its main purpose is to make things compile and give runtime errors
 * if it is actually used.
 */
template <class Scalar, int isLiquidV=-1>
class NullComponent : public Component<Scalar, NullComponent<Scalar> >
{
};

} // end namepace

#endif
