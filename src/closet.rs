//
// Copyright (c) 2025-2026 Varazdat Stepanyan.
//
// This file is part of MetroFerris, an open-source engine for kinetic
// monte carlo (and beyond) simulations. MetroFerris is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License
// version 3 as published by the Free Software Foundation.
//
// MetroFerris is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with MetroFerris.
// If not, see <https://www.gnu.org/licenses/>.
//

// Mark as a physical system, contais a state as well as a measure of time etc
//      get: makes a weak measurement on the state of the system
//      step: makes a monte carlo step during the simulation
//      cond: while true do another step, once false stop the simulation
//      store_cond: if true store the result of get
pub trait IsModel {
    type Obs;
    fn get(&self) -> Self::Obs;
    fn step(&mut self);
    fn cond(&self) -> bool;
    fn store_cond(&mut self) -> bool;
}
