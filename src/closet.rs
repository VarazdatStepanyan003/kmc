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

// Mark as a collection of env info
//      used to initialize the system and passed as a generic to the simulation
//      to retrieve all the types for the simulation
pub trait IsEnv {
    type Model: IsModel;
    fn create(self) -> Self::Model;
}

pub trait IsObs: Copy {} // Mark as a collection of observables

// Mark as a physical state
//      get_obs: makes a weak measurement on the state which returns the observed values
pub trait IsState: Copy {
    type Obs: IsObs;
    fn get_obs(&self) -> Self::Obs;
}

// Represent the result of a measurement at a specific time
#[derive(Clone, Copy, Debug)]
pub struct Result<D: IsObs> {
    pub t: f32,
    pub obs: D,
}

// Mark as a physical system, contais a state as well as a measure of time etc
//      get: makes a weak measurement on the state of the system and returns the Result
//      step: makes a monte carlo step during the simulation
//      cond: if true do another step, if false stop the simulation
//      store_cond: if true store the result of get
pub trait IsModel {
    type State: IsState;
    fn get(&self) -> Result<<Self::State as IsState>::Obs>;
    fn step(&mut self);
    fn cond(&self) -> bool;
    fn store_cond(&mut self) -> bool;
}
