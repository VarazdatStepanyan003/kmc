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

pub trait IsEnv {} // Mark as a collection of env info

pub trait IsObs: Copy {} // Mark as a collection of observables

// Mark as a physical state
//      get_obs: makes a weak measurement on the state which returns the observed values
pub trait IsState: Copy {
    type Obs: IsObs;
    fn get_obs(&self) -> Self::Obs;
}

// Represent the result of a measurement at a specific time
#[derive(Clone, Copy)]
pub struct Result<D: IsObs> {
    pub t: f32,
    pub obs: D,
}

// Mark as a physical system, contais a state as well as a measure of time etc
//      get_obs: makes a weak measurement on the state of the system which returns the observed values
//      suggest: provides a new state suggestion =>
//  =>  decide: processes the suggestionn and decides whether to change the state of the system =>
//  =>  step: processes the decision applying the changes to the state
//      cond: whether the simulation should stop
pub trait IsSystem {
    type State: IsState;
    type Env: IsEnv;

    fn new(e: Option<Self::Env>) -> Self;
    fn get(&self) -> Result<<Self::State as IsState>::Obs>;
    fn step(&mut self);
    fn cond(&self) -> bool;
    fn store_cond(&mut self) -> bool;
}
