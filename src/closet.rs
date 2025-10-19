//
// Copyright (c) 2025 Varazdat Stepanyan.
//
// This file is part of KMC, an open-source engine for kinetic
// monte carlo simulations. KMC is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License
// version 3 as published by the Free Software Foundation.
//
// KMC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with KMC.
// If not, see <https://www.gnu.org/licenses/>.
//

pub trait IsEnv {} // Mark as a collection of env info

pub trait IsObs {} // Mark as a collection of observables

// Mark as a physical state
//      get_obs: makes a weak measurement on the state which returns the observed values
pub trait IsState<D: IsObs + Clone> {
    fn get_obs(&self) -> D;
}

// Represent a decision made by the system
pub enum Decision<S> {
    Skip { dt: f32 },
    Do { dt: f32, dec: S },
}

// Represent the result of a measurement at a specific time
#[derive(Clone)]
pub struct Result<D: IsObs + Clone> {
    pub t: f32,
    pub obs: D,
}

// Mark as a physical system, contais a state as well as a measure of time etc
//      get_obs: makes a weak measurement on the state of the system which returns the observed values
//      suggest: provides a new state suggestion =>
//  =>  decide: processes the suggestionn and decides whether to change the state of the system =>
//  =>  step: processes the decision applying the changes to the state
//      cond: whether the simulation should stop
pub trait IsSystem<D: IsObs + Clone, S: IsState<D> + Clone, E: IsEnv> {
    fn new(e: Option<E>) -> Self;
    fn get(&self) -> Result<D>;
    fn suggest(&self) -> S;
    fn decide(&self, new: S) -> Decision<S>;
    fn step(&mut self, dec: Decision<S>);
    fn cond(&self) -> bool;
}
