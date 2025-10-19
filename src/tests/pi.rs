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

#[cfg(test)]
use crate::engine::simulate;
use rand::{self, Rng};

pub fn test() -> bool {
    let res = simulate(&mut System::new(None));
    let r = res.last().expect("Critical error in testing: Pi");
    if (r.obs.0 / r.t - std::f32::consts::PI).abs() < 0.1 {
        return true;
    }
    println!("{}", r.obs.0 / r.t);
    false
}

use super::super::closet::{Decision, IsEnv, IsObs, IsState, IsSystem, Result};

#[derive(Clone)]
struct Observables(f32);

impl IsObs for Observables {}

#[derive(Clone)]
struct State(f32);

impl State {}

impl IsState<Observables> for State {
    fn get_obs(&self) -> Observables {
        Observables(4.0 * self.0)
    }
}

struct Env(f32);
impl IsEnv for Env {}

#[derive(Clone)]
struct System {
    state: State,
    i: f32,
}

impl IsSystem<Observables, State, Env> for System {
    fn new(_: Option<Env>) -> Self {
        System {
            state: State(0.0),
            i: 0.0,
        }
    }

    fn get(&self) -> Result<Observables> {
        Result {
            t: self.i,
            obs: self.state.get_obs(),
        }
    }

    fn suggest(&self) -> State {
        State(self.state.0 + 1.0)
    }

    fn decide(&self, new: State) -> Decision<State> {
        let x = rand::rng().random::<f32>();
        let y = rand::rng().random::<f32>();
        if x.powi(2) + y.powi(2) > 1.0 {
            return Decision::Skip { dt: 1.0 };
        }
        Decision::Do { dt: 1.0, dec: new }
    }

    fn step(&mut self, dec: Decision<State>) {
        match dec {
            Decision::Skip { dt } => self.i += dt,
            Decision::Do { dt, dec } => {
                self.i += dt;
                self.state = dec
            }
        };
    }

    fn cond(&self) -> bool {
        if self.i < 100000.0 {
            return true;
        }
        false
    }
}
