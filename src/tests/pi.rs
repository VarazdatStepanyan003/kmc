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
use kmc_derive::Observable;
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

use super::super::closet::{IsEnv, IsObs, IsState, IsSystem, Result};

#[derive(Clone, Observable)]
struct Observables(f32);

#[derive(Clone)]
struct State(f32);

impl State {}

impl IsState for State {
    type Obs = Observables;
    fn get_obs(&self) -> Observables {
        Observables(4.0 * self.0)
    }
}

struct Env(f32);
impl IsEnv for Env {}

#[derive(Clone)]
struct System {
    state: State,
    i: usize,
}

impl IsSystem for System {
    type State = State;
    type Env = Env;

    fn new(_: Option<Env>) -> Self {
        System {
            state: State(0.0),
            i: 0,
        }
    }

    fn get(&self) -> Result<Observables> {
        Result {
            t: self.i as f32,
            obs: self.state.get_obs(),
        }
    }

    fn step(&mut self) {
        let x = rand::rng().random::<f32>();
        let y = rand::rng().random::<f32>();
        if x.powi(2) + y.powi(2) > 1.0 {
            self.i += 1;
        } else {
            self.state = State(self.state.0 + 1.0);
            self.i += 1;
        }
    }

    fn cond(&self) -> bool {
        if self.i < 100000 {
            return true;
        }
        false
    }

    fn store_cond(&mut self) -> bool {
        true
    }
}
