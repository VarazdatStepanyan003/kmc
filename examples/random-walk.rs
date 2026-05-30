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

use std::fmt::Display;
use std::ops::AddAssign;

use metroferris::prelude::*;
use rand::RngExt;

pub fn main() {
    let e = Env;
    let r = engine::simulate::<Env>(&mut e.create());
    //for (j, rr) in r.iter().enumerate() {
    //    if j < r.iter().len() - 1 {
    //        print!("{} @ t={}, ", rr.obs, rr.t);
    //    } else {
    //        println!("{} @ t={}", rr.obs, rr.t);
    //    }
    //}
    let last = r.last().unwrap();
    println!("Got to distance of {} in time {}", last.obs.norm(), last.t)
}

struct Env;

impl IsEnv for Env {
    type Model = ThreeD;
    fn create(self) -> Self::Model {
        ThreeD {
            state: Vec3(0, 0, 0),
            t: 0.0,
        }
    }
}

#[derive(Clone, Copy, Observable)]
struct Vec3(isize, isize, isize);

impl Display for Vec3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Vec3({}, {}, {})", self.0, self.1, self.2)
    }
}

impl IsState for Vec3 {
    type Obs = Vec3;
    fn get_obs(&self) -> Self::Obs {
        *self
    }
}

impl AddAssign for Vec3 {
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0;
        self.1 += rhs.1;
        self.2 += rhs.2;
    }
}

impl Vec3 {
    fn norm(&self) -> f32 {
        ((self.0.pow(2) + self.1.pow(2) + self.2.pow(2)) as f32).sqrt()
    }
}

struct ThreeD {
    state: Vec3,
    t: f32,
}

impl IsModel for ThreeD {
    type State = Vec3;
    fn get(&self) -> Result<<Self::State as IsState>::Obs> {
        Result {
            t: self.t,
            obs: self.state,
        }
    }
    fn step(&mut self) {
        let mut rng = rand::rng();

        let d: u8 = rng.random_range(0..6);
        let mut del = Vec3(0, 0, 0);
        match d {
            0 => del.0 = 1,
            1 => del.0 = -1,
            2 => del.1 = 1,
            3 => del.1 = -1,
            4 => del.2 = 1,
            5 => del.2 = -1,
            _ => {}
        }

        self.state += del;
        self.t -= rng.random::<f32>().ln();
    }
    fn cond(&self) -> bool {
        self.t < 1000000.0
    }
    fn store_cond(&mut self) -> bool {
        true
    }
}
