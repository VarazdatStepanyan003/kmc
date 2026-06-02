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

#[cfg(test)]
use crate::prelude::*;
use rand::{self, rngs::ThreadRng, RngExt};

pub fn test() -> bool {
    let res = engine::simulate(&mut Pi {
        hit: 0,
        num: 0,
        rng: rand::rng(),
    });
    let res: Value = *res.last().expect("Critical error in testing: Pi");
    if (res.0 - std::f32::consts::PI).abs() < 0.1 {
        return true;
    }
    println!("{}", res.0);
    false
}

#[derive(Clone, Copy)]
struct Value(f32);

struct Pi {
    hit: usize,
    num: usize,
    rng: ThreadRng,
}

impl IsModel for Pi {
    type Obs = Value;

    fn get(&self) -> Self::Obs {
        Value(4.0 * self.hit as f32 / (self.num as f32))
    }

    fn step(&mut self) {
        let x = self.rng.random::<f32>();
        let y = self.rng.random::<f32>();
        self.num += 1;
        if x.powi(2) + y.powi(2) <= 1.0 {
            self.hit += 1;
        }
    }

    fn cond(&self) -> bool {
        if self.num < 100000 {
            return true;
        }
        false
    }

    fn store_cond(&mut self) -> bool {
        !self.cond()
    }
}
