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

use metroferris::prelude::*;
use rand::{rngs::ThreadRng, RngExt};
use std::io::{self, Write};
use std::ops::{AddAssign, DivAssign};
use std::sync::{Arc, Mutex};
use std::thread;

pub fn main() {
    let t_max: f32;
    let dt: f32;
    let bj: f32;
    let bh: f32;

    read_var!(t_max, "Max Time");
    read_var!(dt, "δt");
    read_var!(bj, "βj");
    read_var!(bh, "βh");

    let e = Env { bj, bh, t_max };

    let rep_num: usize;
    read_var!(rep_num, "Repetitions");

    let res = Arc::new(Mutex::new(Res::new(dt, t_max)));
    let mut handles = vec![];

    for _ in 0..rep_num {
        let resloc = Arc::clone(&res);
        let eclone = e;
        let handle = thread::spawn(move || {
            for v in engine::simulate(&mut eclone.create()) {
                resloc.lock().unwrap().add(v);
            }
        });

        handles.push(handle);
    }

    for handle in handles {
        handle.join().unwrap();
    }

    let mut res = res.lock().unwrap();

    res.ready();
    res.to_str();
    println!("{}", res.to_str());
    //std::fs::write("res.txt", res.to_str().as_str()).expect("did not write");
}

struct Res {
    time: Vec<f32>,
    moms: Vec<Moments>,
    am: Vec<u32>,
}

impl Res {
    fn new(dt: f32, t_max: f32) -> Res {
        let size = (t_max / dt).ceil() as usize;
        let mut time: Vec<f32> = Vec::new();
        let moms: Vec<Moments> = vec![Moments::zero(); size];
        let am: Vec<u32> = vec![0; size];
        for i in 0..=size {
            time.push((i as f32) * dt);
        }
        Res { time, moms, am }
    }

    fn add(&mut self, r: Observables) {
        if let Some(i) = helpers::binary_search(r.t, &self.time) {
            self.moms[i] += r.moms;
            self.am[i] += 1;
        }
    }

    fn ready(&mut self) {
        self.am.iter().enumerate().rev().for_each(|(i, a)| {
            if *a != 0 {
                self.moms[i] /= *a as f32;
            } else {
                self.time.remove(i);
                self.moms.remove(i);
            }
        });
    }

    fn to_str(&self) -> String {
        let mut s = String::new();
        self.time.iter().for_each(|t| {
            s.push_str(&t.to_string());
            s.push(',');
        });
        s.pop();
        s.push('\n');

        self.moms.iter().for_each(|o| {
            s.push_str(&o.avg.to_string());
            s.push(',');
        });
        s.pop();
        s.push('\n');

        self.moms.iter().for_each(|o| {
            s.push_str(&o.corr.to_string());
            s.push(',');
        });
        s.pop();
        s.push('\n');

        s
    }
}

#[derive(Clone, Copy)]
struct Moments {
    avg: f32,
    corr: f32,
}

impl Moments {
    fn zero() -> Moments {
        Moments {
            avg: 0.0,
            corr: 0.0,
        }
    }
}

impl AddAssign for Moments {
    fn add_assign(&mut self, rhs: Self) {
        self.avg += rhs.avg;
        self.corr += rhs.corr;
    }
}

impl DivAssign<f32> for Moments {
    fn div_assign(&mut self, rhs: f32) {
        self.avg /= rhs;
        self.corr /= rhs;
    }
}

struct Observables {
    t: f32,
    moms: Moments,
}

#[derive(Clone, Copy)]
struct Env {
    bj: f32,
    bh: f32,
    t_max: f32,
}

impl Env {
    fn create(self) -> Ising {
        Ising {
            state: 0,
            bj: self.bj,
            bh: self.bh,
            t: 0.0,
            t_max: self.t_max,
            rng: rand::rng(),
        }
    }
}

struct Ising {
    state: u128,
    bj: f32,
    bh: f32,
    t: f32,
    t_max: f32,
    rng: ThreadRng,
}

impl Ising {
    fn energy_diff(&self, j: usize) -> f32 {
        let (l, c, r): (bool, bool, bool);
        if j == 0 {
            (l, c, r) = (
                (self.state >> 127) & 1 > 0,
                self.state & 1 > 0,
                (self.state >> 1) & 1 > 0,
            );
        } else if j == 127 {
            (l, c, r) = (
                (self.state >> 126) & 1 > 0,
                (self.state >> 127) & 1 > 0,
                self.state & 1 > 0,
            );
        } else {
            (l, c, r) = (
                (self.state >> (j - 1)) & 1 > 0,
                (self.state >> j) & 1 > 0,
                (self.state >> (j + 1)) & 1 > 0,
            );
        }
        let mut res = 2.0 * self.bh;
        if l && r {
            res += 4.0 * self.bj
        } else if !(l || r) {
            res -= 4.0 * self.bj
        }
        if !c {
            res = -res;
        }
        res
    }
}

impl IsModel for Ising {
    type Obs = Observables;
    fn get(&self) -> Observables {
        let mut avg: f32 = 0.0;
        let mut corr: f32 = 0.0;
        for i in 0..127 {
            avg += ((self.state >> i) & 1) as f32;
            corr += (((self.state >> i) & 1) * ((self.state >> (i + 1)) & 1)) as f32;
        }
        avg += ((self.state >> 127) & 1) as f32;
        corr += (((self.state >> 127) & 1) * (self.state & 1)) as f32;
        Observables {
            moms: Moments {
                avg: 2.0 * avg - 128.0,
                corr: 4.0 * corr - 4.0 * avg + 128.0,
            },
            t: self.t,
        }
    }

    fn step(&mut self) {
        let j = self.rng.random_range(0..128);
        let r = helpers::sigmoid(-self.energy_diff(j));
        let u: f32 = -(1.0 - self.rng.random::<f32>()).ln();
        let dt = u / 128.0;
        self.t += dt;
        if self.rng.random::<f32>() < r {
            self.state ^= 1 << j;
        }
    }

    fn cond(&self) -> bool {
        self.t < self.t_max
    }

    fn store_cond(&mut self) -> bool {
        true
    }
}
