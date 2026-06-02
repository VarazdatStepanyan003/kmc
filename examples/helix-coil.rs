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
use std::{
    sync::{Arc, Mutex},
    thread,
};

pub fn main() {
    let (t_max, j, bi, bc, bf) = (3000.0, 1.0, 1.8, 1.5, 1.2);

    let res = engine::simulate(&mut ZB::new(j, bi, bc, bf, t_max));

    for r in res {
        println!(
            "Helicity: {}, N of domains: {}, β: {}",
            r.moms.avg,
            r.moms.doms,
            beta(r.t / t_max, bi, bf)
        );
    }
}

fn beta(t: f32, bi: f32, bf: f32) -> f32 {
    (bf - bi) * t + bi
}

struct Moments {
    avg: f32,
    doms: u8,
}

struct Observables {
    t: f32,
    moms: Moments,
}

struct ZB {
    state: u128,
    j: f32,
    bi: f32,
    bc: f32,
    bf: f32,
    t: f32,
    t_st: f32,
    t_max: f32,
    t_st_d: f32,
    rng: ThreadRng,
}

impl ZB {
    fn new(j: f32, bi: f32, bc: f32, bf: f32, t_max: f32) -> ZB {
        ZB {
            state: u128::MAX,
            t: 0.0,
            j,
            bi,
            bc,
            bf,
            t_max,
            t_st: t_max / 30.0,
            t_st_d: t_max / 30.0,
            rng: rand::rng(),
        }
    }
}

impl IsModel for ZB {
    type Obs = Observables;
    fn get(&self) -> Observables {
        let mut avg: u8 = 0;
        let mut corr: u8 = 0;
        for i in 0..127 {
            avg += ((self.state >> i) & 1) as u8;
            corr += (((self.state >> i) & 1) * ((self.state >> (i + 1)) & 1)) as u8;
        }
        avg += ((self.state >> 127) & 1) as u8;
        corr += (((self.state >> 127) & 1) * (self.state & 1)) as u8;

        Observables {
            t: self.t,
            moms: Moments {
                avg: (avg as f32) / 128.0,
                doms: 4 * avg - 4 * corr,
            },
        }
    }

    fn step(&mut self) {
        let u: Vec<f32> = (0..128)
            .map(|_| -((1.0 - self.rng.random::<f32>()).ln()))
            .collect();

        let res = Arc::new(Mutex::new([0.0; 128]));
        let mut handles = vec![];

        for k in 0..128 {
            let resloc = Arc::clone(&res);
            let u = u[k];
            let (sm, s, sp): (bool, bool, bool);
            s = self.state & (1 << k) > 0;
            if k == 0 {
                sp = (self.state & (1 << (k + 1))) > 0;
                sm = (self.state & (1 << 127)) > 0;
            } else {
                sm = (self.state & (1 << (k - 1))) > 0;
                if k == 127 {
                    sp = (self.state & 1) > 0;
                } else {
                    sp = (self.state & (1 << (k + 1))) > 0;
                }
            }
            let (t, t_max, bi, bc, bf, j) = (self.t, self.t_max, self.bi, self.bc, self.bf, self.j);
            let handle = thread::spawn(move || {
                let res = integrate((sm, sp, s), (t, t_max), (bi, bc, bf), j, u);
                match res {
                    Some(del) => resloc.lock().unwrap()[k] = del,
                    None => resloc.lock().unwrap()[k] = f32::MAX,
                }
            });

            handles.push(handle);
        }

        for handle in handles {
            handle.join().unwrap();
        }

        let res = res.lock().unwrap();

        let (res, delta) = res
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap();

        self.t += delta;
        self.state ^= 1 << res;
    }

    fn cond(&self) -> bool {
        self.t < self.t_max
    }

    fn store_cond(&mut self) -> bool {
        if self.t > self.t_st {
            self.t_st += self.t_st_d;
            return true;
        }
        false
    }
}

fn coefs((sm, s, sp): (bool, bool, bool), j: f32) -> (f32, f32) {
    let mut rj;
    if sm && sp {
        rj = 4.0 * j;
    } else if !(sm || sp) {
        rj = -4.0 * j;
    } else {
        rj = 0.0;
    }

    let mut rh = 4.0 * j;

    if s {
        rj = -rj;
        rh = -rh;
    }

    (rj, rh)
}

fn integrate(
    (sm, s, sp): (bool, bool, bool),
    (t_0, t_max): (f32, f32),
    (bi, bc, bf): (f32, f32, f32),
    j: f32,
    u: f32,
) -> Option<f32> {
    let (rj, rh) = coefs((sm, s, sp), j);

    let mut t = 0.0;
    let mut be = beta(t_0 / t_max, bi, bf);
    let mut de = rj * be - rh * (bc - be);

    let mut r = helpers::sigmoid(de);

    let dt = u / r / 10.0;

    if dt > t_max / 100.0 {
        return None;
    }

    let mut v = 0.0;
    while v < u {
        if r < 0.0001 {
            return None;
        }

        be = beta((t_0 + t + dt) / t_max, bi, bf);
        de = rj * be - rh * (bc - be);

        let rr = helpers::sigmoid(de);

        v += dt * (rr + r) / 2.0;
        t += dt;
        r = rr;
    }
    if t > t_max {
        Some(t_max)
    } else {
        Some(t)
    }
}
