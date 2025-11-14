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

use kmc::engine;
use rand::{self, Rng};
use rayon::prelude::*;
use std::io::{self, Write};

pub fn main() {
    let mut tmp = String::new();
    print!("max time: ");
    io::stdout().flush().unwrap();
    io::stdin()
        .read_line(&mut tmp)
        .expect("unable to read line");
    let t_max: f32 = tmp.trim().parse().expect("not a valid max time");

    let mut tmp = String::new();
    print!("dt: ");
    io::stdout().flush().unwrap();
    io::stdin()
        .read_line(&mut tmp)
        .expect("unable to read line");
    let dt: f32 = tmp.trim().parse().expect("not a valid dt");

    let mut tmp = String::new();
    print!("beta*J: ");
    io::stdout().flush().unwrap();
    io::stdin()
        .read_line(&mut tmp)
        .expect("unable to read line");
    let bj = tmp.trim().parse().expect("not a valid bj");

    let mut tmp = String::new();
    print!("beta*h: ");
    io::stdout().flush().unwrap();
    io::stdin()
        .read_line(&mut tmp)
        .expect("unable to read line");
    let bh = tmp.trim().parse().expect("not a valid bh");

    let env = Env { bj, bh, t_max };
    let sys = Ising::new(Some(env));

    let mut tmp = String::new();
    print!("repetitions: ");
    io::stdout().flush().unwrap();
    io::stdin()
        .read_line(&mut tmp)
        .expect("unable to read line");
    let rep_num: usize = tmp
        .trim()
        .parse()
        .expect("not a valid amount of repetitions");

    let mut r: Vec<Vec<Result<Observables>>> = vec![Vec::new(); rep_num];

    r.par_iter_mut().for_each(|rt| {
        *rt = engine::simulate(&mut sys.clone());
    });

    let mut res = Results::new(dt, t_max);

    r.iter().for_each(|tmp| {
        tmp.iter().for_each(|rs| {
            res.add(rs.clone());
        });
    });

    res.ready();

    res.to_str();
    //println!("{}", res.to_str());
    //std::fs::write("res.txt", res.to_str().as_str()).expect("did not write");
}

use kmc::closet::{IsEnv, IsObs, IsState, IsSystem, Result};
use kmc::helpers;
use kmc_derive::Observable;

struct Results {
    time: Vec<f32>,
    obs: Vec<Observables>,
    am: Vec<u32>,
}

impl Results {
    fn new(dt: f32, t_max: f32) -> Results {
        let mut time: Vec<f32> = Vec::new();
        let mut obs: Vec<Observables> = Vec::new();
        let mut am: Vec<u32> = Vec::new();
        for i in 0..=((t_max / dt).ceil() as usize) {
            time.push((i as f32) * dt);
            obs.push(Observables {
                avg: 0.0,
                corr: 0.0,
            });
            am.push(0);
        }
        Results { time, obs, am }
    }

    fn add(&mut self, r: Result<Observables>) {
        if let Some(i) = helpers::binary_search(r.t, &self.time) {
            self.obs[i].avg += r.obs.avg;
            self.obs[i].corr += r.obs.corr;
            self.am[i] += 1;
        }
    }

    fn ready(&mut self) {
        self.am.iter().enumerate().rev().for_each(|(i, a)| {
            if *a != 0 {
                self.obs[i].avg /= *a as f32;
                self.obs[i].corr /= *a as f32;
            } else {
                self.time.remove(i);
                self.obs.remove(i);
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

        self.obs.iter().for_each(|o| {
            s.push_str(&o.avg.to_string());
            s.push(',');
        });
        s.pop();
        s.push('\n');

        self.obs.iter().for_each(|o| {
            s.push_str(&o.corr.to_string());
            s.push(',');
        });
        s.pop();
        s.push('\n');

        s
    }
}

#[derive(Debug, Clone, Observable)]
struct Observables {
    pub avg: f32,
    pub corr: f32,
}

#[derive(Clone)]
struct State {
    state: u128,
}

impl State {
    fn energy(&self, bj: f32, bh: f32) -> f32 {
        let ob = self.get_obs();
        -bj * ob.corr - bh * ob.avg
    }
}

impl IsState for State {
    type Obs = Observables;
    fn get_obs(&self) -> Observables {
        let mut avg: f32 = 0.0;
        let mut corr: f32 = 0.0;
        for i in 0..127 {
            avg += ((self.state >> i) & 1) as f32;
            corr += (((self.state >> i) & 1) * ((self.state >> (i + 1)) & 1)) as f32;
        }
        avg += ((self.state >> 127) & 1) as f32;
        corr += (((self.state >> 127) & 1) * (self.state & 1)) as f32;
        Observables {
            avg: 2.0 * avg - 128.0,
            corr: 4.0 * corr - 4.0 * avg + 128.0,
        }
    }
}

struct Env {
    bj: f32,
    bh: f32,
    t_max: f32,
}

impl IsEnv for Env {}

#[derive(Clone)]
struct Ising {
    state: State,
    bj: f32,
    bh: f32,
    t: f32,
    t_max: f32,
    t_store: f32,
}

impl IsSystem for Ising {
    type State = State;
    type Env = Env;
    fn new(e: Option<Env>) -> Self {
        let env = e.unwrap_or(Env {
            bj: 1.0,
            bh: 0.5,
            t_max: 10.0,
        });
        Ising {
            state: State { state: 0 },
            bj: env.bj,
            bh: env.bh,
            t: 0.0,
            t_max: env.t_max,
            t_store: 0.0,
        }
    }

    fn get(&self) -> Result<Observables> {
        Result {
            t: self.t,
            obs: self.state.get_obs(),
        }
    }

    fn step(&mut self) {
        let new = State {
            state: self.state.state ^ (1 << rand::rng().random_range(0..128)),
        };

        let mut rng = rand::rng();
        let r =
            helpers::sigmoid(self.state.energy(self.bj, self.bh) - new.energy(self.bj, self.bh));
        let u: f32 = 1.0 - rng.random::<f32>();
        let dt = 1.0 / u / 128.0;
        if rng.random::<f32>() > r {
            self.t += dt;
        } else {
            self.t += dt;
            self.state = new;
        }
    }

    fn cond(&self) -> bool {
        self.t < self.t_max
    }

    fn store_cond(&mut self) -> bool {
        true
    }
}
