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
use kmc_derive::read_var;
use rayon::prelude::*;
use std::io::{self, Write};

pub fn main() {
    let lambda: f32;
    read_var!(lambda, "Lambda");

    let del: f32;
    read_var!(del, "Del");

    let a: f32;
    read_var!(a, "A");

    let b: f32;
    read_var!(b, "B");

    let dt: f32;
    read_var!(dt, "dt");

    let t_max = 10.0 / lambda;

    let env = Env {
        lambda,
        del,
        a,
        b,
        t_max,
    };

    let sys = System::new(Some(env));

    let rep_num: usize;
    read_var!(rep_num, "Repetitions");

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
    println!("{}", res.to_str());
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
        let size = (t_max / dt).ceil() as usize;
        let mut time: Vec<f32> = Vec::new();
        let obs: Vec<Observables> = vec![Observables { prop: 0.0 }; size];
        let am: Vec<u32> = vec![0; size];
        for i in 0..=size {
            time.push((i as f32) * dt);
        }
        Results { time, obs, am }
    }

    fn add(&mut self, r: Result<Observables>) {
        if let Some(i) = helpers::binary_search(r.t, &self.time) {
            self.obs[i].prop += r.obs.prop;
            self.am[i] += 1;
        }
    }

    fn ready(&mut self) {
        self.am.iter().enumerate().rev().for_each(|(i, a)| {
            if *a != 0 {
                self.obs[i].prop /= *a as f32;
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
            s.push_str(&o.prop.to_string());
            s.push(',');
        });
        s.pop();
        s.push('\n');

        s
    }
}

#[derive(Debug, Clone, Copy, Observable)]
struct Observables {
    pub prop: f32,
}

#[derive(Clone, Copy)]
struct State {
    prop: f32,
    eps: i8,
}

impl IsState for State {
    type Obs = Observables;
    fn get_obs(&self) -> Observables {
        Observables { prop: self.prop }
    }
}

#[derive(Clone, Copy)]
struct Env {
    lambda: f32,
    del: f32,
    a: f32,
    b: f32,
    t_max: f32,
}

impl IsEnv for Env {}

#[derive(Clone, Copy)]
struct System {
    state: State,
    lambda: f32,
    del: f32,
    a: f32,
    b: f32,
    t: f32,
    t_max: f32,
}

impl IsSystem for System {
    type State = State;
    type Env = Env;
    fn new(e: Option<Env>) -> Self {
        let env = e.unwrap_or(Env {
            lambda: 1.0,
            del: 0.5,
            a: 2.0,
            b: 0.5,
            t_max: 10.0,
        });
        System {
            state: State { prop: 0.5, eps: 1 },
            lambda: env.lambda,
            del: env.del,
            a: env.a,
            b: env.b,
            t: 0.0,
            t_max: env.t_max,
        }
    }

    fn get(&self) -> Result<Observables> {
        Result {
            t: self.t,
            obs: self.state.get_obs(),
        }
    }

    fn step(&mut self) {
        let u: f32 = -(1.0 - rand::random::<f32>()).ln();

        let dt = u / (2.0 * self.lambda);

        self.ode_solv(dt);

        if rand::random::<f32>() > self.lambda + (self.state.eps as f32) * self.del {
            self.t += dt;
        } else {
            self.t += dt;
            self.state.eps *= -1;
        }
    }

    fn cond(&self) -> bool {
        self.t < self.t_max
    }

    fn store_cond(&mut self) -> bool {
        true
    }
}

impl System {
    fn ode_solv(&mut self, t: f32) {
        let h = t / 10.0;
        for _ in 0..10 {
            self.state.prop += self.ode(self.state.prop + self.ode(self.state.prop) * h / 2.0) * h;
        }
    }

    fn ode(&self, x: f32) -> f32 {
        x * (1.0 - x) * ((self.state.eps as f32) + self.a - self.b * x)
    }
}
