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
use std::io::{self, Write};
use std::sync::{Arc, Mutex};
use std::thread;

pub fn main() {
    let lambda: f32;
    let del: f32;
    let a: f32;
    let b: f32;
    let dt: f32;

    read_var!(lambda, "Λ");
    read_var!(del, "Δ");
    read_var!(a, "A");
    read_var!(b, "B");
    read_var!(dt, "δt");

    let t_max = 10.0 / lambda;

    let sys = Env {
        lambda,
        del,
        a,
        b,
        t_max,
    }
    .create();

    let rep_num: usize;
    read_var!(rep_num, "Repetitions");

    let res = Arc::new(Mutex::new(Results::new(dt, t_max)));
    let mut handles = vec![];

    for _ in 0..rep_num {
        let resloc = Arc::clone(&res);
        let mut sysclone = sys;
        let handle = thread::spawn(move || {
            for v in engine::simulate::<Env>(&mut sysclone) {
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

impl IsEnv for Env {
    type Model = Stochastic;
    fn create(self) -> Self::Model {
        Stochastic {
            state: State { prop: 0.5, eps: 1 },
            lambda: self.lambda,
            del: self.del,
            a: self.a,
            b: self.b,
            t: 0.0,
            t_max: self.t_max,
        }
    }
}

#[derive(Clone, Copy)]
struct Stochastic {
    state: State,
    lambda: f32,
    del: f32,
    a: f32,
    b: f32,
    t: f32,
    t_max: f32,
}

impl IsModel for Stochastic {
    type State = State;

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

impl Stochastic {
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
