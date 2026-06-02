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
use std::sync::{Arc, Mutex};
use std::{env, thread};

pub fn main() {
    let tmp: Vec<String> = env::args().skip(1).collect();

    let (lambda, del, a, b, c, dt, rep_num): (f32, f32, f32, f32, f32, f32, usize);

    if tmp.len() < 7 {
        lambda = 1.0;
        del = 0.5;
        a = 1.0;
        b = 0.5;
        c = 0.5;
        dt = 0.1;
        rep_num = 16;
    } else {
        lambda = tmp[0].parse::<f32>().expect("Wrong Lambda");
        del = tmp[1].parse::<f32>().expect("Wrong Delta");
        assert!(del.abs() <= 1.0);
        a = tmp[2].parse::<f32>().expect("Wrong A");
        b = tmp[3].parse::<f32>().expect("Wrong B");
        c = tmp[4].parse::<f32>().expect("Wrong C");
        dt = tmp[5].parse::<f32>().expect("Wrong Time Step");
        rep_num = tmp[6].parse::<usize>().expect("Wrong Repetition Number");
    }

    let t_max = 100.0 / lambda;

    let e = Env {
        lambda,
        del,
        a,
        b,
        c,
        t_max,
    };

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
    props: Vec<f32>,
    am: Vec<u32>,
}

impl Res {
    fn new(dt: f32, t_max: f32) -> Res {
        let size = (t_max / dt).ceil() as usize;
        let mut time: Vec<f32> = Vec::new();
        let props: Vec<f32> = vec![0.0; size];
        let am: Vec<u32> = vec![0; size];
        for i in 0..=size {
            time.push((i as f32) * dt);
        }
        Res { time, props, am }
    }

    fn add(&mut self, r: Observables) {
        if let Some(i) = helpers::binary_search(r.t, &self.time) {
            self.props[i] += r.prop;
            self.am[i] += 1;
        }
    }

    fn ready(&mut self) {
        self.am.iter().enumerate().rev().for_each(|(i, a)| {
            if *a != 0 {
                self.props[i] /= *a as f32;
            } else {
                self.time.remove(i);
                self.props.remove(i);
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

        self.props.iter().for_each(|o| {
            s.push_str(&o.to_string());
            s.push(',');
        });
        s.pop();
        s.push('\n');

        s
    }
}

struct Observables {
    prop: f32,
    t: f32,
}

struct State {
    prop: f32,
    eps: i8,
}

#[derive(Clone, Copy)]
struct Env {
    lambda: f32,
    del: f32,
    a: f32,
    b: f32,
    c: f32,
    t_max: f32,
}

impl Env {
    fn create(self) -> Stochastic {
        Stochastic {
            state: State { prop: 0.5, eps: 1 },
            lambda: self.lambda,
            del: self.del,
            a: self.a,
            b: self.b,
            c: self.c,
            aux: self.c.abs() + self.b.abs() + self.a.abs(),
            t: 0.0,
            t_max: self.t_max,
            rng: rand::rng(),
        }
    }
}

struct Stochastic {
    state: State,
    lambda: f32,
    del: f32,
    a: f32,
    b: f32,
    c: f32,
    aux: f32,
    t: f32,
    t_max: f32,
    rng: ThreadRng,
}

impl IsModel for Stochastic {
    type Obs = Observables;

    fn get(&self) -> Observables {
        Observables {
            t: self.t,
            prop: self.state.prop,
        }
    }

    fn step(&mut self) {
        let u: f32 = -(1.0 - self.rng.random::<f32>()).ln();

        let dt = u / (2.0 * self.lambda);

        self.ode_solv(dt);

        if 2.0 * rand::random::<f32>() < 1.0 + (self.state.eps as f32) * self.del {
            self.state.eps *= -1;
        }
        self.t += dt;
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
        let mut n = 1;
        let mut h = t / 10.0;
        while h * self.aux > 0.01 {
            n += 1;
            h /= 10.0
        }
        for _ in 0..10i32.pow(n) {
            self.state.prop += self.ode(self.state.prop + self.ode(self.state.prop) * h / 2.0) * h;
        }
    }

    fn ode(&self, x: f32) -> f32 {
        x * (1.0 - x) * ((self.state.eps as f32) * self.c + self.a - self.b * x)
    }
}
