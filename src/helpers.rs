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

pub fn binary_search(x: f32, arr: &[f32]) -> Option<usize> {
    let mut a: usize = 0;
    let mut b: usize = arr.len() - 1;
    if x > arr[b] || x < arr[a] {
        return None;
    }
    while a < b {
        let mut m = (a + b) / 2;
        if arr[m + 1] <= x {
            a = m + 1
        } else if arr[m] > x {
            b = m
        } else {
            while m > 0 && arr[m - 1] == x {
                m -= 1;
            }
            return Some(m);
        }
    }
    None
}

pub fn sigmoid(x: f32) -> f32 {
    1.0 / (1.0 + (-x).exp())
}
