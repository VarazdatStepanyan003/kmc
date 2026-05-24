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

pub fn binary_search(x: f32, arr: &[f32]) -> Option<usize> {
    if x == arr[0] {
        return Some(0);
    }
    if x < arr[0] {
        return None;
    }
    let mut b: usize = arr.len() - 1;
    if x > arr[b] {
        return None;
    }
    let mut a: usize = 0;
    while a < b {
        let mut m = (a + b) / 2;
        if arr[m + 1] <= x {
            a = m + 1;
        } else if arr[m] > x {
            b = m;
        } else {
            while m > 0 && arr[m - 1] == x {
                m -= 1;
            }
            return Some(m);
        }
    }
    panic!("This should never happen: binary_search")
}

pub fn sigmoid(x: f32) -> f32 {
    1.0 / (1.0 + (-x).exp())
}
