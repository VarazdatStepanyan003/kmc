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

#![allow(dead_code)]
pub mod closet;
pub mod engine;
pub mod helpers;
pub use kmc_derive::Observable;
mod tests;

#[cfg(test)]
mod lib {
    use super::tests::pi;

    #[test]
    fn calculate_pi() {
        assert!(pi::test());
    }
}
