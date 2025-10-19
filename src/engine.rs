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

use super::closet::{IsEnv, IsObs, IsState, IsSystem, Result};

pub fn simulate<D: IsObs + Clone, S: IsState<D> + Clone, E: IsEnv>(
    sys: &mut (impl IsSystem<D, S, E> + Clone),
) -> Vec<Result<D>> {
    let mut res: Vec<Result<D>> = Vec::new();

    res.push(sys.get());
    while sys.cond() {
        sys.step(sys.decide(sys.suggest()));
        res.push(sys.get());
    }
    res
}
