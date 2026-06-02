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

use crate::closet::IsModel;

pub fn simulate<Model: IsModel>(sys: &mut Model) -> Vec<<Model as IsModel>::Obs> {
    let mut res: Vec<<Model as IsModel>::Obs> = Vec::new();

    res.push(sys.get());
    while sys.cond() {
        sys.step();
        if sys.store_cond() {
            res.push(sys.get());
        }
    }
    res
}
