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

use proc_macro::TokenStream;
use quote::quote;
use syn::{parse, parse_macro_input, Expr, Token};

struct ReadStruct {
    var: Expr,
    name: Expr,
}

impl parse::Parse for ReadStruct {
    fn parse(input: parse::ParseStream) -> syn::Result<Self> {
        let var: Expr = input.parse().expect("macro first parsing failed");
        input
            .parse::<Token![,]>()
            .expect("macro comma parsing failed");
        let name: Expr = input.parse().expect("macro third parsing failed");
        Ok(ReadStruct { var, name })
    }
}

#[proc_macro]
pub fn read_var(input: TokenStream) -> TokenStream {
    let ReadStruct { var, name } = parse_macro_input!(input as ReadStruct);

    quote! {
        let mut tmp = String::new();
        print!("{}: ", stringify!(#name));
        io::stdout().flush().unwrap();
        io::stdin().read_line(&mut tmp).expect("unable to read line");
        #var = tmp.trim().parse().expect("not a valid value")
    }
    .into()
}
