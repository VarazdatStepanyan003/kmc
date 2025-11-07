use proc_macro::TokenStream;
use quote::quote;

#[proc_macro_derive(IsObs)]
pub fn obs_derive(input: TokenStream) -> TokenStream {
    let ast: syn::DeriveInput = syn::parse(input).unwrap();

    let name = ast.ident;
    quote! {
        impl IsObs for #name {

        }
    }
    .into()
}
