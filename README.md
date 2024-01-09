# Installing
1. [Install Cargo](https://doc.rust-lang.org/cargo/getting-started/installation.html) (Rust package manager)
2. Clone this repository
3. Move to your local copy of the repository and execute the command
 ```
cargo build --release
```
4. Execute the program by typing
```
./target/release/rust <path_to_nwk_file> <output_directory> <path_to_transfers_file> <seed>
```

# Example
1. Type
   ```
   ./target/release/gene_tree_sim "./tree4.nwk" "./test_folder" "./transfer_file.txt"
   ```
2. This will simulate 1000 transfers on 1024 gene trees for the species tree tree4.nwk, in test_folder/tree_0. Note that if tree4.nwk contains several trees separated by commas, each of these trees will be treated in a folder test_folder/tree_i.
