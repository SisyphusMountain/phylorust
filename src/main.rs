
/*
1. Define the data structures used to represent the tree parenthood information
2. Use the crate pest to convert the dataframe into a newich string.



3. We can use SORTING on vectors, and then use a binary search to find the element corresponding to the drawn random number
4. Whenever there are two transfers between two nodes, we can just keep the last one.


Gros problème pour des arbres avec énormément de noeuds ou pour le parallélisme : la complexité en mémoire est O(n^2) où n est le nombre de noeuds,
à cause du stockage de la matrice de contemporanéité
*/

// ATTENTION : Le code tel quel ne fonctionne pas si les noeuds intérieurs n'ont pas de nom,
// et si la racine n'a pas de longueur.
#[macro_use]
extern crate pest_derive;
#[allow(dead_code)]
use pest::Parser;
use newick;
use std::fs;
use std::env;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use prettytable::{Table, format, row, cell};
use time::Instant;
extern crate regex;
use regex::Regex;
use rand::{thread_rng, Rng};
use rand::seq::SliceRandom;
use pest::error::Error;
use csv::Writer;

// ATTENTION : les enfants et parents d'un noeud devront être mutables, ainsi que la longueur.

#[derive(Parser)]
#[grammar = "newick.pest"]
struct NewickParser;
#[derive(Clone, Debug)]
struct Node {
    name: String,
    left_child: Option<Box<Node>>,
    right_child: Option<Box<Node>>,
    parent: Option<usize>, // Using index for parent
    depth: Option<f64>,
    length: f64,
}

#[derive(Clone)]
struct Tree {
    nodes: Vec<Node>,
    root: Option<usize>,
}
#[derive(Clone, Debug)]
struct FlatNode {
    name: String,
    left_child: Option<usize>,
    right_child: Option<usize>,
    parent: Option<usize>,
    depth: Option<f64>,
    length: f64,
}
/* I will later create a subdivision, which will be a vector where all 
elements have a duration and a contemporaneity.*/
struct DurationAndContemporaneity {
    duration: f32,
    contemporaneity: Vec<usize>,
}
fn print_flat_node_table(flat_nodes: &[FlatNode]) {
    let mut table = Table::new();
    table.set_format(*format::consts::FORMAT_NO_LINESEP_WITH_TITLE);

    // Adding the header
    table.add_row(row!["Name", "Left Child", "Right Child", "Parent", "Depth", "Length"]);

    // Adding the rows
    for flat_node in flat_nodes {
        table.add_row(row![
            flat_node.name,
            flat_node.left_child.map_or("None".to_string(), |v| v.to_string()),
            flat_node.right_child.map_or("None".to_string(), |v| v.to_string()),
            flat_node.parent.map_or("None".to_string(), |v| v.to_string()),
            flat_node.depth.map_or("None".to_string(), |v| format!("{:.6}", v)),
            format!("{:.6}", flat_node.length),
        ]);
    }

    // Printing the table
    table.printstd();
}
// The following two functions are used to convert a newick string into an arborescent Tree structure, where each "Node" object owns its two children in a Box object.
fn newick_to_tree(pair: pest::iterators::Pair<Rule>) -> Vec<Node> {
    let mut vec_trees:Vec<Node> = Vec::new();
    for inner in pair.into_inner() {
        let tree = handle_pair(inner);
        vec_trees.push(tree.unwrap());
    }
    return vec_trees
}
fn handle_pair(pair: pest::iterators::Pair<Rule>) -> Option<Node> {
    match pair.as_rule() {
        Rule::newick => {
            // newick = { subtree ~ ";" }
            // For now, we only do the implementation for one single tree.
            // This case is trivial : we just call handle_pair on the only tree we found.
            for inner in pair.into_inner() {
                if inner.as_rule() == Rule::subtree {
                    return handle_pair(inner);
                }
            }
            None
        },
        Rule::subtree => {
            // subtree = { leaf | internal }
            // Subtree is like the choice between an inner node and a leaf. Either way, we need to pass it to handle_pair.
            handle_pair(pair.into_inner().next().unwrap())
        },
        Rule::leaf => {
            // Choose default values for the name and length of the leaf.
            // The defaults are an empty string and a 0.0 length.
            let mut name = String::new();
            let mut length = 0.0;

            // leaf = { NAME? ~ ":"? ~ LENGTH? }
            // Therefore, we use a match statement to handle the cases NAME and LENGTH because ":" is unimportant.
            for inner_pair in pair.into_inner() {
                match inner_pair.as_rule() {
                    Rule::NAME => {
                        name = inner_pair.as_str().to_string();
                    }
                    Rule::LENGTH => {
                        let val = inner_pair.as_str();
                        if let Err(_) = val.parse::<f64>() {
                            println!("Failed to parse LENGTH: {}", val);
                        }
                        length = val.parse::<f64>().unwrap_or(0.0);
                    },
                    
                    _ => {} // Ignore other rules
                }
            }
            let node = Node {
                name: name,
                left_child: None,
                right_child: None,
                parent: None,
                depth: None,
                length: length,
            };
            Some(node)
        },
        Rule::internal => {
            // internal = { "(" ~ subtree ~ "," ~ subtree ~ ")" ~ NAME? ~ ":"? ~ LENGTH? }
            
            // Initialize default values for the name and length of the internal node.
            // The defaults are an empty string and a 0.0 length.
            let mut name = String::new();
            let mut length = 0.0;
        
            let mut first_subtree = None;
            let mut second_subtree = None;
        
            // Iterate through the inner rules without assuming their order
            for inner_pair in pair.into_inner() {
                match inner_pair.as_rule() {
                    Rule::subtree => {
                        let subtree = handle_pair(inner_pair).unwrap();
                        if first_subtree.is_none() {
                            first_subtree = Some(subtree);
                        } else {
                            second_subtree = Some(subtree);
                        }
                    },
                    Rule::NAME => {
                        name = inner_pair.as_str().to_string();
                    },
                    Rule::LENGTH => {
                        let val = inner_pair.as_str();
                        if let Err(_) = val.parse::<f64>() {
                            println!("Failed to parse LENGTH: {}", val);
                        }
                        length = val.parse::<f64>().unwrap_or(0.0);
                    },
                    
                    _ => {} // Ignore other rules
                }
            }
        
            let node = Node {
                name,
                left_child: first_subtree.map(Box::new),
                right_child: second_subtree.map(Box::new),
                parent: None,
                depth: None,
                length,
            };
            Some(node)
        }
        Rule::NAME | Rule::LENGTH => {
            // We should never directly handle these outside their containing rules.
            None
        },
    }
}
// Convert a Tree in "Node" form to a Newick string.
fn node_to_newick(node: &Node) -> String {
    /* Takes a node and returns the corresponding subtree in Newick format.
        --------------------------------
        INPUT:
            - node: the node to convert to Newick format.

        OUTPUT:
            - the Newick representation of the subtree rooted at node.
        Warning: rounds the lengths to 6 decimal places.
    */
    if let (Some(left_child), Some(right_child)) = (&node.left_child, &node.right_child) {
        // This is an internal node with both left and right children.
        format!(
            "({},{}){}:{:.6}",
            node_to_newick(left_child),
            node_to_newick(right_child),
            node.name,
            node.length
        )
    } else {
        // This is a leaf node.
        format!("{}:{:.6}", node.name, node.length)
    }
}
// Convert from FlatNode to Node
fn flat_to_node(flat_tree: &[FlatNode], index: usize, parent_index: Option<usize>) -> Option<Node> {
    /*
    This function converts a flat tree into an arborescent tree.

    Input: a flat tree, the index of a node in the flat tree, and the index of its parent.
    Output: the corresponding node in the arborescent tree.

    Warning: can bug if applied directly to a non-root node, which will be mistaken for a root node.
    */
    let flat_node = &flat_tree[index];
    let left_child = flat_node.left_child.and_then(|i| {
        flat_to_node(flat_tree, i, Some(index)).map(Box::new)
    });
    let right_child = flat_node.right_child.and_then(|i| {
        flat_to_node(flat_tree, i, Some(index)).map(Box::new)});
    
    Some(Node {
        name: flat_node.name.clone(),
        left_child: left_child,
        right_child: right_child,
        parent: parent_index,
        depth: flat_node.depth,
        length: flat_node.length,
    })
}
// Convert from Node to FlatNode
fn node_to_flat(node: &Node, flat_tree: &mut Vec<FlatNode>, parent: Option<usize>) -> usize {
    // Transforms the arborescent tree into a "linear tree", which is just a vector of nodes with parent and descendants.

    let index = flat_tree.len();
    flat_tree.push(FlatNode {
        name: node.name.clone(),
        left_child: None,  // Will fill this in a moment
        right_child: None, // Will fill this in a moment
        parent: parent,
        depth: node.depth,
        length: node.length,
    });

    if let Some(left) = &node.left_child {
        let left_index = node_to_flat(left, flat_tree, Some(index));
        flat_tree[index].left_child = Some(left_index);
    }

    if let Some(right) = &node.right_child {
        let right_index = node_to_flat(right, flat_tree, Some(index));
        flat_tree[index].right_child = Some(right_index);
    }

    index
}
fn traverse_tree(node: &Node) {
    println!("Node: {}", node.name);

    if let Some(left_child) = &node.left_child {
        println!("Traversing left child of {}", node.name);
        traverse_tree(left_child);
    }

    if let Some(right_child) = &node.right_child {
        println!("Traversing right child of {}", node.name);
        traverse_tree(right_child);
    }
}
fn total_length_of_flat_tree(flat_tree: &[FlatNode]) -> f64 {
    flat_tree.iter().map(|node| node.length).sum()
}
fn give_depth(node: &mut Node, depth: f64) {
    node.depth = Some(depth);
    if let Some(left_child) = &mut node.left_child {
        give_depth(left_child, depth + left_child.length);
    }
    if let Some(right_child) = &mut node.right_child {
        give_depth(right_child, depth + right_child.length);
    }
}
fn depths_to_lengths(node: &mut Node, parent_depth: f64) {
    let depth = node.depth.unwrap();
    if node.parent.is_some() {
        node.length = depth - parent_depth;
    }
    let parent_depth = node.depth.unwrap();
    if let Some(left_child) = &mut node.left_child {
        depths_to_lengths(left_child, depth);
    }
    if let Some(right_child) = &mut node.right_child {
        depths_to_lengths(right_child, depth);
    }
}
fn make_subdivision(flat_tree: &mut Vec<FlatNode>) -> Vec<f64> {
    /* Input: a flat tree.
    Output: the time subdivision induced by the nodes of the tree.
    Example : If the original nwk tree was ((A:1,B:2)C:1,D:5)R:0, the subdivision should be [0,1,2,3,5]
    */
    let mut depths: Vec<f64> = flat_tree.iter().filter_map(|node| node.depth).collect();
    depths.sort_by(|a, b| a.partial_cmp(b).unwrap());
    depths.dedup();
    return depths;
}
fn make_intervals(depths: &Vec<f64>) -> Vec<f64> {
    /* Input: a vector of depths.
    Output: a vector of intervals, where the i-th interval is the difference between the i-th and (i+1)-th depth.
    Adds a 0 interval at the beginning so the size of the intervals vector is the same as the size of the depths vector.
    Example: if the depths are [0,1,2,3,5], the intervals are [0,1,1,1,2].
    */
    let mut intervals:Vec<f64> = Vec::new();
    intervals.push(0.0 as f64);
    for i in 0..depths.len()-1 {
        intervals.push(depths[i+1] - depths[i]);
    }
    return intervals;
}
fn find_closest_index(contemporaneous_species_vector: &[f64], value: f64) -> usize {
    match contemporaneous_species_vector.binary_search_by(|probe| probe.partial_cmp(&value).unwrap()) {
        Ok(idx) => idx,
        Err(idx) => {
            if idx == 0 {
                0
            } else if idx == contemporaneous_species_vector.len() {
                contemporaneous_species_vector.len() - 1
            } else {
                // Determine which of contemporaneous_species_vector[idx - 1] or contemporaneous_species_vector[idx] is closer to the value
                if (value - contemporaneous_species_vector[idx - 1]).abs() < (value - contemporaneous_species_vector[idx]).abs() {
                    idx - 1
                } else {
                    idx
                }
            }
        }
    }
}
fn find_contemporaneity(flat_tree: &mut Vec<FlatNode>, depths: &Vec<f64>) -> Vec<Vec<usize>> {
    // Returns a vector of vectors, where each vector contains the indices of the nodes that are contemporaneous.
    // Quadratic time complexity for now, maybe we can improve.
    let mut contemporaneity: Vec<Vec<usize>> = vec![Vec::new(); depths.len()];
    for i in 0..flat_tree.len() {
        // Attention à la racine ! L'indice du début et l'indice de la fin pour la racine devraient être identiques...
        let start_time = flat_tree[i].depth.unwrap() - flat_tree[i].length;
        let end_time = flat_tree[i].depth.unwrap();
        let start_index = find_closest_index(&depths, start_time);
        let end_index = find_closest_index(&depths, end_time);
        //println!("Start index: {}, end index: {}", start_index, end_index);
        for j in start_index+1..end_index+1 {
            // We don't count the start index, because the node is not alive on the interval that ends at the start index.
            contemporaneity[j].push(i);
        }
    }
    return contemporaneity;
}
fn number_of_species(contemporaneity: &Vec<Vec<usize>>) -> Vec<f64> {
    let mut number_of_species: Vec<f64> = Vec::new();
    for i in 0..contemporaneity.len() {
        number_of_species.push(contemporaneity[i].len() as f64);
    }
    return number_of_species;
}
fn make_CDF(intervals: Vec<f64>, number_of_species: Vec<f64>) -> Vec<f64> {
    // Returns the CDF of the depth of a random time of a transfer, uniformly distributed on branches.
    let n = intervals.len();
    let mut CDF: Vec<f64> = Vec::new();
    CDF.push(intervals[0]*number_of_species[0]);
    for i in 1..n {
        CDF.push(CDF[i-1] + intervals[i]*number_of_species[i]);
    }
    let total_value = CDF[n-1];
    for i in 0..n {
        CDF[i] = CDF[i]/total_value;
    }
    return CDF;
}
fn choose_from_CDF(CDF: &Vec<f64>, depths: &Vec<f64>) -> (f64, usize) {
    /* The CDF allows us to compute the probability of a given interval being chosen,
    by picking a random number between 0 and 1 and finding the interval it falls into.
    First, we compute the index of the first depth that is larger than the random number.
    
    This function computes the depth that was chosen. To compute this depth, we need to
    know where exactly the random number fell in the interval, and then compute which
    depth this corresponds to, via
    depth = (r - CDF[i-1])/(CDF[i] - CDF[i-1]) * (depths[i] - depths[i-1]) + depths[i-1]

     returns the index of the smallest depth larger than the depth that was picked
    */

    let mut rng = thread_rng();
    let r: f64 = rng.gen_range(0.0..1.0);
    let index = match CDF.binary_search_by(|&probe| probe.partial_cmp(&r).unwrap_or(std::cmp::Ordering::Less)) {
        Ok(index) => index, // Found an exact match, use `index`.
        Err(index) => {
            if index == CDF.len() {
                // `r` is greater than any number in CDF.
                // Return the last index, or handle this case differently depending on your needs.
                panic!("r is greater than any number in CDF, including the last one which should be 1.0. The CDF is: {:?}", CDF);
            } else {
                // No exact match, but `index` is where `r` would be inserted to maintain order.
                // This means `CDF[index]` is the first value greater than or equal to `r`.
                index
            }
        }
    };
    //println!("r: {}, index: {}", r, index);
    let depth = (r - CDF[index-1])/(CDF[index] - CDF[index-1]) * (depths[index] - depths[index-1]) + depths[index-1];
    return (depth, index);
}
fn random_pair(vec: &Vec<usize>) -> Option<(usize, usize)> {
    if vec.len() < 2 {
        return None;
    }

    let mut rng = thread_rng();
    let first_index = rng.gen_range(0..vec.len());
    let mut second_index = rng.gen_range(0..vec.len());
    while second_index == first_index {
        second_index = rng.gen_range(0..vec.len());
    }

    Some((vec[first_index], vec[second_index]))
}
fn generate_transfers(n_transfers: usize, contemporaneity: &Vec<Vec<usize>>, CDF: &Vec<f64>, depths: &Vec<f64>) -> Vec<(usize, usize, f64)> {
    /* Generates an arbitrary number of uniformly distributed transfers.
    */
    let mut transfers: Vec<(usize, usize, f64)> = Vec::new();
    for _ in 0..n_transfers {
        let (depth, index) = choose_from_CDF(CDF, depths);
        let contemporaneous_species_vector = &contemporaneity[index];
        if let Some((first_species, second_species)) = random_pair(contemporaneous_species_vector) {
            transfers.push((first_species, second_species, depth));
        }
    }
    transfers.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap());
    return transfers;
}
fn make_transfers_csv(flat_tree: &Vec<FlatNode>, transfers: &Vec<(usize, usize, f64)>, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
    let mut writer = Writer::from_path(filename)?;

    // Writing the header
    writer.write_record(&["Donor", "Recipient", "Depth"])?;

    for transfer in transfers.iter() {
        let donor_name = &flat_tree[transfer.0].name;
        let recipient_name = &flat_tree[transfer.1].name;
        let depth = transfer.2;

        // Writing the data rows
        writer.write_record(&[donor_name, recipient_name, &depth.to_string()])?;
    }

    writer.flush()?;
    Ok(())
}
fn change_tree(flat_tree: &mut Vec<FlatNode>, transfer: (usize, usize, f64)){
    //println!("{:?}, {:?}, {:?}", flat_tree[transfer.0].name, flat_tree[transfer.1].name, transfer.2);
    /* 
    This function computes the changes that occur in the tree when a transfer occurs.
    --------------------------------
    Input: a flat tree, and a transfer (A, B, depth).
    Output: None. The flat tree is modified in place.
    --------------------------------
    By making a drawing, we can observe that there are two cases. Either A and B have the same parent, or not.
    If they have the same parent, the topology does not change, only the lengths of A and B.
    However, if they do not have the same parent, the change in topology is more complicated. We must do the following:
    --------------------------------
    TOPOLOGY:
    A.parent = B.parent
    FA.child which is A = FB
    SB.parent = FFB (which can be None if B's father is the root)
    if FFB != None:
        FFB.child which is FB = SB
    --------------------------------
    METRIC:
    Remark that several lengths can change, but only one depth can change: the depth of FB.
    We apply this change, knowing that we will recalculate all lengths at the end.
    FB.depth = depth






    Explanation:

    We have to understand what should happen when a transfer occurs.
    Given a node N, denote by SN its sister and FN its parent.
    When a transfers occurs from A to B:

    The nodes which are changed in this operations are therefore: FA (new father, new son replacing SA), FB (New son replacing B), B (New father), SA (new father), FFA (new son replacing FA).
    Remark: the depth of the parent FA of A is also updated, because its new depth is the depth of the transfer.
    Remark: There is a special case if the transfer comes from a node that has the root as a parent. In this case, do not change 
    */
    let (first_species, second_species, depth) = transfer;
    // --------------------------------
    //println!("State of the tree before doing anything");
    //print_flat_node_table(&flat_tree);
    // --------------------------------
    // Node indexes
    let A_index = first_species;
    let B_index = second_species;
    let FA_index = flat_tree[A_index].parent.unwrap();
    let FB_index = flat_tree[B_index].parent.unwrap();
    let SB_index;
    if flat_tree[FB_index].left_child.unwrap() == B_index {
        SB_index = flat_tree[FB_index].right_child.unwrap()} else 
        {SB_index = flat_tree[FB_index].left_child.unwrap()};

    let FFB_index_opt = flat_tree[FB_index].parent;


    // --------------------------------

    //println!("A: {}, B: {}, FA: {}, FB: {}", A_index, B_index, FA_index, FB_index);

    //println!("FA: {}, FB: {}\n", FA_index, FB_index);

    //println!("A: {}", flat_tree[A_index].name);
    //println!("B: {}", flat_tree[B_index].name);
    //println!("FA: {}", flat_tree[FA_index].name);
    //println!("FB: {}", flat_tree[FB_index].name);
    //println!("FFA: {}", flat_tree[FFA_index].name);


    if FA_index == FB_index {
            // The two nodes have the same parent. Do nothing to the topology.

            // Update the depth of FB
            flat_tree[FB_index].depth = Some(depth);
        } else
        // The two nodes have the same parent. Do nothing to the topology.
        {
        // Updating A's parent
        flat_tree[A_index].parent = Some(FB_index);
        // Updating FB's parent
        flat_tree[FB_index].parent = Some(FA_index);

        // Changing the parent of SB to FFB.

        flat_tree[SB_index].parent = FFB_index_opt;

        // Changing the child of FA which is A to FB
        if flat_tree[FA_index].left_child == Some(A_index) {
            flat_tree[FA_index].left_child = Some(FB_index)} else{
                flat_tree[FA_index].right_child = Some(FB_index)};
        
        // Changing the child of FB which is SB to A
        if flat_tree[FB_index].left_child == Some(B_index) {
            flat_tree[FB_index].right_child = Some(A_index)} else{
                flat_tree[FB_index].left_child = Some(A_index)};

        // Changing the child of FFB which is FB to SB
        if let Some(FFB_index) = FFB_index_opt {
            if flat_tree[FFB_index].left_child == Some(FB_index) {
                flat_tree[FFB_index].left_child = Some(SB_index)} else
                {flat_tree[FFB_index].right_child = Some(SB_index)
                };
            }
        }
        // Update the depth of FB
        flat_tree[FB_index].depth = Some(depth);

    }

fn create_new_tree(new_flat_tree: &mut Vec<FlatNode>, transfers: Vec<(usize, usize, f64)>) -> () {
    /* Only applies all transfers to the tree
    */
    for transfer in transfers {
        change_tree(new_flat_tree, transfer);
    }
}
fn find_root_in_flat_tree(flat_tree: &Vec<FlatNode>) -> Option<usize> {
    let mut root_index: Option<usize> = None;
    for i in 0..flat_tree.len() {
        if flat_tree[i].parent.is_none() {
            if root_index.is_some() {
                // There is already a node with no parent
                //print_flat_node_table(&flat_tree);
                panic!("There should be only one node with no parent");
            }
            root_index = Some(i);
        }
    }
    return root_index;
}

fn one_gene_sim_to_string(copied_flat_tree: &mut Vec<FlatNode>, n_transfers: usize, contemporaneity: &Vec<Vec<usize>>, cdf: &Vec<f64>, depths: &Vec<f64>, gene_index: u32) -> String {
    // Create transfers
    let mut transfers = generate_transfers(n_transfers, &contemporaneity, &cdf, &depths);
    // Export them to a CSV file
    let filename = format!("transfers_{}.csv", gene_index);
    make_transfers_csv(copied_flat_tree, &transfers, filename.as_str()).expect("Failed to write CSV file");
    // Generate gene tree as a flat tree
    create_new_tree(copied_flat_tree, transfers);

    let root_of_gene_tree = find_root_in_flat_tree(copied_flat_tree).expect("There should be a root in the gene tree");
    // flat tree to classical arborescent tree
    let mut reconstructed_tree = flat_to_node(copied_flat_tree, root_of_gene_tree, None).unwrap();
    // The lengths of the nodes are not correct, so we need to update them.
    let root_depth = reconstructed_tree.depth.unwrap();
    depths_to_lengths(&mut reconstructed_tree, root_depth);
    // Convert the arborescent tree to a newick string
    let reconstructed_newick = node_to_newick(&reconstructed_tree) + ";";
    let mut gene_file = File::create(format!("gene_{}.nwk", gene_index)).expect("Failed to create file");
    gene_file.write_all(reconstructed_newick.as_bytes()).expect("Failed to write to file");
    return reconstructed_newick;
}

fn create_many_genes(copied_flat_tree: &mut Vec<FlatNode>, n_transfers_vec: Vec<usize>, contemporaneity: &Vec<Vec<usize>>, cdf: &Vec<f64>, depths: &Vec<f64>) -> () {
    for (i, n_transfers) in n_transfers_vec.iter().enumerate(){
        let mut copied_flat_tree2 = copied_flat_tree.clone();
        let n_transfers_usize = *n_transfers; // Convert &usize to usize, then to u32
        one_gene_sim_to_string(&mut copied_flat_tree2, n_transfers_usize, contemporaneity, cdf, depths, i as u32);
        //print_flat_node_table(&copied_flat_tree2);
    }
}


fn tree_length(flat_tree: &Vec<FlatNode>) -> f64 {
    let mut length: f64 = 0.0;
    for i in 0..flat_tree.len() {
        length += flat_tree[i].length;
    }
    return length;
}
fn main() {
    // Read command line arguments
    env::set_var("RUST_BACKTRACE", "1");
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        eprintln!("Usage: {} <path_to_nwk_file>", args[0]);
        return;
    }

    // Read the .nwk file content
    let start = Instant::now();
    let content = fs::read_to_string(&args[1]).expect("Failed to read the file");
    let re = Regex::new(r"\s+").unwrap();  // Matches any whitespace characters
    let sanitized_content = re.replace_all(&content, "");
    // Split the content by ; to get individual trees
    let trees: Vec<String> = sanitized_content.split(';')
    .filter_map(|tree| {
        let trimmed = tree.trim();
        if trimmed.is_empty() {
            None
        } else {
            Some(trimmed.to_string() + ";")
        }
    })
    .collect();


    let mut number_of_trees = 0;

    for tree in trees {
        // Parse each tree to get the Node representation
        let pairs = NewickParser::parse(Rule::newick, &tree).unwrap_or_else(|e| panic!("fail"));
        for pair in pairs {
            number_of_trees += 1;
            let mut vec_trees = newick_to_tree(pair);
            let mut example_tree = vec_trees.pop().unwrap();
            give_depth(&mut example_tree, 0.0);

            // Convert Node representation to FlatNode representation
            let mut flat_tree = Vec::new();
            node_to_flat(&example_tree, &mut flat_tree, None);
            flat_to_node(&flat_tree, 0, None);
            //println!("Flat_tree: {:?}", flat_tree);

            let depths = make_subdivision(&mut flat_tree);
            //println!("Depths: {:?}", depths);
            let intervals = make_intervals(&depths);
            let contemporaneity = find_contemporaneity(&mut flat_tree, &depths);
            //println!("Contemporaneity: {:?}", contemporaneity);
            let n_species = number_of_species(&contemporaneity);
            //println!("Number of species: {:?}", n_species);
            //println!("intervals: {:?}", intervals);
            let cdf = make_CDF(intervals, n_species);
            //println!("CDF in main: {:?}", cdf);
            let choice = choose_from_CDF(&cdf, &depths);



            let mut copied_flat_tree = flat_tree.clone();
            create_many_genes(&mut copied_flat_tree, vec![10000;1000], &contemporaneity, &cdf, &depths);
            //print_flat_node_table(&copied_flat_tree);
            //println!("copied flat tree before: {:?}", copied_flat_tree);
            /*
            let transfers = generate_transfers(1000, &contemporaneity, &cdf, &depths);
            create_new_tree(&mut copied_flat_tree, transfers);*/
            //one_gene_sim_to_string(&mut copied_flat_tree, 1000, &contemporaneity, &cdf, &depths);
            //one_gene_simulation(&mut copied_flat_tree, 1000, &contemporaneity, &cdf, &depths);

            //println!("Root of tree: {}", root_of_tree);


        }
    }

    let duration = start.elapsed();
    println!("Number of trees: {}", number_of_trees);
    println!("Time taken: {} seconds", duration.as_seconds_f64());
}
fn compare_trees(node1: &Node, node2: &Node) -> bool {
    if node1.name != node2.name || (node1.length - node2.length).abs() > 1e-6 {
        return false;
    }

    match (&node1.left_child, &node2.left_child) {
        (Some(l1), Some(l2)) if !compare_trees(l1, l2) => return false,
        (None, None) => {},
        _ => return false,
    }

    match (&node1.right_child, &node2.right_child) {
        (Some(r1), Some(r2)) if !compare_trees(r1, r2) => return false,
        (None, None) => {},
        _ => return false,
    }
    
    true
}