/*
SCRIPT gene_tree_sim
--------------------------------------------------------------------------------
This script takes:
- a tree in a format like such: ((a:1.0,b:1.0)c:1.0,d:2.0)e:0.0; (with necessary ; at the end, and no [&R] at the beginning)
- a path to the output directory
- a path to the "transfer file", which is a text file containing a list
of number of transfers to simulate for each gene tree that the algorithm will produce. example:  23,10,11100,121,
- a seed for the random number generator.
--------------------------------------------------------------------------------
Output:
new files created: 
- output_folder/genes/gene_i.nwk (newick gene tree)
with format (a:1.000000,(b:0.056798,d:0.056798)e:0.943202)c:1.000000;
- output_folder/transfers/transfers_i.csv
with format 
Donor,Recipient,Depth
b,d,1.9432018717076178
--------------------------------------------------------------------------------


TODO: The function computing the time subdivision probably creates too many values, because the terminal contemporaneous nodes
may all have similar but different depths due to approximation errors.


Remark: precomputing the contemporaneity vector may give the fastest results if we have many gene trees,
but it is computationally heavy and takes up a lot of memory (up to about 1GB for 100000 leaves trees).

Let $n$ be the number of leaves in the species tree, $g$ the number of gene trees to be computed, and $t$ the number of transfers.
I think the time complexity is about $O(n^2) + O(t*g)$



Rename contemporanous_vector in find_closest? Appropriate name?
*/

#[macro_use]
extern crate pest_derive;
use pest::Parser;
use std::fs;
use std::io::{self, Write};
use std::path::Path;
use std::env;
use std::fs::File;
use prettytable::{Table, format, row};
use time::Instant;
extern crate regex;
use regex::Regex;
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;
use csv::Writer;



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


#[derive(Clone, Debug)]
struct FlatNode {
    name: String,
    left_child: Option<usize>,
    right_child: Option<usize>,
    parent: Option<usize>,
    depth: Option<f64>,
    length: f64,
}
// --------------------------------
// FUNCTIONS ONLY USED FOR DEBUGGING
fn print_flat_node_table(flat_nodes: &[FlatNode]) {
    /*
    Prints out the current state of all nodes in a flat tree, in a nice table.
    --------------------------------
    Input:
    - An immutable reference to a vector of FlatNode objects
    --------------------------------
    Output:
    None. Print out a table in the terminal for debugging purposes.
    */

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
fn compare_trees(node1: &Node, node2: &Node) -> bool {
    /* Checks that two Nodes representing the subtree of which they are the root, are equal (in name, length and topology)
    --------------------------------
    Input:
    - Two Node objects
    --------------------------------
    Output:
    A boolean value indicating whether the subtrees are equal (equal names)
    */
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
fn traverse_tree(node: &Node) {
    /* Prints out the parenthood relationships of the nodes. Less useful than print_flat_node_table
    */
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
    /*
    Returns the length of the tree (sum of lengths of all branches) in a tree.
    --------------------------------    
    Input:
    - A reference to a vector of FlatNode objects.
    --------------------------------
    Output:
    - An f64 float
    --------------------------------
    */
    flat_tree.iter().map(|node| node.length).sum()
}
// --------------------------------
// The following two functions are used to convert a newick string into an arborescent Tree structure, where each "Node" object owns its two children in a Box object.
fn newick_to_tree(pair: pest::iterators::Pair<Rule>) -> Vec<Node> {
    /*
    The input .nwk file of the script is a concatenation of newick trees, separated by ";".
    This function returns a vector containing each of these trees.
    --------------------------------
    Input:
    - A "pair" object of type "newick" (see newick.pest grammar) which is a newick string representing one single tree.
    ----------------------------------
    Output:

    */
    let mut vec_trees:Vec<Node> = Vec::new();
    for inner in pair.into_inner() {
        let tree = handle_pair(inner);
        vec_trees.push(tree.unwrap());
    }
    return vec_trees
}
fn handle_pair(pair: pest::iterators::Pair<Rule>) -> Option<Node> {
    /* 
    Recursively parses a newick string representing a single newick tree, according to the grammar contained in newick.pest.
    ----------------------------------
    Input:
    - A "pair" object representing a part of this string.
    ----------------------------------
    Output:
    - Either a node if the input string represents a node, or nothing in cases where it does not (e.g. if the string is the length of a node).
    */
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
// --------------------------------
fn tree_length(flat_tree: &Vec<FlatNode>) -> f64 {
    let mut length: f64 = 0.0;
    for i in 0..flat_tree.len() {
        length += flat_tree[i].length;
    }
    return length;
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
    This function converts a flat tree into an arborescent tree recursively.
    To use it, give the flat_tree, as well as the index of the root. The vector will be traversed recursively, 
    following the descendants of each node being examined.
    ----------------------------------
    Input: 
    - A flat tree, the index of a node in the flat tree, and the index of its parent.
    ----------------------------------
    Output:
    - The corresponding node in the arborescent tree.
    ----------------------------------
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
    /* Transforms the arborescent tree into a "linear tree", which is just a vector of nodes with parent and descendants.
    ----------------------------------
    Input:
    - The root node, which contains the whole arborescent tree
    - The flat_tree vector, which is to be filled by the function
    ----------------------------------
    Output:
    - The index of the node currently being added to flat_tree (usize because the indexing of Rust vectors can only be usize).
    */
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
fn give_depth(node: &mut Node, depth: f64) {
    /*
    Recursively gives the depth of nodes in a "Node" arborescent tree.
    ----------------------------------
    Input:
    - A node, and its depth in the tree
    ----------------------------------
    Output:
    - None. The tree is modified in place via a mutable reference (&mut)
    */
    node.depth = Some(depth);
    if let Some(left_child) = &mut node.left_child {
        give_depth(left_child, depth + left_child.length);
    }
    if let Some(right_child) = &mut node.right_child {
        give_depth(right_child, depth + right_child.length);
    }
}
fn depths_to_lengths(node: &mut Node, parent_depth: f64) {
    /*
    Computes the lengths of nodes from the depths.
    Necessary because the gene transfer functions only modify the depth of the
    transfered node and not lengths for convenience.
    ----------------------------------
    Input:
    - A "Node" object representing an arborescent tree.
    - The depth of its parent
    ----------------------------------
    Output:
    - None. The tree is modified in place via a mutable reference.
    */
    let depth = node.depth.unwrap();
    if node.parent.is_some() {
        node.length = depth - parent_depth;
    }
    if let Some(left_child) = &mut node.left_child {
        depths_to_lengths(left_child, depth);
    }
    if let Some(right_child) = &mut node.right_child {
        depths_to_lengths(right_child, depth);
    }
}
fn make_subdivision(flat_tree: &mut Vec<FlatNode>) -> Vec<f64> {
    /* 
    Makes the time subdivision corresponding to all the times of nodes in the tree.
    ----------------------------------    
    Input: 
    - A flat tree.
    ----------------------------------
    Output: 
    - The time subdivision induced by the nodes of the tree.
    Example : If the original nwk tree was ((A:1,B:2)C:1,D:5)R:0, the subdivision should be [0,1,2,3,5]
    */
    let mut depths: Vec<f64> = flat_tree.iter().filter_map(|node| node.depth).collect();
    depths.sort_by(|a, b| a.partial_cmp(b).unwrap());
    depths.dedup();
    return depths;
}
fn make_intervals(depths: &Vec<f64>) -> Vec<f64> {

    /*
    Makes the vector of time intervals, necessary for computing the contemporaneities between all species in the tree.
    ----------------------------------
    Input:
    - A vector of depths.
    ----------------------------------
    Output:
    - A vector of intervals, where the i-th interval is the difference between the i-th and (i+1)-th depth.
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
    /*
    Since the depths are computed with approximation errors, the depth of a node does not necessarily match any value
    in the time subdivision vector. Therefore, we need to find the closest one.
    This should not be a problem because bisection search is O(log(n)).
    ----------------------------------
    Input:
    - The vector of times for the subdivision.
    - Any given time.
    ----------------------------------
    Output:
    - The index of the closest value.
    */
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
    /*
    Constructs a vector of vectors of indexes, containing the indexes of all species present over each interval.

    ----------------------------------
    Input:
    - The flat_tree, containing all the species.
    - The vector of depths (subdivision)
    ----------------------------------
    Output: 
    - A vector representing contemporaneity.

    */
    let mut contemporaneity: Vec<Vec<usize>> = vec![Vec::new(); depths.len()];
    for i in 0..flat_tree.len() {
        let start_time = flat_tree[i].depth.unwrap() - flat_tree[i].length;
        let end_time = flat_tree[i].depth.unwrap();
        // Find the timestep in the subdivision closest to the computed times.
        let start_index = find_closest_index(&depths, start_time);
        let end_index = find_closest_index(&depths, end_time);

        for j in start_index+1..end_index+1 {// Works even if the start index and the end index are identical (for the root), in this case the loop is empty (range 0..0)
            // We don't count the start index, because the node is not alive on the interval that ends at the start index.
            contemporaneity[j].push(i);
        }
    }
    return contemporaneity;
}
fn number_of_species(contemporaneity: &Vec<Vec<usize>>) -> Vec<f64> {
    /*
    Computes the number of species over each time interval from the subdivision.
    ----------------------------------
    Input:
    - The contemporaneity vector, containing the contemporaneous species over any interval.
    ----------------------------------
    Output:
    - A vector containing the number of species over each of these intervals.
    The type is Vec<f64> because we will later need to multiply the length of each interval to compute the
    probability of choosing any given interval.
    */
    let mut number_of_species: Vec<f64> = Vec::new();
    for i in 0..contemporaneity.len() {
        number_of_species.push(contemporaneity[i].len() as f64);
    }
    return number_of_species;
}
fn make_CDF(intervals: Vec<f64>, number_of_species: Vec<f64>) -> Vec<f64> {
    /* Returns the CDF of the depth of a random time of a transfer, uniformly distributed on branches.
    The CDF is an increasing continuous piecewise affine function, which is continuous in our case.
    Therefore, we only need to give its value at the subdivision points.
    We can compute the probability density of choosing any given timepoint for a transfer:
    this probability density is given by a fixed constant * (length of an interval)*(number of species).
    We can therefore integrate the density to obtain the CDF.
    ----------------------------------
    Input:
    - The intervals vector, containing the lengths of intervals in the subdivision.
    - The number_of_species vector, giving the number of species in any interval of the subdivision.
    ----------------------------------
    Output:
    - The values of the probability density function at all timepoints in the subdivision.
    Let b>a>0 be any two times. Let $Y$ be the random variable which gives the time at which a uniformly distributed transfer occurs.
    Let F be the CDF of $Y$, and let $X$ be a uniformly distributed random variable on [0,1].
    Then 
    */
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
fn choose_from_CDF(CDF: &Vec<f64>, depths: &Vec<f64>, rng: &mut StdRng) -> (f64, usize) {
    /* The CDF allows us to compute the probability of a given interval being chosen,
    by picking a random number between 0 and 1 and finding the interval it falls into.
    First, we compute the index of the first depth that is larger than the random number.
    
    Let b>a>0 be any two times. Let $Y$ be the random variable which gives the time at which a uniformly distributed transfer occurs.
    Let F be the CDF of $Y$, and let $X$ be a uniformly distributed random variable on [0,1].
    Then the probability that $Y$ falls between $a$ and $b$ is $F(a)-F(b) = \int_{F(a)}^{F(b)}{dt} = P(X^{-1}([F(a),F(b)]))$
    So the probability that $Y$ falls between $a$ and $b$ is the probability that $X$ falls between $F(a)$ and $F(b)$.



    This function computes the depth that was chosen. To compute this depth, we need to
    know where exactly the random number fell in the interval, and then compute which
    depth this corresponds to, via
    depth = (r - CDF[i-1])/(CDF[i] - CDF[i-1]) * (depths[i] - depths[i-1]) + depths[i-1]
    ----------------------------------
    Input:
    - The CDF
    - The vector of depths
    ----------------------------------
    Output:
    - A choice of a transfer time as well as an interval index.
    */
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
fn random_pair(vec: &Vec<usize>, rng: &mut StdRng) -> Option<(usize, usize)> {
    /* 
    Picks one species, then tries to pick another one that is contemporaneous until the donor and receiver are different.
    ----------------------------------
    Input:
    - A vector of contemporaneous species over an interval
    ----------------------------------
    Output:
    - A pair of indices representing two species.
    */
    if vec.len() < 2 {
        return None;
    }

    let first_index = rng.gen_range(0..vec.len());
    let mut second_index = rng.gen_range(0..vec.len());
    while second_index == first_index {
        second_index = rng.gen_range(0..vec.len());
    }

    Some((vec[first_index], vec[second_index]))
}
fn generate_transfers(n_transfers: usize, contemporaneity: &Vec<Vec<usize>>, CDF: &Vec<f64>, depths: &Vec<f64>, rng: &mut StdRng) -> Vec<(usize, usize, f64)> {
    /* Generates an arbitrary number of uniformly distributed transfers.
    ----------------------------------
    Input:
    - The desired number of transfers
    - The contemporaneity matrix
    - The CDF
    - The depths representing the time subdivision
    ----------------------------------
    Output:
    - A vector containing all the transfers in the format (donor_index, receiver_index, time)
    */
    let mut transfers: Vec<(usize, usize, f64)> = Vec::new();
    for _ in 0..n_transfers {
        let (depth, index) = choose_from_CDF(CDF, depths, rng);
        let contemporaneous_species_vector = &contemporaneity[index];
        if let Some((first_species, second_species)) = random_pair(contemporaneous_species_vector, rng) {
            transfers.push((first_species, second_species, depth));
        }
    }
    transfers.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap());
    return transfers;
}
fn make_transfers_csv(flat_tree: &Vec<FlatNode>, transfers: &Vec<(usize, usize, f64)>, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
    /*
    Writes a .csv file in the path provided, containing columns (donor, receiver, time).
    ----------------------------------
    Input:
    - The flat_tree
    - The transfers vector
    - The filename giving the path to the csv to be created
    ----------------------------------
    Output:
    - None. The file is written directly in the path provided.
    */
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
    Given a node N, denote by SN its sister and FN its parent.
    --------------------------------
    TOPOLOGY:
    A.parent becomes B.parent
    FA's child which is A becomes FB
    SB's parent becomes FFB (which can be None if B's father is the root)
    if FFB != None:
        FFB's child which is FB becomes SB
    --------------------------------
    METRIC:
    Remark that several lengths can change, but only one depth can change: the depth of FB.
    We apply this change, knowing that we will recalculate all lengths at the end.
    FB.depth = depth
    ----------------------------------
    Explanation: easier with a drawing. The operation we are doing is unplugging the receiver's parent and plugging it onto the donor.
    */
    let (first_species, second_species, depth) = transfer;
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
    /* Only applies all transfers to the tree.
    ----------------------------------
    Input:
    - A copy of the species tree, passed by mutable reference
    - The vector of all transfers occuring in the tree
    ----------------------------------
    Output:
    - None. The copy of the tree is modified in place.
    */
    for transfer in transfers {
        change_tree(new_flat_tree, transfer);
    }
}
fn find_root_in_flat_tree(flat_tree: &Vec<FlatNode>) -> Option<usize> {
    /*
    After all transfers, the root of the tree may have changed. This function finds the root of the tree, which should
    be the unique node without a parent.
    ----------------------------------
    Input:
    - The flat_tree.
    ----------------------------------
    Output:
    - The index of the root in the flat_tree.
    */

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
fn one_gene_sim_to_string(
    copied_flat_tree: &mut Vec<FlatNode>,
    n_transfers: usize,
    contemporaneity: &Vec<Vec<usize>>,
    cdf: &Vec<f64>,
    depths: &Vec<f64>,
    gene_index: u32,
    output_dir: &str,
    rng: &mut StdRng,
) -> Result<String, io::Error> {
    /*
    Performs all the steps necessary to create a gene tree from a species tree.
    ----------------------------------
    Input:
    - A copy of the species tree
    - A number of transfers to be placed on the tree
    - The contemporaneity vector
    - The CDF of the tree
    - The depths representing the time subdivision
    - The index of the gene whose gene tree we are simulating
    - The output directory.
    ----------------------------------
    Output:
    - None. The gene tree .nwk as well as the .csv containing the transfers are both written directly in files.
    */

    // Ensure the output directory exists
    fs::create_dir_all(output_dir)?;

    // Paths for transfers and genes directories
    let transfers_dir = Path::new(output_dir).join("transfers");
    let genes_dir = Path::new(output_dir).join("genes");

    // Ensure the transfers and genes directories exist
    fs::create_dir_all(&transfers_dir)?;
    fs::create_dir_all(&genes_dir)?;

    // Create transfers
    let transfers = generate_transfers(n_transfers, contemporaneity, cdf, depths, rng);

    // Export them to a CSV file in the transfers directory
    let transfer_filename = transfers_dir.join(format!("transfers_{}.csv", gene_index));
    make_transfers_csv(copied_flat_tree, &transfers, transfer_filename.to_str().unwrap())
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;

    // Generate gene tree as a flat tree
    create_new_tree(copied_flat_tree, transfers);

    let root_of_gene_tree = find_root_in_flat_tree(copied_flat_tree)
        .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "There should be a root in the gene tree"))?;
    
    // Flat tree to classical arborescent tree
    let mut reconstructed_tree = flat_to_node(copied_flat_tree, root_of_gene_tree, None)
        .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "Failed to convert flat tree to node tree"))?;
    
    // The lengths of the nodes are not correct, so we need to update them.
    let root_depth = reconstructed_tree.depth
        .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "Root depth not found"))?;
    depths_to_lengths(&mut reconstructed_tree, root_depth);
    
    // Convert the arborescent tree to a Newick string
    let reconstructed_newick = node_to_newick(&reconstructed_tree) + ";";
    
    // Save the gene tree as a Newick string in the genes directory
    let gene_filename = genes_dir.join(format!("gene_{}.nwk", gene_index));
    let mut gene_file = File::create(gene_filename)?;
    gene_file.write_all(reconstructed_newick.as_bytes())?;
    
    Ok(reconstructed_newick)
}
fn create_many_genes(copied_flat_tree: &mut Vec<FlatNode>, n_transfers_vec: Vec<usize>, contemporaneity: &Vec<Vec<usize>>, cdf: &Vec<f64>, depths: &Vec<f64>,output_dir: &str, rng: &mut StdRng) -> () {
    /*
    Repeatedly simulates gene trees.
    ----------------------------------
    Input:
    - A copy of the species tree
    - A vector of the number of transfers to be placed on the tree
    - The contemporaneity vector
    - The CDF of the tree
    - The depths representing the time subdivision
    - The output directory.
    ----------------------------------
    Output:
    - None. The .csv and .nwk files are created on the disk.
    */
    for (i, n_transfers) in n_transfers_vec.iter().enumerate() {
        let mut copied_flat_tree2 = copied_flat_tree.clone();
        let n_transfers_usize = *n_transfers; // Convert &usize to usize
        let _ = one_gene_sim_to_string(&mut copied_flat_tree2, n_transfers_usize, contemporaneity, cdf, depths, i as u32, output_dir, rng);// No useful output variable.
    }
}



fn main() {
    // Read command line arguments
    env::set_var("RUST_BACKTRACE", "1");
    let args: Vec<String> = env::args().collect();
    if args.len() != 5 {
        eprintln!("Usage: {} <path_to_nwk_file> <output_directory> <path_to_transfers_file> <rng_seed>", args[0]);
        return;
    }

    let nwk_file_path = &args[1];
    let output_dir = &args[2];
    let transfers_file_path = &args[3];
    let rng_seed_str = &args[4];

    // Convert rng_seed_str to u64
    let seed = rng_seed_str.parse::<u64>().unwrap_or_else(|e| {
        eprintln!("Error parsing RNG seed: {}", e);
        std::process::exit(1);
    });
    
    let mut rng = StdRng::seed_from_u64(seed);
    // Ensure the output directory exists
    fs::create_dir_all(output_dir).expect("Failed to create output directory");

    let start = Instant::now();
    let content = fs::read_to_string(nwk_file_path).expect("Failed to read the file");
    let re = Regex::new(r"\s+").unwrap();  // Matches any whitespace characters
    let sanitized_content = re.replace_all(&content, "");
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

    let transfers_content = fs::read_to_string(transfers_file_path)
        .expect("Failed to read transfers file");
    let n_transfers: Vec<usize> = transfers_content.trim()
        .split(',')
        .filter_map(|n| n.parse().ok())
        .collect();

    let mut number_of_trees = 0;
    for (index, tree) in trees.iter().enumerate() {
        let tree_output_dir = Path::new(output_dir).join(format!("tree_{}", index));
        fs::create_dir_all(&tree_output_dir).expect("Failed to create tree output directory");

        let pairs = NewickParser::parse(Rule::newick, &tree).unwrap_or_else(|e| panic!("fail"));
        for pair in pairs {
            number_of_trees += 1;
            let mut vec_trees = newick_to_tree(pair);
            let mut example_tree = vec_trees.pop().unwrap();
            give_depth(&mut example_tree, 0.0);

            let mut flat_tree = Vec::new();
            node_to_flat(&example_tree, &mut flat_tree, None);

            let depths = make_subdivision(&mut flat_tree);
            let intervals = make_intervals(&depths);
            let contemporaneity = find_contemporaneity(&mut flat_tree, &depths);
            let n_species = number_of_species(&contemporaneity);
            let cdf = make_CDF(intervals, n_species);

            let mut copied_flat_tree = flat_tree.clone();
            create_many_genes(
                &mut copied_flat_tree,
                n_transfers.clone(),
                &contemporaneity,
                &cdf,
                &depths,
                &tree_output_dir.to_string_lossy().to_string(),
                &mut rng
            );
        }
    }

    let duration = start.elapsed();
    //println!("Number of trees: {}", number_of_trees);
    //println!("Time taken: {} seconds", duration.as_seconds_f64());
}
