
/*
1. Define the data structures used to represent the tree parenthood information
2. Use the crate pest to convert the dataframe into a newich string.



3. We can use SORTING on vectors, and then use a binary search to find the element corresponding to the drawn random number
4. Whenever there are two transfers between two nodes, we can just keep the last one.
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
use time::Instant;
extern crate regex;
use regex::Regex;
use rand::{thread_rng, Rng};
// ATTENTION : les enfants et parents d'un noeud devront être mutables, ainsi que la longueur.

#[derive(Parser)]
#[grammar = "newick.pest"]
struct NewickParser;

struct Node {
    name: String,
    left_child: Option<Box<Node>>,
    right_child: Option<Box<Node>>,
    parent: Option<usize>, // Using index for parent
    depth: Option<f64>,
    length: f64,
}


struct Tree {
    nodes: Vec<Node>,
    root: Option<usize>,
}
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


use pest::error::Error;

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

fn node_to_newick(node: &Node) -> String {
    /* Takes a node and returns the corresponding subtree in Newick format.
        ARGUMENTS:
            - node: the node to convert to Newick format.
        RETURNS:
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
        give_depth(left_child, depth + node.length);
    }
    if let Some(right_child) = &mut node.right_child {
        give_depth(right_child, depth + node.length);
    }
}

fn make_subdivision(flat_tree: &mut Vec<FlatNode>) -> Vec<DurationAndContemporaneity> {
    /* Input: a flat tree.
    Output: the time subdivision induced by the nodes of the tree.
    Example : If the original nwk tree was ((A:1,B:2)C:1,D:5)R:0, the subdivision should be [0,1,2,3,5]
    */
    let mut depths: Vec<f64> = flat_tree.iter().filter_map(|node| node.depth).collect();
    depths.sort_by(|a, b| a.partial_cmp(b).unwrap());
    depths.dedup();
    return depths;
}

fn make_intervals(depths: Vec<f64>) -> Vec<f64> {
    /* Input: a vector of depths.
    Output: a vector of intervals, where the i-th interval is the difference between the i-th and (i+1)-th depth.
    Adds a 0 interval at the beginning so the size of the intervals vector is the same as the size of the depths vector.
    Example: if the depths are [0,1,2,3,5], the intervals are [0,1,1,1,2].
    */
    let mut intervals:Vec<f64> = Vec::new();
    intervals.push(0.0 as f64)
    for i in 0..depths.len()-1 {
        intervals.push(depths[i+1] - depths[i]);
    }
    return intervals;
}

fn find_closest_index(vec: &[f64], value: f64) -> usize {
    match vec.binary_search(&value) {
        Ok(idx) => idx,
        Err(idx) => {
            if idx == 0 {
                0
            } else if idx == vec.len() {
                vec.len() - 1
            } else {
                // Determine which of vec[idx - 1] or vec[idx] is closer to the value
                if (value - vec[idx - 1]).abs() < (value - vec[idx]).abs() {
                    idx - 1
                } else {
                    idx
                }
            }
        }
    }
}




fn find_contemporaneity(flat_tree: &mut Vec<FlatNode>, depths: Vec<f64>) -> Vec<Vec<usize>> {
    // Returns a vector of vectors, where each vector contains the indices of the nodes that are contemporaneous.
    // Quadratic time complexity for now, maybe we can improve.
    let mut contemporaneity: Vec<Vec<f64>> = vec![Vec::new(); depths.len()];
    for i in 0..flat_tree.len() {
        // Attention à la racine ! L'indice du début et l'indice de la fin pour la racine devraient être identiques...
        let start_time = flat_tree[i].depth - flat_tree[i].length;
        let end_time = flat_tree[i].depth;
        let start_index = find_closest_index(&depths, start_time);
        let end_index = find_closest_index(&depths, end_time);
        for j in start_index+1..end_index {
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
    CDF.push(0.0 as f64);
    for i in 0..n {
        CDF.push(CDF[i] + intervals[i]*number_of_species[i]);
    }
    let total_value = CDF[n];
    for i in 0..n+1 {
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
    let r: f64 = f64::from(rng.gen());
    let index = match CDF.binary_search_by(|&probe| probe.partial_cmp(&u).unwrap_or(std::cmp::Ordering::Less)) {
        Ok(index) => {
            // Found an exact match, use `index`.
        },
        Err(index) => {
            // No exact match, but `index` tells you where this number could be inserted while maintaining the order.
        }
    }

    let depth = (r - CDF[index-1])/(CDF[index] - CDF[index-1]) * (depths[index] - depths[index-1]) + depths[index-1];
    return (depth, index);
}


fn main() {
    // Read command line arguments
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
            let example_tree = vec_trees.pop().unwrap();

            // Convert Node representation to FlatNode representation
            let mut flat_tree = Vec::new();
            node_to_flat(&example_tree, &mut flat_tree, None);
            
            println!("Length of first flatnode: {}", flat_tree[0].length);
            println!("Total length of flat tree: {}", total_length_of_flat_tree(&flat_tree));
            println!("Number of nodes in flat tree: {}", flat_tree.len());
            // Convert FlatNode representation back to Node representation
            let reconstructed_tree = flat_to_node(&flat_tree, 0, None).unwrap();

            // Convert Node representation back to Newick format and print
            let reconstructed_newick = node_to_newick(&reconstructed_tree) + ";";
        }
    }

    let duration = start.elapsed();
    println!("Number of trees: {}", number_of_trees);
    println!("Time taken: {} seconds and {} nanoseconds", duration.as_seconds_f64(), duration.subsec_nanoseconds());
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