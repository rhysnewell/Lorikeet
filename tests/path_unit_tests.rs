#![allow(
    non_upper_case_globals,
    non_snake_case
)]

#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate approx;

use lorikeet_genome::graphs::base_edge::{BaseEdge, BaseEdgeStruct};
use lorikeet_genome::graphs::path::Path;
use lorikeet_genome::graphs::seq_graph::SeqGraph;
use lorikeet_genome::graphs::seq_vertex::SeqVertex;

#[test]
fn test_make_path() {
    let mut g = SeqGraph::new(3);
    let v1 = SeqVertex::new(b"a".to_vec());
    let v2 = SeqVertex::new(b"b".to_vec());
    let v3 = SeqVertex::new(b"c".to_vec());
    let v4 = SeqVertex::new(b"d".to_vec());

    let v1_index = g.base_graph.add_node(&v1);
    let v2_index = g.base_graph.add_node(&v2);
    let v3_index = g.base_graph.add_node(&v3);
    let v4_index = g.base_graph.add_node(&v4);

    let e1_index =
        g.base_graph
            .graph
            .add_edge(v1_index, v2_index, BaseEdgeStruct::new(false, 1, 0));
    g.base_graph
        .graph
        .edge_weight_mut(e1_index)
        .unwrap()
        .inc_multiplicity(1);

    let e2_index =
        g.base_graph
            .graph
            .add_edge(v2_index, v3_index, BaseEdgeStruct::new(false, 1, 0));
    let _e3_index =
        g.base_graph
            .graph
            .add_edge(v3_index, v4_index, BaseEdgeStruct::new(false, 1, 0));

    let path = Path::new(v2_index, Vec::new());
    let path1 = path.new_add_edge(e2_index, &g.base_graph);
    let path2 = path1.new_prepend_edge(e1_index, &g.base_graph);

    assert_eq!(path.len(), 0);
    assert_eq!(path1.len(), 1);
    assert_eq!(path2.len(), 2);

    assert!(path2.contains_vertex(v1_index, &g.base_graph));
    assert!(path2.contains_vertex(v2_index, &g.base_graph));
    assert!(path2.contains_vertex(v3_index, &g.base_graph));
    assert!(!path2.contains_vertex(v4_index, &g.base_graph));

    assert!(!path1.contains_vertex(v1_index, &g.base_graph));
    assert!(path1.contains_vertex(v2_index, &g.base_graph));
    assert!(path1.contains_vertex(v3_index, &g.base_graph));
    assert!(!path1.contains_vertex(v4_index, &g.base_graph));

    assert_ne!(path2.get_edges(), path1.get_edges());

    assert_eq!(path.get_vertices(&g.base_graph).len(), 1);
    assert_eq!(path.get_last_vertex(), v2_index);
    assert_eq!(path.get_first_vertex(&g.base_graph), v2_index);
    assert_eq!(path1.get_first_vertex(&g.base_graph), v2_index);
    assert_eq!(path1.get_last_vertex(), v3_index);
    assert_eq!(path.get_bases(&g.base_graph), b"b".to_vec());
    assert_eq!(path2.get_bases(&g.base_graph), b"abc".to_vec())
}
