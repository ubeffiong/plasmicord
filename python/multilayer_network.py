"""Optional multilayer network export: one edge per (isolate pair, plasmid unit).

This relabels edges PlasmiCord already computes by two categorical fields it already has
(plasmid unit, chromosomal cluster) -- each plasmid unit is naturally one layer, and every
edge also carries whether it stays within or crosses a known chromosomal cluster. It is
intentionally lightweight and stdlib-only, not a true node-multiplex/projection framework
(e.g. pymnet/multinet); the default network.graphml/network.edges.tsv (which aggregate
across shared units) remain the primary output.
"""
from pathlib import Path
from .build_network import escape

MULTILAYER_EDGE_FIELDS = "source target plasmid_unit cluster_relation minimum_distance threshold_margin".split()


def cluster_relation(meta_by_iso, source, target):
    a = meta_by_iso.get(source, {}).get('chromosomal_cluster', '')
    b = meta_by_iso.get(target, {}).get('chromosomal_cluster', '')
    if not a or not b:
        return 'unknown_cluster'
    return 'same_cluster' if a == b else 'cross_cluster'


def multilayer_edges(edge_details, meta_by_iso):
    return [dict(source=e['source'], target=e['target'], plasmid_unit=e['plasmid_unit'],
                 cluster_relation=cluster_relation(meta_by_iso, e['source'], e['target']),
                 minimum_distance=e['minimum_distance'], threshold_margin=e['threshold_margin'])
            for e in edge_details]


def write_multilayer_graphml(path, edges):
    lines = ['<?xml version="1.0" encoding="UTF-8"?>',
             '<graphml xmlns="http://graphml.graphdrawing.org/xmlns">',
             '  <key id="e_layer" for="edge" attr.name="plasmid_unit" attr.type="string"/>',
             '  <key id="e_relation" for="edge" attr.name="cluster_relation" attr.type="string"/>',
             '  <key id="e_distance" for="edge" attr.name="minimum_distance" attr.type="double"/>',
             '  <key id="e_margin" for="edge" attr.name="threshold_margin" attr.type="double"/>',
             '  <graph edgedefault="undirected">']
    for node in sorted({n for e in edges for n in (e['source'], e['target'])}):
        lines.append(f'    <node id="{escape(node)}"/>')
    for i, e in enumerate(edges):
        lines.append(f'    <edge id="e{i}" source="{escape(e["source"])}" target="{escape(e["target"])}">')
        lines.append(f'      <data key="e_layer">{escape(e["plasmid_unit"])}</data>')
        lines.append(f'      <data key="e_relation">{escape(e["cluster_relation"])}</data>')
        lines.append(f'      <data key="e_distance">{e["minimum_distance"]}</data>')
        lines.append(f'      <data key="e_margin">{e["threshold_margin"]}</data>')
        lines.append('    </edge>')
    lines.append('  </graph>')
    lines.append('</graphml>')
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
