"""Length/similarity heuristic flagging candidate nested-plasmid relationships.

This is NOT alignment-based containment detection (unlike MobMess's MUMmer4 approach).
Symmetric Mash/k-mer distance is a poor proxy for size-disparate containment: Jaccard
distance grows with size disparity even under perfect containment, so this heuristic can
only see near-equal-length, high-similarity pairs -- a narrower "structural variant of the
same backbone" signal, not general containment. Size-disparate containment would need
sequence alignment (MUMmer4/nucmer), which is out of this project's stdlib-only scope.
"""

CONTAINMENT_FIELDS = ("small_plasmid_id large_plasmid_id small_isolate_id large_isolate_id "
                      "small_length large_length length_ratio pairwise_distance "
                      "heuristic_flag interpretation").split()

INTERPRETATION = ("length/similarity heuristic only; NOT alignment-confirmed containment; "
                  "may miss size-disparate containment (see docs)")


def detect_containment(accepted, ids, pos, distances, min_length_ratio=0.5, max_length_ratio=0.95, max_distance=0.05):
    if not 0 < min_length_ratio < max_length_ratio <= 1:
        raise ValueError("Require 0 < containment-min-ratio < containment-max-ratio <= 1")
    if not 0 <= max_distance <= 1:
        raise ValueError("containment-max-distance must be in [0,1]")
    rows = []
    for i in range(len(accepted)):
        for j in range(i + 1, len(accepted)):
            a, b = accepted[i], accepted[j]
            len_a, len_b = int(a["length"]), int(b["length"])
            small, large = (a, b) if len_a <= len_b else (b, a)
            small_len, large_len = min(len_a, len_b), max(len_a, len_b)
            ratio = small_len / large_len
            if not min_length_ratio <= ratio <= max_length_ratio:
                continue
            distance = distances[pos[small["plasmid_id"]]][pos[large["plasmid_id"]]]
            if distance > max_distance:
                continue
            rows.append(dict(small_plasmid_id=small["plasmid_id"], large_plasmid_id=large["plasmid_id"],
                             small_isolate_id=small["isolate_id"], large_isolate_id=large["isolate_id"],
                             small_length=small_len, large_length=large_len, length_ratio=ratio,
                             pairwise_distance=distance, heuristic_flag="heuristic_length_similarity",
                             interpretation=INTERPRETATION))
    return rows
