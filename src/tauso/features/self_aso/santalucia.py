"""The terms the self-structure search needs beyond stacking.

`weights` covers stacking inside a helix. It has no helix end free of fraying and no loop, so
helix initiation and the hairpin loop penalty come from SantaLucia as well, on the same
published kcal/mol scale.

Both are chemistry-blind: a 2'-MOE loop is charged as a DNA loop. Nothing better exists for
the modified sugars, and the alternative -- Turner RNA parameters -- is calibrated to be used
with a terminal mismatch term this search does not have.
"""

from ..hybridization.weights.dna import DNA_HAIRPIN_LOOP_PENALTY, INIT_AT, INIT_GC

HAIRPIN_LOOP_PENALTY = DNA_HAIRPIN_LOOP_PENALTY
"""Cost of closing a hairpin, by the number of unpaired bases in the loop."""

MIN_LOOP = 3
"""A loop needs three bases; the backbone cannot turn tighter."""

MAX_LOOP = len(HAIRPIN_LOOP_PENALTY) - 1
"""Longest loop the table covers. Follows the table so the two cannot drift apart."""

INITIATION_AT = INIT_AT
INITIATION_GC = INIT_GC
"""Helix initiation, charged once per free end. A hairpin has one; a duplex has two."""
