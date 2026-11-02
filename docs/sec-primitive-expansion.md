# SEC Primitive Expansion

The cone builder shared by LEC and SEC turns every combinational Naja
primitive on a path into a Boolean expression. Most cells carry one truth
table per output bit and are expanded directly. Parametric primitives, whose
width would make an explicit table impractical, carry a *generic* table
instead: a type tag plus a little metadata, which the builder expands
structurally. This document lists those generic types and how each is built.

## Generic table types

| Type | Produced for | Expansion |
| --- | --- | --- |
| `AND`, `NAND`, `OR`, `NOR`, `XOR`, `XNOR` | Wide gates | A chain over the inputs, complemented for the negated kinds. |
| `TABLE_SELECT` | Multiplexers and dynamically indexed reads of a signal table | A tree over the address bits, most significant bit first, selecting one data row. Rows at or beyond the table depth read 0. |
| `DIVMOD` | Division and remainder | An unsigned restoring divider; signed operands are negated conditionally on entry and the quotient and remainder on exit. |
| `ADD_SUB` | Addition and subtraction | A ripple-carry adder, or a borrow chain for subtraction. |
| `MULTIPLY` | Multiplication | The low `width` product bits, each partial product added in place. |
| `SHIFT` | Logical and arithmetic shifts | A barrel of multiplexer stages: stage `s` moves the data by `2^s` when count bit `s` is set. Bits shifted in are 0, or the sign bit for an arithmetic right shift, which also covers counts of the width or more. |
| `ROM` | Dynamically indexed constant tables | A `TABLE_SELECT` tree whose rows are constants (section "Read-only memories"). |
| `VECTOR_SELECT` | Dynamically indexed part-select reads, `data[offset +: w]` | A mux tree over the offset bits, most significant bit first, keyed by the offsets that place a data bit at this output bit (section "Dynamic part selects"). |
| `VECTOR_UPDATE` | Dynamically indexed part-select writes, `v[offset +: w] = d` | The same tree keyed by the offsets that place a written bit at this output bit, keeping the base bit otherwise, under the enable. |

A cell with no usable table and no sequential model is opaque, with the reason
`no initialized combinational truth table or usable sequential model`; see the
"Opaque Outputs" section of [sec-sequential-models.md](sec-sequential-models.md).

## Word-level arithmetic

The SystemVerilog frontend lowers RTL arithmetic to word-level primitives,
one instance per operator, rather than to gate-level logic. Each output bit of
such a primitive carries a generic table of type `ADD_SUB`, `MULTIPLY` or
`SHIFT` with the word metadata: the operand width, the count width for a
shift, the index of the output bit, and the subtract, right and arithmetic
flags. The operand bits arrive most significant bit first, the first operand
before the second, as the primitive declares them.

The builder expands the whole word once per instance and caches every output
bit on its iso, so each bit is built once however many readers it has and the
carry, partial-product and shift logic is shared between the bits of one
instance. A table whose arity does not match its metadata is an error, not an
opaque cell.

## Read-only memories

A constant table that RTL indexes dynamically is lowered by Naja to a memory
primitive with read ports but no write port and no clock. Such a memory holds
no state, so SEC keeps it out of the structured-memory model, which assumes
clocked read outputs, and treats each read-data bit as combinational: its
generic table is of type `ROM`, its dependencies are the address bits, most
significant bit first, and the table's metadata column holds that bit of every
row. Rows at or beyond the depth read 0. Memories with write ports keep the
structured-memory model.

## Dynamic part selects

Naja lowers an indexed part-select read `data[offset +: w]` to a
`VECTOR_SELECT` primitive and an indexed write `v[offset +: w] = d` to a
`VECTOR_UPDATE` primitive, which replaces that window of a base vector. Each
output bit carries its own generic table holding the vector width, the window
width, the offset width, the output bit index and whether the offset is
signed. The dependencies arrive in the primitive's term order: for an update
the base bit of this output first, then the data bits, most significant bit
first, then the offset bits, most significant bit first, and for an update the
enable last.

The builder lists the offset values that place a window bit at this output
bit, drops the ones the offset width cannot hold (Naja does not alias offsets
modulo the offset range), and decodes the offset once as a mux tree over its
bits, pruning every prefix no candidate shares. A select whose offset leaves
no data bit at the output reads 0, like a `TABLE_SELECT` row beyond the depth;
Naja's own model calls that bit X. An update with the enable low or with no
window bit at the output keeps the base bit, so a window that runs past the
vector is clipped.
