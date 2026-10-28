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
