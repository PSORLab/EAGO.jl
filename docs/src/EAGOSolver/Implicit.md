## Solving problems with implicitly defined functions

For functions that can be directly rearranged, it's recommended that the user
go through the MathProg base interface and define functions explicitly as
relaxing multiplication can be expensive. If no such functions are present, using
the JuMP interface is recommended.

The implicit solver option embeds the implicit problem in a solver and in turn a model.

## Usage implicit solver
