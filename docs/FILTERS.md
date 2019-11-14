# Filtering Polynomials.

## Where do I add filters?
Add filters by modifying your '.in' file.

## What is the syntax for adding filters?
Add lines to your '.in' file using the syntax described below. The lines may be anywhere
in your '.in' file as long as each filter is on its own line and there is nothing else on
those lines.

Filter syntax:

`add_filter['arg1', 'arg2', 'arg3', ...]`

Each arg must be wrapped in single quotes, even if it is a number.
The args must combine into a valid statement to be a valid filter. A simple
statement has syntax as such: <br>
`'nameOfFilter', 'arg1', 'arg2', ...`

Simple or compound statements can be combined into compound statements by using
one of the conjunctions 'and' or 'or'. Example:" <br>
`'nameOfFilter', 'arg1', 'arg2', 'and', 'nameOfFilter', 'arg1', arg2'`

You may wrap simple or complex statements in parenthesis to group conjunctions: <br>
`'nameOfFilter', 'arg1', 'arg2', 'and', '(', 'nameOfFilter', 'arg1', 'or', 'nameOfFilter', 'arg1', ')'`

Using multiple conjunctions without parenthesis leads to the statements being evaluated
from left to right: <br>
So this: `'nameOfFilter', 'arg1', 'arg2', 'and', 'nameOfFilter', 'arg1', 'or', 'nameOfFilter', 'arg1'` <br>
Is the same as: `'(', 'nameOfFilter', 'arg1', 'arg2', 'and' 'nameOfFilter', 'arg1', ')', 'or', 'nameOfFilter', 'arg1'` <br>
But distict from: `'nameOfFilter', 'arg1', 'arg2', 'and', '(', 'nameOfFilter', 'arg1', 'or', 'nameOfFilter', 'arg1', ')'`

You may prepend 'not' before a statement to invert it: <br>
`'not', 'nameOfFilter', 'arg1', 'arg2', ...`

Unless you use parenthesis, 'not' will affect all following statements:
So this: `'not', 'nameOfFilter', 'arg1', 'arg2', 'and', 'nameOfFilter', 'arg1', 'arg2'` <br>
Is the same as: `'not', '(', 'nameOfFilter', 'arg1', 'arg2', 'and', 'nameOfFilter', 'arg1', 'arg2', ')'` <br>
But disict from: `'(', 'not', 'nameOfFilter', 'arg1', 'arg2', ')', 'and', 'nameOfFilter', 'arg1', 'arg2'`


## What filters are available?

There are currently 4 filters implemented.

Each filter matches a certain set of monomials, each monomial that matches at least
one filter will be removed from the polynomial.

#### Argument Syntax
The filters below require several arguments. Each of these arguments is
a string. Format them like this:

'variable-string': 'x-$a-$b+$c-$d/x-$a-$b+$c-$d/...'
* $a is:
    * 'intra' to match only intramolecular variables.
    * 'inter' to match only intermolecular variables.
    * '*' to match both intramolecular and intermolecular variables.
* $b and $c are each:
    * An atom type to match only variables that are an interaction between this atom type
and the other.
    * '*' to match all variables regardless of this atom type.
    * The order of $b and $c does not matter. $b and $c may be both atom types, both wildcards,
or one of each.
* $c is:
    * This argument refers the the 'fragmentation divergence level'
of a variable. This is equal to the number of fragmentation levels
for which the two atoms in the variable are in the same fragments
before the first layer where they are in different fragments. For
example, in a monomer with no sub-fragments, 'fragmentation divergence' 
of all variables will be 1. In a dimer with no sub-fragments, it
will be 1 for those atoms in the same monomer, and 0 for those in different
monomers. For a dimer with 1 layer of sub-fragmentation in each
monomer, it will be 0 for those atoms in different monomers, 1
for those in the same monomer but different subfragments, and 2
for those in the same monomer and same subfragment.
    * '*' to match all fragmentation divergence levels.
    * 'y' where y is a single number to match levels of that number.
    * 'y+' where y is a single number to match levels of that number or greater.
    * 'y-' where y is a single number to match levels of that number or less.
    * 'y-z' where y and z are single numbers to match levels in the range [y, z]
* You may specify multiple variable-strings by delineating with a '/'. Variables
will be matched as long as they match at least one of them.

'degree-string': '$a/$b/$c/...'
* $a, $b, and $c are each:
    * '*' to match all degrees.
    * 'y' where y is a single number to match degrees of that number.
    * 'y+' where y is a single number to match degrees of that number or greater.
    * 'y-' where y is a single number to match degrees of that number or less.
    * 'y-z' where y and z are single numbers to match degrees in the range [y, z]
(Both inclusive).
    * You may specify multiple degree-strings by delineating with a '/'. Degrees
will be matched as long as they match at least one of them.

'fragments-string': '$a/$b/$c/...'
* $a, $b, and $c are each:
    * '*' to match any number of unique fragment counts.
    * 'y' where y is a single number to match unique fragments counts of that number.
    * 'y+' where y is a single number to match unique fragments counts of that number or greater.
    * 'y-' where y is a single number to match unique fragments counts of that number or less.
    * 'y-z' where y and z are single numbers to match unique fragments counts in the range [y, z]
(Both inclusive).
    * You may specify multiple fragments-strings by delineating with a '/'. Unique fragments
will be matched as long as they match at least one of them.

'term-string': '$a/$b/$c/...'
* $a, $b, and $c are each:
    * '*' to match any total degree.
    * 'y' where y is a single number to match total degrees of that number.
    * 'y+' where y is a single number to match total degrees of that number or greater.
    * 'y-' where y is a single number to match total degrees of that number or less.
    * 'y-z' where y and z are single numbers to match total degrees in the range [y, z]
(Both inclusive).
    * You may specify multiple term-strings by delineating with a '/'. Total degrees
will be matched as long as they match at least one of them.

####Individual Degree Filter

This filter filters out monomials based on the degree of one or more variables.
The degree of each variable is looked at independently: they are not summed.

The monomial will be filtered out if the degree of any variable that matches 'variable-string'
matches 'degree-string'.

The syntax for this filter is as follows:
`'ind-degree', 'variable-string', 'degree-string'`

####Sum Degree Filter

This filter filters out monomials based on the degree of one or more variables.
The degrees of all variables are summed, and then that number is used to filter
the monomials.

The monomial will be filtered out if the sum of the degrees of all variables that match
'variable-string' matches 'degree-string'.

The syntax for this filter is as follows:
`'sum-degree', 'variable-string', 'degree-string'`

####Num Fragments Filter

This filter filters out monomials based on the number of unique fragments
participating in one or more variables in the monomial.

The monomial will be filtered out if the number of unique fragments between all
variables of degree 1 or more that match 'variable-string' in this monomial matches 'fragments-string'.

The syntax for this filter is as follows:
`'num-fragments', 'variable-string', 'fragments-string'`

####Degree Filter

This filter is like the Individual Degree Filter but has an additional argument.
In addition to the individual degrees of variables being checked, the total degree
of the monomial is also checked.

The monomial will be filtered out if the degree of any variable that matches 'variable-string'
matches 'degree-string' and the total degree of all variables in the monomial matches
'term-string'.

The syntax for this filter is as follows:
`'degree', 'variable-string', 'degree-string', 'term-string'`

## Example Filters

Here are some example filters and descriptions of what they do.

`add_filter['ind-degree', '*', '2']` <br>
Filters out all monomials that have degree two in at least one variable.

`add_filter['sum-degree', 'x-intra-*+*-*', '2+']`
Filters out all monomials that have sum degree two or more between all intramolecular variables.

`add_filter['num-fragments', 'x-inter-*+*-*', '2-']`
Filters out all monomials that have 2 or fewer unique fragments between all intermolecular variables.

`add_filter['not', 'ind-degree', 'x-inter-*+*-*', '2+']`
Filters out all monomials that do not have at least degree 2 in at least one intermolecular variable.

`add_filter['degree', 'x-intra-A+*-*', '3/4', '5']`
Filters out all degree 5 monomials that have degree 3 or 4 in at least one intramolecular variable involving
an atom of type A.

Note: 'degree' filters can always be replicated using a 'sum-degree' filter combined with
an 'ind-degree' filter with an 'and'. For example: <br>
This: `add_filter['degree', 'x-intra-A+*-*', '3/4', '5']` <br>
Is the same as: `add_filter['sum-degree', '*', '5', and, 'ind-degree', 'x-intra-A+*-*', '3/4']` <br>
Because of this, 'degree' filters may be depricated at some point in the future.

#### mbpol 2b filters:

`add_filter['sum-degree', '*', '1', 'and', 'sum-degree', 'x-intra-*+*-*', '1+']` <br>
Filters out all degree 1 monomials that have any intramolecular variables. <br>
`add_filter['sum-degree', '*', '2', 'and', 'sum-degree', 'x-intra-*+*-*', '2']` <br>
Filters out all degree 2 monomials that have sum degree 2 between all intramolecular variables. 
This leaves all degree 2 monomials that are at most linear in intramolecular variables. <br>
`add_filter['sum-degree', '*', '3', 'and', 'sum-degree', 'x-intra-*+*-*', '3']` <br>
Filters out all degree 3 monomials that have sum degree 3 between all intramolecular variables.
This leaves all degree 3 monomials that are at most quadratic in intramolecular variables. <br>
`add_filter['sum-degree', '*', '4', 'and', 'sum-degree', 'x-intra-*+*-*', '1-/3+']` <br>
Filters out all degree 4 monomials that have sum degree of 1 or less or 3 or more between all intramolecular variables. <br>

####mbpol 3b filters:

`add_filter[num-fragments', 'x-inter-*+*-*', '2-']` <br>
Filters out all monomials that have two or fewer unique fragments between all intermolecular variables. <br>

`add_filter['sum-degree', '*', '2', 'and', sum-degree', 'x-intra-A+B-*', '1+']` <br>
Filters out all degree 2 monomials with any contribution from a x-intra-A+B variable. <br>
`add_filter['sum-degree', '*', '2', 'and', sum-degree', 'x-intra-B+B-*', '1+']` <br>
Filters out all degree 2 monomials with any contribution from a x-intra-B+B variable. <br>

`add_filter['sum-degree', '*', '3', 'and', sum-degree', 'x-intra-A+B-*', '2+']` <br>
Filters out all degree 3 monomials with sum degree 2 or more between all x-intra-A+B variables. <br>
`add_filter['sum-degree', '*', '3', 'and', sum-degree', 'x-intra-B+B-*', '2+']` <br>
Filters out all degree 3 monomials with sum degree 2 or more between all x-intra-B+B variables. <br>

`add_filter['sum-degree', '*', '4', 'and', sum-degree', 'x-intra-A+B-*', '2+']` <br>
Filters out all degree 4 monomials with sum degree 2 or more between all x-intra-A+B variables. <br>
`add_filter['sum-degree', '*', '4', 'and', sum-degree', 'x-intra-B+B-*', '2+']` <br>
Filters out all degree 4 monomials with sum degree 2 or more between all x-intra-B+B variables. <br>
`add_filter['sum-degree', '*', '4', 'and', sum-degree', 'x-inter-A+A-*', '2+']` <br>
Filters out all degree 4 monomials with sum degree 2 or more between all x-inter-A+A variables. <br>
`add_filter['sum-degree', '*', '4', 'and', sum-degree', 'x-inter-B+B-*', '2+']` <br>
Filters out all degree 4 monomials with sum degree 2 or more between all x-inter-B+B variables. <br>
