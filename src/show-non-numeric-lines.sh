#!/bin/sh

# Prints all lines consisting of a single numeric.
grep -E -v '^[-+]?([0-9]+\.?|\.[0-9])[0-9]*([eE][-+]?[0-9]+)?$' $1



 