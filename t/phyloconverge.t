# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl phyloconverge.t'

#########################

use strict;
use warnings;

use Test::More tests => 1;

#########################

my $cmd = "perl phyloconverge --treefile=examples/converged.tre --tabfile=examples/depths.txt > /dev/null ";
my $test = system $cmd;
is($test, 0, "was able to run example 1");

