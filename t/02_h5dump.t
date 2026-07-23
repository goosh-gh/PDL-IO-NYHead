use strict;
use warnings;
use Test::More;
use FindBin;
use PDL::IO::NYHead::H5Dump;

sub have_h5dump {
    for my $d (split /:/, ($ENV{PATH} // '')) { return 1 if -x "$d/h5dump"; }
    return 0;
}
plan skip_all => "h5dump not found in PATH" unless have_h5dump();

my $mat = "$FindBin::Bin/data/mini_nyhead.mat";
-r $mat or plan skip_all => "fixture $mat not readable";

my $h5 = PDL::IO::NYHead::H5Dump->new;
isa_ok($h5, 'PDL::IO::NYHead::H5Dump');

# MATLAB cell-of-strings を #refs# object-reference から復元。
# 区切りは版非依存(DATASET の objid 有無に依存しない)ことを、
# 複数セルが1本に連結されないことで確認する。
my $ho = $h5->read_cell_strings($mat, '/sa/HO_labels');
is(scalar(@$ho), 3,                       "3 cells, not concatenated into 1");
is($ho->[0], "Left Frontal Pole",         "cell 0 exact");
is($ho->[1], "Right Insular Cortex",      "cell 1 exact");
is($ho->[2], "Subcortical",               "cell 2 exact");

my $clab = $h5->read_cell_strings($mat, '/sa/clab_electrodes');
is(scalar(@$clab), 4,                     "clab_electrodes count");
is($clab->[0], "Fp1",                     "clab_electrodes[0]");

done_testing();
