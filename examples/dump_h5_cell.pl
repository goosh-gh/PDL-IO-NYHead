#!/usr/bin/env perl
# dump_h5_cell.pl — MATLAB cell-of-strings を HDF5(.mat v7.3)から取り出して
# 1 行 1 名で標準出力に書く。PDL::IO::NYHead::H5Dump のラッパ。
#
# 例:
#   perl tools/dump_h5_cell.pl sa_nyhead.mat /sa/HO_labels > HO_labels.txt
#   perl tools/dump_h5_cell.pl sa_nyhead.mat /sa/clab_electrodes > clab.txt
#   perl tools/dump_h5_cell.pl --index sa_nyhead.mat /sa/HO_labels   # "IDX : NAME" 形式
#
# 依存: h5dump(HDF5 tools)。--h5dump で実行ファイルを指定可。
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Getopt::Long;
use PDL::IO::NYHead::H5Dump;

my ($index, $h5dump) = (0, 'h5dump');
GetOptions('index' => \$index, 'h5dump=s' => \$h5dump) or die "bad args\n";
my ($file, $dataset) = @ARGV;
defined $dataset or die "usage: $0 [--index] [--h5dump PATH] <file.mat> <dataset>\n";

my $labels = PDL::IO::NYHead::H5Dump->new(h5dump => $h5dump)
    ->read_cell_strings($file, $dataset);

for my $i (0 .. $#$labels) {
    if ($index) { printf "%d : %s\n", $i, $labels->[$i]; }   # "IDX : NAME"
    else        { print  $labels->[$i], "\n"; }              # 1 行 1 名
}
warn sprintf "dumped %d strings from %s\n", scalar(@$labels), $dataset;
