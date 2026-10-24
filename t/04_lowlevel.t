use strict;
use warnings;
use Test::More;
use FindBin;
use PDL;
use PDL::IO::NYHead;
use PDL::IO::NYHead::H5Dump;

# cell-of-strings を PDL::IO::HDF5 の低水準関数で読む(h5dump 不要)。
# 期待値はフィクスチャを作ったときの文字列そのもの(h5dump の出力とは独立)。
my $mat = "$FindBin::Bin/data/mini_nyhead.mat";
-r $mat or plan skip_all => "fixture $mat not readable";

my @clab = qw(Fp1 Fp2 Cz Oz);
my @ho   = ("Left Frontal Pole", "Right Insular Cortex", "Subcortical");

my $ll = PDL::IO::NYHead::H5Dump::read_cell_strings_lowlevel($mat, '/sa/clab_electrodes');
is_deeply($ll, \@clab, "lowlevel: clab_electrodes, all cells");
is_deeply(PDL::IO::NYHead::H5Dump::read_cell_strings_lowlevel($mat, '/sa/HO_labels'), \@ho,
          "lowlevel: HO_labels, all cells (spaces kept, not concatenated)");

# h5dump が見つからない場所を指定すると auto は lowlevel を選ぶ
my $r = PDL::IO::NYHead::H5Dump->new(h5dump => '/nonexistent/h5dump');
is($r->method, 'lowlevel', "auto falls back to lowlevel without h5dump");
is_deeply($r->read_cell_strings($mat, '/sa/clab_electrodes'), \@clab, "fallback read");

# PDL::IO::NYHead 経由(開いている HDF5 オブジェクトを使う)
my $ny = PDL::IO::NYHead->new($mat, strings => 'lowlevel');
is($ny->strings_method, 'lowlevel', "strings => lowlevel");
is_deeply($ny->electrode_labels, \@clab, "electrode_labels via lowlevel");
is_deeply($ny->ho_labels, \@ho, "ho_labels via lowlevel");
is(($ny->area_of_vertex(0) // ''), "Left Frontal Pole", "area_of_vertex via lowlevel");

my $ny2 = PDL::IO::NYHead->new($mat, h5dump => '/nonexistent/h5dump');
is($ny2->strings_method, 'lowlevel', "NYHead auto without h5dump");
is_deeply($ny2->electrode_labels, \@clab, "electrode_labels, auto without h5dump");

eval { PDL::IO::NYHead::H5Dump->new(method => 'xyz') };
like($@, qr/method must be/, "bad method croaks");

# 否定対照: 数値のデータセットは参照の配列ではないので読めない
eval { PDL::IO::NYHead::H5Dump::read_cell_strings_lowlevel($mat, '/sa/locs_3D_orig') };
ok($@, "numeric dataset is rejected");

SKIP: {
    my $d = PDL::IO::NYHead::H5Dump->new(method => 'h5dump');
    skip "h5dump not found in PATH", 2 unless $d->have_h5dump;
    is_deeply($d->read_cell_strings($mat, '/sa/clab_electrodes'), $ll, "h5dump == lowlevel (clab)");
    is_deeply($d->read_cell_strings($mat, '/sa/HO_labels'), \@ho,       "h5dump == expected (HO)");
}

done_testing();
