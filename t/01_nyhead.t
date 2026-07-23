use strict;
use warnings;
use Test::More;
use FindBin;
use PDL;
use PDL::IO::NYHead;

my $mat = "$FindBin::Bin/data/mini_nyhead.mat";
-r $mat or plan skip_all => "fixture $mat not readable";

# h5dump は cell-of-strings(電極名/HO野名)の読み取りにのみ必要。
sub have_h5dump {
    for my $d (split /:/, ($ENV{PATH} // '')) {
        return 1 if -x "$d/h5dump";
    }
    return 0;
}
my $h5dump = have_h5dump();

my $ny = PDL::IO::NYHead->new($mat);
isa_ok($ny, 'PDL::IO::NYHead');

# ---- 数値系(PDL::IO::HDF5 のみ。h5dump 不要) ------------------------------
my $pos = $ny->electrode_pos;
is(join(",", $pos->dims), "4,3",            "electrode_pos is (Ne,3)");
is($ny->n_electrodes, 4,                    "n_electrodes == 4");
ok($pos->at(0,0)==10 && $pos->at(0,2)==2,   "electrode 0 coords");

my $L = $ny->leadfield;
is(join(",", $L->dims), "4,6",              "leadfield is (Ne,Nv)");
is(join(",", $L->slice("(0),")->list), "0,4,8,12,16,20", "leadfield electrode 0 row");

my $s = $ny->surface('cortex75K');
is(join(",", $s->{vc}->dims), "6,3",        "cortex75K vc is (Nv,3)");
is($s->{tri}->at(0,0), 1,                   "tri is 1-based");
is($s->{tri0}->at(0,0), 0,                  "tri0 is 0-based");

my $c2 = $ny->cortex('cortex2K');
is(join(",", $c2->{in_from}->list), "0,2,4", "cortex2K in_from is 0-based");
is(join(",", $c2->{vc}->slice(",(0)")->list), "0,20,40",
                                            "cortex2K vc reconstructed via fancy index");

my $ai = $ny->atlas_index;
is(join(",", $ai->list), "1,1,2,2,3,0",     "atlas_index (1-based, 0=none)");

# nearest_vertex: 電極0 (10,8,2) の最近傍 head 頂点は index 2 (10,2,9)
my $iv = $ny->nearest_vertex($ny->electrode_pos, 'head');
is(join(",", $iv->dims), "4",               "nearest_vertex returns (Ne)");
is($iv->at(0), 2,                           "electrode 0 nearest head vertex == 2");
is($ny->nearest_vertex(pdl(10,2,9), 'head'), 2, "single-point nearest_vertex scalar");

is(join(",", $ny->idx19->slice("0:2")->list), "156,157,119", "IDX19 constant frozen");

# ---- cell-of-strings 系(h5dump 必要) -------------------------------------
SKIP: {
    skip "h5dump not found in PATH (needed for label fields)", 4 unless $h5dump;

    my $lab = $ny->electrode_labels;
    is(scalar(@$lab), 4,                    "electrode_labels count");
    is($lab->[0], "Fp1",                    "electrode_labels[0] == Fp1");

    my $ho = $ny->ho_labels;
    is($ho->[1], "Right Insular Cortex",    "ho_labels[1] round-trips exactly");
    is(($ny->area_of_vertex(0) // ''), "Left Frontal Pole",
                                            "area_of_vertex(0) maps in_HO->HO_labels");
}

done_testing();
