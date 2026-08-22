use strict; use warnings; use PDL; use lib '/Users/goosh/src/PDL_IO_NYHead/lib'; use PDL::IO::NYHead;
my $ny = PDL::IO::NYHead->new('sa_nyhead.mat');

my $target = pdl([-77.2, -24.2, -59.1]);   # snapped head 頂点(MNI mm)

my $pos    = $ny->electrode_pos;      # (231,3) MNI mm
my $labels = $ny->electrode_labels;   # 231 ラベル(配列参照)

my $diff = $pos - $target->dummy(0);              # (231,3)
my $d    = sqrt( ($diff**2)->mv(1,0)->sumover );  # 座標軸(dim1)を潰す → (231)
my $ord  = $d->qsorti;                            # 近い順 index

printf "target MNI (%.1f, %.1f, %.1f)\n\n", $target->list;
printf "%-8s %8s   %s\n", "elec", "dist_mm", "MNI(x,y,z)";
for my $n (0..7) {
    my $i = $ord->at($n);
    my @p = $pos->slice("($i),:")->list;
    printf "%-8s %8.1f   (%.1f, %.1f, %.1f)\n", $labels->[$i], $d->at($i), @p;
}
