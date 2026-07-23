use strict; use warnings;
use Test::More;
use FindBin;
use File::Temp qw(tempdir);
use PDL;
use PDL::IO::NYHead;

my $mat = "$FindBin::Bin/data/mini_nyhead.mat";
-r $mat or plan skip_all => "fixture $mat not readable";
my $dir = tempdir(CLEANUP => 1);
my $ny  = PDL::IO::NYHead->new($mat);

# mesh() は head / cortex75K / 低解像で {vc,tri(1-based)} を返す
my $m = $ny->mesh('head');
is(join(",", $m->{vc}->dims), "5,3", "mesh(head) vc (Nv,3)");
is($m->{tri}->min, 1,                "mesh tri is 1-based");

# OBJ: 頂点数=v行, 面数=f行, 面 index は 1..Nv
my $obj = $ny->to_obj('head', "$dir/head.obj");
open my $of, '<', $obj or die;
my @lines = <$of>; close $of;
my @v = grep { /^v / } @lines;
my @f = grep { /^f / } @lines;
is(scalar(@v), $m->{vc}->dim(0),   "OBJ vertex count == vc rows");
is(scalar(@f), $m->{tri}->dim(0),  "OBJ face count == tri rows");
my @f0 = (split /\s+/, $f[0])[1..3];
ok((grep { $_ >= 1 && $_ <= $m->{vc}->dim(0) } @f0) == 3, "OBJ faces 1-based in range");

# USDA: faceVertexIndices は 0-based で OBJ-1 と一致
my $usda = $ny->to_usda('head', "$dir/head.usda");
open my $uf, '<', $usda or die; local $/; my $txt = <$uf>; close $uf;
like($txt, qr/#usda 1\.0/,                    "USDA header");
my ($cnt) = $txt =~ /faceVertexCounts = \[([^\]]*)\]/;
is($cnt, join(",", (3) x $m->{tri}->dim(0)),  "faceVertexCounts all 3");
my ($idx) = $txt =~ /faceVertexIndices = \[([^\]]*)\]/;
my @ui = split /,/, $idx;
is_deeply([@ui[0..2]], [map { $_ - 1 } @f0],  "USDA idx == OBJ face - 1 (0-based)");

# write_mesh 分岐
ok($ny->write_mesh('cortex75K', 'obj', "$dir/ctx.obj"),  "write_mesh obj");
ok(-s "$dir/ctx.obj",                                     "cortex75K obj non-empty");
eval { $ny->write_mesh('head', 'stl', "$dir/x") };
like($@, qr/unknown format/,                              "write_mesh rejects bad format");

done_testing();
