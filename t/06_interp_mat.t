use strict;
use warnings;
use Test::More;
use FindBin;
use PDL;
use PDL::IO::NYHead;

my $mat = "$FindBin::Bin/data/mini_nyhead.mat";
-r $mat or plan skip_all => "fixture $mat not readable";
my $ny  = PDL::IO::NYHead->new($mat);
my $pos = $ny->electrode_pos->double;                 # (4,3)
my $L   = $ny->leadfield->double;                     # (4,6)

# 電極 2(Cz)の位置ではその行。4 本なので内挿は 4 点を通る 1 次関数そのもの
my ($R, $W) = $ny->interp_rows($pos->slice('2:2,:'));
is(join(',', $R->dims), '1,6', 'rows (np,ns)');
is(join(',', $W->dims), '4,1', 'weights (ne,np)');
ok(($R - $L->slice('2:2,:'))->abs->max < 1e-9 * $L->abs->max, 'at an electrode: that row');

# 重心では 4 本の平均(1 次関数の重心での値)
my $c = $ny->interp_rows($pos->average->dummy(0), leadfield => $L);
ok(($c - $L->average->dummy(0))->abs->max < 1e-9 * $L->abs->max, 'at the centroid: the mean of the 4 rows');

# 名前で抜く(4 本から 1 本抜くと 3 本で足りない)
my $err = eval { $ny->interp_rows($pos->slice('0:0,:'), exclude => ['Cz']); 1 } ? '' : $@;
like($err, qr/at least 4 electrodes/, 'exclude by name reaches the electrode count check');
$err = eval { $ny->interp_rows($pos->slice('0:0,:'), exclude => ['NoSuch']); 1 } ? '' : $@;
like($err, qr/no electrode named NoSuch/, 'unknown name in exclude croaks');

done_testing;
