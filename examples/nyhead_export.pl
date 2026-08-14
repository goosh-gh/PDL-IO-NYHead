#!/usr/bin/env perl
# ============================================================================
# nyhead_export.pl — New York Head のサーフェスを .obj / .usda に書き出す。
#
# OBJ と USDA は同じ幾何(頂点+面)を吐くだけで、違いは面 index の base
# (OBJ=1-based そのまま / USD=0-based で tri-1)とシリアライズ形式のみ。
# よって1スクリプト + --format 分岐で扱う(PDL::IO::NYHead::write_mesh)。
#
#   perl examples/nyhead_export.pl --format obj  --surf head       sa_nyhead.mat
#   perl examples/nyhead_export.pl --format usda --surf cortex75K  sa_nyhead.mat -o cortex.usda
#   perl examples/nyhead_export.pl --format obj  --surf cortex2K   sa_nyhead.mat
#
# --surf は head / cortex75K / 低解像 cortexNK(cortex1K,2K,5K,10K 等、実在するもの)。
# -o 省略時は "<basename>_<surf>.<ext>" を CWD に作る。
# 注意: これはサーフェス読取の確認であって MRI ボリューム(/sa/mri)の確認ではない。
# ============================================================================
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Getopt::Long;
use File::Basename qw(basename);
use PDL::IO::NYHead;

my ($format, $surf, $out, $h5dump, $axes) = ('obj', 'head', undef, 'h5dump', 0);
GetOptions(
    'format=s' => \$format,
    'surf=s'   => \$surf,
    'o|out=s'  => \$out,
    'h5dump=s' => \$h5dump,
    'axes=f'   => \$axes,     # USDA のみ: 原点から長さ <mm> の xyz 軸を併記(0=なし, 既定)
) or die "bad args\n";
my $file = shift or die <<"USAGE";
usage: $0 [--format obj|usda] [--surf NAME] [-o OUT] <sa_nyhead.mat>
  --surf : head | cortex75K | cortexNK (1K/2K/5K/10K, whichever exist)
USAGE

$format = lc $format;
$format =~ /^(obj|usda)$/ or die "--format must be obj or usda\n";

my $ny = PDL::IO::NYHead->new($file, h5dump => $h5dump);
$out //= sprintf "%s_%s.%s", (basename($file) =~ s/\.[^.]+$//r), $surf, $format;

my $m = $ny->mesh($surf);
printf "surf=%s  vertices=%d  faces=%d  ->  %s (%s)\n",
    $surf, $m->{vc}->dim(0), $m->{tri}->dim(0), $out, $format;

$ny->write_mesh($surf, $format, $out, axes => $axes);
print "done.\n";
