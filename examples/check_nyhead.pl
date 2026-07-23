#!/usr/bin/env perl
# ============================================================================
# check_nyhead.pl — PDL::IO::NYHead の動作確認 & API デモ。
#
# 実 sa_nyhead.mat に対して reader の全経路を叩き、要約を表示しつつ
# 構造整合(サイズ・index 範囲・軸)を PASS/FAIL で検査する。実データ用の
# スモークテスト。231ch/74382頂点に依存する検査は自動でスキップ判定するので
# 縮小フィクスチャでもそのまま走る。
#
#   perl examples/check_nyhead.pl sa_nyhead.mat
#   perl examples/check_nyhead.pl --h5dump /opt/local/bin/h5dump sa_nyhead.mat
# ============================================================================
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Getopt::Long;
use PDL;
use PDL::IO::NYHead;

my $h5dump = 'h5dump';
GetOptions('h5dump=s' => \$h5dump) or die "bad args\n";
my $file = shift @ARGV or die "usage: $0 [--h5dump PATH] <sa_nyhead.mat>\n";

my $ny = PDL::IO::NYHead->new($file, h5dump => $h5dump);

my ($pass, $fail) = (0, 0);
sub check { my ($cond, $msg) = @_;
    printf "  [%s] %s\n", $cond ? 'PASS' : 'FAIL', $msg;
    $cond ? $pass++ : $fail++;
}
sub section { print "\n== $_[0] ==\n" }

# ---- 電極 -----------------------------------------------------------------
section("Electrodes");
my $labels = $ny->electrode_labels;
my $pos    = $ny->electrode_pos;
my $ne     = $ny->n_electrodes;
printf "  n_electrodes = %d   labels[0..2] = %s ...\n", $ne, join(", ", @{$labels}[0..($ne>3?2:$ne-1)]);
printf "  electrode_pos dims = (%s)\n", join(",", $pos->dims);
check(scalar(@$labels) == $ne,          "label count == n_electrodes");
check($pos->dim(1) == 3,                "electrode_pos is (N,3)");
check($pos->dim(0) == $ne,              "electrode_pos rows == n_electrodes");

# ---- Leadfield ------------------------------------------------------------
section("Leadfield");
my $L  = $ny->leadfield;                       # (Ne, Nv)
my $nv = $ny->surface('cortex75K')->{vc}->dim(0);
printf "  leadfield dims = (%s)   [expect (%d, cortex75K vertices)]\n", join(",", $L->dims), $ne;
check($L->dim(0) == $ne,                "leadfield dim0 == n_electrodes");
check($L->dim(1) == $nv,                "leadfield dim1 == cortex75K vertex count");

# ---- Harvard-Oxford atlas -------------------------------------------------
section("Atlas (Harvard-Oxford)");
my $ho = $ny->ho_labels;
my $ai = $ny->atlas_index;
printf "  ho_labels = %d areas   e.g. [0]=%s  [last]=%s\n", scalar(@$ho), $ho->[0], $ho->[-1];
printf "  atlas_index range = %d .. %d over %d vertices\n", $ai->min, $ai->max, $ai->nelem;
check($ai->max <= scalar(@$ho),         "max atlas index <= number of HO areas");
check($ai->min >= 0,                    "atlas index >= 0 (0 = none)");

# ---- Surfaces -------------------------------------------------------------
section("Surfaces");
for my $s (qw(head cortex75K)) {
    my $surf = $ny->surface($s);
    printf "  %-9s vc=(%s)  tri=(%s)  tri range %d..%d (1-based)\n",
        $s, join(",",$surf->{vc}->dims), join(",",$surf->{tri}->dims),
        $surf->{tri}->min, $surf->{tri}->max;
    check($surf->{vc}->dim(1) == 3,                    "$s vc is (Nv,3)");
    check($surf->{tri}->min >= 1,                      "$s tri is 1-based");
    check($surf->{tri}->max <= $surf->{vc}->dim(0),    "$s tri indices within vertex range");
    check($surf->{tri0}->min == $surf->{tri}->min - 1, "$s tri0 == tri-1");
}

# ---- 低解像 cortex の fancy-index 復元 ------------------------------------
section("Low-res cortex reconstruction");
for my $res ($ny->cortex_resolutions) {
    next if $res eq 'cortex75K';
    my $c = $ny->cortex($res);
    printf "  %-9s in_from=%d  vc=(%s)\n", $res, $c->{in_from}->nelem, join(",",$c->{vc}->dims);
    check($c->{vc}->dim(0) == $c->{in_from}->nelem, "$res vc rows == in_from count");
    check($c->{in_from}->max < $nv,                 "$res in_from indices within 75K range (0-based)");
}

# ---- nearest_vertex(19ch を頭皮に投影)-- 実データ(231ch)のみ ------------
section("Nearest scalp vertex (19ch) — requires full 231ch model");
if ($ne >= 224) {
    my $p19  = $ny->electrode_pos_19;                 # (19,3)
    my $iv   = $ny->nearest_vertex($p19, 'head');     # (19)
    my $hvc  = $ny->surface('head')->{vc};
    my $near = $hvc->dice_axis(0, $iv);               # (19,3)
    my $resid = (($p19 - $near)**2)->xchg(0,1)->sumover->sqrt;  # (19) mm
    printf "  19ch -> nearest scalp vertex, residual mean %.2f mm, max %.2f mm\n",
        $resid->avg, $resid->max;
    my $lab19 = $ny->electrode_labels_19;
    printf "  labels_19[0..4] = %s\n", join(", ", @{$lab19}[0..4]);
    check($iv->nelem == 19,        "nearest_vertex returned 19 indices");
    check($resid->max < 20,        "all 19ch within 20mm of a scalp vertex");
} else {
    print "  (skipped: fixture has $ne electrodes, IDX19 needs >=224)\n";
}

# ---- 集計 -----------------------------------------------------------------
printf "\n%s  (%d pass, %d fail)\n", ($fail ? "*** FAILURES ***" : "ALL CHECKS PASSED"), $pass, $fail;
exit($fail ? 1 : 0);
