use strict;
use warnings;
use Test::More;
use PDL;
use PDL::IO::NYHead::Interp qw(interp_weights interp_rows row_errors loo);

# 合成のデータ(実装とは独立の式):
#   電極: 半径 90 mm の球の z > -30 の部分に黄金螺旋で 80 点(間隔 約 30 mm)。
#   行: 無限一様媒質の電流双極子の電位 V(p) = m·(p - q) / |p - q|^3 を電流源ごとに。
#       電流源は半径 55 mm の球の中。並びは決まった式で作る(乱数を使わない)。
my $R  = 90;
my $NE = 80;
my @e;
for my $i (0 .. $NE - 1) {
    my $z   = $R * (1 - ($i + 0.5) / $NE * (1 + 30 / $R));     # R .. -30
    my $rxy = sqrt($R**2 - $z**2);
    my $th  = $i * 3.14159265358979 * (3 - sqrt(5));
    push @e, [ $rxy * cos($th), $rxy * sin($th), $z ];
}
my $E = pdl(\@e)->transpose->copy;                             # (80,3)
is(join(',', $E->dims), '80,3', 'electrodes (ne,3)');

sub dipole_rows {                                              # ($P (np,3)) -> (np,ns)
    my ($P, $ns) = @_;
    my @rows;
    for my $s (0 .. $ns - 1) {
        my $q = pdl(55 * 0.9 * sin(1.3 * $s) * cos(0.7 * $s), 55 * 0.9 * sin(1.3 * $s) * sin(0.7 * $s), 50 * cos(1.3 * $s) - 10);
        my $m = pdl(cos(2.1 * $s), sin(2.1 * $s) * cos(0.4 * $s), sin(0.4 * $s) + 0.3);
        my $dx = $P->slice(':,(0)') - $q->at(0);
        my $dy = $P->slice(':,(1)') - $q->at(1);
        my $dz = $P->slice(':,(2)') - $q->at(2);
        my $r3 = ($dx**2 + $dy**2 + $dz**2)**1.5;
        push @rows, ($m->at(0) * $dx + $m->at(1) * $dy + $m->at(2) * $dz) / $r3;
    }
    return cat(@rows);                                          # (np,ns)
}
sub linear_rows {                                              # 電流源ごとに a + b·p
    my ($P, $ns) = @_;
    cat(map { 0.3 * $_ - 1 + sin($_) * $P->slice(':,(0)') + cos(0.5 * $_) * $P->slice(':,(1)') - 0.2 * $_ * $P->slice(':,(2)') } 0 .. $ns - 1);
}
sub hash_rows {                                                # 位置と無関係の値(否定対照)
    my ($np, $ns) = @_;
    my $i = xvals(double, $np, $ns); my $j = yvals(double, $np, $ns);
    my $h = sin(12.9898 * ($i + 1) + 78.233 * ($j + 1)) * 43758.5453;
    return $h - floor($h) - 0.5;
}

my $NS = 40;
my $P  = pdl([[0, -20, 85], [30, 40, 70], [-60, -50, 20]])->transpose->copy;   # (3,3): 電極の間の点
$P = $P * ($R / sqrt(($P**2)->mv(1, 0)->sumover))->dummy(1);                     # 球面の上へ

# ---- 1 次関数は正確に再現する、重みの性質 -------------------------------------
for my $kern (qw(r r3)) {
    for my $k (undef, 12) {
        my $tag = "kernel $kern, k " . ($k // 'all');
        my $W = interp_weights($E, $P, kernel => $kern, k => $k);
        is(join(',', $W->dims), '80,3', "$tag: weights (ne,np)");
        ok(all(abs($W->sumover - 1) < 1e-10), "$tag: weights sum to 1");
        ok((interp_rows($E, $W) - $P)->abs->max < 1e-9, "$tag: positions reproduced (sum w_i e_i = p)");
        my $est = interp_rows(linear_rows($E, $NS), $W);
        is(join(',', $est->dims), "3,$NS", "$tag: rows (np,ns)");
        my $err = row_errors($est, linear_rows($P, $NS));
        ok($err->{rel}->max < 1e-9, "$tag: linear rows reproduced (rel " . sprintf('%.1e', $err->{rel}->max) . ')');
        ok(($W->slice(':,(0)') != 0)->sum == ($k // 80), "$tag: number of electrodes used");
    }
}

# ---- 電極そのものの位置では、その電極の行 -------------------------------------
{
    my $Ld = dipole_rows($E, $NS);               # (ne,ns)
    my $W  = interp_weights($E, $E->slice('17:17,:'));
    ok(abs($W->at(17, 0) - 1) < 1e-9 && ($W->abs->sum - 1) < 1e-8, 'at an electrode the weight is one-hot');
    ok((interp_rows($Ld, $W) - $Ld->slice('17:17,:'))->abs->max < 1e-9 * $Ld->abs->max, 'and the row is that electrode');

    # 平均基準: 平均基準の行から内挿しても、内挿してから引いても同じ
    my $mean = $Ld->average;                                       # (ns): 電極の軸(dim 0)の平均
    my $Wp   = interp_weights($E, $P);
    my $d    = (interp_rows($Ld - $mean->dummy(0), $Wp) - (interp_rows($Ld, $Wp) - $mean->dummy(0)))->abs->max;
    ok($d < 1e-9 * $Ld->abs->max, 'interpolation commutes with the average reference');
}

# ---- 1 個抜き検証: なめらかな場は戻り、位置と無関係の値は戻らない -------------------
my @idx = map { 5 * $_ + 2 } 0 .. 13;                              # 球の上から下まで 14 本
{
    my $Ld = dipole_rows($E, $NS);
    for my $kern (qw(r r3)) {
        my $cv = loo($E, $Ld, idx => \@idx, kernel => $kern);
        is($cv->{corr}->nelem, 14, "loo kernel $kern: one result per electrode");
        note(sprintf "dipole field, kernel %s: corr median %.3f (min %.3f), rel median %.3f (max %.3f), nearest %.0f..%.0f mm",
            $kern, $cv->{corr}->median, $cv->{corr}->min, $cv->{rel}->median, $cv->{rel}->max, $cv->{nn}->min, $cv->{nn}->max);
        ok($cv->{corr}->median > 0.97, "loo kernel $kern: smooth (dipole) rows come back (corr median > 0.97)");
        ok($cv->{rel}->median < 0.25,  "loo kernel $kern: rel error median < 0.25");
    }
    my $Lh = hash_rows($NE, $NS);
    my $cv = loo($E, $Lh, idx => \@idx);
    note(sprintf "position-free values: corr median %.3f, rel median %.3f", $cv->{corr}->median, $cv->{rel}->median);
    ok(abs($cv->{corr}->median) < 0.4, 'negative control: position-free rows do not come back (|corr| median < 0.4)');
    ok($cv->{rel}->median > 0.8,       'negative control: rel error median > 0.8');

    # 遠いほど悪くなる: 間隔の 2 倍の球(20 点)で抜くと、80 点より誤差が大きい
    my @sub = map { 4 * $_ } 0 .. 19;
    my $E20 = $E->dice_axis(0, pdl(long, \@sub));
    my $cv20 = loo($E20, $Ld->dice_axis(0, pdl(long, \@sub)), idx => [ 0 .. 19 ]);
    my $cv80 = loo($E, $Ld, idx => \@sub);
    note(sprintf "rel median: 80 electrodes %.3f, 20 electrodes %.3f", $cv80->{rel}->median, $cv20->{rel}->median);
    ok($cv20->{rel}->median > $cv80->{rel}->median, 'sparser electrodes give larger errors');
}

# ---- 1 個抜きの閉じた式(Rippa)は、1 本ずつ解き直したものと一致する -------------
{
    my $Ld = dipole_rows($E, $NS);
    for my $kern (qw(r r3)) {
        for my $ex ([], [3, 40]) {
            my $fast = loo($E, $Ld, idx => \@idx, kernel => $kern, exclude => $ex);
            my $slow = loo($E, $Ld, idx => \@idx, kernel => $kern, exclude => $ex, explicit => 1);
            my $d  = ($fast->{rel} - $slow->{rel})->abs->max;
            my $dc = ($fast->{corr} - $slow->{corr})->abs->max;
            ok($d < 1e-8 && $dc < 1e-8 && all($fast->{nn} == $slow->{nn}),
               sprintf "loo kernel %s, exclude [%s]: closed form = explicit (rel %.1e, corr %.1e)", $kern, "@$ex", $d, $dc);
        }
    }
    for my $kern (qw(r r3)) {                                      # 近くの電極もいっしょに抜く
        my $fast = loo($E, $Ld, idx => \@idx, kernel => $kern, gap => 45, exclude => [40]);
        my $slow = loo($E, $Ld, idx => \@idx, kernel => $kern, gap => 45, exclude => [40], explicit => 1);
        my $d = ($fast->{rel} - $slow->{rel})->abs->max;
        ok($d < 1e-8 && ($fast->{abs} - $slow->{abs})->abs->max < 1e-8 * $Ld->abs->max && all($fast->{nn} == $slow->{nn}),
           sprintf "loo kernel %s, gap 45 mm: closed form = explicit (rel %.1e)", $kern, $d);
        ok(all($fast->{nn} >= 45), sprintf 'gap 45 mm: nearest remaining electrode %.1f .. %.1f mm', $fast->{nn}->min, $fast->{nn}->max);
        my $plain = loo($E, $Ld, idx => \@idx, kernel => $kern, exclude => [40]);
        ok($fast->{rel}->median > $plain->{rel}->median, sprintf 'gap 45 mm gives larger errors (rel median %.3f > %.3f)', $fast->{rel}->median, $plain->{rel}->median);
    }
    eval { loo($E, $Ld, idx => [3], exclude => [3]) };
    like($@, qr/also in exclude/, 'loo croaks when idx and exclude overlap');
}

# ---- 入力の誤り --------------------------------------------------------------------
eval { interp_weights($E->transpose, $P) };
like($@, qr/electrodes must be \(ne,3\)/, 'croaks on (3,ne) electrodes');
eval { interp_weights($E, $P, kernel => 'gauss') };
like($@, qr/kernel must be/, 'croaks on an unknown kernel');

done_testing;
