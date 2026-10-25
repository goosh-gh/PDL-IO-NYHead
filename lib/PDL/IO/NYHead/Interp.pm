package PDL::IO::NYHead::Interp;
# ============================================================================
# PDL::IO::NYHead::Interp — 電極の位置から、任意の点の導出行列の行を内挿する。
#
#   電位は場所についてなめらかなので、点 p の行を電極の行の重み付きの和
#   L_p = Σ_i w_i(p) L_i で作る。重みは電極と点の位置だけで決まり、全電流源に
#   同じ重みを使えるので、行の内挿は行列の積 1 回(L x Wᵀ)で済む。
#
#   重みは 3 次元の多重調和スプライン(φ(r) = r または r³)に 1 次の多項式
#   (1, x, y, z)を足した放射基底関数の補間から作る。1 次の項があるので、
#   定数と 1 次関数は正確に再現され(重みの和 = 1、Σ w_i e_i = p)、
#   平均基準の行列から内挿した行も同じ平均基準のままになる。
#
# ---- 次元の約束(PDL の (列, 行)) --------------------------------------
#   電極の座標 $E (ne,3)、点 $P (np,3)、重み $W (ne,np)、
#   導出行列 $L (ne,ns)(PDL::IO::NYHead の leadfield と同じ)、
#   内挿した行 (np,ns)。
#
# 由来: PDL::IO::NYHead05(HBrBP 008 の Cv12 の行のため)。
# ============================================================================
use strict;
use warnings;
use Carp;
use PDL;
use Exporter 'import';

our @EXPORT_OK = qw(interp_weights interp_rows row_errors loo distances);
our $VERSION   = '0.04';

my %KERNEL = (
    r  => sub { $_[0] },
    r3 => sub { $_[0]**3 },
);

# distances($A (na,3), $B (nb,3)) -> (na,nb) mm: 座標の軸ごとに書き下す(sumover の軸を取り違えない)
sub distances {
    my ($A, $B) = @_;
    my $d2 = zeroes(double, $A->dim(0), $B->dim(0));
    for my $c (0 .. 2) {
        $d2 += ($A->slice(":,($c)")->dummy(1) - $B->slice(":,($c)")->dummy(0))**2;
    }
    return sqrt($d2);
}

# interp_weights($E (ne,3), $P (np,3), %opt) -> $W (ne,np)
#   kernel  => 'r' | 'r3'            φ(r)                                   ('r')
#   k       => N                     点ごとに近い N 本の電極だけを使う        (全部)
#   exclude => [i, ...]              使わない電極(0 始まり。1 個抜き検証用)
#   使わない電極の重みは 0。
sub interp_weights {
    my ($E, $P, %opt) = @_;
    my $phi = $KERNEL{ $opt{kernel} // 'r' } or croak "interp_weights: kernel must be one of: ", join(', ', sort keys %KERNEL);
    $E = $E->double; $P = $P->double;
    $P = $P->dummy(0) if $P->ndims == 1;
    $E->ndims == 2 && $E->dim(1) == 3 or croak "interp_weights: electrodes must be (ne,3), got (", join(',', $E->dims), ")";
    $P->ndims == 2 && $P->dim(1) == 3 or croak "interp_weights: points must be (np,3), got (", join(',', $P->dims), ")";
    my ($ne, $np) = ($E->dim(0), $P->dim(0));
    my %ex = map { $_ => 1 } @{ $opt{exclude} // [] };
    my @use = grep { !$ex{$_} } 0 .. $ne - 1;
    @use >= 4 or croak "interp_weights: needs at least 4 electrodes, have ", scalar @use;

    my $W = zeroes(double, $ne, $np);
    if (!defined $opt{k} || $opt{k} >= @use) {
        my $w = _solve($E->dice_axis(0, pdl(long, \@use)), $P, $phi);   # (nu,np)
        $W->dice_axis(0, pdl(long, \@use)) .= $w;
        return $W;
    }
    $opt{k} >= 5 or croak "interp_weights: k must be at least 5";
    my $Eu = $E->dice_axis(0, pdl(long, \@use));
    my $D  = distances($Eu, $P);                                        # (nu,np)
    for my $p (0 .. $np - 1) {
        my $near = $D->slice(":,($p)")->qsorti->slice('0:' . ($opt{k} - 1));
        my $sel  = pdl(long, \@use)->index($near);
        my $w    = _solve($E->dice_axis(0, $sel), $P->slice("$p:$p,:"), $phi);   # (k,1)
        $W->slice(":,($p)")->index($sel) .= $w->slice(':,(0)');
    }
    return $W;
}

# 放射基底関数 + 1 次の項の連立方程式の行列 A (n+4,n+4)。座標は使う電極の重心で中心化し、
# 大きさ(重心からの距離の RMS)で割ってから組む(条件数のため。補間の結果は変わらない)。
# 返り値: ($A, $center (3), $scale)
sub _system {
    my ($Eu, $phi) = @_;
    my $n   = $Eu->dim(0);
    my $c   = $Eu->average;                                         # (3): 電極の軸(dim 0)を畳む
    my $s   = sqrt((($Eu - $c->dummy(0))**2)->mv(1, 0)->sumover->avg);
    $s > 0 or croak "interp_weights: the electrodes coincide";
    my $Es  = ($Eu - $c->dummy(0)) / $s;                            # (n,3)
    my $A   = zeroes(double, $n + 4, $n + 4);
    $A->slice("0:" . ($n - 1) . ",0:" . ($n - 1)) .= $phi->(distances($Es, $Es));
    my $Q   = ones(double, $n, 1)->glue(1, $Es);                    # (n,4): 1, x, y, z
    $A->slice("0:" . ($n - 1) . ",$n:" . ($n + 3)) .= $Q;
    $A->slice("$n:" . ($n + 3) . ",0:" . ($n - 1)) .= $Q->transpose;
    return ($A, $c, $s);
}

# _lu_solve($A, $B (n+4,nrhs)) -> X (n+4,nrhs): 右辺ごとに A x = b を解き、残差を確かめる
sub _lu_solve {
    my ($A, $B) = @_;
    my ($lu, $perm, $par) = lu_decomp($A);
    my $X = lu_backsub($lu, $perm, $par, $B);
    my $res = (($A x $X->transpose) - $B->transpose)->abs->max;
    $res < 1e-6 * (1 + $B->abs->max) * (1 + $X->abs->max) or croak sprintf "interp_weights: the system is ill-conditioned (residual %.2g)", $res;
    return $X;
}

sub _solve {
    my ($Eu, $P, $phi) = @_;
    my $n  = $Eu->dim(0);
    my ($A, $c, $s) = _system($Eu, $phi);
    my $Es = ($Eu - $c->dummy(0)) / $s;
    my $Ps = ($P  - $c->dummy(0)) / $s;                             # (np,3)
    my $np = $Ps->dim(0);
    my $B  = $phi->(distances($Es, $Ps))->glue(0, ones(double, 1, $np)->glue(0, $Ps->transpose));   # (n+4,np)
    return _lu_solve($A, $B)->slice('0:' . ($n - 1) . ',:')->sever;  # (n,np)
}

# interp_rows($L (ne,ns), $W (ne,np)) -> (np,ns): 行の重み付きの和を行列の積 1 回で
sub interp_rows {
    my ($L, $W) = @_;
    $L->dim(0) == $W->dim(0) or croak "interp_rows: leadfield has ", $L->dim(0), " electrodes, weights ", $W->dim(0);
    return $L x $W->transpose;
}

# row_errors($est (np,ns), $true (np,ns)) -> 各点の
#   corr(電流源にわたる相関)と rel(差の RMS / 本物の行の RMS)
sub row_errors {
    my ($est, $true) = @_;
    $est  = $est->dummy(0)  if $est->ndims == 1;
    $true = $true->dummy(0) if $true->ndims == 1;
    my $te = $est->mv(1, 0);  my $tt = $true->mv(1, 0);            # (ns,np): 電流源の軸を先頭に
    my $ce = $te - $te->average->dummy(0);
    my $ct = $tt - $tt->average->dummy(0);
    my $corr = ($ce * $ct)->sumover / sqrt(($ce**2)->sumover * ($ct**2)->sumover);
    my $rel  = sqrt((($te - $tt)**2)->average) / sqrt(($tt**2)->average);
    return { corr => $corr, rel => $rel };
}

# loo($E (ne,3), $L (ne,ns), idx => [i,...], %opt) -> 1 個抜き検証
#   各 i を抜いて残りから i の位置の行を内挿し、本物の行と比べる。
#   gap => G を付けると、i から G mm 未満の電極もいっしょに抜く(最も近い電極まで G mm 以上の
#   点を内挿するときの誤差を、どの電極でも作れる)。
#   返り値: { corr (ni), rel (ni): 差の RMS / 本物の行の RMS, abs (ni): 差の RMS,
#            nn (ni): 残りの電極で最も近いものまでの距離 mm }
#   k を付けない(電極を全部使う)ときは、連立方程式の逆行列 1 回から閉じた式で出す
#   (Rippa 1999 を集合に広げたもの: 集合 S を抜いたときの S での誤差 = (B_SS)^-1 c_S、
#   B = A^-1 の電極の部分、c = B f)。explicit => 1 で、抜くたびに解き直す(テストの照合用。遅い)。
sub loo {
    my ($E, $L, %opt) = @_;
    my @idx = @{ delete $opt{idx} // [ 0 .. $E->dim(0) - 1 ] };
    my @ex  = @{ delete $opt{exclude} // [] };
    my $explicit = delete $opt{explicit};
    my $gap = delete $opt{gap};
    my $E2  = $E->double;
    my $ne  = $E2->dim(0);
    my %exh = map { $_ => 1 } @ex;
    !grep { $exh{$_} } @idx or croak "loo: an electrode in idx is also in exclude";
    my $D   = distances($E2, $E2->dice_axis(0, pdl(long, \@idx)));  # (ne,ni)
    my @set;                                                        # 抜く集合(i が先頭、ex は含めない)
    for my $j (0 .. $#idx) {
        my $i = $idx[$j];
        my @near = defined $gap ? grep { $_ != $i && !$exh{$_} } which($D->slice(":,($j)") < $gap)->list : ();
        push @set, [ $i, @near ];
        $D->slice("($_),($j)") .= 9e9 for @{ $set[-1] }, @ex;
    }
    my $truth = $L->dice_axis(0, pdl(long, \@idx))->double;          # (ni,ns)
    my @use = grep { !$exh{$_} } 0 .. $ne - 1;
    my $est;
    if (!$explicit && (!defined $opt{k} || $opt{k} >= @use - 1)) {
        my $phi = $KERNEL{ $opt{kernel} // 'r' } or croak "loo: kernel must be one of: ", join(', ', sort keys %KERNEL);
        my $n   = @use;
        my ($A) = _system($E2->dice_axis(0, pdl(long, \@use)), $phi);
        my $B   = _lu_solve($A, identity($n + 4))->slice('0:' . ($n - 1) . ',0:' . ($n - 1));   # (n,n), symmetric
        my %pos; @pos{@use} = 0 .. $n - 1;
        my $C   = $L->dice_axis(0, pdl(long, \@use))->double x $B;   # (n,ns): c = B f for every source
        $est = $truth->copy;
        for my $j (0 .. $#idx) {
            my $S   = pdl(long, [ map { $pos{$_} } @{ $set[$j] } ]);
            my $Bss = $B->dice_axis(0, $S)->dice_axis(1, $S);
            my $v   = inv($Bss)->slice('(0),')->dummy(0);          # (1,|S|): column of (B_SS)^-1 for i (symmetric)
            $est->slice("$j:$j,:") -= $C->dice_axis(0, $S) x $v;   # (|S|,ns) x (1,|S|) -> (1,ns): e_i over the sources
        }
    } else {
        my $W = zeroes(double, $ne, scalar @idx);                  # 抜いた電極ごとの重みを列に集める
        for my $j (0 .. $#idx) {
            my $i = $idx[$j];
            $W->slice(":,($j)") .= interp_weights($E2, $E2->slice("$i:$i,:"), %opt, exclude => [ @{ $set[$j] }, @ex ])->slice(':,(0)');
        }
        $est = interp_rows($L, $W);                                # 積は 1 回
    }
    my $err = row_errors($est, $truth);
    my $abs = sqrt((($est - $truth)**2)->mv(1, 0)->average);       # (ni): over the sources (dim 1)
    return { corr => $err->{corr}, rel => $err->{rel}, abs => $abs, nn => $D->minimum };
}

1;

__END__

=head1 NAME

PDL::IO::NYHead::Interp - interpolate leadfield rows at arbitrary points from the electrode rows

=head1 SYNOPSIS

  use PDL::IO::NYHead::Interp qw(interp_weights interp_rows row_errors loo);
  my $W    = interp_weights($E, $P, kernel => 'r');   # $E (ne,3), $P (np,3) -> (ne,np)
  my $rows = interp_rows($L, $W);                     # $L (ne,ns) -> (np,ns)
  my $err  = row_errors($rows, $true);                # {corr, rel}, one per point
  my $cv   = loo($E, $L, idx => [5, 9]);              # leave-one-out: {corr, rel, nn}

=head1 DESCRIPTION

Polyharmonic radial basis function interpolation in 3-D (phi(r) = r or r^3) with a linear
polynomial term, so constants and linear functions are reproduced exactly (the weights sum
to 1, and an average-referenced leadfield stays average-referenced). The weights depend
only on the positions, so all sources are interpolated with one matrix product.
Option C<k> uses only the k nearest electrodes of each point; C<exclude> leaves electrodes
out (leave-one-out).

=cut
