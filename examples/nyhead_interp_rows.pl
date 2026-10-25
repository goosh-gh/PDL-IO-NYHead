#!/usr/bin/env perl
#
# nyhead_interp_rows.pl — leadfield rows of the New York Head at own points (neck electrodes
# such as Cv12), interpolated from the rows of the 231 NYHead electrodes, with the error of
# the interpolation estimated by leave-one-out, and a table of row differences between any
# pairs of NYHead electrodes and own points (how much brain activity a bipolar lead carries).
#
#   perl examples/nyhead_interp_rows.pl --points neck_points.tsv --pairs 'Cv12-Cv45,Cv12-Nk2,Cz-Cv12,Exxz-Nk2'
#   perl examples/nyhead_interp_rows.pl --point 'Cv12=0,-113.0,-114.5' --point 'Cv45=0,-103.5,-152.9'
#   perl examples/nyhead_interp_rows.pl --selftest          # synthetic sphere, no .mat needed
#
#   --mat FILE        sa_nyhead.mat                                         (~/src/NYHead/sa_nyhead.mat)
#   --point N=X,Y,Z   an own point in MNI mm (repeatable)
#   --points TSV      own points from a file: columns name, x, y, z[, kind]; with a kind column
#                     only the rows of kind 'point' are used (the --points-out file of
#                     nyhead_neck_electrodes.pl)
#   --method K:N      how the rows of the own points are interpolated: kernel r or r3 (phi(r)),
#                     N = all or the number of nearest electrodes                      (r3:all)
#   --methods LIST    methods checked by leave-one-out                   (r:all,r3:all,r3:12)
#   --show LIST       electrodes listed one by one in the leave-one-out table
#                                                                 (Exxz,Exx1,Exx2,Nk2,Iz,Oz,Cz)
#   --pairs LIST      A-B,...  pairs of NYHead electrodes and own points for the difference table
#   --strings S       auto | h5dump | lowlevel (how the names are read)                  (auto)
#
# Printed:
#  1. the own points: the nearest NYHead electrodes, a flag (inside: within 45 mm of an electrode
#     and not below the lowest electrode; otherwise 'far' or 'below the lowest electrode',
#     where the interpolation becomes an extrapolation), and the electrodes with the largest
#     weights for --method;
#  2. leave-one-out: every electrode is left out in turn and its row is interpolated from the
#     others; correlation over the 74382 sources and relative RMS error (|est - true| RMS /
#     true RMS), summarised by the distance from the left-out electrode to its nearest
#     remaining one, and listed for the --show electrodes.  Then, for each own point, the
#     same with a gap: every electrode is left out together with the electrodes closer to it
#     than the point is to its nearest electrode, so that every electrode is interpolated
#     from the same spacing as the point; the RMS error of the row relative to the median
#     row RMS (median and 90th percentile) is the expected error of the point's row, over
#     all electrodes and over those lower than Iz and behind y = -60 mm (the lower back of
#     the head, around the neck points).  For a point marked extrapolated this is a lower
#     bound;
#  3. with --pairs: distance, correlation of the rows, and the RMS of the row difference
#     relative to the median row RMS of the 231 electrodes, with the expected error of the
#     difference (the errors of the interpolated rows of section 2, median and 90th
#     percentile, combined as independent when both ends are own points) and the range of
#     the difference over the --methods; for comparison each non-neck electrode minus its
#     nearest non-neck neighbour (median and range).
# The rows are V_fem_normal as stored (average reference over the 231 electrodes); the
# interpolation keeps that reference.  The model's neck ends below Nk1 (z -179): its effect
# remains in interpolated rows too.  Row signs are arbitrary: compare correlations by size.

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Getopt::Long;
use PDL;
use PDL::IO::NYHead;
use PDL::IO::NYHead::Interp qw(interp_weights interp_rows row_errors loo distances);

my %o = (mat => "$ENV{HOME}/src/NYHead/sa_nyhead.mat", point => [], points => undef,
         method => 'r3:all', methods => 'r:all,r3:all,r3:12', show => 'Exxz,Exx1,Exx2,Nk2,Iz,Oz,Cz',
         pairs => '', strings => 'auto', selftest => 0);
GetOptions('mat=s' => \$o{mat}, 'point=s@' => $o{point}, 'points=s' => \$o{points},
           'method=s' => \$o{method}, 'methods=s' => \$o{methods}, 'show=s' => \$o{show},
           'pairs=s' => \$o{pairs}, 'strings=s' => \$o{strings}, 'selftest' => \$o{selftest})
    or die "bad options (see the header of $0)\n";
if ($o{selftest}) { exit(selftest() ? 0 : 1) }
-f $o{mat} or die "not found: $o{mat} (give it with --mat)\n";

my $ny = PDL::IO::NYHead->new($o{mat}, strings => $o{strings});
printf "%s: strings read by %s\n", $o{mat}, $ny->strings_method;
my @pts = read_points(\%o);
run(\%o, [ @{ $ny->electrode_labels } ], $ny->electrode_pos->double, $ny->leadfield->double, \@pts);
exit 0;

# =====================================================================================
sub read_points {
    my ($opt) = @_;
    my @p;
    for my $s (@{ $opt->{point} }) {
        my ($n, $xyz) = $s =~ /^([^=]+)=(.+)$/ or die "--point NAME=X,Y,Z, not '$s'\n";
        my @c = split /,/, $xyz; @c == 3 or die "--point $n: three coordinates\n";
        push @p, [ $n, @c ];
    }
    if (defined $opt->{points}) {
        open my $fh, '<', $opt->{points} or die "$opt->{points}: $!\n";
        my @col;
        while (my $l = <$fh>) {
            chomp $l;
            next if $l =~ /^\s*(#|$)/;
            my @f = split /\t/, $l;
            if (!@col) { @col = map { lc } @f; next }
            my %r; @r{@col} = @f;
            defined $r{name} && defined $r{z} or die "$opt->{points}: needs columns name, x, y, z\n";
            next if defined $r{kind} && $r{kind} ne 'point';
            push @p, [ @r{qw(name x y z)} ];
        }
    }
    return @p;
}

sub parse_method {
    my ($m) = @_;
    my ($k, $n) = $m =~ /^(r3?):(all|\d+)$/ or die "a method is KERNEL:N with KERNEL r or r3 and N all or a number, not '$m'\n";
    return (kernel => $k, ($n eq 'all' ? () : (k => $n)));
}

sub run {
    my ($opt, $lab, $E, $L, $pts) = @_;
    my $ne  = $E->dim(0);
    my %at; @at{@$lab} = 0 .. $#$lab;
    my $rms = sqrt(($L**2)->mv(1, 0)->average);                  # (ne): over the sources (dim 1)
    my $med = $rms->median;
    my $zlow = $E->slice(':,(2)')->min;
    my %meth = parse_method($opt->{method});
    printf "%d electrodes, %d sources; median row RMS %.4g; lowest electrode z %.1f mm\n", $ne, $L->dim(1), $med, $zlow;

    # ---- 1. own points ------------------------------------------------------------------
    my ($P, $R, %prow, %ppos, %pnear, %mrow, %perr);
    if (@$pts) {
        $P = pdl(map { [ @$_[1 .. 3] ] } @$pts)->transpose->copy;  # (np,3)
        $P = $P->dummy(0) if $P->ndims == 1;
        my $W = interp_weights($E, $P, %meth);                     # (ne,np)
        $R = interp_rows($L, $W);                                  # (np,ns)
        for my $m (split /,/, $opt->{methods}) {                   # rows by every method, for the range of the differences
            my $Rm = $m eq $opt->{method} ? $R : interp_rows($L, interp_weights($E, $P, parse_method($m)));
            $mrow{$m}{ $pts->[$_][0] } = $Rm->slice("($_),") for 0 .. $#$pts;
        }
        my $D = distances($E, $P);                                 # (ne,np)
        print "\n-- 1. own points: nearest NYHead electrodes, and the largest weights ($opt->{method})\n";
        for my $p (0 .. $#$pts) {
            my $n  = $pts->[$p][0];
            my $dp = $D->slice(":,($p)");
            my @o3 = ($dp->qsorti->list)[0 .. 2];
            my $dn = $dp->at($o3[0]);
            my $z  = $P->at($p, 2);
            my $flag = $z < $zlow ? sprintf('below the lowest electrode by %.0f mm: extrapolated', $zlow - $z)
                     : $dn > 45  ? 'far from the electrodes (> 45 mm): extrapolated'
                     :             'inside';
            printf "  %-6s (%.1f %.1f %.1f)  %s   %s\n", $n, $P->slice("($p),")->list,
                join('  ', map { sprintf '%s %.1f mm', $lab->[$_], $dp->at($_) } @o3), $flag;
            my $wp = $W->slice(":,($p)");
            my @top = ($wp->abs->qsorti->list)[-5 .. -1];
            printf "         weights: %s  (sum %.3f)\n", join('  ', map { sprintf '%s %+.3f', $lab->[$_], $wp->at($_) } reverse @top), $wp->sum;
            $prow{$n} = $R->slice("($p),"); $ppos{$n} = $P->slice("($p),"); $pnear{$n} = $dn;
            $pts->[$p][4] = $flag ne 'inside';
        }
    }

    # ---- 2. leave-one-out -------------------------------------------------------------
    print "\n-- 2. leave-one-out: each electrode's row from the others (corr over sources, rel = RMS error / RMS of the row)\n";
    my @show = grep { exists $at{$_} } split /,/, $opt->{show};
    my %cv;
    for my $m (split /,/, $opt->{methods}) {
        $cv{$m} = loo($E, $L, parse_method($m));
        my ($c, $r, $nn) = @{ $cv{$m} }{qw(corr rel nn)};
        printf "  %-7s all %d: corr median %.3f, rel median %.3f (max %.3f)\n", $m, $ne, $c->median, $r->median, $r->max;
        for my $bin ([0, 20], [20, 30], [30, 45], [45, 1e9]) {
            my $w = which(($nn >= $bin->[0]) * ($nn < $bin->[1]));
            next if $w->isempty;
            printf "          nearest %s mm: %3d electrodes, corr median %.3f (min %.3f), rel median %.3f (max %.3f)\n",
                $bin->[1] > 1e8 ? ">= $bin->[0]" : "$bin->[0]-$bin->[1]", $w->nelem,
                $c->index($w)->median, $c->index($w)->min, $r->index($w)->median, $r->index($w)->max;
        }
    }
    if (@show) {
        my @m = split /,/, $opt->{methods};
        printf "  %-6s %9s%s\n", 'left', 'nearest', join('', map { sprintf '  %15s', "$_ corr/rel" } @m);
        for my $n (@show) {
            my $i = $at{$n};
            printf "  %-6s %6.1f mm%s\n", $n, $cv{ $m[0] }{nn}->at($i),
                join('', map { sprintf '  %7.3f %7.3f', $cv{$_}{corr}->at($i), $cv{$_}{rel}->at($i) } @m);
        }
    }
    if (@$pts) {
        my @low = exists $at{Iz}
            ? grep { $E->at($_, 2) < $E->at($at{Iz}, 2) && $E->at($_, 1) < -60 } 0 .. $ne - 1 : ();
        print "  expected error of the own points' rows ($opt->{method}): each electrode left out with the electrodes\n",
              "  closer to it than the point's nearest electrode (err = RMS error / median row RMS)\n";
        my %bygap;
        for my $p (@$pts) {
            my $g  = sprintf '%.1f', $pnear{ $p->[0] };
            my $cg = $bygap{$g} //= loo($E, $L, %meth, gap => $g);
            my $err = $cg->{abs} / $med;
            my $w  = which($cg->{nn} < 1e8);
            my $line = sub {
                my ($tag, $ix) = @_;
                return sprintf "%s %d: corr median %.3f, rel median %.3f, err median %.3f (90%% %.3f), nearest remaining %.0f..%.0f mm",
                    $tag, $ix->nelem, $cg->{corr}->index($ix)->median, $cg->{rel}->index($ix)->median,
                    $err->index($ix)->median, $err->index($ix)->pct(0.9), $cg->{nn}->index($ix)->min, $cg->{nn}->index($ix)->max;
            };
            printf "  %-6s gap %5.1f mm: %s\n", $p->[0], $g, $line->('all', $w);
            printf "  %-6s %15s %s\n", '', '', $line->('lower back', pdl(long, \@low)) if @low;
            printf "  %-6s %15s %s\n", '', '', 'extrapolated point: these errors are a lower bound' if $p->[4];
            my $ref = @low ? pdl(long, \@low) : $w;
            $perr{ $p->[0] } = [ $err->index($ref)->median, $err->index($ref)->pct(0.9) ];
        }
    }

    # ---- 3. pairs -----------------------------------------------------------------------
    if ($opt->{pairs} ne '') {
        my $row = sub { my $n = shift; exists $prow{$n} ? $prow{$n} : exists $at{$n} ? $L->slice("($at{$n}),") : die "--pairs: no electrode or point named $n\n" };
        my $pos = sub { my $n = shift; exists $ppos{$n} ? $ppos{$n} : $E->slice("($at{$n}),") };
        print "\n-- 3. pairs: distance, correlation of the rows, RMS of the row difference / median row RMS\n";
        for my $pr (split /,/, $opt->{pairs}) {
            my ($na, $nb) = $pr =~ /^([^-]+)-([^-]+)$/ or die "--pairs A-B, not '$pr'\n";
            my ($ra, $rb) = ($row->($na), $row->($nb));
            my $ca = $ra - $ra->avg; my $cb = $rb - $rb->avg;
            my $corr = ($ca * $cb)->sum / sqrt(($ca**2)->sum * ($cb**2)->sum);
            my $diff = sqrt((($ra - $rb)**2)->avg) / $med;
            my $dist = sqrt((($pos->($na) - $pos->($nb))**2)->sum);
            my @ip = grep { exists $prow{$_} } $na, $nb;
            my $extra = '';
            if (@ip) {
                my ($e50, $e90) = (0, 0);
                for (@ip) { $e50 += $perr{$_}[0]**2; $e90 += $perr{$_}[1]**2 }
                my @dm = map {
                    my $m = $_;
                    my ($xa, $xb) = map { exists $mrow{$m}{$_} ? $mrow{$m}{$_} : $L->slice("($at{$_}),") } $na, $nb;
                    (sqrt((($xa - $xb)**2)->avg) / $med)->sclr;
                } sort keys %mrow;
                $extra = sprintf '   +- %.3f (90%% %.3f)   methods %.3f..%.3f', sqrt($e50), sqrt($e90), (sort { $a <=> $b } @dm)[0, -1];
            }
            printf "  %-14s %6.1f mm   corr %6.3f   diff %.3f%s\n", "$na - $nb", $dist->sclr, $corr->sclr, $diff->sclr, $extra;
        }
        my @nn = grep { $lab->[$_] !~ /^Nk\d/ } 0 .. $ne - 1;
        my $En = $E->dice_axis(0, pdl(long, \@nn));
        my $Dn = distances($En, $En);
        $Dn->diagonal(0, 1) .= 9e9;
        my $nbr = $Dn->minimum_ind;                                # (nn): nearest of each column
        my $Ln  = $L->dice_axis(0, pdl(long, \@nn));
        my $dr  = sqrt((($Ln - $Ln->dice_axis(0, $nbr))**2)->mv(1, 0)->average) / $med;   # (nn)
        printf "  for comparison, each non-neck electrode - its nearest neighbour (median %.1f mm apart): diff median %.3f (%.3f .. %.3f)\n",
            $Dn->minimum->median, $dr->median, $dr->min, $dr->max;
    }
    print "\nThe rows are cortical sources only (V_fem_normal, average reference over the 231 electrodes).\n",
          "The model's neck ends below Nk1; interpolated rows near it carry that truncation too.\n";
}

# synthetic check: electrodes on a sphere, rows from current dipoles (V = m.(p-q)/|p-q|^3)
sub selftest {
    my ($R0, $NE, $NS) = (90, 60, 30);
    my (@e, @lab);
    for my $i (0 .. $NE - 1) {
        my $z = $R0 * (1 - ($i + 0.5) / $NE * 1.4); my $rxy = sqrt($R0**2 - $z**2);
        my $th = $i * 3.14159265358979 * (3 - sqrt(5));
        push @e, [ $rxy * cos($th), $rxy * sin($th), $z ]; push @lab, sprintf 'E%02d', $i;
    }
    $lab[10] = 'Nk1';                                              # one 'neck' name, kept out of the neighbour baseline
    my $E = pdl(\@e)->transpose->copy;
    my @rows;
    for my $s (0 .. $NS - 1) {
        my $q = pdl(40 * sin(1.3 * $s), 40 * cos(0.7 * $s), 30 * sin(0.9 * $s));
        my $m = pdl(cos(2.1 * $s), sin(2.1 * $s), 0.5);
        my $d = $E - $q->dummy(0);
        push @rows, ($d * $m->dummy(0))->mv(1, 0)->sumover / (($d**2)->mv(1, 0)->sumover)**1.5;
    }
    my $L = cat(@rows);                                            # (ne,ns)
    $L -= $L->average->dummy(0);                                   # average reference like V_fem_normal
    my $ok = 1;
    my $chk = sub { my ($c, $msg) = @_; printf "  [%s] %s\n", $c ? 'PASS' : 'FAIL', $msg; $ok &&= $c };
    # a point placed exactly at E05 must give E05's row: the pair E05 - P05 has zero difference
    my @pts = ([ 'P05', $E->slice('(5),')->list ], [ 'Pmid', (($E->slice('(5),') + $E->slice('(6),')) / 2)->list ]);
    my $out = '';
    {
        open my $fh, '>', \$out or die; my $old = select $fh;
        run({ method => 'r3:all', methods => 'r3:all', show => 'E05', pairs => 'E05-P05,E05-E06,Pmid-E05' }, \@lab, $E, $L, \@pts);
        select $old; close $fh;
    }
    print $out;
    $chk->($out =~ /E05 - P05\s+0\.0 mm\s+corr\s+1\.000\s+diff 0\.000/, 'a point at an electrode has that electrode\'s row');
    $chk->($out =~ /P05 .*inside/, 'a point at an electrode is flagged inside');
    my ($rel) = $out =~ /r3:all\s+all 60: corr median [\d.]+, rel median ([\d.]+)/;
    $chk->(defined $rel && $rel < 0.3, sprintf 'leave-one-out on a dipole field: rel median %s < 0.3', $rel // 'n/a');
    my ($g0) = $out =~ /Pmid\s+gap\s+16\.1 mm: all 60: .*err median ([\d.]+)/;
    $chk->(defined $g0, 'expected error of a point from the leave-one-out with a gap');
    $chk->($out =~ /Pmid - E05 .*\+- [\d.]+ \(90% [\d.]+\)/, 'the pair with an own point carries an error');
    print $ok ? "selftest passed\n" : "selftest FAILED\n";
    return $ok;
}
