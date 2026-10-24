#!/usr/bin/env perl
#
# nyhead_neck_electrodes.pl — where are the neck electrodes of the New York Head (Nk1..Nk4),
# how close are they to the lower edge of the model, how much of the cortical leadfield do
# they carry, and which NYHead electrode is nearest to a given point (own neck / chest sites).
#
#   perl examples/nyhead_neck_electrodes.pl --mat ~/src/NYHead/sa_nyhead.mat
#   perl examples/nyhead_neck_electrodes.pl --mat ... --point 'Cv12=0,-95,-60' --point 'Th12=0,-110,-200'
#   perl examples/nyhead_neck_electrodes.pl --mat ... --point ... --rename YJ0394VB_rename.tsv --out new.tsv
#   perl examples/nyhead_neck_electrodes.pl --mat ... --arc 'Cv12=85,Cv45=125,BN1=165,Th12=200'
#   perl examples/nyhead_neck_electrodes.pl --selftest        # synthetic mesh, no .mat needed
#
#   --mat FILE       sa_nyhead.mat                                        (~/src/NYHead/sa_nyhead.mat)
#   --neck REGEX     electrodes to examine                                (^Nk\d)
#   --compare LIST   electrodes to compare the leadfield rows with        (Iz,Oz,Cz,T7,T8,Fpz)
#   --point N=X,Y,Z  a point in MNI mm (repeatable): nearest NYHead electrodes and the
#                    distance to the scalp mesh
#   --rename TSV     a channel table (tab separated, 2nd column = the new channel name);
#                    with --out, written again with the columns nyhead_nearest,
#                    nyhead_dist_mm and nyhead_point (how the point was placed: given,
#                    arc on the mesh, arc extended) for the rows given with --point or --arc
#   --out TSV        output of --rename
#   --arc LIST       NAME=MM,...  points on the back midline of the head mesh, MM measured along
#                    the skin downwards from the --from electrode (a tape measure along the
#                    neck); each is then treated like a --point
#   --from NAME      start of --arc                                       (Iz)
#   --points-out TSV write the --point / --arc points and the back midline profile (on the
#                    mesh, and its extension down to 20 mm below the deepest --arc point) as
#                    name, x, y, z, kind (point | profile_mesh | profile_extended), for
#                    nyhead_points_usda.pl
#   --strings S      auto | h5dump | lowlevel (how the names are read)    (auto)
#
# Printed:
#  1. the neck electrodes: locs_3D_orig, the nearest head-mesh vertex and its distance, and
#     the height above the lowest vertex of the head mesh (/sa/head);
#  2. the head mesh: z range, and the number of vertices below each neck electrode;
#  3. leadfield rows (V_fem_normal, as stored = average reference over the 231 electrodes):
#     RMS over the 74382 cortical sources relative to the median RMS of all electrodes, and
#     the correlation of each neck row with the --compare rows and with its nearest
#     electrode that is not a neck electrode (how much of the same brain activity it sees).
#     Then every pair of neck electrodes: distance, correlation of the rows, and the RMS of
#     the difference of the rows (the brain activity a bipolar neck lead picks up) relative
#     to the same median; for comparison the RMS of the difference between every non-neck
#     electrode and its nearest non-neck neighbour (median and range), and Cz - each neck
#     electrode (a scalp channel against a neck reference).
#     The leadfield holds cortical sources only: nothing here says how the cardiac field
#     reaches these sites.
#  4. with --arc: the back midline profile of the head mesh (the cut of /sa/head by the plane
#     x = 0, the most posterior point in every 1 mm of height below the --from electrode).
#     /sa/head ends in the upper neck and is closed by a cap: the profile is cut where it
#     first runs forward at less than 45 deg from the horizontal (the start of the cap), and
#     below that it is extended straight along the direction of its last 20 mm.  Printed:
#     the length of the profile on the mesh, the points at the given lengths (marked
#     'extended' below the cut), and for every neck electrode its length along the profile
#     (extension included) and its distance from it.  /sa/head is not the true outer skin,
#     the NYHead neck posture is not that of a seated subject, and the extension is a
#     straight line, so the points are approximate.

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Getopt::Long;
use PDL;
use PDL::IO::NYHead;

my %o = (mat => "$ENV{HOME}/src/NYHead/sa_nyhead.mat", neck => '^Nk\d', compare => 'Iz,Oz,Cz,T7,T8,Fpz',
         point => [], rename => undef, out => undef, strings => 'auto', arc => '', from => 'Iz', selftest => 0,
         points_out => undef);
GetOptions('mat=s' => \$o{mat}, 'neck=s' => \$o{neck}, 'compare=s' => \$o{compare},
           'point=s@' => $o{point}, 'rename=s' => \$o{rename}, 'out=s' => \$o{out},
           'strings=s' => \$o{strings}, 'arc=s' => \$o{arc}, 'from=s' => \$o{from}, 'selftest' => \$o{selftest},
           'points-out=s' => \$o{points_out})
    or die "bad options (see the header of $0)\n";
if ($o{selftest}) { selftest(); exit 0 }
-f $o{mat} or die "not found: $o{mat} (give it with --mat)\n";
defined $o{rename} == defined $o{out} or die "--rename and --out go together\n";

my $ny  = PDL::IO::NYHead->new($o{mat}, strings => $o{strings});
my @lab = @{ $ny->electrode_labels };
my $pos = $ny->electrode_pos->double;                       # (ne, 3) MNI mm
my %at; @at{@lab} = 0 .. $#lab;
my $hd  = $ny->surface('head');
my $vc  = $hd->{vc}->double;                                # (nv, 3)
my $vz  = $vc->slice(':,(2)');
my $zmin = $vz->min;
printf "%s: %d electrodes, head mesh %d vertices, z %.1f .. %.1f mm; strings read by %s\n",
    $o{mat}, scalar @lab, $vc->dim(0), $zmin, $vz->max, $ny->strings_method;

my @neck = grep { $lab[$_] =~ /$o{neck}/ } 0 .. $#lab;
@neck or print "no electrode matches /$o{neck}/\n";

# ---- 1. neck electrodes and the mesh -----------------------------------------
if (@neck) {
    my $np = $pos->dice_axis(0, pdl(long, \@neck));
    my $iv = $ny->nearest_vertex($np, 'head');
    $iv = pdl(long, [$iv]) unless ref $iv;
    print "\n-- 1. neck electrodes (MNI mm): position, nearest head-mesh vertex, distance, height above the mesh bottom\n";
    for my $k (0 .. $#neck) {
        my $p = $np->slice("($k),");
        my $v = $vc->slice("(" . $iv->at($k) . "),");
        printf "  %-5s (%6.1f %6.1f %6.1f)  vertex %6d (%6.1f %6.1f %6.1f)  %5.1f mm  %+6.1f mm\n",
            $lab[ $neck[$k] ], $p->list, $iv->at($k), $v->list, sqrt((($p - $v)**2)->sum), $p->at(2) - $zmin;
    }
    # ---- 2. mesh below them
    print "\n-- 2. head-mesh vertices below each neck electrode (of ", $vc->dim(0), ")\n";
    printf "  %-5s %7d vertices below z = %.1f; within 30 mm horizontally of it: %d below\n",
        $lab[ $neck[$_] ], ($vz < $np->at($_, 2))->sum, $np->at($_, 2),
        (($vz < $np->at($_, 2)) * (((($vc->slice(':,0:1') - $np->slice("($_),0:1")->dummy(0))**2)->mv(1, 0)->sumover) < 900))->sum
        for 0 .. $#neck;
}

# ---- 3. leadfield rows -----------------------------------------------------------
if (@neck) {
    my $L   = $ny->leadfield;                                # (ne, ns)
    my $rms = sqrt(($L**2)->mv(1, 0)->average);              # (ne): mean over the sources (dim 1)
    my $med = $rms->median;
    my @cmp = grep { exists $at{$_} } split /,/, $o{compare};
    my @nn  = grep { $lab[$_] !~ /$o{neck}/ } 0 .. $#lab;   # non-neck electrodes
    printf "\n-- 3. leadfield rows: RMS / median RMS of all %d electrodes, and correlation across sources\n", scalar @lab;
    printf "  %-5s %8s  %s  %s\n", '', 'RMS/med', join(' ', map { sprintf '%6s', $_ } @cmp), 'nearest non-neck electrode';
    my $Lc  = $L - $L->mv(1, 0)->average->dummy(1);         # each row centred over sources
    my $nrm = sqrt(($Lc**2)->mv(1, 0)->sumover);             # (ne)
    my $corr = sub { my ($i, $j) = @_; (($Lc->slice("($i),") * $Lc->slice("($j),"))->sum / ($nrm->at($i) * $nrm->at($j)))->sclr };
    for my $i (@neck, map { $at{$_} } @cmp) {
        my $d2 = (($pos->dice_axis(0, pdl(long, \@nn)) - $pos->slice("($i),")->dummy(0))**2)->mv(1, 0)->sumover;
        my ($j) = grep { $_ != $i } map { $nn[$_] } $d2->qsorti->list;   # nearest, not itself
        printf "  %-5s %8.3f  %s  %s (%.1f mm) r %.3f\n", $lab[$i], $rms->at($i) / $med,
            join(' ', map { sprintf '%6.3f', $corr->($i, $at{$_}) } @cmp),
            $lab[$j], sqrt((($pos->slice("($i),") - $pos->slice("($j),"))**2)->sum), $corr->($i, $j);
    }
    # differences of rows = bipolar leads
    my $drms = sub { my ($i, $j) = @_; (sqrt((($L->slice("($i),") - $L->slice("($j),"))**2)->avg) / $med)->sclr };
    my $dist = sub { my ($i, $j) = @_; sqrt((($pos->slice("($i),") - $pos->slice("($j),"))**2)->sum)->sclr };
    print "\n  pairs of neck electrodes: distance, correlation of the rows, RMS of the difference / median RMS\n";
    for my $x (0 .. $#neck) {
        for my $y ($x + 1 .. $#neck) {
            my ($i, $j) = ($neck[$x], $neck[$y]);
            printf "  %-5s - %-5s %6.1f mm  r %6.3f  diff %.3f\n", $lab[$i], $lab[$j], $dist->($i, $j), $corr->($i, $j), $drms->($i, $j);
        }
    }
    my (@nd, @ndist);                                        # every non-neck electrode minus its nearest non-neck neighbour
    my $pn = $pos->dice_axis(0, pdl(long, \@nn));
    for my $a (0 .. $#nn) {
        my $d2 = (($pn - $pn->slice("($a),")->dummy(0))**2)->mv(1, 0)->sumover;
        $d2->set($a, 9e9);
        my $b = $d2->minimum_ind->sclr;
        push @nd, $drms->($nn[$a], $nn[$b]); push @ndist, sqrt($d2->at($b));
    }
    my $ndp = pdl(\@nd);
    printf "  for comparison, each non-neck electrode - its nearest neighbour (median %.1f mm apart): diff median %.3f (%.3f .. %.3f)\n",
        pdl(\@ndist)->median, $ndp->median, $ndp->min, $ndp->max;
    printf "  Cz - %-5s diff %.3f   (a scalp channel against that neck electrode)\n", $lab[$_], $drms->($at{Cz}, $_)
        for exists $at{Cz} ? @neck : ();
}

# ---- 4. along the back midline -----------------------------------------------------
my %how;                                                         # point name -> how it was placed
my @profile_rows;                                                # for --points-out
if ($o{arc} ne '') {
    exists $at{ $o{from} } or die "--from: no electrode named $o{from}\n";
    my $p0 = $pos->slice("($at{$o{from}}),");
    my ($full, $flen) = back_profile($vc, $hd->{tri0}, $p0);     # (np, 3), (np)
    my ($prof, $len, $lmesh, $u) = extend_profile($full, $flen);
    printf "\n-- 4. back midline of the head mesh from %s (%.1f %.1f %.1f)\n", $o{from}, $p0->list;
    printf "  %s is %.1f mm off the mesh profile start (%.1f %.1f %.1f)\n", $o{from},
        sqrt((($prof->slice('(0),') - $p0)**2)->sum), $prof->slice('(0),')->list;
    printf "  on the mesh: %.1f mm down to (y %.1f, z %.1f); below that the mesh cap (the mesh profile went on for %.1f mm, to z %.1f)\n",
        $lmesh, _at_len($prof, $len, $lmesh)->at(1), _at_len($prof, $len, $lmesh)->at(2), $flen->at(-1) - $lmesh, $full->at(-1, 2);
    printf "  extended straight below that, direction (y %+.2f, z %+.2f) per mm (%.0f deg from vertical)\n",
        $u->at(1), $u->at(2), 57.29578 * atan2(abs($u->at(1)), abs($u->at(2)));
    print "  the profile every 20 mm: ", join('  ', map { my $q = _at_len($prof, $len, $_); sprintf '%d%s (%.0f %.0f)', $_, $_ > $lmesh ? '*' : '', $q->at(1), $q->at(2) }
        map { 20 * $_ } 0 .. 11), "  [mm: (y z)], * extended\n";
    for my $kv (split /,/, $o{arc}) {
        my ($n, $mm) = $kv =~ /^([^=]+)=([\d.]+)$/ or die "--arc NAME=MM, not '$kv'\n";
        my $q = _at_len($prof, $len, $mm);
        printf "  %-6s %5.0f mm: (%.1f %.1f %.1f)%s\n", $n, $mm, $q->list, $mm > $lmesh ? sprintf('  extended (%.0f mm below the mesh)', $mm - $lmesh) : '';
        push @{ $o{point} }, sprintf '%s=%.2f,%.2f,%.2f', $n, $q->list;
        $how{$n} = $mm > $lmesh ? sprintf('arc %g mm, extended', $mm) : sprintf('arc %g mm, on the mesh', $mm);
    }
    my $deep = 0;
    for my $kv (split /,/, $o{arc}) { my (undef, $mm) = split /=/, $kv; $deep = $mm if $mm > $deep }
    for my $k (0 .. $prof->dim(0) - 2) {                         # the mesh part (the last row is the extension end)
        push @profile_rows, [ sprintf('profile%03d', $k), $prof->slice("($k),")->list, 'profile_mesh' ];
    }
    push @profile_rows, [ 'profile_end', _at_len($prof, $len, $deep + 20)->list, 'profile_extended' ] if $deep + 20 > $lmesh;
    for my $i (@neck) {                                          # neck electrodes against the profile
        my ($s, $d) = project_on_profile($prof, $len, $pos->slice("($i),"));
        printf "  %-6s at %.1f mm along the profile from %s%s, %.1f mm off it\n", $lab[$i], $s, $o{from},
            $s > $lmesh ? ' (extended part)' : '', $d;
    }
}

# ---- points -> nearest electrodes --------------------------------------------------
my %pt;
if (@{ $o{point} }) {
    print "\n-- points (MNI mm): nearest NYHead electrodes, distance to the nearest head-mesh vertex\n";
    for my $s (@{ $o{point} }) {
        my ($n, $xyz) = $s =~ /^([^=]+)=(.+)$/ or die "--point NAME=X,Y,Z, not '$s'\n";
        my @c = split /,/, $xyz; @c == 3 or die "--point $n: three coordinates\n";
        my $p  = pdl(@c);
        my $d  = sqrt((($pos - $p->dummy(0))**2)->mv(1, 0)->sumover);   # (ne): sum over x,y,z (dim 1)
        my @o3 = ($d->qsorti->list)[0 .. 2];
        my $iv = $ny->nearest_vertex($p, 'head');
        my $dv = sqrt((($vc->slice("($iv),") - $p)**2)->sum);
        printf "  %-6s (%s)  %s   head mesh %.1f mm away\n", $n, join(' ', map { sprintf '%.1f', $_ } @c),
            join('  ', map { sprintf '%s %.1f mm', $lab[$_], $d->at($_) } @o3), $dv;
        $pt{$n} = [ $lab[ $o3[0] ], $d->at($o3[0]), [@c] ];
        push @{ $o{_order} }, $n;
    }
}

if (defined $o{points_out}) {
    open my $fo, '>', $o{points_out} or die "$o{points_out}: $!";
    print $fo "# points and back midline profile from nyhead_neck_electrodes.pl (MNI mm)\n";
    print $fo join("\t", qw(name x y z kind)), "\n";
    printf $fo "%s\t%.2f\t%.2f\t%.2f\t%s\n", $_, @{ $pt{$_}[2] }, 'point' for @{ $o{_order} // [] };
    printf $fo "%s\t%.2f\t%.2f\t%.2f\t%s\n", @$_ for @profile_rows;
    close $fo;
    print "\nwrote $o{points_out}\n";
}

if (defined $o{rename}) {
    open my $in,  '<', $o{rename} or die "$o{rename}: $!";
    open my $out, '>', $o{out}    or die "$o{out}: $!";
    my $seen = 0;                                            # header = the first line that is not a comment
    while (my $l = <$in>) {
        chomp $l;
        if ($l =~ /^#/ || $l eq '') { print $out "$l\n"; next }
        my @f = split /\t/, $l;
        if (!$seen++ && $f[0] =~ /^recorded$/i) { print $out join("\t", @f, 'nyhead_nearest', 'nyhead_dist_mm', 'nyhead_point'), "\n"; next }
        my $p = $pt{ $f[1] // '' };
        print $out join("\t", @f, $p ? ($p->[0], sprintf('%.1f', $p->[1]), $how{ $f[1] } // 'given') : ('', '', '')), "\n";
    }
    close $out;
    print "\nwrote $o{out}\n";
}

# =====================================================================================
# back_profile($vc (nv,3), $tri0 (nf,3) 0-based, $p0 (3)): the cut of the mesh by the plane
# x = p0.x, keeping for every 1 mm of height below p0.z the most posterior point (smallest y).
# Returns the profile (np, 3), top first, and the length along it from the height of p0 (np).
sub back_profile {
    my ($vc, $tri, $p0) = @_;
    my $x0 = $p0->at(0);
    my @pts;
    for my $e ([0, 1], [1, 2], [2, 0]) {                     # the three edges of every triangle
        my $a = $vc->dice_axis(0, $tri->slice(":,($e->[0])"));    # (nf, 3): first vertex of the edge
        my $b = $vc->dice_axis(0, $tri->slice(":,($e->[1])"));
        my ($xa, $xb) = ($a->slice(':,(0)') - $x0, $b->slice(':,(0)') - $x0);
        my $w = which($xa * $xb < 0);                         # edges that cross the plane
        next if $w->isempty;
        my $t = ($xa->index($w) / ($xa->index($w) - $xb->index($w)))->dummy(1);
        push @pts, $a->dice_axis(0, $w) + $t * ($b->dice_axis(0, $w) - $a->dice_axis(0, $w));
    }
    my $on = which(abs($vc->slice(':,(0)') - $x0) < 1e-9);   # vertices lying on the plane
    push @pts, $vc->dice_axis(0, $on) unless $on->isempty;
    @pts or die "back_profile: the plane x = $x0 does not cut the mesh\n";
    my $P = $pts[0]; $P = $P->glue(0, $_) for @pts[1 .. $#pts];   # (nc, 3): points appended along dim 0
    my ($y, $z) = ($P->slice(':,(1)'), $P->slice(':,(2)'));
    my $sel = which(($z <= $p0->at(2)) & ($y < $p0->at(1) + 60));     # below the start, back half
    my $bin = long(($p0->at(2) - $z->index($sel)));           # 1 mm bins of depth below p0
    my (%best);
    my @yl = $y->index($sel)->list; my @bl = $bin->list;
    for my $k (0 .. $#bl) {                                   # most posterior point per bin
        $best{ $bl[$k] } = $k if !exists $best{ $bl[$k] } || $yl[$k] < $yl[ $best{ $bl[$k] } ];
    }
    my @ord = map { $sel->at($best{$_}) } sort { $a <=> $b } keys %best;
    my $prof = $P->dice_axis(0, pdl(long, \@ord))->sever;    # (np, 3), top first
    # start at the height of p0: the first point moved straight up to z = p0.z (the length
    # is measured from the start electrode's height, not from the first cut point below it)
    my $top = $prof->slice('0:0,:')->copy;
    $top->slice('(0),(2)') .= $p0->at(2);
    $prof = $top->glue(0, $prof);
    my $seg  = sqrt((($prof->slice('1:-1,:') - $prof->slice('0:-2,:'))**2)->mv(1, 0)->sumover);
    my $len  = zeroes($prof->dim(0));
    $len->slice('1:-1') .= $seg->cumusumover;
    return ($prof, $len);
}

# extend_profile($prof, $len): cut the profile where the mesh cap starts (the first segment
# that runs forward, dy > 0, at less than 45 deg from the horizontal, |dy| > |dz|) and append a
# straight extension of 300 mm along the direction of the last 20 mm before the cut.
# Returns (profile with the extension, its lengths, the length on the mesh, the unit direction).
sub extend_profile {
    my ($prof, $len) = @_;
    my $d   = $prof->slice('1:-1,:') - $prof->slice('0:-2,:');            # (np-1, 3)
    my ($dy, $dz) = ($d->slice(':,(1)'), $d->slice(':,(2)'));
    my $cap = which(($dy > 0) & (abs($dy) > abs($dz)));
    my $k   = $cap->isempty ? $prof->dim(0) - 1 : $cap->at(0);             # last point on the mesh
    $k >= 1 or die "extend_profile: the profile turns into the cap at its start\n";
    my $pk  = $prof->slice("($k),");
    my $lk  = $len->at($k);
    my $pb  = _at_len($prof->slice("0:$k,:"), $len->slice("0:$k"), $lk > 20 ? $lk - 20 : 0);
    my $u   = $pk - $pb; $u /= sqrt(($u**2)->sum);
    my $ext = $prof->slice("0:$k,:")->glue(0, ($pk + 300 * $u)->dummy(0));
    my $el  = $len->slice("0:$k")->append(pdl($lk + 300));
    return ($ext->sever, $el->sever, $lk, $u);
}

# project_on_profile($prof, $len, $p): the nearest point of the polyline to $p (on a segment,
# not only at a vertex).  Returns (length along the profile, distance).
sub project_on_profile {
    my ($prof, $len, $p) = @_;
    my $A  = $prof->slice('0:-2,:'); my $B = $prof->slice('1:-1,:');
    my $AB = $B - $A;
    my $l2 = ($AB**2)->mv(1, 0)->sumover;                                  # (ns)
    my $t  = ((($p->dummy(0) - $A) * $AB)->mv(1, 0)->sumover / $l2->clip(1e-12, undef))->clip(0, 1);
    my $q  = $A + $t->dummy(1) * $AB;                                      # (ns, 3)
    my $dd = sqrt((($q - $p->dummy(0))**2)->mv(1, 0)->sumover);
    my $k  = $dd->minimum_ind->sclr;
    return ($len->at($k) + $t->at($k) * sqrt($l2->at($k)), $dd->at($k));
}

sub _at_len {                                                 # point at length $s on the profile
    my ($prof, $len, $s) = @_;
    my $k = which($len <= $s)->max;
    return $prof->slice("($k),")->copy if $k >= $len->nelem - 1;
    my $f = ($s - $len->at($k)) / ($len->at($k + 1) - $len->at($k));
    return $prof->slice("($k),") + $f * ($prof->slice("(" . ($k + 1) . "),") - $prof->slice("($k),"));
}

# synthetic mesh: a cone around the z axis, radius 50 mm at z = 0 growing to 100 mm at z = -100
# (back line y = -(50 + z/-2)), so a length s down the back from (0, -50, 0) is at depth
# s / sqrt(1 + 0.25).
sub selftest {
    my $ok = 1;
    for my $shift (0, 0.5) {          # 0: vertices lie on the plane x = 0;  0.5: only edges cross it
        print "  cone with the vertex rings turned by $shift step\n";
        my $r = _selftest_cone($shift);
        $ok &&= $r;
    }
    print $ok ? "selftest passed\n" : "selftest FAILED\n";
    exit($ok ? 0 : 1);
}
sub _selftest_cone {
    my $shift = shift;
    my ($nth, $nz) = (72, 101);
    my $th = (sequence($nth) + $shift) * (6.283185307179586 / $nth);
    my $zz = -sequence($nz);                                  # 0 .. -100
    my $r  = 50 - 0.5 * $zz;
    my $vc = zeroes($nth * $nz, 3);
    my $ri = $r->dummy(0, $nth)->flat; my $ti = $th->dummy(1, $nz)->flat;
    $vc->slice(':,(0)') .= $ri * sin($ti);
    $vc->slice(':,(1)') .= -$ri * cos($ti);                   # theta = 0 at the back (y < 0)
    $vc->slice(':,(2)') .= $zz->dummy(0, $nth)->flat;
    my @f;
    for my $j (0 .. $nz - 2) { for my $i (0 .. $nth - 1) {
        my ($a0, $a1) = ($j * $nth + $i, $j * $nth + ($i + 1) % $nth);
        my ($b0, $b1) = ($a0 + $nth, $a1 + $nth);
        push @f, [$a0, $a1, $b0], [$a1, $b1, $b0];
    } }
    # bottom cap: a flat cone from the last ring (z = -100) to a centre point at (0, 0, -110),
    # so the cut runs forward under the neck like the cap of /sa/head
    my $c = $nth * $nz;
    $vc = $vc->glue(0, pdl([[0], [0], [-110]]));              # (1, 3)
    my $j = $nz - 1;
    push @f, [$j * $nth + $_, $j * $nth + ($_ + 1) % $nth, $c] for 0 .. $nth - 1;
    my $tri = pdl(long, \@f);                                  # (3, nf) -> want (nf, 3)
    $tri = $tri->transpose if $tri->dim(0) == 3 && $tri->dim(1) != 3;
    my ($full, $flen) = back_profile($vc, $tri, pdl(0, -50, 0));
    my ($prof, $len, $lmesh) = extend_profile($full, $flen);
    my $ok = 1;
    my $chk = sub { my ($c, $m) = @_; printf "  [%s] %s\n", $c ? 'PASS' : 'FAIL', $m; $ok &&= $c };
    my $k = sqrt(1.25);
    $chk->($full->at(-1, 2) < -105, sprintf 'the raw profile runs into the cap (to z %.1f)', $full->at(-1, 2));
    $chk->(abs($lmesh - 100 * $k) < 0.2, sprintf 'cut at the cap: %.2f mm on the mesh (expected %.2f)', $lmesh, 100 * $k);
    my $qe = _at_len($prof, $len, 130);                       # extended: the cone line goes on
    my ($ez, $ey) = (-130 / $k, -(50 + 0.5 * 130 / $k));
    $chk->(abs($qe->at(2) - $ez) < 0.5 && abs($qe->at(1) - $ey) < 0.5,
        sprintf 'extended point at 130 mm: (%.2f %.2f %.2f), expected (0 %.2f %.2f)', $qe->list, $ey, $ez);
    # 10 mm off the line along its normal (y, z) = (-1, 0.5)/sqrt(1.25)
    my ($sp, $dp) = project_on_profile($prof, $len, pdl(0, $ey - 10 / $k, $ez + 5 / $k));
    $chk->(abs($sp - 130) < 0.5 && abs($dp - 10) < 0.3, sprintf 'projection of a point 10 mm off the extension: %.2f mm along, %.2f mm off (expected 130, 10)', $sp, $dp);
    my $q = _at_len($prof, $len, 85);
    $chk->(abs($q->at(2) + 85 / $k) < 0.2 && abs($q->at(1) + 50 + 0.5 * 85 / $k) < 0.3 && abs($q->at(0)) < 1e-9,
        sprintf 'point at 85 mm: (%.2f %.2f %.2f), expected (0 %.2f %.2f)', $q->list, -(50 + 0.5 * 85 / $k), -85 / $k);
    $chk->(all($prof->slice(':,(1)') < 0), 'the profile stays on the back (y < 0)');
    return $ok;
}
