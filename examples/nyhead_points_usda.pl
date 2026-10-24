#!/usr/bin/env perl
#
# nyhead_points_usda.pl — write one .usda with the New York Head scalp (wireframe), the cortex,
# the NYHead electrodes and your own points (e.g. neck / chest electrode sites), so that the
# position of a point relative to the model electrodes (Nk1..Nk4 ...) can be seen in usdview,
# Keynote or Blender.
#
#   perl examples/nyhead_points_usda.pl --mat ~/src/NYHead/sa_nyhead.mat --points neck_points.tsv -o neck.usda
#   perl examples/nyhead_points_usda.pl --mat ... --point 'BN2=0,40,-215' --points neck_points.tsv -o neck.usda
#   perl examples/nyhead_points_usda.pl --self-test -o /tmp/st.usda       # synthetic geometry, no .mat
#
#   --mat FILE       sa_nyhead.mat                                     (~/src/NYHead/sa_nyhead.mat)
#   --points TSV     points to draw: name, x, y, z[, kind] in MNI mm, tab separated, '#' comments
#                    (nyhead_neck_electrodes.pl --points-out writes it). kind 'point' (default)
#                    is drawn as a coloured sphere; 'profile_mesh' / 'profile_extended' rows are
#                    joined in file order into the back-midline line (green on the mesh, orange
#                    for the straight extension below it)
#   --point N=X,Y,Z  one more point (repeatable)
#   --cortex NAME    cortex mesh: cortex75K, cortex10K, cortex5K, ... or none     (cortex10K)
#   --mark REGEX     NYHead electrodes drawn as opaque yellow spheres            (^Nk\d$|^Iz$)
#   --no-electrodes  leave the other NYHead electrodes out (by default small translucent spheres)
#   --link           a thin line from every point to its nearest NYHead electrode (on by default;
#                    --no-link to leave out); the prim name carries the distance, e.g.
#                    link_Cv45_Nk2_9p7mm
#   --axes MM        X/Y/Z axis triad at the MNI origin (red/green/blue), 0 = none  (60)
#   -o FILE          output                                             (nyhead_points.usda)
#   --strings S      auto | h5dump | lowlevel (how the names are read)  (auto)
#
# The file is Z-up in millimetres (upAxis = "Z", metersPerUnit = 0.001).  Every prim has the
# name of what it shows (electrode label, point name), so the outline of usdview identifies it.
# A colour legend for the points is printed.  /sa/head ends in the upper neck (a closed cap
# below), so points lower than that are drawn outside the scalp wireframe.

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Getopt::Long;
use PDL;

my %o = (mat => "$ENV{HOME}/src/NYHead/sa_nyhead.mat", points => [], point => [], cortex => 'cortex10K',
         mark => '^Nk\d$|^Iz$', electrodes => 1, link => 1, axes => 60, out => 'nyhead_points.usda',
         strings => 'auto', self_test => 0);
GetOptions('mat=s' => \$o{mat}, 'points=s@' => $o{points}, 'point=s@' => $o{point}, 'cortex=s' => \$o{cortex},
           'mark=s' => \$o{mark}, 'electrodes!' => \$o{electrodes}, 'link!' => \$o{link}, 'axes=f' => \$o{axes},
           'out|o=s' => \$o{out}, 'strings=s' => \$o{strings}, 'self-test' => \$o{self_test})
    or die "bad options (see the header of $0)\n";

my %G = $o{self_test} ? synthetic() : from_nyhead();

# ---- points ---------------------------------------------------------------------
my (@pts, @prof);                                   # [name, x, y, z], profile rows [x, y, z, kind]
for my $f (@{ $o{points} }) {
    open my $fh, '<', $f or die "$f: $!";
    while (<$fh>) {
        chomp; next if /^\s*(#|$)/;
        my ($n, $x, $y, $z, $k) = split /\t/;
        next if $n eq 'name';                       # header
        $k //= 'point';
        if ($k eq 'point') { push @pts, [$n, $x, $y, $z] }
        elsif ($k =~ /^profile_(mesh|extended)$/) { push @prof, [$x, $y, $z, $1] }
        else { die "$f: unknown kind '$k'\n" }
    }
}
for my $s (@{ $o{point} }) {
    my ($n, $xyz) = $s =~ /^([^=]+)=(.+)$/ or die "--point NAME=X,Y,Z, not '$s'\n";
    my @c = split /,/, $xyz; @c == 3 or die "--point $n: three coordinates\n";
    push @pts, [$n, @c];
}

# ---- write ------------------------------------------------------------------------
my @pal = ([0.90, 0.10, 0.60], [0.10, 0.75, 0.95], [0.20, 0.80, 0.20], [0.95, 0.45, 0.10],
           [0.55, 0.30, 0.95], [0.95, 0.85, 0.20], [0.10, 0.45, 0.85], [0.85, 0.20, 0.20]);
open my $u, '>', $o{out} or die "$o{out}: $!";
print $u qq~#usda 1.0\n(\n    defaultPrim = "NYHead"\n    upAxis = "Z"\n    metersPerUnit = 0.001\n)\n\n~;
print $u qq~def Xform "NYHead"\n{\n~;

# scalp: every triangle edge once, as a wireframe
{
    my ($vc, $tri) = ($G{skin}{vc}, $G{skin}{tri0});
    my $nv = $vc->dim(0);
    my $e  = cat(map { my ($a0, $a1) = @$_;
                       my ($p, $q) = ($tri->slice(":,($a0)"), $tri->slice(":,($a1)"));
                       my $lo = ($p < $q) * $p + ($p >= $q) * $q;
                       my $hi = ($p < $q) * $q + ($p >= $q) * $p;
                       $lo->longlong * $nv + $hi } [0, 1], [1, 2], [2, 0])->flat->uniq;    # edge codes
    my ($ea, $eb) = (($e / $nv)->long, ($e % $nv)->long);
    my $A = $vc->dice_axis(0, $ea); my $B = $vc->dice_axis(0, $eb);
    my $ne = $e->nelem;
    print $u qq~    def BasisCurves "scalp_wire"\n    {\n        uniform token type = "linear"\n~;
    print $u "        int[] curveVertexCounts = [", join(',', (2) x $ne), "]\n";
    print $u "        point3f[] points = [", join(',', map { sprintf '(%.2f,%.2f,%.2f),(%.2f,%.2f,%.2f)', $A->slice("($_),")->list, $B->slice("($_),")->list } 0 .. $ne - 1), "]\n";
    print $u qq~        float[] widths = [0.35] (interpolation = "constant")\n~;
    print $u qq~        color3f[] primvars:displayColor = [(0.80, 0.72, 0.62)]\n    }\n~;
    printf "scalp /sa/head: %d vertices, %d edges as a wireframe, z %.1f .. %.1f mm\n",
        $nv, $ne, $vc->slice(':,(2)')->minmax;
}

# cortex: a solid mesh
if ($G{cortex}) {
    my ($vc, $tri) = ($G{cortex}{vc}, $G{cortex}{tri0});
    my $nv = $vc->dim(0); my $nf = $tri->dim(0);
    print $u qq~    def Mesh "cortex"\n    {\n        uniform bool doubleSided = 1\n~;
    print $u "        int[] faceVertexCounts = [", join(',', (3) x $nf), "]\n";
    print $u "        int[] faceVertexIndices = [", join(',', $tri->transpose->flat->list), "]\n";
    print $u "        point3f[] points = [", join(',', map { sprintf '(%.2f,%.2f,%.2f)', $vc->slice("($_),")->list } 0 .. $nv - 1), "]\n";
    print $u qq~        color3f[] primvars:displayColor = [(0.92, 0.66, 0.66)]\n    }\n~;
    printf "cortex %s: %d vertices, %d faces\n", $G{cortex}{name}, $nv, $nf;
}

# NYHead electrodes
my @lab = @{ $G{labels} };
my $pos = $G{pos};
my %used;
print $u qq~    def Xform "electrodes"\n    {\n~;
for my $i (0 .. $#lab) {
    my $marked = $lab[$i] =~ /$o{mark}/;
    next unless $marked || $o{electrodes};
    print $u _sphere(_prim("e_$lab[$i]", \%used), [ $pos->slice("($i),")->list ],
                     $marked ? (4.0, [1.0, 0.92, 0.10], 1.0) : (2.0, [0.75, 0.75, 0.80], 0.35), '        ');
}
print $u "    }\n";

# own points, links to the nearest electrode
print $u qq~    def Xform "points"\n    {\n~;
print "points (sphere colour, nearest NYHead electrode):\n";
for my $k (0 .. $#pts) {
    my ($n, @c) = @{ $pts[$k] };
    my $col = $pal[ $k % @pal ];
    print $u _sphere(_prim("pt_$n", \%used), \@c, 5.0, $col, 1.0, '        ');
    my $d  = sqrt((($pos - pdl(@c)->dummy(0))**2)->mv(1, 0)->sumover);
    my $j  = $d->minimum_ind->sclr;
    printf "  %-6s (%7.1f %7.1f %7.1f)  rgb(%.2f %.2f %.2f)  nearest %s %.1f mm\n", $n, @c, @$col, $lab[$j], $d->at($j);
    next unless $o{link};
    (my $dm = sprintf '%.1f', $d->at($j)) =~ s/\./p/;
    print $u _line(_prim("link_${n}_$lab[$j]_${dm}mm", \%used), [ [@c], [ $pos->slice("($j),")->list ] ], 0.8, $col, '        ');
}
print $u "    }\n";

# back-midline profile
if (@prof) {
    my @mesh = map { [ @$_[0 .. 2] ] } grep { $_->[3] eq 'mesh' } @prof;
    my @ext  = map { [ @$_[0 .. 2] ] } grep { $_->[3] eq 'extended' } @prof;
    print $u _line('profile_on_mesh', \@mesh, 1.2, [0.10, 0.80, 0.25], '    ') if @mesh >= 2;
    print $u _line('profile_extended', [ $mesh[-1], @ext ], 1.2, [1.00, 0.55, 0.10], '    ') if @mesh && @ext;
    printf "back-midline profile: %d points on the mesh (green)%s\n", scalar @mesh, @ext ? ', straight extension (orange)' : '';
}

if ($o{axes} > 0) {
    my $L = $o{axes};
    print $u qq~    def BasisCurves "axes"\n    {\n        uniform token type = "linear"\n        int[] curveVertexCounts = [2, 2, 2]\n~;
    print $u "        point3f[] points = [(0,0,0),($L,0,0),(0,0,0),(0,$L,0),(0,0,0),(0,0,$L)]\n";
    print $u qq~        color3f[] primvars:displayColor = [(1,0,0),(0,1,0),(0,0,1)] (interpolation = "uniform")\n~;
    print $u qq~        float[] widths = [1.0] (interpolation = "constant")\n    }\n~;
}
print $u "}\n";
close $u;
print "wrote $o{out}\n";

# =====================================================================================
sub _prim {                                         # a valid, unique USD prim name
    my ($n, $used) = @_;
    $n =~ s/[^A-Za-z0-9_]/_/g;
    $n = "_$n" if $n =~ /^\d/;
    my $m = $n; my $k = 1;
    $m = $n . '_' . $k++ while $used->{$m}++;
    return $m;
}
sub _sphere {
    my ($name, $c, $r, $col, $op, $ind) = @_;
    return sprintf qq~%sdef Sphere "%s"\n%s{\n%s    double radius = %g\n%s    double3 xformOp:translate = (%.2f, %.2f, %.2f)\n~
                 . qq~%s    uniform token[] xformOpOrder = ["xformOp:translate"]\n%s    color3f[] primvars:displayColor = [(%.3f, %.3f, %.3f)]\n~
                 . qq~%s    float[] primvars:displayOpacity = [%g]\n%s}\n~,
        $ind, $name, $ind, $ind, $r, $ind, @$c, $ind, $ind, @$col, $ind, $op, $ind;
}
sub _line {
    my ($name, $p, $w, $col, $ind) = @_;
    return sprintf qq~%sdef BasisCurves "%s"\n%s{\n%s    uniform token type = "linear"\n%s    int[] curveVertexCounts = [%d]\n~
                 . qq~%s    point3f[] points = [%s]\n%s    float[] widths = [%g] (interpolation = "constant")\n~
                 . qq~%s    color3f[] primvars:displayColor = [(%.3f, %.3f, %.3f)]\n%s}\n~,
        $ind, $name, $ind, $ind, $ind, scalar @$p, $ind, join(',', map { sprintf '(%.2f,%.2f,%.2f)', @$_ } @$p),
        $ind, $w, $ind, @$col, $ind;
}

sub from_nyhead {
    -f $o{mat} or die "not found: $o{mat} (give it with --mat)\n";
    require PDL::IO::NYHead;
    my $ny = PDL::IO::NYHead->new($o{mat}, strings => $o{strings});
    my %g  = (labels => $ny->electrode_labels, pos => $ny->electrode_pos->double);
    my $s  = $ny->surface('head');
    $g{skin} = { vc => $s->{vc}->double, tri0 => long($s->{tri0}) };
    if ($o{cortex} ne 'none') {
        my $c = $ny->cortex($o{cortex});
        $g{cortex} = { vc => $c->{vc}->double, tri0 => long($c->{tri0}), name => $o{cortex} };
    }
    return %g;
}

# synthetic geometry for --self-test: an icosphere-like UV sphere of radius 90 (scalp),
# radius 70 (cortex), electrodes on the scalp sphere
sub synthetic {
    my $uv = sub {
        my ($r, $nt, $np) = @_;
        my @v; my @f;
        for my $i (0 .. $nt) { my $th = 3.14159265 * $i / $nt;
            for my $j (0 .. $np - 1) { my $ph = 6.2831853 * $j / $np;
                push @v, [ $r * sin($th) * cos($ph), $r * sin($th) * sin($ph), $r * cos($th) ] } }
        for my $i (0 .. $nt - 1) { for my $j (0 .. $np - 1) {
            my ($a0, $a1) = ($i * $np + $j, $i * $np + ($j + 1) % $np);
            push @f, [$a0, $a1, $a0 + $np], [$a1, $a1 + $np, $a0 + $np] } }
        return { vc => pdl(\@v)->transpose->sever, tri0 => long(pdl(\@f))->transpose->sever };   # (N,3), (nf,3)
    };
    my %g = (skin => $uv->(90, 12, 16), labels => [qw(Iz Nk2 Cz 1Bad)],
             pos => pdl([[0, -90, 0], [0, -80, -60], [0, 0, 90], [60, 0, 60]])->transpose->sever);   # (4,3)
    $g{cortex} = { %{ $uv->(70, 8, 12) }, name => 'synthetic' } unless $o{cortex} eq 'none';
    return %g;
}
