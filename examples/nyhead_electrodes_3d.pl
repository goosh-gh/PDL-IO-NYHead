#!/usr/bin/env perl
#
# nyhead_electrodes_3d.pl — New York Head(sa_nyhead.mat)の 231 電極を、名前つきで 3D 表示する。
#
#   perl examples/nyhead_electrodes_3d.pl --mat ~/src/NYHead/sa_nyhead.mat                 # GS3D(giza-server)
#   perl examples/nyhead_electrodes_3d.pl --mat ... --backend trid                          # PDL::Graphics::TriD
#   perl examples/nyhead_electrodes_3d.pl --mat ... --mark 'Iz|^I|O9|O10|^N'                # 名前で強調
#   perl examples/nyhead_electrodes_3d.pl --mat ... --below 10 --list-only                  # 表だけ(窓なし)
#
#   --mat FILE       sa_nyhead.mat(MATLAB v7.3 = HDF5)                    (~/src/NYHead/sa_nyhead.mat)
#   --backend B      gs3d | trid                                              (gs3d)
#   --mark REGEX     名前がこれに合う電極を黄色で描く                         (なし)
#   --below MM       Iz の高さ + MM より下で後ろ半分(y < 0)の電極を表にし、マゼンタで描く (0)
#   --no-head        GS3D で頭皮メッシュ(/sa/head)を描かない
#   --labels / --no-labels   ラベル表示(GS3D は既定 on、TriD は既定 off: 下の注意)
#   --list-only      表だけ出して窓は開かない
#   --only REGEX     窓に描く NYHead 電極をこれに合う名前だけにする(ラベルの重なりを減らす) (全部)
#   --points TSV     自前の点を加える: name, x, y, z[, kind](MNI mm、タブ区切り。
#                    nyhead_neck_electrodes.pl --points-out の出力)。kind が point の行は
#                    名前つきの点(ピンク系)として、profile_mesh / profile_extended の行は
#                    首の後ろの正中の線(緑 = メッシュ上、橙 = まっすぐの延長)として描く
#   --point N=X,Y,Z  自前の点をもう 1 つ(複数可)
#   --link           自前の点と最も近い NYHead 電極を細線で結ぶ(既定 on、--no-link で外す)
#   --surface NAME   GS3D の面: head | cortex75K | cortex10K | ... (head)。面は 1 枚だけ
#   --strings S      電極名の読み方 auto | h5dump | lowlevel                  (auto)
#
# 読み込みは PDL::IO::NYHead に任せる:
#   electrode_labels(231 名。h5dump があれば h5dump、無ければ PDL::IO::HDF5 の低水準関数)、
#   electrode_pos(/sa/locs_3D_orig、(231,3) MNI mm)、surface('head')(頭皮メッシュ、tri0 = 0 始まり)。
#
# 色: 左=オレンジ、右=シアン、正中(|x| < 5 mm)=白、--below の対象=マゼンタ、--mark=黄、
#     自前の点(--points / --point)=ピンク系。
# GS3D の窓の File → Export 3D as OBJ/USDA… で、表示中のシーン(面、点、ラベルの文字)を
# .usda に書き出せる(PDL::Graphics::Cairo::SceneExport。Cocoa の giza-server のみ)。
# 表示の流儀は PDL-EEG の examples/show_electrodes_3d.pl と同じ(TriD = Advent 2024 Day 12 の
# points3d + Labels、GS3D = PDL::Graphics::Cairo::Driver::GS3D)。
#
use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Getopt::Long;
use PDL;
use PDL::IO::NYHead;

my %o = (mat => "$ENV{HOME}/src/NYHead/sa_nyhead.mat", backend => 'gs3d',
         mark => undef, below => 0, head => 1, labels => undef, list_only => 0, strings => 'auto',
         only => undef, points => [], point => [], link => 1, surface => 'head');
GetOptions('mat=s' => \$o{mat}, 'backend=s' => \$o{backend},
           'mark=s' => \$o{mark}, 'below=f' => \$o{below}, 'head!' => \$o{head},
           'labels!' => \$o{labels}, 'list-only' => \$o{list_only}, 'strings=s' => \$o{strings},
           'only=s' => \$o{only}, 'points=s@' => $o{points}, 'point=s@' => $o{point},
           'link!' => \$o{link}, 'surface=s' => \$o{surface})
    or die "bad options (see the header of $0)\n";
-f $o{mat} or die "not found: $o{mat} (give it with --mat)\n";
$o{backend} = lc $o{backend};
$o{backend} =~ /^(gs3d|trid)$/ or die "--backend must be gs3d or trid\n";

my $ny  = PDL::IO::NYHead->new($o{mat}, strings => $o{strings});

# ---- labels and positions ------------------------------------------------------------
my @lab = @{ $ny->electrode_labels };
my $xyz = $ny->electrode_pos->double->sever;                 # (n, 3) MNI mm
my $n   = @lab;
$xyz->dim(0) == $n or die sprintf "%d labels but %d positions\n", $n, $xyz->dim(0);
my ($X, $Y, $Z) = map { $xyz->slice(":,($_)") } 0 .. 2;
printf "%s: %d electrodes, positions /sa/locs_3D_orig (MNI mm)\n", $o{mat}, $n;

# ---- colours ---------------------------------------------------------------------------
my %at; @at{@lab} = 0 .. $#lab;
my $rgb = zeroes(double, $n, 3);
$rgb .= pdl(1.0, 1.0, 1.0)->dummy(0, $n);                       # midline: white
my $left  = which($X < -5);
my $right = which($X >  5);
$rgb->dice_axis(0, $left)  .= pdl(0.9, 0.55, 0.25)->dummy(0, $left->nelem)  if $left->nelem;
$rgb->dice_axis(0, $right) .= pdl(0.3, 0.7,  0.85)->dummy(0, $right->nelem) if $right->nelem;

my @low;
if ($o{below}) {
    exists $at{Iz} or die "--below needs an electrode named Iz in clab_electrodes\n";
    my $zl = $Z->at($at{Iz}) + $o{below};
    @low = sort { $Z->at($a) <=> $Z->at($b) } which(($Z < $zl) & ($Y < 0))->list;      # back half only (neck, not face)
    printf "\n%d electrode(s) behind y = 0 and below z = %.1f mm (Iz %.1f mm %+g mm), lowest first:\n",
        scalar @low, $zl, $Z->at($at{Iz}), $o{below};
    printf "  %-8s %7.1f %7.1f %7.1f\n", $lab[$_], $X->at($_), $Y->at($_), $Z->at($_) for @low;
    $rgb->dice_axis(0, pdl(long, \@low)) .= pdl(1.0, 0.2, 0.9)->dummy(0, scalar @low) if @low;
}
if (defined $o{mark}) {
    my @m = grep { $lab[$_] =~ /$o{mark}/ } 0 .. $n - 1;
    printf "\n%d electrode(s) matching /%s/:\n", scalar @m, $o{mark};
    printf "  %-8s %7.1f %7.1f %7.1f\n", $lab[$_], $X->at($_), $Y->at($_), $Z->at($_) for @m;
    $rgb->dice_axis(0, pdl(long, \@m)) .= pdl(1.0, 0.95, 0.1)->dummy(0, scalar @m) if @m;
}
print "\nall labels:\n";
for (my $i = 0; $i < $n; $i += 12) {
    my $j = $i + 11 < $n ? $i + 11 : $n - 1;
    print '  ', join(' ', map { sprintf '%-7s', $_ } @lab[$i .. $j]), "\n";
}

# ---- own points, the back-midline profile, links -----------------------------------------
my (@own, @prof);                                  # [name, x, y, z]; [x, y, z, mesh|extended]
for my $f (@{ $o{points} }) {
    open my $fh, '<', $f or die "$f: $!";
    while (<$fh>) {
        chomp; next if /^\s*(#|$)/;
        my ($nm, $x, $y, $z, $k) = split /\t/;
        next if $nm eq 'name';
        $k //= 'point';
        if    ($k eq 'point')                     { push @own,  [$nm, $x, $y, $z] }
        elsif ($k =~ /^profile_(mesh|extended)$/) { push @prof, [$x, $y, $z, $1] }
        else  { die "$f: unknown kind '$k'\n" }
    }
}
for my $s (@{ $o{point} }) {
    my ($nm, $c) = $s =~ /^([^=]+)=(.+)$/ or die "--point NAME=X,Y,Z, not '$s'\n";
    my @c = split /,/, $c; @c == 3 or die "--point $nm: three coordinates\n";
    push @own, [$nm, @c];
}
my @pal = ([1.0, 0.25, 0.70], [0.95, 0.10, 0.45], [1.0, 0.55, 0.85], [0.80, 0.20, 0.90], [1.0, 0.40, 0.40]);
my @lines;                                         # rows (x0,y0,z0, x1,y1,z1, r,g,b)
if (@own) {
    print "\nown points: nearest NYHead electrode\n";
    for my $k (0 .. $#own) {
        my ($nm, @c) = @{ $own[$k] };
        my $d = sqrt((($xyz - pdl(@c)->dummy(0))**2)->mv(1, 0)->sumover);
        my $j = $d->minimum_ind->sclr;
        printf "  %-6s (%7.1f %7.1f %7.1f)  %s %.1f mm\n", $nm, @c, $lab[$j], $d->at($j);
        push @lines, [ @c, $xyz->slice("($j),")->list, @{ $pal[ $k % @pal ] } ] if $o{link};
    }
}
for my $kind (['mesh', [0.10, 0.85, 0.25]], ['extended', [1.0, 0.55, 0.10]]) {
    my @p = map { [ @$_[0 .. 2] ] } grep { $_->[3] eq $kind->[0] } @prof;
    unshift @p, [ map { $_ } @{ (grep { $_->[3] eq 'mesh' } @prof)[-1] }[0 .. 2] ] if $kind->[0] eq 'extended' && @p && grep { $_->[3] eq 'mesh' } @prof;
    push @lines, [ @{ $p[$_] }, @{ $p[$_ + 1] }, @{ $kind->[1] } ] for 0 .. $#p - 1;
}
exit 0 if $o{list_only};

# what the window shows: the NYHead electrodes (all, or --only), then the own points
my @show = defined $o{only} ? grep { $lab[$_] =~ /$o{only}/ } 0 .. $n - 1 : 0 .. $n - 1;
@show or die "--only /$o{only}/ matches no electrode\n";
my @plab = (@lab[@show], map { $_->[0] } @own);
my $pxyz = $xyz->dice_axis(0, pdl(long, \@show));
my $prgb = $rgb->dice_axis(0, pdl(long, \@show));
if (@own) {
    $pxyz = $pxyz->glue(0, pdl([ map { [ @$_[1 .. 3] ] } @own ])->transpose);
    $prgb = $prgb->glue(0, pdl([ map { $pal[ $_ % @pal ] } 0 .. $#own ])->transpose);
}
my $plines = @lines ? pdl(\@lines)->transpose->sever : undef;     # (M, 9)

$o{backend} eq 'gs3d' ? show_gs3d() : show_trid();

# =========================================================================================
# GS3D — PDL::Graphics::Cairo::Driver::GS3D (giza-server 3D viewer), with the scalp mesh
# =========================================================================================
sub show_gs3d {
    eval { require PDL::Graphics::Cairo::Driver::GS3D; 1 }
        or die "GS3D backend needs PDL::Graphics::Cairo (Driver::GS3D) in \@INC:\n$@";
    my %scene = (points => $pxyz, colors => $prgb,
                 labels => (($o{labels} // 1) ? \@plab : undef));
    $scene{lines} = $plines if defined $plines;
    my @pivot;
    if ($o{head} && $o{surface} ne 'none') {
        my $hd = $o{surface} eq 'head' ? $ny->surface('head') : $ny->cortex($o{surface});
        my ($vc, $tri) = ($hd->{vc}, $hd->{tri0});
        $scene{mesh} = { verts => $vc->double, faces => $tri->long,     # 0-based
                         color => [0.86, 0.78, 0.70], light => [0.3, 0.45, 1.0], ambient => 0.28,
                         winding => 'auto', fill_step => 1 };
        @pivot = (pivot => [ $vc->double->average->list ]);
        printf "mesh %s: %d vertices, %d faces\n", $o{surface}, $vc->dim(0), $tri->dim(0);
    }
    my $drv = PDL::Graphics::Cairo::Driver::GS3D->new(scene => \%scene, width => 760, height => 760, @pivot);
    print "GS3D: drag to rotate, r = reset view, close the window to quit\n";
    $drv->run;
}

# =========================================================================================
# TriD — Advent Calendar 2024 Day 12 path (points3d + Labels)
# =========================================================================================
sub show_trid {
    eval "use PDL::Graphics::TriD; 1" or die "TriD backend needs PDL::Graphics::TriD:\n$@";
    # Labels: the current TriD master ships no Labels.pm, and an old classic-OpenGL
    # Labels.pm collides with OpenGL::Modern (see PDL-EEG examples/show_electrodes_3d.pl), so only on --labels.
    my $want = $o{labels} // 0;
    if ($want && !eval "use PDL::Graphics::TriD::Labels; 1") {
        warn "TriD labels unavailable, rendering without them:\n$@";
        $want = 0;
    }
    my @c = map { $prgb->slice(":,($_)")->sever } 0 .. 2;
    my ($PX, $PY, $PZ) = map { $pxyz->slice(":,($_)")->sever } 0 .. 2;
    print "TriD: drag to rotate; press 'q' in the window to close\n";
    PDL::Graphics::TriD::points3d([ $PX, $PY, $PZ ], \@c, { PointSize => 8 });
    if ($want) {
        PDL::Graphics::TriD::hold3d();
        PDL::Graphics::TriD::graph_object(
            PDL::Graphics::TriD::Labels->new([ $PX, $PY, $PZ ],
                                             { Strings => [ map { "  $_" } @plab ] }));
        PDL::Graphics::TriD::release3d();
    }
}
