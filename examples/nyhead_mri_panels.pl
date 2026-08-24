#!/usr/bin/perl
# nyhead_mri_panels.pl — NYHead MRI を axial/coronal/sagittal の3断面で描き、
# 脳(brain window の灰階調)の上に骨(T1 の暗い側の window)を色付きで重畳する
# (A案 = ITK-SNAP 流の intensity-window overlay)。断面には電極(EEC/Ex19/LPA 等)を
# 面外距離つきで載せる。
#
# 使い方(実データ):
#   perl nyhead_mri_panels.pl --mat sa_nyhead.mat \
#        --eec -65,-25,-60 --nyhead Ex19,LPA --out panels.png
#   perl nyhead_mri_panels.pl --mat sa_nyhead.mat --center 0,0,0 \
#        --elec 'Ex19=-77,-24,-59' --elec 'MYPT=-65,-25,-60'
#
# 使い方(合成 self-test, .mat 不要):
#   perl nyhead_mri_panels.pl --self-test --out selftest.png
#
# 座標系: NYHead の MRI は MNI 枠(nyhead_mri_to_nii.pl で確定, sform_code=4)。
#   vol は PDL 順 (Ni,Nj,Nk)=(394,466,378), 軸意味 i=x(L->R,+i=+R) /
#   j=y(P->A,+j=+A) / k=z(I->S,+k=+S)。電極も MNI mm なので world->voxel は
#   mni2mri の一発変換。

use strict; use warnings; use PDL;
use lib '/Users/goosh/src/PDL_IO_NYHead/lib/';
use PDL::Image2D;   # cc8compt(内部空気の分離)
use Getopt::Long;
# PDL::IO::HDF5 は実データ .mat 読取時のみ require する

# ---------------------------------------------------------------- opts
my %o = (
    mat          => undef,
    selftest     => 0,
    center       => undef,          # "x,y,z" MNI mm
    eec          => undef,          # "x,y,z" MNI mm (外耳道電極点)
    nyhead       => undef,          # "Ex19,LPA" NYHead 電極名(.mat から座標取得)
    elec         => [],             # 'LABEL=x,y,z' を複数
    brain_window => undef,          # "lo,hi"
    bone_window  => undef,          # "lo,hi"
    bone_color   => "1,0.5,0",      # 骨帯オーバレイ色 0..1
    bone_opacity => 0.55,
    overlay      => "air",          # オーバレイ種別: air|bone|none
    air_side     => "auto",         # 空気の輝度側: auto|low|high(auto=背景の極性)
    air_color    => "0,0.9,1",      # 内部空気の色(シアン)
    radiological => 0,              # 既定 neurological(患者右=画像右)
    slab         => 4.0,            # 面外 |off|<=slab mm で「面上」扱い=塗り
    max_dim      => 340,            # 表示スライスの最大辺(subsample で軽量化)
    elec_orig    => 1,              # NYHead 電極: 1=locs_3D_orig / 0=登録済 locs_3D
    bg           => "auto",         # 背景の扱い: auto|low|high|none(高輝度背景を黒に)
    invert       => 0,              # grayscale 反転(組織を明るく)
    out          => "nyhead_mri_panels.png",
    show         => 0,
);
GetOptions(
    "mat=s"          => \$o{mat},
    "self-test!"     => \$o{selftest},
    "center=s"       => \$o{center},
    "eec=s"          => \$o{eec},
    "nyhead=s"       => \$o{nyhead},
    "elec=s@"        => $o{elec},
    "brain-window=s" => \$o{brain_window},
    "bone-window=s"  => \$o{bone_window},
    "bone-color=s"   => \$o{bone_color},
    "bone-opacity=f" => \$o{bone_opacity},
    "overlay=s"      => \$o{overlay},
    "air-side=s"     => \$o{air_side},
    "air-color=s"    => \$o{air_color},
    "radiological!"  => \$o{radiological},
    "slab=f"         => \$o{slab},
    "max-dim=i"      => \$o{max_dim},
    "elec-orig!"     => \$o{elec_orig},
    "bg=s"           => \$o{bg},
    "invert!"        => \$o{invert},
    "out=s"          => \$o{out},
    "show!"          => \$o{show},
) or die "bad options\n";

die "give --mat FILE or --self-test\n" unless $o{mat} || $o{selftest};

# ---------------------------------------------------------------- helpers
sub parse_xyz {
    my ($s) = @_;
    my @v = split /\s*,\s*/, $s;
    die "expected x,y,z got '$s'\n" unless @v == 3;
    return map { 0 + $_ } @v;
}

# 数値配列(1D piddle)の percentile。PDL::Stats に依存しない。
sub pctl {
    my ($p1d, $p) = @_;
    my $s = $p1d->flat->qsort;
    my $n = $s->nelem;
    return 0 if $n == 0;
    my $i = int($p * ($n - 1) + 0.5);
    $i = 0        if $i < 0;
    $i = $n - 1   if $i > $n - 1;
    return $s->at($i);
}

# math 規約の 4x4 affine を点(x,y,z)に適用 → (o0,o1,o2)
sub apply_affine {
    my ($M, $x, $y, $z) = @_;
    my @o;
    for my $r (0 .. 2) {
        $o[$r] = $M->at($r,0)*$x + $M->at($r,1)*$y + $M->at($r,2)*$z + $M->at($r,3);
    }
    return @o;
}

# 数学規約(最下行 [0,0,0,1])を満たす向きを自動選択して返す。
# PDL::IO::HDF5 の軸反転で転置されて読まれても、平行移動列が正しい方を選ぶ。
sub to_math_affine {
    my ($A) = @_;
    my ($best, $bestdev) = (undef, 1e30);
    for my $M ($A->sever, $A->transpose->sever) {
        my $bottom = abs($M->at(3,0)) + abs($M->at(3,1)) + abs($M->at(3,2));
        my $one    = abs($M->at(3,3) - 1);
        return $M if $bottom < 1e-3 && $one < 1e-3;
        my $dev = $bottom + $one;
        ($best, $bestdev) = ($M, $dev) if $dev < $bestdev;
    }
    warn "to_math_affine: no clean [0,0,0,1] bottom row (dev=$bestdev); using closest\n";
    return $best;
}

# HDF5 の group をたどって dataset を get
sub h5get {
    my ($h5, $path) = @_;
    my @p = grep { length } split m{/}, $path;
    my $ds = pop @p;
    my $node = $h5;
    $node = $node->group($_) for @p;
    return $node->dataset($ds)->get;
}

# 2D スライスの空気マスクを外部背景/内部空洞に分離。
# is_high=1 なら高輝度側が空気, 0 なら低輝度側。thr は境界しきい値。
# 返り値: ($external, $internal_air)  (どちらも byte 2D)
sub air_components {
    my ($S, $is_high, $thr) = @_;
    my $air = ($is_high ? ($S > $thr) : ($S < $thr))->byte;
    return (zeroes(byte,$S->dims), zeroes(byte,$S->dims)) if $air->sum == 0;
    my $lab = cc8compt($air);
    my ($H,$W) = $lab->dims;
    my $border = $lab->slice("0,:")->flat
               ->append($lab->slice("-1,:")->flat)
               ->append($lab->slice(":,0")->flat)
               ->append($lab->slice(":,-1")->flat);
    my $ext = zeroes(byte,$H,$W);
    for my $L (grep { $_>0 } $border->uniq->list) { $ext = $ext | ($lab==$L); }
    my $intair = ($air & ($ext==0))->byte;
    return ($ext->byte, $intair);
}

# 電極 voxel から最寄りの「空気」ボクセルまでの距離(mm)。+-35mm の箱内で探す。
sub nearest_air_mm {
    my ($ci,$cj,$ck, $is_high, $thr, $Ni,$Nj,$Nk, $vsize_ref, $vol) = @_;
    my @c = ($ci,$cj,$ck); my @N = ($Ni,$Nj,$Nk); my @vs = @$vsize_ref;
    my $rv = 70;   # +-35mm @0.5mm
    my @lo = map { my $x=$c[$_]-$rv; $x<0?0:$x } (0,1,2);
    my @hi = map { my $x=$c[$_]+$rv; $x>$N[$_]-1?$N[$_]-1:$x } (0,1,2);
    my $box = $vol->slice("$lo[0]:$hi[0],$lo[1]:$hi[1],$lo[2]:$hi[2]");
    my $air = $is_high ? ($box > $thr) : ($box < $thr);
    return undef if $air->sum == 0;
    my ($bi,$bj,$bk) = $box->dims;
    my $dx = (xvals($bi,$bj,$bk) + $lo[0] - $ci) * $vs[0];
    my $dy = (yvals($bi,$bj,$bk) + $lo[1] - $cj) * $vs[1];
    my $dz = (zvals($bi,$bj,$bk) + $lo[2] - $ck) * $vs[2];
    my $d2 = $dx*$dx + $dy*$dy + $dz*$dz;
    return sqrt( ($d2 + (1 - $air)*1e15)->min );
}

# ---------------------------------------------------------------- load volume + affine
my ($vol, $mri2mni, $mni2mri);

if ($o{selftest}) {
    # 合成: 脳(楕円体, 明)+ 骨(その外側の殻, 暗)+ 空気(0)。
    my ($Ni,$Nj,$Nk) = (80, 90, 70);
    my @ctr = (40, 45, 35);
    my $ii = xvals($Ni,$Nj,$Nk);
    my $jj = yvals($Ni,$Nj,$Nk);
    my $kk = zvals($Ni,$Nj,$Nk);
    my $rx = ($ii - $ctr[0]) / 30.0;
    my $ry = ($jj - $ctr[1]) / 34.0;
    my $rz = ($kk - $ctr[2]) / 26.0;
    my $rr = sqrt($rx*$rx + $ry*$ry + $rz*$rz);
    $vol = zeroes(double, $Ni,$Nj,$Nk);
    my $brain = ($rr < 0.86);                       # 脳(明)
    my $skull = (($rr >= 0.86) & ($rr < 1.0));       # 骨(暗)
    my $scalp = (($rr >= 1.0)  & ($rr < 1.08));      # 頭皮(中)
    $vol = $vol + $brain*800 + $skull*150 + $scalp*600;
    # ちょっとした灰白/白質のむら
    $vol += ($brain * (100 * sin($rr*12)));
    # 内部空気ポケット(外耳道アナログ): 頭部内の左右側方に暗い小球を彫る
    for my $cx (14, 66) {
        my $da = sqrt(($ii-$cx)**2 + ($jj-45)**2 + ($kk-32)**2);
        my $pocket = ($da < 4);
        $vol = $vol * (1 - $pocket);   # =0(空気)
    }
    # affine(math 規約): voxel->MNI = 0.5*(v - ctr)
    $mri2mni = zeroes(double,4,4);
    $mri2mni->set(3,3,1);
    for my $r (0..2) { $mri2mni->set($r,$r,0.5); $mri2mni->set($r,3, -0.5*$ctr[$r]); }
    # MNI->voxel = 2*MNI + ctr
    $mni2mri = zeroes(double,4,4);
    $mni2mri->set(3,3,1);
    for my $r (0..2) { $mni2mri->set($r,$r,2.0); $mni2mri->set($r,3, $ctr[$r]); }
    # 既定電極(self-test 用): EEC っぽい左点 + 右の骨上 2点(1つは面外)
    $o{eec}    //= "-15,0,0";
    push @{$o{elec}}, "R_on=15,0,0"   unless @{$o{elec}};
    push @{$o{elec}}, "R_off=15,0,6";
} else {
    print STDERR "reading $o{mat} ...\n";
    require PDL::IO::HDF5;
    my $h5 = PDL::IO::HDF5->new($o{mat}) or die "cannot open $o{mat}\n";
    $vol = h5get($h5, '/sa/mri/data')->double;      # PDL (394,466,378)
    # affine は向きを自動判定(最下行 [0,0,0,1] 不変条件)して math 規約へ
    $mri2mni = to_math_affine( h5get($h5, '/sa/mri2mni')->double );  # voxel->MNI
    $mni2mri = to_math_affine( h5get($h5, '/sa/mni2mri')->double );  # MNI->voxel
    printf STDERR "vol dims = %s\n", join("x", $vol->dims);
}

my ($Ni,$Nj,$Nk) = $vol->dims;

# voxel サイズ(mm)= mri2mni 各列のノルム(面外 mm 換算用)
my @vsize = map {
    my $c = $_;
    sqrt($mri2mni->at(0,$c)**2 + $mri2mni->at(1,$c)**2 + $mri2mni->at(2,$c)**2)
} (0,1,2);

# --------------------------- sanity: MNI(0,0,0) がどの voxel に落ちるか
{
    my @v0 = apply_affine($mni2mri, 0,0,0);
    printf STDERR "MNI(0,0,0) -> voxel (%.1f, %.1f, %.1f)  [dims %dx%dx%d, vsize %.3f/%.3f/%.3f mm]\n",
        @v0, $Ni,$Nj,$Nk, @vsize;
}

# ---------------------------------------------------------------- center
my @center_mni = defined $o{center} ? parse_xyz($o{center})
               : defined $o{eec}    ? parse_xyz($o{eec})
               :                       (0,0,0);
my @cvox = apply_affine($mni2mri, @center_mni);   # (i0,j0,k0) float
my @cidx = map { my $v = int($_ + 0.5); $v } @cvox;
$cidx[0] = 0 if $cidx[0] < 0; $cidx[0] = $Ni-1 if $cidx[0] > $Ni-1;
$cidx[1] = 0 if $cidx[1] < 0; $cidx[1] = $Nj-1 if $cidx[1] > $Nj-1;
$cidx[2] = 0 if $cidx[2] < 0; $cidx[2] = $Nk-1 if $cidx[2] > $Nk-1;
printf STDERR "center MNI (%.1f,%.1f,%.1f) -> voxel (%d,%d,%d)\n", @center_mni, @cidx;

# ---------------------------------------------------------------- electrodes
# 各要素: { label, x,y,z(MNI), type }  type: eec|nyhead|user
my @elecs;
if (defined $o{eec}) {
    my @p = parse_xyz($o{eec});
    push @elecs, { label => "EEC", x=>$p[0], y=>$p[1], z=>$p[2], type=>"eec" };
}
for my $spec (@{$o{elec}}) {
    my ($lab, $rest) = split /=/, $spec, 2;
    die "bad --elec '$spec' (want LABEL=x,y,z)\n" unless defined $rest;
    my @p = parse_xyz($rest);
    push @elecs, { label => $lab, x=>$p[0], y=>$p[1], z=>$p[2], type=>"user" };
}
if (defined $o{nyhead} && !$o{selftest}) {
    my @want = split /\s*,\s*/, $o{nyhead};
    my $ok = eval {
        require PDL::IO::HDF5;
        require PDL::IO::NYHead;
        my $ny  = PDL::IO::NYHead->new($o{mat});
        my $lab = $ny->electrode_labels;            # arrayref of 231 names
        # 位置: 既定 locs_3D_orig(テンプレ標準), --no-elec-orig で登録済 locs_3D(cols0-2)
        my $pos;
        if ($o{elec_orig}) {
            $pos = $ny->electrode_pos;              # (231,3) MNI mm
        } else {
            $pos = h5get(PDL::IO::HDF5->new($o{mat}), '/sa/locs_3D')->slice(':,0:2')->sever;
        }
        my %idx; $idx{ $lab->[$_] } = $_ for 0 .. $#$lab;
        for my $name (@want) {
            unless (exists $idx{$name}) {
                warn "  --nyhead: label '$name' not found, skipping\n";
                next;
            }
            my $r = $idx{$name};
            push @elecs, {
                label => $name,
                x => $pos->at($r,0), y => $pos->at($r,1), z => $pos->at($r,2),
                type => "nyhead",
            };
        }
        1;
    };
    warn "  --nyhead read failed ($@); pass positions via --elec instead\n" unless $ok;
}

# ---------------------------------------------------------------- windows + background
my @bcol = parse_xyz($o{bone_color});
my @acol = parse_xyz($o{air_color});
my ($bw_lo, $bw_hi, $bo_lo, $bo_hi);
my $bg_mode = $o{bg};
my ($vmin, $vmax, $range, $bg_thr, $air_high, $air_thr);
my $headvals;
{
    $vmin  = $vol->min;
    $vmax  = $vol->max;
    $range = ($vmax - $vmin) || 1;

    # 背景の極性判定(min 側 / max 側 どちらに多くのボクセルが溜まっているか)
    if ($bg_mode eq 'auto') {
        my $lo_cnt = ($vol < ($vmin + 0.01*$range))->sum;
        my $hi_cnt = ($vol > ($vmax - 0.01*$range))->sum;
        $bg_mode = ($hi_cnt > $lo_cnt) ? 'high' : 'low';
        printf STDERR "background = %s (auto; lo=%d hi=%d voxels at extremes)\n",
            $bg_mode, $lo_cnt, $hi_cnt;
    }
    $bg_thr = ($bg_mode eq 'high') ? $vmax - 0.03*$range
            : ($bg_mode eq 'low')  ? $vmin + 0.03*$range : undef;

    # 空気の輝度側(既定=背景と同じ極性)。ear canal 等の内部空気の色付けに使う。
    $air_high = ($o{air_side} eq 'high') ? 1
              : ($o{air_side} eq 'low')  ? 0
              : ($bg_mode eq 'high' ? 1 : 0);
    $air_thr  = $air_high ? $vmax - 0.03*$range : $vmin + 0.03*$range;

    # 背景を除いた「頭部ボクセル」で窓を決める
    my $mask;
    if    ($bg_mode eq 'high') { $mask = ($vol < ($vmax - 0.05*$range)); }
    elsif ($bg_mode eq 'low')  { $mask = ($vol > ($vmin + 0.05*$range)); }
    else                       { $mask = ones(byte, $vol->dims); }
    my $vals = $vol->flat->where($mask->flat);
    $vals = $vol->flat if $vals->nelem < 100;
    $headvals = $vals;

    printf STDERR "intensity pct(head): p1=%.0f p5=%.0f p25=%.0f p50=%.0f p75=%.0f p95=%.0f p99=%.0f  (min=%.0f max=%.0f)\n",
        (map { pctl($vals,$_) } (0.01,0.05,0.25,0.50,0.75,0.95,0.99)), $vmin, $vmax;

    if (defined $o{brain_window}) { ($bw_lo,$bw_hi) = map { 0+$_ } split /\s*,\s*/, $o{brain_window}; }
    else { $bw_lo = pctl($vals, 0.10); $bw_hi = pctl($vals, 0.98); }
    if (defined $o{bone_window})  { ($bo_lo,$bo_hi) = map { 0+$_ } split /\s*,\s*/, $o{bone_window}; }
    else { $bo_lo = pctl($vals, 0.02); $bo_hi = pctl($vals, 0.40); }
    printf STDERR "brain window = [%.1f, %.1f]   overlay=%s   invert=%d\n",
        $bw_lo,$bw_hi,$o{overlay},$o{invert};
}

# ---------------------------------------------------------------- electrode intensity readout
# クリック点が素の MRI 上でどんな輝度(=空気/骨/軟部/脳)に落ちているかを数値で出す。
{
    print STDERR "electrode intensities on raw /sa/mri/data (cube = +-3 voxel neighborhood):\n";
    for my $e (@elecs) {
        my @v  = apply_affine($mni2mri, $e->{x}, $e->{y}, $e->{z});
        my @iv = map { int($_ + 0.5) } @v;
        my @N  = ($Ni,$Nj,$Nk);
        my @c  = map { my $x=$iv[$_]; $x=0 if $x<0; $x=$N[$_]-1 if $x>$N[$_]-1; $x } (0,1,2);
        my $val = $vol->at(@c);
        my @lo = map { my $x=$c[$_]-3; $x<0?0:$x } (0,1,2);
        my @hi = map { my $x=$c[$_]+3; $x>$N[$_]-1?$N[$_]-1:$x } (0,1,2);
        my $cube = $vol->slice("$lo[0]:$hi[0],$lo[1]:$hi[1],$lo[2]:$hi[2]");
        my $rank = 100 * ($headvals < $val)->sum / $headvals->nelem;
        my $is_air = $air_high ? ($val > $air_thr) : ($val < $air_thr);
        my $tag = $is_air        ? "AIR"
                : $rank > 85     ? "bright tissue (scalp/skin/fat)"
                : $rank < 15     ? "dark tissue (brain/CSF)"
                :                  "mid tissue";
        my $dair = nearest_air_mm(@c, $air_high, $air_thr, $Ni,$Nj,$Nk, \@vsize, $vol);
        printf STDERR "  %-6s vox(%d,%d,%d) I=%.0f  rank=%.0f%% head  -> %s;  nearest air = %s\n",
            $e->{label}, @c, $val, $rank, $tag,
            (defined $dair ? sprintf("%.1f mm", $dair) : ">35mm");
    }
}

# ---------------------------------------------------------------- plane specs
# 各 plane: fixed 軸, 面内 (col 軸,方向) (row 軸,方向), コーナーラベル。
# neurological 既定(患者右=画像右)。radiological は L-R を反転。
my $lr_flip = $o{radiological} ? 1 : 0;
my %PLANE = (
    axial => {                       # k 固定: col=i(x), row=j(y) anterior up
        fixed => 2, cidx => $cidx[2],
        col_axis => 0, col_flip => $lr_flip,
        row_axis => 1, row_flip => 1,          # +j=anterior を上へ
        title => sub { sprintf("Axial  (z = %.0f mm)", $center_mni[2]) },
        labels => $o{radiological}
            ? { top=>"A", bottom=>"P", left=>"R", right=>"L" }
            : { top=>"A", bottom=>"P", left=>"L", right=>"R" },
    },
    coronal => {                     # j 固定: col=i(x), row=k(z) superior up
        fixed => 1, cidx => $cidx[1],
        col_axis => 0, col_flip => $lr_flip,
        row_axis => 2, row_flip => 1,
        title => sub { sprintf("Coronal  (y = %.0f mm)", $center_mni[1]) },
        labels => $o{radiological}
            ? { top=>"S", bottom=>"I", left=>"R", right=>"L" }
            : { top=>"S", bottom=>"I", left=>"L", right=>"R" },
    },
    sagittal => {                    # i 固定: col=j(y) anterior left, row=k(z) sup up
        fixed => 0, cidx => $cidx[0],
        col_axis => 1, col_flip => 1,          # anterior(+j) を左へ
        row_axis => 2, row_flip => 1,
        title => sub { sprintf("Sagittal  (x = %.0f mm)", $center_mni[0]) },
        labels => { top=>"S", bottom=>"I", left=>"A", right=>"P" },
    },
);
my @order = qw(axial coronal sagittal);

# plane の 2D スライス(full-res, [row,col])を取り出す
sub plane_slice {
    my ($spec) = @_;
    my @idx = (undef,undef,undef);
    $idx[$spec->{fixed}] = $spec->{cidx};
    # slice string
    my @ss = map { defined $idx[$_] ? "($idx[$_])" : ":" } (0,1,2);
    my $s2 = $vol->slice(join(",",@ss))->squeeze;   # 面内2軸(元の軸順)
    # s2 の dim0,dim1 = 面内で残った小さい方/大きい方の軸番号
    my @rem = grep { $_ != $spec->{fixed} } (0,1,2);   # 昇順
    # s2 dims 対応: dim0=rem[0], dim1=rem[1]
    # col 軸 / row 軸 を s2 の dim に割り当てて (row,col) 並びに整える
    my ($cd) = grep { $rem[$_] == $spec->{col_axis} } (0,1);
    my ($rd) = grep { $rem[$_] == $spec->{row_axis} } (0,1);
    # 目標 (row,col) = (rd, cd) の順に転置
    my $img = ($rd == 0) ? $s2 : $s2->transpose;   # -> dim0=row 軸, dim1=col 軸
    $img = $img->sever;
    $img = $img->slice("-1:0,:")->sever if $spec->{row_flip};
    $img = $img->slice(":,-1:0")->sever if $spec->{col_flip};
    return $img;   # [Hrow, Wcol]
}

# voxel(i,j,k) -> この plane の full-res (row,col,off_vox)
sub map_vox {
    my ($spec, $i,$j,$k) = @_;
    my @v = ($i,$j,$k);
    my @N = ($Ni,$Nj,$Nk);
    my $ca = $spec->{col_axis};
    my $ra = $spec->{row_axis};
    my $col = $v[$ca]; $col = ($N[$ca]-1) - $col if $spec->{col_flip};
    my $row = $v[$ra]; $row = ($N[$ra]-1) - $row if $spec->{row_flip};
    my $off = $v[$spec->{fixed}] - $spec->{cidx};   # voxel
    return ($row, $col, $off);
}

# ---------------------------------------------------------------- render
require PDL::Graphics::Cairo;
PDL::Graphics::Cairo->import(qw(subplots));

my ($fig, @ax) = subplots(1, 3, figsize => [16.5, 6.2]);

for my $pi (0 .. $#order) {
    my $spec = $PLANE{ $order[$pi] };
    my $ax   = $ax[$pi];
    my $img  = plane_slice($spec);                  # [H,W] full-res
    my ($Hf, $Wf) = $img->dims;

    # subsample(軽量化)。step は整数ストライド。
    my $step = 1;
    my $md = $o{max_dim};
    $step = int( (($Hf > $Wf ? $Hf : $Wf) + $md - 1) / $md );
    $step = 1 if $step < 1;
    my $S = ($step > 1) ? $img->slice("0:-1:$step,0:-1:$step")->sever : $img;
    my ($Hs, $Ws) = $S->dims;

    # --- 合成: 組織 grayscale + オーバレイ(air=内部空気 / bone=輝度帯 / none)
    my $gn = (($S - $bw_lo) / (($bw_hi - $bw_lo) || 1))->clip(0,1);
    $gn = 1 - $gn if $o{invert};

    # 外部背景(padding)の連結成分だけ黒に落とす → 内部空洞(ear canal 等)は残す
    my $ext = zeroes(byte, $Hs, $Ws);
    if ($bg_mode ne 'none') { ($ext, undef) = air_components($S, $bg_mode eq 'high', $bg_thr); }
    my $g = $gn * (1 - $ext);

    my ($a, @ocol);
    if ($o{overlay} eq 'bone') {
        my $edge = (($bo_hi - $bo_lo) * 0.15) || 1e-6;
        my $up   = (($S - $bo_lo) / $edge)->clip(0,1);
        my $dn   = (($bo_hi - $S) / $edge)->clip(0,1);
        $a    = ($up * $dn)->clip(0,1) * $o{bone_opacity} * (1 - $ext);
        @ocol = @bcol;
    } elsif ($o{overlay} eq 'air') {
        my (undef, $intair) = air_components($S, $air_high, $air_thr);
        $a    = $intair * $o{bone_opacity};
        @ocol = @acol;
    } else {
        $a    = zeroes(float, $Hs, $Ws);
        @ocol = (0,0,0);
    }

    my $rgb  = zeroes(float, $Hs, $Ws, 3);
    $rgb->slice(":,:,(0)") .= $g*(1-$a) + $ocol[0]*$a;
    $rgb->slice(":,:,(1)") .= $g*(1-$a) + $ocol[1]*$a;
    $rgb->slice(":,:,(2)") .= $g*(1-$a) + $ocol[2]*$a;

    $ax->imshow($rgb, origin => 'upper');
    $ax->axis('off');
    $ax->set_aspect('equal');
    $ax->set_title($spec->{title}->());

    # 表示座標: imshow は xmax=Ws, ymax=Hs, origin upper。
    # マーカーは (x=col_s, y=Hs-row_s) で画像と一致する(origin upper 補正)。
    my $to_xy = sub { my ($row,$col)=@_; ($col/$step, $Hs - $row/$step) };

    my %col_of = (eec=>[1,0,0], nyhead=>[0,0.8,1], user=>[1,0.85,0]);

    # 電極の表示情報を先に集める(描画は「マーカー全部→テキスト全部」の順で行う。
    # P:G:C の 'o' マーカーは arc 前に new_path しないため、直前の text が残した
    # current point から arc 始点へ迷い線が引かれる。マーカーを連続描画(各 stroke
    # が path をクリア)し、text を後回しにすることでこれを避ける)。
    my @draw;
    for my $e (@elecs) {
        my @ev = apply_affine($mni2mri, $e->{x}, $e->{y}, $e->{z});
        my ($row,$col,$offv) = map_vox($spec, @ev);
        next if $col < -2 || $col > $Wf+2 || $row < -2 || $row > $Hf+2;
        my ($ex,$ey) = $to_xy->($row,$col);
        my $off_mm = abs($offv) * $vsize[ $spec->{fixed} ];
        my $c = $col_of{ $e->{type} } // [1,0.85,0];
        my $on = ($off_mm <= $o{slab});
        my $tag = $on ? $e->{label}
                      : sprintf("%s %+.0fmm", $e->{label}, $offv>0 ? $off_mm : -$off_mm);
        push @draw, { ex=>$ex, ey=>$ey, c=>$c, on=>$on, tag=>$tag };
    }

    # (1) center crosshair — 最後の描画は axhline の stroke なので path はクリア済
    my ($crow,$ccol,undef) = map_vox($spec, @cvox);
    my ($cx,$cy) = $to_xy->($crow,$ccol);
    $ax->axvline($cx, color=>[0.2,1,0.4], lw=>0.8, ls=>'dashed');
    $ax->axhline($cy, color=>[0.2,1,0.4], lw=>0.8, ls=>'dashed');

    # (2) 電極マーカー(text を挟まず連続で)
    for my $d (@draw) {
        if ($d->{on}) {
            $ax->scatter(pdl($d->{ex}), pdl($d->{ey}), s=>7, color=>$d->{c},
                         marker=>'o', alpha=>0.95, fillstyle=>'full');
        } else {
            $ax->scatter(pdl($d->{ex}), pdl($d->{ey}), s=>7, color=>$d->{c},
                         marker=>'o', alpha=>0.45, fillstyle=>'none');
        }
    }

    # (3) テキスト(電極ラベル + コーナーの方位ラベル)= すべてマーカーの後
    for my $d (@draw) {
        $ax->text($d->{ex} + $Ws*0.015, $d->{ey}, $d->{tag},
                  color=>$d->{c}, fontsize=>9, ha=>'left', va=>'middle');
    }
    my $L = $spec->{labels};
    $ax->text($Ws*0.5, $Hs*0.98, $L->{top},    color=>'white', fontsize=>13, ha=>'center', va=>'top');
    $ax->text($Ws*0.5, $Hs*0.02, $L->{bottom}, color=>'white', fontsize=>13, ha=>'center', va=>'bottom');
    $ax->text($Ws*0.02,$Hs*0.5,  $L->{left},   color=>'white', fontsize=>13, ha=>'left',   va=>'middle');
    $ax->text($Ws*0.98,$Hs*0.5,  $L->{right},  color=>'white', fontsize=>13, ha=>'right',  va=>'middle');
}

$fig->tight_layout;
$fig->save($o{out});
print STDERR "wrote $o{out}\n";

if ($o{show}) {
    eval { $fig->show(backend=>'osx'); 1 } or warn "show failed: $@";
}
