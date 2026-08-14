package PDL::IO::NYHead;
# ============================================================================
# PDL::IO::NYHead  —  New York Head (sa_nyhead.mat, MATLAB v7.3=HDF5) を
#                     PDL::EEG に持ち込むためのトップレベル reader。
#
# 数値データセットは PDL::IO::HDF5 で、MATLAB cell-of-strings(電極名・HO野名)は
# PDL::IO::NYHead::H5Dump(h5dump CLI 経由)で読む。両者を1つの OO API に統合する。
#
# ---- 軸規約(実測確定) ----------------------------------------------------
#   HDF5/h5dump のシェイプ {A,B,...} は PDL では逆順 (...,B,A) で読まれる。
#   NYHead は座標系を {3,N}(=MATLAB 列優先の (N,3))で持つため、PDL では
#   vc/normals/locs/tri がそのまま (N,3) / (Nf,3) になる(転置不要)。
#     /sa/cortex75K/vc  {3,74382}         -> PDL (74382,3)      各行=頂点xyz
#     /sa/cortex75K/tri {3,148756}        -> PDL (148756,3)     1-based
#     /sa/cortex75K/V_fem_normal {74382,231} -> PDL (231,74382) (電極,源) CAR
#     /sa/cortex75K/V_fem {3,74382,231}   -> PDL (231,74382,3)  (電極,源,xyz)
#     /sa/cortex75K/in_HO {1,74382}       -> PDL (74382,1)      flat で (74382)
#     /sa/locs_3D_orig {3,231}            -> PDL (231,3)
#
# ---- PDL 最適化 -----------------------------------------------------------
#   * 低解像 cortex 頂点の復元は Perl ループでなく fancy index(dice_axis)。
#   * tri / in_from の 1-based→0-based は PDL 演算($x-1)で一括。
#   * nearest_vertex は電極ごとに内ループのみベクトル化(巨大な全対全配列を
#     作らない):素朴フル・ブロードキャストより速くメモリも軽い。
#
# 由来: goosh の NYHead 取り込み作業(P:EEG14〜, H5Dump は下位モジュール)。
# ============================================================================
use strict;
use warnings;
use Carp;
use PDL;
use PDL::IO::HDF5;
use PDL::IO::NYHead::H5Dump;

our $VERSION = '0.01';

# 19ch 対応(10-20)の 231ch モデル内 0-based index。P:EEG14 で sed 相互検証し凍結。
# 並びは Fp1 Fp2 F3 F4 C3 C4 P3 P4 O1 O2 F7 F8 T7 T8 P7 P8 Fz Cz Pz。
our @IDX19 = (156,157,119,115,159,116,120,222,19,48,20,223,182,178,220,179,183,168,169);

# new($matfile, h5dump => 'h5dump')
sub new {
    my ($class, $file, %arg) = @_;
    defined $file && -r $file or croak "PDL::IO::NYHead: cannot read '$file'";
    my $self = {
        file   => $file,
        h5     => PDL::IO::HDF5->new($file),
        h5dump => PDL::IO::NYHead::H5Dump->new(h5dump => $arg{h5dump} || 'h5dump'),
        cache  => {},
    };
    return bless $self, $class;
}

# ---- 低レベル: 数値データセット(キャッシュ付き) ----------------------------
sub _get {
    my ($self, $path) = @_;
    return $self->{cache}{$path} //= $self->{h5}->dataset($path)->get;
}

# ---- 電極 -----------------------------------------------------------------
sub electrode_labels {                       # 231 個(cell-of-strings)
    my $self = shift;
    return $self->{cache}{clab} //=
        $self->{h5dump}->read_cell_strings($self->{file}, '/sa/clab_electrodes');
}
sub electrode_pos { $_[0]->_get('/sa/locs_3D_orig') }   # (231,3) mm, MNI
sub n_electrodes  { $_[0]->electrode_pos->dim(0) }

sub idx19 { return pdl(long, \@IDX19) }
sub electrode_pos_19 {                        # (19,3)
    my $self = shift;
    return $self->electrode_pos->dice_axis(0, $self->idx19);
}
sub electrode_labels_19 {                     # 19 個
    my $self = shift;
    my $all = $self->electrode_labels;
    return [ @{$all}[@IDX19] ];
}

# ---- リードフィールド ------------------------------------------------------
sub leadfield { $_[0]->_get('/sa/cortex75K/V_fem_normal') }   # (231,74382) 法線射影 CAR
sub leadfield_full { $_[0]->_get('/sa/cortex75K/V_fem') }     # (231,74382,3) 3成分(巨大)
sub leadfield_19 {                                            # (19,74382)
    my $self = shift;
    return $self->leadfield->dice_axis(0, $self->idx19);
}

# ---- サーフェス(head / cortex75K) -----------------------------------------
#   {vc=>(Nv,3), tri=>(Nf,3) 1-based(MATLAB/OBJ そのまま), tri0=>(Nf,3) 0-based}
sub surface {
    my ($self, $name) = @_;
    $name //= 'head';
    my $g = "/sa/$name";
    my $vc  = $self->_get("$g/vc");
    my $tri = $self->_get("$g/tri");
    return { vc => $vc, tri => $tri, tri0 => long($tri) - 1 };
}

# ---- 低解像 cortex(cortex1K/2K/5K/10K)を 75K 頂点から fancy index 復元 -----
#   {vc=>(M,3), tri=>(Mf,3) 1-based(低解像頂点内), tri0=>0-based,
#    in_from=>(M) 0-based の 75K 頂点 index}
sub cortex {
    my ($self, $res) = @_;
    $res //= 'cortex75K';
    if ($res eq 'cortex75K') {
        my $s = $self->surface('cortex75K');
        $s->{in_from} = sequence(long, $s->{vc}->dim(0));
        return $s;
    }
    my $in_from = long($self->_get("/sa/$res/in_from_cortex75K")->flat) - 1;  # 0-based
    my $vc75    = $self->_get('/sa/cortex75K/vc');
    my $vc      = $vc75->dice_axis(0, $in_from);          # (M,3) fancy index
    my $tri     = $self->_get("/sa/$res/tri");            # (Mf,3) 1-based(低解像内)
    return { vc => $vc, tri => $tri, tri0 => long($tri) - 1, in_from => $in_from };
}

# ---- Harvard-Oxford atlas -------------------------------------------------
sub ho_labels {                              # 97 個
    my $self = shift;
    return $self->{cache}{ho} //=
        $self->{h5dump}->read_cell_strings($self->{file}, '/sa/HO_labels');
}
sub atlas_index {                            # (74382) 1-based(0=none)
    my $self = shift;
    return long($self->_get('/sa/cortex75K/in_HO')->flat);
}
sub area_of_vertex {                         # 75K 頂点 index(0-based)→野名 or undef
    my ($self, $v) = @_;
    my $ho = $self->atlas_index->at($v);     # 1-based
    return undef if $ho <= 0;
    return $self->ho_labels->[$ho - 1];
}

# ---- 座標変換(4x4 affine) --------------------------------------------------
#   注意: PDL は {4,4} を転置で読む。MATLAB 側 M に対し PDL は M' を返すので、
#   列ベクトル p(4) に MATLAB の M を掛けるには (p x $T) を使う(下記 apply_xform)。
sub transform { my ($self, $name) = @_; $self->_get("/sa/$name") }   # (4,4) = MATLAB の転置
sub apply_xform {                            # $pts (N,3) を 4x4 affine で変換 -> (N,3)
    my ($self, $pts, $name) = @_;
    my $T = $self->transform($name);                 # PDL(4,4) = MATLAB M の転置
    my $N = $pts->dim(0);
    my $hom = ones(4, $N);
    $hom->slice("0:2,:") .= $pts->transpose;         # (4,N): 上3行=xyz, 4行目=1
    my $out = $T x $hom;                             # MATLAB: M*[x;y;z;1]
    return $out->slice("0:2,:")->transpose->sever;   # (N,3)
}

# ---- 最近傍頂点(ベクトル化。全対全の巨大配列を作らない) --------------------
#   $pts: (P,3) または (3)。$which: サーフェス名(既定 cortex75K)。
#   戻り: 0-based 頂点 index の (P) piddle。
sub nearest_vertex {
    my ($self, $pts, $which) = @_;
    $which //= 'cortex75K';
    my $vT = $self->_get("/sa/$which/vc")->transpose;  # (3,Nvert) 座標優先
    $pts = $pts->dummy(0) if $pts->ndims == 1;         # (3)->(1,3)
    my $P = $pts->dim(0);
    my $out = zeroes(long, $P);
    for my $i (0 .. $P-1) {
        my $d2 = (($vT - $pts->slice("($i),"))**2)->sumover;  # (Nvert)
        $out->set($i, $d2->minimum_ind);
    }
    return $P == 1 ? $out->at(0) : $out;
}

# ---- メッシュ取得(head / cortex75K / 低解像 cortexNK を一元化) -------------
#   { vc=>(Nv,3), tri=>(Nf,3) 1-based }。低解像は cortex() で 75K から復元。
sub mesh {
    my ($self, $name) = @_;
    $name //= 'head';
    my $g = ($name eq 'head' || $name eq 'cortex75K')
          ? $self->surface($name) : $self->cortex($name);
    return { vc => $g->{vc}, tri => $g->{tri} };
}

# /sa 直下に実在する cortex 解像度グループ(cortex1K.. など)を昇順で返す
sub cortex_resolutions {
    my $self = shift;
    my @g = eval { $self->{h5}->group('/sa')->groups };
    return sort grep { /^cortex\d+K$/ } @g;
}

# ---- .obj 書き出し(OBJ は 1-based。tri をそのまま使う) --------------------
sub to_obj {
    my ($self, $name, $out) = @_;
    my $m  = $self->mesh($name);
    my @v  = $m->{vc}->xchg(0,1)->flat->list;    # [x0,y0,z0, x1,...] 行優先
    my @t  = $m->{tri}->xchg(0,1)->flat->list;   # [a0,b0,c0, ...] 1-based
    my $fh = _open_out($out);
    printf $fh "# New York Head '%s' mesh (via PDL::IO::NYHead)\n", $name;
    printf $fh "v %.6f %.6f %.6f\n", @v[3*$_, 3*$_+1, 3*$_+2] for 0 .. @v/3 - 1;
    printf $fh "f %d %d %d\n",       @t[3*$_, 3*$_+1, 3*$_+2] for 0 .. @t/3 - 1;
    close $fh;
    return $out;
}

# ---- .usda 書き出し(USD は 0-based。faceVertexIndices は tri-1) -----------
sub to_usda {
    my ($self, $name, $out, %opt) = @_;
    my $axes = $opt{axes} // 0;   # >0 で原点から長さ axes(mm) の xyz 軸(赤=X/緑=Y/青=Z)を併記。0=軸なし(既定)
    my $m   = $self->mesh($name);
    my $nv  = $m->{vc}->dim(0);
    my $nf  = $m->{tri}->dim(0);
    my @v   = $m->{vc}->xchg(0,1)->flat->list;
    my @idx = (long($m->{tri}) - 1)->xchg(0,1)->flat->list;   # 0-based
    my $fh  = _open_out($out);
    # MNI/NY Head は Z-up・mm。USD 既定(Y-up, m)のままだと頭が横倒し・軸が変になるので明示。
    print  $fh "#usda 1.0\n(\n    upAxis = \"Z\"\n    metersPerUnit = 0.001\n)\n\ndef Xform \"NYHead\"\n{\n    def Mesh \"$name\"\n    {\n";
    print  $fh "        int[] faceVertexCounts = [", join(",", (3) x $nf), "]\n";
    print  $fh "        int[] faceVertexIndices = [", join(",", @idx), "]\n";
    print  $fh "        point3f[] points = [",
               join(",", map { sprintf "(%.6f,%.6f,%.6f)", @v[3*$_,3*$_+1,3*$_+2] } 0..$nv-1),
               "]\n";
    print  $fh "    }\n";
    if ($axes > 0) {
        my ($L, $W) = ($axes, $axes / 50);   # MNI mm フレーム。head 半径相当は ~90
        print $fh <<"USDA";
    def BasisCurves "Axes"
    {
        uniform token type = "linear"
        int[] curveVertexCounts = [2, 2, 2]
        point3f[] points = [(0,0,0),($L,0,0), (0,0,0),(0,$L,0), (0,0,0),(0,0,$L)]
        color3f[] primvars:displayColor = [(1,0,0),(0,1,0),(0,0,1)] (interpolation = "uniform")
        float[] widths = [$W] (interpolation = "constant")
    }
USDA
    }
    print  $fh "}\n";
    close $fh;
    return $out;
}

# フォーマット分岐(スクリプトの --format はこれを呼ぶ)
sub write_mesh {
    my ($self, $name, $format, $out, %opt) = @_;
    $format = lc $format;
    return $self->to_obj($name, $out)        if $format eq 'obj';
    return $self->to_usda($name, $out, %opt) if $format eq 'usda';
    croak "write_mesh: unknown format '$format' (obj|usda)";
}

sub _open_out {
    my $out = shift;
    return $out if ref $out;                       # 既に filehandle
    open(my $fh, '>', $out) or croak "cannot write '$out': $!";
    return $fh;
}

1;

__END__

=encoding utf-8

=head1 NAME

PDL::IO::NYHead - read the New York Head (sa_nyhead.mat) into PDL for PDL::EEG

=head1 SYNOPSIS

    use PDL::IO::NYHead;
    my $ny = PDL::IO::NYHead->new('sa_nyhead.mat');

    my $labels = $ny->electrode_labels;       # 231 names
    my $pos    = $ny->electrode_pos;          # (231,3) MNI mm
    my $L      = $ny->leadfield;              # (231,74382) normal-oriented CAR

    my $ctx    = $ny->surface('cortex75K');   # {vc,tri,tri0}
    my $lo     = $ny->cortex('cortex2K');     # reconstructed via in_from_cortex75K

    my $areas  = $ny->ho_labels;              # 97 Harvard-Oxford area names
    print $ny->area_of_vertex(1234);          # e.g. 'Left Frontal Pole'

    my $iv = $ny->nearest_vertex($pos, 'head'); # each electrode -> nearest scalp vertex

=head1 DESCRIPTION

Loads numeric datasets via L<PDL::IO::HDF5> and MATLAB cell-of-string fields
(electrode / atlas labels) via L<PDL::IO::NYHead::H5Dump>. Axis order is the
reverse of the h5dump shape; helpers return the natural PDL orientation.

=head1 REQUIREMENTS

L<PDL::IO::HDF5>, and the C<h5dump> executable for label fields.

=head1 SEE ALSO

L<PDL::IO::NYHead::H5Dump>, L<PDL::IO::HDF5>, L<PDL::EEG>

=cut
