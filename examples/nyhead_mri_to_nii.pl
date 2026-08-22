#!/usr/bin/env perl
# nyhead_mri_to_nii.pl — small prototype
#
# NYHead sa_nyhead.mat の /sa/mri/data (378x466x394 f64) を NIfTI-1 (.nii) に書き出す。
# 目的は「NYHead -> .nii -> PDL::IO::BIDS で読む」道の真ん中(NIfTI writer)を通すこと。
#
# 設計上の要点:
#  - affine(sform) を完全制御したいので PDL::IO::Nifti->write_nii は使わず、自前の
#    最小 NIfTI-1 ライタ write_nii1() を持つ。sform に NYHead の voxel->MNI 変換を入れる。
#  - ボリューム本体は get_dataref + (host が big-endian の時のみ) bswap + syswrite の
#    直流し込み。piddle -> Perl list -> pack は通らない(6900万要素のリストを作らない)。
#  - HDF5 読取は PDL::IO::HDF5 が bulk で piddle に入れる(リストを経ない)ので既に速い。
#
# 実データの校正(--mat 経路)で Mac 上で確認すべき2点はコード内の CALIBRATE コメント参照。

use strict; use warnings;
use PDL;
use PDL::IO::Misc;      # bswap*
use Getopt::Long;

my ($mat, $out, $xform, $self_test, $reorder, $dtype);
$xform = 'mri2mni';     # sform に使う NYHead 変換名(voxel index -> MNI mm 想定)
$out   = 'nyhead_mri.nii';
$dtype = 'keep';        # keep|float32 : NYHead MRI は f64(過剰)。float32 推奨(半分/警告消/視認上ロス無)
GetOptions(
  'mat=s'     => \$mat,
  'out|o=s'   => \$out,
  'xform=s'   => \$xform,
  'reorder!'  => \$reorder,   # PDL(394,466,378) を (378,466,394) に軸反転して書く
  'dtype=s'   => \$dtype,     # keep(既定) | float32
  'self-test' => \$self_test,
) or die "bad args\n";

# ---- NIfTI-1 データ型コード ----
my %NIFTI_DT = (   # PDL type name => [datatype code, bitpix]
  'byte'   => [2,   8],   # DT_UINT8
  'short'  => [4,  16],   # DT_INT16
  'long'   => [8,  32],   # DT_INT32
  'float'  => [16, 32],   # DT_FLOAT32
  'double' => [64, 64],   # DT_FLOAT64
);

# ---- 最小 NIfTI-1 ライタ ----
# $vol    : 3-D piddle (i,j,k) 順(=NIfTI voxel 順)
# $affine : 4x4 piddle, voxel index [i,j,k,1] -> world [x,y,z,1] mm (RAS)
sub write_nii1 {
  my ($vol, $affine, $path) = @_;
  my ($ni,$nj,$nk) = $vol->dims;  $nk //= 1;
  my $tname = "".$vol->type;
  my $dt = $NIFTI_DT{$tname} or die "unsupported PDL type '$tname' for NIfTI\n";
  my ($datatype,$bitpix) = @$dt;

  # affine の各列ノルム = voxel spacing(pixdim)
  my @a = list $affine;             # row-major 16 要素
  my @A = ([@a[0..3]],[@a[4..7]],[@a[8..11]],[@a[12..15]]);
  my $col = sub { my $c=shift; sqrt($A[0][$c]**2+$A[1][$c]**2+$A[2][$c]**2) };
  my ($dx,$dy,$dz) = ($col->(0)||1, $col->(1)||1, $col->(2)||1);

  # ---- 348 byte ヘッダを LE で組む ----
  my $h = "\0" x 348;
  substr($h,  0,4) = pack('V', 348);                 # sizeof_hdr
  substr($h, 40,2) = pack('v', 3);                   # dim[0] = ndim
  substr($h, 42,2) = pack('v', $ni);                 # dim[1]
  substr($h, 44,2) = pack('v', $nj);                 # dim[2]
  substr($h, 46,2) = pack('v', $nk);                 # dim[3]
  substr($h, 48,10)= pack('v5', 1,1,1,1,1);          # dim[4..7]=1
  substr($h, 70,2) = pack('v', $datatype);
  substr($h, 72,2) = pack('v', $bitpix);
  # pixdim[0]=qfac(未使用,1), pixdim[1..3]=spacing
  substr($h, 76,4) = pack('f', 1.0);
  substr($h, 80,4) = pack('f', $dx);
  substr($h, 84,4) = pack('f', $dy);
  substr($h, 88,4) = pack('f', $dz);
  substr($h,108,4) = pack('f', 352.0);               # vox_offset
  substr($h,112,4) = pack('f', 1.0);                 # scl_slope
  substr($h,116,4) = pack('f', 0.0);                 # scl_inter
  substr($h,252,2) = pack('v', 0);                   # qform_code = 0 (sform のみ使用)
  substr($h,254,2) = pack('v', 4);                   # sform_code = 4 (MNI152)
  # srow_x/y/z = affine の上3行
  substr($h,280,16) = pack('f4', @{$A[0]});
  substr($h,296,16) = pack('f4', @{$A[1]});
  substr($h,312,16) = pack('f4', @{$A[2]});
  substr($h,344,4)  = "n+1\0";                        # magic (single-file .nii)

  # ---- 本体 = get_dataref 直流し込み ----
  my $v = $vol->copy->sever;         # contiguous 保証
  my $host_be = unpack('C', pack('S',1)) == 0;   # host big-endian?
  if ($host_be) {                    # ファイルは LE 宣言なので BE host なら swap
    my %sw = (short=>'bswap2',long=>'bswap4',float=>'bswap4',double=>'bswap8');
    if (my $m = $sw{$tname}) { $v->$m; }
  }
  my $ref = $v->get_dataref;         # \ raw bytes, list を経ない
  my $expect = $ni*$nj*$nk*($bitpix/8);
  length($$ref)==$expect or die sprintf "payload %d != expected %d\n",length($$ref),$expect;

  open my $fh,'>:raw',$path or die "open $path: $!";
  print $fh $h;
  print $fh pack('V',0);             # 4-byte extender flag = 0
  print $fh $$ref;                   # payload
  close $fh;
  return { dims=>[$ni,$nj,$nk], datatype=>$datatype, bytes=>352+$expect };
}

# ---- 実データ経路 ----
sub load_from_mat {
  my ($mat,$xform) = @_;
  require PDL::IO::HDF5;
  my $h5  = PDL::IO::HDF5->new($mat);
  my $vol = $h5->dataset('/sa/mri/data')->get;    # HDF5 {378,466,394} -> PDL (394,466,378)
  my $T   = $h5->dataset("/sa/$xform")->get;      # 4x4
  # FIX(BIDS02): HDF5/MATLAB は column-major。PDL は数学行列の転置で読む。
  # write_nii1 は「行=各世界軸(srow)」の数学規約前提なので、ここで転置して戻す。
  # 検証: mri2mni の並進(-98.5,-134.5,-72.5)は as-read では最終行、transpose で
  # 最終列(NIfTI srow の並進位置)に来る。転置を怠ると並進が落ちて diag(0.5) だけ残る。
  $T = $T->transpose->sever;
  # CALIBRATE-1: 軸順。HDF5 逆順で PDL は (394,466,378)。affine が MATLAB voxel 順
  #   (378,466,394)=(i,j,k) 前提なら --reorder で軸反転してから書く。
  if ($reorder) { $vol = $vol->xchg(0,2)->sever; }  # (394,466,378)->(378,466,394)
  # CALIBRATE-2: affine の転置と index base。NYHead 4x4 が PDL で転置して入る/1-based
  #   voxel の可能性(BIDS01 で指摘の転置トラップ)。実データで既知ランドマークが正しい
  #   MNI mm に落ちるかを nibabel/BIDS voxel_to_world で1度だけ確認して確定する。
  return ($vol, $T);
}

# ---- self-test 用の合成データ(非対称 dims + 回転/異方 spacing/offset) ----
sub synth {
  my $vol = (sequence(double, 5,7,9) % 13);        # 非対称 = 軸取り違えを検出
  # 異方 spacing(2,3,4) + 90度回転混ぜ + offset。voxel[i,j,k]->mm
  my $aff = pdl([[ 0,  3, 0, -10],
                 [ 2,  0, 0,  20],
                 [ 0,  0, 4,  -5],
                 [ 0,  0, 0,   1]])->double;
  return ($vol, $aff);
}

# ================= main =================
my ($vol, $aff);
if ($self_test) {
  ($vol,$aff) = synth();
} elsif ($mat) {
  ($vol,$aff) = load_from_mat($mat,$xform);
} else {
  die "usage: $0 --self-test | --mat sa_nyhead.mat [--xform mri2mni] [--reorder] -o out.nii\n";
}

if ($dtype eq 'float32') { $vol = $vol->float; }
elsif ($dtype ne 'keep') { die "unknown --dtype '$dtype' (keep|float32)\n"; }

printf STDERR "vol dims = (%s)  type=%s  sform xform=%s\n",
  join(",",$vol->dims), $vol->type, ($mat?$xform:'synthetic');
my $info = write_nii1($vol, $aff, $out);
printf STDERR "wrote %s  (%d bytes, datatype=%d)\n", $out, $info->{bytes}, $info->{datatype};
