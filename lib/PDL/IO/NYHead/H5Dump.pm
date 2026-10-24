package PDL::IO::NYHead::H5Dump;
# ============================================================================
# PDL::IO::NYHead::H5Dump — MATLAB v7.3(HDF5)の cell-of-strings を読む。
#
# cell-of-strings は「オブジェクト参照の配列」で、各参照の先が uint16 の文字コード列。
# 読み方は 2 通り:
#   h5dump   … h5dump(HDF5 tools)を呼び、出力を解析する。
#   lowlevel … PDL::IO::HDF5 の低水準関数で直接読む(h5dump 不要)。
#              H5Dread(H5T_STD_REF_OBJ) で参照の列を読み、各参照を
#              H5Rdereference(id, 0 = H5R_OBJECT, ref) で開いて H5T_NATIVE_UINT16 で読む。
#
# read_cell_strings は、h5dump が PATH(または h5dump => で与えた場所)にあれば h5dump、
# 無ければ lowlevel を使う。method => 'h5dump' | 'lowlevel' で固定できる。
# ============================================================================
use strict; use warnings; use Carp;
use File::Spec;
our $VERSION = '0.05';

sub new {
    my ($c, %a) = @_;
    my $m = $a{method} // 'auto';
    $m =~ /^(auto|h5dump|lowlevel)$/ or croak "method must be auto, h5dump or lowlevel";
    bless { h5dump => $a{h5dump} || 'h5dump', method => $m }, $c;
}

# h5dump の実行ファイルが見つかるか(パスを含むならそのファイル、名前だけなら PATH を探す)
sub have_h5dump {
    my $self = shift;
    my $exe = $self->{h5dump};
    return(-x $exe ? 1 : 0) if $exe =~ m{/};
    for my $d (File::Spec->path) { return 1 if -x File::Spec->catfile($d, $exe); }
    return 0;
}

# 実際に使う読み方('h5dump' か 'lowlevel')
sub method {
    my $self = shift;
    return $self->{method} if $self->{method} ne 'auto';
    return $self->have_h5dump ? 'h5dump' : 'lowlevel';
}

sub read_cell_strings {
    my ($self, $file, $dataset) = @_;
    -r $file or croak "cannot read '$file'";
    return $self->method eq 'h5dump'
         ? $self->_read_h5dump($file, $dataset)
         : read_cell_strings_lowlevel($file, $dataset);
}

sub _read_h5dump {
    my ($self,$file,$dataset)=@_;
    my $cmd=sprintf(q{%s -d "%s" "%s"},$self->{h5dump},$dataset,$file);
    open(my $fh,'-|',$cmd) or croak "cannot run '$cmd': $!";
    my (@labels,@codes); my $reading=0;
    while (<$fh>) {
        # 参照先 (/#refs#/...) が現れたら直前セルを確定。
        # h5dump の版で `DATASET <id> "..."` / `DATASET "..."` の差があるが
        # どちらも #refs#/ を含むので、その出現だけで区切れる(版非依存)。
        if (m{#refs#/}) {
            push @labels, join('',map{chr($_)}@codes) if @codes;
            @codes=(); $reading=0; next;
        }
        if ($reading==0 && /^\s*DATA\s*\{/) { $reading=1; next; }
        next unless $reading;
        if (/^\s*\}/) { $reading=0; next; }
        while (/\(\d+,\d+\):\s*(\d+)/g) { push @codes,$1; }
    }
    push @labels, join('',map{chr($_)}@codes) if @codes;
    close($fh);
    return \@labels;
}

# ---- 低水準(h5dump 不要)---------------------------------------------------
# 引数: ファイル名(または PDL::IO::HDF5 オブジェクト)、データセットのパス('/sa/clab_electrodes')。
# group・dataset のオブジェクトは捨てた時点で ID を閉じるので、読み終わるまで変数に持つ。
sub read_cell_strings_lowlevel {
    my ($file, $dataset) = @_;
    require PDL::IO::HDF5;
    my $h5 = ref $file ? $file : PDL::IO::HDF5->new($file);
    my ($gpath, $name) = $dataset =~ m{^(.*)/([^/]+)$} or croak "bad dataset path '$dataset'";
    my $grp = $gpath eq '' ? $h5 : $h5->group($gpath);
    my $ds  = $grp->dataset($name);
    my $id  = $ds->IDget;
    defined $id && $id >= 0 or croak "cannot open $dataset";
    my $tid = PDL::IO::HDF5::H5Dget_type($id);
    my $cls = PDL::IO::HDF5::H5Tget_class($tid);
    PDL::IO::HDF5::H5Tclose($tid);
    $cls == PDL::IO::HDF5::H5T_REFERENCE()
        or croak "$dataset: not an array of object references (not a MATLAB cell array)";

    my $n    = _npoints($id);
    my $rtyp = PDL::IO::HDF5::H5T_STD_REF_OBJ();
    my $rsz  = PDL::IO::HDF5::H5Tget_size($rtyp);
    my $buf  = "\0" x ($rsz * $n);
    PDL::IO::HDF5::H5Dread($id, $rtyp, PDL::IO::HDF5::H5S_ALL(), PDL::IO::HDF5::H5S_ALL(),
                           PDL::IO::HDF5::H5P_DEFAULT(), $buf) >= 0
        or croak "$dataset: cannot read object references (not a cell array?)";

    my $u16 = PDL::IO::HDF5::H5T_NATIVE_UINT16();
    my @out;
    for my $i (0 .. $n - 1) {
        # 第 2 引数 0 = H5R_OBJECT(モジュールが定数で出すのは H5R_DATASET_REGION = 1 だけ)
        my $oid = PDL::IO::HDF5::H5Rdereference($id, 0, substr($buf, $i * $rsz, $rsz));
        $oid >= 0 or croak "$dataset: cannot dereference cell $i";
        my $m  = _npoints($oid);
        my $cb = "\0" x (2 * $m);
        my $rc = $m ? PDL::IO::HDF5::H5Dread($oid, $u16, PDL::IO::HDF5::H5S_ALL(), PDL::IO::HDF5::H5S_ALL(),
                                             PDL::IO::HDF5::H5P_DEFAULT(), $cb) : 0;
        PDL::IO::HDF5::H5Dclose($oid);
        $rc >= 0 or croak "$dataset: cannot read cell $i";
        # 空の文字列は MATLAB では [0 0] として保存される(MATLAB_empty 属性)。0 を落とす。
        push @out, pack('U*', grep { $_ } unpack('S*', $cb));
    }
    return \@out;
}

# データセット ID の要素数(全次元の積。スカラーは 1)
sub _npoints {
    my $id = shift;
    my $sp = PDL::IO::HDF5::H5Dget_space($id);
    my $nd = PDL::IO::HDF5::H5Sget_simple_extent_ndims($sp);
    my $n  = 1;
    if ($nd > 0) {
        my $d  = PDL::IO::HDF5::packList((0) x $nd);
        my $d2 = PDL::IO::HDF5::packList((0) x $nd);
        PDL::IO::HDF5::H5Sget_simple_extent_dims($sp, $d, $d2) == $nd
            or croak "cannot get dataspace dims";
        $n *= $_ for PDL::IO::HDF5::unpackList($d);
    }
    PDL::IO::HDF5::H5Sclose($sp);
    return $n;
}

1;

__END__

=encoding utf-8

=head1 NAME

PDL::IO::NYHead::H5Dump - read MATLAB v7.3 cell-of-strings (HDF5 object references)

=head1 SYNOPSIS

    my $r = PDL::IO::NYHead::H5Dump->new;                 # h5dump if found, else low-level
    my $labels = $r->read_cell_strings('sa_nyhead.mat', '/sa/clab_electrodes');
    print $r->method, "\n";                               # 'h5dump' or 'lowlevel'

    PDL::IO::NYHead::H5Dump->new(method => 'lowlevel');   # never call h5dump
    PDL::IO::NYHead::H5Dump->new(method => 'h5dump', h5dump => '/opt/local/bin/h5dump');

    my $l = PDL::IO::NYHead::H5Dump::read_cell_strings_lowlevel($file_or_h5obj, $dataset);

=head1 DESCRIPTION

C<read_cell_strings> uses the C<h5dump> executable when it is found (on C<PATH>, or
at the path given with C<< h5dump => >>) and otherwise reads the object references
directly with the low-level functions of L<PDL::IO::HDF5>.

=cut
