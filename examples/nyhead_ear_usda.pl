#!/usr/bin/env perl
# nyhead_ear_usda.pl — 外耳道電極の位置検討用に、NY Head の関係ジオメトリを1本の .usda に
# 書き出す(Keynote / usdview / Blender でグルグル)。
#
#   skin (/sa/head)        : 三角メッシュの稜線だけ(BasisCurves ワイヤ)=塗らない・透ける籠
#   BEM shells (#refs#/b,c,d): しっかり塗る(不透明・濃い白/グレー、doubleSided)
#   電極群                 : 透明球(既定 透明度 0.3)
#   snapped / 指示点 / 候補 : 塗りつぶし球(不透明・色分け)
#   軸                     : 赤=X 緑=Y 青=Z の BasisCurves(--axes mm)
#
# **重要な限界**: NY Head の「骨」は順問題用に単純化された 3 層 BEM シェル(滑らかな殻)で、
# 外耳道の穴は開いていない。よってこの図で見えるのは「候補点が skin/頭蓋シェルに対して
# 外か内か・どの深さか」まで。骨の外耳道そのものとの近さは MRI ボリューム(.nii)にしかない。
#
# 使い方:
#   perl nyhead_ear_usda.pl --mat ~/src/NYHead/sa_nyhead.mat -o ear.usda
#   perl nyhead_ear_usda.pl --self-test -o /tmp/st.usda      # 合成データで書式確認(.mat不要)
#   perl nyhead_ear_usda.pl ... --ear "-65,-25,-60" --candidate LPA,Ex19
#   perl nyhead_ear_usda.pl ... --no-bem-shells              # BEM 読めない時は skin+球だけ
#
# 依存: PDL, PDL::IO::HDF5(--self-test 時は不要)。
use strict;
use warnings;
use Getopt::Long;
use PDL;

my %o = (
    mat       => "$ENV{HOME}/src/NYHead/sa_nyhead.mat",
    skin      => '/sa/head',
    bem       => '/#refs#/b,/#refs#/c,/#refs#/d',   # 3層 BEM シェル(vc/tri)。h5ls -r で要確認
    out       => 'nyhead_ear.usda',
    ear       => '-65,-25,-60',            # 私(元)の指示点(MNI mm)
    target    => '-77.2,-24.2,-59.1',      # snapped head 頂点
    candidate => 'LPA',                    # 候補電極(塗り球で強調)。カンマ区切りで複数可
    elec_r    => 3.0,                      # 電極球 半径 mm
    mark_r    => 4.0,                      # 塗り球 半径 mm
    elec_op   => 0.30,                     # 電極群の透明度(displayOpacity)
    axes      => 60,                       # 軸長 mm(0=描かない)
    bem_on    => 1,
    self_test => 0,
);
GetOptions(
    'mat=s'=>\$o{mat}, 'skin=s'=>\$o{skin}, 'bem=s'=>\$o{bem}, 'out|o=s'=>\$o{out},
    'ear=s'=>\$o{ear}, 'target=s'=>\$o{target}, 'candidate=s'=>\$o{candidate},
    'elec-r=f'=>\$o{elec_r}, 'mark-r=f'=>\$o{mark_r}, 'elec-op=f'=>\$o{elec_op},
    'axes=f'=>\$o{axes}, 'bem-shells!'=>\$o{bem_on}, 'self-test'=>\$o{self_test},
) or die "bad args\n";

# 近傍電極(名前 => MNI mm)。seek_ext_ear.pl の距離順上位。透明球で全部出す。
my @ELEC = (
    [ 'Ex19', -82.5, -18.5, -62.0 ],
    [ 'Ex21', -81.5,  -4.0, -63.0 ],
    [ 'LPA',  -83.5, -17.5, -38.5 ],
    [ 'Ex13', -74.5, -56.5, -60.5 ],
    [ 'Ex23', -79.5,   9.0, -64.5 ],
);

# ---------------- HDF5 サーフェスローダ(demo_gs3d_nyhead.pl 流用) ----------------
sub _as_Nx3 {
    my ($p,$n)=@_; my @d=$p->dims;
    die "$n: expected 2-D (@d)\n" unless @d==2;
    return $p if $d[1]==3; return $p->transpose->copy if $d[0]==3;
    die "$n: no len-3 axis (@d)\n";
}
sub _face_base { my ($tri,$nv)=@_; return 0 if $tri->min < 0.5; return 1; }
sub load_surface {
    my ($file,$path)=@_; $path=~s{/$}{};
    require PDL::IO::HDF5; PDL::IO::HDF5->import;
    my $h=PDL::IO::HDF5->new($file);
    my ($vc,$tri);
    eval { $vc=_as_Nx3($h->dataset("$path/vc")->get,"$path/vc");
           $tri=_as_Nx3($h->dataset("$path/tri")->get,"$path/tri"); 1 }
        or die "could not read $path/vc + $path/tri: $@";
    my $nv=$vc->dim(0);
    my $verts=$vc->copy;                                  # (Nv,3)
    my $based = _face_base($tri,$nv);                     # 0=0-based, 1=1-based
    my $faces = ($tri+0.5)->long;                         # 読値を整数化
    $faces = $faces - 1 if $based == 1;                   # 1-based -> 0-based
    printf STDERR "surface %s: %d verts, %d faces (%s)\n",$path,$nv,$tri->dim(0),$based?'1-based':'0-based';
    return ($verts, $faces);   # faces は 0-based (Nf,3)
}

# ---------------- USDA 断片ジェネレータ ----------------
sub _san { my $s=shift; $s=~s/[^A-Za-z0-9_]/_/g; $s='_'.$s if $s=~/^[0-9]/; $s }
sub _pts3 {   # (Nv,3) piddle -> "(x,y,z), (x,y,z), ..."
    my $v=shift; my @x=$v->slice(':,(0)')->list; my @y=$v->slice(':,(1)')->list; my @z=$v->slice(':,(2)')->list;
    join(", ", map { sprintf "(%.4f, %.4f, %.4f)",$x[$_],$y[$_],$z[$_] } 0..$#x);
}
sub _extent {
    my $v=shift;
    sprintf "[(%.4f, %.4f, %.4f), (%.4f, %.4f, %.4f)]",
        $v->slice(':,(0)')->min,$v->slice(':,(1)')->min,$v->slice(':,(2)')->min,
        $v->slice(':,(0)')->max,$v->slice(':,(1)')->max,$v->slice(':,(2)')->max;
}
sub mesh_solid {   # 塗りメッシュ(不透明 or 半透明)
    my (%a)=@_;
    my ($name,$V,$F,$c,$op)=@a{qw(name verts faces color opacity)};
    my @fc = ('3') x $F->dim(0);
    my @fi = ($F->transpose->flat->list);
    my $col = sprintf "(%.3f, %.3f, %.3f)", @$c;
    my $out = qq!    def Mesh "${\_san($name)}"\n    {\n!;
    $out .= "        point3f[] points = [". _pts3($V) ."]\n";
    $out .= "        int[] faceVertexCounts = [". join(", ",@fc) ."]\n";
    $out .= "        int[] faceVertexIndices = [". join(", ",@fi) ."]\n";
    $out .= "        color3f[] primvars:displayColor = [$col] (interpolation = \"constant\")\n";
    $out .= sprintf "        float[] primvars:displayOpacity = [%.3f] (interpolation = \"constant\")\n",$op;
    $out .= "        uniform bool doubleSided = true\n";
    $out .= "        uniform token subdivisionScheme = \"none\"\n";
    $out .= "        float3[] extent = "._extent($V)."\n";
    $out .= "    }\n";
    return $out;
}
sub wireframe {   # 三角メッシュの稜線を BasisCurves で(=塗らない籠)
    my (%a)=@_;
    my ($name,$V,$F,$c,$w)=@a{qw(name verts faces color width)};
    my %seen; my @edges;
    my @fi = $F->transpose->flat->list;
    for (my $t=0; $t<@fi; $t+=3) {
        my ($i,$j,$k)=@fi[$t,$t+1,$t+2];
        for my $e ([$i,$j],[$j,$k],[$k,$i]) {
            my ($a,$b)= $e->[0]<$e->[1] ? @$e : ($e->[1],$e->[0]);
            next if $seen{"$a-$b"}++; push @edges,[$a,$b];
        }
    }
    my @x=$V->slice(':,(0)')->list; my @y=$V->slice(':,(1)')->list; my @z=$V->slice(':,(2)')->list;
    my @cvc = ('2') x scalar(@edges);
    my @pts; for my $e (@edges) { push @pts, sprintf("(%.4f, %.4f, %.4f)",$x[$e->[0]],$y[$e->[0]],$z[$e->[0]]),
                                              sprintf("(%.4f, %.4f, %.4f)",$x[$e->[1]],$y[$e->[1]],$z[$e->[1]]); }
    my $col = sprintf "(%.3f, %.3f, %.3f)", @$c;
    my $out = qq!    def BasisCurves "${\_san($name)}"\n    {\n!;
    $out .= "        uniform token type = \"linear\"\n";
    $out .= "        int[] curveVertexCounts = [". join(", ",@cvc) ."]\n";
    $out .= "        point3f[] points = [". join(", ",@pts) ."]\n";
    $out .= sprintf "        float[] widths = [%.3f] (interpolation = \"constant\")\n",$w;
    $out .= "        color3f[] primvars:displayColor = [$col] (interpolation = \"constant\")\n";
    $out .= "    }\n";
    printf STDERR "wireframe %s: %d edges\n",$name,scalar(@edges);
    return $out;
}
sub sphere {   # 球(電極=透明 / マーク=塗り)
    my (%a)=@_;
    my ($name,$x,$y,$z,$r,$c,$op)=@a{qw(name x y z r color opacity)};
    my $col=sprintf "(%.3f, %.3f, %.3f)",@$c;
    my $out = qq!        def Sphere "${\_san($name)}"\n        {\n!;
    $out .= sprintf "            double radius = %.3f\n",$r;
    $out .= "            color3f[] primvars:displayColor = [$col] (interpolation = \"constant\")\n";
    $out .= sprintf "            float[] primvars:displayOpacity = [%.3f] (interpolation = \"constant\")\n",$op;
    $out .= sprintf "            double3 xformOp:translate = (%.4f, %.4f, %.4f)\n",$x,$y,$z;
    $out .= "            uniform token[] xformOpOrder = [\"xformOp:translate\"]\n";
    $out .= "        }\n";
    return $out;
}
sub axes_curves {
    my $L=shift; return "" if $L<=0;
    my @seg=([[0,0,0],[$L,0,0],[1,0,0],'X'],[[0,0,0],[0,$L,0],[0,1,0],'Y'],[[0,0,0],[0,0,$L],[0,0,1],'Z']);
    my $out="    def BasisCurves \"axes\"\n    {\n        uniform token type = \"linear\"\n";
    $out.="        int[] curveVertexCounts = [2, 2, 2]\n";
    my @pts; for my $s (@seg){ push @pts, sprintf("(%.1f, %.1f, %.1f)",@{$s->[0]}), sprintf("(%.1f, %.1f, %.1f)",@{$s->[1]}); }
    $out.="        point3f[] points = [".join(", ",@pts)."]\n";
    $out.="        float[] widths = [1.5] (interpolation = \"constant\")\n";
    my @col; for my $s (@seg){ push @col, sprintf("(%.1f, %.1f, %.1f)",@{$s->[2]}) for 1..2; }
    $out.="        color3f[] primvars:displayColor = [".join(", ",@col)."] (interpolation = \"vertex\")\n";
    $out.="    }\n";
    return $out;
}

# ---------------- ジオメトリ取得 ----------------
sub _mkmesh {   # 合成: 四面体 verts(4,3), faces0(4,3) 0-based
    my ($cx,$cy,$cz,$s)=@_;
    my $V=pdl([[ $cx+$s,$cy,$cz],[ $cx-$s,$cy+$s,$cz],[ $cx-$s,$cy-$s,$cz+$s],[ $cx-$s,$cy-$s,$cz-$s]])->transpose->sever;  # (4,3)
    my $F=pdl([[0,1,2],[0,2,3],[0,3,1],[1,2,3]])->transpose->long->sever;   # (4,3) 面ごと
    return ($V,$F);
}

my ($skinV,$skinF, @shells, @elec, @marks);
if ($o{self_test}) {
    ($skinV,$skinF)=_mkmesh(0,0,0, 80);
    my ($bV,$bF)=_mkmesh(0,0,0, 70); push @shells,['shell_b',$bV,$bF];
    @elec = map { [$_->[0],$_->[1],$_->[2],$_->[3]] } @ELEC;
} else {
    -f $o{mat} or die "not found: $o{mat}\n";
    ($skinV,$skinF)=load_surface($o{mat},$o{skin});
    if ($o{bem_on}) {
        for my $p (split /,/,$o{bem}) {
            my ($V,$F);
            my $ok = eval { ($V,$F)=load_surface($o{mat},$p); 1 };
            if ($ok) { (my $nm=$p)=~s{.*/}{}; push @shells,["shell_$nm",$V,$F]; }
            else { warn "BEM shell '$p' 読めず(skip): $@  → h5ls -r $o{mat} | grep -iE 'refs|vc|tri' で実パス確認、--bem で指定\n"; }
        }
    }
    @elec = map { [$_->[0],$_->[1],$_->[2],$_->[3]] } @ELEC;
}

# マーク(塗り球): snapped(赤) / 指示点(マゼンタ) / 候補電極(橙)
sub _xyz { my @p=split /,/,$_[0]; $_+=0 for @p; die "bad coord '$_[0]'\n" unless @p==3; @p }
{
    my @t=_xyz($o{target}); push @marks,['target',@t, [1.0,0.15,0.15]];
    my @e=_xyz($o{ear});    push @marks,['instructed',@e, [1.0,0.10,0.90]];
    my %cand = map { $_=>1 } split /,/,$o{candidate};
    for my $el (@ELEC) { push @marks,["cand_$el->[0]",$el->[1],$el->[2],$el->[3],[1.0,0.55,0.0]] if $cand{$el->[0]}; }
}

# ---------------- 出力組み立て ----------------
my $body = "";
$body .= wireframe(name=>'skin', verts=>$skinV, faces=>$skinF, color=>[0.70,0.70,0.78], width=>0.6);
my @greys=([0.92,0.92,0.92],[0.80,0.80,0.82],[0.66,0.66,0.70]);
for my $i (0..$#shells) {
    my ($nm,$V,$F)=@{$shells[$i]};
    $body .= mesh_solid(name=>$nm, verts=>$V, faces=>$F, color=>$greys[$i % @greys], opacity=>1.0);
}
$body .= axes_curves($o{axes});

# 電極群(透明球)
my $eblock = "    def Xform \"electrodes\"\n    {\n";
for my $el (@elec) { $eblock .= sphere(name=>$el->[0], x=>$el->[1],y=>$el->[2],z=>$el->[3], r=>$o{elec_r}, color=>[1.0,0.85,0.10], opacity=>$o{elec_op}); }
$eblock .= "    }\n";
$body .= $eblock;

# マーク(塗り球)
my $mblock = "    def Xform \"marks\"\n    {\n";
for my $m (@marks) { $mblock .= sphere(name=>$m->[0], x=>$m->[1],y=>$m->[2],z=>$m->[3], r=>$o{mark_r}, color=>$m->[4], opacity=>1.0); }
$mblock .= "    }\n";
$body .= $mblock;

my $usda = qq{#usda 1.0\n(\n    upAxis = "Z"\n    metersPerUnit = 0.001\n    defaultPrim = "NYHeadEar"\n)\n\n};
$usda .= qq{def Xform "NYHeadEar"\n{\n$body}\n};

open my $fh,'>',$o{out} or die "open $o{out}: $!";
print $fh $usda; close $fh;
printf STDERR "wrote %s (%d bytes): skin wire + %d shell(s) + %d elec + %d mark\n",
    $o{out}, -s $o{out}, scalar(@shells), scalar(@elec), scalar(@marks);
