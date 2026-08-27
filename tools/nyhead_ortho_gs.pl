#!/usr/bin/perl
# nyhead_ortho_gs.pl — giza 対話 2x2 ortho ビューア(+ズーム)。
#   左上 axial / 右上 coronal / 左下 sagittal / 右下 ヒストグラム or オーバービュー。
#   * クリック → その断面の点へ全ヘアライン移動(=パン) + 端末に `--eec X,Y,Z` 出力
#   * スライダは下=id0(level=窓中心)・右=id1(width=窓幅)の2本。ズーム率は起動時 --zoom で
#     指定する(現行ビューアにズーム用スライダは無く、セッション中は --zoom の値で固定)。
#   * 内部空気(外耳道等)をシアン重畳、電極を各断面に表示
#   * ズーム>1 のとき: 各断面は十字中心の周りを full-res で切り出して拡大表示、
#     右下ペインは3面オーバービュー(全体像+現在のズーム枠)に切替。
# giza-server 再ビルド不要(PICK/CURSOR/SLIDER は既存配線)。実機は Mac。
#
#   perl nyhead_ortho_gs.pl --mat sa_nyhead.mat --center -65,-25,-60 --eec -65,-25,-60 --nyhead Ex19,LPA
#   perl nyhead_ortho_gs.pl --self-test --probe coronal,0.60,0.52 --zoom 4   # 対話なしで pick+zoom を検証

use strict; use warnings; use PDL; use PDL::Image2D; use Getopt::Long;
use lib '/Users/goosh/src/PDL_IO_NYHead/lib/';

# ------------------------------------------------------------------ options
my %o = (
    mat=>undef, selftest=>0, center=>undef, eec=>undef, elec=>[], nyhead=>undef,
    overlay=>'air', air_side=>'auto', air_color=>"0,0.9,1", bg=>'auto',
    radiological=>0, max_dim=>260, elec_orig=>1, zoom=>1,
    probe=>undef, out=>"nyhead_ortho.png", debug=>1,
);
GetOptions(
    "mat=s"=>\$o{mat}, "self-test!"=>\$o{selftest}, "center=s"=>\$o{center},
    "eec=s"=>\$o{eec}, "elec=s@"=>$o{elec}, "nyhead=s"=>\$o{nyhead},
    "overlay=s"=>\$o{overlay}, "air-side=s"=>\$o{air_side}, "air-color=s"=>\$o{air_color},
    "bg=s"=>\$o{bg}, "radiological!"=>\$o{radiological}, "max-dim=i"=>\$o{max_dim},
    "elec-orig!"=>\$o{elec_orig}, "zoom=f"=>\$o{zoom},
    "probe=s"=>\$o{probe}, "out=s"=>\$o{out}, "debug!"=>\$o{debug},
) or die "bad options\n";
die "give --mat or --self-test\n" unless $o{mat} || $o{selftest};

# ズーム率の範囲(スライダ id2 は [0,1] を対数で ZMIN..ZMAX に写す)
use constant { ZMIN=>1.0, ZMAX=>12.0 };
sub zoom_from_s2 { my($s)=@_; $s=0 if !defined $s || $s<0; $s=1 if $s>1;
                   ZMIN*(ZMAX/ZMIN)**$s }
sub s2_from_zoom { my($z)=@_; $z=ZMIN if $z<ZMIN; $z=ZMAX if $z>ZMAX;
                   log($z/ZMIN)/log(ZMAX/ZMIN) }

# ------------------------------------------------------------------ helpers (engine)
sub parse_xyz { my @v=split /\s*,\s*/,$_[0]; die "want x,y,z\n" unless @v==3; map {0+$_} @v }
sub pctl { my($p,$q)=@_; my $s=$p->flat->qsort; my $n=$s->nelem; return 0 unless $n;
           my $i=int($q*($n-1)+0.5); $i=0 if $i<0; $i=$n-1 if $i>$n-1; $s->at($i) }
sub apply_affine { my($M,$x,$y,$z)=@_; map { $M->at($_,0)*$x+$M->at($_,1)*$y+$M->at($_,2)*$z+$M->at($_,3) } (0,1,2) }
sub to_math_affine { my($A)=@_; my($best,$bd)=(undef,1e30);
    for my $M ($A->sever,$A->transpose->sever){ my $b=abs($M->at(3,0))+abs($M->at(3,1))+abs($M->at(3,2));
        my $one=abs($M->at(3,3)-1); return $M if $b<1e-3&&$one<1e-3; ($best,$bd)=($M,$b+$one) if $b+$one<$bd } $best }
sub h5get { my($h5,$path)=@_; my @p=grep{length}split m{/},$path; my $ds=pop @p;
            my $n=$h5; $n=$n->group($_) for @p; $n->dataset($ds)->get }
sub air_components { my($S,$hi,$thr)=@_;
    my $air=($hi?($S>$thr):($S<$thr))->byte;
    return (zeroes(byte,$S->dims),zeroes(byte,$S->dims)) if $air->sum==0;
    my $lab=cc8compt($air); my($H,$W)=$lab->dims;
    my $b=$lab->slice("0,:")->flat->append($lab->slice("-1,:")->flat)
         ->append($lab->slice(":,0")->flat)->append($lab->slice(":,-1")->flat);
    my $ext=zeroes(byte,$H,$W); for my $L (grep{$_>0} $b->uniq->list){ $ext=$ext|($lab==$L) }
    ($ext->byte, ($air&($ext==0))->byte) }
sub auto_window { my($vals)=@_;
    my @los=map{pctl($vals,$_)} (0.005,0.01,0.02,0.05);
    my @his=map{pctl($vals,$_)} (0.95,0.98,0.99,0.995);
    my($best,$bE)=([$los[2],$his[1]],-1);
    for my $lo (@los){ for my $hi (@his){ next unless $hi>$lo;
        my $g=(($vals-$lo)/($hi-$lo))->clip(0,1); my $h=histogram($g,1/64,0,64)->double;
        my $p=$h/($h->sum||1); my $E=-($p*($p+1e-12)->log)->sum; ($best,$bE)=([$lo,$hi],$E) if $E>$bE }}
    @$best }

# ------------------------------------------------------------------ load
my ($vol,$mri2mni,$mni2mri);
if ($o{selftest}) {
    my ($Ni,$Nj,$Nk)=(80,90,70); my @ctr=(40,45,35);
    my $ii=xvals($Ni,$Nj,$Nk); my $jj=yvals($Ni,$Nj,$Nk); my $kk=zvals($Ni,$Nj,$Nk);
    my $rr=sqrt((($ii-$ctr[0])/30)**2+(($jj-$ctr[1])/34)**2+(($kk-$ctr[2])/26)**2);
    my $brain=($rr<0.86); my $skull=(($rr>=0.86)&($rr<1.0)); my $scalp=(($rr>=1.0)&($rr<1.08));
    $vol=$brain*800+$skull*150+$scalp*600+$brain*(100*sin($rr*12));
    for my $cx (14,66){ my $da=sqrt(($ii-$cx)**2+($jj-45)**2+($kk-32)**2); $vol=$vol*($da>=4) }
    $mri2mni=zeroes(double,4,4); $mri2mni->set(3,3,1);
    for my $r (0..2){ $mri2mni->set($r,$r,0.5); $mri2mni->set($r,3,-0.5*$ctr[$r]) }
    $mni2mri=zeroes(double,4,4); $mni2mri->set(3,3,1);
    for my $r (0..2){ $mni2mri->set($r,$r,2); $mni2mri->set($r,3,$ctr[$r]) }
    $o{eec}//="-14,0,-2"; $o{center}//="0,0,0";
} else {
    print STDERR "reading $o{mat} ...\n"; require PDL::IO::HDF5;
    my $h5=PDL::IO::HDF5->new($o{mat}) or die "cannot open $o{mat}\n";
    $vol=h5get($h5,'/sa/mri/data')->double;
    $mri2mni=to_math_affine(h5get($h5,'/sa/mri2mni')->double);
    $mni2mri=to_math_affine(h5get($h5,'/sa/mni2mri')->double);
}
my ($Ni,$Nj,$Nk)=$vol->dims;
my @vsize=map{ my $c=$_; sqrt($mri2mni->at(0,$c)**2+$mri2mni->at(1,$c)**2+$mri2mni->at(2,$c)**2) } (0,1,2);

# ------------------------------------------------------------------ bg / air / window params
my ($vmin,$vmax,$range)=($vol->min,$vol->max,($vol->max-$vol->min)||1);
my $bg_mode=$o{bg};
if ($bg_mode eq 'auto'){ my $lo=($vol<($vmin+0.01*$range))->sum; my $hi=($vol>($vmax-0.01*$range))->sum;
    $bg_mode=($hi>$lo)?'high':'low' }
my $bg_thr=($bg_mode eq 'high')?$vmax-0.03*$range:($bg_mode eq 'low')?$vmin+0.03*$range:undef;
my $air_high=($o{air_side} eq 'high')?1:($o{air_side} eq 'low')?0:($bg_mode eq 'high'?1:0);
my $air_thr=$air_high?$vmax-0.03*$range:$vmin+0.03*$range;
my $mask=($bg_mode eq 'high')?($vol<($vmax-0.05*$range)):($bg_mode eq 'low')?($vol>($vmin+0.05*$range)):ones(byte,$vol->dims);
my $headvals=$vol->flat->where($mask->flat); $headvals=$vol->flat if $headvals->nelem<100;
my @acol=parse_xyz($o{air_color});

# ヒストグラム(頭部ボクセル分布)は不変 → 起動時に1回だけ計算
my ($HLO,$HHI,$HX,$HY);
{
    $HLO=pctl($headvals,0.005); $HHI=pctl($headvals,0.995);
    my $nb=80; my $hstep=(($HHI-$HLO)||1)/$nb;
    my $hh=histogram($headvals,$hstep,$HLO,$nb)->double;
    $HX=$HLO+$hstep*(sequence($nb)+0.5); $HY=($hh+1)->log;
}

# 初期窓(auto-contrast) → スライダ正規化 [0,1]
my ($aw_lo,$aw_hi)=auto_window($headvals);
my $init_level=($aw_lo+$aw_hi)/2; my $init_width=($aw_hi-$aw_lo);
my $s0_init=($init_level-$vmin)/$range; my $s1_init=$init_width/$range;
$s0_init=0 if $s0_init<0; $s0_init=1 if $s0_init>1; $s1_init=0.02 if $s1_init<0.02; $s1_init=1 if $s1_init>1;
my $s2_init=s2_from_zoom($o{zoom});
my $ZOOM=$o{zoom};   # 直近レンダのズーム率(pick 解決が同じ窓を使うため global)

# ------------------------------------------------------------------ electrodes
my @elecs;
if (defined $o{eec}){ my @p=parse_xyz($o{eec}); push @elecs,{label=>"EEC",x=>$p[0],y=>$p[1],z=>$p[2],type=>"eec"} }
for my $spec (@{$o{elec}}){ my($l,$r)=split/=/,$spec,2; my @p=parse_xyz($r);
    push @elecs,{label=>$l,x=>$p[0],y=>$p[1],z=>$p[2],type=>"user"} }
if (defined $o{nyhead} && !$o{selftest}){ my @want=split/\s*,\s*/,$o{nyhead};
    eval { require PDL::IO::HDF5; require PDL::IO::NYHead; my $ny=PDL::IO::NYHead->new($o{mat});
        my $lab=$ny->electrode_labels; my $pos=$o{elec_orig}?$ny->electrode_pos
              :h5get(PDL::IO::HDF5->new($o{mat}),'/sa/locs_3D')->slice(':,0:2')->sever;
        my %idx; $idx{$lab->[$_]}=$_ for 0..$#$lab;
        for my $nm (@want){ exists $idx{$nm} or next; my $r=$idx{$nm};
            push @elecs,{label=>$nm,x=>$pos->at($r,0),y=>$pos->at($r,1),z=>$pos->at($r,2),type=>"nyhead"} } 1 }
    or warn "  --nyhead failed ($@)\n" }

# ------------------------------------------------------------------ geometry (center-dependent)
my $lr=$o{radiological}?1:0;
my @center_mni = defined $o{center}?parse_xyz($o{center}) : defined $o{eec}?parse_xyz($o{eec}) : (0,0,0);
my @anchor_mni = defined $o{eec} ? parse_xyz($o{eec}) : ();   # original(起動時 --eec)。クリックで動かさない
sub dist3 { my($a,$b)=@_; sqrt(($a->[0]-$b->[0])**2+($a->[1]-$b->[1])**2+($a->[2]-$b->[2])**2) }
my @cidx;
sub set_center_mni { @center_mni=@_;
    my @v=apply_affine($mni2mri,@center_mni);
    @cidx=map{ my $x=int($v[$_]+0.5); $x=0 if $x<0; my @N=($Ni,$Nj,$Nk); $x=$N[$_]-1 if $x>$N[$_]-1; $x } (0,1,2) }
set_center_mni(@center_mni);
sub plane_spec {
    my ($name)=@_;
    my %P=(
      axial   =>{fixed=>2,cidx=>$cidx[2],col_axis=>0,col_flip=>$lr,row_axis=>1,row_flip=>1,
                 title=>sprintf("Axial  z=%.0f",$center_mni[2]),
                 labels=>$o{radiological}?{top=>"A",bottom=>"P",left=>"R",right=>"L"}:{top=>"A",bottom=>"P",left=>"L",right=>"R"}},
      coronal =>{fixed=>1,cidx=>$cidx[1],col_axis=>0,col_flip=>$lr,row_axis=>2,row_flip=>1,
                 title=>sprintf("Coronal  y=%.0f",$center_mni[1]),
                 labels=>$o{radiological}?{top=>"S",bottom=>"I",left=>"R",right=>"L"}:{top=>"S",bottom=>"I",left=>"L",right=>"R"}},
      sagittal=>{fixed=>0,cidx=>$cidx[0],col_axis=>1,col_flip=>1,row_axis=>2,row_flip=>1,
                 title=>sprintf("Sagittal  x=%.0f",$center_mni[0]),
                 labels=>{top=>"S",bottom=>"I",left=>"A",right=>"P"}},
    );
    $P{$name};
}
sub plane_slice { my($spec)=@_; my @idx=(undef,undef,undef); $idx[$spec->{fixed}]=$spec->{cidx};
    my @ss=map{ defined $idx[$_]?"($idx[$_])":":" } (0,1,2);
    my $s2=$vol->slice(join(",",@ss))->squeeze; my @rem=grep{ $_!=$spec->{fixed} } (0,1,2);
    my ($rd)=grep{ $rem[$_]==$spec->{row_axis} } (0,1);
    my $img=($rd==0)?$s2:$s2->transpose; $img=$img->sever;
    $img=$img->slice("-1:0,:")->sever if $spec->{row_flip};
    $img=$img->slice(":,-1:0")->sever if $spec->{col_flip}; $img }
sub map_vox { my($spec,$i,$j,$k)=@_; my @v=($i,$j,$k); my @N=($Ni,$Nj,$Nk);
    my $col=$v[$spec->{col_axis}]; $col=($N[$spec->{col_axis}]-1)-$col if $spec->{col_flip};
    my $row=$v[$spec->{row_axis}]; $row=($N[$spec->{row_axis}]-1)-$row if $spec->{row_flip};
    ($row,$col,$v[$spec->{fixed}]-$spec->{cidx}) }
sub rc_to_vox { my($spec,$row,$col)=@_; my @N=($Ni,$Nj,$Nk); my @v; $v[$spec->{fixed}]=$spec->{cidx};
    $v[$spec->{col_axis}]=$spec->{col_flip}?($N[$spec->{col_axis}]-1-$col):$col;
    $v[$spec->{row_axis}]=$spec->{row_flip}?($N[$spec->{row_axis}]-1-$row):$row; @v }

# 幅を保ったまま [0,N-1] に収める(窓が断面以上なら全面)
sub clamp_range { my($a,$b,$N)=@_; my $w=$b-$a;
    return (0,$N-1) if $w>=$N-1;
    if ($a<0){ $b-=$a; $a=0 } if ($b>$N-1){ $a-=($b-($N-1)); $b=$N-1 } $a=0 if $a<0; ($a,$b) }

# 断面の full-res 画像 + 全面オーバレイ(cidx が変わったときだけ再計算)
my %SLICE;
sub slice_cache {
    my ($name)=@_; my $spec=plane_spec($name);
    my $c=$SLICE{$name}; return $c if $c && $c->{cidx}==$spec->{cidx};
    my $img=plane_slice($spec); my ($Hf,$Wf)=$img->dims;
    my $ext=zeroes(byte,$Hf,$Wf);
    ($ext,undef)=air_components($img,$bg_mode eq 'high',$bg_thr) if $bg_mode ne 'none';
    my $intair=zeroes(byte,$Hf,$Wf);
    (undef,$intair)=air_components($img,$air_high,$air_thr) if $o{overlay} eq 'air';
    $SLICE{$name}={cidx=>$spec->{cidx},img=>$img,ext=>$ext,intair=>$intair,Hf=>$Hf,Wf=>$Wf};
}

# 十字中心の (row,col) を中心に、zoom 率で窓 [r0..r1]x[c0..c1] を決める(full-res px)
sub window_box {
    my ($name,$zoom)=@_; my $spec=plane_spec($name); my $sc=slice_cache($name);
    my ($Hf,$Wf)=($sc->{Hf},$sc->{Wf});
    $zoom=ZMIN if $zoom<ZMIN;
    my ($crow,$ccol)=map_vox($spec, apply_affine($mni2mri,@center_mni));
    my $hh=0.5*$Hf/$zoom; my $hw=0.5*$Wf/$zoom;
    my ($r0,$r1)=clamp_range(int($crow-$hh+0.5),int($crow+$hh+0.5),$Hf);
    my ($c0,$c1)=clamp_range(int($ccol-$hw+0.5),int($ccol+$hw+0.5),$Wf);
    ($r0,$r1,$c0,$c1);
}

# 窓を切り出し → subsample した表示用パネル(窓が変わったときだけ再計算)
my %CACHE;
sub panel_cache {
    my ($name,$zoom)=@_; $zoom//=1; my $spec=plane_spec($name); my $sc=slice_cache($name);
    my ($r0,$r1,$c0,$c1)=window_box($name,$zoom);
    my $key="$spec->{cidx}:$r0:$r1:$c0:$c1";
    my $c=$CACHE{$name}; return $c if $c && $c->{key} eq $key;
    my $winimg=$sc->{img}->slice("$r0:$r1,$c0:$c1")->sever;
    my ($Hw,$Ww)=$winimg->dims;
    my $step=int((($Hw>$Ww?$Hw:$Ww)+$o{max_dim}-1)/$o{max_dim}); $step=1 if $step<1;
    my $S=($step>1)?$winimg->slice("0:-1:$step,0:-1:$step")->sever:$winimg;
    my ($Hs,$Ws)=$S->dims;
    my $ext   =$sc->{ext}   ->slice("$r0:$r1,$c0:$c1"); my $intair=$sc->{intair}->slice("$r0:$r1,$c0:$c1");
    $ext   =($step>1)?$ext   ->slice("0:-1:$step,0:-1:$step")->sever:$ext->sever;
    $intair=($step>1)?$intair->slice("0:-1:$step,0:-1:$step")->sever:$intair->sever;
    $CACHE{$name}={key=>$key,S=>$S,Hs=>$Hs,Ws=>$Ws,step=>$step,ext=>$ext,intair=>$intair,
                   r0=>$r0,c0=>$c0,r1=>$r1,c1=>$c1,Hf=>$sc->{Hf},Wf=>$sc->{Wf}};
}

# ------------------------------------------------------------------ draw one panel into $ax
sub draw_panel {
    my ($ax,$name,$w_lo,$w_hi,$zoom)=@_;
    my $spec=plane_spec($name); my $pc=panel_cache($name,$zoom);
    my ($S,$Hs,$Ws,$step,$ext,$r0,$c0)=@{$pc}{qw(S Hs Ws step ext r0 c0)};
    my $gn=(($S-$w_lo)/(($w_hi-$w_lo)||1))->clip(0,1);
    my $g=$gn*(1-$ext);
    my ($alpha,@ocol)=(zeroes(float,$Hs,$Ws),(0,0,0));
    if ($o{overlay} eq 'air'){ $alpha=$pc->{intair}*0.55; @ocol=@acol; }
    my $rgb=zeroes(float,$Hs,$Ws,3);
    $rgb->slice(":,:,($_)").= $g*(1-$alpha)+$ocol[$_]*$alpha for (0,1,2);
    $ax->imshow($rgb,origin=>'upper'); $ax->axis('off'); $ax->set_aspect('equal');
    # タイトル: ズーム時は窓の MNI 範囲も表示
    my $title=$spec->{title};
    if ($zoom>1.01) {
        my @a_mni=apply_affine($mri2mni, rc_to_vox($spec,$r0,$c0));
        my @b_mni=apply_affine($mri2mni, rc_to_vox($spec,$pc->{r1},$pc->{c1}));
        my $ca=$spec->{col_axis}; my $ra=$spec->{row_axis};
        my ($cA,$cB)=($a_mni[$ca],$b_mni[$ca]); ($cA,$cB)=($cB,$cA) if $cA>$cB;
        my ($rA,$rB)=($a_mni[$ra],$b_mni[$ra]); ($rA,$rB)=($rB,$rA) if $rA>$rB;
        $title.=sprintf("  x%.1f  [%.0f..%.0f / %.0f..%.0f]",$zoom,$cA,$cB,$rA,$rB);
    }
    $ax->set_title($title);
    my $to_xy=sub{ my($row,$col)=@_; (($col-$c0)/$step,$Hs-($row-$r0)/$step) };
    # 電極(窓外はカリング)
    my @draw;
    for my $e (@elecs){ my @ev=apply_affine($mni2mri,$e->{x},$e->{y},$e->{z});
        my ($row,$col,$off)=map_vox($spec,@ev);
        my ($ex,$ey)=$to_xy->($row,$col);
        next if $ex< -$Ws*0.06 || $ex>$Ws*1.06 || $ey< -$Hs*0.06 || $ey>$Hs*1.06;
        my $off_mm=$off*$vsize[$spec->{fixed}];         # 符号つき: +=このスライスより奥
        my %cof=(eec=>[1,0,0],nyhead=>[0,0.8,1],user=>[1,0.85,0]);
        push @draw,{ex=>$ex,ey=>$ey,c=>($cof{$e->{type}}//[1,0.85,0]),
                    offabs=>abs($off_mm),on=>(abs($off_mm)<=4),
                    tag=>(abs($off_mm)<=4?sprintf("%s  (in plane)",$e->{label})
                                         :sprintf("%s  %+.0fmm",$e->{label},$off_mm))} }
    my ($crow,$ccol)=map_vox($spec,apply_affine($mni2mri,@center_mni));
    my ($cx,$cy)=$to_xy->($crow,$ccol);
    $ax->axvline($cx,color=>[0.2,1,0.4],lw=>0.8,ls=>'dashed');
    $ax->axhline($cy,color=>[0.2,1,0.4],lw=>0.8,ls=>'dashed');
    # 電極/最寄り点: 色線オープンサークル(近いほど不透明)
    for my $d (@draw){ $ax->scatter(pdl($d->{ex}),pdl($d->{ey}),s=>8,color=>$d->{c},marker=>'o',
                        alpha=>($d->{on}?0.95:0.5),fillstyle=>'none') }
    # 現在クリック点のマーカ(緑オープンサークル, 十字の交点)
    $ax->scatter(pdl($cx),pdl($cy),s=>9,color=>[0.2,1,0.4],marker=>'o',fillstyle=>'none',alpha=>0.9);
    # 右下: orig(起動 --eec)と cur(現在クリック点)から各参照電極への直線距離(3D, mm)。
    #        クリックで cur 側だけ書き換わり、orig の場所と計測は残る。
    {
        my @ref = grep { $_->{type} ne 'eec' } @elecs;
        my @lines;   # [text,[r,g,b]] を下から上へ
        for my $e (@ref) {
            my @p=($e->{x},$e->{y},$e->{z});
            my $col=($e->{type} eq 'nyhead')?[0,0.8,1]:[1,0.85,0];
            my $dc=dist3(\@center_mni,\@p);
            push @lines,[ (@anchor_mni ? sprintf("%s  o:%.1f c:%.1f",$e->{label},dist3(\@anchor_mni,\@p),$dc)
                                       : sprintf("%s  c:%.1f",$e->{label},$dc)), $col ];
        }
        push @lines,[ (@anchor_mni ? sprintf("cur (%.0f,%.0f,%.0f)  \x{0394}%.1f",@center_mni,dist3(\@anchor_mni,\@center_mni))
                                   : sprintf("cur (%.0f,%.0f,%.0f)",@center_mni)), [0.2,1,0.4] ];
        push @lines,[ sprintf("orig (%.0f,%.0f,%.0f)",@anchor_mni), [1,0.35,0.35] ] if @anchor_mni;
        my $y=$Hs*0.05;
        for my $ln (@lines){ $ax->text($Ws*0.985,$y,$ln->[0],color=>$ln->[1],fontsize=>9,ha=>'right',va=>'bottom');
                             $y+=$Hs*0.048; }
    }
    # スケールバー(ズーム時, col 方向 10mm)
    if ($zoom>1.01) {
        my $mm_per_px=$step*$vsize[$spec->{col_axis}];
        if ($mm_per_px>0) {
            my ($mm,$len)=(0,0);
            for my $cand (10,5,2,1){ my $l=$cand/$mm_per_px; if ($l>4 && $l<$Ws*0.45){ ($mm,$len)=($cand,$l); last } }
            if ($len>0) {
                my $x0=$Ws*0.06; my $y0=$Hs*0.06;
                $ax->line(pdl($x0,$x0+$len),pdl($y0,$y0),color=>'yellow',lw=>2.0);
                $ax->text($x0+$len/2,$y0+$Hs*0.02,"${mm}mm",color=>'yellow',fontsize=>8,ha=>'center',va=>'bottom');
            }
        }
    }
    my $L=$spec->{labels};
    $ax->text($Ws*0.5,$Hs*0.98,$L->{top},color=>'white',fontsize=>12,ha=>'center',va=>'top');
    $ax->text($Ws*0.5,$Hs*0.02,$L->{bottom},color=>'white',fontsize=>12,ha=>'center',va=>'bottom');
    $ax->text($Ws*0.02,$Hs*0.5,$L->{left},color=>'white',fontsize=>12,ha=>'left',va=>'middle');
    $ax->text($Ws*0.98,$Hs*0.5,$L->{right},color=>'white',fontsize=>12,ha=>'right',va=>'middle');
}

# ------------------------------------------------------------------ overview (右下, ズーム時)
# 3面の全体像を1枚にタイルして imshow し、現在のズーム枠と十字を重ねる。
sub draw_overview {
    my ($ax,$w_lo,$w_hi,$zoom)=@_;
    my $THUMB=120; my $GAP=10; my $TOP=16;
    my (@thumbs,@place); my $x=0; my $maxh=0;
    for my $name (qw(axial coronal sagittal)) {
        my $sc=slice_cache($name);
        my ($Hf,$Wf)=($sc->{Hf},$sc->{Wf});
        my $ts=int((($Hf>$Wf?$Hf:$Wf)+$THUMB-1)/$THUMB); $ts=1 if $ts<1;
        my $t=$sc->{img}->slice("0:-1:$ts,0:-1:$ts")->sever;
        my ($th,$tw)=$t->dims;
        my $gn=(($t-$w_lo)/(($w_hi-$w_lo)||1))->clip(0,1);
        push @thumbs,{name=>$name,g=>$gn,th=>$th,tw=>$tw,ts=>$ts,x=>$x};
        $x+=$tw+$GAP; $maxh=$th if $th>$maxh;
    }
    my $Wc=$x-$GAP; $Wc=1 if $Wc<1; my $Hc=$maxh+$TOP;
    my $comp=zeroes(float,$Hc,$Wc,3)+0.12;
    for my $t (@thumbs) {
        my ($th,$tw,$xo)=@{$t}{qw(th tw x)};
        $comp->slice("$TOP:".($TOP+$th-1).",$xo:".($xo+$tw-1).",($_)") .= $t->{g} for (0,1,2);
    }
    $ax->imshow($comp,origin=>'upper'); $ax->axis('off'); $ax->set_aspect('equal');
    $ax->set_title(sprintf("overview  (zoom x%.1f)",$zoom));
    my $to_xy=sub{ my($r,$c)=@_; ($c,$Hc-$r) };   # comp 座標(row,col)→表示(x,y)
    for my $t (@thumbs) {
        my $name=$t->{name}; my $spec=plane_spec($name); my $ts=$t->{ts}; my $xo=$t->{x};
        my $pc=panel_cache($name,$zoom);
        # ズーム枠(full-res box → thumb → comp)
        my @rr=($pc->{r0},$pc->{r1}); my @cc=($pc->{c0},$pc->{c1});
        my (@xs,@ys);
        for my $cn ([0,0],[0,1],[1,1],[1,0],[0,0]) {
            my ($bxp,$byp)=$to_xy->($TOP+$rr[$cn->[0]]/$ts, $xo+$cc[$cn->[1]]/$ts);
            push @xs,$bxp; push @ys,$byp;
        }
        $ax->line(pdl(@xs),pdl(@ys),color=>[1,0.85,0],lw=>1.2);
        # 十字中心
        my ($crow,$ccol)=map_vox($spec,apply_affine($mni2mri,@center_mni));
        my ($px,$py)=$to_xy->($TOP+$crow/$ts,$xo+$ccol/$ts);
        $ax->scatter(pdl($px),pdl($py),s=>10,color=>[0.2,1,0.4],marker=>'+');
        $ax->text($xo+$t->{tw}/2,$Hc-($TOP-3),substr(ucfirst($name),0,3),
                  color=>'white',fontsize=>9,ha=>'center',va=>'bottom');
    }
}

# ------------------------------------------------------------------ pick 座標 → center 更新
my @PLANE_BY_IDX=qw(axial coronal sagittal);
my $LAST_FIG; my $BUILT_KEY;
sub pick_to_mni {
    my ($ax_idx,$cur_x,$cur_y)=@_;
    return () unless defined $ax_idx && $ax_idx>=0 && $ax_idx<=2;
    return () unless defined $cur_x && defined $cur_y;
    my $name=$PLANE_BY_IDX[$ax_idx]; my $pc=panel_cache($name,$ZOOM); my $spec=plane_spec($name);
    my $col_s=$cur_x; my $row_s=$pc->{Hs}-$cur_y;         # imshow origin upper の y 反転を戻す
    my @v=rc_to_vox($spec, $pc->{r0}+$row_s*$pc->{step}, $pc->{c0}+$col_s*$pc->{step});
    my @mni=apply_affine($mri2mni,@v);
    (\@mni,\@v,$name);
}
sub pick_fallback {
    my ($fx,$fy)=@_;
    return () unless defined $fx && defined $fy;
    my ($name,$lx,$ly);
    if ($fy<0.5){ $ly=$fy/0.5; ($name,$lx)= $fx<0.5 ? ('axial',$fx/0.5) : ('coronal',($fx-0.5)/0.5); }
    else       { return () if $fx>=0.5; $name='sagittal'; $lx=$fx/0.5; $ly=($fy-0.5)/0.5; }
    my $pc=panel_cache($name,$ZOOM); my $spec=plane_spec($name);
    my @v=rc_to_vox($spec, $pc->{r0}+($ly*$pc->{Hs})*$pc->{step}, $pc->{c0}+($lx*$pc->{Ws})*$pc->{step});
    my @mni=apply_affine($mri2mni,@v);
    (\@mni,\@v,"$name(approx)");
}

# ------------------------------------------------------------------ render callback (2x2)
sub build_window { my ($state)=@_;
    my $s0=defined $state->{0}?$state->{0}:$s0_init;
    my $s1=defined $state->{1}?$state->{1}:$s1_init;
    my $level=$vmin+$s0*$range; my $width=$s1*$range; $width=1 if $width<1;
    ($level-$width/2, $level+$width/2) }
sub cur_zoom { my ($state)=@_; my $s2=defined $state->{2}?$state->{2}:$s2_init; zoom_from_s2($s2) }

sub apply_pick_fxfy {
    my ($fx,$fy)=@_;
    return 0 unless $LAST_FIG;
    my ($idx,$xd,$yd)=(-1,undef,undef); my $i=0;
    for my $ax (@{ $LAST_FIG->axes_list }) {
        last if $i>2;                       # メイン3面だけで解決(オーバービュー/ヒストは無視)
        my ($x,$y)=eval { $ax->image_frac_to_data($fx,$fy) };
        if (defined $x) { ($idx,$xd,$yd)=($i,$x,$y); last }
        $i++;
    }
    warn sprintf("pick resolve: fx=%.3f fy=%.3f -> ax_idx=%d\n",$fx,$fy,$idx) if $o{debug};
    my @res=pick_to_mni($idx,$xd,$yd);
    @res=pick_fallback($fx,$fy) unless @res && ref $res[0];
    if (@res && ref $res[0]) {
        set_center_mni(@{$res[0]});
        printf STDERR "  --eec %.1f,%.1f,%.1f   (vox %.0f,%.0f,%.0f, %s)\n",
            @{$res[0]}, @{$res[1]}, $res[2];
        return 1;
    }
    return 0;
}

sub render {
    my ($state,$w,$h)=@_;
    my ($w_lo,$w_hi)=build_window($state);
    my $zoom=cur_zoom($state); $ZOOM=$zoom;
    my $key=join(",",@cidx)."|".sprintf("%.1f,%.1f,%.1f",@center_mni)
           ."|".sprintf("%.0f,%.0f",$w_lo,$w_hi)."|z".sprintf("%.4f",$zoom);
    return $LAST_FIG if $LAST_FIG && defined $BUILT_KEY && $BUILT_KEY eq $key;

    require PDL::Graphics::Cairo; PDL::Graphics::Cairo->import(qw(subplots));
    my ($fig,@rows)=subplots(2,2, width=>($w||1100), height=>($h||1000));
    draw_panel($rows[0][0],'axial',   $w_lo,$w_hi,$zoom);
    draw_panel($rows[0][1],'coronal', $w_lo,$w_hi,$zoom);
    draw_panel($rows[1][0],'sagittal',$w_lo,$w_hi,$zoom);

    my $ax=$rows[1][1];
    if ($zoom>1.01) {
        draw_overview($ax,$w_lo,$w_hi,$zoom);
    } else {
        $ax->line($HX,$HY,color=>[0.7,0.7,0.7],lw=>1.2);
        $ax->axvspan($w_lo,$w_hi,color=>[0.2,0.7,1],alpha=>0.20) if $ax->can('axvspan');
        $ax->axvline($w_lo,color=>[0.2,0.7,1],lw=>1.2); $ax->axvline($w_hi,color=>[0.2,0.7,1],lw=>1.2);
        $ax->set_title(sprintf("hist + window [%.0f, %.0f]",$w_lo,$w_hi));
        $ax->set_xlabel("intensity"); $ax->set_ylabel("log count");
    }

    if ($fig->can('suptitle')) {
        my @cv=map{int($_+0.5)} apply_affine($mni2mri,@center_mni);
        my $ci=$vol->at(map{ my $x=$cv[$_]; $x=0 if $x<0; my @N=($Ni,$Nj,$Nk); $x=$N[$_]-1 if $x>$N[$_]-1; $x } (0,1,2));
        $fig->suptitle(sprintf("center MNI(%.0f,%.0f,%.0f)  vox(%d,%d,%d)  I=%.0f  zoom x%.1f   [click=move  ·  sliders: level/width]",
            @center_mni,@cv,$ci,$zoom));
    }
    $fig->tight_layout if $fig->can('tight_layout');
    $LAST_FIG=$fig; $BUILT_KEY=$key;
    return $fig;
}

# ------------------------------------------------------------------ probe(対話なしで pick 経路を検証)
if (defined $o{probe}) {
    my ($pname,$fx,$fy)=split /\s*,\s*/,$o{probe};
    if ($pname eq 'canvas') {
        my $f0=render({}, 1100,1000); $f0->save("/tmp/_probe_warm.png");
        apply_pick_fxfy($fx,$fy);
        my $fig=render({}, 1100,1000);
        $fig->save($o{out}); print STDERR "probe-canvas: wrote $o{out}\n"; exit 0;
    }
    my %idx=(axial=>0,coronal=>1,sagittal=>2); die "bad panel\n" unless defined $idx{$pname};
    $ZOOM=$o{zoom};
    my $pc=panel_cache($pname,$ZOOM);
    my @res=pick_to_mni($idx{$pname}, $fx*$pc->{Ws}, $pc->{Hs}-$fy*$pc->{Hs});
    if (@res && ref $res[0]) { set_center_mni(@{$res[0]});
        printf STDERR "  --eec %.1f,%.1f,%.1f   (vox %.0f,%.0f,%.0f, %s)\n",@{$res[0]},@{$res[1]},$res[2]; }
    my $fig=render({}, 1100,1000);
    $fig->save($o{out}); print STDERR "probe: wrote $o{out}\n"; exit 0;
}

# ------------------------------------------------------------------ interactive (Mac / giza)
require PDL::Graphics::Cairo;
require PDL::Graphics::Cairo::Driver::GS;
my $drv = PDL::Graphics::Cairo::Driver::GS->new(width=>1100, height=>1000, title=>"NYHead ortho");
$drv->show_interactive(
    render         => \&render,
    cursor_overlay => 1,
    init           => { 0=>$s0_init, 1=>$s1_init },   # 下=level, 右=width の2本(ズームは --zoom)
    on_pick   => sub { my ($fx,$fy,$btn)=@_;
                       warn sprintf("PICK   fx=%.3f fy=%.3f btn=%d\n",$fx,$fy,$btn) if $o{debug};
                       apply_pick_fxfy($fx,$fy); 1 },
    on_cursor => sub { 1 },
);
