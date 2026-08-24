#!/usr/bin/perl
# nyhead_ortho_gs.pl — giza 対話 2x2 ortho ビューア。
#   左上 axial / 右上 coronal / 左下 sagittal / 右下 ヒストグラム(+窓+読み出し)。
#   * クリック → その断面の点へ全ヘアライン移動 + 端末に `--eec X,Y,Z` 出力(コピペ用)
#   * スライダ2本: id0=level(窓中心), id1=width(窓幅) で grayscale コントラスト調整
#   * 内部空気(外耳道等)をシアン重畳、電極を各断面に表示
# giza-server 再ビルド不要(PICK/CURSOR/SLIDER は既存配線)。実機は Mac。
#
#   perl nyhead_ortho_gs.pl --mat sa_nyhead.mat --center -65,-25,-60 --eec -65,-25,-60 --nyhead Ex19,LPA
#   perl nyhead_ortho_gs.pl --self-test --probe coronal,0.60,0.52   # 対話なしで pick 経路を検証

use strict; use warnings; use PDL; use PDL::Image2D; use Getopt::Long;
use lib '/Users/goosh/src/PDL_IO_NYHead/lib/';

# ------------------------------------------------------------------ options
my %o = (
    mat=>undef, selftest=>0, center=>undef, eec=>undef, elec=>[], nyhead=>undef,
    overlay=>'air', air_side=>'auto', air_color=>"0,0.9,1", bg=>'auto',
    radiological=>0, max_dim=>260, elec_orig=>1,
    probe=>undef, out=>"nyhead_ortho.png", debug=>1,
);
GetOptions(
    "mat=s"=>\$o{mat}, "self-test!"=>\$o{selftest}, "center=s"=>\$o{center},
    "eec=s"=>\$o{eec}, "elec=s@"=>$o{elec}, "nyhead=s"=>\$o{nyhead},
    "overlay=s"=>\$o{overlay}, "air-side=s"=>\$o{air_side}, "air-color=s"=>\$o{air_color},
    "bg=s"=>\$o{bg}, "radiological!"=>\$o{radiological}, "max-dim=i"=>\$o{max_dim},
    "elec-orig!"=>\$o{elec_orig}, "probe=s"=>\$o{probe}, "out=s"=>\$o{out}, "debug!"=>\$o{debug},
) or die "bad options\n";
die "give --mat or --self-test\n" unless $o{mat} || $o{selftest};

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

# ヒストグラム(頭部ボクセル分布)は不変 → 起動時に1回だけ計算(毎フレームの qsort を回避)
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

# slice/air キャッシュ(center 変更時のみ再計算)
my %CACHE;   # plane => {cidx, S, Hs, Ws, step, ext}
sub panel_cache {
    my ($name)=@_; my $spec=plane_spec($name);
    my $c=$CACHE{$name};
    return $c if $c && $c->{cidx}==$spec->{cidx};
    my $img=plane_slice($spec); my ($Hf,$Wf)=$img->dims;
    my $step=int((($Hf>$Wf?$Hf:$Wf)+$o{max_dim}-1)/$o{max_dim}); $step=1 if $step<1;
    my $S=($step>1)?$img->slice("0:-1:$step,0:-1:$step")->sever:$img;
    my ($Hs,$Ws)=$S->dims;
    my $ext=zeroes(byte,$Hs,$Ws);
    ($ext,undef)=air_components($S,$bg_mode eq 'high',$bg_thr) if $bg_mode ne 'none';
    my $intair=zeroes(byte,$Hs,$Ws);
    (undef,$intair)=air_components($S,$air_high,$air_thr) if $o{overlay} eq 'air';
    $CACHE{$name}={cidx=>$spec->{cidx},S=>$S,Hs=>$Hs,Ws=>$Ws,step=>$step,ext=>$ext,intair=>$intair,Hf=>$Hf,Wf=>$Wf};
}

# ------------------------------------------------------------------ draw one panel into $ax
sub draw_panel {
    my ($ax,$name,$w_lo,$w_hi)=@_;
    my $spec=plane_spec($name); my $pc=panel_cache($name);
    my ($S,$Hs,$Ws,$step,$ext)=@{$pc}{qw(S Hs Ws step ext)};
    my $gn=(($S-$w_lo)/(($w_hi-$w_lo)||1))->clip(0,1);
    my $g=$gn*(1-$ext);
    my ($a,@ocol)=(zeroes(float,$Hs,$Ws),(0,0,0));
    if ($o{overlay} eq 'air'){ $a=$pc->{intair}*0.55; @ocol=@acol; }
    my $rgb=zeroes(float,$Hs,$Ws,3);
    $rgb->slice(":,:,($_)").= $g*(1-$a)+$ocol[$_]*$a for (0,1,2);
    $ax->imshow($rgb,origin=>'upper'); $ax->axis('off'); $ax->set_aspect('equal'); $ax->set_title($spec->{title});
    my $to_xy=sub{ my($row,$col)=@_; ($col/$step,$Hs-$row/$step) };
    my @draw;
    for my $e (@elecs){ my @ev=apply_affine($mni2mri,$e->{x},$e->{y},$e->{z});
        my ($row,$col,$off)=map_vox($spec,@ev);
        next if $col<-2||$col>$pc->{Wf}+2||$row<-2||$row>$pc->{Hf}+2;
        my ($ex,$ey)=$to_xy->($row,$col); my $off_mm=abs($off)*$vsize[$spec->{fixed}];
        my %cof=(eec=>[1,0,0],nyhead=>[0,0.8,1],user=>[1,0.85,0]);
        push @draw,{ex=>$ex,ey=>$ey,c=>($cof{$e->{type}}//[1,0.85,0]),on=>($off_mm<=4),
                    tag=>($off_mm<=4?$e->{label}:sprintf("%s %+.0fmm",$e->{label},$off>0?$off_mm:-$off_mm))} }
    my ($crow,$ccol)=map_vox($spec,apply_affine($mni2mri,@center_mni));
    my ($cx,$cy)=$to_xy->($crow,$ccol);
    $ax->axvline($cx,color=>[0.2,1,0.4],lw=>0.8,ls=>'dashed');
    $ax->axhline($cy,color=>[0.2,1,0.4],lw=>0.8,ls=>'dashed');
    for my $d (@draw){ $ax->scatter(pdl($d->{ex}),pdl($d->{ey}),s=>7,color=>$d->{c},marker=>'o',
                        alpha=>($d->{on}?0.95:0.45),fillstyle=>($d->{on}?'full':'none')) }
    for my $d (@draw){ $ax->text($d->{ex}+$Ws*0.015,$d->{ey},$d->{tag},color=>$d->{c},fontsize=>9,ha=>'left',va=>'middle') }
    my $L=$spec->{labels};
    $ax->text($Ws*0.5,$Hs*0.98,$L->{top},color=>'white',fontsize=>12,ha=>'center',va=>'top');
    $ax->text($Ws*0.5,$Hs*0.02,$L->{bottom},color=>'white',fontsize=>12,ha=>'center',va=>'bottom');
    $ax->text($Ws*0.02,$Hs*0.5,$L->{left},color=>'white',fontsize=>12,ha=>'left',va=>'middle');
    $ax->text($Ws*0.98,$Hs*0.5,$L->{right},color=>'white',fontsize=>12,ha=>'right',va=>'middle');
}

# ------------------------------------------------------------------ pick 座標 → center 更新
my @PLANE_BY_IDX=qw(axial coronal sagittal);   # axes_list の並び(2x2 行優先)
my $LAST_FIG;                                  # 直近に描いた Figure(軸箱で pick を正確解決)
my $BUILT_KEY;                                 # 直近に構築した (center|window) キー
sub pick_to_mni {
    my ($ax_idx,$cur_x,$cur_y)=@_;
    return () unless defined $ax_idx && $ax_idx>=0 && $ax_idx<=2;
    return () unless defined $cur_x && defined $cur_y;
    my $name=$PLANE_BY_IDX[$ax_idx]; my $pc=panel_cache($name); my $spec=plane_spec($name);
    my $col_s=$cur_x; my $row_s=$pc->{Hs}-$cur_y;         # imshow origin upper の y 反転を戻す
    my @v=rc_to_vox($spec,$row_s*$pc->{step},$col_s*$pc->{step});
    my @mni=apply_affine($mri2mni,@v);
    (\@mni,\@v,$name);
}
# _cursor_ax_idx が取れないとき用: 生画像分率(fx,fy)を 2x2 象限で割ってパネル決定(近似)
sub pick_fallback {
    my ($fx,$fy)=@_;
    return () unless defined $fx && defined $fy;
    my ($name,$lx,$ly);
    if ($fy<0.5){ $ly=$fy/0.5; ($name,$lx)= $fx<0.5 ? ('axial',$fx/0.5) : ('coronal',($fx-0.5)/0.5); }
    else       { return () if $fx>=0.5; $name='sagittal'; $lx=$fx/0.5; $ly=($fy-0.5)/0.5; }
    my $pc=panel_cache($name); my $spec=plane_spec($name);
    my @v=rc_to_vox($spec, ($ly*$pc->{Hs})*$pc->{step}, ($lx*$pc->{Ws})*$pc->{step});
    my @mni=apply_affine($mri2mni,@v);
    (\@mni,\@v,"$name(approx)");
}

# ------------------------------------------------------------------ render callback (2x2)
sub build_window { my ($state)=@_;
    my $s0=defined $state->{0}?$state->{0}:$s0_init;
    my $s1=defined $state->{1}?$state->{1}:$s1_init;
    my $level=$vmin+$s0*$range; my $width=$s1*$range; $width=1 if $width<1;
    ($level-$width/2, $level+$width/2) }

# クリックの生 (fx,fy) を、直近に描いた図の軸箱で解決 → center 更新 + --eec 出力。
# on_pick から即座に呼ぶ(再描画より前に center を更新するので1フレーム遅れが出ない)。
sub apply_pick_fxfy {
    my ($fx,$fy)=@_;
    return 0 unless $LAST_FIG;
    my ($idx,$xd,$yd)=(-1,undef,undef); my $i=0;
    for my $ax (@{ $LAST_FIG->axes_list }) {
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
    # center も window も変わっていなければ(hover 等)再構築しない
    my $key=join(",",@cidx)."|".sprintf("%.0f,%.0f",$w_lo,$w_hi);
    return $LAST_FIG if $LAST_FIG && defined $BUILT_KEY && $BUILT_KEY eq $key;

    require PDL::Graphics::Cairo;
    PDL::Graphics::Cairo->import(qw(subplots));
    my ($fig,@rows)=subplots(2,2, width=>($w||1100), height=>($h||1000));
    draw_panel($rows[0][0],'axial',   $w_lo,$w_hi);
    draw_panel($rows[0][1],'coronal', $w_lo,$w_hi);
    draw_panel($rows[1][0],'sagittal',$w_lo,$w_hi);

    # 右下: ヒストグラム(事前計算済) + 窓
    my $ax=$rows[1][1];
    $ax->line($HX,$HY,color=>[0.7,0.7,0.7],lw=>1.2);
    $ax->axvspan($w_lo,$w_hi,color=>[0.2,0.7,1],alpha=>0.20) if $ax->can('axvspan');
    $ax->axvline($w_lo,color=>[0.2,0.7,1],lw=>1.2); $ax->axvline($w_hi,color=>[0.2,0.7,1],lw=>1.2);
    $ax->set_title(sprintf("hist + window [%.0f, %.0f]",$w_lo,$w_hi));
    $ax->set_xlabel("intensity"); $ax->set_ylabel("log count");

    if ($fig->can('suptitle')) {
        my @cv=map{int($_+0.5)} apply_affine($mni2mri,@center_mni);
        my $ci=$vol->at(map{ my $x=$cv[$_]; $x=0 if $x<0; my @N=($Ni,$Nj,$Nk); $x=$N[$_]-1 if $x>$N[$_]-1; $x } (0,1,2));
        $fig->suptitle(sprintf("center MNI(%.0f,%.0f,%.0f)  vox(%d,%d,%d)  I=%.0f   [click=move crosshair, sliders=level/width]",
            @center_mni,@cv,$ci));
    }
    $fig->tight_layout if $fig->can('tight_layout');   # レイアウト確定(+ image_frac_to_data 用の枠)
    $LAST_FIG=$fig; $BUILT_KEY=$key;
    return $fig;
}

# ------------------------------------------------------------------ probe(対話なしで pick 経路を検証)
if (defined $o{probe}) {
    my ($pname,$fx,$fy)=split /\s*,\s*/,$o{probe};
    if ($pname eq 'canvas') {
        # path A(実対話と同経路): 初期フレームを描画 → 生 (fx,fy) を軸箱で解決 → 再描画
        my $f0=render({}, 1100,1000); $f0->save("/tmp/_probe_warm.png");   # draw → 枠が付く(driver 相当)
        apply_pick_fxfy($fx,$fy);
        my $fig=render({}, 1100,1000);
        $fig->save($o{out}); print STDERR "probe-canvas: wrote $o{out}\n"; exit 0;
    }
    # パネル指定(ax_idx 直接): 座標を直接 MNI 化して center 更新
    my %idx=(axial=>0,coronal=>1,sagittal=>2); die "bad panel\n" unless defined $idx{$pname};
    my $pc=panel_cache($pname);
    my @res=pick_to_mni($idx{$pname}, $fx*$pc->{Ws}, $pc->{Hs}-$fy*$pc->{Hs});
    if (@res && ref $res[0]) { set_center_mni(@{$res[0]});
        printf STDERR "  --eec %.1f,%.1f,%.1f   (vox %.0f,%.0f,%.0f, %s)\n",@{$res[0]},@{$res[1]},$res[2]; }
    my $fig=render({}, 1100,1000);
    $fig->save($o{out});
    print STDERR "probe: wrote $o{out}\n";
    exit 0;
}

# ------------------------------------------------------------------ interactive (Mac / giza)
require PDL::Graphics::Cairo;
require PDL::Graphics::Cairo::Driver::GS;
my $drv = PDL::Graphics::Cairo::Driver::GS->new(width=>1100, height=>1000, title=>"NYHead ortho");
$drv->show_interactive(
    render         => \&render,
    cursor_overlay => 1,                       # PICK で自動再描画(+ _cursor_ax_idx/_x/_y を供給)
    init           => { 0=>$s0_init, 1=>$s1_init },
    on_pick   => sub { my ($fx,$fy,$btn)=@_;
                       warn sprintf("PICK   fx=%.3f fy=%.3f btn=%d\n",$fx,$fy,$btn) if $o{debug};
                       apply_pick_fxfy($fx,$fy); 1 },
    on_cursor => sub { 1 },   # hover は移動しない(再描画は cursor_overlay が担当)
);
