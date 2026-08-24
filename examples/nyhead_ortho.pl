#!/usr/bin/perl
# nyhead_ortho.pl — NYHead MRI を 2x2 の ortho ビューで描く。
#   左上 axial / 右上 coronal / 左下 sagittal / 右下 ヒストグラム(+現在窓)。
# grayscale 窓は min/max から始めて auto-contrast(エントロピー最大の窓)を探せる。
# 内部空気(外耳道・副鼻腔)をシアンで重畳可。クリック点→MNI の逆写像を持ち、
# --click で座標出力を検証する(将来 giza 対話ビューアの pick 配線がこれを呼ぶ)。
#
#   perl nyhead_ortho.pl --mat sa_nyhead.mat --center -65,-25,-60 --eec -65,-25,-60
#   perl nyhead_ortho.pl --self-test --out ortho.png
#   perl nyhead_ortho.pl --self-test --click coronal,0.62,0.55   # クリック模擬→--eec 出力
#
# 逆写像(クリック→MNI)は giza の PICK(fx,fy=画像比率)をそのまま食える形にしてある。

use strict; use warnings; use PDL;
use PDL::Image2D;   # cc8compt
use Getopt::Long;
use lib '/Users/goosh/src/PDL_IO_NYHead/lib/';

# ------------------------------------------------------------------ options
my %o = (
    mat => undef, selftest => 0,
    center => undef, eec => undef, elec => [], nyhead => undef,
    window => undef,                 # "lo,hi" 明示窓
    window_mode => 'auto',           # auto|minmax|manual
    overlay => 'air',                # air|none
    air_side => 'auto', air_color => "0,0.9,1",
    bg => 'auto',
    radiological => 0,
    max_dim => 300,
    elec_orig => 1,
    click => undef,                  # "PANEL,fx,fy" (fx,fy=0..1 画像比率, 左上原点)
    out => "nyhead_ortho.png",
);
GetOptions(
    "mat=s"=>\$o{mat}, "self-test!"=>\$o{selftest},
    "center=s"=>\$o{center}, "eec=s"=>\$o{eec}, "elec=s@"=>$o{elec}, "nyhead=s"=>\$o{nyhead},
    "window=s"=>\$o{window}, "window-mode=s"=>\$o{window_mode},
    "overlay=s"=>\$o{overlay}, "air-side=s"=>\$o{air_side}, "air-color=s"=>\$o{air_color},
    "bg=s"=>\$o{bg}, "radiological!"=>\$o{radiological}, "max-dim=i"=>\$o{max_dim},
    "elec-orig!"=>\$o{elec_orig}, "click=s"=>\$o{click}, "out=s"=>\$o{out},
) or die "bad options\n";
die "give --mat FILE or --self-test\n" unless $o{mat} || $o{selftest};

# ------------------------------------------------------------------ helpers
sub parse_xyz { my @v=split /\s*,\s*/,$_[0]; die "want x,y,z got '$_[0]'\n" unless @v==3; map {0+$_} @v }
sub pctl { my($p,$q)=@_; my $s=$p->flat->qsort; my $n=$s->nelem; return 0 unless $n;
           my $i=int($q*($n-1)+0.5); $i=0 if $i<0; $i=$n-1 if $i>$n-1; $s->at($i) }
sub apply_affine { my($M,$x,$y,$z)=@_; map { $M->at($_,0)*$x+$M->at($_,1)*$y+$M->at($_,2)*$z+$M->at($_,3) } (0,1,2) }
sub to_math_affine {
    my ($A)=@_; my ($best,$bd)=(undef,1e30);
    for my $M ($A->sever, $A->transpose->sever) {
        my $b=abs($M->at(3,0))+abs($M->at(3,1))+abs($M->at(3,2)); my $one=abs($M->at(3,3)-1);
        return $M if $b<1e-3 && $one<1e-3;
        ($best,$bd)=($M,$b+$one) if $b+$one<$bd;
    }
    warn "to_math_affine: no clean bottom row; using closest\n"; $best;
}
sub h5get { my($h5,$path)=@_; my @p=grep{length}split m{/},$path; my $ds=pop @p;
            my $n=$h5; $n=$n->group($_) for @p; $n->dataset($ds)->get }

sub air_components {
    my ($S,$is_high,$thr)=@_;
    my $air=($is_high?($S>$thr):($S<$thr))->byte;
    return (zeroes(byte,$S->dims),zeroes(byte,$S->dims)) if $air->sum==0;
    my $lab=cc8compt($air); my ($H,$W)=$lab->dims;
    my $border=$lab->slice("0,:")->flat->append($lab->slice("-1,:")->flat)
              ->append($lab->slice(":,0")->flat)->append($lab->slice(":,-1")->flat);
    my $ext=zeroes(byte,$H,$W);
    for my $L (grep{$_>0} $border->uniq->list){ $ext=$ext|($lab==$L) }
    ($ext->byte, ($air&($ext==0))->byte);
}

# auto-contrast: 頭部ボクセルの分布で、窓内に写した表示ヒストの entropy を最大化する (lo,hi)
sub auto_window {
    my ($vals)=@_;
    my @los=map{pctl($vals,$_)} (0.005,0.01,0.02,0.05);
    my @his=map{pctl($vals,$_)} (0.95,0.98,0.99,0.995);
    my ($best,$bestE)=([$los[2],$his[1]],-1);
    for my $lo (@los){ for my $hi (@his){
        next unless $hi>$lo;
        my $g=(($vals-$lo)/($hi-$lo))->clip(0,1);
        my $h=histogram($g,1/64,0,64)->double;
        my $p=$h/($h->sum||1);
        my $E=-($p*($p+1e-12)->log)->sum;
        ($best,$bestE)=([$lo,$hi],$E) if $E>$bestE;
    }}
    @$best;
}

# ------------------------------------------------------------------ load
my ($vol,$mri2mni,$mni2mri);
if ($o{selftest}) {
    my ($Ni,$Nj,$Nk)=(80,90,70); my @ctr=(40,45,35);
    my $ii=xvals($Ni,$Nj,$Nk); my $jj=yvals($Ni,$Nj,$Nk); my $kk=zvals($Ni,$Nj,$Nk);
    my $rr=sqrt((($ii-$ctr[0])/30)**2+(($jj-$ctr[1])/34)**2+(($kk-$ctr[2])/26)**2);
    my $brain=($rr<0.86); my $skull=(($rr>=0.86)&($rr<1.0)); my $scalp=(($rr>=1.0)&($rr<1.08));
    $vol=$brain*800+$skull*150+$scalp*600 + $brain*(100*sin($rr*12));
    for my $cx (14,66){ my $da=sqrt(($ii-$cx)**2+($jj-45)**2+($kk-32)**2); $vol=$vol*($da>=4) }
    $mri2mni=zeroes(double,4,4); $mri2mni->set(3,3,1);
    for my $r (0..2){ $mri2mni->set($r,$r,0.5); $mri2mni->set($r,3,-0.5*$ctr[$r]) }
    $mni2mri=zeroes(double,4,4); $mni2mri->set(3,3,1);
    for my $r (0..2){ $mni2mri->set($r,$r,2); $mni2mri->set($r,3,$ctr[$r]) }
    $o{eec}    //= "-14,0,-2";
    $o{center} //= "0,0,0";
} else {
    print STDERR "reading $o{mat} ...\n";
    require PDL::IO::HDF5;
    my $h5=PDL::IO::HDF5->new($o{mat}) or die "cannot open $o{mat}\n";
    $vol=h5get($h5,'/sa/mri/data')->double;
    $mri2mni=to_math_affine(h5get($h5,'/sa/mri2mni')->double);
    $mni2mri=to_math_affine(h5get($h5,'/sa/mni2mri')->double);
    printf STDERR "vol dims = %s\n", join("x",$vol->dims);
}
my ($Ni,$Nj,$Nk)=$vol->dims;
my @vsize=map{ my $c=$_; sqrt($mri2mni->at(0,$c)**2+$mri2mni->at(1,$c)**2+$mri2mni->at(2,$c)**2) } (0,1,2);

# ------------------------------------------------------------------ background / window / air side
my ($vmin,$vmax,$range)=($vol->min,$vol->max,($vol->max-$vol->min)||1);
my $bg_mode=$o{bg};
if ($bg_mode eq 'auto'){
    my $lo=($vol<($vmin+0.01*$range))->sum; my $hi=($vol>($vmax-0.01*$range))->sum;
    $bg_mode=($hi>$lo)?'high':'low';
}
my $bg_thr=($bg_mode eq 'high')?$vmax-0.03*$range:($bg_mode eq 'low')?$vmin+0.03*$range:undef;
my $air_high=($o{air_side} eq 'high')?1:($o{air_side} eq 'low')?0:($bg_mode eq 'high'?1:0);
my $air_thr=$air_high?$vmax-0.03*$range:$vmin+0.03*$range;

my $mask=($bg_mode eq 'high')?($vol<($vmax-0.05*$range))
        :($bg_mode eq 'low') ?($vol>($vmin+0.05*$range)):ones(byte,$vol->dims);
my $headvals=$vol->flat->where($mask->flat); $headvals=$vol->flat if $headvals->nelem<100;

my ($w_lo,$w_hi);
if (defined $o{window}) { ($w_lo,$w_hi)=map{0+$_} split/\s*,\s*/,$o{window}; }
elsif ($o{window_mode} eq 'minmax') { ($w_lo,$w_hi)=($vmin,$vmax); }
else { ($w_lo,$w_hi)=auto_window($headvals); }   # auto(既定)
printf STDERR "background=%s  window[%.0f,%.0f] (mode=%s)  head p2/p50/p98=%.0f/%.0f/%.0f  min/max=%.0f/%.0f\n",
    $bg_mode,$w_lo,$w_hi,$o{window_mode},pctl($headvals,0.02),pctl($headvals,0.5),pctl($headvals,0.98),$vmin,$vmax;

# ------------------------------------------------------------------ center + electrodes
my @center_mni = defined $o{center}?parse_xyz($o{center})
               : defined $o{eec}   ?parse_xyz($o{eec}) : (0,0,0);
my @cvox=apply_affine($mni2mri,@center_mni);
my @cidx=map{ my $v=int($_+0.5); $v } @cvox;
$cidx[0]=0 if $cidx[0]<0; $cidx[0]=$Ni-1 if $cidx[0]>$Ni-1;
$cidx[1]=0 if $cidx[1]<0; $cidx[1]=$Nj-1 if $cidx[1]>$Nj-1;
$cidx[2]=0 if $cidx[2]<0; $cidx[2]=$Nk-1 if $cidx[2]>$Nk-1;

my @elecs;
if (defined $o{eec}){ my @p=parse_xyz($o{eec}); push @elecs,{label=>"EEC",x=>$p[0],y=>$p[1],z=>$p[2],type=>"eec"} }
for my $spec (@{$o{elec}}){ my($l,$r)=split/=/,$spec,2; die "bad --elec '$spec'\n" unless defined $r;
    my @p=parse_xyz($r); push @elecs,{label=>$l,x=>$p[0],y=>$p[1],z=>$p[2],type=>"user"} }
if (defined $o{nyhead} && !$o{selftest}){
    my @want=split/\s*,\s*/,$o{nyhead};
    eval { require PDL::IO::HDF5; require PDL::IO::NYHead;
        my $ny=PDL::IO::NYHead->new($o{mat}); my $lab=$ny->electrode_labels;
        my $pos=$o{elec_orig}?$ny->electrode_pos
              : h5get(PDL::IO::HDF5->new($o{mat}),'/sa/locs_3D')->slice(':,0:2')->sever;
        my %idx; $idx{$lab->[$_]}=$_ for 0..$#$lab;
        for my $nm (@want){ exists $idx{$nm} or (warn("  --nyhead: '$nm' not found\n"),next);
            my $r=$idx{$nm}; push @elecs,{label=>$nm,x=>$pos->at($r,0),y=>$pos->at($r,1),z=>$pos->at($r,2),type=>"nyhead"} }
        1 } or warn "  --nyhead failed ($@); use --elec\n";
}

# ------------------------------------------------------------------ plane specs
my $lr=$o{radiological}?1:0;
my %PLANE=(
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

sub plane_slice {
    my ($spec)=@_; my @idx=(undef,undef,undef); $idx[$spec->{fixed}]=$spec->{cidx};
    my @ss=map{ defined $idx[$_]?"($idx[$_])":":" } (0,1,2);
    my $s2=$vol->slice(join(",",@ss))->squeeze;
    my @rem=grep{ $_!=$spec->{fixed} } (0,1,2);
    my ($cd)=grep{ $rem[$_]==$spec->{col_axis} } (0,1);
    my ($rd)=grep{ $rem[$_]==$spec->{row_axis} } (0,1);
    my $img=($rd==0)?$s2:$s2->transpose; $img=$img->sever;
    $img=$img->slice("-1:0,:")->sever if $spec->{row_flip};
    $img=$img->slice(":,-1:0")->sever if $spec->{col_flip};
    $img;
}
sub map_vox {
    my ($spec,$i,$j,$k)=@_; my @v=($i,$j,$k); my @N=($Ni,$Nj,$Nk);
    my $col=$v[$spec->{col_axis}]; $col=($N[$spec->{col_axis}]-1)-$col if $spec->{col_flip};
    my $row=$v[$spec->{row_axis}]; $row=($N[$spec->{row_axis}]-1)-$row if $spec->{row_flip};
    ($row,$col,$v[$spec->{fixed}]-$spec->{cidx});
}
# 逆写像: plane 上の full-res (row,col) → voxel(i,j,k)
sub rc_to_vox {
    my ($spec,$row,$col)=@_; my @N=($Ni,$Nj,$Nk); my @v; $v[$spec->{fixed}]=$spec->{cidx};
    $v[$spec->{col_axis}]=$spec->{col_flip}?($N[$spec->{col_axis}]-1-$col):$col;
    $v[$spec->{row_axis}]=$spec->{row_flip}?($N[$spec->{row_axis}]-1-$row):$row;
    @v;
}

# ------------------------------------------------------------------ click -> MNI (giza PICK fx,fy を食う形)
if (defined $o{click}) {
    my ($pname,$fx,$fy)=split /\s*,\s*/,$o{click};
    die "bad --click (want PANEL,fx,fy)\n" unless $PLANE{$pname} && defined $fy;
    my $spec=$PLANE{$pname};
    my $img=plane_slice($spec); my ($Hf,$Wf)=$img->dims;
    # fx=左→右(col比率), fy=上→下(row比率, 画像上端=row0)
    my $col=$fx*$Wf; my $row=$fy*$Hf;
    my @v=rc_to_vox($spec,$row,$col);
    my @mni=apply_affine($mri2mni,@v);
    my $val=$vol->at(map{ my $x=int($_+0.5); $x=0 if $x<0; $x } @v);
    printf "click %s (fx=%.3f,fy=%.3f) -> vox(%.0f,%.0f,%.0f)  I=%.0f\n", $pname,$fx,$fy,@v,$val;
    printf "  --eec %.1f,%.1f,%.1f\n", @mni;    # コピペ用
    exit 0;
}

# ------------------------------------------------------------------ render 2x2
require PDL::Graphics::Cairo;
PDL::Graphics::Cairo->import(qw(subplots));
my ($fig,@rows)=subplots(2,2,figsize=>[12,11]);
my %CELL=(axial=>[0,0], coronal=>[0,1], sagittal=>[1,0]);   # 右下[1][1]=ヒストグラム

for my $pname (qw(axial coronal sagittal)) {
    my $spec=$PLANE{$pname}; my ($r,$c)=@{$CELL{$pname}}; my $ax=$rows[$r][$c];
    my $img=plane_slice($spec); my ($Hf,$Wf)=$img->dims;
    my $step=int((($Hf>$Wf?$Hf:$Wf)+$o{max_dim}-1)/$o{max_dim}); $step=1 if $step<1;
    my $S=($step>1)?$img->slice("0:-1:$step,0:-1:$step")->sever:$img;
    my ($Hs,$Ws)=$S->dims;

    my $gn=(($S-$w_lo)/(($w_hi-$w_lo)||1))->clip(0,1);
    my $ext=zeroes(byte,$Hs,$Ws);
    ($ext,undef)=air_components($S,$bg_mode eq 'high',$bg_thr) if $bg_mode ne 'none';
    my $g=$gn*(1-$ext);
    my ($a,@ocol)=(zeroes(float,$Hs,$Ws),(0,0,0));
    if ($o{overlay} eq 'air'){ my(undef,$ia)=air_components($S,$air_high,$air_thr); $a=$ia*0.55; @ocol=parse_xyz($o{air_color}); }
    my $rgb=zeroes(float,$Hs,$Ws,3);
    $rgb->slice(":,:,($_)").= $g*(1-$a)+$ocol[$_]*$a for (0,1,2);

    $ax->imshow($rgb,origin=>'upper'); $ax->axis('off'); $ax->set_aspect('equal'); $ax->set_title($spec->{title});
    my $to_xy=sub{ my($row,$col)=@_; ($col/$step, $Hs-$row/$step) };

    my @draw;
    for my $e (@elecs){
        my @ev=apply_affine($mni2mri,$e->{x},$e->{y},$e->{z});
        my ($row,$col,$off)=map_vox($spec,@ev);
        next if $col<-2||$col>$Wf+2||$row<-2||$row>$Hf+2;
        my ($ex,$ey)=$to_xy->($row,$col);
        my $off_mm=abs($off)*$vsize[$spec->{fixed}];
        my %cof=(eec=>[1,0,0],nyhead=>[0,0.8,1],user=>[1,0.85,0]);
        push @draw,{ex=>$ex,ey=>$ey,c=>($cof{$e->{type}}//[1,0.85,0]),on=>($off_mm<=4),
                    tag=>($off_mm<=4?$e->{label}:sprintf("%s %+.0fmm",$e->{label},$off>0?$off_mm:-$off_mm))};
    }
    my ($crow,$ccol)=map_vox($spec,@cvox); my ($cx,$cy)=$to_xy->($crow,$ccol);
    $ax->axvline($cx,color=>[0.2,1,0.4],lw=>0.8,ls=>'dashed');
    $ax->axhline($cy,color=>[0.2,1,0.4],lw=>0.8,ls=>'dashed');
    for my $d (@draw){ $ax->scatter(pdl($d->{ex}),pdl($d->{ey}),s=>7,color=>$d->{c},marker=>'o',
                        alpha=>($d->{on}?0.95:0.45),fillstyle=>($d->{on}?'full':'none')) }
    for my $d (@draw){ $ax->text($d->{ex}+$Ws*0.015,$d->{ey},$d->{tag},color=>$d->{c},fontsize=>9,ha=>'left',va=>'middle') }
    my $L=$spec->{labels};
    $ax->text($Ws*0.5,$Hs*0.98,$L->{top},color=>'white',fontsize=>13,ha=>'center',va=>'top');
    $ax->text($Ws*0.5,$Hs*0.02,$L->{bottom},color=>'white',fontsize=>13,ha=>'center',va=>'bottom');
    $ax->text($Ws*0.02,$Hs*0.5,$L->{left},color=>'white',fontsize=>13,ha=>'left',va=>'middle');
    $ax->text($Ws*0.98,$Hs*0.5,$L->{right},color=>'white',fontsize=>13,ha=>'right',va=>'middle');
}

# 右下: ヒストグラム(頭部ボクセル) + 現在窓
{
    my $ax=$rows[1][1];
    my $lo=pctl($headvals,0.005); my $hi=pctl($headvals,0.995);
    my $nb=80; my $hstep=(($hi-$lo)||1)/$nb;
    my $h=histogram($headvals,$hstep,$lo,$nb)->double;
    my $x=$lo + $hstep*(sequence($nb)+0.5);
    my $hy=($h+1)->log;                              # log スケールで裾も見えるように
    $ax->line($x,$hy,color=>[0.7,0.7,0.7],lw=>1.2);
    $ax->axvspan($w_lo,$w_hi,color=>[0.2,0.7,1],alpha=>0.20) if $ax->can('axvspan');
    $ax->axvline($w_lo,color=>[0.2,0.7,1],lw=>1.2);
    $ax->axvline($w_hi,color=>[0.2,0.7,1],lw=>1.2);
    $ax->set_title(sprintf("intensity hist + window [%.0f, %.0f]",$w_lo,$w_hi));
    $ax->set_xlabel("intensity"); $ax->set_ylabel("log count");
}

$fig->tight_layout;
$fig->save($o{out});
print STDERR "wrote $o{out}\n";
