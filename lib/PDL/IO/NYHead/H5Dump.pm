package PDL::IO::NYHead::H5Dump;
use strict; use warnings; use Carp;
our $VERSION = '0.04';
sub new { my ($c,%a)=@_; bless { h5dump=>$a{h5dump}||'h5dump' }, $c }
sub read_cell_strings {
    my ($self,$file,$dataset)=@_;
    -r $file or croak "cannot read '$file'";
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
1;
