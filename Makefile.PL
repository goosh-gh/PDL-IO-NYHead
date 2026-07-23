use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME             => 'PDL::IO::NYHead',
    VERSION_FROM     => 'lib/PDL/IO/NYHead.pm',
    ABSTRACT         => 'Read the New York Head (sa_nyhead.mat, MATLAB v7.3/HDF5) into PDL',
    AUTHOR           => 'goosh',
    LICENSE          => 'perl_5',
    MIN_PERL_VERSION => '5.016',
    PREREQ_PM        => {
        'PDL'           => '2.000',
        'PDL::IO::HDF5' => '0',
        'Carp'          => '0',
    },
    TEST_REQUIRES    => {
        'Test::More' => '0.88',
    },
    META_MERGE => {
        'meta-spec' => { version => 2 },
        resources   => {
            repository => {
                type => 'git',
                url  => 'https://github.com/goosh-gh/PDL-IO-NYHead.git',
                web  => 'https://github.com/goosh-gh/PDL-IO-NYHead',
            },
        },
        # h5dump は CPAN 依存ではなく実行時の外部ツール(HDF5 tools)。
        # 電極名/atlas 名(cell-of-strings)の読み取りにのみ必要。
        prereqs => {
            runtime => {
                requires   => {},
                recommends => {},
            },
        },
    },
    dist  => { COMPRESS => 'gzip -9f', SUFFIX => 'gz' },
    clean => { FILES => 'PDL-IO-NYHead-*' },
);
