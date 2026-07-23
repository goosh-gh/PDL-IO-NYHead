use strict;
use warnings;
use Test::More tests => 3;

use_ok('PDL::IO::NYHead');
use_ok('PDL::IO::NYHead::H5Dump');
ok(defined $PDL::IO::NYHead::VERSION, "PDL::IO::NYHead sets \$VERSION ($PDL::IO::NYHead::VERSION)");
