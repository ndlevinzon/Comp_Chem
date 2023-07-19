#!/usr/bin/perl

use strict;
use warnings;

# Get a list of all .tar.gz files in the current directory
my @tarballs = glob('*.tar.gz');

# Loop through each tarball and extract its contents
foreach my $tarball (@tarballs) {
    # Run the 'tar xf' command for each tarball
    system("tar xf $tarball");
}
