#!/usr/bin/perl

use strict;
use warnings;

# Get the current directory
my $current_directory = `pwd`;
chomp($current_directory);

# Get a list of all .tar.gz files in the current directory
my @tarballs = glob('*.tar.gz');

# Loop through each tarball and extract its contents
foreach my $tarball (@tarballs) {
    # Run the 'tar xzvf' command for each tarball
    system("tar xf $tarball");
}

# Get the parent directory (repack_directory is located in the directory directly above the current directory)
my $repack_directory = `cd .. && pwd`;
chomp($repack_directory);
$repack_directory .= '/main/repack';

# After extracting, run the 'make_tarballs.bash' script with the specified arguments
my $make_tarballs_script = 'make_tarballs.bash';
system("bash $make_tarballs_script $current_directory $repack_directory");
