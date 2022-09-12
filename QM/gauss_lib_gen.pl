#!/usr/bin/perl
#gauss_lib_geb.pl: Generate AMBER FF Parameters from Ligands in Seperate .PDB Files
use warnings;
use 5.010;
use Module::Load;
use File::Basename;

# Modules
require "$ENV{\"LMOD_PKG\"}/init/perl";
module("load openbabel/2.4.1"); # OpenBabel
module("load gaussian16/SSE4.C01"); # Gaussian16

# Set Working Directory Here
$home = '/uufs/chpc.utah.edu/common/home/cheatham-group7/nate/barrios22/SHP_1/QM';
opendir(DIR, $home) or die "Could not open $home\n";

# Make Directory For Each File
foreach my $filename (<*.pdb>) 
{
  	(my $without_extension = $filename) =~ s/\.[^.]+$//;
	if (! -e '$without_extension')
	{
  		system("mkdir $without_extension");
  		print "Making Directory: $without_extension\n";
  		system("mv $filename $without_extension");
  		print "Moving File: $filename\n";
	}
}

foreach my $newdirectory ( <$home/*> ) 
{
	if (-d $newdirectory)
	{
		chdir( $newdirectory ) or die "Couldn't go inside $newdirectory directory, $!";
		print "Current Working Directory is $newdirectory\n";
		foreach my $pdbfile ( glob  "$newdirectory/*.pdb" ) 
		{
			($file,$dir,$ext) = fileparse($pdbfile, qr/\.[^.]*/);
			@extensions = (".pdb", ".com", "_opt.com", "_HF631Gs_charge_nosym.com", "_opt.chk", "_HF631Gs_charge_nosym.chk");
			my $pdb    = $file.$extensions[0];
			my $com    = $file.$extensions[1];
			my $opt    = $file.$extensions[2];
			my $charge = $file.$extensions[3];

			# Convert .PDB to .COM
			if (! -f $com)
			{
				system("obabel -ipdb $pdb -ocom $com -m");
			}
			# Create GeomOpt Gaussian Input
			open my $in,  '<', $com  or die "Can't read COM file: $!";
			open my $out, '>', $opt  or die "Can't write OPT file: $!";
			print $out "%chk=$file$extensions[4]\n#P B3LYP/6-31G* Opt Freq=NoRaman\n$opt\n\n";
			while( <$in> )
    			{
				next if $. <= 4;
   				print $out $_;
   			}
			close $out;

			# Create Charge Gaussian Input
			open my $in,  '<', $com    or die "Can't read COM file: $!";
			open my $out, '>', $charge or die "Can't write CHARGE file: $!";
			print $out "%chk=$file$extensions[5]\n#P HF/6-31G* Pop=MK iop(6/33=2,6/41=10,6/42=6) nosym Test\n\n$charge\n\n";	
			while( <$in> )
    			{
				next if $. <= 4;
   				print $out $_;
   			}
			close $out;
			
			# Run Gaussian
			system("g16 $opt");
			system("g16 $charge");

		}
	}
	else
	{
		print "$newdirectory is not a directory, skipping...\n";	
	}
}

exit;

#      ::::    ::: :::::::::  :::            ::::::::   :::::::   ::::::::   :::::::: #
#     :+:+:   :+: :+:    :+: :+:           :+:    :+: :+:   :+: :+:    :+: :+:    :+: #
#    :+:+:+  +:+ +:+    +:+ +:+                 +:+  +:+   +:+       +:+        +:+   #
#   +#+ +:+ +#+ +#+    +:+ +#+               +#+    +#+   +:+     +#+        +#+      #
#  +#+  +#+#+# +#+    +#+ +#+             +#+      +#+   +#+   +#+        +#+         #
# #+#   #+#+# #+#    #+# #+#            #+#       #+#   #+#  #+#        #+#           #
####    #### #########  ##########    ##########  #######  ########## ##########      #  
