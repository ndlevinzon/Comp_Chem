#!/usr/bin/perl
#LibGen.pl: Generate AMBER FF Parameters from Ligands in Seperate .PDB Files
use warnings;
use 5.010;
use Module::Load;
use File::Basename;

# Load Modules
require "$ENV{\"LMOD_PKG\"}/init/perl";
module("load openbabel/2.4.1");                                 # OpenBabel
module("load gaussian16/SSE4.C01");                             # GAUSSIAN
module("load gcc/8.5.0 intel-oneapi-mpi/2021.4.0 amber/20.20"); # AMBER


# Set Working Directory Here
$home = '';
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
			@extensions = (".pdb", ".com", "_opt.com", "_HF631Gs_charge_nosym.com", "_opt.chk", "_HF631Gs_charge_nosym.chk", "_HF631Gs_charge_nosym.log", "_HF631Gs_charge_nosym.mol2", ".frcmod", ".tleap", ".lib");
			my $pdb             = $file.$extensions[0];
			my $com             = $file.$extensions[1];
			my $opt             = $file.$extensions[2];
			my $charge          = $file.$extensions[3];
			my $opt_check       = $file.$extensions[4];
			my $charge_check    = $file.$extensions[5];
			my $log             = $file.$extensions[6];
			my $mol2            = $file.$extensions[7];
			my $frcmod          = $file.$extensions[8];
			my $tleap           = $file.$extensions[9];
			my $lib             = $file.$extensions[10];

			# Convert .PDB to .COM
			if (! -e $com)
			{
				system("obabel -ipdb $pdb -ocom $com -m");
			}

			# Create GeomOpt Gaussian Input
			if (! -e $opt_check)
			{
				open my $in,  '<', $com  or die "Can't read COM file: $!";
				open my $out, '>', $opt  or die "Can't write OPT file: $!";
				print $out "%chk=$opt_check\n#P B3LYP/6-31G* Opt Freq=NoRaman\n\n$opt\n\n";
				while( <$in> )
    				{
					next if $. <= 4;
   					print $out $_;
   				}
				close $out;
			}

			# Create Charge Gaussian Input
			if (! -e $charge_check)
			{
				open my $in,  '<', $com    or die "Can't read COM file: $!";
				open my $out, '>', $charge or die "Can't write CHARGE file: $!";
				print $out "%chk=$charge_check\n#P HF/6-31G* Pop=MK iop(6/33=2,6/42=6) pop=mk scf=tight Test\n\n$charge\n\n";	
				while( <$in> )
    				{
					next if $. <= 4;
   					print $out $_;
   				}
				close $out;
			}
			
			# Run Gaussian
			if (! -e $log)
			{
				system("g16 $opt");
				system("g16 $charge");
			}

			# Get Charge from Gaussian .COM
			open my $in,  '<', $com  or die "Can't read COM file: $!";
			while( <$in> )
    			{
				next unless 5 .. undef;
				my @spin_charge = split('', $_, length($_))			
   				}
			close $in;
			
			# Get Charge from Gaussian .COM
			open my $in,  '<', $com  or die "Can't read COM file: $!";
			my $start_line = 5;
			while( <$in> )
    			{
				if( $. == $start_line ) 
				{ 
       					$line = $_;
					my @charge_spin = split('', $line, length($line));			
					$charge = $charge_spin[0].$charge_spin[1];
        				last;
				}
   			}
			close $in;

			# Run AMBER Antechamber
			if (! -e $lib)
			{
				system("antechamber -i $log -fi gout -o $mol2 -fo mol2 -c resp -nc $charge");
				system("parmchk2 -i $mol2 -f mol2 -o $frcmod");
				open my $tLEAP, '>', $tleap or die "Can't write tLEAP input file: $!";
				print $tLEAP "source leaprc.gaff\nloadamberparams $frcmod\nlig = loadmol2 $mol2\ncheck lig\nsaveoff lig $lib\nquit";
				close $tLEAP;
				system("tleap -f $tleap");
			}
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
