#!/usr/bin/perl -w
#Compactness of a protein conformation through
#radius of gyration and moment of inertia.
#Higher the value-higher the compactness of the conformation.

use strict;
use warnings;
use List::Util qw(sum);

unless($#ARGV==1){
	print "USAGE: perl RogMoi.pl file.pdb chainID\n";
	exit;
}

my $id = $ARGV[0];
open(IN, $id);
my $code = substr($id, 0,4);
my $chain = $ARGV[1];
my $s1 = 0;
my $s2 = 0;
my $s3 = 0;

my $cx=0;
my $cy=0;
my $cz=0;

my @pdb;
my @new;
my @new1;

while(<IN>){
	my @col = split;
	next unless $col[0] eq 'ATOM' and $col[2] eq 'CA' and $col[4] eq $chain;
	push @pdb, [@col[2,3,5,6,7,8]];
}

for (my $i=0;$i<=$#pdb;$i++){
	my($a, $r, $n, $x, $y, $z) = @{$pdb[$i]};
	
	$s1 = $s1+$x;
	$cx++;
	$s2 = $s2+$y;
	$cy++;
	$s3 = $s3+$z;
	$cz++;
}

#Coordinates of COM

my $X = sprintf "%0.3f", $s1/$cx;
my $Y = sprintf "%0.3f", $s2/$cy;
my $Z = sprintf "%0.3f", $s3/$cz;

#distance of each residue from COM

for my $j(0..$#pdb){
	my($a1, $r1, $n1, $x1, $y1, $z1) = @{$pdb[$j]};
	if($a1 eq 'CA'){
	my $dist = sprintf "%0.3f", sqrt(($X-$x1)**2 + ($Y-$y1)**2 + ($Z-$z1)**2);
	my $dis_rog = sprintf "%0.2f", (($dist)**2)*12.01; 	#for Radius of gyration
	my $dist_moi = sprintf "%0.2f", ($dist)**2;		#for Moment of inertia
	push(@new1,"$dist_moi\n");
	push(@new, "$dis_rog\n");
	}
}
my $length = @pdb;

#calculation - radius of gyration
my $sum_mass = $length*12.01;
my $sum_distances_rog = sum @new;
my $rog = sqrt($sum_distances_rog/$sum_mass);

#calculation - Moment of inertia
my $sum_distances_moi = sum @new1;
my $moi = ($sum_distances_moi/$length);
print "$code: Radius of gyration: $rog\n";
print "$code: Moment of inertia: $moi\n";

#End of program
