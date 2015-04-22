#!/usr/bin/perl -w
#comparision of compactness of N and C terminal
#useful as a parameter in determination of protein cotranslational folding
# log(N-terminal)/(C-termminal)
# if the above value is negative N-terminal is more compact 
# otherwise the C-terminal

use strict;
use warnings;
use List::Util qw(sum min);
use POSIX;

unless($#ARGV==1){
	print "USAGE: perl rmin.pl file.pdb chainID\n";
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
	my $dis_sq = sprintf "%0.2f", (($dist)**2)*12.01; 	#for Radius of gyration
	my $dist_moi = sprintf "%0.2f", ($dist)**2;		#for Moment of inertia
	push(@new1,"$dist_moi\n");
	push(@new, "$dis_sq\n");
	}
}

# Moment of inertia of N and C terminal
my $length = @pdb;
my $N = ceil ($length/2);
my $C = $length-$N;
my $n_sum_moi = sum @new1[0..$N];
my $c_sum_moi = sum @new1[-$C..-1];
my $n_moi = $n_sum_moi/$N;
my $c_moi = $c_sum_moi/$C;
my $pre_moi = $n_moi/$c_moi;
my $MOI = log($pre_moi)/log(10);

# Radius of gyration of N and C terminal
my $N_term_mass =sprintf "%0.2f", ($N*12.01);
my $C_term_mass =sprintf "%0.2f", ($C*12.01);
my $N_dist = sum @new[0..$N];
my $C_dist = sum @new[-$C..-1];
my $N_rog = sprintf "%0.2f", sqrt($N_dist/$N_term_mass);
my $C_rog = sprintf "%0.2f", sqrt($C_dist/$C_term_mass);
my $NbyC = ($N_rog/$C_rog);
my $ROG = log($NbyC)/log(10);
print "ID\tROG\t\t\tMOI\n";
print "$code\t$ROG\t$MOI\n";

#End of program
