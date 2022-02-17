#!/usr/bin/perl

# hypergeom_tester.pl
# Some driver code to test the Hypergeometric.pm package

use strict;
use warnings;

# @INC libraries for my PC and Biolinux or HS custom classes 
use lib '/media/hunter/Data/scripts';
use lib 'D:\scripts';
use lib '/scicomp/home/ngr8/Biolinux/scripts';
use Hypergeometric;

my @ns = qw(100 200 300);
my @Ks = qw(50 60 70);
my @ks = qw(0 1 2 3 4 5 6 7 8 9 10);
my $N = 500;

foreach my $n ( @ns )	{
	foreach my $K ( @Ks )	{
		print "X ~ Hypergeom($N, $K, $n)\n";
		my $Hypergeom = Hypergeometric->new( "N" => $N, "K" => $K, "n" => $n );
		foreach my $k ( @ks )	{
			my $prob = $Hypergeom->probability($k);
			print "k = $k    --- HypergeometricPr(X = k) = $prob\n";
		}
		print "\n\n";
	}
}

# Draw random samples from a Hypergeometric distribution
my $Hypergeom = Hypergeometric->new( "N" => $N, "K" => 70, "n" => 300 );
my $random_sample = $Hypergeom->rhypergeom(75);
print "\n\nRandomly sample a Hypergeometric distribution with parameters N=500, K=70, n=300\n";
print join(" ", @{$random_sample}), "\n\n";

exit;