#!/usr/bin/perl

# negbinom_tester.pl
# Some driver code to test the NegativeBinomial.pm package

use strict;
use warnings;

# @INC libraries for my PC and Biolinux or HS custom classes 
use lib '/media/hunter/Data/scripts';
use lib 'D:\scripts';
use lib '/scicomp/home/ngr8/Biolinux/scripts';
use NegativeBinomial;

my @ks = qw(5 15 85);
my @ps = qw(0.2 0.5 0.91);
my @rs = qw(0 1 2 3 4 5 6 7 8 9 10);

foreach my $r ( @rs )	{
	foreach my $p ( @ps )	{
		print "X ~ NB($r, $p)\n";
		my $NB = NegativeBinomial->new( "r" => $r, "p" => $p );
		foreach my $k ( @ks )	{
			my $prob = $NB->probability($k);
			print "k = $k    --- NB(X = k) = $prob\n";
		}
		print "\n\n";
	}
}

# Draw random samples from a Hypergeometric distribution
my $NB2 = NegativeBinomial->new( "p" => 0.22, "r" => 5 );
my $random_sample = $NB2->rnegbinom(75);
print "\n\nRandomly sample a Negative Binomial distribution with parameters p=0.22, r=5\n";
print join(" ", @{$random_sample}), "\n\n";

exit;