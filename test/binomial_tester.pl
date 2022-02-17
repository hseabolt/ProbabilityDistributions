#!/usr/bin/perl

# binomial_tester.pl
# Some driver code to test the Bionomial.pm package

use strict;
use warnings;
use Binomial;

my @ps = qw(0.3 0.5 0.7);
my @ns = qw(10 20 40 50);

for ( my $i=0; $i < scalar @ps; $i++ )	{
	for ( my $j=0; $j < scalar @ns; $j++ )	{
		print "P = $ps[$i]  --- N = $ns[$j]\n";
		my $Binom = Binomial->new( "p" => $ps[$i], "n" => $ns[$j] );
		
		# Binomial probability for each number of successes K from 0 -> N
		foreach my $k ( 0 .. $ns[$j] )	{	
			my $prob = $Binom->probability($k);
			print "k = $k    --- Binom(X=k) = $prob\n";
		}
	}
}

# Draw random samples from a Binomial distribution
print "\n\nRandomly sample a Binomial distribution with parameter N=50 and p=0.6.  The numbers returned are the number of successes out of N trials.\n";
my $Binom = Binomial->new( "p" => 0.6, "n" => 50 );
my $random_sample = $Binom->rbinom(75);
print join(" ", @{$random_sample}), "\n\n";

# Maximum Likelihood Estimation of p from a list of k-sampled trials
print "\n\nMLE estimation of p from a sampled set of Bernolli trials.\n";	
my @samples = qw(0 0 0 1 1 0 1 1 0 0 1 0 1 0 1 1 1 1 0 0 0);
my $p = Binomial->binomial_mle( \@samples );
print " === P =~ $p === \n\n";

exit;