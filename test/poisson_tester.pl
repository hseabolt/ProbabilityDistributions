#!/usr/bin/perl

# binomial_tester.pl
# Some driver code to test the Bionomial.pm package

use strict;
use warnings;
use Poisson;

my @lambdas = qw(1 4 10);
my @ks = qw(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20);

foreach my $lambda ( @lambdas )	{
	print "lambda = $lambda\n";
	my $Pois = Poisson->new( "lambda" => $lambda );
	foreach my $k ( @ks )	{
		my $prob = $Pois->probability($k);
		print "k = $k    --- PoissonPr(X = k) = $prob\n";
	}
	print "\n\n";
}

# Estimate K from a PoissonPr and a given lambda
print "\n\nEstimate K --> PoissonPr = 0.156, lambda=4\n";		# K should equal 1
my $Pois = Poisson->new( "lambda" => 4 );
my $estk = $Pois->estimate_k(0.156);
print " === K = $estk === \n\n";

# Maximum Likelihood Estimation of p from a list of k-sampled trials
print "\n\nMLE estimation of lambda from a list of samples.\n";	
my @samples = qw(1 1 1 1 1 2);
my $lambda = Poisson->lambda_mle( \@samples );
print " === lambda =~ $lambda === \n\n";

# Draw random samples from a Poisson distribution
$Pois->{_lambda} = 4;
my $random_sample = $Pois->rpoisson(75);
print "\n\nRandomly sample a Poisson distribution with parameter lambda=4\n";
print join(" ", @{$random_sample}), "\n\n";

###########################################################################
$Pois->{_lambda} = 50;
my $cdf_prob = $Pois->probability( 400 );
print "CDF Pois(l=100, k=20) = $cdf_prob\n";


exit;