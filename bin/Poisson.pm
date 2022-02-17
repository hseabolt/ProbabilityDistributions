#!/usr/bin/perl

# package Poisson.pm
# An HS custom Perl package for Poisson distributions, properties, and sampling methods.

# Author: MH Seabolt
# Last Updated: 9-27-2020	

package Poisson; 

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(max min sum lg log10 pow round);			 # Standard import from other packages

use strict;
use warnings;
use List::Util qw( max min sum first shuffle uniq );
use Scalar::Util;
use Storable qw(dclone);
use Data::Dumper;
use Carp;
our $AUTOLOAD;

# @INC libraries for my PC and Biolinux or HS custom classes 
use lib '/media/hunter/Data/scripts';
use lib 'D:\scripts';
use lib '/scicomp/home/ngr8/Biolinux/scripts';

####################################################################################################################
# The Poisson distribution is popular for modeling the number of times a discrete event occurs in an interval of time or space.
# The Poisson distribution may be useful to model events such as:
# -- The number of meteorites greater than 1 meter diameter that strike Earth in a year
# -- The number of patients arriving in an emergency room between 10 and 11 pm
# -- The number of laser photons hitting a detector in a particular time interval
#
# Assumptions and validity
# The Poisson distribution is an appropriate model if the following assumptions are true:
# -- k is the number of times an event occurs in an interval and k can take values 0, 1, 2, ....
# -- The occurrence of one event does not affect the probability that a second event will occur. That is, events occur independently.
# -- The average rate at which events occur is independent of any occurrences. 
# 	 For simplicity, this is usually assumed to be constant, but may in practice vary with time.
# -- Two events cannot occur at exactly the same instant; instead, at each very small sub-interval
#	 exactly one event either occurs or does not occur.
#
# If these conditions are true, then k is a Poisson random variable, and the distribution of k is a Poisson distribution.
# The Poisson distribution is also the limit of a binomial distribution, for which the probability 
# of success for each trial equals λ divided by the number of trials, as the number of trials approaches infinity.
# In cases where the set of Bernoulli trials is suitably large and the number of successes is small, 
# the Poisson distribution is a suitable approximation of the Binomial distribution.
#
# Examples which violate the Poisson assumptions:
# -- The number of students who arrive at the student union per minute will likely not follow a Poisson distribution,
# 	 because the rate is not constant (low rate during class time, high rate between class times) and the arrivals
# 	 of individual students are not independent (students tend to come in groups).
# -- The number of magnitude 5 earthquakes per year in a country may not follow a Poisson distribution
# 	 if one large earthquake increases the probability of aftershocks of similar magnitude.
####################################################################################################################

####################################################################################################################
# Class data and methods 
# Attributes 
{
	my %_attribute_properties = (
		_lambda => '',  						# Must be positive real number greater than or equal to 0
	);
	
	# Global variable counter
	my $_count = 0;
	
	# Return a list of all attributes
	sub _all_attributes	{
		keys %_attribute_properties;
	}
	
	# Return the default value for a given attribute
    	sub _attribute_default 	{
     	my( $self, $attribute ) = @_;
        	$_attribute_properties{$attribute};
    	}
    
	# Manage the count of existing objects
	sub get_count	{
		$_count;
	}
	sub _incr_count	{
		++$_count;
	}
	sub _decr_count	{
		--$_count;
	}	
}

############################################################
#                       CONSTRUCTORS                       #
############################################################

# The contructor method
# Construct a new graph (my $node = Markov->new() );
# Returns a scalar reference to a
sub new				{
	my ( $class, %arg ) = @_;
	
	# Create the new object
	my $self = bless {}, $class;

	foreach my $attribute ( $self->_all_attributes() ) {
		# E.g. attribute = "_name",  argument = "name"
		my ($argument) = ( $attribute =~ /^_(.*)/ );
		# If explicitly given
		if (exists $arg{$argument}) 	{
			$self->{$attribute} = $arg{$argument};
		}
		else	{
			$self->{$attribute} = $self->_attribute_default($attribute);
		}
   	}
	
	# Construction code specific to this class
	######################################################################################
  
	# Just check that the value of lambda is a positive real number  
	return if ( $self->{_lambda} < 0 ); 
	
	######################################################################################
	
    $class->_incr_count();
	return $self;
}

# The clone method
# All attributes are copied from the calling object, unless specifically overriden
# Called from an existing object ( Syntax: $cloned_obj = $obj->clone(); )
sub clone	{
	my ( $caller, %arg ) = @_;
	# Extract the class name from the calling object
	my $class = ref($caller);
		
	# Create a new object
	my $self = bless {}, $class;
		
	foreach my $attribute ( $self->_all_attributes() )	{
		my ($argument) = ( $attribute =~ /^_(.*)/ );
			
		# If explicitly given
		if ( exists $arg{$argument} )	{
			$self->{$attribute} = $arg{$argument};
		}
			
		# Otherwise, copy attribute of new object from the calling object
		else	{
			$self->{$attribute} = $caller->{$attribute};
		}
	}
	$self->_incr_count();
	return $self;
}

# When an object is no longer being used, garbage collect it and adjust count of existing objects
sub DESTROY	{
	my ( $self ) = @_;
	$self->_decr_count();
}

# Get and set the value of lambda
sub get_lambda	{	return $_[0]->{_lambda};	}
sub set_lambda	{	$_[0]->{_lambda} = $_[1];	}
sub lambda 		{	return $_[0]->{_lambda};	}

####################################################################################################################

# Note: many of the below routines use $_[0] in place of defining $self in @_ --> this is useful for speed.

################   SUBROUTINES     ###################
# Probability mass function and various aliases.  
# Requires additional input value K for the distribution.
# Poisson probability --> P(X=k)

sub poisson {
 my ($x, $a) = @_;
 return unless $a >= 0 && $x >= 0 && $x == int($x);
 return ($a ** $x) * exp(-$a) / factorial($x);
}
sub dpoisson	{	
	my ( $self, $k ) = @_;
	my $l = $self->{_lambda};
	return ($l**$k) * exp(-$l) / $self->_factorial($k);
}
sub poisson		{	
	my ( $self, $k ) = @_;
	my $l = $self->{_lambda};
	return ($l**$k) * exp(-$l) / $self->_factorial($k);
}
sub probability		{	
	my ( $self, $k ) = @_;
	my $l = $self->{_lambda};
	return ($l**$k) * exp(-$l) / $self->_factorial($k);
}
sub pmf		{	
	my ( $self, $k ) = @_;
	my $l = $self->{_lambda};
	return ($l**$k) * exp(-$l) / $self->_factorial($k);
}
sub pois		{	
	my ( $self, $k ) = @_;
	my $l = $self->{_lambda};
	return ($l**$k) * exp(-$l) / $self->_factorial($k);
}

# Mean and variance, plus aliases
sub expected	{	return $_[0]->{_lambda};		}
sub mean		{	return $_[0]->{_lambda};		}
sub variance	{	return $_[0]->{_lambda};		}
sub stdev		{	return $_[0]->{_lambda}**0.5;	}
sub sd			{	return $_[0]->{_lambda}**0.5;	}

# Median, mode
sub median		{	return int( ($_[0]->{_lambda}+(1/3)) - (0.02/$_[0]->{_lambda}) );		}
sub mode	{
	if ( $_[0]->{_lambda} == int($_[0]->{_lambda}) )		{	return ( $_[0]->{_lambda}, $_[0]->{_lambda} - 1 );	}
	else	{	return int($_[0]->{_lambda});	}
}

# Skewness and kurtosis
sub skewness	{	return $_[0]->{_lambda}**(-1/2);	}
sub kurtosis	{	return $_[0]->{_lambda}**-1;		}
sub coefficient_of_variation	{	return $_[0]->{_lambda}**(-1/2);	}

# Fisher information
sub fisher_information	{	return 1/$_[0]->{_lambda};		}

# Cumulative distribution function and aliases
sub cumulative_distribution		{
	my ( $self, $k ) = @_;
	my $bigsum;
	for ( my $i=0; $i <= int($k); $i++ )	{
		$bigsum += ( $self->{_lambda}**$i / $self->_factorial($i) );
	}
	return ( (exp($self->{_lambda} * -1)) * $bigsum);
}
sub cdf		{
	my ( $self, $k ) = @_;
	my $bigsum;
	return 0 if ( $k == 0 );
	for ( my $i=1; $i <= int($k); $i++ )	{
		$bigsum += ( $self->{_lambda}**$i / $self->_factorial($i) )
	}
	return ( (exp($self->{_lambda} * -1)) * $bigsum);
}

########################################
#    ESTIMATE LAMBDA FROM A SAMPLE     #
########################################

# Maximum likelihood estimation of lambda, and alias
# Given a sample of n measured values k_sub{i} in {0,1,...}, for i = 1, ..., n, 
# we wish to estimate the value of the parameter λ of the Poisson population from which the sample was drawn.
sub lambda_mle	{
	my @sample = @{$_[1]};
	my $bigsum;
	$bigsum += $_ foreach ( @sample );
	return (1/scalar @sample) * $bigsum;
}
sub estimate_lambda_from_sample		{
	my @sample = @{$_[1]};
	my $bigsum;
	$bigsum += $_ foreach ( @sample );
	return (1/scalar @sample) * $bigsum;
}

########################################
# ESTIMATE K FROM POISSON PROBABILITY  #
########################################

# Given a Poisson probability and the distribution it came from,
# Estimate the value of K.  This solution is crude and follows the CDF loop, 
# breaking when we exceed the probability value given.
sub estimate_k	{
	my ( $self, $p ) = @_;
	return 0 if ( $p == 0 );
	my $cdf = 0;
	my $k = 0;
	while ( $cdf < $p )	{
		my $bigsum = 0;
		for ( my $i=0; $i < $k; $i++ )	{
			$bigsum += ( $self->{_lambda}**$i / $self->_factorial($i) );
		}
		$cdf = exp(-$self->{_lambda}) * $bigsum;
		$k++;
	}
	return $k - 1;
}

########################################
#     SAMPLE A POISSON DISTRIBUTION    #
########################################

# Draw a sample of N random real numbers from a Poisson distribution characterized by lambda
# Follows the algorithm described by Knuth.  rpoisson is an alias for sample().
sub sample 		{
	my @sample;
	while ( scalar @sample <= $_[1] )	{
		push @sample, $_[0]->_rand_poisson();
	}
	return \@sample;
}
sub rpoisson		{
	my @sample;
	while ( scalar @sample <= $_[1] )	{
		push @sample, $_[0]->_rand_poisson();
	}
	return \@sample;
}

# (Internal subroutine) Generate a random Poisson-distributed number, following Knuth's algorithm.
sub _rand_poisson		{
	my $L = exp(-$_[0]->{_lambda});
	my $p = 1.0;
	my $k = 0;
	while ( $p > $L )	{
		$k++;
		$p *= rand();
	}
	return $k - 1;
}

########################################
#       BASIC UTILITY SUBROUTINES      #
########################################

# Computes the factorial of an integer value.
# Use the Gamma function for non-integers!
sub _factorial	{
	my ( $self, $n ) = @_;
	return undef unless $n >= 0 and $n == int($n);			# Note: Non-integers require the Gamma function!
	my $res = 1;
	$res *= $n-- while ( $n > 1 );
	return $res;
}

# Normalizes a list of numerical values	-- here probabilities	
sub _normalize	{
	my ( @list ) = @_;
	my $probabilities;
	my $total = sum @list;	
	foreach my $probability ( @list )	{
		push @{$probabilities}, ($probability/$total);
	}
	return $probabilities;
}
####################################################################################################################


# Poisson probability of success here: 1.0 :-)
1;