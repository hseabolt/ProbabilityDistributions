#!/usr/bin/perl

# package Geometric.pm
# An HS custom Perl package for Geometric distributions, properties, and sampling methods.

# Author: MH Seabolt
# Last Updated: 9-27-2020	

package Geometric; 

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
# Consider a sequence of trials, where each trial has only two possible outcomes (designated failure and success). 
# The probability of success is assumed to be the same for each trial. In such a sequence of trials, the geometric 
# distribution is useful to model the number of failures before the first success. The distribution gives the 
# probability that there are zero failures before the first success, one failure before the first success,
# two failures before the first success, and so on.
#
# Assumptions: When is the Geometric distribution an appropriate model?
# The geometric distribution is an appropriate model if the following assumptions are true.
# -- The phenomenon being modeled is a sequence of independent trials.
# -- There are only two possible outcomes for each trial, often designated success or failure.
# -- The probability of success, p, is the same for every trial.
# -- If these conditions are true, then the geometric random variable Y is the count of the number of failures 
# 	 before the first success. The possible number of failures before the first success is 0, 1, 2, 3, and so on.
#
# An alternative formulation is that the geometric random variable X is the total number of trials up to and including
#  the first success, and the number of failures is X − 1. In the graphs above, this formulation is shown on the left.
####################################################################################################################

####################################################################################################################
# Class data and methods 
# Attributes 
{
	my %_attribute_properties = (
		_p => '',  						# Must be positive real number greater than 0 and less than or equal to 1.
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
  
	# Just check that the value of p is a positive real number greater than 0 and less than or equal to 1.
	return if ( $self->{_p} <= 0 && $self->{_p} > 1.0 ); 
	
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

# Get and set the value of p
sub get_p	{	return $_[0]->{_p};	}
sub set_p	{	$_[0]->{_p} = $_[1];}
sub p 		{	return $_[0]->{_p};	}

####################################################################################################################

# Note: many of the below routines use $_[0] in place of defining $self in @_ --> this is useful for algorithmic speed.

################   SUBROUTINES     ###################
# Probability mass function and various aliases.  
# Requires additional input value K for the distribution.
# Geometric probability --> P(X=k) = p*(1-p)^(k)
# Calculates the probability that there are k failures before the first success, where the argument "p" is the probability of success on each trial.
sub dgeometric	{	
	my ( $self, $k ) = @_;
	return $self->{_p} * (1 - $self->{_p})**($k);
}
sub geometric		{	
	my ( $self, $k ) = @_;
	return $self->{_p} * (1 - $self->{_p})**($k);
}
sub probability		{	
	my ( $self, $k ) = @_;
	return $self->{_p} * (1 - $self->{_p})**($k);
}
sub pmf		{	
	my ( $self, $k ) = @_;
	return $self->{_p} * (1 - $self->{_p})**($k);
}
sub geom	{	
	my ( $self, $k ) = @_;
	return $self->{_p} * (1 - $self->{_p})**($k);
}

# Mean and variance, plus aliases
sub expected	{	return (1 - $_[0]->{_p}) / $_[0]->{_p};		}
sub mean		{	return (1 - $_[0]->{_p}) / $_[0]->{_p};		}
sub variance	{	return (1 - $_[0]->{_p}) / $_[0]->{_p}**2;		}
sub stdev		{	return ((1 - $_[0]->{_p}) / $_[0]->{_p}**2)**0.5;	}
sub sd			{	return ((1 - $_[0]->{_p}) / $_[0]->{_p}**2)**0.5;	}

# Median, mode
sub median		{	return int( -1 / $_[0]->log2(1 - $_[0]->{_p})) - 1;		}
sub mode		{ 	return 0;		}

# Skewness and kurtosis
sub skewness	{	return (2 - $_[0]->{_p}) / (1 - $_[0]->{_p})**0.5;		}
sub kurtosis	{	return 6 + (($_[0]->{_p}**2) / (1 - $_[0]->{_p}));		}

# Cumulative distribution function and aliases
sub cumulative_distribution		{
	my ( $self, $k ) = @_;
	return 1 - (1 - $self->{_p})**$k;
}
sub cdf		{
	my ( $self, $k ) = @_;
	return 1 - (1 - $self->{_p})**($k+1);
}

########################################
#       ESTIMATE P FROM A SAMPLE       #
########################################

# Maximum likelihood estimation of p, and alias
# The parameter p can be estimated by equating the expected value with the sample mean. 
# This is the method of moments, which in this case happens to yield maximum likelihood estimates of p.
sub geom_mle	{
	my @sample = @{$_[1]};
	my $bigsum;
	$bigsum += $_ foreach ( @sample );
	return scalar @sample / $bigsum;
}
sub estimate_p_from_sample		{
	my @sample = @{$_[1]};
	my $bigsum;
	$bigsum += $_ foreach ( @sample );
	return scalar @sample / $bigsum;
}

##########################################
#     SAMPLE A GEOMETRIC DISTRIBUTION    #
##########################################

# Draw a sample of N random integers from a Geometric distribution characterized by p.
# The returned list represents the number of failures (k) before the first success in each trial.
# rgeom() is an alias for sample().
sub sample 		{
	my @sample;
	while ( scalar @sample <= $_[1] )	{
		push @sample, $_[0]->_rand_geometric();
	}
	return \@sample;
}
sub rgeom		{
	my @sample;
	while ( scalar @sample <= $_[1] )	{
		push @sample, $_[0]->_rand_geometric();
	}
	return \@sample;
}

# (Internal subroutine) Generate a random Geometric-distributed integer.
sub _rand_geometric		{
	my $L = $_[0]->{_p};
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

sub log2 {
    return log($_[1])/log(2);
}

####################################################################################################################


# Geometric probability of success here: 1.0 :-)
1;