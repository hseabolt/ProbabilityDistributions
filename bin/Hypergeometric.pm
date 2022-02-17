#!/usr/bin/perl

# package Hypergeometric.pm
# An HS custom Perl package for Hypergeometric distributions, properties, and sampling methods.

# Author: MH Seabolt
# Last Updated: 5-18-2021	

package Hypergeometric; 

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
use POSIX;				# Using this for ceil() and floor() functions
our $AUTOLOAD;

# @INC libraries for my PC and Biolinux or HS custom classes 
use lib '/media/hunter/Data/scripts';
use lib 'D:\scripts';
use lib '/scicomp/home/ngr8/Biolinux/scripts';

########################################################################################################################
# Describes the probability of k successes (random draws for which the object drawn has a specified feature) in 
# n draws, without replacement, from a finite population of size N that contains exactly K objects with that feature, 
# wherein each draw is either a success or a failure. 
#
# In contrast, the binomial distribution describes the probability of k successes in n draws with replacement.
########################################################################################################################

####################################################################################################################
# Class data and methods 
# Attributes 
{
	my %_attribute_properties = (
		_N => '',  						# Total (finite) population size (e.g. a bag of colored marbles)
		_K => '',						# Number of objects of a specific type (e.g. green marbles)
		_n => '',						# number of "draws" / Bernoulli trials
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
  
	# Just check that the values of N, K, and n are are positive integer numbers (allowing any of them to be zero)
	return if ( $self->{_N} <= 0 && ($self->{_N} %% 1 != 0) ); 
	return if ( $self->{_K} <= 0 && ($self->{_K} %% 1 != 0) ); 
	return if ( $self->{_n} <= 0 && ($self->{_n} %% 1 != 0) ); 
	
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

# Getters and setters for N, K, and n since we are not using AUTOLOAD here
sub get_N	{	return $_[0]->{_N};	}
sub set_N	{	$_[0]->{_N} = $_[1];}
sub p 		{	return $_[0]->{_N};	}

sub get_K	{	return $_[0]->{_K};	}
sub set_K	{	$_[0]->{_K} = $_[1];}
sub K 		{	return $_[0]->{_K};	}

sub get_n	{	return $_[0]->{_n};	}
sub set_n	{	$_[0]->{_n} = $_[1];}
sub n 		{	return $_[0]->{_n};	}

####################################################################################################################

# Note: many of the below routines use $_[0] in place of defining $self in @_ --> this is useful for algorithmic speed.

################   SUBROUTINES     ###################
# Probability mass function and various aliases.  
# Requires additional input value k for the distribution.
# Hypergeometric probability --> P(X=k) = choose(K,k)*choose(N-K,n-k) / choose(N,n)
# Calculates the probability that there are k failures before the first success,
sub dhypergeometric	{	
	my ( $self, $k ) = @_;
	return (choose($self->{_K},$k) * choose($self->{_N}-$self->{_K},$self->{_n}-$k))/choose($self->{_N},$self->{_n});
}
sub hypergeometric		{	
	my ( $self, $k ) = @_;
	return (choose($self->{_K},$k) * choose($self->{_N}-$self->{_K},$self->{_n}-$k))/choose($self->{_N},$self->{_n});
}
sub probability		{	
	my ( $self, $k ) = @_;
	return (choose($self->{_K},$k) * choose($self->{_N}-$self->{_K},$self->{_n}-$k))/choose($self->{_N},$self->{_n});
}
sub pmf		{	
	my ( $self, $k ) = @_;
	return (choose($self->{_K},$k) * choose($self->{_N}-$self->{_K},$self->{_n}-$k))/choose($self->{_N},$self->{_n});
}
sub hypergeom	{	
	my ( $self, $k ) = @_;
	return (choose($self->{_K},$k) * choose($self->{_N}-$self->{_K},$self->{_n}-$k))/choose($self->{_N},$self->{_n});
}

# Mean and variance, plus aliases
sub expected	{	return $_[0]->{_n}*$_[0]->{_K}/$_[0]->{_N};		}
sub mean		{	return $_[0]->{_n}*$_[0]->{_K}/$_[0]->{_N};		}
sub variance	{	return ($_[0]->{_n}*($_[0]->{_K}/$_[0]->{_N}))*(($_[0]->{_N}-$_[0]->{_K})/$_[0]->{_N})*(($_[0]->{_N}-$_[0]->{_n})/($_[0]->{_N}-1));		}
sub stdev		{	return (($_[0]->{_n}*($_[0]->{_K}/$_[0]->{_N}))*(($_[0]->{_N}-$_[0]->{_K})/$_[0]->{_N})*(($_[0]->{_N}-$_[0]->{_n})/($_[0]->{_N}-1)))**0.5;	}
sub sd			{	return (($_[0]->{_n}*($_[0]->{_K}/$_[0]->{_N}))*(($_[0]->{_N}-$_[0]->{_K})/$_[0]->{_N})*(($_[0]->{_N}-$_[0]->{_n})/($_[0]->{_N}-1)))**0.5;	}

# Mode --> Hypergeometric has 2 modes (upper and lower)
sub mode		{ 	
	my ( $self ) = @_;
	my $upper = ceil((($self->{_n}+1)*($self->{_K}+1)/($self->{_N}+2))-1);
	my $lower = floor((($self->{_n}+1)*($self->{_K}+1)/($self->{_N}+2)));
	return ( $lower, $upper);	
}

# Note: The Hypergeometric distribution does not have a median

# Skewness and kurtosis
sub skewness	{	
	my ( $self ) = @_;
	my $num = ($self->{_N}-2*$self->{_K})*(($self->{_N}-1)**0.5)*($self->{_N}-2*$self->{_n});
	my $denom = (($self->{_n}*$self->{_K}*($self->{_N}-$self->{_K})*($self->{_N}-$self->{_n}))**0.5)*($self->{_N}-2);
	return $num/$denom;	
}
sub kurtosis	{	
	return 6 + (($_[0]->{_p}**2) / (1 - $_[0]->{_p}));		
}

# Cumulative distribution function and aliases
sub cumulative_distribution		{
	my ( $self, $k ) = @_;
	return 1 - ((choose($self->{_n},$k+1) * choose($self->{_N}-$self->{_n},$self->{_K}-$k-1))/choose($self->{_N},$self->{_K}));
}
sub cdf		{
	my ( $self, $k ) = @_;
	return 1 - ((choose($self->{_n},$k+1) * choose($self->{_N}-$self->{_n},$self->{_K}-$k-1))/choose($self->{_N},$self->{_K}));
}

###############################################
#     SAMPLE A HYPERGEOMETRIC DISTRIBUTION    #
###############################################

# Draw a sample of H random integers from a Hypereometric distribution characterized by N, K, and n.
# The returned list represents the number of successes (k) drawn from the given distribution (parameters N, K, and n are already characterized).
# rhypergeom() is an alias for sample().
sub sample 		{
	my @sample;
	while ( scalar @sample <= $_[1] )	{
		push @sample, $_[0]->_rand_hypergeometric();
	}
	return \@sample;
}
sub rhypergeom		{
	my @sample;
	while ( scalar @sample <= $_[1] )	{
		push @sample, $_[0]->_rand_hypergeometric();
	}
	return \@sample;
}

# (Internal subroutine) Generate a random Hypergeometric-distributed integer.
# This integer represents the number k of successes drawn from the given distribution.
# Unlike other distributions that have a given probability p associated with 
# successes, here we will use the mean as a stand-in.
sub _rand_hypergeometric		{
	my $L = $_[0]->{_K}/$_[0]->{_N};
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

# choose(n,k) is the number of ways to choose k elements from a set of n elements,
# when the order of selection is irrelevant.
# Given by the formula for n choose k =
#   (n - k)!
# -------------
#  k!(n - k)!
sub choose	{
	my ( $n, $k ) = @_;
	my ( $result, $j ) = ( 1, 1 );
	
	return 0 if $k > $n || $k < 0;
	$k = ($n-$k) if ( $n-$k ) < $k;
	
	while ( $j <= $k )	{
		$result *= $n--;
		$result /= $j++;
	}
	return $result;
}
####################################################################################################################


# Geometric probability of success here: 1.0 :-)
1;