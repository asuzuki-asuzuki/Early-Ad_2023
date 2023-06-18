#!/usr/bin/perl

use strict ;
use warnings ;

#ARGV[0] --- ID

#input: output file of extract_HP_MNV_qname.pl
open (IN, "${ARGV[0]}_T_HP_MNV_qname.txt") or die ;

my $line ;
my @LINE ;
my @flag = ("REF_1", "REF_2", "REF_un", "REF_otherPS", "ALT_1", "ALT_2", "ALT_un", "ALT_otherPS", "OTHERS_1", "OTHERS_2", "OTHERS_un", "OTHERS_otherPS") ;
my @QNAME ;
my @qname ;
my %count ;
my $i ;

my $title ;

while (<IN>) {
    $line = $_ ;
    chomp ($line) ;
    if ($line =~ /^Chr\t/) {
	$title = $line ;
	next ;
    }
    @LINE = () ;
    @LINE = split (/\t/, $line) ;
    if ($line =~ /^[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t/) {
	@QNAME = () ;
	@QNAME = split (/\t/, "$'") ;
    } else {
	die ;
    }
    unless (@QNAME == @flag) {
	die ;
    }
    for ($i = 0 ; $i < @flag ; $i++) {
	if ($QNAME[$i] eq "0") {
	    next ;
	}
	@qname = () ;
	@qname = split (/,/, $QNAME[$i]) ;
	foreach my $item (@qname) {
	    $count{$LINE[0]}{$LINE[1]}{$LINE[2]}{$LINE[3]}{$flag[$i]} .= $item."," ;
	} 
    }
}

close (IN) ;

open (OUT, ">${ARGV[0]}_T_HP_MNV.txt") or die ;

print OUT "$title\n" ;

#input: output file of list_phased.pl
open (IN, "${ARGV[0]}_T.phased.list") or die ;

my $len ;
my %ok_qname ;
my %count_qname ;

my $total ; 
my %ok ;
my %OK ;
my $HP ;

my @refref ;
my @altalt ;
my $p ;

while (<IN>) {
    $line = $_ ;
    chomp ($line) ;
    @LINE = () ;
    @LINE = split (/\t/, $line) ;
    if (length ($LINE[2]) == 1 and length ($LINE[3]) == 1) { # SNV
        next ;
    } elsif (length ($LINE[2]) == length ($LINE[3])) { # MNV
        unless ($LINE[4] eq "MNV") {
            die ;
        }
	$len = length ($LINE[2]) ;
	%count_qname = () ;
	%ok_qname = () ;
	@refref = () ;
	@refref = split (//, $LINE[2]) ;
	@altalt = () ;
        @altalt = split (//, $LINE[3]) ;
        for ($i = 0 ; $i < $len ; $i++) {
	    $p = $LINE[1] + $i ;
	    foreach my $item (@flag) {
		if (exists ($count{$LINE[0]}{$p}{$refref[$i]}{$altalt[$i]}{$item})) {
		    while ($count{$LINE[0]}{$p}{$refref[$i]}{$altalt[$i]}{$item} =~ /[^,]+/g) {
			$ok_qname{$&}++ ;
			$count_qname{$item}{$&}++ ;
		    }
		}
	    }
	}
	$total = 0 ;
	%ok = () ;
	foreach my $item (keys %ok_qname) {
	    if ($ok_qname{$item} == $len) { # count qname == length of MNV (cover!)
		$total++ ;
		%OK = () ;
		foreach my $item2 (@flag) {
		    if ($item2 =~ /_([^_]+)$/) {
			$HP = $1 ;
		    } else {
			die ;
		    }
		    if (exists ($count_qname{$item2}{$item})) {
			if ($count_qname{$item2}{$item} == $len) { # all qname in the same hp and category
			    $OK{$item2} = "kekeke" ;
			} else {
			    $OK{"OTHERS_".$HP} = "kekeke" ; 
			}
		    }
		}
		unless (keys %OK == 1) {
		    die ;
		}
		foreach my $item2 (keys %OK) {
		    $ok{$item2}++ ;
		}
	    } else {
		next ;
	    }
	}
	print OUT "$LINE[0]\t$LINE[1]\t$LINE[2]\t$LINE[3]\t$total\t$LINE[5]" ;
	foreach my $item (@flag) {
	    $ok{$item} += 0 ;
	    print OUT "\t$ok{$item}" ;
	}
	print OUT "\n" ;
    } else {
	next ; #indel
    }
}

close (OUT) ;
