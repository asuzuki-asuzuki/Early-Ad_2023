#!/usr/bin/perl

use strict ;
use warnings ;

#ARGV[0] --- ID

#input: output file of list_phased.pl
open (IN, "${ARGV[0]}.phased.list") or die ; 

my $line ;
my @LINE ;
my %ref ;
my %alt ;
my %phased ;
my @P ;

my @refref ;
my @altalt ;
my $p ;
my $i ;

while (<IN>) {
    $line = $_ ;
    chomp ($line) ;
    if ($line =~ /^\#/) {
	next ;
    }
    @LINE = () ;
    @LINE = split (/\t/, $line) ;
    if (length ($LINE[2]) == 1 and length ($LINE[3]) == 1) { # SNV
	next ;
    } elsif (length ($LINE[2]) == length ($LINE[3])) { # MNV
	unless ($LINE[4] eq "MNV") {
	    die ;
	}
	@refref = () ;
 	@refref = split (//, $LINE[2]) ;
	@altalt = () ;
        @altalt = split (//, $LINE[3]) ;
	for ($i = 0 ; $i < @refref ; $i++) {
	    $p = $LINE[1] + $i ;
	    $ref{$LINE[0]}{$p} = $refref[$i] ;
	    $alt{$LINE[0]}{$p} = $altalt[$i] ;
	    $phased{$LINE[0]}{$p} = $LINE[5] ;
	    push (@P, "$LINE[0]\t$p\t$refref[$i]\t$altalt[$i]") ;
	}
    }
}

close (IN) ;

my @HP = (1, 2, "*", "otherPS") ; # HP1, HP2, unphased, other phased block

#input: output file of samtools mpileup
open (IN, "${ARGV[0]}.MNV.pileup") or die ;

my @hp ;
my @ps ;
my @qname ;
my $k ;
my $pileup ;
my $base ;
my %REF ;
my %ALT ;
my %OTHERS ;

my %ok ;

while (<IN>) {
    $line = $_ ;
    chomp ($line) ;
    @LINE = () ;
    @LINE = split (/\t/, $line) ;
    #HP tag
    unless ($LINE[7] =~ /^[12,\*]+$/) {
	print "$line\n" ;
	die ;
    }
    @hp = () ;
    @hp = split (/,/, $LINE[7]) ;
    #PS tag
    unless ($LINE[8] =~ /^[\d,\*]+$/) {
	print "$line\n" ;
	die ;
    }
    @ps = () ;
    @ps = split (/,/, $LINE[8]) ;
    #read name
    @qname = () ;
    @qname = split (/,/, $LINE[6]) ;
    #pileup
    $k = -1 ;
    $pileup = $LINE[4] ;
    %REF = () ;
    %ALT = () ;
    %OTHERS = () ;
    while ($pileup =~ /^./g) {
	$base = $& ;
	$pileup = "$'" ;
	if ($base =~ /^[\.,]$/) { # skip (N)
	    print "$line\n" ;
	    die ;
	} elsif ($base =~ /^[ATCGNatcgn]$/) {
	    if ($base =~ /^$ref{$LINE[0]}{$LINE[1]}$/i) { # REF
		$k++ ;
		if ($ps[$k] eq $phased{$LINE[0]}{$LINE[1]} or $ps[$k] eq "*") { # PS check
		    $REF{$hp[$k]} .= $qname[$k]."," ;
		} else {
		    $REF{"otherPS"} .= $qname[$k]."," ;
		}
            } elsif ($base =~ /^$alt{$LINE[0]}{$LINE[1]}$/i) { # ALT
		$k++ ;
		if ($ps[$k] eq $phased{$LINE[0]}{$LINE[1]} or $ps[$k] eq "*") { # PS check
		    $ALT{$hp[$k]} .= $qname[$k]."," ;
		} else {
		    $ALT{"otherPS"} .= $qname[$k]."," ;
		}
	    } else {
		$k++ ;
		if ($ps[$k] eq $phased{$LINE[0]}{$LINE[1]} or $ps[$k] eq "*") { # PS check
		    $OTHERS{$hp[$k]} .= $qname[$k]."," ;
		} else {
		    $OTHERS{"otherPS"} .= $qname[$k]."," ;
		}
	    }
        } elsif ($base =~ /^[\*]$/) {
	    $k++ ;
	    if ($ps[$k] eq $phased{$LINE[0]}{$LINE[1]} or $ps[$k] eq "*") { # PS check
		$OTHERS{$hp[$k]} .= $qname[$k]."," ;
	    } else {
		$OTHERS{"otherPS"} .= $qname[$k]."," ;
	    }
	} elsif ($base eq "+") { # information insertion after a base
	    if ($pileup =~ /^(\d+)/) {
		$pileup = substr ("$'", $&) ;
	    } else {
		print "$line\n" ;
		die ;
	    }
	} elsif ($base eq "-") { # information insertion after a base
            if ($pileup =~ /^(\d+)/) {
                $pileup = substr ("$'", $&) ;
	    } else {
                print "$line\n" ;
                die ;
            }
	} elsif ($base =~ /^[\^]$/) { # read head
	    $pileup = substr ($pileup, 1) ;
	} elsif ($base =~ /^[\$]$/) { # read tail 
	    next ;
	} else {
	    print "$line\n" ;
	    die ;
	}
    }
    unless (@hp == $k + 1) {
	print "$line\n" ;
	die ;
    }
    $ok{$LINE[0]}{$LINE[1]}{$ref{$LINE[0]}{$LINE[1]}}{$alt{$LINE[0]}{$LINE[1]}} = $LINE[3]."\t".$phased{$LINE[0]}{$LINE[1]} ;
    foreach my $item (@HP) {
	if (exists ($REF{$item})) {
	    $ok{$LINE[0]}{$LINE[1]}{$ref{$LINE[0]}{$LINE[1]}}{$alt{$LINE[0]}{$LINE[1]}} .= "\t$REF{$item}" ;
	} else {
	    $ok{$LINE[0]}{$LINE[1]}{$ref{$LINE[0]}{$LINE[1]}}{$alt{$LINE[0]}{$LINE[1]}} .= "\t0" ;
	}
    }
    foreach my $item (@HP) {
	if (exists ($ALT{$item})) {
	    $ok{$LINE[0]}{$LINE[1]}{$ref{$LINE[0]}{$LINE[1]}}{$alt{$LINE[0]}{$LINE[1]}} .= "\t$ALT{$item}" ;
	} else {
	    $ok{$LINE[0]}{$LINE[1]}{$ref{$LINE[0]}{$LINE[1]}}{$alt{$LINE[0]}{$LINE[1]}} .= "\t0" ;
	}
    }
    foreach my $item (@HP) {
	if (exists ($OTHERS{$item})) {
	    $ok{$LINE[0]}{$LINE[1]}{$ref{$LINE[0]}{$LINE[1]}}{$alt{$LINE[0]}{$LINE[1]}} .= "\t$OTHERS{$item}" ;
	} else {
	    $ok{$LINE[0]}{$LINE[1]}{$ref{$LINE[0]}{$LINE[1]}}{$alt{$LINE[0]}{$LINE[1]}} .= "\t0" ;
	}
    }
}

close (IN) ;


open (OUT, ">${ARGV[0]}_HP_MNV_qname.txt") or die ;

print OUT "Chr\tPosition\tRef\tAlt\tTotal\tPS\t" ;
print OUT "REF_1\tREF_2\tREF_un\tREF_otherPS\t" ;
print OUT "ALT_1\tALT_2\tALT_un\tALT_otherPS\t" ;
print OUT "OTHERS_1\tOTHERS_2\tOTHERS_un\tOTHERS_otherPS\n" ;

foreach my $item (@P) {
    @LINE = split (/\t/, $item) ;
    if (exists ($ok{$LINE[0]}{$LINE[1]}{$LINE[2]}{$LINE[3]})) {
	print OUT "$LINE[0]\t$LINE[1]\t$LINE[2]\t$LINE[3]\t" ;
	print OUT "$ok{$LINE[0]}{$LINE[1]}{$LINE[2]}{$LINE[3]}\n" ;
    } else {
	print OUT "$LINE[0]\t$LINE[1]\t$LINE[2]\t$LINE[3]\t0\t$phased{$LINE[0]}{$LINE[1]}\t" ;
	print OUT "0\t0\t0\t0\t" ;
        print OUT "0\t0\t0\t0\t" ;
	print OUT "0\t0\t0\t0\n" ;
    }
}

close (OUT) ;
