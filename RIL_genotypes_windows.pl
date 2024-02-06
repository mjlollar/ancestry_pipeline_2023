#!/usr/bin/perl -w
#use strict;

#my $OutputFile = '48old_418old_50kb_genotypes_121923.txt';
my $OutputFile = 'test.txt';

my $WinFilePrefix = 'windows_50kb_Chr';

my @chrs = ('X','2L','2R','3L','3R');

my @RILS = ("F1");
#my @RILs = ("PO-4", "PH-15");
#my @RILs = ("F1", "Q-14", "CF-6");

my @handles = @RILs;
for ($i = 0; $i < @handles; $i++){
  $handles[$i] =~ s/1/A/;
  $handles[$i] =~ s/2/B/;
  $handles[$i] =~ s/3/C/;
  $handles[$i] =~ s/4/D/;
  $handles[$i] =~ s/5/E/;
  $handles[$i] =~ s/6/F/;
  $handles[$i] =~ s/7/G/;
  $handles[$i] =~ s/8/H/;
  $handles[$i] =~ s/9/I/;
  $handles[$i] =~ s/0/J/;
}

my $c = 0;
my $f = 0;
my $i = 0;
my $j = 0;
my $r = 0;
my $s = 0;
my $w = 0;
my $chr = '';
my $file = '';
my $site = 0;
my $Win2Sums = 0;
my $Win1Sums = 0;
my $Win0Sums = 0;

my @line = ();
my @WinStarts = ();
my @WinStops = ();
my @InputAoA = ();
my @ChrGenoAoA = ();

#open output file
open O, ">$OutputFile";
print O "chr\tWinStart\tWinStop";
for ($r = 0; $r < @RILs; $r++){
  print O "\t$RILs[$r]";
}
print O "\n";

#for each chr arm
for ($c = 0; $c < @chrs; $c++){
  $chr = $chrs[$c];
  @WinStarts = ();
  @WinStops = ();
  @ChrGenoAoA = ();
  @line = ();
  for ($w = 0; $w < @WinStops; $w++){
    push @ChrGenoAoA, [ @line ];
  }
  
#open window file
  $file = $WinFilePrefix . $chr . '.txt';
  open W, "<$file" or die;
  while (<W>){
    chomp;
    last if m/^$/;
    @line = split;
    push @WinStarts, $line[0];
    push @WinStops, $line[1];
  }
  close W;
  $w = 0;
  
#open and analyze each RIL file
  for ($r = 0; $r < @RILs; $r++){
    @InputAoA = ();
    $file = $RILs[$r] . '.posterior';
    open $handles[$r], "<./$chr/$file" or die "can not open ./$chr/$file\n";
    $file = $handles[$r];
    scalar (<$file>);
    while (<$file>){
      chomp;
      last if m/^$/;
      @line = split;
      push @InputAoA, [ @line ];
    }
    close $file;
    $site = 0;
    for ($w = 0; $w < @WinStops; $w++){
      $Win2Sums = 0;
      $Win1Sums = 0;
      $Win0Sums = 0;
      for ($s = $site; $s < @InputAoA; $s++){
	$site = $s;
	last if ($InputAoA[$s][1] > $WinStops[$w]);
	$Win2Sums += $InputAoA[$s][2];
	$Win1Sums += $InputAoA[$s][3];
	$Win0Sums += $InputAoA[$s][4];
      }
      if (($Win2Sums + $Win1Sums + $Win0Sums) == 0){
	push @{$ChrGenoAoA[$w]}, -999;
	next;
      }
      if ($Win2Sums > $Win1Sums){
	if ($Win2Sums > $Win0Sums){
	  push @{$ChrGenoAoA[$w]}, 2;
	}
	else{
	  push @{$ChrGenoAoA[$w]}, 0;
	}
      }
      else{
	if ($Win1Sums > $Win0Sums){
	  push @{$ChrGenoAoA[$w]}, 1;
	}
	else{
	  push @{$ChrGenoAoA[$w]}, 0;
	}
      }
###
#      print "$RILs[$r] $site $ChrGenoAoA[$w][-1] $InputAoA[$s][2] $InputAoA[$s][3] $InputAoA[$s][4]\n";
#      last;
###      
    }
    print "Finished chr $chr, RIL $RILs[$r]\n";
  }

###
#  die;
  $j = @ChrGenoAoA;
  print "$j x ";
  $j =  @{$ChrGenoAoA[0]};
  print "$j\n";
#  for ($w = 0; $w < @ChrGenoAoA; $w++){
#    for ($r = 0; $r < @{$ChrGenoAoA[0]}; $r++
###
    
#add chr genotypes to output
  for ($w = 0; $w < @WinStops; $w++){
    print O "$chr\t$WinStarts[$w]\t$WinStops[$w]";
    for ($r = 0; $r < @{$ChrGenoAoA[0]}; $r++){
      print O "\t$ChrGenoAoA[$w][$r]"
    }
    print O "\n";
  }
}
close O;



    
