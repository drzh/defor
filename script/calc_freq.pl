#!/usr/bin/perl

# Usage: prog -i FILE_mpileup -d INT_min_depth -f FLOAT_min_freq -F FLOAT_max_freq -s INT_column_start -k BINARY_flag_filter -t BINARY_flag_head

# Usage: -k represent whether : 0) remove the whole line, 1) keep to the filtered out items as '.' or 2) keep it as is

use strict;
use warnings;
use Getopt::Long;

my $file = '/dev/stdin';
my $para_dmin = 30;
my $para_fmin = 0.01;
my $para_fmax = 0.99;
my $para_start = 3;
my $para_k = 0;
my $para_t = 0;

Getopt::Long::Configure("bundling");

GetOptions(
           'i|input:s' => \$file,
           'd:i' => \$para_dmin,
           'f:f' => \$para_fmin,
           'F:f' => \$para_fmax,
           's:i' => \$para_start,
	   'k' => \$para_k,
 	   't' => \$para_t,
          );

while (<>) {
  chomp;
  my @e = split /\t/, $_, -1;
  my $i = 4;
  my %stat;
  my %sum;
  my %sumsamp;
  if ($e[$i] ne '') {
    my $len = length($e[$i]);
    my $j = 0;
    while ($j < $len) {
      my $b = uc(substr($e[$i], $j, 1));
      if ($b eq '.' || $b eq ',') {
        $stat{$i}{$e[2]}++;
        $sumsamp{$i}++;
        $j++;
      }
      elsif ($b eq 'A' ||
             $b eq 'C' ||
             $b eq 'G' ||
             $b eq 'T' ||
             $b eq 'N') {
        $stat{$i}{$b}++;
        $sum{$b}++;
        $sumsamp{$i}++;
        $j++;
      }
      elsif ($b eq '+' || $b eq '-') {
        my $k = $j;
        $j += 2;
        while (substr($e[$i], $j, 1) ge '0' && substr($e[$i], $j, 1) le '9') {
          $j++;
        }
        my $n = substr($e[$i], $k + 1, $j - 1 - $k);
        $b = uc(substr($e[$i], $k, $j - $k + $n));
        $stat{$i}{$b}++;
        $sum{$b}++;
        $stat{$i}{$e[2]}--;
        $j += $n;
      }
      elsif ($b eq '^') {
        $j += 2;
      }
      else {
        $j++;
      }
    }
  }
  my @alt = sort {$sum{$b} <=> $sum{$a}} keys %sum;
  unshift @alt, $e[2];
  if (exists $sumsamp{$i}) {
    my $flag = 0;
    my $ale = '';
    my $cnt = '';
    my $per = '';
    foreach my $a (@alt) {
      if ($ale ne '') {
        $ale .= ',';
        $cnt .= ',';
        $per .= ',';
      }
      $ale .= $a;
      if (exists $stat{$i}{$a}) {
        $cnt .= $stat{$i}{$a};
        $per .= sprintf("%.3g", $stat{$i}{$a} / $sumsamp{$i});
      }
      else {
        $cnt .= '0';
        $per .= '0';
      }
    }
    my $dep = 0;
    foreach my $d (split /,/, $cnt) {
      $dep += $d;
    }
    my $frq = (split /,/, $per)[0];
    if ($dep >= $para_dmin && $frq >= $para_fmin && $frq <= $para_fmax) {
      print join("\t", $e[0], $e[1], sprintf("%.3f", $frq)), "\n";
    }
  }
}
