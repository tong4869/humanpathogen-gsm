#!/usr/bin/perl
use strict;
use warnings;

if (@ARGV != 3) {
    die "Wrong: insufficient parameters";
}

my ($sample_name, $r1_file, $r2_file) = @ARGV;

my $r1_avg_len = calculate_avg_len($r1_file);
my $r2_avg_len = calculate_avg_len($r2_file);
my $r1_count = calculate_seq_count($r1_file);
my $r2_count = calculate_seq_count($r2_file);
my $stdnum = (( $r1_avg_len - 50 ) * $r1_count + ( $r2_avg_len - 50 ) * $r2_count) / (10**10);


print "$sample_name\t$r1_avg_len\t$r1_count\t$r2_avg_len\t$r2_count\t$stdnum\n";

sub calculate_avg_len {
    my ($filename) = @_;
    open(my $fh, '<', $filename) or die "Could not open file '$filename' $!";
    my $total_len = 0;
    my $count = 0;
    while (my $header = <$fh>) {
        my $seq = <$fh>;
        my $plus = <$fh>;
        my $qual = <$fh>;
        $total_len += length($seq) - 1; 
        $count++;
    }
    close $fh;
    return $count ? $total_len / $count : 0;
}

sub calculate_seq_count {
    my ($filename) = @_;
    open(my $fh, '<', $filename) or die "Could not open file '$filename' $!";
    my $count = 0;
    while (<$fh>) {
        $count++ if $. % 4 == 1; # 每4行计数一次
    }
    close $fh;
    return $count;
}
