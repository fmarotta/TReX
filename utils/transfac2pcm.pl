#!/usr/bin/env perl
use warnings;
# This program converts PCM in transfac format in the PCM format accepted by matrix_rider/vcf_rider
# with counts and no zeroes. We convert 0 to 1 directly, i.e. we do not add a pseudocount to all entries.
if ($#ARGV != -1) {
	die "Usage: $0 < transfac_pcm  > pcm_mrformat_add1\n";
}

my $status = 0; # Waiting to see an id.
my $id = "";
my @counts = ();
while(<STDIN>) {
	chomp $_;
	my @line = split(/\s+/, $_);
	if ($line[0] eq 'ID') {
		if ($status == 0) {
			$id = $line[1];
			$status = 1; # Waiting to see P0.
		} else {
			die "Malformed line at $., I saw too ID lines without counts after the first one!\n$_\n";
		}
	} elsif ($status == 1 && $line[0] eq 'P0') {
		$status = 2; # Waiting to get counts.
		if ($line[1] ne 'A' || $line[2] ne 'C' || $line[3] ne 'G' || $line[4] ne 'T') {
			die "Saw at line $. a not expected nucleotide order!\n$_\n";
		}
	} elsif ($status == 2) {
		if ($line[0] =~ /\d+/) {
			for (my $i = 1; $i <= 4; $i++) {
				if ($line[$i] != 0) {
					push(@counts, $line[$i]);
				} else {
					push(@counts, 1);
				}
			}
		} elsif ($line[0] eq "XX") {
			if (&print_pcm($id, \@counts)) {
				die "Malformed line at $., the previous pcm had something off!\n$_\n!";
			}
			$status = 0;
			$id = "";
			@counts = ();
		} else {
			die "Malformed line at $., I was expecting counts and saw a line with no counts or XX\n$_\n!";
		}
	}
	#print $status . "\t" . $id . "\n";
	# Not all erroneous combos of status and seen line are checked but should be sufficient.
}

sub print_pcm {
	my $id = shift;
	my $r_counts = shift;
	my $m_len = scalar( @{ $r_counts } );
	if ($id eq '' or $m_len % 4 != 0) {
		print STDERR $id . "\t" . $m_len . "\n";
		return 1;
	} else {
		$m_len = $m_len / 4;
		my $j = 0;
		# $j is the index on the flattened array, $i the position inside the pwm and $k the nucleotide.
		for (my $i = 0; $i < $m_len; $i++) {
			print $id . "\t" . $i;
			for (my $k = 0; $k < 4; $k++) {
				print "\t" . $r_counts->[$j];
				$j++;
			}
			print "\n";
		}
		return 0;
	}
}
