#!/usr/bin/perl -w

use strict;
use DBI;
#use Bio::Tools::Run::Ensembl;

my (@liste_gene, $raw_data, $treated, $taille);

#my $nb_nt = 0;
if (!$ARGV[0]) {
	print "\n\nProgramme pour rechercher les coordonnees hg19 des exons pour un ou des genes donnes en argument\n\nUsage:\n\n$./ngs_designer.pl nom_gene\n\nou en mode batch\n\n$./ngs_designer.pl nom_fichier.txt\n\navec nom_fichier.txt un fichier contenant une liste de genes (nom HGNC)\n\nMettre le fichier dans le meme repertoire que le programme\n\nPour modifier le nombre de nucleotides a ajouter avant et apres chaque exon (25 par defaut), le specifier en deuxieme argument; ex:\n\n$./ngs_designer.pl nom_gene 300\n\n";
	exit;
}
elsif ($ARGV[0]) {
	if ($ARGV[0] =~ /^(\w+\.txt)$/o) {
		open F, $1;
		while (<F>) {
			/^(\w+)$/o;
			push @liste_gene, $1;
		}
		close F;
	}
	elsif ($ARGV[0] =~ /(\w+)/o) {
		push @liste_gene, $1;
	}
	else {exit}
}
my $ajout = 25;
if ($ARGV[1] && $ARGV[1] =~ /(\d+)/o) {$ajout = $1}

my $dbh = DBI->connect(    "DBI:mysql:database=hg19;host=genome-mysql.cse.ucsc.edu;",
                        "genome",
                        "",
                        {'RaiseError' => 1}
                ) or die $!;




my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date =  $mday."_".($mon+1)."_".(1900+$year);


open I, ">Nimblegen_".$date."_$ajout.txt";
open J, ">Illumina_".$date."_$ajout.csv";
print J "Chromosome,Start,Stop,TargetType,DesiredprobeSpacing,Include3UTR,Include5UTR,Labels\n";
open K, ">Asked_".$date."_$ajout.bed";
print K "track name=\"Design Capture $ARGV[0] $ajout\" description=\"Asked design $ajout\" visibility=1\n";
open L, ">Agilent_".$date."_$ajout.txt";


my $num_region = 1;

foreach (@liste_gene) {
	my ($raw_data, $treated, @exon_s, @exon_e, $long_t, $chr, $s, @temp, $elu, $nb_ex, $ex_elu);
	my $nom = $_;
	my $query = "SELECT name, chrom, strand, txStart, txEnd, exonCount, exonStarts, exonEnds FROM refGene WHERE name2 = '$nom';";# AND cdsStartStat != 'none';";# AND name like 'NM_%';";
	my $sth = $dbh->prepare($query);
	my $res = $sth->execute();
	if ($res ne '0E0') {
		while (my $result = $sth->fetchrow_hashref()) {
			if ($result->{'chrom'} !~ /^chr[\dX]+$/og) {next}
			$raw_data->{$result->{'name'}} = [$result->{'chrom'}, $result->{'strand'}, $result->{'txStart'}+1, $result->{'txEnd'}, $result->{'exonCount'}, $result->{'exonStarts'}, $result->{'exonEnds'}];
			($chr, $s) = ($result->{'chrom'}, $result->{'strand'});
		}
		$ex_elu = 0;

		foreach my $keys (keys(%{$raw_data})) {
			#undef @temp;
			#undef @exon_s;
			#undef @exon_e;
			##print $raw_data->{$keys}->[5]."\n";
			#@temp = split(/,/, $raw_data->{$keys}->[5]);
			##print $raw_data->{$keys}->[5];exit;
			#foreach (@temp) {
			#	if ($raw_data->{$keys}->[1] eq '+') {push @exon_s, $_+1}
			#	else {unshift @exon_s, $_+1}
			#}
			#@temp = split(/,/, $raw_data->{$keys}->[6]);
			#foreach (@temp) {
			#	if ($raw_data->{$keys}->[1] eq '+') {push @exon_e, $_}
			#	else {unshift @exon_e, $_}
			#}
			#calcul nb exon puis taille transcrit puis $elu
			#$elu => pas NR_
			if ($keys !~ /NR/) {
				$nb_ex = $raw_data->{$keys}->[4];
				#print $nb_ex;exit;
				if ($nb_ex > $ex_elu) {$elu = $keys;$ex_elu = $nb_ex;}
				elsif ($nb_ex == $ex_elu) {
					my $long_ancien = $raw_data->{$elu}->[3] - $raw_data->{$elu}->[2];
					my $long_nouveau = $raw_data->{$keys}->[3] - $raw_data->{$keys}->[2];
					if ($long_ancien > $long_nouveau) {$elu = $keys;}
				}
			}
		}
		undef @temp;
		undef @exon_s;
		undef @exon_e;
		if ($nom eq 'MYO7A') {$elu = "NM_000260"}
		elsif ($nom eq 'PCDH15') {$elu = "NM_033056"}
		elsif ($nom eq 'DFNB31') {$elu = "NM_015404"}
		elsif ($nom eq 'VEZT') {$elu = "NM_017599"}
		elsif ($nom eq 'WFS1') {$elu = "NM_006005"}
		@temp = split(/,/, $raw_data->{$elu}->[5]);
		foreach (@temp) {
			if ($raw_data->{$elu}->[1] eq '+') {push @exon_s, $_+1}
			else {unshift @exon_s, $_+1}
		}
		@temp = split(/,/, $raw_data->{$elu}->[6]);
		foreach (@temp) {
			if ($raw_data->{$elu}->[1] eq '+') {push @exon_e, $_}
			else {unshift @exon_e, $_}
		}
		for (my $i = 1; $i <= $#temp+1; $i++) {
			$treated->{$elu}->{$i} = [shift @exon_s, shift @exon_e];
		}

		foreach my $keys(keys %{$treated}) {
			foreach my $exon (sort {$a <=> $b} keys %{$treated->{$keys}}) {
				my ($start, $end);
				if ($s eq '+') {
					my $taille_intron5 = $ajout*2+1;
					my $taille_intron3 = $ajout*2+1;
					if ($exon > 1 && $exon < $#temp+1) {
						$taille_intron5 = $treated->{$keys}->{$exon}->[0]-$treated->{$keys}->{$exon-1}->[1]-1;
						$taille_intron3 = $treated->{$keys}->{$exon+1}->[0]-$treated->{$keys}->{$exon}->[1]-1;
					}
					elsif ($exon == 1 && $#temp+1 > 1) {$taille_intron3 = $treated->{$keys}->{$exon+1}->[0]-$treated->{$keys}->{$exon}->[1]}
					elsif ($exon == 1 && $#temp+1 == 1) {$taille_intron3 = 0}
					elsif ($exon == $#temp+1) {$taille_intron5 = $treated->{$keys}->{$exon}->[0]-$treated->{$keys}->{$exon-1}->[1]}

					if ($taille_intron5 > $ajout*2 && $taille_intron3 > $ajout*2) {
						($start, $end) = ($treated->{$keys}->{$exon}->[0]-$ajout, $treated->{$keys}->{$exon}->[1]+$ajout);
					}
					elsif ($taille_intron5 <= $ajout*2 && $taille_intron3 <= $ajout*2) {
						($start, $end) = ($treated->{$keys}->{$exon}->[0]-(sprintf("%.0f", $taille_intron5/2)), $treated->{$keys}->{$exon}->[1]+(sprintf("%.0f", $taille_intron3/2)));
					}
					elsif ($taille_intron5 <= $ajout*2) {
						($start, $end) = ($treated->{$keys}->{$exon}->[0]-(sprintf("%.0f", $taille_intron5/2)), $treated->{$keys}->{$exon}->[1]+$ajout);
					}
					elsif ($taille_intron3 <= $ajout*2) {
						($start, $end) = ($treated->{$keys}->{$exon}->[0]-$ajout, $treated->{$keys}->{$exon}->[1]+(sprintf("%.0f", $taille_intron3/2)));
					}

				}
				else {
					my $taille_intron5 = $ajout*2+1;
					my $taille_intron3 = $ajout*2+1;
					if ($exon > 1 && $exon < $#temp+1) {
						$taille_intron5 = $treated->{$keys}->{$exon-1}->[0]-$treated->{$keys}->{$exon}->[1]-1;
						$taille_intron3 = $treated->{$keys}->{$exon}->[0]-$treated->{$keys}->{$exon+1}->[1]-1;
					}
					elsif ($exon == 1 && $#temp+1 > 1) {$taille_intron3 = $treated->{$keys}->{$exon}->[0]-$treated->{$keys}->{$exon+1}->[1]-1}
					elsif ($exon == 1 && $#temp+1 == 1) {$taille_intron3 = 0}
					elsif ($exon == $#temp+1) {$taille_intron5 = $treated->{$keys}->{$exon-1}->[0]-$treated->{$keys}->{$exon}->[1]-1}

					if ($taille_intron5 > $ajout*2 && $taille_intron3 > $ajout*2) {
						($start, $end) = ($treated->{$keys}->{$exon}->[0]-$ajout, $treated->{$keys}->{$exon}->[1]+$ajout);
					}
					elsif ($taille_intron5 <= $ajout*2 && $taille_intron3 <= $ajout*2) {
						($start, $end) = ($treated->{$keys}->{$exon}->[0]-(sprintf("%.0f", $taille_intron3/2)), $treated->{$keys}->{$exon}->[1]+(sprintf("%.0f", $taille_intron5/2)));
					}
					elsif ($taille_intron5 <= $ajout*2) {
						($start, $end) = ($treated->{$keys}->{$exon}->[0]-$ajout, $treated->{$keys}->{$exon}->[1]+(sprintf("%.0f", $taille_intron5/2)));
					}
					elsif ($taille_intron3 <= $ajout*2) {
						($start, $end) = ($treated->{$keys}->{$exon}->[0]-(sprintf("%.0f", $taille_intron3/2)), $treated->{$keys}->{$exon}->[1]+$ajout);
					}
					if (!$start) {
						print "\n\n$taille_intron5\t$taille_intron3\t$exon\n\n"
					}

				}
				my $curr_size = $end-$start+1;
				&write_file($num_region, $nom, $keys, $exon, $start, $end, $curr_size, $chr, 'MAIN');
				$taille += $end-$start+1;
				$num_region++;
			}
		}




		#we add vega genes in addition to refgene works but chiant
		#my $query_vega = "SELECT name, chrom, strand, txStart, txEnd, exonCount, exonStarts, exonEnds FROM vegaGene WHERE name2 = '$nom' AND cdsStartStat != 'none';";# AND name like 'NM_%';";
		#my $sth_vega = $dbh->prepare($query_vega);
		#my $res_vega = $sth_vega->execute();
		#if ($res_vega ne '0E0') {
		#	while (my $result_vega = $sth_vega->fetchrow_hashref()) {
		#		if ($result_vega->{'chrom'} !~ /^chr[\dX]+$/og) {next}
		#		$raw_data->{$result_vega->{'name'}} = [$result_vega->{'chrom'}, $result_vega->{'strand'}, $result_vega->{'txStart'}+1, $result_vega->{'txEnd'}, $result_vega->{'exonCount'}, $result_vega->{'exonStarts'}, $result_vega->{'exonEnds'}];
		#	}
		#}
		#my $query_ens = "SELECT a.name, chrom, strand, txStart, txEnd, exonCount, exonStarts, exonEnds FROM ensGene a, ensemblToGeneName b WHERE a.name = b.name AND b.value = '$nom';";# AND cdsStartStat != 'none';";# AND name like 'NM_%';";
		#my $sth_ens = $dbh->prepare($query_ens);
		#my $res_ens = $sth_ens->execute();
		#if ($res_ens ne '0E0') {
		#	while (my $result_ens = $sth_ens->fetchrow_hashref()) {
		#		if ($result_ens->{'chrom'} !~ /^chr[\dX]+$/og) {next}
		#		$raw_data->{$result_ens->{'name'}} = [$result_ens->{'chrom'}, $result_ens->{'strand'}, $result_ens->{'txStart'}+1, $result_ens->{'txEnd'}, $result_ens->{'exonCount'}, $result_ens->{'exonStarts'}, $result_ens->{'exonEnds'}];
		#	}
		#}

		###COMMENT BELOW TO REMOVE ALT TRANSCRIPT
		my $retenu = '';
		foreach my $keys (keys(%{$raw_data})) {
			if ($keys eq $elu) {next}
			undef @temp;
			undef @exon_s;
			undef @exon_e;
			@temp = split(/,/, $raw_data->{$keys}->[5]);
			foreach (@temp) {
				if ($raw_data->{$keys}->[1] eq '+') {push @exon_s, $_+1}
				else {unshift @exon_s, $_+1 }
			}
			@temp = split(/,/, $raw_data->{$keys}->[6]);
			foreach (@temp) {
				if ($raw_data->{$keys}->[1] eq '+') {push @exon_e, $_}
				else {unshift @exon_e, $_}
			}
			#my $retenu = '';
			for (my $i = 0; $i <= $#exon_s ; $i++) {
				my ($start, $end);
				my $semaph = 0;
				foreach my $exon (sort {$a <=> $b} keys %{$treated->{$elu}}) {
					if ($exon eq '0') {print "\n\n$treated->{$elu}->{$exon}->[0]"}
					##select regions not already included we use a padding of 25nts
					#print "\n\n1: $exon_s[$i]-2: $treated->{$elu}->{$exon}->[0]-3:$exon_e[$i]-4:$treated->{$elu}->{$exon}->[1]\n\n";
					if (($exon_s[$i] == $treated->{$elu}->{$exon}->[0] && $exon_e[$i] == $treated->{$elu}->{$exon}->[1])) {$semaph = 1;last;}
					elsif ($exon_s[$i] >= $treated->{$elu}->{$exon}->[0]-25 && $exon_e[$i] <= $treated->{$elu}->{$exon}->[1]+25) {$semaph = 1;last;}
				}
				if ($semaph == 0) {
					if ($retenu !~ /.*$exon_s[$i]-$exon_e[$i].*/) {
						#my @selected = split(/_/, $retenu); useless
						#foreach (@selected) {
						#	/(\d+)-(\d+)/o;
						#	if ($exon_s[$i] == $1 && $exon_e[$i] == $2) {next}
						#	elsif ($exon_s[$i] >= $1-25 && $exon_e[$i] >= $2+25) {next}
						#}


						$retenu .= $exon_s[$i]."-".$exon_e[$i]."_";

						if ($s eq '+') {

							my $taille_intron5 = $ajout*2+1;
							my $taille_intron3 = $ajout*2+1;
							#print "\n$i - $#exon_s - $exon_s[$i] a $exon_e[$i-1] b $exon_s[$i+1] c $exon_e[$i]d\n";
							#print "\n$nom - $i - $retenu\n";
							if ($i > 0 && $i < $#exon_s) {
								$taille_intron5 = $exon_s[$i]-$exon_e[$i-1];
								$taille_intron3 = $exon_s[$i+1]-$exon_e[$i];

							}
							elsif ($i == 0) {$taille_intron3 = $exon_s[$i+1]-$exon_e[$i]}
							elsif ($i == $#exon_s) {$taille_intron5 = $exon_s[$i]-$exon_e[$i-1]}

							if ($taille_intron5 > $ajout*2 && $taille_intron3 > $ajout*2) {
								($start, $end) = ($exon_s[$i]-$ajout, $exon_e[$i]+$ajout);
							}
							elsif ($taille_intron5 <= $ajout*2 && $taille_intron3 <= $ajout*2) {
								($start, $end) = ($exon_s[$i]-(sprintf("%.0f", $taille_intron5/2)), $exon_e[$i]+(sprintf("%.0f", $taille_intron3/2)));
							}
							elsif ($taille_intron5 <= $ajout*2) {
								($start, $end) = ($exon_s[$i]-(sprintf("%.0f", $taille_intron5/2)), $exon_e[$i]+$ajout);
							}
							elsif ($taille_intron3 <= $ajout*2) {
								($start, $end) = ($exon_s[$i]-$ajout, $exon_e[$i]+(sprintf("%.0f", $taille_intron3/2)));
							}
						}
						else {
							my $taille_intron5 = $ajout*2+1;
							my $taille_intron3 = $ajout*2+1;
							if ($i > 0 && $i < $#exon_s) {
								$taille_intron5 = $exon_s[$i-1]-$exon_e[$i]-1;
								$taille_intron3 = $exon_s[$i]-$exon_e[$i+1]-1;
							}
							elsif ($i == 0) {$taille_intron3 = $exon_s[$i]-$exon_e[$i+1]-1}
							elsif ($i == $#exon_s) {$taille_intron5 = $exon_s[$i-1]-$exon_e[$i]-1}
							if ($taille_intron5 > $ajout*2 && $taille_intron3 > $ajout*2) {
								($start, $end) = ($exon_s[$i]-$ajout, $exon_e[$i]+$ajout);
							}
							elsif ($taille_intron5 <= $ajout*2 && $taille_intron3 <= $ajout*2) {
								($start, $end) = ($exon_s[$i]-(sprintf("%.0f", $taille_intron3/2)), $exon_e[$i]+(sprintf("%.0f", $taille_intron5/2)));
							}
							elsif ($taille_intron5 <= $ajout*2) {
								($start, $end) = ($exon_s[$i]-$ajout, $exon_e[$i]+(sprintf("%.0f", $taille_intron5/2)));
							}
							elsif ($taille_intron3 <= $ajout*2) {
								($start, $end) = ($exon_s[$i]-(sprintf("%.0f", $taille_intron3/2)), $exon_e[$i]+$ajout);
							}
						}
						my $curr_size = $end-$start+1;
						&write_file($num_region, $nom, $keys, $i+1, $start, $end, $curr_size, $chr, 'ALT');
						$taille += $end-$start+1;
						$num_region++;
						#$taille += ($exon_e[$i]-$exon_s[$i]+1) + (2*$ajout);
					}
				}

				#}
			}
			#if ($semaph == 1) {$class++;$semaph = 0;}
		}

		### END COMMENT BELOW

		print "\n\nTaille totale:".$taille."\n\n";



	}
	else {print "\n\nle nom $nom ne semble pas etre le nom HGNC\n\n"}


}
close I;
close J;
close K;
close L;

#new feature to sort and merge data into a single bed file 08/2016 - bedtools muste be in the PATH
my $bed = "Asked_".$date."_$ajout";
system "sort -k1,1 -k2,2n $bed.bed > $bed.sorted.bed; bedtools merge -i $bed.sorted.bed -c 4 -o collapse > $bed.sorted.merged.bed";


sub write_file {
	my ($num_region, $nom, $keys, $exon, $start, $end, $curr_size, $chr, $type) = @_;
	print I "$chr\t$start\t$end\t$nom:$keys-$exon\n";
	print J "$chr,$start,$end,FullRegion,Overlapping,FALSE,FALSE,$nom:$keys-$exon\n";
	print K "$chr\t$start\t$end\t$nom:$keys-$exon\n";
	print L "$chr:$start-$end\t$nom:$keys-$exon\n";
}

exit;
