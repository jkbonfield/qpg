#!/usr/bin/perl -w
use strict;

# Estimates the quality of a candidate solution by analysing the bam alignment stats
# Usage: estimate_qual.pl run/seq_base ...
# eg: estimate_qual.pl seq_1050-0003-#1#1

# Output is the estimated quality and the actual quality (eval_cons)

foreach my $base (@ARGV) {
    foreach my $fn (sort glob("$base.bam.*")) {
	# Remap again?
	# We're evaluating the *path_cons* files, so we remap the reads to
	# this to get BAM stats from.
	#
	# This slows down the estimator considerably, but it ups correlation
	# from R=0.919 to 0.926
#	`minimap2 -x lr:hq -a $base.path_cons.0.0 $base.shred.fq | samtools sort -o $base.bam2.0.0 2>/dev/null`;
#	$fn = "$base.bam2.0.0";
	
	# Read mapping stats
	my ($nmapped,$mapped,$supp);
	open(FH, "samtools flagstat $fn |") || die;
	while (<FH>) {
	    chomp($_);
	    my @F = split(/\s+/, $_);
	    $mapped  = $F[5] if /primary mapped/;
	    $nmapped = $F[0] if /primary mapped/;
	    $supp    = $F[0] if /supplementary/;
	}
	close(FH);

	# Read mapping stats
	my ($mismatch, $err, $indel) = (0,0,0);
	open(FH, "samtools stats $fn|") || die;
	while (<FH>) {
	    chomp($_);
	    my @F = split(/\s+/, $_);
	    $mismatch = $F[2] if /^SN\tmismatches:/;
	    $err      = $F[3] if /^SN\terror rate:/;

	    # May not need as indel is included in NM which provides "mismatch"
	    $indel   += $F[1]*$F[2] if /^ID\t/;
	    $indel   += $F[1]*$F[3] if /^ID\t/;
	}
	close(FH);

	# Read depth variation
	open(FH, "samtools depth -a $fn|") || die;
	my $sum = 0;
	my $sum_sq = 0;
	my $count = 0;
	while (<FH>) {
	    chomp($_);
	    my @F = split(/\s+/, $_);
	    $sum    += $F[-1];
	    $sum_sq += $F[-1]*$F[-1];
	    $count++;
	}
	my $var = $sum_sq/$count - ($sum/$count)*($sum/$count);
	close(FH);

	my $path_seq = $fn;
	$path_seq =~ s/\.bam\d*\./.path_seq./;
	my $path_cons = $fn;
	$path_cons =~ s/\.bam\d+\./.path_cons./;

	my $seq_len  = `wc -c < $path_seq`;
	my $cons_len = `wc -c < $path_cons`;
	my $len_diff = abs($seq_len - $cons_len);
	my $indel_frac = $indel / $seq_len;

	# Optimiser stats (diff between seq and cons)
	my $opts="-xlr:hq -N0 -m100 -O 2,10 -E 2,1 -B 2";
	open(FH, "minimap2 $opts -a $path_seq $path_cons 2>/dev/null|") || die;
	my ($nalign,$nm) = (0,0);
	while (<FH>) {
	    next if (/^@/);
	    $nalign++;
	    $nm += $1 if (/\tNM:i:(\d+)/);
	}
	close(FH);

	# Measured eval_cons stats
	my $supp_perc = int($supp/($supp+$nmapped) * 100 * 1000 + 0.5)/1000;
	$mapped=~tr/%(//d;
	$fn =~ s/bam\d*/eval_cons/;
	open(FH, "<$fn") || die "$fn";
	my ($cov,$used,$contigs,$breaks);
	while (<FH>) {
	    chomp($_);
	    next if /contig/;
	    s/%//g;
	    my @F = split("\t", $_);
	    $cov     = $F[3];
	    $used    = $F[4];
	    $contigs = $F[5];
	    $breaks  = $F[6];
	}
	my $f1 = int(2*$cov*$used/($cov+$used)*1000+.5)/1000;
	close(FH);

	# An empirically derived formula for a score that's correlated against
	# the harmonic mean of used% and covered%.
	my $score = $mapped**1.4 * 0.1614;
	$score -=   $supp_perc**0.5 * 2.74;
	$score -=   $nalign**0.7 * 2.22;
	$score -=   $nm/100;
	$score -=   $var/30;
	$score = 0.1 if $score < 0.1;
	$score = $score**0.9 + 42;
	$score = 100 if $score > 100;

	#        0   1        2           3     X4    X5     X6   X7        X8       9        10   11
	print "$fn\t$mapped\t$supp_perc\t$var\t$cov\t$used\t$f1\t$contigs\t$breaks\t$nalign\t$nm\t$score\n";
    }
}

__END__;

$mapped -2*$supp_perc;
/tmp/_km  	R = 0.879078390806322
/tmp/_mg  	R = 0.763659182787955
/tmp/_ga  	R = 0.796459953546781
/tmp/_vg  	R = 0.825282663698217
/tmp/_km2 	R = 0.786158839376276
/tmp/_mg2 	R = 0.71195936773364
/tmp/_ga2 	R = 0.78127604750198
/tmp/_vg2 	R = 0.782713334690815
/tmp/_mqkm	R = 0.759969892742108
/tmp/_mqmg	R = 0.700206215557556
/tmp/_mqga	R = 0.67581079108933
/tmp/_mqvg	R = 0.740773274738386
0.766946 0.591056

$mapped               - 2*($nalign-1) - $nm/30;                                   
/tmp/_km  	R = 0.778651293344621
/tmp/_mg  	R = 0.720748050442272
/tmp/_ga  	R = 0.753433877065854
/tmp/_vg  	R = 0.782056208054651
/tmp/_km2 	R = 0.692740432336673
/tmp/_mg2 	R = 0.699916498839288
/tmp/_ga2 	R = 0.737483730233721
/tmp/_vg2 	R = 0.747848399970645
/tmp/_mqkm	R = 0.713149592857818
/tmp/_mqmg	R = 0.746215030076823
/tmp/_mqga	R = 0.789657981373488
/tmp/_mqvg	R = 0.763899000689913
0.743817 0.554209

$mapped -2*$supp_perc - 2*($nalign-1) - $nm/30;                                   
/tmp/_km  	R = 0.80244581549214
/tmp/_mg  	R = 0.769225924035529
/tmp/_ga  	R = 0.801364507131914
/tmp/_vg  	R = 0.819135011938861
/tmp/_km2 	R = 0.755490972559358
/tmp/_mg2 	R = 0.757003672092635
/tmp/_ga2 	R = 0.791625548746423
/tmp/_vg2 	R = 0.795115100436295
/tmp/_mqkm	R = 0.734564003290637
/tmp/_mqmg	R = 0.760651775972985
/tmp/_mqga	R = 0.804338080574051
/tmp/_mqvg	R = 0.791610325069013
0.781881 0.61194

$mapped -2*$supp_perc - 2*($nalign-1) - $nm/100;
/tmp/_km  	R = 0.885280606326711
/tmp/_mg  	R = 0.822883119165361
/tmp/_ga  	R = 0.839540039279067
/tmp/_vg  	R = 0.856027833132445
/tmp/_km2 	R = 0.814063775025119
/tmp/_mg2 	R = 0.798790906997559
/tmp/_ga2 	R = 0.832385424833525
/tmp/_vg2 	R = 0.840337012957649
/tmp/_mqkm	R = 0.79539155490957
/tmp/_mqmg	R = 0.795872447645792
/tmp/_mqga	R = 0.829562396765195
/tmp/_mqvg	R = 0.857489095201401
0.830635 0.690649  <<< best

$mapped -1.5*$supp_perc - 2*($nalign-1) - $nm/100;
0.827064 0.684701

$mapped -2.5*$supp_perc - 2*($nalign-1) - $nm/100;
0.823328 0.678603

$mapped -2*$supp_perc - 3*($nalign-1) - $nm/100;
0.816087 0.666688  (neg on all)

$mapped -2*$supp_perc - 1*($nalign-1) - $nm/100;
/tmp/_km  	R = 0.890148930721132 +
/tmp/_mg  	R = 0.816936603731278 -
/tmp/_ga  	R = 0.83980236220396  ~
/tmp/_vg  	R = 0.859242589463312 ~
/tmp/_km2 	R = 0.813462255247315 ~
/tmp/_mg2 	R = 0.788871703418606 -
/tmp/_ga2 	R = 0.832261318845892 ~
/tmp/_vg2 	R = 0.835285300171234 -
/tmp/_mqkm	R = 0.792075660685925 -
/tmp/_mqmg	R = 0.780831205247754 --
/tmp/_mqga	R = 0.806390313095143 --
/tmp/_mqvg	R = 0.840012701283525 --
0.82461 0.68089

$mapped -2*$supp_perc - 1.5*($nalign-1) - $nm/100;
/tmp/_km  	R = 0.889294138676572 +
/tmp/_mg  	R = 0.823487928733569 +
/tmp/_ga  	R = 0.841174028470496 +
/tmp/_vg  	R = 0.858403821880511 +
/tmp/_km2 	R = 0.815972394618732 +
/tmp/_mg2 	R = 0.797867382897333 ~
/tmp/_ga2 	R = 0.833840994417062 ~
/tmp/_vg2 	R = 0.840090691867397 ~
/tmp/_mqkm	R = 0.795881242405055 ~
/tmp/_mqmg	R = 0.792372344356714 ~
/tmp/_mqga	R = 0.806478958600503 --
/tmp/_mqvg	R = 0.840890879543538 --
0.82798 0.686303


$mapped -2*$supp_perc - 2*($nalign-1) - $nm/200;
/tmp/_km  	R = 0.886260245704367 +
/tmp/_mg  	R = 0.821740124064794 -
/tmp/_ga  	R = 0.836864456673975 -
/tmp/_vg  	R = 0.852698091741069 -
/tmp/_km2 	R = 0.809847567760025 --
/tmp/_mg2 	R = 0.793003104178669 -
/tmp/_ga2 	R = 0.828946507610176 -
/tmp/_vg2 	R = 0.835500314978508 -
/tmp/_mqkm	R = 0.7923618928533 -
/tmp/_mqmg	R = 0.789212581973958 -
/tmp/_mqga	R = 0.746829176580446 ---
/tmp/_mqvg	R = 0.797642178316014 ---
0.815909 0.666896

$mapped -2*$supp_perc - 2*($nalign-1) - $nm/75;
/tmp/_km  	R = 0.880242550034959 -
/tmp/_mg  	R = 0.819952119413151 -
/tmp/_ga  	R = 0.837335462714399 -
/tmp/_vg  	R = 0.85419753252847  -
/tmp/_km2 	R = 0.811499161480452 -
/tmp/_mg2 	R = 0.798275114123371 -
/tmp/_ga2 	R = 0.830261612491023 -
/tmp/_vg2 	R = 0.838624106953412 -
/tmp/_mqkm	R = 0.792356613014944 -
/tmp/_mqmg	R = 0.79624885378827 +
/tmp/_mqga	R = 0.818623422697463 -
/tmp/_mqvg	R = 0.848921339544384 --
0.827211 0.684915

$mapped -2*$supp_perc - 2*($nalign-1) - $nm/125;
/tmp/_km  	R = 0.886670260493407 +
/tmp/_mg  	R = 0.82330370739173  +
/tmp/_ga  	R = 0.839510847915176 ~
/tmp/_vg  	R = 0.855746076998629 ~
/tmp/_km2 	R = 0.813690232213426 -
/tmp/_mg2 	R = 0.797521251445785 -
/tmp/_ga2 	R = 0.832154705418196 ~
/tmp/_vg2 	R = 0.83964519252817  -
/tmp/_mqkm	R = 0.79540479494367  ~
/tmp/_mqmg	R = 0.794169175736228 -
/tmp/_mqga	R = 0.786658660561881 --
/tmp/_mqvg	R = 0.82812876664509  --
0.824384 0.680393


=============================================================================

Final:

vs eval_seq 
$mapped -2*$supp_perc - 2*($nalign-1) - $nm/100;
/tmp/_km  	R = 0.885280606326711
/tmp/_mg  	R = 0.822883119165361
/tmp/_ga  	R = 0.839540039279067
/tmp/_vg  	R = 0.856027833132445
/tmp/_km2 	R = 0.814063775025119
/tmp/_mg2 	R = 0.798790906997559
/tmp/_ga2 	R = 0.832385424833525
/tmp/_vg2 	R = 0.840337012957649
/tmp/_mqkm	R = 0.79539155490957
/tmp/_mqmg	R = 0.795872447645792
/tmp/_mqga	R = 0.829562396765195
/tmp/_mqvg	R = 0.857489095201401
0.830635 0.690649

vs eval_cons
$mapped -2*$supp_perc - 2*($nalign-1) - $nm/100;
/tmp/_km  	R = 0.872029576080218
/tmp/_mg  	R = 0.825288639329911
/tmp/_ga  	R = 0.821629073029486
/tmp/_vg  	R = 0.841679596516064
/tmp/_km2 	R = 0.787710539919411
/tmp/_mg2 	R = 0.799846346113465
/tmp/_ga2 	R = 0.812023764631774
/tmp/_vg2 	R = 0.821796914760769
/tmp/_mqkm	R = 0.768170830907575
/tmp/_mqmg	R = 0.78783346783834
/tmp/_mqga	R = 0.789848113570591
/tmp/_mqvg	R = 0.831884380828404
0.813312 0.662223

vs eval_cons, revised minimap2 options
/tmp/_km  	R = 0.876156800172059
/tmp/_mg  	R = 0.82052795177086
/tmp/_ga  	R = 0.83063492712825
/tmp/_vg  	R = 0.85656581906909
/tmp/_km2 	R = 0.787188123805958
/tmp/_mg2 	R = 0.79221769934815
/tmp/_ga2 	R = 0.822651494703104
/tmp/_vg2 	R = 0.828804710970885
/tmp/_mqkm	R = 0.768416987294951
/tmp/_mqmg	R = 0.775518006380906
/tmp/_mqga	R = 0.796354702682633
/tmp/_mqvg	R = 0.834704273577235
0.815812 0.666527


-$var/50
/tmp/_km  	R = 0.920241233358204++
/tmp/_mg  	R = 0.907016811132933+++
/tmp/_ga  	R = 0.893857154722622+++
/tmp/_vg  	R = 0.912667543890098+++
/tmp/_km2 	R = 0.880657454403117++
/tmp/_mg2 	R = 0.893459676679665+++
/tmp/_ga2 	R = 0.886615459567163++
/tmp/_vg2 	R = 0.897984501533859++
/tmp/_mqkm	R = 0.860119186193415+++
/tmp/_mqmg	R = 0.880892332947483+++
/tmp/_mqga	R = 0.800489785082632+
/tmp/_mqvg	R = 0.837359748568438+
0.880947 0.777124

-$var/20
0.866062 0.753132

-$var/60
0.877426 0.770772

-$var/40
/tmp/_km  	R = 0.920120040014206
/tmp/_mg  	R = 0.90978229210557
/tmp/_ga  	R = 0.898128818922253
/tmp/_vg  	R = 0.915931329118107
/tmp/_km2 	R = 0.887497428371653
/tmp/_mg2 	R = 0.899960589482393
/tmp/_ga2 	R = 0.890869858276209
/tmp/_vg2 	R = 0.90286515117268
/tmp/_mqkm	R = 0.867323082577373
/tmp/_mqmg	R = 0.885984888830719
/tmp/_mqga	R = 0.789783454575931
/tmp/_mqvg	R = 0.828283318308358
0.883044 0.78111

-$var/30
0.880784 0.777671

-var/40 score**0.9
/tmp/_km  	R = 0.922220258257691+
/tmp/_mg  	R = 0.909522760907347~
/tmp/_ga  	R = 0.897118511690500~
/tmp/_vg  	R = 0.915710381934181~
/tmp/_km2 	R = 0.887100357224534~
/tmp/_mg2 	R = 0.899856232385883~
/tmp/_ga2 	R = 0.889965816681326-
/tmp/_vg2 	R = 0.902365354591502~
/tmp/_mqkm	R = 0.866327198600395-
/tmp/_mqmg	R = 0.884707631768855-
/tmp/_mqga	R = 0.786150547244411-
/tmp/_mqvg	R = 0.822305387013423-
0.881946 0.779292

-var/40, score**0.5
/tmp/_km  	R = 0.928075399127473++
/tmp/_mg  	R = 0.878565625647535-
/tmp/_ga  	R = 0.890996962399642-
/tmp/_vg  	R = 0.911475907192114-
/tmp/_km2 	R = 0.883474205571470~
/tmp/_mg2 	R = 0.895549629528005-
/tmp/_ga2 	R = 0.884183003880789-
/tmp/_vg2 	R = 0.898449129478462-
/tmp/_mqkm	R = 0.859599451937629-
/tmp/_mqmg	R = 0.871662265391011--
/tmp/_mqga	R = 0.766979506221110--
/tmp/_mqvg	R = 0.777516655339857---
0.870544 0.760074


sqrt(var)/1.5 (vs *1)
/tmp/_km  	R = 0.91315797636784
/tmp/_mg  	R = 0.897798840170256
/tmp/_ga  	R = 0.881940677392581
/tmp/_vg  	R = 0.902873236430702
/tmp/_km2 	R = 0.869325151777437
/tmp/_mg2 	R = 0.882746084653387
/tmp/_ga2 	R = 0.873970948302093
/tmp/_vg2 	R = 0.884584998668781
/tmp/_mqkm	R = 0.852013384691823
/tmp/_mqmg	R = 0.86930081574845
/tmp/_mqga	R = 0.803969865580464
/tmp/_mqvg	R = 0.843527598403691
0.872934 0.762808

sqrt(var)  vs var/40
/tmp/_km  	R = 0.912847549145163--
/tmp/_mg  	R = 0.908582610707660-
/tmp/_ga  	R = 0.888513264954981-
/tmp/_vg  	R = 0.908069210789187--
/tmp/_km2 	R = 0.881908676502458-
/tmp/_mg2 	R = 0.897039843843536~
/tmp/_ga2 	R = 0.880121833241239-
/tmp/_vg2 	R = 0.891824468647944-
/tmp/_mqkm	R = 0.865690856104352--
/tmp/_mqmg	R = 0.884030123567996-
/tmp/_mqga	R = 0.791497517099595+
/tmp/_mqvg	R = 0.832879080279079+
0.878584 0.773029

sqrt(var, score**0.5
/tmp/_km  	R = 0.925430952946734++
/tmp/_mg  	R = 0.907013526136008-
/tmp/_ga  	R = 0.883161146428275-
/tmp/_vg  	R = 0.905471292056314-
/tmp/_km2 	R = 0.879370161261689-
/tmp/_mg2 	R = 0.897254153711969+
/tmp/_ga2 	R = 0.875546629478285-
/tmp/_vg2 	R = 0.889814030845963-
/tmp/_mqkm	R = 0.858733742924395-
/tmp/_mqmg	R = 0.880042955429537-
/tmp/_mqga	R = 0.775693172187874--
/tmp/_mqvg	R = 0.777904238736052---
0.871286 0.761199

sqrt(var, score**1.2
/tmp/_km  	R = 0.905299500924453
/tmp/_mg  	R = 0.905803155853624
/tmp/_ga  	R = 0.888497338064308
/tmp/_vg  	R = 0.905594058452952
/tmp/_km2 	R = 0.880934634945166
/tmp/_mg2 	R = 0.894360882657123
/tmp/_ga2 	R = 0.879748430420042
/tmp/_vg2 	R = 0.890689863951212
/tmp/_mqkm	R = 0.865680453947826
/tmp/_mqmg	R = 0.882031264316343
/tmp/_mqga	R = 0.795206289054026
/tmp/_mqvg	R = 0.837820699816309
0.877639 0.771205

sqrt(var)*1.25 vs var/40
/tmp/_km  	R = 0.908498387726912--
/tmp/_mg  	R = 0.910724886212387~
/tmp/_ga  	R = 0.889719533835701-
/tmp/_vg  	R = 0.908286394218983-
/tmp/_km2 	R = 0.885505318449292~
/tmp/_mg2 	R = 0.901173154121354+
/tmp/_ga2 	R = 0.880914219816004-
/tmp/_vg2 	R = 0.893106128339798-
/tmp/_mqkm	R = 0.869908030363461~
/tmp/_mqmg	R = 0.888095595242678~
/tmp/_mqga	R = 0.779309258662264-
/tmp/_mqvg	R = 0.821953480521084~
0.8781 0.77247


# final (for now)
score = $mapped -2*$supp_perc - 2*($nalign-1) - $nm/100 - $var/40;
clipped with >= 1 <= 100 (max 100 helps R for F1 percentage with some funcs)
/tmp/_km  	R = 0.920120040014206
/tmp/_mg  	R = 0.910931097329467
/tmp/_ga  	R = 0.898128818922253
/tmp/_vg  	R = 0.915931329118107
/tmp/_km2 	R = 0.887497428371653
/tmp/_mg2 	R = 0.899960589482393
/tmp/_ga2 	R = 0.890869858276209
/tmp/_vg2 	R = 0.90286515117268
/tmp/_mqkm	R = 0.867323082577373
/tmp/_mqmg	R = 0.885984888830719
/tmp/_mqga	R = 0.789783454575931
/tmp/_mqvg	R = 0.828283318308358
0.88314 0.781284 <<<
