# /usr/local/bin/perl -wT

####################
# All proteins MDM
####################

$sql= qq{ SELECT f.ensp_id, ensg_id, aa_start, aa_stop
					FROM feature f
					JOIN (SELECT ensp_id, COUNT(DISTINCT feature_type_id) num_type
								FROM feature
								WHERE segment IN ('t', 'r')
								GROUP BY ensp_id
								HAVING num_type>3) AS tm USING (ensp_id)
					WHERE segment IN ('t', 'r') AND feature_type_id NOT IN (3,4,7,11) };
$sth = $atlas_dbh->prepare($sql);
$sth->execute();

#Analyze one ensp_id/TM region at a time
while (($ensp_id, $ensg_id, $start, $stop) = $sth->fetchrow()) {
	$ensg_hash{$ensp_id} = $ensg_id;
	my $pos = $start;
	#Check Signal Peptide (SP) for positions < 50
	if ($start < 50) {
		if (exists $sp_hash{$ensp_id} && $sp_hash{$ensp_id} > $start) {
			#Skip TM region of overlap >= 5
			if ($sp_hash{$ensp_id}-$start >= 5) {
				next;
			}
			#Else start count from end of SP
			else {
				$pos = $sp_hash{$ensp_id} + 1;
			}
		}
	}
	while ($pos <= $stop) {
		if (exists $tm_hash{$ensp_id}{$pos}) {
			$tm_hash{$ensp_id}{$pos}++;
		}
		else {
			$tm_hash{$ensp_id}{$pos} = 1;
		}
		$pos++;
	}
}
$time=localtime;
print $time."  Finished hashing\n";

$num = 0;
while (($ensp_id, $ensg_id) = each %ensg_hash) {
#foreach $ensp_id (sort keys %ensg_hash) {
	$num++;
	if ($num % 500 == 0 ) {
		$time=localtime;
		print "$time \t $num ensp_ids processed\n";
	}

	$in_TM = 0;
	$last_pos = 0;
	$tm_start = 0;
	$tm_stop = 0;
	
	foreach $pos (sort {$a<=>$b} keys %{$tm_hash{$ensp_id}}) {
		#If less then 4 methods say TM, skip position
		if ($tm_hash{$ensp_id}{$pos} < 4) {
			#print "next\n";
			next;
		}
		
		#Save segment to database when entering a new segment
		if ($tm_start && ($pos-$last_pos) > 1) {
			&save_TM($atlas_dbh, $ensp_id, $ensg_id, $tm_start, $last_pos);
			$tm_start = 0;
		}
		if (!$tm_start) {
			$tm_start = $pos;
		}
		$last_pos = $pos;
	}
	if ($tm_start) {
		&save_TM($atlas_dbh, $ensp_id, $ensg_id, $tm_start, $last_pos);
	}
}

# create proteinclass

my $sub_source="MDM";
my $sub_code="Md";
my $sub_url="http://www.ncbi.nlm.nih.gov/pubmed/20175080";
my $sub_name='Membrane proteins predicted by MDM';
my $sub_description="Transmembrane proteins predictions based on a majority decision method";
my $sub_id = &insert_protein_class($sub_name, $sub_code, $public, $sub_description, $sub_source, $external_db, $institute, $sub_url, $reference, $external_version);
&set_parent($sub_id, $class_id);
$sql= qq{ 	SELECT distinct et.et_id, uniprot_id 
						FROM prest_atlas_db.feature tm, lims.ensembl_transcripts et, lims.ensembl_genes eg 
						WHERE et.eg_id=eg.eg_id AND et.ensp_id=tm.ensp_id AND feature_type_id=11 };
$sth = $lims_dbh->prepare($sql);
$sth->execute();
$total=0;
$num_et=0;
while (my ($et_id, $uniprot_id) = $sth->fetchrow_array()) {
	$total++;
	$num_et++;
	&insert_et_ids($et_id, $parent_class_id, $uniprot_id, '');
	&insert_et_ids($et_id, $class_id, $uniprot_id, '');
	&insert_et_ids($et_id, $sub_id, $uniprot_id, '');
}
print LOG "Total MDM $total\nTotal enst_ids found: $num_et\n";

######################
# All proteins MDSEC
######################

my $sub_source="HPA";
my $sub_code="Sa";
my $sub_url="http://www.proteinatlas.org/humanproteome/secretome#prediction";
my $sub_name='Secreted proteins predicted by MDSEC';
my $sub_description="Secreted protein predictions based on a majority decision method";
my $sub_id = &insert_protein_class($sub_name, $sub_code, $public, $sub_description, $sub_source, $external_db, $institute, $sub_url, $reference, $external_version);
&set_parent($sub_id, $class_id);

$sql= qq{ 	SELECT distinct et.et_id, uniprot_id, if(count(distinct f1.feature_type_id)>=1,count(distinct f1.feature_type_id),0) as SP
							FROM lims.ensembl_genes eg 
							JOIN lims.ensembl_transcripts et USING (eg_id) 
							JOIN prest_atlas_db.feature f1 USING (ensp_id)
							WHERE f1.segment='s' AND f1.feature_type_id IN (1, 10, 22) AND aa_stop>10 AND f1.ensp_id NOT IN (select distinct ensp_id from prest_atlas_db.feature where segment='c' )
							GROUP BY et.et_id 
							HAVING SP>1 };
$sth = $lims_dbh->prepare($sql);
$sth->execute();	
$total=0;
$num_et=0;
while (my ($et_id, $uniprot_id, $num_sp) = $sth->fetchrow_array()) {
	$total++;
	$num_et++;
	&insert_et_ids($et_id, $class_id, $uniprot_id, '');
	&insert_et_ids($et_id, $sub_id, $uniprot_id, '');
}
print LOG "Total MDSEC $total\nTotal enst_ids found: $num_et\n";

###############
# Intracellular
###############

$name="Predicted intracellular proteins";
$code="Za";
$public=1;
$description="Predicted intracellular proteins";
$source="HPA";
$external_db="";
$institute="KTH";
$url="http://www.proteinatlas.org/humanproteome/secretome#prediction";
$reference="Linn Fagerberg";
$external_version=localtime;

my $sql= qq{ 	SELECT DISTINCT et_id, uniprot_id
								FROM lims.ensembl_transcripts
								LEFT JOIN (
									SELECT DISTINCT et_id
									FROM lims.proteinclass_transcripts
									JOIN lims.proteinclass USING (class_id)
									JOIN lims.ensembl_transcripts USING (et_id)
									WHERE name IN ('Secreted proteins predicted by MDSEC', 'Membrane proteins predicted by MDM')
								) AS tmp USING (et_id)
								WHERE tmp.et_id IS NULL
							};
	my $sth = $lims_dbh->prepare($sql);
	$sth->execute();
	
	#Analyze one resultrow at a time
	my $total=0;
	while (my ($et_id, $uniprot_id) = $sth->fetchrow_array()) {
		$total++;
		$num_et++;
		&insert_et_ids($et_id, $class_id, $uniprot_id, '');	
	}
	if($notfound==$total) {
		print ERR "$name Error - no transcripts found!\n\n\n";
	}
	print LOG "Total $name: $total\nTotal enst_ids found: $num_et\n";
}

#########################################
# Subroutines used by the other modules
########################################

sub save_TM {
	my $feature_type_id = 11;
	my $segment = 'c';
	my ($db, $ensp_id, $ensg_id, $tm_start, $tm_stop) = @_;
	$sql = qq{ INSERT INTO feature (ensp_id, ensg_id, aa_start, aa_stop, segment, feature_type_id)
						 VALUES (?,?,?,?,?,?) };
	my $sth5 = $db->prepare($sql);
	$sth5->execute($ensp_id, $ensg_id, $tm_start, $tm_stop, $segment, $feature_type_id);
}

sub insert_et_ids {
	my ($et_id, $class_id, $uniprot_id, $external_id) = @_;
	#Insert the class_id
	my $sql= qq{ 	INSERT IGNORE INTO lims.proteinclass_transcripts (class_id, et_id, uniprot_id, external_id) VALUES (?, ?, ?, ?) };
	my $sth = $lims_dbh->prepare($sql);
	$sth->execute($class_id, $et_id, $uniprot_id, $external_id);
}

sub set_parent {
	my ($sub_id, $class_id) = @_;
	my $sql= qq{ 	UPDATE lims.proteinclass SET parent_class_id=? WHERE class_id=? };
	my $sth = $lims_dbh->prepare($sql);
	$sth->execute($class_id, $sub_id);
}

sub set_class_type {
	my ($class_id, $type_id) = @_;
	my $sql= qq{ 	UPDATE lims.proteinclass SET proteinclass_type_id=? WHERE class_id=? };
	my $sth = $lims_dbh->prepare($sql);
	$sth->execute($type_id, $class_id);
}