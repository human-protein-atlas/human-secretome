<?
# Classify gene
function classify_gene($eg_id) {
	$cl['comment'] = '';
	$cl['location'] = '';
	$cl['automatic_source'] = '';
	$up_ids = get_up_ids($eg_id);
	$up = get_up_data($up_ids);
	$eg = get_eg_data($eg_id);
	$hpa_secreted = has_secreted_transcript($eg_id);
	$up_secreted = in_array('Secreted',$up['keywords']);
	# Is hpa secreted or uniprot secreted	
	if ($hpa_secreted) {
		if ($up_secreted) {
			$cl['automatic_source'] = 'HPA_UP_SEC';
		} else {
			$cl['automatic_source'] = 'HPA_SEC';
		} 
	} else if ($up_secreted) {
		$cl['automatic_source'] = 'UP_SEC';
	} else {
		# not secreted
		$cl['automatic_source'] = 'NOT_SECRETED';
		$cl['comment'] = 'Not secreted in uniprot or ensembl';
		$cl['location'] = 'NOT_SECRETED';
	}
	return $cl;
}
	
function has_secreted_transcript($eg_id) {
	$sql = "SELECT 1
					FROM ensembl_transcripts
					JOIN proteinclass_transcripts USING (et_id)
					JOIN proteinclass USING (class_id)
					WHERE eg_id=$eg_id AND code='Sa' AND tsl=1";
	$res = sql_query($sql);
	return sql_num_rows($res);		
}

# get all uniprot ids from xref
function get_up_ids($eg_id) {
	$up_ids = [];
	# Get all uniprot ids for gene
	$sql = "SELECT DISTINCT external_id						
					FROM ensembl_transcripts
					JOIN ensembl_transcripts_xref USING (et_id)
					JOIN ensembl_xref x USING (xref_id)
					JOIN ensembl_xref_db xd USING (xref_db_id)
					WHERE eg_id=$eg_id AND xd.db_name='Uniprot/SWISSPROT'
					ORDER BY external_id";
	$res = sql_query($sql);
	while	(list($uniprot_id) = sql_fetch($res)) {
		$up_ids[] = $uniprot_id; 
	}
	return $up_ids;
}

# uniprot data
function get_up_data($up_ids) {	
	$up['locations'] = [];
	$up['keywords'] = [];	
	if (empty($up_ids)) return $up;
	foreach($up_ids as $up_id) {
		try {
			$entry_version = \HPA\Uniprot\Uniprot::get_latest_entry_version_from_db($up_id);
			$entry = \HPA\Uniprot\Uniprot::get_decoded_parsed_entry_from_db($up_id, $entry_version);		
		} catch (Exception $e) {
			echo "Uniprot entry not found: $up_id\r\n";
			continue;
		}
		# up subcellular location
		if (isset($entry->subcellular_location) && !empty((array)$entry->subcellular_location)) {
			$locs = [];
			$molecules = [];
			foreach ($entry->subcellular_location as $molecule => $molvalues) {						
				foreach ($molvalues->locations as $location_group) {
					foreach ($location_group->locations as $locName) {
						$up['locations'][] = $locName;
					}
				}
			}
		}
		# up keywords
		if (isset($entry->keyword) && !empty((array)$entry->keyword)) {
			$up['keywords'] = $entry->keyword;
		}
	}
	return $up;
}
function get_eg_data($eg_id) {
	$sql = "SELECT display_name,description,ensg_id FROM ensembl_genes WHERE eg_id=$eg_id";
	list($eg['gene_name'],$eg['description'],$eg['ensg_id']) = sql_query_fetch($sql);	
	return $eg;	
}