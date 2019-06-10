<?
function get_sa_protein_function($ensg_id) {
	# use uniprot keywords and some added rules to annotate one function per protein
	$up_prio = ['Blood coagulation','Complement pathway','Acute phase','Cytokine','Hormone','Neuropeptide','Growth factor','Receptor','Transport','Developmental protein','Defense','Enzyme','Enzymeinhibitor','Transcription','Immunity','Cell adhesion'];
	$types['Enzyme'] = ['Acyltransferase','Allosteric enzyme','Aminoacyl-tRNA synthetase','Aminopeptidase','Aminotransferase','Aspartyl protease','Bacteriolytic enzyme','Carboxypeptidase','Decarboxylase','Dioxygenase','Dipeptidase','DNA-directed DNA polymerase','Endonuclease','Excision nuclease','Exonuclease','Glycosidase','Glycosyltransferase','Helicase','Hydrolase','Isomerase','Kinase','Ligase','Lyase','Metalloprotease','Methyltransferase','Monooxygenase','Multifunctional enzyme','Nuclease','Nucleotidyltransferase','Oxidoreductase','Peroxidase','Prenyltransferase','Protease','Protein phosphatase','RNA-directed DNA polymerase','RNA-directed RNA polymerase','Rotamase','Serine esterase','Serine protease','Serine/threonine-protein kinase','Thiol protease','Threonine protease','Topoisomerase','Transferase','Tyrosine-protein kinase'];
	$types['Enzyme inhibitor'] = ["Protease inhibitor","Serine protease inhibitor","Thiol protease inhibitor","Phospholipase A2 inhibitor","Metalloenzyme inhibitor","Metalloprotease inhibitor","Protein phosphatase inhibitor","Aspartic protease inhibitor","Protein kinase inhibitor"];
	$types['Defense'] =  ['Defensin','Antiviral protein','Antimicrobial','Antibiotic','Fungicide',"Antiviral defense"];
	$types['Complement pathway'] = ["Complement pathway","Complement alternate pathway","Complement activation lectin pathway"];
		
	$kws = [];
	$uniprot_ids = get_swissprot_ids_for_ensgid($ensg_id);
	$sql = "SELECT display_name FROM lims.ensembl_genes WHERE ensg_id='$ensg_id'";
	list($gene_name) = sql_query_fetch($sql);	
	
	foreach($uniprot_ids as $uniprot_id) {
		$entry = \HPA\Uniprot\Uniprot::get_decoded_parsed_entry_from_db($uniprot_id);	
		$up_kws = get_uniprot_keywords($entry);
		$kws = array_merge($kws, @$up_kws['Biological process']?:[],@$up_kws['Molecular function']?:[]);
	}
	$match = array_intersect($up_prio, $kws);
	$function = array_shift($match);
	if (empty($function)) {		
		foreach($types as $type=>$type_kws) {			
			if ($match = array_intersect($kws, $type_kws)) {
				$function = $type;
				break;
			}
		}
	}
	if ($function == 'Cytokine') {
		if (substr($gene_name,0,2) == "IL") {
			$function = 'Interleukin';			
		}
		if (substr($gene_name,0,3) == "IFN") {
			$function = 'Interferon';			
		}
		if (in_array(substr($gene_name,0,3),["CCL","CXC","CX3","XCL"])) {
			$function = 'Chemokine';
		}
	}
	else if (substr($gene_name,0,3) == "APO") {
			$function = 'Apolipoprotein';
	}
	else if (empty($function) AND in_array('Angiogenesis',$kws)) {
		$function = 'Developmental protein';
	}
	else if (substr($gene_name,0,3) == "C1Q") {
		$function = 'Complement pathway';
	}
	else if (substr($gene_name,0,4) == "CFHR") {
		$function = 'Complement pathway';
	}
	else if (in_array($gene_name,["ALB","SHBG","SERPINA7"])) {
		$function = 'Transport';
	}
	else if (substr($gene_name,0,4) == "SEMA") {
		$function = 'Developmental protein';
	}	
	
	if ($function) {
		return $function;
	} else {
		if (!empty($kws)) {
			return 'Other';
		} else {
			return 'No annotated function';
		}
	}
}

function get_swissprot_ids_for_ensgid($ensg_id) {
	$sql = "SELECT DISTINCT external_id
					FROM ensembl_genes
					JOIN ensembl_transcripts USING (eg_id)
					JOIN ensembl_transcripts_xref USING (et_id)
					JOIN ensembl_xref x USING (xref_id)
					JOIN ensembl_xref_db xd ON xd.xref_db_id=x.xref_db_id AND xd.db_name='Uniprot/SWISSPROT'
					WHERE ensg_id='$ensg_id'
					ORDER BY ensg_id, external_id";
	$res = sql_query($sql);
	$out = array();
	while (list($external_id) = sql_fetch($res)) {
		$out[] = $external_id;
	}
	return $out;
}

function get_uniprot_keywords($entry) {
	$keywords = [];
	foreach ($entry->keyword as $keyword) {
		$sql = "SELECT category FROM uniprot_keyword WHERE keyword='$keyword'";
		list($category) = sql_query_fetch($sql);
		$keywords[$category][] = $keyword;
	}
	return $keywords;
}