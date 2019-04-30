# OrthoDB

## API Commands

### /tree

> Retreive the full tree of species

[Json Tree](https://www.orthodb.org/tree)

### /search

> Find all cluster id's matching a given query.

#### Arguments

- query : full query string
- ncbi : flag - if 0, then generic search, if 1 the query is assumed to be a NCBI gene id
- level : NCBI taxon id of the clade
- skip : number of hits to skip
- limit : maximum nr of hits (cluster ids) to return - default is 1000
- universal : phyloprofile filter, present in 1.0, 0.9, 0.8 of all species in the clade
- singlecopy :  phyloprofile filter, singlecopy in 1.0, 0.9, 0.8 of all species in the clade

#### Returns

- A list of clusters, the maximum number of clusters is defined by 'limit'

#### Example 1

```sh
https://www.orthodb.org/search?query=p450&limit=2&level=33208&singlecopy=0.8
```

```json
{"status": "ok", "message": null, "data": ["101276at33208", "104147at33208", "102241at33208"], "count": 552, "skip": 0, "limit": 3, "query": "p450", "level": 33208, "url": "https://www.orthodb.org//search?query=p450&limit=3&level=33208", "universal": null, "singlecopy": null, "inclusive": 1}
```

---

### /blast

> This finds all cluster id's with genes matching the given sequence. The list is sorted with the best matching cluster first.

#### Argument

- level : NCBI taxon id of the clade.
- skip : number of hits to skip.
- limit : maximum nr of hits (cluster ids) to return - default is 1000.
- universal : phyloprofile filter, present in 1.0, 0.9, 0.8 of all species in the clade.
- singlecopy :  phyloprofile filter, singlecopy in 1.0, 0.9, 0.8 of all species in the clade.
- seq : sequence string, without fasta-header.
- species : comma separated list of NCBI numerical taxonomy ids.
- inclusive : flag; 0 - return clusters containing at least one of given species, 1 - return all matches ignoring the species list (default).

#### Returns

- list of OrthoDB clusters

#### Example 2

```sh
https://www.orthodb.org/blast?level=33208&seq=MGDSHEDTSATVPEAVAEEVSLFSTTDIVLF
```

```json
{"status": "ok", "message": null, "data": ["405462at33208"], "count": 1, "skip": 0, "limit": 10000, "query": null, "level": 33208, "url": "https://www.orthodb.org//blast?level=33208&seq=MGDSHEDTSATVPEAVAEEVSLFSTTDIVLF", "seq": "MGDSHEDTSATVPEAVAEEVSLFSTTDIVLF", "species": null, "inclusive": 1, "universal": null, "singlecopy": null
```

---

### /group

> Retrieve detailed annotation information on the given cluster.

#### Argument

- id : OrthoDB cluster id.

#### Returns

> Annotation details on the given cluster id

#### Example 3

```sh
https://www.orthodb.org/group?id=716at7742
```

```json
{"status": "ok", "message": null, "data": {"cellular_component": [{"count": 73, "description": "dynein complex", "id": "GO:0030286", "name": "GO:0030286", "type": "GeneOntology"}], "interpro_domains": [{"count": 47, "description": "Dynein heavy chain", "id": "IPR026983", "name": "IPR026983", "type": "interpro"}, {"count": 45, "description": "P-loop containing nucleoside triphosphate hydrolase", "id": "IPR027417", "name": "IPR027417", "type": "interpro"}, {"count": 44, "description": "Dynein heavy chain, hydrolytic ATP-binding dynein motor region D1", "id": "IPR035699", "name": "IPR035699", "type": "interpro"}, {"count": 44, "description": "Dynein heavy chain, domain-2", "id": "IPR013602", "name": "IPR013602", "type": "interpro"}, {"count": 44, "description": "Dynein heavy chain, AAA module D4", "id": "IPR024317", "name": "IPR024317", "type": "interpro"}, {"count": 44, "description": "ATPase, dynein-related, AAA domain", "id": "IPR011704", "name": "IPR011704", "type": "interpro"}, {"count": 43, "description": "Dynein heavy chain, coiled coil stalk", "id": "IPR024743", "name": "IPR024743", "type": "interpro"}, {"count": 43, "description": "Dynein heavy chain, ATP-binding dynein motor region D5", "id": "IPR035706", "name": "IPR035706", "type": "interpro"}, {"count": 43, "description": "Dynein heavy chain domain", "id": "IPR004273", "name": "IPR004273", "type": "interpro"}, {"count": 42, "description": "Dynein heavy chain, domain-1", "id": "IPR013594", "name": "IPR013594", "type": "interpro"}, {"count": 41, "description": "AAA+ ATPase domain", "id": "IPR003593", "name": "IPR003593", "type": "interpro"}], "id": "716at7742", "public_id": "716at7742", "tax_id": 7742, "level_name": "Vertebrata", "name": "dynein heavy chain 8, axonemal", "gene_architecture": {"protein_median_length": 4604, "protein_stdev_length": 1106.9, "exon_median_counts": 87, "exon_stdev_counts": 12.77}, "biological_process": [{"count": 79, "description": "microtubule-based movement", "id": "GO:0007018", "name": "GO:0007018", "type": "GeneOntology"}], "molecular_function": [{"count": 87, "description": "microtubule motor activity", "id": "GO:0003777", "name": "GO:0003777", "type": "GeneOntology"}, {"count": 75, "description": "ATP binding", "id": "GO:0005524", "name": "GO:0005524", "type": "GeneOntology"}], "phyletic_profile": {"gene_count": 609, "species_count": 243, "present_in": 243, "multi_copy": 229, "single_copy": 14}, "evolutionary_rate": 0.9204, "functional_category": [{"description": "U: Intracellular trafficking, secretion, and vesicular transport", "id": "U", "type": "text"}, {"description": "L: Replication, recombination and repair", "id": "L", "type": "text"}, {"description": "V: Defense mechanisms", "id": "V", "type": "text"}]}, "url": "https://www.orthodb.org//group?id=716at7742"}
```

---

### /orthologs

> Retrieve all genes in a given cluster, possibly filtered wrt species.

#### Arguments

- id : OrthoDB cluster id.
- species : optional comma-separated list of species taxid's.

#### Returns

- A dictionary of tax id's, each contain a list of OrthoDB gene id's

#### Example 4

```sh
https://www.orthodb.org/orthologs?id=716at7742&species=9606_0
```

```json
{"status": "ok", "message": null, "data": [{"organism": {"description": "man", "param": "species in Primates", "id": "9606_0", "name": "Homo sapiens", "type": "taxonomy"}, "species_statistics": {"in_clusters_count": 156629, "genes_count": 21416, "mapping_type": "C"}, "genes": [{"gene_id": {"id": "DNAH5", "param": "9606_0:0017fc"}, "exons": "86", "interpro": [{"description": "Dynein heavy chain", "param": "107..4624", "id": "IPR026983", "type": "interpro"}, {"description": "Dynein heavy chain, domain-1", "param": "252..802", "id": "IPR013594", "type": "interpro"}, {"description": "Dynein heavy chain, domain-2", "param": "1404..1807", "id": "IPR013602", "type": "interpro"}, {"description": "Dynein heavy chain, hydrolytic ATP-binding dynein motor region D1", "param": "1942..2173", "id": "IPR035699", "type": "interpro"}, {"description": "P-loop containing nucleoside triphosphate hydrolase", "param": "1959..2174", "id": "IPR027417", "type": "interpro"}, {"description": "AAA+ ATPase domain", "param": "1975..2111", "id": "IPR003593", "type": "interpro"}, {"description": "ATPase, dynein-related, AAA domain", "param": "2260..2392", "id": "IPR011704", "type": "interpro"}, {"description": "Dynein heavy chain, AAA module D4", "param": "2924..3186", "id": "IPR024317", "type": "interpro"}, {"description": "Dynein heavy chain, coiled coil stalk", "param": "3202..3544", "id": "IPR024743", "type": "interpro"}, {"description": "Dynein heavy chain, ATP-binding dynein motor region D5", "param": "3575..3795", "id": "IPR035706", "type": "interpro"}, {"description": "Dynein heavy chain domain", "param": "3935..4621", "id": "IPR004273", "type": "interpro"}], "aas": "4660", "how_much_more_info": 3, "description": "Dynein heavy chain 5, axonemal", "uniprot": {"id": "Q8TE73", "name": "Dynein heavy chain 5, axonemal", "type": "uniprot"}, "more_info": true}, {"gene_id": {"id": "DNAH8", "param": "9606_0:0019b4"}, "exons": "101", "interpro": [{"description": "Dynein heavy chain", "param": "22..4490", "id": "IPR026983", "type": "interpro"}, {"description": "Dynein heavy chain, domain-1", "param": "139..693", "id": "IPR013594", "type": "interpro"}, {"description": "Dynein heavy chain, domain-2", "param": "1268..1673", "id": "IPR013602", "type": "interpro"}, {"description": "Dynein heavy chain, hydrolytic ATP-binding dynein motor region D1", "param": "1808..2039", "id": "IPR035699", "type": "interpro"}, {"description": "P-loop containing nucleoside triphosphate hydrolase", "param": "1825..2034", "id": "IPR027417", "type": "interpro"}, {"description": "AAA+ ATPase domain", "param": "1841..1985", "id": "IPR003593", "type": "interpro"}, {"description": "ATPase, dynein-related, AAA domain", "param": "2125..2259", "id": "IPR011704", "type": "interpro"}, {"description": "Dynein heavy chain, AAA module D4", "param": "2788..3052", "id": "IPR024317", "type": "interpro"}, {"description": "Dynein heavy chain, coiled coil stalk", "param": "3065..3412", "id": "IPR024743", "type": "interpro"}, {"description": "Dynein heavy chain, ATP-binding dynein motor region D5", "param": "3440..3660", "id": "IPR035706", "type": "interpro"}, {"description": "Dynein heavy chain domain", "param": "3801..4487", "id": "IPR004273", "type": "interpro"}], "aas": "4707", "how_much_more_info": 3, "description": "Dynein heavy chain 8, axonemal", "uniprot": {"id": "A0A075B6F3", "name": "Dynein heavy chain 8, axonemal", "type": "uniprot"}, "more_info": true}]}], "url": "https://www.orthodb.org//orthologs?id=716at7742&species=9606_0"}
```

---

### /ogdetails

> Retrieve further details on a given gene id.

#### Argument

- id : detailed information on the given gene id.

#### Returns

- detailed information on the given gene id

#### Example 5

```sh
https://www.orthodb.org/ogdetails?id=9606_0:0017fc
```

---

### /siblings

> Retrieve all siblings to the given cluster.

#### Argument

- id : detailed information on the given gene id.
- limit : max nr of returned siblings

#### Returns

- a list of OrthoDB cluster id's

#### Example 6

```sh
https://www.orthodb.org/siblings?id=716at7742&limit=2
```

```json
{"status": "ok", "message": null, "data": [{"param": "68%", "interpro": [{"description": "Dynein heavy chain", "id": "IPR026983", "name": "IPR026983", "type": "interpro"}, {"description": "P-loop containing nucleoside triphosphate hydrolase", "id": "IPR027417", "name": "IPR027417", "type": "interpro"}, {"description": "Dynein heavy chain, domain-2", "id": "IPR013602", "name": "IPR013602", "type": "interpro"}, {"description": "ATPase, dynein-related, AAA domain", "id": "IPR011704", "name": "IPR011704", "type": "interpro"}, {"description": "Dynein heavy chain, AAA module D4", "id": "IPR024317", "name": "IPR024317", "type": "interpro"}, {"description": "Dynein heavy chain, domain-1", "id": "IPR013594", "name": "IPR013594", "type": "interpro"}, {"description": "Dynein heavy chain, hydrolytic ATP-binding dynein motor region D1", "id": "IPR035699", "name": "IPR035699", "type": "interpro"}, {"description": "Dynein heavy chain, coiled coil stalk", "id": "IPR024743", "name": "IPR024743", "type": "interpro"}, {"description": "Dynein heavy chain, ATP-binding dynein motor region D5", "id": "IPR035706", "name": "IPR035706", "type": "interpro"}, {"description": "Dynein heavy chain domain", "id": "IPR004273", "name": "IPR004273", "type": "interpro"}, {"description": "AAA+ ATPase domain", "id": "IPR003593", "name": "IPR003593", "type": "interpro"}], "description": "68% overlap between 716at7742 and 1570at7742 by 11 Interpro ID(s)", "type": "ODBcluster", "name": "1570at7742", "id": "1570at7742"}, {"param": "57%", "interpro": [{"description": "Dynein heavy chain", "id": "IPR026983", "name": "IPR026983", "type": "interpro"}, {"description": "P-loop containing nucleoside triphosphate hydrolase", "id": "IPR027417", "name": "IPR027417", "type": "interpro"}, {"description": "Dynein heavy chain, domain-2", "id": "IPR013602", "name": "IPR013602", "type": "interpro"}, {"description": "Dynein heavy chain, AAA module D4", "id": "IPR024317", "name": "IPR024317", "type": "interpro"}, {"description": "ATPase, dynein-related, AAA domain", "id": "IPR011704", "name": "IPR011704", "type": "interpro"}, {"description": "Dynein heavy chain domain", "id": "IPR004273", "name": "IPR004273", "type": "interpro"}, {"description": "Dynein heavy chain, hydrolytic ATP-binding dynein motor region D1", "id": "IPR035699", "name": "IPR035699", "type": "interpro"}, {"description": "Dynein heavy chain, domain-1", "id": "IPR013594", "name": "IPR013594", "type": "interpro"}, {"description": "Dynein heavy chain, ATP-binding dynein motor region D5", "id": "IPR035706", "name": "IPR035706", "type": "interpro"}, {"description": "Dynein heavy chain, coiled coil stalk", "id": "IPR024743", "name": "IPR024743", "type": "interpro"}, {"description": "AAA+ ATPase domain", "id": "IPR003593", "name": "IPR003593", "type": "interpro"}], "description": "57% overlap between 716at7742 and 5176at7742 by 11 Interpro ID(s)", "type": "ODBcluster", "name": "5176at7742", "id": "5176at7742"}], "url": "https://www.orthodb.org//siblings?id=716at7742&limit=2"}
```

---

### /fasta

> Retrieve FASTA sequences of the given clusters.


#### Argument 1

- id : OrthoDB cluster id.
- species : list of NCBI species taxonomy id's.

#### Argument 2

- all arguments for /search.
- species : list of NCBI species taxonomy id's.

#### Returns

- sequences in fasta format.

#### Example 7

```sh
https://www.orthodb.org/fasta?id=716at7742
```

---

### /tab

#### Argument

- same arguments as for /fasta
- long : 0 (default) -> without sequence ; 1 -> include sequence.

#### Returns

- tab-separated table of gene annotations

#### Example 7

```sh
https://www.orthodb.org/tab?id=716at7742&species=9606_0&long=0
```

```json
pub_og_id	og_name	level_taxid	organism_taxid	organism_name	int_prot_id	pub_gene_id	description
716at7742	dynein heavy chain 8, axonemal	7742	9606_0	Homo sapiens	9606_0:0017fc	DNAH5	Dynein heavy chain 5, axonemal
716at7742	dynein heavy chain 8, axonemal	7742	9606_0	Homo sapiens	9606_0:0019b4	DNAH8	Dynein heavy chain 8, axonemal
```

