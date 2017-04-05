import pymongo
import utils
import exac


def load_analyse(cohort_name,analyse_name):
    db = exac.get_db()


def filterForFrequency(variant, key):
    af = float(variant[key]) if key in variant and variant[key] else 0
    if af > 0.01:
        return False
    else:
        return True
        
def filterForAllFrequencies(variant):
    if filterForFrequency(variant,'ExAC_AF') and filterForFrequency(variant,'g1k') and filterForFrequency(variant,'Kaviar_AF') and filterForFrequency(variant,'HRC_AF') and filterForFrequency(variant,'ESP') :
        return True
    else:
        return False
        
def filterForSpliceSite(variant):
    key = 'dbscSNV_ADA_SCORE'
    ada = float(variant[key]) if key in variant and variant[key] else 0
    if ada > 0.6:
        return True
    else:
        return False

def filterByFilter(cohort):
    if cohort['filter'] == 'PASS':
        return True
    else:
        return False
        
def worst_csq_from_csq(csq):
    """
    Input possibly &-filled csq string (e.g. 'non_coding_exon_variant&nc_transcript_variant')
    Return the worst consequence (In this case, 'non_coding_exon_variant')
    """
    return rev_csq_order_dict[worst_csq_index(csq.split('&'))]
        

def worst_csq_index(csq_list):
    """
    Input list of consequences (e.g. ['frameshift_variant', 'missense_variant'])
    Return index of the worst consequence (In this case, index of 'frameshift_variant', so 4)
    Works well with worst_csq_index('non_coding_exon_variant&nc_transcript_variant'.split('&'))
    """
    return min([csq_order_dict[csq] for csq in csq_list])
    
def filterByCsq(csq,variant):
    conseq = list()
    for ann in variant['vep_annotations']:
        conseq.append(worst_csq_from_csq(ann['Consequence']))
    conseq =  sorted(conseq, key=(lambda ann:csq_order_dict[ann]))
    if len(conseq) > 0 and csq_order_dict[conseq[0]] <= csq_order_dict[csq]:
        variant['major_consequence'] = conseq[0]
        return True
    else:
        return False
        
def filterByGenes(variant):
    intersection = set.intersection(set(variant['genes']), set(ibd_genes))
    disease_genes = list()
    for ensembl in list(intersection):
        uniprot = ensembl2uniprot[ensembl]
        disease_gene = dict()
        disease_gene['ensembl'] = ensembl
        disease_gene['uniprot'] = uniprot
        disease_genes.append(disease_gene)

    variant['disease_genes'] = disease_genes
    return bool(intersection)
    
        
def filterByHomCount(variant,cohortName):
    for cohort in variant['cohorts']:
        if cohort['name'] == cohortName:
            return float(cohort['hom_count']) >= 1
            break
    return False
    
def prettyPrint(element):
    if isinstance(element, list):
        print element
        return ",".join(element)
    else:
        return element
    
    
# Note that this is the current as of v81 with some included for backwards compatibility (VEP <= 75)
csq_order = ["transcript_ablation",
"splice_acceptor_variant",
"splice_donor_variant",
"stop_gained",
"frameshift_variant",
"stop_lost",
"start_lost",  # new in v81
"initiator_codon_variant",  # deprecated
"transcript_amplification",
"inframe_insertion",
"inframe_deletion",
"missense_variant",
"protein_altering_variant",  # new in v79
"splice_region_variant",
"incomplete_terminal_codon_variant",
"stop_retained_variant",
"synonymous_variant",
"coding_sequence_variant",
"mature_miRNA_variant",
"5_prime_UTR_variant",
"3_prime_UTR_variant",
"non_coding_transcript_exon_variant",
"non_coding_exon_variant",  # deprecated
"intron_variant",
"NMD_transcript_variant",
"non_coding_transcript_variant",
"nc_transcript_variant",  # deprecated
"upstream_gene_variant",
"downstream_gene_variant",
"TFBS_ablation",
"TFBS_amplification",
"TF_binding_site_variant",
"regulatory_region_ablation",
"regulatory_region_amplification",
"feature_elongation",
"regulatory_region_variant",
"feature_truncation",
"intergenic_variant",
"other_variant",
"",
"?"]
csq_order_dict = {csq:i for i,csq in enumerate(csq_order)}
rev_csq_order_dict = dict(enumerate(csq_order))

#init genes
gene_mapping_file = open('gene_mapping.txt','r')
uniprot2ensembl = dict()
ensembl2uniprot = dict()
for line in gene_mapping_file:
    ensembl = line.strip().split(",")[0]
    uniprot = line.strip().split(",")[1]
    uniprot2ensembl[uniprot] = ensembl
    ensembl2uniprot[ensembl] = uniprot
ibd_gene_file = open('ibd_genes.txt','r')
ibd_genes = list()
for line in ibd_gene_file:
    gene = line.strip()
    ibd_genes.append(uniprot2ensembl[gene])


client = pymongo.MongoClient(host='localhost', port=27017)
db = client['exac']
db.analyses.drop()
#db.analyses.drop_indexes()
#db.variants.ensure_index('cohorts.name')
#db.variants.ensure_index('ExAC_AF')
#db.variants.ensure_index('g1k')
#db.variants.ensure_index('Kaviar_AF')
#db.variants.ensure_index('HRC_AF')
#db.variants.ensure_index('ESP')
#db.variants.ensure_index('dbscSNV_ADA_SCORE')
print db.variants.index_information()
#import sys
#sys.exit()
variants = db.variants
query = dict()
cohort_name = 'pediatricIBD'
exac_or = {'$or' : [{'ExAC_AF': {'$exists': False}}, {'ExAC_AF': {'$lt': 0.01}},{'ExAC_AF': {'$eq': ''}}]}
g1k_or = {'$or' : [{'g1k': {'$exists': False}}, {'g1k': {'$lt': 0.01}},{'g1k': {'$eq': ''}}]}
kaviar = {'$or' : [{'Kaviar_AF': {'$exists': False}}, {'Kaviar_AF': {'$lt': 0.01}},{'Kaviar_AF': {'$eq': ''}}]}
hrc = {'$or' : [{'HRC_AF': {'$exists': False}}, {'HRC_AF': {'$lt': 0.01}},{'HRC_AF': {'$eq': ''}}]}
esp = {'$or' : [{'ESP': {'$exists': False}}, {'ESP': {'$lt': 0.01}},{'ESP': {'$eq': ''}}]}
cohort = {'cohorts.name': {'$eq': cohort_name}}
query['$and'] = [exac_or,g1k_or,kaviar,hrc,esp,cohort]
variant_all = variants.find(query)

result_variants = list()
nr_results = 0
nr_all = 0
nr_bygene = 0
nr_bycsq = 0
print "Start filtering"
for variant in variant_all:
    nr_all += 1
    if filterByCsq('missense_variant',variant) or filterForSpliceSite(variant):
        nr_bycsq += 1
        if filterByGenes(variant):
            for cohort in variant['cohorts']:
                if cohort['name'] == cohort_name:
                    if filterByFilter(cohort):
                        nr_bygene += 1
                        #if filterByHomCount(variant,'IBM_ped'):
                        variant['analyse_name'] = cohort_name
                        result_variants.append(variant)
                        nr_results += 1
                            #if nr_results % 100 == 0:
                            #print nr_results
                        break;

print nr_all
print nr_bycsq
print nr_bygene
print nr_results

for variant in result_variants:
    utils.add_consequence_to_variant(variant)
    variant.pop('vep_annotations', None)
    for cohort in variant['cohorts']:
        if cohort['name'] == cohort_name:
            variant['hom_count'] = cohort['hom_count']
            variant['allele_count'] = cohort['allele_count']
            variant['allele_freq'] = cohort['allele_freq']
            variant['allele_num'] = cohort['allele_num']
            samples = list()
            for track in cohort['tracks']:
                for sample in track['samples']:
                    samples.append(sample)
            variant['samples'] = samples
            break
    variant.pop('cohorts', None)

    db.analyses.insert(variant)

#header = ['chrom','pos','ref','alt','genes','disease_genes','allele_num','hom_count']
#with open( 'filtered_variants.txt', 'wb' ) as out_file:
#    csv_w = csv.writer( out_file ,delimiter=";" )
#    csv_w.writerow( header)
#
#    for i_r in result_variants:
#        print i_r
#        csv_w.writerow( map( lambda x: prettyPrint(i_r.get( x, "" )), header ))