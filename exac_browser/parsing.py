"""
Utils for reading flat files that are loaded into database
"""
import re
import traceback
from utils import *
import time
from decimal import *
import copy

POPS = {
    'AFR': 'African',
    'AMR': 'Latino',
    'EAS': 'East Asian',
    'FIN': 'European (Finnish)',
    'NFE': 'European (Non-Finnish)',
    'SAS': 'South Asian',
    'OTH': 'Other'
}


def get_base_coverage_from_file(base_coverage_file):
    """
    Read a base coverage file and return iter of dicts that look like:
    {
        'xpos': 1e9+1,
        'mean': 0.0,
        'median': 0.0,
        '1': 0.0,
        '5': 0.0,
        '10': 0.0,
        '15': 0.0,
        '20': 0.0,
        '25': 0.0,
        '30': 0.0,
        '50': 0.0,
        '100': 0.0,
    }
    """

    float_header_fields = ['mean', 'median', '1', '5', '10', '15', '20', '25', '30', '50', '100']
    for line in base_coverage_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')
        d = {
            'xpos': get_xpos(fields[0], int(fields[1])),
            'pos': int(fields[1]),
        }
        for i, k in enumerate(float_header_fields):
            d[k] = float(fields[i+2])
        yield d


def get_variants_from_sites_vcf(sites_vcf):
    """
    Parse exac sites VCF file and return iter of variant dicts
    sites_vcf is a file (gzipped), not file path
    """
    vep_field_names = None
    for line in sites_vcf:
        try:
            line = line.strip('\n')
            if line.startswith('##INFO=<ID=CSQ'):
                vep_field_names = line.split('Format: ')[-1].strip('">').split('|')
            if line.startswith('#'):
                continue

            # If we get here, it's a variant line
            if vep_field_names is None:
                raise Exception("VEP_field_names is None. Make sure VCF header is present.")
            # This elegant parsing code below is copied from https://github.com/konradjk/loftee
            fields = line.split('\t')
            info_field = dict([(x.split('=', 1)) if '=' in x else (x, x) for x in re.split(';(?=\w)', fields[7])])
            consequence_array = info_field['CSQ'].split(',') if 'CSQ' in info_field else []
            annotations = [dict(zip(vep_field_names, x.split('|'))) for x in consequence_array if len(vep_field_names) == len(x.split('|'))]
            coding_annotations = [ann for ann in annotations if ann['Feature'].startswith('ENST')]

            alt_alleles = fields[4].split(',')

            # different variant for each alt allele
            for i, alt_allele in enumerate(alt_alleles):

                vep_annotations = [ann for ann in coding_annotations if int(ann['ALLELE_NUM']) == i + 1]

                # Variant is just a dict
                # Make a copy of the info_field dict - so all the original data remains
                # Add some new keys that are allele-specific
                pos, ref, alt = get_minimal_representation(fields[1], fields[3], alt_allele)

                exacvariant = {}
                exacvariant['chrom'] = fields[0]
                exacvariant['pos'] = pos
                exacvariant['rsid'] = fields[2]
                exacvariant['xpos'] = get_xpos(exacvariant['chrom'], exacvariant['pos'])
                exacvariant['ref'] = ref
                exacvariant['alt'] = alt
                exacvariant['xstart'] = exacvariant['xpos']
                exacvariant['xstop'] = exacvariant['xpos'] + len(exacvariant['alt']) - len(exacvariant['ref'])
                exacvariant['variant_id'] = '{}-{}-{}-{}'.format(exacvariant['chrom'], exacvariant['pos'], exacvariant['ref'], exacvariant['alt'])
                exacvariant['orig_alt_alleles'] = [
                    '{}-{}-{}-{}'.format(exacvariant['chrom'], *get_minimal_representation(fields[1], fields[3], x))
                    for x in alt_alleles
                ]

                exacvariant['site_quality'] = float(fields[5])
                exacvariant['filter'] = 'PASS' if fields[6] == 'PASS' or ('filter' in exacvariant and exacvariant['filter'] == 'PASS') else fields[6]
                exacvariant['vep_annotations'] = vep_annotations

                old_allele_count = exacvariant['allele_count'] if 'allele_count' in exacvariant else 0
                exacvariant['allele_count'] = int(info_field['AC'].split(',')[i]) + old_allele_count #changed from AC_Adj to AC
                if not exacvariant['allele_count'] and exacvariant['filter'] == 'PASS': exacvariant['filter'] = 'AC_Adj0' # Temporary filter
                old_allele_num = exacvariant['allele_num'] if 'allele_num' in exacvariant else 0
                exacvariant['allele_num'] = int(info_field['AN']) + old_allele_num #changed from AN_adj to AN

                if exacvariant['allele_num'] > 0:
                    exacvariant['allele_freq'] = exacvariant['allele_count']/float(exacvariant['allele_num']) #changed from AC_Adj to AN
                else:
                    exacvariant['allele_freq'] = None

                exacvariant['pop_acs'] = dict([(POPS[x], int(info_field['AC_%s' % x].split(',')[i])) for x in POPS])
                exacvariant['pop_ans'] = dict([(POPS[x], int(info_field['AN_%s' % x])) for x in POPS])
                exacvariant['pop_homs'] = dict([(POPS[x], int(info_field['Hom_%s' % x].split(',')[i])) for x in POPS])
                exacvariant['an_male'] = info_field['AN_MALE']
                exacvariant['an_female'] = info_field['AN_FEMALE']
                exacvariant['hom_count'] = sum(exacvariant['pop_homs'].values())
                if exacvariant['chrom'] in ('X', 'Y'):
                    exacvariant['pop_hemis'] = dict([(POPS[x], int(info_field['Hemi_%s' % x].split(',')[i])) for x in POPS])
                    exacvariant['hemi_count'] = sum(exacvariant['pop_hemis'].values())
                exacvariant['quality_metrics'] = dict([(x, info_field[x]) for x in METRICS if x in info_field])

                exacvariant['genes'] = list({annotation['Gene'] for annotation in vep_annotations})
                exacvariant['transcripts'] = list({annotation['Feature'] for annotation in vep_annotations})

                yield exacvariant
        except Exception:
            print("Error parsing vcf line: " + line)
            traceback.print_exc()
            break
            
def get_variants_from_sites_vcf_ikmb(sites_vcf,cohort_name):
    """
    Parse exac sites VCF file and return iter of variant dicts
    sites_vcf is a file (gzipped), not file path
    """
    vep_field_names = None
    sample_names = None
    getcontext().prec = 5


    for line in sites_vcf:
        try:
            line = line.strip('\n')
            if line.startswith('##INFO=<ID=CSQ'):
                vep_field_names = line.split('Format: ')[-1].strip('">').split('|')
            if line.startswith('##INFO=<ID=DP_HIST'):
                dp_mids = map(float, line.split('Mids: ')[-1].strip('">').split('|'))
            if line.startswith('##INFO=<ID=GQ_HIST'):
                gq_mids = map(float, line.split('Mids: ')[-1].strip('">').split('|'))
            if line.startswith('#CHROM'):
                sample_names = line.split('\t')[9::]
            if line.startswith('#'):
                continue

            # If we get here, it's a variant line
            if vep_field_names is None:
                raise Exception("VEP_field_names is None. Make sure VCF header is present.")
            # This elegant parsing code below is copied from https://github.com/konradjk/loftee
            fields = line.split('\t')
            info_field = dict([(x.split('=', 1)) if '=' in x else (x, x) for x in re.split(';(?=\w)', fields[7])])
            consequence_array = info_field['CSQ'].split(',') if 'CSQ' in info_field else []
            annotations = [dict(zip(vep_field_names, x.split('|'))) for x in consequence_array if len(vep_field_names) == len(x.split('|'))]
            coding_annotations = [ann for ann in annotations if ann['Feature'].startswith('ENST')]

            alt_alleles = fields[4].split(',')

            # different variant for each alt allele
            for i, alt_allele in enumerate(alt_alleles):

                vep_annotations = [ann for ann in coding_annotations if int(ann['ALLELE_NUM']) == i + 1]

                # Variant is just a dict
                # Make a copy of the info_field dict - so all the original data remains
                # Add some new keys that are allele-specific
                pos, ref, alt = get_minimal_representation(fields[1], fields[3], alt_allele)
                chrom = fields[0].replace("chr","")
                variant = get_variant_and_replace(chrom,pos,ref,alt)
                if variant == None:
                    variant = {}
                #variant = {}

                variant['chrom'] = chrom
                variant['pos'] = pos
                variant['rsid'] = fields[2]
                variant['xpos'] = get_xpos(variant['chrom'], variant['pos'])
                variant['ref'] = ref
                variant['alt'] = alt
                variant['xstart'] = variant['xpos']
                variant['xstop'] = variant['xpos'] + len(variant['alt']) - len(variant['ref'])
                variant['variant_id'] = '{}-{}-{}-{}'.format(variant['chrom'], variant['pos'], variant['ref'], variant['alt'])
                variant['orig_alt_alleles'] = [
                    '{}-{}-{}-{}'.format(variant['chrom'], *get_minimal_representation(fields[1], fields[3], x))
                    for x in alt_alleles
                ]
                variant['site_quality'] = float(fields[5])
                variant['filter'] = 'PASS' if fields[6] == 'PASS' or ('filter' in variant and variant['filter'] == 'PASS') else fields[6]
                variant['vep_annotations'] = vep_annotations

                old_allele_count = variant['allele_count'] if 'allele_count' in variant else 0
                variant['allele_count'] = int(info_field['AC'].split(',')[i]) + old_allele_count #changed from AC_Adj to AC
                if not variant['allele_count'] and variant['filter'] == 'PASS': variant['filter'] = 'AC_Adj0' # Temporary filter
                old_allele_num = variant['allele_num'] if 'allele_num' in variant else 0
                variant['allele_num'] = int(info_field['AN']) + old_allele_num #changed from AN_adj to AN

                if variant['allele_num'] > 0:
                    variant['allele_freq'] = variant['allele_count']/float(variant['allele_num']) #changed from AC_Adj to AN
                else:
                    variant['allele_freq'] = None
                    
                #for key, value in info_field.iteritems():
                 #   variant[key] = value

                #variant['pop_acs'] = dict([(POPS[x], int(info_field['AC_%s' % x].split(',')[i])) for x in POPS])
                #variant['pop_ans'] = dict([(POPS[x], int(info_field['AN_%s' % x])) for x in POPS])
                #variant['pop_homs'] = dict([(POPS[x], int(info_field['Hom_%s' % x].split(',')[i])) for x in POPS])
                #variant['ac_male'] = info_field['AC_MALE']
                #variant['ac_female'] = info_field['AC_FEMALE']
                #variant['an_male'] = info_field['AN_MALE']
                #variant['an_female'] = info_field['AN_FEMALE']


                sample_infos = get_genotype(i,line,sample_names) #Todo split the ith allele?
                samples = sample_infos[0]
                hom_count = sample_infos[1]
                het_count = sample_infos[2]

                old_hom_count = variant['hom_count'] if 'hom_count' in variant else 0
                variant['hom_count'] = hom_count + old_hom_count
                old_het_count = variant['hemi_count'] if'hemi_count' in variant else 0
                variant['hemi_count'] = het_count + old_het_count
                variant['quality_metrics'] = dict([(x, info_field[x]) for x in METRICS if x in info_field])

                variant['genes'] = list({annotation['Gene'] for annotation in vep_annotations})
                variant['transcripts'] = list({annotation['Feature'] for annotation in vep_annotations})

                #custom annotations
                exac_dict = dict()
                #domains_dict = dict()
                for single_annotation in vep_annotations:
                    exac_dict.update([(key, value) for key, value in single_annotation.iteritems() if key.startswith("ExAC")])
                #    domains_dict.update([(key, value) for key, value in single_annotation.iteritems() if key.startswith("DOMAINS")])

                for key,value in exac_dict.iteritems():
                        value = value.split("&")[min(i,len(value.split("&")) -1)]
                        value = value.split(":")[len(value.split(":"))-1]
                        if value == '':
                            variant[key] = value
                        else:
                            variant[key] = float(value)

                #population annotations
                cohort = get_cohort(variant,cohort_name)
                cohort['name'] = cohort_name
                old_co_hom_count = cohort['hom_count'] if 'hom_count' in cohort else 0
                cohort['hom_count'] = hom_count + old_co_hom_count
                old_co_hemi_count = cohort['hemi_count'] if 'hemi_count' in cohort else 0
                cohort['hemi_count'] = het_count + old_co_hemi_count

                cohort['filter'] = 'PASS' if fields[6] == 'PASS' or ('filter' in cohort and cohort['filter'] == 'PASS') else fields[6]
                old_co_allele_count = cohort['allele_count'] if 'allele_count' in cohort else 0
                cohort['allele_count'] = int(info_field['AC'].split(',')[i]) + old_co_allele_count
                if not cohort['allele_count'] and cohort['filter'] == 'PASS': cohort['filter'] = 'AC_Adj0' # Temporary filter
                old_co_allele_num = cohort['allele_num'] if 'allele_num' in cohort else 0
                cohort['allele_num'] = int(info_field['AN']) + old_co_allele_num

                if cohort['allele_num'] > 0:
                    cohort['allele_freq'] = cohort['allele_count']/ float(cohort['allele_num'])
                else:
                    cohort['allele_freq'] = None

                track = {}
                track['date'] = time.strftime("%d/%m/%Y")
                track['samples'] = samples
                tracks = cohort['tracks'] if 'tracks' in cohort else list()
                tracks.append(track)
                cohort['tracks'] = tracks
                cohorts = variant['cohorts'] if 'cohorts' in variant else list()
                cohorts.append(cohort)
                variant['cohorts'] = cohorts
                add_consequence_to_variant(variant);

                if 'DP_HIST' in info_field:
                    hists_all = [info_field['DP_HIST'].split(',')[0], info_field['DP_HIST'].split(',')[i+1]]
                    variant['genotype_depths'] = [zip(dp_mids, map(int, x.split('|'))) for x in hists_all]
                if 'GQ_HIST' in info_field:
                    hists_all = [info_field['GQ_HIST'].split(',')[0], info_field['GQ_HIST'].split(',')[i+1]]
                    variant['genotype_qualities'] = [zip(gq_mids, map(int, x.split('|'))) for x in hists_all]

                yield variant
        except Exception:
            print("Error parsing vcf line: " + line)
            traceback.print_exc()
            break


def get_variant_and_replace(chrom,pos,ref,alt):
    import exac
    db = exac.get_db()
    variants = db.variants
    return variants.find_one_and_delete({'chrom':chrom,'pos':pos,'ref':ref,'alt':alt})

def get_variant(chrom,pos,ref,alt):
    import exac
    db = exac.get_db()
    variants = db.variants
    return variants.find_one({'chrom':chrom,'pos':pos,'ref':ref,'alt':alt})

def get_cohort(variant,cohort_name):
    if 'cohorts' not in variant:
        return {}
    else:
        cohorts = list()
        current_cohort = {}
        for cohort in variant['cohorts']:
            if cohort['name'] == cohort_name:
                current_cohort =  cohort
            else:
                cohorts.append(cohort)
        variant['cohorts'] = cohorts
        return current_cohort


def get_genotype(i,line,sample_names):
    format_keys = line.split('\t')[8]
    sample_values = line.split('\t')[9::]
    samples = list()
    hom_count = 0
    het_count = 0
    for key, value in zip(sample_names, sample_values):
        sample = {}
        sample['name'] = key
        me = re.search('_exome_(.*)_hg19',key)
        library_number = key
        if me:
            library_number = me.group(1)
        sample['library_number'] = library_number
        format_values = value.split(":")
        for formatkey, formatvalue in zip(format_keys.split(":"),format_values):
            sample[formatkey] = formatvalue
            if formatkey == 'GT':
                if formatvalue == "1/1":
                    formatvalue = "homozygous"
                    hom_count += 1
                elif formatvalue == "0/1" or value == "1/0":
                    formatvalue = "heterozygous"
                    het_count += 1
                elif formatvalue == "0/0":
                    formatvalue = "both_ref"
                else:
                    formatvalue = "no_call"
                sample[formatkey+'_format'] = formatvalue
        samples.append(sample)

    return (samples,hom_count,het_count)


def get_annotations_vcf_ikmb(sites_vcf):
    """
    Parse exac sites VCF file and return iter of variant dicts
    sites_vcf is a file (gzipped), not file path
    """
    getcontext().prec = 5
    for line in sites_vcf:
        try:
            line = line.strip('\n')
            if line.startswith('#'):
                continue

            fields = line.split('\t')
            alt_alleles = fields[4].split(',')
            chrom = fields[0].replace("chr","")
            info_field = dict([(x.split('=', 1)) if '=' in x else (x, x) for x in re.split(';(?=\w)', fields[7])])

            # different variant for each alt allele
            for i, alt_allele in enumerate(alt_alleles):
                # If we get here, it's a variant line
                pos, ref, alt = get_minimal_representation(fields[1], fields[3], alt_allele)
                variant = get_variant_and_replace(chrom,pos,ref,alt)


#1. Frequenzen: HRC, ESP, ExAC, Kaviar, 1000G (und rs-Nummern aus dbSNP)
#2. prediction tools: DANN, FATHMM, CADD (zum Teil auch genomweit verfuegbar)
#3. conservation scores: GERP, PhyloP
#Ausserdem noch dbscSNV (zur splice-site prediction), interpro-domain (Info ueber betroffene Proteindomaenen, genaue Ressource auf die ANNOVAR zugreift muesste ich nochmal nachschauen) und CLINVAR.

                #parse frequencies
                hrc_dict = dict()
                hrc_dict.update([(key, value) for key, value in info_field.iteritems() if key.startswith("HRC")])
                for key,value in hrc_dict.iteritems():
                    value = value if value != '.' else ""
                    if value == '':
                        try:
                            variant[key] = value
                        except TypeError:
                            print variant
                            print key
                            print value
                    else:
                        variant[key] = float(value)

                esp6500siv2_all = info_field['esp6500siv2_all'] if info_field['esp6500siv2_all'] != '.' else ""
                variant['ESP'] = float(esp6500siv2_all) if esp6500siv2_all!='' else esp6500siv2_all

                #exac03nontcga = info_field['exac03nontcga'] if info_field['exac03nontcga'] != '.' else ""
                #variant['ExAC'] = exac03nontcga

                kav_dict = dict()
                kav_dict.update([(key, value) for key, value in info_field.iteritems() if key.startswith("Kaviar")])
                for key,value in kav_dict.iteritems():
                    value = value if value != '.' else ""
                    if value == '':
                        variant[key] = value
                    else:
                        variant[key] = float(value)

                kgenomes = info_field['1000g2014oct_all'] if info_field['1000g2014oct_all'] != '.' else ""
                variant['g1k'] = float(kgenomes) if kgenomes!='' else kgenomes

                #predictions scores
                dann_gw = info_field['DANN_score'] if info_field['DANN_score'] != '.' else ""
                variant['DANN'] = float(dann_gw) if dann_gw!='' else dann_gw

                fathmm_gw = info_field['FATHMM_score'] if info_field['FATHMM_score'] != '.' else ""
                variant['FATHMM'] = float(fathmm_gw) if fathmm_gw!='' else fathmm_gw

                cadd_gw = info_field['CADD_raw'] if info_field['CADD_raw'] != '.' else ""
                variant['CADD'] = float(cadd_gw) if cadd_gw!='' else cadd_gw

                #conservation scores
                #GERP++_RS=2.31;phyloP46way_placental=0.267;phyloP100way_vertebrate=1.636;SiPhy_29way_logOdds=7.538
                gerp = info_field['GERP++_RS'] if info_field['GERP++_RS'] != '.' else ""
                variant['GERP'] = float(gerp) if gerp!='' else gerp

                phylo_placental = info_field['phyloP46way_placental'] if info_field['phyloP46way_placental'] != '.' else ""
                variant['phylo_placental'] = float(phylo_placental) if phylo_placental!='' else phylo_placental

                pyhlo_vertebrate = info_field['phyloP100way_vertebrate'] if info_field['phyloP100way_vertebrate'] != '.' else ""
                variant['pyhlo_vertebrate'] = float(pyhlo_vertebrate) if pyhlo_vertebrate!='' else pyhlo_vertebrate

                siPhy = info_field['SiPhy_29way_logOdds'] if info_field['SiPhy_29way_logOdds'] != '.' else ""
                variant['SiPhy'] = float(siPhy) if siPhy!='' else siPhy

                #additional
                dbsnv_dict = dict()
                dbsnv_dict.update([(key, value) for key, value in info_field.iteritems() if key.startswith("dbscSNV")])
                for key,value in dbsnv_dict.iteritems():
                    value = value if value != '.' else ""
                    variant[key] = float(value) if value!='' else value

                dbnsfp31a_interpro = info_field['Interpro_domain'] if info_field['Interpro_domain'] != '.' else ""
                variant['interpro_domain'] = dbnsfp31a_interpro

                #CLINSIG=.;CLNDBN=.;CLNACC=.;CLNDSDB=.;CLNDSDBID=
                CLINSIG = info_field['CLINSIG'] if info_field['CLINSIG'] != '.' else ""
                variant['CLINSIG'] = CLINSIG

                CLNDBN = info_field['CLNDBN'] if info_field['CLNDBN'] != '.' else ""
                variant['CLNDBN'] = CLNDBN

                CLNACC = info_field['CLNACC'] if info_field['CLNACC'] != '.' else ""
                variant['CLNACC'] = CLNACC

                CLNDSDB = info_field['CLNDSDB'] if info_field['CLNDSDB'] != '.' else ""
                variant['CLNDSDB'] = CLNDSDB

                CLNDSDBID = info_field['CLNDSDBID'] if info_field['CLNDSDBID'] != '.' else ""
                variant['CLNDSDBID'] = CLNDSDBID

                yield variant

        except Exception:
            print("Error parsing vcf line: " + line)
            traceback.print_exc()
            break



def get_mnp_data(mnp_file):
    header = mnp_file.readline().strip().split('\t')
    for line in mnp_file:
        data = dict(zip(header, line.split('\t')))
        if any(map(lambda x: x == 'True', data['QUESTIONABLE_PHASING'])): continue
        chroms = data['CHROM'].split(',')
        chrom = chroms[0]
        sites = data['SITES'].split(',')
        refs = data['REF'].split(',')
        alts = data['ALT'].split(',')
        for i, site in enumerate(sites):
            all_sites = zip(chroms, sites, refs, alts)
            all_sites.remove(all_sites[i])
            mnp = {}
            mnp['xpos'] = get_xpos(chrom, site)
            mnp['ref'] = refs[i]
            mnp['alt'] = alts[i]
            mnp['site2'] = '-'.join(all_sites[0])
            if len(all_sites) > 1:
                mnp['site3'] = all_sites[1]
            mnp['combined_codon_change'] = data['COMBINED_CODON_CHANGE']
            mnp['category'] = data['CATEGORY']
            mnp['number_samples'] = data['NSAMPS']
            yield mnp


def get_constraint_information(constraint_file):
    _, _, _, header = constraint_file.readline().strip().split(None, 3)
    header = header.split()
    for line in constraint_file:
        transcript, gene, chrom, info = line.strip().split(None, 3)
        transcript_info = dict(zip(header, map(float, info.split())))
        transcript_info['transcript'] = transcript.split('.')[0]
        yield transcript_info


def get_canonical_transcripts(canonical_transcript_file):
    for line in canonical_transcript_file:
        gene, transcript = line.strip().split()
        yield gene, transcript


def get_omim_associations(omim_file):
    for line in omim_file:
        fields = line.strip().split('\t')
        if len(fields) == 4:
            yield fields
        else:
            yield None


def get_genes_from_gencode_gtf(gtf_file):
    """
    Parse gencode GTF file;
    Returns iter of gene dicts
    """
    for line in gtf_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')

        if fields[2] != 'gene':
            continue

        chrom = fields[0][3:]
        start = int(fields[3]) + 1  # bed files are 0-indexed
        stop = int(fields[4]) + 1
        info = dict(x.strip().split() for x in fields[8].split(';') if x != '')
        info = {k: v.strip('"') for k, v in info.items()}
        gene_id = info['gene_id'].split('.')[0]

        gene = {
            'gene_id': gene_id,
            'gene_name': info['gene_name'],
            'gene_name_upper': info['gene_name'].upper(),
            'chrom': chrom,
            'start': start,
            'stop': stop,
            'strand': fields[6],
            'xstart': get_xpos(chrom, start),
            'xstop': get_xpos(chrom, stop),
        }
        yield gene


def get_transcripts_from_gencode_gtf(gtf_file):
    """
    Parse gencode GTF file;
    Returns iter of transcript dicts
    """
    for line in gtf_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')

        if fields[2] != 'transcript':
            continue

        chrom = fields[0][3:]
        start = int(fields[3]) + 1  # bed files are 0-indexed
        stop = int(fields[4]) + 1
        info = dict(x.strip().split() for x in fields[8].split(';') if x != '')
        info = {k: v.strip('"') for k, v in info.items()}
        transcript_id = info['transcript_id'].split('.')[0]
        gene_id = info['gene_id'].split('.')[0]

        gene = {
            'transcript_id': transcript_id,
            'gene_id': gene_id,
            'chrom': chrom,
            'start': start,
            'stop': stop,
            'strand': fields[6],
            'xstart': get_xpos(chrom, start),
            'xstop': get_xpos(chrom, stop),
        }
        yield gene


def get_exons_from_gencode_gtf(gtf_file):
    """
    Parse gencode GTF file;
    Returns iter of transcript dicts
    """
    for line in gtf_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')

        if fields[2] not in ['exon', 'CDS', 'UTR']:
            continue

        chrom = fields[0][3:]
        feature_type = fields[2]
        start = int(fields[3]) + 1  # bed files are 0-indexed
        stop = int(fields[4]) + 1
        info = dict(x.strip().split() for x in fields[8].split(';') if x != '')
        info = {k: v.strip('"') for k, v in info.items()}
        transcript_id = info['transcript_id'].split('.')[0]
        gene_id = info['gene_id'].split('.')[0]

        exon = {
            'feature_type': feature_type,
            'transcript_id': transcript_id,
            'gene_id': gene_id,
            'chrom': chrom,
            'start': start,
            'stop': stop,
            'strand': fields[6],
            'xstart': get_xpos(chrom, start),
            'xstop': get_xpos(chrom, stop),
        }
        yield exon


def get_cnvs_from_txt(cnv_txt_file):
    """                                                                                                                                                                                                                                  
    Parse gencode txt file;                                                                                                                                                                                                              
    Returns iter of gene dicts                                                                                                                                                                                                           
    """
    header = cnv_txt_file.next() # gets rid of the header                                                                                                                                                                                
    #print header                                                                                                                                                                                                                        
    for line in cnv_txt_file:

        fields = line.rsplit()
        transcript = fields[0]
        gene = fields[1]
        chrom = fields[2]
        start = int(fields[3])
        stop = int(fields[4])
        del0 = int(fields[5])
        del60 = int(fields[6])
        dup0 = int(fields[7])
        dup60 = int(fields[8])
        delpop0 = fields[9]
        delpop60 = fields[10]
        duppop0 = fields[11]
        duppop60 = fields[12]
        

        #find gene from DB.genes, get ID                                                                                                                                                                                                 
        #find exon of that gene that this CNV referes to from db.exons, get ID                                                                                                                                                           
        #add the object reference to the cnv dict.                                                                                                                                                                                       
        cnv = {
            'transcript': transcript,
            'gene': gene,
            'chrom': chrom,
            'start': start,
            'stop': stop,
            'del0': del0,
            'dup0': dup0,
            'dup60': dup60,
            'del60' : del60,
            'delpop0' : delpop0,
            'delpop60' : delpop60,
            'duppop0' : duppop0,
            'duppop60' : duppop60,
            'xstart': get_xpos(chrom, start),
            'xstop': get_xpos(chrom, stop),
        }
        yield cnv


def get_cnvs_per_gene(cnv_gene_file):
    header = cnv_gene_file.next() # gets rid of the header                                                                                                                                                                               
    for line in cnv_gene_file:

        fields = line.rsplit()
        gene = fields[0]
        symbol = fields[1]
        del0 = int(fields[2])
        dup0 = int(fields[3])
        cnv0 = int(fields[4])
        del60 = int(fields[5])
        dup60 = int(fields[6])
        cnv60 = int(fields[7])
        del_score = float(fields[8])
        dup_score = float(fields[9])
        cnv_score = float(fields[10])
        rank = int(fields[11])

        cnv_gene = {
            'gene': gene,
            'symbol': symbol,
            'del0': del0,
            'dup0': dup0,
            'cnv0': cnv0,
            'del60': del60,
            'dup60': dup60,
            'cnv60' : cnv60,
            'del_score': del_score,
            'dup_score': dup_score,
            'cnv_score': cnv_score,
            'rank': rank,
            }
        yield cnv_gene




def get_dbnsfp_info(dbnsfp_file):
    """
    Parse dbNSFP_gene file;
    Returns iter of transcript dicts
    """
    header = dbnsfp_file.next().split('\t')
    fields = dict(zip(header, range(len(header))))
    for line in dbnsfp_file:
        line = line.split('\t')
        other_names = line[fields["Gene_old_names"]].split(';') if line[fields["Gene_old_names"]] != '.' else []
        if line[fields["Gene_other_names"]] != '.':
            other_names.extend(line[fields["Gene_other_names"]].split(';'))
        gene_info = {
            'gene_name': line[fields["Gene_name"]],
            'ensembl_gene': line[fields["Ensembl_gene"]],
            'gene_full_name': line[fields["Gene_full_name"]],
            'gene_other_names': other_names
        }
        yield gene_info


def get_snp_from_dbsnp_file(dbsnp_file):
    for line in dbsnp_file:
        fields = line.split('\t')
        if len(fields) < 3: continue
        rsid = int(fields[0])
        chrom = fields[1].rstrip('T')
        if chrom == 'PAR': continue
        start = int(fields[2]) + 1
        snp = {
            'xpos': get_xpos(chrom, start),
            'rsid': rsid
        }
        yield snp
