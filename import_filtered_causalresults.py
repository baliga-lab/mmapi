#!/usr/bin/env python3

import argparse
import MySQLdb
import pandas as pd
import numpy as np
from collections import defaultdict


DESCRIPTION = "import_filtered_causalresults - import filtered causal results to MM backend"

"""
Tables in mm_api_v2:


ignored for now (because patients are needed):

patient_tf
patients
bicluster_patients
bicluster_boxplot_data

bc_tf_roles: 1 = activates, 2 = represses
bc_mutation_tf_roles: 1 = down-regulates 2 = up-regulates
"""

def fill_roles(conn):
    with conn.cursor() as cur:
        cur.execute('select count(*) from bc_tf_roles')
        # only insert if not exist
        if cur.fetchone()[0] == 0:
            cur.execute("insert into bc_tf_roles (id, name) values (1, 'activates')")
            cur.execute("insert into bc_tf_roles (id, name) values (2, 'represses')")
            cur.execute("insert into bc_mutation_tf_roles (id, name) values (1, 'down-regulates')")
            cur.execute("insert into bc_mutation_tf_roles (id, name) values (2, 'up-regulates')")
            conn.commit()

def import_tfs(conn, df):
    """tfs: id: int, name: string, cox_hazard_ratio: float
    In the import document this is the Regulator"""
    tfs = {}
    for index, row in df.iterrows():
        tf = row['Regulator']
        cox_hr = row['HazardRatio']
        tfs[tf] = cox_hr
    result = {}
    with conn.cursor() as cur:
        cur.execute('select count(*) from tfs')
        num_tfs = cur.fetchone()[0]
        if num_tfs > 0:
            print('TFs found, reading from database')
            cur.execute('select id,name,cox_hazard_ratio from tfs')
            for pk, name, cox_hr in cur.fetchall():
                result[name] = (pk, cox_hr)
        else:
            print('TFs not found, inserting into database')
            for tf, cox_hr in tfs.items():
                cur.execute('insert into tfs (name,cox_hazard_ratio) values (%s,%s)',
                            [tf, cox_hr])
                result[tf] = (cur.lastrowid, cox_hr)
            conn.commit()
    return result


def import_mutations(conn, df):
    """mutations: id: int, name: string"""
    mutations = sorted(set(df["Mutation"]))

    result = {}
    with conn.cursor() as cur:
        cur.execute('select count(*) from mutations')
        if cur.fetchone()[0] > 0:
            print("Mutations found, reading from database")
            cur.execute('select id,name from mutations')
            for pk, name in cur.fetchall():
                result[name] = pk
        else:
            print("Mutations not found, inserting into database")
            for name in mutations:
                cur.execute('insert into mutations (name) values (%s)', [name])
                result[name] = cur.lastrowid
            conn.commit()
    return result


def import_genes(conn, df, ens2pref, ens2entrez):
    """genes: id: int, ensembl_id: string, entrez_id: string, preferred: string
    In the import document this is the Regulator
    """
    genes = set()
    for index, row in df.iterrows():
        row_genes = row['RegulonGenes'].split(',')
        genes.update(row_genes)
    print("# unique genes found: %d" % len(genes))

    result = {}
    with conn.cursor() as cur:
        cur.execute('select count(*) from genes')
        if cur.fetchone()[0] > 0:
            print("Genes found reading from database")
            cur.execute('select id,ensembl_id,entrez_id,preferred from genes')
            for pk, ens, entrez, pref in cur.fetchall():
                result[ens] = (pk, entrez, pref)
        else:
            print("Genes not found, inserting into database")
            for ens in sorted(genes):
                try:
                    entrez = ens2entrez[ens]
                except:
                    entrez = None
                try:
                    pref = ens2pref[ens]
                except:
                    pref = None
                cur.execute('insert into genes (ensembl_id,entrez_id,preferred) values (%s,%s,%s)',
                            [ens, entrez, pref])
                result[ens] = (cur.lastrowid, entrez, pref)
            conn.commit()
    return result


def import_regulons(conn, df, mutations):
    """biclusters: id, name, cox_hazard_ratio"""
    result = {}
    with conn.cursor() as cur:
        cur.execute('select count(*) from biclusters')
        if cur.fetchone()[0] > 0:
            print('Regulons found, reading from database')
            cur.execute('select b.id,b.name,b.cox_hazard_ratio from biclusters b')
            for pk, regulon, hr in cur.fetchall():
                result[regulon] = (pk, hr)
        else:
            print('Regulons not found, inserting into database')
            for index, row in df.iterrows():
                regulon = row['Regulon']
                cox_hr = row['HazardRatio']
                cur.execute('select count(*) from biclusters where name=%s', [regulon])
                if cur.fetchone()[0] == 0:
                    cur.execute('insert into biclusters (name,cox_hazard_ratio) values (%s,%s)',
                                [regulon, cox_hr])
                    result[regulon] = (cur.lastrowid, cox_hr)
            conn.commit()
    return result


def import_mutation_regulator(conn, df, regulons, mutations, tfs):
    """table: bc_mutation_tf, field: MutationRegulatorEdge
    id, bicluster_id, mutation_id, tf_id, role"""
    """bc_tf_roles: 1 = activates, 2 = represses
    bc_mutation_tf_roles: 1 = down-regulates 2 = up-regulates"""
    with conn.cursor() as cur:
        cur.execute('select count(*) from bc_mutation_tf')
        if cur.fetchone()[0] == 0:
            print('no mutation_regulator edges found, insert into database')

            for index, row in df.iterrows():
                regulon = row['Regulon']
                regulator = row['Regulator']
                mutation = row['Mutation']
                edge = row['MutationRegulatorEdge']
                regulon_id, cox_hr1  = regulons[regulon]
                regulator_id, cox_hr2 = tfs[regulator]
                mutation_id = mutations[mutation]
                if edge < 0:
                    role = 1  # down-regulates
                else:
                    role = 2  # up-regulates
                cur.execute('insert into bc_mutation_tf (bicluster_id,mutation_id,tf_id,role) values (%s,%s,%s,%s)',
                            [regulon_id, mutation_id, regulator_id, role])
            conn.commit()
        else:
            print('skip inserting mutation regulator edges')

"""
Index(['Unnamed: 0', 'Mutation', 'Regulator', 'Regulon',
       'MutationRegulatorEdge', '-log10(p)_MutationRegulatorEdge',
       'RegulatorRegulon_Spearman_R', 'RegulatorRegulon_Spearman_p-value',
       'Regulon_stratification_t-statistic',
       '-log10(p)_Regulon_stratification',
       'Fraction_of_edges_correctly_aligned', 'TranscriptionalProgram',
       'HazardRatio', 'HazardRatioPval', 'RegulonGenes', 'DrugEnrichment',
       'HallmarksEnrichment', 'LinHallmarks', 'TargetClassEnrichment',
       'TargetClass_p-value', 'MechanismOfActionEnrichment',
       'MechanismOfAction_p-value', 'mirv_miRNA', 'PITA_miRNA', 'PITA_pval',
       'TargetScan_miRNA', 'TargetScan_pval'],
      dtype='object')

"""

def import_regulon_regulator(conn, df, regulons, tfs):
    """table: bc_tf (bicluster_id,tf_id,role)"""
    with conn.cursor() as cur:
        cur.execute('select count(*) from bc_tf')
        if cur.fetchone()[0] == 0:
            print('insert regulon-regulator into database')
            for index, row in df.iterrows():
                regulon = row['Regulon']
                regulator = row['Regulator']
                regulon_id = regulons[regulon][0]
                regulator_id = tfs[regulator][0]
                spearman_r = row['RegulatorRegulon_Spearman_R']
                if spearman_r > 0:
                    role = 1  # 1 activates
                else:
                    role = 2  # 2 represses
                #print(row['RegulatorRegulon_Spearman_p-value'])
                cur.execute('select count(*) from bc_tf where bicluster_id=%s and tf_id=%s',
                            [regulon_id, regulator_id])
                if cur.fetchone()[0] == 0:
                    # only insert once
                    cur.execute('insert into bc_tf (bicluster_id,tf_id,role) values (%s,%s,%s)',
                                [regulon_id, regulator_id, role])
            conn.commit()
        else:
            print('regulon regulator relations exist, skip')


def import_regulon_genes(conn, df, regulons, genes):
    """table: bicluster_genes (bicluster_id, gene_id)"""
    with conn.cursor() as cur:
        cur.execute('select count(*) from bicluster_genes')
        if cur.fetchone()[0] == 0:
            print('insert regulon genes into database')
            for index, row in df.iterrows():
                regulon = row['Regulon']
                regulon_genes = set(row['RegulonGenes'].strip().split(','))
                regulon_id = regulons[regulon][0]
                for gene in regulon_genes:
                    gene_id = genes[gene][0]
                    cur.execute('select count(*) from bicluster_genes where bicluster_id=%s and gene_id=%s',
                                [regulon_id, gene_id])
                    if cur.fetchone()[0] == 0:  # prevent double insertions in the same transaction
                        cur.execute('insert into bicluster_genes(bicluster_id,gene_id) values (%s,%s)',
                                    [regulon_id, gene_id])
            conn.commit()
        else:
            print('regulon genes exist, skip')


def import_transcriptional_programs(conn, df, regulons):
    print('import transcriptional programs')
    with conn.cursor() as cur:
        for index, row in df.iterrows():
            regulon = row['Regulon']
            regulon_id = regulons[regulon][0]
            prog = row['TranscriptionalProgram']
            prog = int(prog)
            cur.execute('update biclusters set trans_program=%s where id=%s', [prog, regulon_id])
        conn.commit()


def import_drug_enrichment(conn, df, regulons):
    print('import drug enrichment')
    regulon2drugs = defaultdict(set)
    all_drugs = set()
    with conn.cursor() as cur:
        cur.execute('select count(*) from regulon_drugs')
        if cur.fetchone()[0] > 0:
            print("drug enrichment already exists, skipping")
            return
        for index, row in df.iterrows():
            regulon = row['Regulon']
            regulon_id = regulons[regulon][0]
            drugs = row['DrugEnrichment']
            try:
                if not np.isnan(drugs):
                    # drugs are text
                    drugs = drugs.split(',')
                else:
                    drugs = set()
            except:
                drugs = drugs.split(',')

            regulon2drugs[regulon_id].update(drugs)
            all_drugs.update(drugs)

        # insert all the drugs first
        drug_ids = {}
        for drug in all_drugs:
            cur.execute('insert into drugs (name) values (%s)', [drug])
            pk = cur.lastrowid
            drug_ids[drug] = pk

        for regulon_id, drugs in regulon2drugs.items():
            for drug in drugs:
                drug_id = drug_ids[drug]
                cur.execute('insert into regulon_drugs (regulon_id, drug_id) values (%s,%s)',
                            [regulon_id, drug_id])
        conn.commit()


def import_target_class_enrichment(conn, df, regulons):
    print('import target class enrichment')

    regulon2targetclass = defaultdict(set)
    all_targetclasses = set()
    with conn.cursor() as cur:
        cur.execute('select count(*) from regulon_target_class')
        if cur.fetchone()[0] > 0:
            print('regulon target class enrichment found, skipping')
            return

        for index, row in df.iterrows():
            regulon = row['Regulon']
            regulon_id = regulons[regulon][0]
            enrich = row['TargetClassEnrichment']
            pval = row['TargetClass_p-value']
            try:
                if np.isnan(enrich):
                    enrich = None
            except:
                pass
            if enrich is not None:
                all_targetclasses.add(enrich)
                regulon2targetclass[regulon_id].add((enrich, pval))

        tc_ids = {}
        for tc in all_targetclasses:
            cur.execute('insert into target_classes (name) values (%s)', [tc])
            tc_ids[tc] = cur.lastrowid
        for regulon_id, target_classes in regulon2targetclass.items():
            for tc, pval in target_classes:
                tc_id = tc_ids[tc]
                cur.execute('insert into regulon_target_class (regulon_id,target_class_id,pval) values (%s,%s,%s)',
                            [regulon_id, tc_id, pval])
        conn.commit()


def parse_list(s):
    try:
        if not np.isnan(s):
            result = s.split(',')
        else:
            result = []
    except:
        result = s.split(',')
    return result


def import_hallmarks_enrichment(conn, df, regulons):
    print('import hallmarks enrichment')
    reg_hallmarks = defaultdict(set)
    reg_lin = defaultdict(set)
    all_hallmarks = set()
    all_lins = set()
    hallmark_ids = {}
    lin_hallmark_ids = {}

    with conn.cursor() as cur:
        cur.execute('select count(*) from regulon_hallmarks')
        if cur.fetchone()[0] > 0:
            print('entries exist, skipping')
            return

        for index, row in df.iterrows():
            regulon = row['Regulon']
            regulon_id = regulons[regulon][0]
            enrich = row['HallmarksEnrichment']
            lin_enrich = row['LinHallmarks']
            enrich = parse_list(enrich)
            lin_enrich = parse_list(lin_enrich)

            all_hallmarks.update(enrich)
            all_lins.update(lin_enrich)

            reg_hallmarks[regulon_id].update(enrich)
            reg_lin[regulon_id].update(lin_enrich)

        for hm in all_hallmarks:
            cur.execute('insert into hallmarks (name) values (%s)', [hm])
            hallmark_ids[hm] = cur.lastrowid

        for hm in all_lins:
            cur.execute('insert into lin_hallmarks (name) values (%s)', [hm])
            lin_hallmark_ids[hm] = cur.lastrowid

        for regulon_id, hallmarks in reg_hallmarks.items():
            for hm in hallmarks:
                hm_id = hallmark_ids[hm]
                cur.execute('insert into regulon_hallmarks (regulon_id,hallmark_id) values (%s,%s)',
                            [regulon_id, hm_id])

        for regulon_id, hallmarks in reg_lin.items():
            for hm in hallmarks:
                hm_id = lin_hallmark_ids[hm]
                cur.execute('insert into regulon_lin_hallmarks (regulon_id,hallmark_id) values (%s,%s)',
                            [regulon_id, hm_id])
        conn.commit()


def import_mechanism_of_action_enrichment(conn, df, regulons):
    print('import mechanism of action enrichment')

    all_moas = set()
    regulon_moas = defaultdict(set)
    moa_ids = {}

    with conn.cursor() as cur:
        cur.execute('select count(*) from regulon_mechanism_of_action')
        if cur.fetchone()[0] > 0:
            print('entries found, skipping')
            return

        for index, row in df.iterrows():
            regulon = row['Regulon']
            regulon_id = regulons[regulon][0]
            enrich = row['MechanismOfActionEnrichment']
            enrich = parse_list(enrich)
            all_moas.update(enrich)
            regulon_moas[regulon_id].update(enrich)

        for moa in all_moas:
            cur.execute('insert into mechanisms_of_action (name) values (%s)', [moa])
            moa_ids[moa] = cur.lastrowid

        for regulon_id, moas in regulon_moas.items():
            for moa in moas:
                moa_id = moa_ids[moa]
                cur.execute('insert into regulon_mechanism_of_action (regulon_id,mechanism_of_action_id) values (%s,%s)',
                            [regulon_id, moa_id])

        conn.commit()


DBNAME = 'mm_api_v3'

HEADERS = [
    "Mutation", "Regulator", "Regulon", "MutationRegulatorEdge", "-log10(p)_MutationRegulatorEdge",
    "RegulatorRegulon_Spearman_R", "RegulatorRegulon_Spearman_p-value",
    "Regulon_stratification_t-statistic-log10(p)_Regulon_stratification	Fraction_of_edges_correctly_aligned",
    "TranscriptionalProgram", "HazardRatio", "HazardRatioPval", "RegulonGenes",
    "DrugEnrichment",
    "HallmarksEnrichment", "LinHallmarks",
    "TargetClassEnrichment", "TargetClass_p-value",
    "MechanismOfActionEnrichment", "MechanismOfAction_p-valuemirv_miRNA",
    "PITA_miRNA", "PITA_pval", "TargetScan_miRNA", "TargetScan_pval"
]

def dbconn():
    return MySQLdb.connect(host="localhost", user="root", passwd="root",
                           db=DBNAME)


def read_synonyms(idconv_file):
    df = pd.read_csv(idconv_file, sep='\t')
    ens2pref = {}
    ens2entrez = {}
    for index, row in df.iterrows():
        ensembl = row['ensembl']
        try:
            entrez = int(row['entrez'].strip())
            ens2entrez[ensembl] = entrez
        except:
            pass
        try:
            preferred = row['preferred']
            ens2pref[ensembl] = preferred
        except:
            pass
    return ens2pref, ens2entrez


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('infile', help="input file in CSV format")
    parser.add_argument('idconv', help="id conversion file")
    args = parser.parse_args()
    conn = dbconn()

    df = pd.read_csv(args.infile, sep='\t')
    #print(df['MutationRegulatorEdge'])

    ens2pref, ens2entrez = read_synonyms(args.idconv)
    #print(df.columns)
    # Step 1: collect the Regulator field
    #print(set(df['Regulator']))
    #print(set(df['RegulonGenes']))
    fill_roles(conn)
    tfs = import_tfs(conn, df)
    mutations = import_mutations(conn, df)
    genes = import_genes(conn, df, ens2pref, ens2entrez)
    regulons = import_regulons(conn, df, mutations)
    import_mutation_regulator(conn, df, regulons, mutations, tfs)
    import_regulon_regulator(conn, df, regulons, tfs)
    import_regulon_genes(conn, df, regulons, genes)
    import_transcriptional_programs(conn, df, regulons)
    import_drug_enrichment(conn, df, regulons)
    import_target_class_enrichment(conn, df, regulons)
    import_hallmarks_enrichment(conn, df, regulons)
    import_mechanism_of_action_enrichment(conn, df, regulons)
