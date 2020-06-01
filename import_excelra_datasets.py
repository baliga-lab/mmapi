#!/usr/bin/env python3
import os, sys
from collections import defaultdict
import MySQLdb
import traceback

DBNAME = 'mm_api_v3'

def dbconn():
    return MySQLdb.connect(host="localhost", user="root", passwd="root",
                           db=DBNAME)

def add_or_get_cancer(cur, cancer):
    cur.execute('select count(*) from exc_cancers where name=%s', [cancer])
    if cur.fetchone()[0] == 0:
        cur.execute('insert into exc_cancers (name) values (%s)', [cancer])
        return cur.lastrowid
    else:
        cur.execute('select id from exc_cancers where name=%s', [cancer])
        return cur.fetchone()[0]


def add_or_get_mutation(cur, mutation):
    cur.execute('select count(*) from exc_mutations where name=%s', [mutation])
    if cur.fetchone()[0] == 0:
        cur.execute('insert into exc_mutations (name) values (%s)', [mutation])
        return cur.lastrowid
    else:
        cur.execute('select id from exc_mutations where name=%s', [mutation])
        return cur.fetchone()[0]


def add_or_get_disease(cur, disease):
    cur.execute('select count(*) from exc_diseases where name=%s', [disease])
    if cur.fetchone()[0] == 0:
        cur.execute('insert into exc_diseases (name) values (%s)', [disease])
        return cur.lastrowid
    else:
        cur.execute('select id from exc_diseases where name=%s', [disease])
        return cur.fetchone()[0]


def add_or_get_regulator(cur, regulator):
    cur.execute('select count(*) from exc_regulators where name=%s', [regulator])
    if cur.fetchone()[0] == 0:
        cur.execute('insert into exc_regulators (name) values (%s)', [regulator])
        return cur.lastrowid
    else:
        cur.execute('select id from exc_regulators where name=%s', [regulator])
        return cur.fetchone()[0]


def add_or_get_regulon(cur, regulon):
    cur.execute('select count(*) from exc_regulons where name=%s', [regulon])
    if cur.fetchone()[0] == 0:
        cur.execute('insert into exc_regulons (name) values (%s)', [regulon])
        return cur.lastrowid
    else:
        cur.execute('select id from exc_regulons where name=%s', [regulon])
        return cur.fetchone()[0]

def add_or_get_drug(cur, drug):
    cur.execute('select count(*) from exc_drugs where name=%s', [drug])
    if cur.fetchone()[0] == 0:
        cur.execute('insert into exc_drugs (name) values (%s)', [drug])
        return cur.lastrowid
    else:
        cur.execute('select id from exc_drugs where name=%s', [drug])
        return cur.fetchone()[0]


def read_cancer_mutation(datadir, conn, hr):
    """Return a PMID dataset that is indexed for both cancers and mutations"""

    with open(os.path.join(datadir, 'Mut_Cancer_HR_%d.txt' % hr)) as infile:
        header = infile.readline().strip().split('\t')
        cancer_col = header.index('Disease')
        mutation_col = header.index('Mutation')
        pmids_col = header.index('PMIDs')

        with conn.cursor() as cur:
            for line in infile:
                row = line.strip().split('\t')
                if len(row) > 1:
                    # synonyms for cancer
                    cancers = [c for c in row[cancer_col].strip('"').split(' or ')
                               if len(c.strip()) > 0]
                    mutation = row[mutation_col].strip()
                    cancer_pks = []
                    for cancer in cancers:
                        cancer_pk = add_or_get_cancer(cur, cancer)
                        cancer_pks.append(cancer_pk)

                    if len(mutation) > 0:
                        mutation_pk = add_or_get_mutation(cur, mutation)
                    else:
                        continue

                    pmids = [s.strip() for s in row[pmids_col].split('|') if len(s.strip()) > 0]
                    for pmid in pmids:
                        for cancer_pk in cancer_pks:
                            cur.execute('select count(*) from exc_cancer_mutation where hr=%s and cancer_id=%s and mutation_id=%s and pmid=%s',
                                        [hr, cancer_pk, mutation_pk, pmid])
                            if cur.fetchone()[0] == 0:
                                cur.execute('insert into exc_cancer_mutation (hr,cancer_id,mutation_id,pmid) values (%s,%s,%s,%s)',
                                            [hr, cancer_pk, mutation_pk, pmid])
            conn.commit()

def read_mm_mutation(datadir, conn, hr):
    """Also known as disease-mutation"""
    with open(os.path.join(datadir, 'Mut_MM_HR_%d.txt' % hr)) as infile:
        header = infile.readline().strip().split('\t')
        disease_col = header.index('Disease')
        pmids_col = header.index('PMIDs')
        try:
            mutation_col = header.index('Mutation')
        except:
            print("WARNING: Dataset 'Mut_MM_HR_%d.txt' has no 'Mutation' column" % hr)
            mutation_col = header.index('Mutations')

        with conn.cursor() as cur:
            for line in infile:
                row = line.split('\t')
                diseases = [d for d in row[disease_col].strip('"').split(' or ') if len(d.strip()) > 0]
                disease_pks = []
                for disease in diseases:
                    disease_pk = add_or_get_disease(cur, disease)
                    disease_pks.append(disease_pk)

                mutations = [m.strip() for m in row[mutation_col].strip().split(' or ')
                             if len(m.strip()) > 0]
                mutation_pks = []
                for mutation in mutations:
                    mutation_pk = add_or_get_mutation(cur, mutation)
                    mutation_pks.append(mutation_pk)

                pmids = [s.strip() for s in row[pmids_col].split('|') if len(s.strip()) > 0]
                for pmid in pmids:
                    for disease_pk in disease_pks:
                        for mutation_pk in mutation_pks:
                            cur.execute('select count(*) from exc_disease_mutation where hr=%s and disease_id=%s and mutation_id=%s and pmid=%s',
                                        [hr, disease_pk, mutation_pk, pmid])
                            if cur.fetchone()[0] == 0:
                                cur.execute('insert into exc_disease_mutation (hr,disease_id,mutation_id,pmid) values (%s,%s,%s,%s)',
                                            [hr, disease_pk, mutation_pk, pmid])
            conn.commit()


def read_mm_regulator(datadir, conn, hr):
    """also known as disease regulator"""
    with open(os.path.join(datadir, 'Regulator_MM_HR_%d.txt' % hr)) as infile:
        header = infile.readline().strip().split('\t')
        disease_col = header.index('Disease')
        regulator_col = header.index('Regulator')
        pmids_col = header.index('PMIDs')

        with conn.cursor() as cur:
            for line in infile:
                row = line.split('\t')
                if len(row) > 1:
                    diseases = [d for d in row[disease_col].strip('"').split(' or ')
                                if len(d.strip()) > 0]
                    disease_pks = []
                    for disease in diseases:
                        disease_pk = add_or_get_disease(cur, disease)
                        disease_pks.append(disease_pk)

                    regulator = row[regulator_col].strip()
                    if len(regulator) > 0:
                        regulator_pk = add_or_get_regulator(cur, regulator)

                    pmids = [s.strip() for s in row[pmids_col].split('|') if len(s.strip()) > 0]
                    for pmid in pmids:
                        for disease_pk in disease_pks:
                            cur.execute('select count(*) from exc_disease_regulator where hr=%s and disease_id=%s and regulator_id=%s and pmid=%s',
                                        [hr, disease_pk, regulator_pk, pmid])
                            if cur.fetchone()[0] == 0:
                                cur.execute('insert into exc_disease_regulator (hr,disease_id,regulator_id,pmid) values (%s,%s,%s,%s)',
                                            [hr, disease_pk, regulator_pk, pmid])
            conn.commit()


def read_mm_regulon(datadir, conn, hr):
    """also known as disease regulon"""
    with open(os.path.join(datadir, 'Regulon_MM_HR_%d.txt' % hr)) as infile:
        header = infile.readline().strip().split('\t')
        disease_col = header.index('Disease')
        regulon_col = header.index('Regulon')
        pmids_col = header.index('PMIDs')

        with conn.cursor() as cur:
            for line in infile:
                row = line.split('\t')
                if len(row) > 1:
                    diseases = [d for d in row[disease_col].strip('"').split(' or ')
                                if len(d.strip()) > 0]
                    disease_pks = []
                    for disease in diseases:
                        disease_pk = add_or_get_disease(cur, disease)
                        disease_pks.append(disease_pk)

                    regulon = row[regulon_col].strip()
                    if len(regulon) > 0:
                        regulon_pk = add_or_get_regulon(cur, regulon)
                    else:
                        continue

                    pmids = [s.strip() for s in row[pmids_col].split('|') if len(s.strip()) > 0]
                    for pmid in pmids:
                        for disease_pk in disease_pks:
                            cur.execute('select count(*) from exc_disease_regulon where hr=%s and disease_id=%s and regulon_id=%s and pmid=%s',
                                        [hr, disease_pk, regulon_pk, pmid])
                            if cur.fetchone()[0] == 0:
                                cur.execute('insert into exc_disease_regulon (hr,disease_id,regulon_id,pmid) values (%s,%s,%s,%s)',
                                            [hr, disease_pk, regulon_pk, pmid])
            conn.commit()


def read_mutation_regulator(datadir, conn, hr):
    """The input data sets are inconsistent !!! The Regulator column can either
    be anywhere or not there at all !"""
    with open(os.path.join(datadir, 'Mut_Regulators_Raw_HR_%d.txt' % hr)) as infile:
        header = infile.readline().strip().split('\t')
        mutation_col = header.index('Mutation')
        pmids_col = header.index('PMIDs')
        try:
            regulator_col = header.index('Regulator')
        except:
            # there is a multivalued Regulators column
            print("WARNING: no 'Regulator' column found in 'Mut_Regulators_Raw_HR_%d.txt'" % hr)
            regulator_col = header.index('Regulators')

        with conn.cursor() as cur:
            for line in infile:
                row = line.strip().split('\t')
                if len(row) > 1:
                    mutation = row[mutation_col].strip()
                    if len(mutation.strip()) > 0:
                        mutation_pk = add_or_get_mutation(cur, mutation)
                    else:
                        continue

                    regulators = [r.strip() for r in row[regulator_col].strip().split(' or ')
                                  if len(r.strip()) > 0]

                    regulator_pks = []
                    if len(regulators) > 0:
                        for regulator in regulators:
                            regulator_pk = add_or_get_regulator(cur, regulator)
                            regulator_pks.append(regulator_pk)
                    else:
                        continue

                    pmids = [s.strip() for s in row[pmids_col].split('|') if len(s.strip()) > 0]
                    for pmid in pmids:
                        for regulator_pk in regulator_pks:
                            cur.execute('select count(*) from exc_mutation_regulator where hr=%s and mutation_id=%s and regulator_id=%s and pmid=%s',
                                        [hr, mutation_pk, regulator_pk, pmid])
                            if cur.fetchone()[0] == 0:
                                cur.execute('insert into exc_mutation_regulator (hr,mutation_id,regulator_id,pmid) values (%s,%s,%s,%s)',
                                            [hr, mutation_pk, regulator_pk, pmid])
            conn.commit()


def read_mutation_drugs(datadir, conn, hr):
    with open(os.path.join(datadir, 'Mut_Drugs_Raw_HR_%d.txt' % hr)) as infile:
        header = infile.readline().strip().split('\t')
        mutation_col = header.index('Mutation')
        drugs_col = header.index('Drugs')
        pmids_col = header.index('PMIDs')
        with conn.cursor() as cur:
            for line in infile:
                row = line.split('\t')
                if len(row) > 1:
                    mutation = row[mutation_col].strip()
                    if len(mutation) > 0:
                        mutation_pk = add_or_get_mutation(cur, mutation)
                    else:
                        continue

                    drug_pks = []
                    drugs = [d for d in row[drugs_col].strip('"').split(' or ') if len(d.strip()) > 0]
                    for drug in drugs:
                        drug_pk = add_or_get_drug(cur, drug)
                        drug_pks.append(drug_pk)

                    pmids = [s.strip() for s in row[pmids_col].split('|') if len(s.strip()) > 0]
                    for pmid in pmids:
                        for drug_pk in drug_pks:
                            cur.execute('select count(*) from exc_mutation_drug where hr=%s and mutation_id=%s and drug_id=%s and pmid=%s',
                                        [hr, mutation_pk, drug_pk, pmid])
                            if cur.fetchone()[0] == 0:
                                try:
                                    pmid = pmid.replace("']", "")  # hack
                                    cur.execute('insert into exc_mutation_drug (hr,mutation_id,drug_id,pmid) values (%s,%s,%s,%s)',
                                                [hr, mutation_pk, drug_pk, pmid])
                                except Exception as e:
                                    print('Exception in PMID: [%s]'% pmid)
                                    traceback.print_exc()
                                    raise

            conn.commit()


def read_all_datasets(datadir):
    print('Reading in all datasets')
    conn = dbconn()
    for hr in range(2, 7):
        print('Hazard ratio %d...' % hr)
        #read_cancer_mutation(datadir, conn, hr)
        #read_mm_mutation(datadir, conn, hr)
        #read_mm_regulator(datadir, conn, hr)
        #read_mm_regulon(datadir, conn, hr)
        #read_mutation_regulator(datadir, conn, hr)
        ## TODO: Check why there are no entries !!!!
        read_mutation_drugs(datadir, conn, hr)

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print('usage: import_excelra_datasets.py <datadir>')
    else:
        read_all_datasets(sys.argv[1])
