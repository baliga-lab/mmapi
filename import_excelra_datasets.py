#!/usr/bin/env python3
import os, sys
from collections import defaultdict
import MySQLdb

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
                    for cancer in cancers:
                        cancer_pk = add_or_get_cancer(cur, cancer)

                    if len(mutation) > 0:
                        mutation_pk = add_or_get_mutation(cur, mutation)
                    """
                    all_cancers.update(cancers)
                    all_mutations.add(mutation)
                    pmids = [s.strip() for s in row[pmids_col].split('|') if len(s.strip()) > 0]
                    for c in cancers:
                        result['cancer'][c].update(pmids)
                    for pmid in pmids:
                        hr_result['pmid_cancer'][pmid].update(cancers)
                        hr_result['pmid_mutation'][pmid].add(mutation)
                    result['mutation'][mutation].update(pmids)
                    """
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
                for disease in diseases:
                    disease_pk = add_or_get_disease(cur, disease)

                mutations = [m.strip() for m in row[mutation_col].strip().split(' or ')
                             if len(m.strip()) > 0]
                for mutation in mutations:
                    mutation_pk = add_or_get_mutation(cur, mutation)

                pmids = [s.strip() for s in row[pmids_col].split('|') if len(s.strip()) > 0]
                """
                for d in diseases:
                    result['disease'][d].update(pmids)
                for mutation in mutations:
                    result['mutation'][mutation].update(pmids)

                for pmid in pmids:
                    hr_result['pmid_disease'][pmid].update(diseases)
                    hr_result['pmid_mutation'][pmid].update(mutations)"""
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
                    for disease in diseases:
                        disease_pk = add_or_get_disease(cur, disease)
                    regulator = row[regulator_col].strip()
                    if len(regulator) > 0:
                        regulator_pk = add_or_get_regulator(cur, regulator)
                    """
                    if len(regulator.strip()) > 0:
                        all_regulators.add(regulator)"""

                    regulator_ensg = row[2].strip()
                    pmids = [s.strip() for s in row[pmids_col].split('|') if len(s.strip()) > 0]
                    """
                    for d in diseases:
                        result['disease'][d].update(pmids)
                    result['regulator'][regulator].update(pmids)

                    for pmid in pmids:
                        hr_result['pmid_disease'][pmid].update(diseases)
                        hr_result['pmid_regulator'][pmid].add(regulator)"""
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
                    for disease in diseases:
                        disease_pk = add_or_get_disease(cur, disease)
                    regulon = row[regulon_col].strip()
                    if len(regulon) > 0:
                        regulon_pk = add_or_get_regulon(cur, regulon)

                    regulon_ensg = row[2]
                    pmids = [s.strip() for s in row[pmids_col].split('|') if len(s.strip()) > 0]
                    """
                    for d in diseases:
                        result['disease'][d].update(pmids)
                    result['regulon'][regulon].update(pmids)

                    for pmid in pmids:
                        hr_result['pmid_disease'][pmid].update(diseases)
                        hr_result['pmid_regulon'][pmid].add(regulon)"""
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

                    if len(regulators) > 0:
                        for regulator in regulators:
                            regulator_pk = add_or_get_regulator(cur, regulator)
                    else:
                        continue

                    pmids = [s.strip() for s in row[pmids_col].split('|') if len(s.strip()) > 0]

                    """
                    result['mutation'][mutation].update(pmids)
                    for regulator in regulators:
                        result['regulator'][regulator].update(pmids)

                    for pmid in pmids:
                        hr_result['pmid_regulator'][pmid].update(regulators)
                        hr_result['pmid_mutation'][pmid].add(mutation)"""
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

                    drugs = [d for d in row[drugs_col].strip('"').split(' or ') if len(d.strip()) > 0]
                    for drug in drugs:
                        drug_pk = add_or_get_drug(cur, drug)

                    pmids = [s.strip() for s in row[pmids_col].split('|') if len(s.strip()) > 0]
                    """
                    if len(mutation.strip()) > 0:
                        all_mutations.add(mutation)
                        result['mutation'][mutation].update(pmids)
                    for drug in drugs:
                        result['drug'][drug].update(pmids)

                    for pmid in pmids:
                        hr_result['pmid_drug'][pmid].update(drugs)
                        hr_result['pmid_mutation'][pmid].add(mutation)"""
            conn.commit()


def read_all_datasets(datadir):
    print('Reading in all datasets')
    conn = dbconn()
    for hr in range(2, 7):
        print('Hazard ratio %d...' % hr)
        read_cancer_mutation(datadir, conn, hr)
        read_mm_mutation(datadir, conn, hr)
        read_mm_regulator(datadir, conn, hr)
        read_mm_regulon(datadir, conn, hr)
        read_mutation_regulator(datadir, conn, hr)
        read_mutation_drugs(datadir, conn, hr)

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print('usage: import_excelra_datasets.py <datadir>')
    else:
        read_all_datasets(sys.argv[1])
