import os
from collections import defaultdict

def read_cancer_mutation(datadir, hr):
    """Return a PMID dataset that is indexed for both cancers and mutations"""
    result = {}
    result['cancer'] = defaultdict(set)
    result['mutation'] = defaultdict(set)

    with open(os.path.join(datadir, 'Mut_Cancer_HR_%d.txt' % hr)) as infile:
        header = infile.readline()
        for line in infile:
            row = line.strip().split('\t')
            if len(row) > 1:
                # synonyms for cancer
                diseases = [d for d in row[0].strip('"').split(' or ') if len(d.strip()) > 0]
                mutation = row[1].strip()
                pmids = [s.strip() for s in row[6].split('|') if len(s.strip()) > 0]
                for d in diseases:
                    result['cancer'][d].update(pmids)
                result['mutation'][mutation].update(pmids)
    return result


def read_mm_mutation(datadir, hr):
    """Also known as disease-mutation"""
    result = {}
    result['disease'] = defaultdict(set)
    result['mutation'] = defaultdict(set)
    all_diseases = set()

    with open(os.path.join(datadir, 'Mut_MM_HR_%d.txt' % hr)) as infile:
        header = infile.readline()
        for line in infile:
            row = line.split('\t')
            diseases = [d for d in row[0].strip('"').split(' or ') if len(d.strip()) > 0]
            all_diseases.update(diseases)
            mutation = row[1].strip()
            pmids = [s.strip() for s in row[6].split('|') if len(s.strip()) > 0]
            for d in diseases:
                result['disease'][d].update(pmids)
            result['mutation'][mutation].update(pmids)
    return result, all_diseases


def read_mm_regulator(datadir, hr):
    """also known as disease regulator"""
    all_regulators = set()
    result = {}
    result['disease'] = defaultdict(set)
    result['regulator'] = defaultdict(set)

    with open(os.path.join(datadir, 'Regulator_MM_HR_%d.txt' % hr)) as infile:
        header = infile.readline()
        for line in infile:
            row = line.split('\t')
            if len(row) > 1:
                diseases = [d for d in row[0].strip('"').split(' or ') if len(d.strip()) > 0]
                regulator = row[1].strip()
                regulator_ensg = row[2].strip()
                pmids = [s.strip() for s in row[7].split('|') if len(s.strip()) > 0]
                for d in diseases:
                    result['disease'][d].update(pmids)
                result['regulator'][regulator].update(pmids)

                if len(regulator.strip()) > 0:
                    all_regulators.add(regulator)

    return result, all_regulators


def read_mm_regulon(datadir, hr):
    """also known as disease regulon"""
    all_regulons = set()
    result = {}
    result['disease'] = defaultdict(set)
    result['regulon'] = defaultdict(set)

    with open(os.path.join(datadir, 'Regulon_MM_HR_%d.txt' % hr)) as infile:
        header = infile.readline()
        for line in infile:
            row = line.split('\t')
            if len(row) > 1:
                diseases = [d for d in row[0].strip('"').split(' or ') if len(d.strip()) > 0]
                regulon = row[1]
                regulon_ensg = row[2]
                pmids = [s.strip() for s in row[7].split('|') if len(s.strip()) > 0]
                for d in diseases:
                    result['disease'][d].update(pmids)
                result['regulon'][regulon].update(pmids)

                if len(regulon.strip()) > 0:
                    all_regulons.add(regulon)

    return result, all_regulons


def read_mutation_regulator(datadir, hr):
    result = {}
    result['mutation'] = defaultdict(set)
    result['regulator'] = defaultdict(set)

    with open(os.path.join(datadir, 'Mut_Regulators_Raw_HR_%d.txt' % hr)) as infile:
        header = infile.readline()
        for line in infile:
            row = line.strip().split('\t')
            if len(row) > 1:
                mutation = row[0]
                regulator = row[1]
                pmids = [s.strip() for s in row[9].split('|') if len(s.strip()) > 0]
                result['mutation'][mutation].update(pmids)
                result['regulator'][regulator].update(pmids)
    return result


def read_mutation_drugs(datadir, hr):
    result = {}
    result['mutation'] = defaultdict(set)
    result['drug'] = defaultdict(set)

    all_drugs = set()
    with open(os.path.join(datadir, 'Mut_Drugs_Raw_HR_%d.txt' % hr)) as infile:
        header = infile.readline()
        for line in infile:
            row = line.split('\t')
            if len(row) > 1:
                mutation = row[0]
                drugs = row[1].strip('"').split(' or ')
                pmids = [s.strip() for s in row[8].split('|') if len(s.strip()) > 0]
                result['mutation'][mutation].update(pmids)
                if len(drugs) > 0:
                    for drug in drugs:
                        if len(drug.strip()) > 0:
                            all_drugs.add(drug)
                            result['drug'][drug].update(pmids)
    return result, all_drugs


def read_all_datasets(datadir):
    print('Reading in all datasets')
    result = {}
    all_cancers = set()
    all_diseases = set()
    all_mutations = set()
    all_regulators = set()
    all_regulons = set()
    all_drugs = set()

    for hr in range(2, 7):
        result[hr] = {}
        print('Hazard ratio %d...' % hr)
        # add the cancer dataset for this hr
        ds = read_cancer_mutation(datadir, hr)
        cancers = ds['cancer'].keys()
        result[hr]['cancers'] = ds
        for c in cancers:
            all_cancers.add(c)

        ds, tmp_all_diseases = read_mm_mutation(datadir, hr)
        result[hr]['mm_mutation'] = ds
        all_diseases.update(tmp_all_diseases)

        ds, tmp_all_regulators = read_mm_regulator(datadir, hr)
        result[hr]['mm_regulator'] = ds
        all_regulators.update(tmp_all_regulators)

        ds, tmp_all_regulons = read_mm_regulon(datadir, hr)
        result[hr]['mm_regulon'] = ds
        all_regulons.update(tmp_all_regulons)

        ds = read_mutation_regulator(datadir, hr)
        result[hr]['mutation_regulator'] = ds
        mutations = ds['mutation'].keys()
        for m in mutations:
            all_mutations.add(m)

        ds, tmp_all_drugs = read_mutation_drugs(datadir, hr)
        result[hr]['mutation_drugs'] = ds
        all_drugs.update(tmp_all_drugs)

    result['cancers'] = sorted(all_cancers)
    result['diseases'] = sorted(all_diseases)
    result['mutations'] = sorted(all_mutations)
    result['regulators'] = sorted(all_regulators)
    result['regulons'] = sorted(all_regulons)
    result['drugs'] = sorted(all_drugs)
    return result


"""
Queries
"""

def cancer_mutation(all_datasets, hr, cancer, mutation):
    mm_mutation = all_datasets[int(hr)]['cancers']['mutation']
    mm_cancer = all_datasets[int(hr)]['cancers']['cancer']
    mutation_pmids = set()
    cancer_pmids = set()
    if mutation == 'All':
        for tmp_pmids in mm_mutation.values():
            mutation_pmids.update(tmp_pmids)
    else:
        mutation_pmids.update(mm_mutation[mutation])

    if cancer == 'All':
        for tmp_pmids in mm_cancer.values():
            cancer_pmids.update(tmp_pmids)
    else:
        cancer_pmids.update(mm_cancer[cancer])
    pmids = mutation_pmids.intersection(cancer_pmids)
    return list(pmids)


def disease_mutation(all_datasets, hr, disease, mutation):
    mm_mutation = all_datasets[int(hr)]['mm_mutation']['mutation']
    mm_disease = all_datasets[int(hr)]['mm_mutation']['disease']
    disease_pmids = set()
    mutation_pmids = set()
    if mutation == 'All':
        for tmp_pmids in mm_mutation.values():
            mutation_pmids.update(tmp_pmids)
    else:
        mutation_pmids.update(mm_mutation[mutation])

    if disease == 'All':
        for tmp_pmids in mm_disease.values():
            disease_pmids.update(tmp_pmids)
    else:
        disease_pmids.update(mm_disease[disease])

    pmids = mutation_pmids.intersection(disease_pmids)
    return list(pmids)


def disease_regulator(all_datasets, hr, disease, regulator):
    disease_regulator = all_datasets[int(hr)]['mm_regulator']['regulator']
    mm_disease = all_datasets[int(hr)]['mm_mutation']['disease']
    disease_pmids = set()
    regulator_pmids = set()

    if regulator == 'All':
        for tmp_pmids in disease_regulator.values():
            regulator_pmids.update(tmp_pmids)
    else:
        regulator_pmids.update(disease_regulator[regulator])

    if disease == 'All':
        for tmp_pmids in mm_disease.values():
            disease_pmids.update(tmp_pmids)
    else:
        disease_pmids.update(mm_disease[disease])

    pmids = regulator_pmids.intersection(disease_pmids)
    return list(pmids)


def disease_regulon(all_datasets, hr, disease, regulon):
    disease_regulon = all_datasets[int(hr)]['mm_regulon']['regulon']
    mm_disease = all_datasets[int(hr)]['mm_mutation']['disease']
    disease_pmids = set()
    regulon_pmids = set()

    if regulon == 'All':
        for tmp_pmids in disease_regulon.values():
            regulon_pmids.update(tmp_pmids)
    else:
        regulon_pmids.update(disease_regulon[regulon])

    if disease == 'All':
        for tmp_pmids in mm_disease.values():
            disease_pmids.update(tmp_pmids)
    else:
        disease_pmids.update(mm_disease[disease])

    pmids = regulon_pmids.intersection(disease_pmids)
    return list(pmids)


def mutation_regulator(all_datasets, hr, mutation, regulator):
    mutation_regulator_m = all_datasets[int(hr)]['mutation_regulator']['mutation']
    mutation_regulator_r = all_datasets[int(hr)]['mutation_regulator']['regulator']
    mutation_pmids = set()
    regulator_pmids = set()
    if mutation == 'All':
        for tmp_pmids in mutation_regulator_m.values():
            mutation_pmids.update(tmp_pmids)
    else:
        mutation_pmids.update(mutation_regulator_m[mutation])

    if regulator == 'All':
        for tmp_pmids in mutation_regulator_r.values():
            regulator_pmids.update(tmp_pmids)
    else:
        regulator_pmids.update(mutation_regulator_r[regulator])

    pmids = mutation_pmids.intersection(regulator_pmids)
    return list(pmids)


def mutation_drug(all_datasets, hr, mutation, drug):
    mutation_drug_m = all_datasets[int(hr)]['mutation_drugs']['mutation']
    mutation_drug_d = all_datasets[int(hr)]['mutation_drugs']['drug']
    mutation_pmids = set()
    drug_pmids = set()

    if mutation == 'All':
        for tmp_pmids in mutation_drug_m.values():
            mutation_pmids.update(tmp_pmids)
    else:
        mutation_pmids.update(mutation_drug_m[mutation])

    if drug == 'All':
        for tmp_pmids in mutation_drug_d.values():
            drug_pmids.update(tmp_pmids)
    else:
        drug_pmids.update(mutation_drug_d[drug])

    pmids = mutation_pmids.intersection(drug_pmids)
    return list(pmids)
