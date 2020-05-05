#!/usr/bin/env python

import logging
import json
import os
import traceback
from collections import defaultdict
import mysql.connector

# for mockup
import random

from flask import Flask, Response, url_for, redirect, render_template, request, session, flash, jsonify
import flask
from sqlalchemy import and_
# because we have an API, we need to allow cross-origin here
from flask_cors import CORS, cross_origin

from werkzeug import secure_filename
import requests
import xml.etree.ElementTree as ET
import re
import datasets

app = Flask(__name__)
CORS(app)

app.config.from_envvar('APP_SETTINGS')

MUTATION_TF_ROLES = { 1: 'down-regulates', 2: 'up-regulates'}
TF_BC_ROLES = { 1: 'activates', 2: 'represses' }


all_datasets = datasets.read_all_datasets(app.config['EXCELRA_DATASETS'])


def dbconn():
    return mysql.connector.connect(user=app.config['DATABASE_USER'],
                                   password=app.config['DATABASE_PASSWORD'],
                                   database=app.config['DATABASE_NAME'])

@app.route('/api/v1.0.0/bicluster/<cluster_id>')
def bicluster_info(cluster_id):
    """all genes in the specified bicluster"""
    conn = dbconn()
    cursor = conn.cursor(buffered=True)
    try:
        cursor.execute('select cox_hazard_ratio,trans_program from biclusters where name=%s', [cluster_id])
        hazard_ratio, trans_program = cursor.fetchone()

        # mutation role -> transcription factors
        """
        cursor.execute('select m.name,tfs.name,g.preferred,role from biclusters bc join bc_mutation_tf bmt on bc.id=bmt.bicluster_id join mutations m on m.id=bmt.mutation_id join tfs on tfs.id=bmt.tf_id left join genes g on tfs.name=g.ensembl_id where bc.name=%s', [cluster_id])
        mutations_tfs = [{"mutation": mutation,
                          "tf": tf, "tf_preferred": tf_preferred if tf_preferred is not None else tf,
                          "role": MUTATION_TF_ROLES[role]}
                         for mutation, tf, tf_preferred, role in cursor.fetchall()]"""
        cursor.execute('select mut.name as mutation,tfs.name as regulator,g.preferred as reg_preferred,bmt.role as mutation_role,bc_tf.role as regulator_role,bc.cox_hazard_ratio from bc_mutation_tf bmt join biclusters bc on bmt.bicluster_id=bc.id join mutations mut on bmt.mutation_id=mut.id join tfs on bmt.tf_id=tfs.id join bc_tf on bc.id=bc_tf.bicluster_id and tfs.id=bc_tf.tf_id left join genes g on g.ensembl_id=tfs.name where bc.name=%s',
                       [cluster_id])
        mutations_tfs = [{"mutation": mutation,
                          "tf": tf, "tf_preferred": tf_preferred if tf_preferred is not None else tf,
                          "mutation_role": MUTATION_TF_ROLES[mut_role],
                          'regulator_role': TF_BC_ROLES[bc_role],
                          'hazard_ratio': hazard_ratio}
                         for mutation, tf, tf_preferred, mut_role, bc_role, hazard_ratio in cursor.fetchall()]
        # transcription factor -> bicluster
        cursor.execute('select tfs.name,g.preferred,role,tfs.cox_hazard_ratio from biclusters bc join bc_tf bt on bc.id=bt.bicluster_id join tfs on tfs.id=bt.tf_id left join genes g on tfs.name=g.ensembl_id  where bc.name=%s', [cluster_id])
        tfs_bc = [{"tf": tf, "tf_preferred": tf_preferred if tf_preferred is not None else tf, "role": TF_BC_ROLES[role], "hazard_ratio": tf_hazard_ratio}
                  for tf, tf_preferred, role, tf_hazard_ratio in cursor.fetchall()]

        # bicluster genes
        cursor.execute('select g.preferred from biclusters bc join bicluster_genes bcg on bc.id=bcg.bicluster_id join genes g on g.id=bcg.gene_id where bc.name=%s', [cluster_id])
        genes = [row[0] for row in cursor.fetchall()]

        # regulon drugs
        cursor.execute('select d.name from biclusters bc join regulon_drugs rd on rd.regulon_id=bc.id join drugs d on rd.drug_id=d.id where bc.name=%s',
                       [cluster_id])
        drugs = [row[0] for row in cursor.fetchall()]

        # mechanism of action
        cursor.execute('select m.name from biclusters bc join regulon_mechanism_of_action rm on rm.regulon_id=bc.id join mechanisms_of_action m on rm.mechanism_of_action_id=m.id where bc.name=%s',
                       [cluster_id])
        moas = [row[0] for row in cursor.fetchall()]

        # hallmarks
        cursor.execute('select h.name from biclusters bc join regulon_hallmarks rh on rh.regulon_id=bc.id join hallmarks h on rh.hallmark_id=h.id where bc.name=%s',
                       [cluster_id])
        hallmarks = [row[0] for row in cursor.fetchall()]

        # target classes
        cursor.execute('select tc.name,rtc.pval from biclusters bc join regulon_target_class rtc on rtc.regulon_id=bc.id join target_classes tc on rtc.target_class_id=tc.id where bc.name=%s',
                       [cluster_id])
        target_classes = [{'name': name, 'pval': pval} for name, pval in cursor.fetchall()]

        return jsonify(cluster=cluster_id, trans_program=trans_program,
                       hazard_ratio=hazard_ratio,
                       mutations_tfs=mutations_tfs,
                       tfs_bc=tfs_bc,
                       genes=genes,
                       drugs=drugs,
                       mechanism_of_action=moas,
                       hallmarks=hallmarks,
                       target_classes=target_classes,
                       num_causal_flows=len(mutations_tfs))
    except:
        traceback.print_exc()
    finally:
        cursor.close()
        conn.close()


@app.route('/api/v1.0.0/search/<term>')
def simple_search(term):
    """retrieves the type of the term or empty result if not found"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute('select count(*) from tfs where tfs.name=%s', [term])
        num_results = cursor.fetchone()[0]
        if num_results > 0:
            return jsonify(found="yes", data_type="regulator")
        cursor.execute('select count(*) from mutations where mutations.name=%s', [term])
        num_results = cursor.fetchone()[0]
        if num_results > 0:
            return jsonify(found="yes", data_type="mutation")
        cursor.execute('select count(*) from biclusters bc where bc.name=%s', [term])
        num_results = cursor.fetchone()[0]
        if num_results > 0:
            return jsonify(found="yes", data_type="bicluster")
        cursor.execute('select count(*) from genes where ensembl_id=%s or preferred=%s',
                       [term, term])
        num_results = cursor.fetchone()[0]
        if num_results > 0:
            return jsonify(found="yes", data_type="gene")

        # transcription factor -> bicluster
        return jsonify(found="no", data_type="NA")
    finally:
        cursor.close()
        conn.close()


@app.route('/api/v1.0.0/cfsearch/<term>')
def causal_flow_search(term):
    """retrieves all causal flows related to the search term"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        # search causal flows by mutation or regulator
        cursor.execute("""select bc.name,mut.name,tfs.name,g.preferred,bmt.role,bc_tf.role,bc.cox_hazard_ratio,bgg.num_genes,bc.trans_program from bc_mutation_tf bmt join biclusters bc on bmt.bicluster_id=bc.id join mutations mut on bmt.mutation_id=mut.id join tfs on bmt.tf_id=tfs.id join bc_tf on bc.id=bc_tf.bicluster_id and tfs.id=bc_tf.tf_id join (select bc.id,count(bg.gene_id) as num_genes from biclusters bc join bicluster_genes bg on bc.id=bg.bicluster_id group by bc.id) as bgg on bc.id=bgg.id left join genes g on g.ensembl_id=tfs.name where g.preferred=%s or mut.name=%s""", [term, term])
        by_mutation = []
        by_regulator = []
        for bc,mut,tf,tf_preferred,mut_role,tf_role,hratio,ngenes,trans_program in cursor.fetchall():
            entry = {
                'bicluster': bc,
                'mutation': mut,
                'regulator': tf,
                'regulator_preferred': tf_preferred if tf_preferred is not None else tf,
                'mutation_role': MUTATION_ROLES[mut_role],
                'regulator_role': REGULATOR_ROLES[tf_role],
                'hazard_ratio': hratio,
                'num_genes': ngenes,
                'trans_program': trans_program
            }
            if tf_preferred == term:
                by_regulator.append(entry)

            if mut == term:
                by_mutation.append(entry)

        # search causal flows by regulon genes
        cursor.execute("""select bc.name,mut.name,tfs.name,g.preferred,bmt.role,bc_tf.role,bc.cox_hazard_ratio,bgg.num_genes,bc.trans_program from bc_mutation_tf bmt join biclusters bc on bmt.bicluster_id=bc.id join mutations mut on bmt.mutation_id=mut.id join tfs on bmt.tf_id=tfs.id join bc_tf on bc.id=bc_tf.bicluster_id and tfs.id=bc_tf.tf_id join (select bc.id,count(bg.gene_id) as num_genes from biclusters bc join bicluster_genes bg on bc.id=bg.bicluster_id group by bc.id) as bgg on bc.id=bgg.id left join genes g on g.ensembl_id=tfs.name where bc.id in (select bicluster_id from bicluster_genes bg join genes g on bg.gene_id=g.id where g.preferred=%s)""", [term])
        by_reggenes = [{
            'bicluster': bc,
            'mutation': mut,
            'regulator': tf,
            'regulator_preferred': tf_preferred if tf_preferred is not None else tf,
            'mutation_role': MUTATION_ROLES[mut_role],
            'regulator_role': REGULATOR_ROLES[tf_role],
            'hazard_ratio': hratio,
            'num_genes': ngenes,
            'trans_program': trans_program
        } for bc,mut,tf,tf_preferred,mut_role,tf_role,hratio,ngenes,trans_program in cursor.fetchall()]

        if len(by_mutation) == 0 and len(by_regulator) == 0 and len(by_reggenes) == 0:
            return jsonify(found='no')

        return jsonify(found="yes",
            by_mutation=by_mutation,
            by_regulator=by_regulator,
            by_reggenes=by_reggenes)
    finally:
        cursor.close()
        conn.close()


@app.route('/api/v1.0.0/causal_flows_with_program/<program>')
def causal_flows_with_program(program):
    """retrieves all causal flows containing this program"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        # search causal flows by regulon genes
        cursor.execute("""select bc.name,mut.name,tfs.name,g.preferred,bmt.role,bc_tf.role,bc.cox_hazard_ratio,bgg.num_genes,bc.trans_program from bc_mutation_tf bmt join biclusters bc on bmt.bicluster_id=bc.id join mutations mut on bmt.mutation_id=mut.id join tfs on bmt.tf_id=tfs.id join bc_tf on bc.id=bc_tf.bicluster_id and tfs.id=bc_tf.tf_id join (select bc.id,count(bg.gene_id) as num_genes from biclusters bc join bicluster_genes bg on bc.id=bg.bicluster_id group by bc.id) as bgg on bc.id=bgg.id left join genes g on g.ensembl_id=tfs.name where bc.trans_program=%s""", [program])
        result = [{
            'bicluster': bc,
            'mutation': mut,
            'regulator': tf,
            'regulator_preferred': tf_preferred if tf_preferred is not None else tf,
            'mutation_role': MUTATION_ROLES[mut_role],
            'regulator_role': REGULATOR_ROLES[tf_role],
            'hazard_ratio': hratio,
            'num_genes': ngenes,
            'trans_program': trans_program
        } for bc,mut,tf,tf_preferred,mut_role,tf_role,hratio,ngenes, trans_program in cursor.fetchall()]

        return jsonify(entries=result)
    finally:
        cursor.close()
        conn.close()


@app.route('/api/v1.0.0/completions/<term>')
def completions(term):
    """this is a function to serve the jquery autocomplete box"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute('select name from tfs where name like %s', ["%s%%" % term])
        terms = [{"id": row[0], "label": row[0], "value": row[0]}
                 for row in cursor.fetchall()]
        cursor.execute('select name from mutations where name like %s',
                       ["%s%%" % term])
        mutations = [{"id": row[0], "label": row[0], "value": row[0]}
                     for row in cursor.fetchall()]
        cursor.execute('select preferred from genes where preferred like %s',
                       ["%s%%" % term])
        gene_preferred = [{"id": row[0], "label": row[0], "value": row[0]}
                          for row in cursor.fetchall()]
        cursor.execute('select ensembl_id from genes where ensembl_id like %s',
                       ["%s%%" % term])
        gene_ensembl = [{"id": row[0], "label": row[0], "value": row[0]}
                        for row in cursor.fetchall()]

        completions = terms + mutations + gene_preferred + gene_ensembl
        return jsonify(completions=completions)
    finally:
        cursor.close()
        conn.close()


@app.route('/api/v1.0.0/biclusters')
def biclusters():
    """all biclusters in the system"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute('select name, cox_hazard_ratio from biclusters')
        biclusters = [{'cluster_id': row[0], 'hazard_ratio': row[1]}
                      for row in cursor.fetchall()]
        return jsonify(biclusters=biclusters)
    finally:
        cursor.close()
        conn.close()


@app.route('/api/v1.0.0/mutation/<mutation_name>')
def mutation(mutation_name):
    """information for the specified mutation"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute('select tfs.name,g.preferred,bc.name,bmt.role,bc.cox_hazard_ratio from bc_mutation_tf bmt join biclusters bc on bmt.bicluster_id=bc.id join mutations m on m.id=bmt.mutation_id join tfs on tfs.id=bmt.tf_id left join genes g on tfs.name=g.ensembl_id where m.name=%s',
                       [mutation_name])
        result = [{"regulator": tf, "regulator_preferred": tf_preferred if tf_preferred is not None else tf, "bicluster": bc,
                   "role": MUTATION_TF_ROLES[role],
                   "bc_cox_hazard_ratio": bc_cox_hazard_ratio}
                  for tf, tf_preferred, bc, role,bc_cox_hazard_ratio in cursor.fetchall()]
        return jsonify(mutation=mutation_name, entries=result)
    finally:
        cursor.close()
        conn.close()


@app.route('/api/v1.0.0/regulator/<tf_name>')
def regulator(tf_name):
    """information for the specified mutation"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute('select cox_hazard_ratio from tfs where name=%s', [tf_name])
        hazard_ratio = cursor.fetchone()[0]
        cursor.execute('select bc.name,bt.role,bc.cox_hazard_ratio,mut.name from bc_tf bt join biclusters bc on bt.bicluster_id=bc.id join tfs on bt.tf_id=tfs.id join bc_mutation_tf bmt on bmt.bicluster_id=bt.bicluster_id and bmt.tf_id=bt.tf_id join mutations mut on mut.id=bmt.mutation_id where tfs.name=%s',
                       [tf_name])
        result = [{"bicluster": bc, "role": TF_BC_ROLES[role],
                   "hazard_ratio": bc_hazard_ratio,
                   "mutation": mut}
                  for bc, role, bc_hazard_ratio, mut in cursor.fetchall()]
        cursor.execute('select preferred from genes where ensembl_id=%s', [tf_name])
        row = cursor.fetchone()
        if row is not None:
            reg_preferred = row[0]
        else:
            reg_preferred = tf_name
        return jsonify(regulator=tf_name, regulator_preferred=reg_preferred,
                       hazard_ratio=hazard_ratio, entries=result)
    finally:
        cursor.close()
        conn.close()

@app.route('/api/v1.0.0/summary')
def summary():
    """model summary"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute('select count(*) from biclusters')
        num_biclusters = cursor.fetchone()[0]
        cursor.execute('select count(*) from tfs')
        num_regulators = cursor.fetchone()[0]
        cursor.execute('select count(*) from mutations')
        num_mutations = cursor.fetchone()[0]
        cursor.execute('select count(*) from patients')
        num_patients = cursor.fetchone()[0]
        cursor.execute('select count(*) from bc_mutation_tf bmt join biclusters bc on bmt.bicluster_id=bc.id join mutations mut on bmt.mutation_id=mut.id join tfs on bmt.tf_id=tfs.id join bc_tf on bc.id=bc_tf.bicluster_id and tfs.id=bc_tf.tf_id left join genes g on g.ensembl_id=tfs.name')
        num_causal_flows = cursor.fetchone()[0]
        cursor.execute('select count(distinct trans_program) from biclusters')
        num_trans_programs = cursor.fetchone()[0]
        return jsonify(num_biclusters=num_biclusters,
                       num_regulators=num_regulators,
                       num_mutations=num_mutations,
                       num_patients=num_patients,
                       num_causal_flows=num_causal_flows,
                       num_trans_programs=num_trans_programs)
    finally:
        cursor.close()
        conn.close()


@app.route('/api/v1.0.0/gene_info/<gene_name>')
def gene_info(gene_name):
    """return all biclusters that contain this gene"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute('select entrez_id,ensembl_id,preferred from genes where entrez_id=%s or ensembl_id=%s or preferred=%s',
                       [gene_name, gene_name, gene_name])
        results = [(entrez_id, ensembl_id, preferred)
                   for entrez_id, ensembl_id, preferred in cursor.fetchall()]
        if len(results) > 0:
            ensembl_id = results[0][1]
            r = requests.get('https://rest.ensembl.org/lookup/id/%s?content-type=application/json;expand=1' % ensembl_id)
            if r.status_code == 200:
                ensdata = r.json()
                print(ensdata['description'])
            r = requests.get('https://rest.ensembl.org/xrefs/id/' + ensembl_id +
                             '?content-type=application/json;external_db=uniprot%;all_levels=1')
            if r.status_code == 200:
                xrefdata = r.json()
                uniprot_ids = [entry['primary_id'] for entry in xrefdata]
                if len(uniprot_ids) > 0:
                    uniprot_id = uniprot_ids[0]
                else:
                    uniprot_id = '-'
            if uniprot_id != '-':
                # retrieve the function information from UniProt
                r = requests.get('https://www.uniprot.org/uniprot/%s.xml' % uniprot_id)
                doc = ET.fromstring(r.text)
                function = '-'
                for child in doc:
                    localname = re.sub(r'{.*}', '', child.tag)
                    if localname == 'entry':
                        for c in child:
                            localname = re.sub(r'{.*}', '', c.tag)
                            if localname == 'comment' and c.attrib['type'] == 'function':
                                function = ""
                                for node in c:
                                    function += node.text
            return jsonify(entrez_id=results[0][0], ensembl_id=ensembl_id,
                           preferred=results[0][2], description=ensdata['description'],
                           uniprot_id=uniprot_id, function=function)
        else:
            return jsonify(entrez_id="NA", ensembl_id="NA", preferred="NA", description='NA',
                           uniprot_id='NA', function='NA')
    finally:
        cursor.close()
        conn.close()


@app.route('/api/v1.0.0/biclusters_for_gene/<gene_name>')
def biclusters_for_gene(gene_name):
    """return all biclusters that contain this gene"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute('select bc.name,cox_hazard_ratio from bicluster_genes bg join genes g on bg.gene_id=g.id join biclusters bc on bc.id=bg.bicluster_id where g.entrez_id=%s or g.ensembl_id=%s or g.preferred=%s',
                       [gene_name, gene_name, gene_name])
        return jsonify(biclusters=[{'cluster_id': row[0],
                                    'hazard_ratio': row[1]} for row in cursor.fetchall()])
    finally:
        cursor.close()
        conn.close()


@app.route('/api/v1.0.0/bicluster_network/<cluster_id>')
def bicluster_network(cluster_id):
    """all genes in the specified bicluster"""
    conn = dbconn()
    cursor = conn.cursor()
    elements = [{"data": {"id": cluster_id}, "classes": "bicluster"}]
    try:
        edge_count = 0
        # retrieve the mutation and transcription factor nodes/edges
        # transcription factor -> bicluster
        cursor.execute('select tfs.name,g.preferred,role from biclusters bc join bc_tf bt on bc.id=bt.bicluster_id join tfs on tfs.id=bt.tf_id left join genes g on tfs.name=g.ensembl_id where bc.name=%s', [cluster_id])
        for tf, tf_preferred, role in cursor.fetchall():
            tf = tf_preferred if tf_preferred is not None else tf
            elements.append({"data": {"id": tf}, "classes": "tf"})
            elements.append({"data": {"id": str(edge_count), "source": tf, "target": cluster_id},
                                      "classes": TF_BC_ROLES[role]})
            edge_count += 1

        # mutation role -> transcription factors
        cursor.execute('select m.name,tfs.name,g.preferred,role from biclusters bc join bc_mutation_tf bmt on bc.id=bmt.bicluster_id join mutations m on m.id=bmt.mutation_id join tfs on tfs.id=bmt.tf_id left join genes g on tfs.name=g.ensembl_id where bc.name=%s', [cluster_id])
        for mutation, tf, tf_preferred, role in cursor.fetchall():
            tf = tf_preferred if tf_preferred is not None else tf
            elements.append({"data": {"id": mutation}, "classes": "mutation"})
            elements.append({"data": {"id": str(edge_count), "source": mutation, "target": tf},
                                      "classes": MUTATION_TF_ROLES[role].replace('-', '_') })
            edge_count += 1


        return jsonify(elements=elements)
    finally:
        cursor.close()
        conn.close()


@app.route('/api/v1.0.0/bicluster_expressions/<cluster_id>')
def bicluster_expression_data(cluster_id):
    """returns data plot data in Highcharts format for bicluster expressions"""
    """this is actually box plot data, so the series needs to be a list of
    six tuples [condition_id, min, lower quartile, mean, upper quartile, max]]
    """
    conn = dbconn()
    cursor = conn.cursor()
    data = []
    try:
        cursor.execute('select p.name,median,min_value,max_value,lower_quartile,upper_quartile from bicluster_boxplot_data bbd join patients p on bbd.patient_id=p.id join biclusters bc on bbd.bicluster_id=bc.id where bc.name=%s order by median',
                       [cluster_id])
        for patient, median, minval, maxval, lower_quart, upper_quart in cursor.fetchall():
            data.append([patient, minval, lower_quart, median, upper_quart, maxval])
        return jsonify(data=data)
    finally:
        cursor.close()
        conn.close()


@app.route('/api/v1.0.0/bicluster_enrichment/<cluster_id>')
def bicluster_enrichment(cluster_id):
    """returns barplot enrichment data for tumor subtypes in quitiles
    for the given bicluster"""
    conn = dbconn()
    cursor = conn.cursor()
    subtypes = ['g_cimp', 'proneural', 'neural', 'classical', 'mesenchymal', 'control']
    # series is gene -> list of values
    series  = defaultdict(list)
    # mockup some data for now (3 conditions)
    for s in subtypes:
        series[s] = [random.uniform(-10.0, 10.0) for i in range(5)]
    conds = ['All', 'All', 'All', 'All', 'All']
    return jsonify(expressions=series, conditions=conds)


@app.route('/api/v1.0.0/patient/<patient_id>')
def patient_info(patient_id):
    """patient information"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute('select pfs_survival, pfs_status, age, sex from patients where name=%s',
                       [patient_id])
        surv, status, age, sex = cursor.fetchone()
        cursor.execute('select t.name,pt.tf_activity from patient_tf pt join patients p on pt.patient_id=p.id join tfs t on pt.tf_id=t.id where p.name=%s', [patient_id])
        tf_activity = [{'tf': tf, 'activity': activity} for tf, activity in cursor.fetchall()]
        return jsonify(patient=patient_id, pfs_survival=surv, pfs_status=status,
                       age=age, sex=sex,
                       tf_activity=tf_activity)
    finally:
        cursor.close()
        conn.close()


MUTATION_ROLES = {1: 'down-regulates', 2: 'up-regulates'}
REGULATOR_ROLES = {1: 'activates', 2: 'represses'}

# TODO: This should actually be simply a pregenerated JSON file
@app.route('/api/v1.0.0/causal_flow')
def causal_flow():
    """causal flow"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute("""select bc.name,mut.name,tfs.name,g.preferred,bmt.role,bc_tf.role,bc.cox_hazard_ratio,bgg.num_genes,bc.trans_program from bc_mutation_tf bmt join biclusters bc on bmt.bicluster_id=bc.id join mutations mut on bmt.mutation_id=mut.id join tfs on bmt.tf_id=tfs.id join bc_tf on bc.id=bc_tf.bicluster_id and tfs.id=bc_tf.tf_id join (select bc.id,count(bg.gene_id) as num_genes from biclusters bc join bicluster_genes bg on bc.id=bg.bicluster_id group by bc.id) as bgg on bc.id=bgg.id left join genes g on g.ensembl_id=tfs.name""")
        return jsonify(entries=[{
            'bicluster': bc,
            'mutation': mut,
            'regulator': tf,
            'regulator_preferred': tf_preferred if tf_preferred is not None else tf,
            'mutation_role': MUTATION_ROLES[mut_role],
            'regulator_role': REGULATOR_ROLES[tf_role],
            'hazard_ratio': hratio,
            'num_genes': ngenes,
            'trans_program': trans_program
        } for bc,mut,tf,tf_preferred,mut_role,tf_role,hratio,ngenes,trans_program in cursor.fetchall()])
    finally:
        cursor.close()
        conn.close()


@app.route('/api/v1.0.0/bicluster_patients/<cluster_id>')
def bicluster_patients(cluster_id):
    """returns the information for all the patients in the specified bicluster"""
    conn = dbconn()
    cursor = conn.cursor()
    data = []
    try:
        cursor.execute('select p.name,p.pfs_survival,p.pfs_status,p.sex,p.age from bicluster_patients bp join patients p on bp.patient_id=p.id join biclusters bc on bp.bicluster_id=bc.id where bc.name=%s', [cluster_id])
        for patient, survival, survival_status, sex, age in cursor.fetchall():
            data.append({'name': patient, 'survival': survival,
                         'survival_status': survival_status, 'sex': sex, 'age': age})
        return jsonify(data=data)
    finally:
        cursor.close()
        conn.close()


@app.route('/api/v1.0.0/bicluster_patient_survival/<cluster_id>')
def bicluster_patient_survival(cluster_id):
    """returns the information for all the patients in the specified bicluster"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute('select p.pfs_survival from bicluster_patients bp join patients p on bp.patient_id=p.id join biclusters bc on bp.bicluster_id=bc.id where bc.name=%s and p.pfs_survival is not null', [cluster_id])
        values = [row[0] for row in cursor.fetchall()]
        return jsonify(data=values)
    finally:
        cursor.close()
        conn.close()


@app.route('/api/v1.0.0/bicluster_patient_ages/<cluster_id>')
def bicluster_patient_ages(cluster_id):
    """returns the information for all the patients in the specified bicluster"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute('select p.age from bicluster_patients bp join patients p on bp.patient_id=p.id join biclusters bc on bp.bicluster_id=bc.id where bc.name=%s and p.age is not null', [cluster_id])
        values = [int(row[0]) for row in cursor.fetchall()]
        return jsonify(data=values)
    finally:
        cursor.close()
        conn.close()


@app.route('/api/v1.0.0/bicluster_patient_status/<cluster_id>')
def bicluster_patient_status(cluster_id):
    """returns the information for all the patients in the specified bicluster"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute('select p.age,count(p.age) from bicluster_patients bp join patients p on bp.patient_id=p.id join biclusters bc on bp.bicluster_id=bc.id where bc.name=%s and p.age is not null group by p.age;', [cluster_id])
        ages = [(age, count) for age, count in cursor.fetchall()]
        total = sum([t[1] for t in ages])
        buckets = defaultdict(int)
        age_values = []
        keys = ['<= 45', '46-65', '> 65']
        for age, count in ages:
            if age <= 45:
                buckets['<= 45'] += count
            elif age <= 65:
                buckets['46-65'] += count
            else:
                buckets['> 65'] += count
        for key, count in buckets.items():
            age_values.append({'name': key, 'y': float(count) / float(total)})

        cursor.execute('select p.sex,count(p.sex) from bicluster_patients bp join patients p on bp.patient_id=p.id join biclusters bc on bp.bicluster_id=bc.id where bc.name=%s and p.sex is not null group by p.sex;', [cluster_id])
        sex_map = {sex: count for sex, count in cursor.fetchall()}
        total = sum(sex_map.values())
        sexvals = []
        if total > 0:
            if 1 in sex_map:
                sexvals.append({'name': 'male', 'y': float(sex_map[1]) / float(total)})
            if 2 in sex_map:
                sexvals.append({'name': 'female', 'y': float(sex_map[2]) / float(total)})

        cursor.execute('select p.pfs_survival,count(p.pfs_survival) from bicluster_patients bp join patients p on bp.patient_id=p.id join biclusters bc on bp.bicluster_id=bc.id where bc.name=%s and p.pfs_survival is not null group by p.pfs_survival;', [cluster_id])
        surv_buckets = defaultdict(int)
        surv_values = []
        survs = [(survival, count) for survival, count in cursor.fetchall()]
        total = sum([t[1] for t in survs])
        keys = ['< 500', '500-1000', '> 1000']
        for surv, count in survs:
            if surv < 500:
                surv_buckets['< 500'] += count
            elif surv <= 1000:
                surv_buckets['500-1000'] += count
            else:
                surv_buckets['> 1000'] += count
        for key, count in surv_buckets.items():
            surv_values.append({'name': key, 'y': float(count) / float(total)})

        return jsonify(data={'age': age_values, 'sex': sexvals, 'survival': surv_values})

    finally:
        cursor.close()
        conn.close()


"""
EXCELRA ACCESS FUNCTIONS
"""

@app.route('/cancers')
def cancers():
    return jsonify(status='ok', cancers=all_datasets['cancers'])


@app.route('/diseases')
def diseases():
    return jsonify(status='ok', diseases=all_datasets['diseases'])


@app.route('/mutations')
def mutations():
    return jsonify(status='ok', mutations=all_datasets['mutations'])


@app.route('/regulators')
def regulators():
    return jsonify(status='ok', regulators=all_datasets['regulators'])


@app.route('/regulons')
def regulons():
    return jsonify(status='ok', regulons=all_datasets['regulons'])

@app.route('/drugs')
def drugs():
    return jsonify(status='ok', drugs=all_datasets['drugs'])

@app.route('/pmid_counts/<hr>/<cancer>/<disease>/<mutation>/<regulator>/<regulon>/<drug>')
def pmid_counts(hr, cancer, disease, mutation, regulator, regulon, drug):
    return jsonify(status='ok',
                   num_cancer_mutation_pmids=len(datasets.cancer_mutation(all_datasets, hr, cancer, mutation)),
                   num_disease_mutation_pmids=len(datasets.disease_mutation(all_datasets, hr, disease, mutation)),
                   num_disease_regulator_pmids=len(datasets.disease_regulator(all_datasets, hr, disease, regulator)),
                   num_disease_regulon_pmids=len(datasets.disease_regulon(all_datasets, hr, disease, regulon)),
                   num_mutation_regulator_pmids=len(datasets.mutation_regulator(all_datasets, hr, mutation, regulator)),
                   num_mutation_drug_pmids=len(datasets.mutation_drug(all_datasets, hr, mutation, drug)))


@app.route('/search_pmid_counts/<hr>', methods=['POST'])
def search_pmid_counts(hr):
    reqdata = request.get_json()
    search_term = reqdata['search']
    print(search_term)
    return jsonify(status="ok",
                   num_cancer_mutation_pmids=len(datasets.search_cancer_mutation(all_datasets, hr, search_term)),
                   num_disease_mutation_pmids=len(datasets.search_disease_mutation(all_datasets, hr, search_term)),
                   num_disease_regulator_pmids=len(datasets.search_disease_regulator(all_datasets, hr, search_term)),
                   num_disease_regulon_pmids=len(datasets.search_disease_regulon(all_datasets, hr, search_term)),
                   num_mutation_regulator_pmids=len(datasets.search_mutation_regulator(all_datasets, hr, search_term)),
                   num_mutation_drug_pmids=len(datasets.search_mutation_drug(all_datasets, hr, search_term)))


def fetch_articles(pmids):
    pmid_str = ','.join(pmids)
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=%s&retmode=xml" % pmid_str

    r = requests.get(url)
    try:
        root = ET.fromstring(r.content)
        result = []
        for item in root.findall('./PubmedArticle/MedlineCitation'):
            pmid = item.find('PMID').text
            article = item.find('Article')
            title = article.find('ArticleTitle').text
            try:
                abstract = article.find('Abstract/AbstractText').text
            except:
                abstract = ''
            pubdate_str = ''
            try:
                for pubdate in article.findall('.//PubDate'):
                    pubyear = pubdate.find('Year').text
                    pubmonth = pubdate.find('Month').text
                    pubdate_str = '%s/%s' % (pubyear, pubmonth)
            except:
                pass

            result.append({'pmid': pmid, 'title': title, 'abstract': abstract,
                           'pubdate': pubdate_str})

        return result
    except:
        raise


def batch_results(payload, pmids):
    items_per_page = payload['itemsPerPage']
    page = payload['page']

    result = []
    if len(pmids) > 0:
        total = len(pmids)
        if len(pmids) > items_per_page:
            offset = (page - 1) * items_per_page
            pmids = pmids[offset:offset + items_per_page]
        return fetch_articles(pmids)
    return result


def markup_assocs(content, assocs):
    for assoc in assocs:
        regexp = re.compile(re.escape(assoc), re.IGNORECASE)
        content = regexp.sub('<span class="marked">' + assoc + '</span>', content)
    return content

@app.route('/cancer_mutation_docs/<hr>/<cancer>/<mutation>', methods=['POST'])
def cancer_mutation_docs(hr, cancer, mutation):
    pmids = datasets.cancer_mutation(all_datasets, hr, cancer, mutation)
    total = len(pmids)
    result = batch_results(request.get_json(), pmids)
    for entry in result:
        #assocs = all_datasets[int(hr)]['pmid_mutation'][entry['pmid']]
        assocs = []
        if cancer != 'All':
            assocs.append(cancer)
        if mutation != 'All':
            assocs.append(mutation)

        entry['abstract'] = markup_assocs(entry['abstract'], assocs)
        entry['assocs'] = '->'.join(list(assocs))
    return jsonify(status='ok', total=total, data=result)


def __make_search_results(pmids, search_term):
    total = len(pmids)
    result = batch_results(request.get_json(), pmids)
    assocs = [search_term]
    for entry in result:
        entry['assocs'] = search_term
        entry['abstract'] = markup_assocs(entry['abstract'], assocs)
    return result, total


@app.route('/cancer_mutation_search/<hr>/<search_term>', methods=['POST'])
def cancer_mutation_search(hr, search_term):
    pmids = datasets.search_cancer_mutation(all_datasets, hr, search_term)
    result, total = __make_search_results(pmids, search_term)
    return jsonify(status='ok', total=total, data=result)


@app.route('/disease_mutation_docs/<hr>/<disease>/<mutation>', methods=['POST'])
def disease_mutation_docs(hr, disease, mutation):
    pmids = datasets.disease_mutation(all_datasets, hr, disease, mutation)
    total = len(pmids)
    result = batch_results(request.get_json(), pmids)
    for entry in result:
        #assocs = all_datasets[int(hr)]['pmid_mutation'][entry['pmid']]
        assocs = []
        if disease != 'All':
            assocs.append(disease)
        if mutation != 'All':
            assocs.append(mutation)

        entry['assocs'] = '->'.join(list(assocs))
        entry['abstract'] = markup_assocs(entry['abstract'], assocs)
    return jsonify(status='ok', total=total, data=result)


@app.route('/disease_mutation_search/<hr>/<search_term>', methods=['POST'])
def disease_mutation_search(hr, search_term):
    pmids = datasets.search_disease_mutation(all_datasets, hr, search_term)
    result, total = __make_search_results(pmids, search_term)
    return jsonify(status='ok', total=total, data=result)


@app.route('/disease_regulator_docs/<hr>/<disease>/<regulator>', methods=['POST'])
def disease_regulator_docs(hr, disease, regulator):
    pmids = datasets.disease_regulator(all_datasets, hr, disease, regulator)
    total = len(pmids)
    result = batch_results(request.get_json(), pmids)
    for entry in result:
        #assocs = all_datasets[int(hr)]['pmid_regulator'][entry['pmid']]
        assocs = []
        if disease != 'All':
            assocs.append(disease)
        if regulator != 'All':
            assocs.append(regulator)

        entry['assocs'] = '->'.join(list(assocs))
        try:
            entry['abstract'] = markup_assocs(entry['abstract'], assocs)
        except:
            traceback.print_exc()
    return jsonify(status='ok', total=total, data=result)


@app.route('/disease_regulator_search/<hr>/<search_term>', methods=['POST'])
def disease_regulator_search(hr, search_term):
    pmids = datasets.search_disease_regulator(all_datasets, hr, search_term)
    result, total = __make_search_results(pmids, search_term)
    return jsonify(status='ok', total=total, data=result)


@app.route('/disease_regulon_docs/<hr>/<disease>/<regulon>', methods=['POST'])
def disease_regulon_docs(hr, disease, regulon):
    pmids = datasets.disease_regulon(all_datasets, hr, disease, regulon)
    total = len(pmids)
    result = batch_results(request.get_json(), pmids)
    for entry in result:
        #assocs = all_datasets[int(hr)]['pmid_regulon'][entry['pmid']]
        assocs = []
        if disease != 'All':
            assocs.append(disease)
        if regulon != 'All':
            assocs.append(regulon)

        entry['assocs'] = '->'.join(list(assocs))
        entry['abstract'] = markup_assocs(entry['abstract'], assocs)
    return jsonify(status='ok', total=total, data=result)


@app.route('/disease_regulon_search/<hr>/<search_term>', methods=['POST'])
def disease_regulon_search(hr, search_term):
    pmids = datasets.search_disease_regulon(all_datasets, hr, search_term)
    result, total = __make_search_results(pmids, search_term)
    return jsonify(status='ok', total=total, data=result)


@app.route('/mutation_regulator_docs/<hr>/<mutation>/<regulator>', methods=['POST'])
def mutation_regulator_docs(hr, mutation, regulator):
    pmids = datasets.mutation_regulator(all_datasets, hr, mutation, regulator)
    total = len(pmids)
    result = batch_results(request.get_json(), pmids)
    for entry in result:
        #assocs = all_datasets[int(hr)]['pmid_regulator'][entry['pmid']]
        assocs = []
        if mutation != 'All':
            assocs.append(mutation)
        if regulator != 'All':
            assocs.append(regulator)

        entry['assocs'] = '->'.join(list(assocs))
        entry['abstract'] = markup_assocs(entry['abstract'], assocs)
    return jsonify(status='ok', total=total, data=result)


@app.route('/mutation_regulator_search/<hr>/<search_term>', methods=['POST'])
def mutation_regulator_search(hr, search_term):
    pmids = datasets.search_mutation_regulator(all_datasets, hr, search_term)
    result, total = __make_search_results(pmids, search_term)
    return jsonify(status='ok', total=total, data=result)


@app.route('/mutation_drug_docs/<hr>/<mutation>/<drug>', methods=['POST'])
def mutation_drug_docs(hr, mutation, drug):
    pmids = datasets.mutation_drug(all_datasets, hr, mutation, drug)
    total = len(pmids)
    result = batch_results(request.get_json(), pmids)
    for entry in result:
        assocs = []
        if mutation != 'All':
            #mutation_assocs = all_datasets[int(hr)]['pmid_mutation'][entry['pmid']]
            assocs.append(mutation)
        if drug != 'All':
            #drug_assocs = all_datasets[int(hr)]['pmid_drug'][entry['pmid']]
            assocs.append(drug)

        entry['assocs'] = '->'.join(list(assocs))
        entry['abstract'] = markup_assocs(entry['abstract'], assocs)
    return jsonify(status='ok', total=total, data=result)



@app.route('/mutation_drug_search/<hr>/<search_term>', methods=['POST'])
def mutation_drug_search(hr, search_term):
    pmids = datasets.search_mutation_drug(all_datasets, hr, search_term)
    result, total = __make_search_results(pmids, search_term)
    return jsonify(status='ok', total=total, data=result)


if __name__ == '__main__':
    app.debug = True
    app.secret_key = 'trstrestnorgp654g'
    app.run(host='0.0.0.0', debug=True)
