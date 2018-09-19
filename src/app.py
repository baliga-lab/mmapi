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

app = Flask(__name__)
CORS(app)

app.config.from_envvar('APP_SETTINGS')

MUTATION_TF_ROLES = { 1: 'down-regulates', 2: 'up-regulates'}
TF_BC_ROLES = { 1: 'activates', 2: 'represses' }

def dbconn():
    return mysql.connector.connect(user=app.config['DATABASE_USER'],
                                   password=app.config['DATABASE_PASSWORD'],
                                   database=app.config['DATABASE_NAME'])

@app.route('/api/v1.0.0/bicluster/<cluster_id>')
def bicluster_info(cluster_id):
    """all genes in the specified bicluster"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute('select cox_hazard_ratio from biclusters where name=%s', [cluster_id])
        hazard_ratio = cursor.fetchone()[0]

        # mutation role -> transcription factors
        cursor.execute('select m.name,tfs.name,role from biclusters bc join bc_mutation_tf bmt on bc.id=bmt.bicluster_id join mutations m on m.id=bmt.mutation_id join tfs on tfs.id=bmt.tf_id where bc.name=%s', [cluster_id])
        mutations_tfs = [{"mutation": mutation, "tf": tf, "role": MUTATION_TF_ROLES[role]}
                         for mutation, tf, role in cursor.fetchall()]

        # transcription factor -> bicluster
        cursor.execute('select tfs.name,role,tfs.cox_hazard_ratio from biclusters bc join bc_tf bt on bc.id=bt.bicluster_id join tfs on tfs.id=bt.tf_id where bc.name=%s', [cluster_id])
        tfs_bc = [{"tf": tf, "role": TF_BC_ROLES[role], "hazard_ratio": tf_hazard_ratio}
                  for tf, role, tf_hazard_ratio in cursor.fetchall()]

        # bicluster genes
        cursor.execute('select g.preferred from biclusters bc join bicluster_genes bcg on bc.id=bcg.bicluster_id join genes g on g.id=bcg.gene_id where bc.name=%s', [cluster_id])
        genes = [row[0] for row in cursor.fetchall()]

        # TODO: number of causal flows, which is the number of entries in the
        # causal flow table
        return jsonify(cluster=cluster_id,
                       hazard_ratio=hazard_ratio,
                       mutations_tfs=mutations_tfs,
                       tfs_bc=tfs_bc, genes=genes,
                       num_causal_flows=1)
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
        cursor.execute('select tfs.name,bc.name,bmt.role,bc.cox_hazard_ratio from bc_mutation_tf bmt join biclusters bc on bmt.bicluster_id=bc.id join mutations m on m.id=bmt.mutation_id join tfs on tfs.id=bmt.tf_id where m.name=%s',
                       [mutation_name])
        result = [{"regulator": tf, "bicluster": bc,
                   "role": MUTATION_TF_ROLES[role],
                   "bc_cox_hazard_ratio": bc_cox_hazard_ratio}
                  for tf, bc, role,bc_cox_hazard_ratio in cursor.fetchall()]
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
        cursor.execute('select bc.name,bt.role,bc.cox_hazard_ratio from bc_tf bt join biclusters bc on bt.bicluster_id=bc.id join tfs on bt.tf_id=tfs.id where tfs.name=%s',
                       [tf_name])
        result = [{"bicluster": bc, "role": TF_BC_ROLES[role], "hazard_ratio": bc_hazard_ratio}
                  for bc, role, bc_hazard_ratio in cursor.fetchall()]
        return jsonify(regulator=tf_name, hazard_ratio=hazard_ratio, entries=result)
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
        return jsonify(num_biclusters=num_biclusters,
                       num_regulators=num_regulators,
                       num_mutations=num_mutations,
                       num_patients=num_patients)
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
        cursor.execute('select tfs.name,role from biclusters bc join bc_tf bt on bc.id=bt.bicluster_id join tfs on tfs.id=bt.tf_id where bc.name=%s', [cluster_id])
        for tf, role in cursor.fetchall():
            elements.append({"data": {"id": tf}, "classes": "tf"})
            elements.append({"data": {"id": str(edge_count), "source": tf, "target": cluster_id, "role": TF_BC_ROLES[role]}})
            edge_count += 1

        # mutation role -> transcription factors
        cursor.execute('select m.name,tfs.name,role from biclusters bc join bc_mutation_tf bmt on bc.id=bmt.bicluster_id join mutations m on m.id=bmt.mutation_id join tfs on tfs.id=bmt.tf_id where bc.name=%s', [cluster_id])
        for mutation, tf, role in cursor.fetchall():
            elements.append({"data": {"id": mutation}, "classes": "mutation"})
            elements.append({"data": {"id": str(edge_count), "source": mutation, "target": tf, "role": MUTATION_TF_ROLES[role] }})
            edge_count += 1


        return jsonify(elements=elements)
    finally:
        cursor.close()
        conn.close()


#
# TODO APIs for bicluster details
#
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

@app.route('/api/v1.0.0/causal_flow')
def causal_flow():
    """causal flow"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute("""select bc.name,mut.name,tfs.name,bmt.role,bc_tf.role,bc.cox_hazard_ratio,bgg.num_genes from bc_mutation_tf bmt join biclusters bc on bmt.bicluster_id=bc.id join mutations mut on bmt.mutation_id=mut.id join tfs on bmt.tf_id=tfs.id join bc_tf on bc.id=bc_tf.bicluster_id and tfs.id=bc_tf.tf_id join (select bc.id,count(bg.gene_id) as num_genes from biclusters bc join bicluster_genes bg on bc.id=bg.bicluster_id group by bc.id) as bgg on bc.id=bgg.id""")
        return jsonify(entries=[{
            'bicluster': bc,
            'mutation': mut,
            'regulator': tf,
            'mutation_role': MUTATION_ROLES[mut_role],
            'regulator_role': REGULATOR_ROLES[tf_role],
            'hazard_ratio': hratio,
            'num_genes': ngenes
        } for bc,mut,tf,mut_role,tf_role,hratio,ngenes in cursor.fetchall()])
    finally:
        cursor.close()
        conn.close()


if __name__ == '__main__':
    app.debug = True
    app.secret_key = 'trstrestnorgp654g'
    app.run(host='0.0.0.0', debug=True)
