
#!/usr/bin/env python

import logging
import json
import os
import traceback
import mysql.connector

from flask import Flask, Response, url_for, redirect, render_template, request, session, flash, jsonify
import flask
from sqlalchemy import and_
# because we have an API, we need to allow cross-origin here
from flask_cors import CORS, cross_origin

from werkzeug import secure_filename

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
        # mutation role -> transcription factors
        cursor.execute('select m.name,tfs.name,role from biclusters bc join bc_mutation_tf bmt on bc.id=bmt.bicluster_id join mutations m on m.id=bmt.mutation_id join tfs on tfs.id=bmt.tf_id where bc.name=%s', [cluster_id])
        mutations_tfs = [{"mutation": mutation, "tf": tf, "role": MUTATION_TF_ROLES[role]}
                         for mutation, tf, role in cursor.fetchall()]

        # transcription factor -> bicluster
        cursor.execute('select tfs.name,role from biclusters bc join bc_tf bt on bc.id=bt.bicluster_id join tfs on tfs.id=bt.tf_id where bc.name=%s', [cluster_id])
        tfs_bc = [{"tf": tf, "role": TF_BC_ROLES[role]}
                  for tf, role in cursor.fetchall()]

        # bicluster genes
        cursor.execute('select g.preferred from biclusters bc join bicluster_genes bcg on bc.id=bcg.bicluster_id join genes g on g.id=bcg.gene_id where bc.name=%s', [cluster_id])
        genes = [row[0] for row in cursor.fetchall()]
        return jsonify(cluster=cluster_id, mutations_tfs=mutations_tfs,
                       tfs_bc=tfs_bc, genes=genes)
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
        cursor.execute('select name from tfs where tfs.name like %s', ["%s%%" % term])
        terms = [{"id": row[0], "label": row[0], "value": row[0]}
                 for row in cursor.fetchall()]
        cursor.execute('select name from mutations where mutations.name like %s',
                       ["%s%%" % term])
        mutations = [{"id": row[0], "label": row[0], "value": row[0]}
                     for row in cursor.fetchall()]
        completions = terms + mutations
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
        cursor.execute('select name from biclusters')
        biclusters = [row[0] for row in cursor.fetchall()]
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
        cursor.execute('select tfs.name,bc.name,bmt.role from bc_mutation_tf bmt join biclusters bc on bmt.bicluster_id=bc.id join mutations m on m.id=bmt.mutation_id join tfs on tfs.id=bmt.tf_id where m.name=%s',
                       [mutation_name])
        result = [{"regulator": tf, "bicluster": bc, "role": MUTATION_TF_ROLES[role]}
                  for tf, bc, role in cursor.fetchall()]
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
        cursor.execute('select bc.name,bt.role from bc_tf bt join biclusters bc on bt.bicluster_id=bc.id join tfs on bt.tf_id=tfs.id where tfs.name=%s',
                       [tf_name])
        result = [{"bicluster": bc, "role": TF_BC_ROLES[role]}
                  for bc, role in cursor.fetchall()]
        return jsonify(regulator=tf_name, entries=result)
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
        return jsonify(num_biclusters=num_biclusters,
                       num_regulators=num_regulators,
                       num_mutations=num_mutations)
    finally:
        cursor.close()
        conn.close()

if __name__ == '__main__':
    app.debug = True
    app.secret_key = 'trstrestnorgp654g'
    app.run(host='0.0.0.0', debug=True)
