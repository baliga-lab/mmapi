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

@app.route('/api/v1.0.0/bicluster_genes/<cluster_id>')
def bicluster_genes(cluster_id):
    """all genes in the specified bicluster"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute('select g.name from biclusters bc join bicluster_genes bcg on bc.id=bcg.bicluster_id join genes g on g.id=bcg.gene_id where bc.name=%s', [cluster_id])
        genes = [row[0] for row in cursor.fetchall()]
        return jsonify(cluster=cluster_id, genes=genes)
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
