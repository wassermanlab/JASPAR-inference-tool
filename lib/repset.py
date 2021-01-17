#!/usr/bin/env python

"""
Source: https://github.com/mlibbrecht/submodular_sequence_repset
"""

import argparse
import copy
import heapq
import logging
import math
import os
import subprocess
import sys

from Bio import SeqIO

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", required=True, help="output directory")
    parser.add_argument("--seqs", required=True, help="input sequences (in FASTA format)")
    parser.add_argument("--mixture", type=float, default=0.5, help="mixture parameter determining the relative weight of facility-location relative to sum-redundancy (default=0.5)")

    args = parser.parse_args()
    workdir = os.path.abspath(args.outdir)

    assert args.mixture >= 0.0
    assert args.mixture <= 1.0

    if not os.path.isdir(workdir):
        os.makedirs(workdir)

    # Logging
    logging.basicConfig(format="%(asctime)s %(levelname)s:%(message)s")
    logger = logging.getLogger("log")
    logger.setLevel(logging.DEBUG)
    # create file handler which logs even debug messages
    fh = logging.FileHandler(os.path.join(workdir, "stdout.txt"))
    fh.setLevel(logging.DEBUG) # >> this determines the file level
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)# >> this determines the output level
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)
    # add the handlers to logger
    logger.addHandler(ch)
    logger.addHandler(fh)
else:
    logger = logging.getLogger("log")

###############################################################
###############################################################
# run_psiblast()
# ---------------
# Runs psiblast to get a similarity matrix between sequences.
# Input:
#   - workdir: directory for working files
#   - seqs: fasta file with sequences
# Output:
#   - db with sequence and similarity information.
###############################################################
###############################################################
def run_psiblast(workdir, seqs):
    # Create psiblast db
    psiblast_db = os.path.join(workdir, "db")
    if not os.path.exists(psiblast_db):
        cmd = ["makeblastdb",
          "-in", seqs,
          "-input_type", "fasta",
          "-out", psiblast_db,
          "-dbtype", "prot"]
        logger.info(" ".join(cmd))
        subprocess.check_call(cmd)
    # Run psiblast
    psiblast_out = os.path.join(workdir, "psiblast_result.tab")
    if not os.path.exists(psiblast_out):
        cmd = ["psiblast",
          "-query", seqs,
          "-db", psiblast_db,
          "-num_iterations", "6",
          "-outfmt", "6 qseqid sseqid pident length mismatch evalue bitscore",
          "-comp_based_stats", "1",
          "-seg", "yes",
          "-out", psiblast_out
        ]
        logger.info(" ".join(cmd))
        subprocess.check_call(cmd)
    # Read psiblast output
    db = {}
    fasta_sequences = SeqIO.parse(open(seqs),'fasta')
    for seq in fasta_sequences:
        seq_id = seq.id
        db[seq_id] = {"neighbors": {}, "in_neighbors": {}, "seq": seq.seq}
    with open(psiblast_out, "r") as f:
        for line in f:
            if line.strip() == "": continue
            if line.startswith("Search has CONVERGED!"): continue
            line = line.split()
            seq_id1 = line[0]
            seq_id2 = line[1]
            pident = float(line[2])
            evalue = line[5]
            try:
                log10_e = math.log10(float(evalue))
            except:
                # fix 0 error
                log10_e = math.log10(sys.float_info.min)
            if float(evalue) <= 1e-2:
                db[seq_id2]["neighbors"][seq_id1] = {"log10_e": log10_e, "pct_identical": pident}
                db[seq_id1]["in_neighbors"][seq_id2] = {"log10_e": log10_e, "pct_identical": pident}
    return(db)

###############################################################
###############################################################
# Similarity functions
###############################################################
###############################################################

def sim_from_neighbor(sim, d):
    return(sim(d["log10_e"], d["pct_identical"]))

def fraciden(log10_e, pct_identical):
    return(float(pct_identical) / 100)

###############################################################
###############################################################
# summaxacross
# AKA facility location
###############################################################
###############################################################

def summaxacross_eval(db, seq_ids, sim):
    max_sim = {seq_id:0 for seq_id in db}
    for chosen_seq_id in seq_ids:
        for neighbor_seq_id, d in db[chosen_seq_id]["in_neighbors"].items():
            if neighbor_seq_id in max_sim:
                sim_val = sim(d["log10_e"], d["pct_identical"])
                if sim_val > max_sim[neighbor_seq_id]:
                    max_sim[neighbor_seq_id] = sim_val
            else:
                pass

    return(sum(max_sim.values()))

# summaxacross data:
# Who each example is represented by
# "examples": {seq_id: (representive, val) }
# Who each represntative represents
# "representatives" {seq_id: {example: val}}

summaxacross_base_data = lambda db, sim: {"examples": {seq_id: (None, 0) for seq_id in db},
                                              "representatives": {}}

summaxacross_full_data = lambda db, sim: {"examples": {seq_id: (seq_id, sim_from_db(db, sim, seq_id, seq_id)) for seq_id in db},
                                              "representatives": {seq_id: {seq_id: sim_from_db(db, sim, seq_id, seq_id)} for seq_id in db}}

def summaxacross_diff(db, seq_id, sim, data):
    diff = 0
    for neighbor_seq_id, d in db[seq_id]["in_neighbors"].items():
        if neighbor_seq_id in data["examples"]:
            sim_val = sim(d["log10_e"], d["pct_identical"])
            if sim_val > data["examples"][neighbor_seq_id][1]:
                diff += sim_val - data["examples"][neighbor_seq_id][1]
        else:
            pass
            #raise Exception("Found node with neighbor not in set")
    return(diff)

def summaxacross_update(db, seq_id, sim, data):
    data = copy.deepcopy(data)
    data["representatives"][seq_id] = {}
    for neighbor_seq_id, d in db[seq_id]["in_neighbors"].items():
        if neighbor_seq_id in data["examples"]:
            sim_val = sim(d["log10_e"], d["pct_identical"])
            if sim_val > data["examples"][neighbor_seq_id][1]:
                data["examples"][neighbor_seq_id] = (seq_id, sim_val)
                data["representatives"][seq_id][neighbor_seq_id] = sim_val
        else:
            pass
    return(data)

# O(D^2)
def summaxacross_negdiff(db, seq_id, sim, data):
    diff = 0
    # For each neighbor_seq_id that was previously represented by seq_id
    new_representatives = (set(data["representatives"].keys()) - set([seq_id]))
    for neighbor_seq_id, d in data["representatives"][seq_id].items():
        # Find the next-best representative for neighbor_seq_id
        candidate_ids = set(db[neighbor_seq_id]["neighbors"].keys()) & new_representatives
        if len(candidate_ids) == 0:
            diff += -d
        else:
            best_id = max(candidate_ids, key=lambda x: sim_from_db(db, sim, neighbor_seq_id, x))
            diff += sim_from_db(db, sim, neighbor_seq_id, best_id) - d
    return(diff)

# O(D^2)
def summaxacross_negupdate(db, seq_id, sim, data):
    data = copy.deepcopy(data)
    new_representatives = (set(data["representatives"].keys()) - set([seq_id]))
    for neighbor_seq_id, d in data["representatives"][seq_id].items():
        # Find the next-best representative for neighbor_seq_id
        candidate_ids = set(db[neighbor_seq_id]["neighbors"].keys()) & new_representatives
        if len(candidate_ids) == 0:
            data["examples"][neighbor_seq_id] = (None, 0)
        else:
            best_id = max(candidate_ids, key=lambda x: sim_from_db(db, sim, neighbor_seq_id, x))
            data["examples"][neighbor_seq_id] = (best_id, sim_from_db(db, sim, neighbor_seq_id, best_id))
            data["representatives"][best_id][neighbor_seq_id] = sim_from_db(db, sim, neighbor_seq_id, best_id)
    del data["representatives"][seq_id]
    return(data)

summaxacross = {
    "eval": summaxacross_eval,
    "diff": summaxacross_diff,
    "negdiff": summaxacross_negdiff,
    "update": summaxacross_update,
    "negupdate": summaxacross_negupdate,
    "base_data": summaxacross_base_data,
    "full_data": summaxacross_full_data,
    "name": "summaxacross"
}

###############################################################
###############################################################
# sumsumwithin
###############################################################
###############################################################

def bisim(db, sim, seq_id1, seq_id2):
    ret = 0
    if seq_id2 in db[seq_id1]["neighbors"]:
        d = db[seq_id1]["neighbors"][seq_id2]
        ret += sim(d["log10_e"], d["pct_identical"])
    if seq_id1 in db[seq_id2]["neighbors"]:
        d = db[seq_id2]["neighbors"][seq_id1]
        ret += sim(d["log10_e"], d["pct_identical"])
    return(ret)

def sumsumwithin_eval(db, seq_ids, sim):
    seq_ids = set(seq_ids)
    s = 0
    for chosen_id in seq_ids:
        for neighbor, d in db[chosen_id]["neighbors"].items():
            if chosen_id == neighbor: continue
            if neighbor in seq_ids:
                s += -sim(d["log10_e"], d["pct_identical"])
    return(s)

sumsumwithin_base_data = lambda db, sim: set()
sumsumwithin_full_data = lambda db, sim: set(db.keys())

def sumsumwithin_diff(db, seq_id, sim, data):
    diff = 0
    data = data | set([seq_id])
    for neighbor, d in db[seq_id]["neighbors"].items():
        if seq_id == neighbor: continue
        if not (neighbor in data): continue
        diff += -sim_from_neighbor(sim, d)
        #neighbor_bisim = bisim(db, sim, seq_id, neighbor)
        #diff += -neighbor_bisim
    for neighbor, d in db[seq_id]["in_neighbors"].items():
        if seq_id == neighbor: continue
        if not (neighbor in data): continue
        diff += -sim_from_neighbor(sim, d)
    return(diff)

def sumsumwithin_update(db, seq_id, sim, data):
    data.add(seq_id)
    return(data)

def sumsumwithin_negdiff(db, seq_id, sim, data):
    diff = 0
    #data = data - set([seq_id])
    for neighbor, d in db[seq_id]["neighbors"].items():
        if seq_id == neighbor: continue
        if not (neighbor in data): continue
        #neighbor_bisim = bisim(db, sim, seq_id, neighbor)
        #diff -= -neighbor_bisim
        diff += sim_from_neighbor(sim, d) # removing a penalty
    for neighbor, d in db[seq_id]["in_neighbors"].items():
        if seq_id == neighbor: continue
        if not (neighbor in data): continue
        diff += sim_from_neighbor(sim, d) # removing a penalty
    return(diff)

def sumsumwithin_negupdate(db, seq_id, sim, data):
    data.remove(seq_id)
    return(data)

sumsumwithin = {
    "eval": sumsumwithin_eval,
    "diff": sumsumwithin_diff,
    "negdiff": sumsumwithin_negdiff,
    "update": sumsumwithin_update,
    "negupdate": sumsumwithin_negupdate,
    "base_data": sumsumwithin_base_data,
    "full_data": sumsumwithin_full_data,
    "name": "sumsumwithin"
}

###############################################################
###############################################################
# MixtureObjective
# ------------------------
# Create a mixture objective with:
# MixtureObjective([summaxacross, sumsumwithin], [0.1, 1.2])
# Must be used with a sim of the form
# [sim1, sim2]
# (same number of sims as objectives)
###############################################################
###############################################################

class MixtureObjective(object):
    def __init__(self, objectives, weights):
        self.objectives = objectives
        self.weights = weights
        self.name = "mix-" + "-".join(["{0}({1})".format(objective["name"], self.weights[i]) for i,objective in enumerate(self.objectives)])

    def __getitem__(self, key):
        return(self.__getattribute__(key))

    def __contains__(self, item):
        all_contain = True
        for i, objective in enumerate(self.objectives):
            all_contain = all_contain and (item in objective)
        return(all_contain)

    def eval(self, db, seq_ids, sims):
        s = 0
        for i, objective in enumerate(self.objectives):
            s += self.weights[i]*objective["eval"](db, seq_ids, sims[i])
        return(s)

    def diff(self, db, seq_id, sims, datas):
        s = 0
        for i, objective in enumerate(self.objectives):
            s += self.weights[i]*objective["diff"](db, seq_id, sims[i], datas[i])
        return(s)

    def negdiff(self, db, seq_id, sims, datas):
        s = 0
        for i, objective in enumerate(self.objectives):
            s += self.weights[i]*objective["negdiff"](db, seq_id, sims[i], datas[i])
        return(s)

    def update(self, db, seq_id, sims, datas):
        new_datas = []
        for i, objective in enumerate(self.objectives):
            new_datas.append(objective["update"](db, seq_id, sims[i], datas[i]))
        return(new_datas)

    def negupdate(self, db, seq_id, sims, datas):
        new_datas = []
        for i, objective in enumerate(self.objectives):
            new_datas.append(objective["negupdate"](db, seq_id, sims[i], datas[i]))
        return(new_datas)

    def base_data(self, db, sims):
        datas = []
        for i, objective in enumerate(self.objectives):
            datas.append(objective["base_data"](db, sims[i]))
        return(datas)

    def full_data(self, db, sims):
        datas = []
        for i, objective in enumerate(self.objectives):
            datas.append(objective["full_data"](db, sims[i]))
        return(datas)

###############################################################
###############################################################
# Optimization algorithms
# -------------------------
# Each returns either a specific
# subset or an order.
###############################################################
###############################################################

def accelerated_greedy_selection(
    db, objective, sim, max_evals=float("inf"), diff_approx_ratio=1.0,
    repset_size=float("inf"), target_obj_val=float("inf")
):
    """
    Returns an order.
    """
    assert diff_approx_ratio <= 1.0
    repset = []
    pq = [(-float("inf"), seq_id) for seq_id in db]
    objective_data = objective["base_data"](db, sim)
    cur_objective = 0
    num_evals = 0
    while (len(repset) < repset_size) and (len(pq) > 1) and (cur_objective < target_obj_val):
        possible_diff, seq_id = heapq.heappop(pq)
        diff = objective["diff"](db, seq_id, sim, objective_data)
        next_diff = -pq[0][0]
        num_evals += 1
        if (num_evals >= max_evals) or (((diff - next_diff) / (abs(diff)+0.01)) >= (diff_approx_ratio - 1.0)):
            repset.append(seq_id)
            objective_data = objective["update"](db, seq_id, sim, objective_data)
            cur_objective += diff
            #assert(abs(cur_objective - objective["eval"](db, repset, sim)) < 1e-3)
            #if (len(repset) % 100) == 0: logger.debug("Accelerated greedy iteration: %s", len(repset))
            #print len(repset), num_evals # XXX
            #print len(repset), cur_objective
            num_evals = 0
        else:
            heapq.heappush(pq, (-diff, seq_id))
    if len(pq) == 1:
        repset.append(pq[0][1])
    return(repset)

###############################################################
###############################################################
# Run optimization and output
###############################################################
###############################################################

if __name__ == "__main__":
    db = run_psiblast(workdir, os.path.abspath(args.seqs))
    objective = MixtureObjective([summaxacross, sumsumwithin], [args.mixture, 1.0-args.mixture])
    logger.info("-----------------------")
    logger.info("Starting mixture of summaxacross and sumsumwithin with weight %s...", args.mixture)
    sim, sim_name = ([fraciden, fraciden], "fraciden-fraciden")
    repset_order = accelerated_greedy_selection(db, objective, sim)

    with open(os.path.join(workdir, "repset.txt"), "w") as f:
        for seq_id in repset_order:
            f.write(seq_id)
            f.write("\n")
