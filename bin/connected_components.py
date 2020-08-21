#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#usage: python connected_components.py fichier_sommets.txt fichier_arretes.txt
# >connected_components.txt
#Form groupe between lists.

import csv
import argparse as ap

class node(object):

    def __init__(self, label):
        # label is not necessary for the algorithm but makes it easier to print
        self.label = label
        self.p = None
        self.rank = 0

    def __str__(self):
        # enabling to print simply the label when calling print(u) or str(u)
        return self.label


def make_set(e):
    """
    Create a new node and make it point to itself
    We do this in two steps :
      - construct a node
      - make it point to itself
    because pointing to itself within the constructor is not possible
    """
    n = node(e)
    n.p = n
    return n


def link(x, y):
    """
    Link two nodes choosing the higher in the "tree", or the older (of higher
    rank) as the representative. We cannot make an arbitrary choice for
    the representative because we need a single one for each connected
    component
    """
    if x.rank > y.rank:
        y.p = x
    else:
        x.p = y
        if x.rank == y.rank:
            # nodes of equal rank, y becomes the father of x, hence y+=1
            y.rank = y.rank + 1


def find_set(x):
    """
    Find the set of a node, in fact the representative
    Each time a node is crossed it's representative is updated
    :param x: The node
    :return: the representative of it's set
    """
    if x is not x.p:
        x.p = find_set(x.p)
    return x.p


def union(x, y):
    """
    Join two elements in fact their corresponding set
    :param x: a node
    :param y: a second node
    """
    link(find_set(x), find_set(y))


def connected_components(sommets, adjacencies):
    """
    Compute the connected components
    :param adjacencies: a dictionnary of adjacencies
         keys are the nodes
         value is a set of node
    :return the connected components as a dictionnary
       key : the representative of a set
       value: the list of nodes
    """
    nodes = dict()
    # for each node u, a set is constructed
    for e in sommets:
        nodes[e[0]] = make_set(e[0])
    # for each node u
    for u in adjacencies:
        # for each adjacency v join u and v (if they are not already joined)
        u1 = u[0]
        u2 = u[1]
        if find_set(nodes[u1]) != find_set(nodes[u2]):
            union(nodes[u1], nodes[u2])
    # Now we seek for the different representatives,
    # hence the different connected components
    components = dict()
    for u in nodes:
        representative = find_set(nodes[u])
        if representative not in components:
            components[representative] = [u]
        else:
            components[representative].append(u)
    return components



parser = ap.ArgumentParser()
parser.add_argument("fichier_sommets")
parser.add_argument("fichier_arretes")
args = parser.parse_args()

sommets = []
with open (args.fichier_sommets,"r") as fichier:
    sommet = csv.reader (fichier, delimiter ='\t')
    for line in sommet:
        sommets.append(line)

adjacencies = []
with open (args.fichier_arretes, "r") as fichier2:
    arretes = csv.reader (fichier2, delimiter = '\t')
    for line in arretes:
        adjacencies.append(line)


# Finding the connected components
components = connected_components(sommets, adjacencies)

# Printing them
for c in components:
    i=0
    output = ''
    while i<len(components[c]):
        if output == '':
            output = components[c][i]
        else:
            output = '\t'.join((output,components[c][i]))
        i+=1
    print(output)
