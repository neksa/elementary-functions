#!/usr/bin/env python
"""
CDD 3.11
Network generator: v3.0
"""

import string
import math
import re
from collections import defaultdict

CDD_dir = "../input/CDD/3.11/"
SCOP_dir = "../input/SCOP/1.75/"

net_path = "../output/"
debug_ef = False


#CHL00003       164415  cl08220 174986
fam_supfam_file = open(CDD_dir + "family_superfamily_links")
fam_supfam = {}
supfam_size = {}  # in number of families in the sf
for line in fam_supfam_file:
    fam_name, fam_id, sf_name, sf_id = line.strip().split("\t", 3)
    sf_id = int(sf_id)
    fam_id = int(fam_id)
    fam_supfam[fam_id] = sf_id

    if sf_id not in supfam_size:
        supfam_size[sf_id] = 0
    supfam_size[sf_id] += 1

fam_supfam_file.close()

#164413 CHL00001        rpoB    RNA polymerase beta subunit     1071
cddid_tbl = open(CDD_dir + "cddid.tbl")
cdd_descr = {}
cdd_codes = {}
for line in cddid_tbl:
    (cdd_id, cdd_code, short_name, description, pssm_length) = line.strip().split("\t", 4)
    cdd_descr[int(cdd_id)] = {'code': cdd_code, 'short': short_name, 'description': description}
    cdd_codes[cdd_code] = int(cdd_id)
cddid_tbl.close()

##################################################################################################
import gzip

seq_per_sf = {}
with gzip.open("../input/all_CDD-3.11.fasta.gz") as fasta:
    for line in fasta:
        if line[0] == '>':
            cdd_id = int(line.split(" ", 1)[0].split("|", 2)[2])
            sf = fam_supfam[cdd_id]
            if sf not in seq_per_sf:
                seq_per_sf[sf] = 0
            seq_per_sf[sf] += 1

            #if sf == 198671:
            #    print sf, cdd_id, seq_per_sf[sf]

#for sf, seq in seq_per_sf.iteritems():
#    print sf, seq

#######################################################################################################
#
# load CDD site information
#

site_keywords = defaultdict(set)
with open("../input/site_keywords.tab") as f:
    for line in f:
        ef, keyword = line.strip().split(' ', 1)
        keyword = keyword.replace(' ', '_')
        site_keywords[ef].add(keyword)

# OLD FORMAT:
# site_keywords = {}
# f = open("site_keywords.tab")
# for line in f:
#     site_class, keywords = line.strip().split("\t", 1)
#     site_keywords[site_class] = []
#     for kw in keywords.split(";"):
#         site_keywords[site_class].append(kw.strip())
# f.close()


def detect_elementary_functions(site_type, text):
    default_ef = {
        0: "GENERIC_OTHER",
        1: "GENERIC_ACTIVE_SITE",
        2: "GENERIC_PPI",
        3: "GENERIC_DNA_RNA",
        4: "GENERIC_ION",
        5: "GENERIC_CHEMICAL",
        6: "GENERIC_PTM"}

    detected_ef = set()
    text = text.strip()
    for ef, list_keywords in site_keywords.iteritems():
        detected = False
        for kw in list_keywords:
            # ATTENTION: CASE SENSITIVE STRING COMPARISON
            if -1 != text.find(kw):
                detected = True
                break
        if detected:
            detected_ef.add(ef)

    if len(detected_ef) == 0:
        detected_ef.add(default_ef.get(site_type, "UNKNOWN_EF_TYPE"))
    return detected_ef

#########################
cdd_sites = defaultdict(dict)
elementary_functions = set()

with open(CDD_dir + "cddannot.dat") as f:
    for line in f:
        (cdd, cdd_code, short_name, feature_number, feature_descr, str_based, ref_based, comments, addresses, site_type) = line.strip().split("\t", 9)

        positions = addresses.split(",")
        pos_list = []
        for position in positions:
            pos_list.append(int(position))

        site = int(feature_number)
        site_type = int(site_type)
        cdd = int(cdd)

        # classify elementary function of sites based on type/description here or load external classification
        ef = detect_elementary_functions(site_type, feature_descr.replace(" ", "_"))
        # ef = set()

        # ATTENTION!
        # POTENTIALLY DANGEROUS: if we do not know the EF of the site -- skip this site
        # Changed behaviour -- now we set the EF to uncategorized original feature description
        if len(ef) == 0:
            if debug_ef:
                ef.add(feature_descr.replace(" ", "_"))
            else:
                continue

        elementary_functions |= ef

        cdd_sites[cdd][site] = {
            "cdd_code": cdd_code, "name": short_name, "description": feature_descr,
            "positions": pos_list, "type": site_type, "ef": ef}

##################################################################################################
#
# Load profile -> CDD domain sequence coordinate mapping (calculated by model_alignment.py script):
#
"""
From profile_to_pssm_mapping.tab:

147     12      cd05229 29604089        2
293     1       cd05229 29604089        -9
235     0       cd05229 29604089        -10
"""
cdd_profile_coordinates = {}

f = open("../input/profile_to_pssm_mapping.tab", "r")
for line in f:
    profile, offset, cdd_code, gi, coordinate = line.strip().split("\t", 4)
    profile = int(profile)
    offset = int(offset)
    gi = int(gi)
    coordinate = int(coordinate)

    cdd = cdd_codes[cdd_code]

    if cdd not in cdd_profile_coordinates:
        cdd_profile_coordinates[cdd] = {}
    if profile not in cdd_profile_coordinates[cdd]:
        cdd_profile_coordinates[cdd][profile] = []

    # cdd_profile_coordinates[cdd][profile].append(coordinate)
    cdd_profile_coordinates[cdd][profile].append(offset - 10)

f.close()


##################################################################################################

#176    1.27    1.00    24      gnl|CDD|164534 CHL00141, rpl24, ribosomal protein L24; Validated
supfam_edges = {}
supfam_fam_covered = {}
sf_seq_hits = {1: {}, 10: {}, 100: {}, 1000: {}}

matches = {}
search_matches = gzip.open("../input/search_matches.E1000.tab.gz", "r")
for match in search_matches:
    (matrix, score, evalue, position, protein) = match.split("\t", 4)
    matrix = int(matrix)
    _gln, _cdd, domain = protein.split("|", 2)
    cdd_id, domain_description = domain.split(" ", 1)
    cdd = int(cdd_id)
    evalue = float(evalue)
    protein = domain_description.split(" ", 1)[1].split(",", 2)[1]
    idx = (matrix, protein)
    # we need unique hits matrix -> protein (in one position) with the best hit (according to evalue)
    if idx not in matches:
        matches[idx] = {}
        matches[idx]['evalue'] = evalue
        matches[idx]['cdd'] = cdd
        matches[idx]['position'] = int(position)
    else:
        if evalue < matches[idx]['evalue']:
            matches[idx]['evalue'] = evalue
            matches[idx]['cdd'] = cdd
            matches[idx]['position'] = int(position)

cdd_profile_offsets = {}

for (matrix, protein), match in matches.iteritems():
    cdd = match['cdd']
    evalue = match['evalue']

    # has a superfamily
    if cdd in fam_supfam:
        sf = fam_supfam[cdd]

        if (evalue == 1.0):
            if (matrix, sf) not in supfam_edges:
                supfam_edges[(matrix, sf)] = 0
            supfam_edges[(matrix, sf)] += 1

            if sf not in supfam_fam_covered:
                supfam_fam_covered[sf] = set()
            supfam_fam_covered[sf] |= set([cdd])

            if cdd not in cdd_profile_offsets:
                cdd_profile_offsets[cdd] = {}
            if matrix not in cdd_profile_offsets[cdd]:
                cdd_profile_offsets[cdd][matrix] = []
            cdd_profile_offsets[cdd][matrix].append(match['position'])

        for e in (1, 10, 100, 1000):
            if evalue > e:
                continue
            if (matrix, sf) not in sf_seq_hits[e]:
                sf_seq_hits[e][(matrix, sf)] = 0
            sf_seq_hits[e][(matrix, sf)] += 1

    else:
        print "SF UNKNOWN for FAM", cdd
#######################################################################################################
#
# calculate coverage of cdd sites (as a percent of covered site positions) by profile search
#
with open("../input/profile_signatures.csv") as sign_file:
    first = True
    profile_signature = {}
    profile_lengths = {}
    for line in sign_file:
        if first:
            first = not first
            continue
        profile, original_signature = line.strip().split("\t", 1)
        signature = re.sub("^[\-]*", "", original_signature)
        signature = re.sub("[\-]*$", "", signature)
        signature = re.sub("[\-]", "x", signature)
        profile = int(profile)
        signature = signature.strip()
        profile_signature[profile] = signature
        profile_lengths[profile] = len(original_signature)

cdd_site_profile_coverage = {}
for cdd, sites in cdd_sites.iteritems():
    cdd_site_profile_coverage[cdd] = {}
    for n, site in sites.iteritems():
        cdd_site_profile_coverage[cdd][n] = {}
        if cdd not in cdd_profile_offsets:
            continue
        for profile, offsets in cdd_profile_offsets[cdd].iteritems():
            cdd_site_profile_coverage[cdd][n][profile] = 0.0
            covered_pos = 0.0
            sum_pos = 0.0
            for pos in site["positions"]:
                # here we should use correct mapping of site <-> profile hits
                if cdd in cdd_profile_coordinates:
                    if profile in cdd_profile_coordinates[cdd]:
                        for x in cdd_profile_coordinates[cdd][profile]:
                            if x <= pos <= x + profile_lengths[profile]:
                                covered_pos += 1.0
                                break
                sum_pos += 1.0
            cdd_site_profile_coverage[cdd][n][profile] = covered_pos / sum_pos

#######################################################################################################
"""
site_types = {
    0 : "other",\
    1: "active site",\
    2: "polypeptide binding",\
    3: "nucleic acid binding",\
    4: "ion binding",\
    5: "chemical binding",\
    6: "posttranslational modification" }

features = {}
for cdd, sites in cdd_sites.iteritems():
    for n, site in sites.iteritems():

        hit_site = False
        for profile, coverage in cdd_site_profile_coverage[cdd][n].iteritems():
            if coverage > 0.3:
                hit_site = True
                break

        if not hit_site:
            continue

        site_type = site["type"]
        feature_descr = site["description"]

        if site_type not in features:
            features[site_type] = {}
        if feature_descr not in features[site_type]:
            features[site_type][feature_descr] = 0

        features[site_type][feature_descr] += 1

print "Number of different FOUND sites in each type group:"
for site_type, sites in features.iteritems():
    print site_type, site_types[site_type], len(sites)
print

print "Listing of FOUND sites:"
for site_type, sites in features.iteritems():
    print site_type, site_types[site_type]
    for site, counter in sites.iteritems():
        print "\t", site, counter
"""
#######################################################################################################
#s = float(seq_per_sf[199311])
#print seq_per_sf[199311], sf_seq_hits[1][199311],   sf_seq_hits[10][199311],  sf_seq_hits[100][199311], sf_seq_hits[1000][199311], supfam_size[199311]
#print sf_seq_hits[1][199311]/s ,   sf_seq_hits[10][199311]/ s,  sf_seq_hits[100][199311]/ s, sf_seq_hits[1000][199311]/ s, supfam_size[199311]
# Profile attributes
"""
PROFILE 4
BEGIN
ORIGIN binding_Ca_Zn; PDOC00018; PS00018; 1
MATRIX ID=0 K=450 L=30
"""

profile_origin = {}
origin = ""
matrix_id = 0
with open("../input/extracted.matrix", "r") as mfile:
    for line in mfile:
        if line.startswith("BEGIN"):
            origin = ""
            matrix_id = 0

        if line.startswith("ORIGIN"):
            origin = line[7:].strip()

        if line.startswith("MATRIX"):
            mline = line[7:].strip()
            (id_line, _, _) = mline.split(" ", 3)
            (ID, matrix_id) = id_line.split("=", 1)
            matrix_id = int(matrix_id)
            profile_origin[matrix_id] = origin


#######################################################################################################
# LOAD SCOP #

filename_des = SCOP_dir + "dir.des.scop.txt_1.75"
scan_des = re.compile('(.*?)\t(.*?)\t(.*?)\t.*?\t(.*)')
scan_class = re.compile('(.+?)\.?(.+)?\.?(.+)?\.?(.+)?')
with open(filename_des, 'r') as file_des:
    scop_attr = {}
    scop_fam_fold = {}
    scop_fold_fam = {}
    scop_fold = 0
    for line in file_des:
        m = scan_des.match(line)
        if m is not None:
            #print m.group(0)
            sunid = int(m.group(1))
            level = m.group(2)
            classification = m.group(3)
            description = m.group(4)

            if level == 'cf':
                g = scan_class.match(classification)
                if g.group(1) is not None and g.group(2) is not None:
                    scop_fold = sunid  # g.group(1) + "." + g.group(2)
                    #print "cf", scop_fold, description
                    scop_attr[scop_fold] = g.group(1) + "." + g.group(2) + " " + description
                scop_fold_fam[sunid] = 0
            if level == 'fa':
                if sunid not in scop_fam_fold:
                    scop_fam_fold[sunid] = scop_fold
                scop_fold_fam[scop_fold] += 1
            #print sunid, level, classification, description

#CDD_masters_cd_only.ass
"""
gnl|CDD|28986   0044234 1-89    2.15e-23        5       PVITSISPSSGPvsGGTEVTITGSNFGSGSN--LRVTFGGgVPCSVL--SVSSTAIVCTTPPYANPGPGPVEVTVDRGnggITSSPLTFTYVP   0.1     66214   81279
gnl|CDD|28997   0049520 1-107   2.18e-45        3       KMRLRPWLVEQVDSGTYPGLIWLDEEKTIFRIPWKHAARHDVQEADAKIFKAWAVERGIYQPGGTP-DPAEWKARLLCALRSSRGFEEVKDKSKdTPGDPHRVYRLLP    0.1     123768  46859
"""
cdd_sf_scop_folds = defaultdict(lambda: defaultdict(lambda: 1))
with open("../input/cddmasters_cd_only.ass") as ass_file:
    for line in ass_file:
        fields = line.split("\t", 8)
        if fields[7] == '-':
            continue
        evalue = float(fields[3])

        cdd = int(fields[0].split("|", 2)[2])
        sf = fam_supfam[cdd]
        scop_domain = int(fields[7])
        scop_family = int(fields[8])
        fold = scop_fam_fold[scop_family]

        # defaultdict:
        # if sf not in cdd_sf_scop_folds:
        #     cdd_sf_scop_folds[sf] = {}
        # if fold not in cdd_sf_scop_folds[sf]:
        #     cdd_sf_scop_folds[sf][fold] = 1

        # record the minimal HMMER E-value
        if evalue < cdd_sf_scop_folds[sf][fold]:
            cdd_sf_scop_folds[sf][fold] = evalue

        #print cdd, sf, scop_fam_fold[scop_family], scop_attr[scop_fam_fold[scop_family]]

#print cdd_sf_scop_folds
#######################################################################################################
scop_folds_in_net = set()
SITE_COVERAGE_THRESHOLD = 0.3

# Create network
with open(net_path+"net.sif", "w") as supfam_net:
    for profile, sf in supfam_edges.iterkeys():
        supfam_net.write("profile_"+str(profile) + " EFL " + "cdd_"+str(sf) + "\n")

    sf_list = [sf for (profile, sf) in supfam_edges.iterkeys()]
    unique_sf = set(sf_list)

    for sf in unique_sf:
        for fold in cdd_sf_scop_folds[sf].iterkeys():
            supfam_net.write("cdd_"+str(sf) + " SCOP " + "scop_"+str(fold) + "\n")
            scop_folds_in_net.add(fold)

    for cdd, sites in cdd_sites.iteritems():
        sf = fam_supfam[cdd]
        if sf not in unique_sf:
            continue

        for n, site in sites.iteritems():
            # is there (SITE_COVERAGE_THRESHOLD) coverage of this particular site by any profile?
            is_covered = False
            efs = set()
            for ef in site["ef"]:
                if ef in ('active_site', 'allosteric_site', 'PPI'):
                    ef = "{}_{}_{}".format(ef, sf, n)
                efs.add(ef)
            elementary_functions |= efs

            for profile, coverage in cdd_site_profile_coverage[cdd][n].iteritems():
                if coverage > SITE_COVERAGE_THRESHOLD:
                    is_covered = True
                    #supfam_net.write("profile_%d SITE_HIT cdd_site_%d_%d\n" % (profile, cdd, n))
                    for ef in efs:
                        supfam_net.write("profile_%d EF ef_%s\n" % (profile, ef))
            # don't include sites that have no coverage
            if is_covered:
                #supfam_net.write("cdd_%d CDD_FEATURE cdd_site_%d_%d\n" % (sf, cdd, n))
                # for ef in site["ef"]:
                for ef in efs:
                    supfam_net.write("cdd_%d CDD_EF ef_%s\n" % (sf, ef))

    #for profile in profile_ef.iterkeys():
    #    supfam_net.write("profile_"+str(profile) + " EF " + "ef_"+ profile_ef[profile] + "\n")

# Create edge attributes
with open(net_path + "edge_attributes.tab", "w") as o:
    o.write("edge\tweight\tE1_cov\tE1_hits\tE10_cov\tE10_hits\tE100_cov\tE100_hits\tE1000_cov\tE1000_hits\n")
    for (profile, sf) in supfam_edges.iterkeys():
        hits = {}
        cov = {}
        for e in (1, 10, 100, 1000):
            hits[e] = sf_seq_hits[e][(profile, sf)]
            cov[e] = float(sf_seq_hits[e][(profile, sf)]) / seq_per_sf[sf]
        weight = float(sf_seq_hits[100][(profile, sf)]) / seq_per_sf[sf]
        o.write("profile_{} (EFL) cdd_{}\t{:.2}\t{:.2}\t{}\t{:.2}\t{}\t{:.2}\t{}\t{:.2}\t{}\n".format(
                profile, sf, weight, cov[1], hits[1], cov[10], hits[10], cov[100], hits[100], cov[1000], hits[1000]))

    # SCOP connection between cdd superfamily and SCOP fold
    for sf in unique_sf:
        for fold, evalue in cdd_sf_scop_folds[sf].iteritems():
            #print sf, fold, evalue,
            if evalue > 0:
                weight = 0 - math.log10(evalue)
                if weight < 1:
                    weight = 1.0
                if weight > 100:
                    weight = 100.0
                weight = weight / 100.0
            else:
                weight = 1.0
            #print weight
            o.write("cdd_{} (SCOP) scop_{}\t{:.2}\t0.0\t0\t0.0\t0\t0.0\t0\t0.0\t0\n".format(sf, fold, weight))

    """
    # SITE_HIT connection between profile and site
    for cdd, sites in cdd_sites.iteritems():
        sf = fam_supfam[cdd]
        if sf not in unique_sf:
            continue

        for n, site in sites.iteritems():
            for profile, coverage in cdd_site_profile_coverage[cdd][n].iteritems():
                if coverage > SITE_COVERAGE_THRESHOLD:
                    #FIXME THIS IS INCORRECT!
                    #eda.write("profile_%d (SITE_HIT) cdd_site_%d_%d = %.4f\n" % (profile, cdd, n, coverage))
                    pass

    for profile in profile_ef.iterkeys():
        eda.write("profile_%s (EF) ef_%s = 0.3\n" % (profile, profile_ef[profile]))
    """

#######################################################################################################
# Common node attribute names:
#
# ID = cdd_number, ef_name, profile_number, scop_sunid
# type = cdd, ef, profile, scop
# label = a short node label
# popup = longer name or description on mouse over
# size = size of the node, normalized from 1 to 100, either linear or logarithmic scale, depending on the distribution shape
# URL = http://URL?id

node_attributes = []

# Superfamily attributes
"""

attr = {}
attr['ID'] =
attr['type'] =
attr['label'] =
attr['popup'] =
attr['size'] =
attr['URL'] =
attr['Nfamilies'] =
attr['short_name'] =
attr['long_name'] =
attr['description'] =
"""

max_seq = float(max(seq_per_sf.values()))

for sf in unique_sf:
    attr = {}
    attr['ID'] = 'cdd_' + str(sf)
    attr['type'] = 'cdd'
    attr['label'] = cdd_descr[sf]['short']
    attr['popup'] = cdd_descr[sf]['description'].split(".", 1)[0]
    # min size = 40, max size = 100
    attr['size'] = "%.4f" % (40.0 + (60.0 * math.log(seq_per_sf[sf]) / math.log(max_seq)))
    #attr['URL'] = "http://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid="+sf
    attr['URL'] = "http://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid="+cdd_descr[sf]['code']
    attr['Nfamilies'] = supfam_size[sf]
    attr['Nsequences'] = seq_per_sf[sf]
    attr['short_name'] = cdd_descr[sf]['short']
    attr['cdd_code'] = cdd_descr[sf]['code']
    #attr['long_name'] = ""
    attr['description'] = cdd_descr[sf]['description']
    node_attributes.append(attr)

max_families = float(max(scop_fold_fam.values()))

for fold in scop_folds_in_net:
    attr = {}
    attr['ID'] = 'scop_' + str(fold)
    attr['type'] = 'scop'
    attr['label'] = scop_attr[fold].split(" ", 1)[0]
    attr['popup'] = scop_attr[fold]
    # min size = 40, max size = 100
    attr['size'] = "%.4f" % (40.0 + (60.0 * math.log(scop_fold_fam[fold]) / math.log(max_families)))
    attr['URL'] = "http://scop.mrc-lmb.cam.ac.uk/scop/search.cgi?ver=1.75&key=" + str(fold)
    attr['Nfamilies'] = scop_fold_fam[fold]
    #attr['short_name'] =
    #attr['long_name'] =
    attr['description'] = scop_attr[fold].split(" ", 1)[1]
    node_attributes.append(attr)

for profile in profile_origin.iterkeys():
    attr = {}
    attr['ID'] = 'profile_' + str(profile)
    attr['type'] = 'profile'
    attr['label'] = profile
    attr['popup'] = profile_signature[profile]
    attr['size'] = 30.0
    attr['URL'] = "http://neksa.net/EF/profile/" + str(profile) + ".html"
    attr['description'] = profile_origin[profile]
    # attr['ef_group'] = profile_ef[profile].replace(" ", "_").replace(":", "_")
    node_attributes.append(attr)

for ef in elementary_functions:
    attr = {}
    attr['ID'] = 'ef_' + str(ef).replace(" ", "_").replace(":", "_")
    attr['type'] = 'ef'
    attr['label'] = ef
    attr['popup'] = "Group of elementary functions " + ef
    attr['size'] = 45.0
    attr['URL'] = "http://neksa.net/EF/ef/" + str(ef) + ".html"
    node_attributes.append(attr)

#######################################################################################################
# export attributes
attribute_names = ['ID', 'type', 'label', 'popup', 'size', 'URL']
for attr_row in node_attributes:
    for attr in attr_row.iterkeys():
        if attr not in attribute_names:
            attribute_names.append(attr)

noa = open(net_path+"node_atributes.tab", "w")
noa.write(string.join(attribute_names, "\t") + "\n")
for node in node_attributes:
    export_attributes = []
    for attr in attribute_names:
        if attr in node:
            export_attributes.append(str(node[attr]))
        else:
            export_attributes.append("")
    noa.write(string.join(export_attributes, "\t") + "\n")
noa.close()
