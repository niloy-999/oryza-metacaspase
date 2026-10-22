"""Canonical rice metacaspase nomenclature (Huang et al. 2015, Int J Mol Sci 16:16216) and the manuscript's ad-hoc 'oryza-N' names."""
from collections import OrderedDict
# OsMC name -> (MSU locus, type)  [Type I = OsMC1-3, Type II = OsMC4-8; Huang et al. 2015; Wang & Zhang 2014]
OSMC = OrderedDict([
 ("OsMC1", ("LOC_Os03g27120", "I")), ("OsMC2", ("LOC_Os03g27210", "I")), ("OsMC3", ("LOC_Os03g27170", "I")),
 ("OsMC4", ("LOC_Os05g41660", "II")), ("OsMC5", ("LOC_Os05g41670", "II")), ("OsMC6", ("LOC_Os01g58580", "II")),
 ("OsMC7", ("LOC_Os11g04010", "II")), ("OsMC8", ("LOC_Os03g27190", "II")),
])
LOCUS2OSMC = {v[0]: k for k, v in OSMC.items()}
# manuscript's provisional names ("metacaspase_1_oryza-N") -> locus
MS_NAME = OrderedDict([
 ("oryza-1", "LOC_Os01g58580"), ("oryza-2", "LOC_Os03g27190"), ("oryza-3", "LOC_Os03g27120"),
 ("oryza-4", "LOC_Os03g27210"), ("oryza-5", "LOC_Os03g27170"), ("oryza-6", "LOC_Os09g36640"),
 ("oryza-7", "LOC_Os05g41660"), ("oryza-8", "LOC_Os05g41670"), ("oryza-9", "LOC_Os11g04010")])
LOCUS2MS = {v: k for k, v in MS_NAME.items()}
# values printed in the manuscript's Table 2 (PDF, Nov 2025): chr, start, end, strand, CDS nt, 'protein length', MW kDa, pI, GRAVY, intron:exon
MS_TABLE2 = {
 "oryza-1": (1, 33869141, 33871840, "+", 1254, 418, 45.89, 5.51, -0.595, "1:02"),
 "oryza-2": (3, 15567898, 15572135, "-", 909, 303, 33.52, 6.46, -0.312, "4:05"),
 "oryza-3": (3, 15530612, 15532459, "-", 1107, 369, 39.66, 7.85, -0.188, "3:04"),
 "oryza-4": (3, 15584825, 15587152, "-", 1110, 370, 39.05, 7.91, 0.03, "2:03"),
 "oryza-5": (3, 15563468, 15564927, "-", 1203, 401, 42.17, 7.01, -0.156, "3:04"),
 "oryza-6": (9, 21142247, 21143387, "-", 1140, 380, 38.49, 6.17, -0.34, "0:01"),
 "oryza-7": (5, 24383179, 24384628, "-", 1230, 410, 44.55, 6.35, -0.565, "1:02"),
 "oryza-8": (5, 24384900, 24387468, "-", 1263, 421, 46.02, 6.04, -0.565, "1:02"),
 "oryza-9": (11, 1616647, 1618300, "-", 1026, 342, 35.66, 5.24, -0.074, "0:01"),
}
AT_MC = OrderedDict([("AtMC1","AT1G02170"),("AtMC2","AT4G25110"),("AtMC3","AT5G64240"),("AtMC4","AT1G79340"),("AtMC5","AT1G79330"),
                     ("AtMC6","AT1G79320"),("AtMC7","AT1G79310"),("AtMC8","AT1G16420"),("AtMC9","AT5G04200")])
AT_TYPE = {"AtMC1":"I","AtMC2":"I","AtMC3":"I","AtMC4":"II","AtMC5":"II","AtMC6":"II","AtMC7":"II","AtMC8":"II","AtMC9":"II"}
