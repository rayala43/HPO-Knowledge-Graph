"""
data/hpo_data.py
----------------
Curated HPO gene–phenotype–disease dataset.
32 conditions · 108 gene–disease links · 96 unique genes

Sources (all peer-reviewed):
  - HPO database  : hpo.jax.org
  - OMIM          : omim.org
  - Orphanet      : orpha.net
  - GWAS Catalog  : ebi.ac.uk/gwas
  - ClinVar       : clinvar.ncbi.nlm.nih.gov
"""

# HPO term labels (used for display in network tooltips)
HPO_TERM_LABELS = {
    "HP:0000819": "Diabetes mellitus",
    "HP:0001513": "Obesity",
    "HP:0000855": "Insulin resistance",
    "HP:0003362": "Insulin receptor defect",
    "HP:0000825": "Hyperinsulinemia",
    "HP:0005978": "Type 2 diabetes",
    "HP:0003074": "Hyperglycemia",
    "HP:0001508": "Failure to thrive",
    "HP:0003124": "Hypercholesterolemia",
    "HP:0001695": "Cardiac arrest",
    "HP:0002621": "Atherosclerosis",
    "HP:0001297": "Stroke",
    "HP:0001659": "Aortic regurgitation",
    "HP:0000822": "Hypertension",
    "HP:0001873": "Thrombocytopenia",
    "HP:0001974": "Leukocytosis",
    "HP:0003002": "Breast carcinoma",
    "HP:0100615": "Ovarian neoplasm",
    "HP:0002664": "Neoplasm",
    "HP:0002671": "Basal cell carcinoma",
    "HP:0001903": "Anemia",
    "HP:0000818": "Abnormality of endocrine system",
    "HP:0003003": "Colon carcinoma",
    "HP:0005275": "Hamartoma of the intestine",
    "HP:0200008": "Intestinal polyposis",
    "HP:0002027": "Abdominal pain",
    "HP:0001300": "Parkinsonism",
    "HP:0001288": "Gait disturbance",
    "HP:0002072": "Chorea",
    "HP:0002074": "Rigidity",
    "HP:0000726": "Dementia",
    "HP:0002527": "Falls",
    "HP:0002511": "Alzheimer disease",
    "HP:0001260": "Dysarthria",
    "HP:0001268": "Mental deterioration",
    "HP:0002099": "Asthma",
    "HP:0002110": "Bronchiectasis",
    "HP:0001696": "Situs inversus",
    "HP:0007760": "Sickle cell disease",
    "HP:0001744": "Splenomegaly",
    "HP:0100651": "Type 1 diabetes",
    "HP:0001250": "Seizures",
    "HP:0002960": "Autoimmunity",
    "HP:0001370": "Rheumatoid arthritis",
    "HP:0001369": "Arthritis",
    "HP:0001327": "Psychosis",
    "HP:0000729": "Autistic behavior",
    "HP:0001256": "Intellectual disability",
    "HP:0000718": "Aggressive behavior",
    "HP:0000752": "Hyperactivity",
    "HP:0100753": "Schizophrenia",
    "HP:0000708": "Behavioral abnormality",
    "HP:0000716": "Depression",
    "HP:0001249": "Intellectual disability",
    "HP:0001836": "Hemophilia",
    "HP:0003645": "Prolonged bleeding time",
    "HP:0000978": "Easy bruising",
    "HP:0001166": "Arachnodactyly",
    "HP:0000518": "Ectopia lentis",
    "HP:0001083": "Ectopia lentis",
    "HP:0002616": "Aortic root dilatation",
    "HP:0006746": "Neurofibromatosis",
    "HP:0009732": "Plexiform neurofibroma",
    "HP:0009733": "Hamartoma",
    "HP:0001263": "Global developmental delay",
    "HP:0003560": "Muscular dystrophy",
    "HP:0003701": "Proximal muscle weakness",
    "HP:0003236": "Elevated CK",
    "HP:0003557": "Increased variability in muscle fiber diameter",
    "HP:0001360": "Holoprosencephaly",
    "HP:0003128": "Lactic acidosis",
    "HP:0001410": "Hepatic failure",
    "HP:0010278": "Kayser-Fleischer ring",
    "HP:0001394": "Cirrhosis",
    "HP:0003429": "CNS demyelination",
    "HP:0001371": "Flexion contracture",
    "HP:0000786": "Primary amenorrhea",
    "HP:0004322": "Short stature",
    "HP:0000823": "Delayed puberty",
    "HP:0000744": "Azoospermia",
    "HP:0008669": "Abnormal spermatogenesis",
    "HP:0002725": "Systemic lupus erythematosus",
    "HP:0002608": "Celiac disease",
    "HP:0002013": "Vomiting",
    "HP:0002243": "Protein-losing enteropathy",
    "HP:0007354": "Amyotrophic lateral sclerosis",
    "HP:0001257": "Spasticity",
    "HP:0003326": "Myalgia",
    "HP:0002500": "Abnormal pyramidal signs",
    "HP:0002354": "Memory impairment",
    "HP:0004378": "Abnormality of the kidney",
    "HP:0001659": "Aortic regurgitation",
}

# Main dataset: 32 conditions with genes, HPO terms, evidence, PubMed IDs
HPO_DATA = {
    "Type 2 Diabetes": {
        "hpo_id": "HP:0005978", "omim": "125853", "category": "Metabolic",
        "description": "Impaired insulin secretion and/or insulin resistance causing hyperglycemia",
        "genes": {
            "TCF7L2":  {"hpo_terms": ["HP:0000819","HP:0001513","HP:0000855"], "evidence":"PCS", "pubmed":17293463, "effect":"risk"},
            "PPARG":   {"hpo_terms": ["HP:0000819","HP:0001513","HP:0003362"], "evidence":"PCS", "pubmed":10545951, "effect":"risk"},
            "KCNJ11":  {"hpo_terms": ["HP:0000819","HP:0000825","HP:0005978"], "evidence":"PCS", "pubmed":15781101, "effect":"risk"},
            "GCK":     {"hpo_terms": ["HP:0000819","HP:0003074","HP:0001508"], "evidence":"PCS", "pubmed":9039267,  "effect":"causal"},
            "HNF1A":   {"hpo_terms": ["HP:0000819","HP:0005978","HP:0001508"], "evidence":"PCS", "pubmed":17855067, "effect":"causal"},
            "SLC30A8": {"hpo_terms": ["HP:0000819","HP:0001513","HP:0000855"], "evidence":"PCS", "pubmed":17293464, "effect":"risk"},
            "CDKAL1":  {"hpo_terms": ["HP:0000819","HP:0001513","HP:0000855"], "evidence":"PCS", "pubmed":17460697, "effect":"risk"},
            "MTNR1B":  {"hpo_terms": ["HP:0000819","HP:0003074","HP:0000855"], "evidence":"PCS", "pubmed":19060910, "effect":"risk"},
            "ADCY5":   {"hpo_terms": ["HP:0000819","HP:0001508","HP:0001513"], "evidence":"PCS", "pubmed":20581827, "effect":"risk"},
        }
    },
    "Coronary Artery Disease": {
        "hpo_id": "HP:0001695", "omim": "608320", "category": "Cardiovascular",
        "description": "Atherosclerotic narrowing of coronary arteries leading to myocardial ischemia",
        "genes": {
            "PCSK9":  {"hpo_terms": ["HP:0003124","HP:0001695","HP:0001513"], "evidence":"PCS", "pubmed":15654334, "effect":"causal"},
            "APOE":   {"hpo_terms": ["HP:0003124","HP:0001695","HP:0002621"], "evidence":"PCS", "pubmed":21738483, "effect":"risk"},
            "LPA":    {"hpo_terms": ["HP:0003124","HP:0001695","HP:0001297"], "evidence":"PCS", "pubmed":19106084, "effect":"causal"},
            "LDLR":   {"hpo_terms": ["HP:0003124","HP:0001695","HP:0001659"], "evidence":"PCS", "pubmed":9390937,  "effect":"causal"},
            "APOB":   {"hpo_terms": ["HP:0003124","HP:0001695","HP:0002621"], "evidence":"PCS", "pubmed":11073740, "effect":"causal"},
            "CETP":   {"hpo_terms": ["HP:0003124","HP:0001695","HP:0000822"], "evidence":"PCS", "pubmed":18193043, "effect":"protective"},
            "HMGCR":  {"hpo_terms": ["HP:0003124","HP:0001695","HP:0001659"], "evidence":"PCS", "pubmed":17634449, "effect":"risk"},
            "SH2B3":  {"hpo_terms": ["HP:0001695","HP:0001873","HP:0001974"], "evidence":"PCS", "pubmed":19198612, "effect":"risk"},
        }
    },
    "Breast Cancer": {
        "hpo_id": "HP:0003002", "omim": "114480", "category": "Oncology",
        "description": "Malignant neoplasm of breast epithelial tissue",
        "genes": {
            "BRCA1": {"hpo_terms": ["HP:0003002","HP:0100615","HP:0002664"], "evidence":"PCS", "pubmed":7545954,  "effect":"causal"},
            "BRCA2": {"hpo_terms": ["HP:0003002","HP:0002671","HP:0002664"], "evidence":"PCS", "pubmed":8524414,  "effect":"causal"},
            "TP53":  {"hpo_terms": ["HP:0003002","HP:0002671","HP:0002664"], "evidence":"PCS", "pubmed":1565144,  "effect":"causal"},
            "PALB2": {"hpo_terms": ["HP:0003002","HP:0002671","HP:0002664"], "evidence":"PCS", "pubmed":17200671, "effect":"causal"},
            "ATM":   {"hpo_terms": ["HP:0003002","HP:0002664","HP:0001903"], "evidence":"PCS", "pubmed":16116422, "effect":"risk"},
            "CHEK2": {"hpo_terms": ["HP:0003002","HP:0002664","HP:0000818"], "evidence":"PCS", "pubmed":12610780, "effect":"risk"},
        }
    },
    "Colorectal Cancer": {
        "hpo_id": "HP:0003003", "omim": "114500", "category": "Oncology",
        "description": "Malignant neoplasm of the colon or rectum",
        "genes": {
            "APC":   {"hpo_terms": ["HP:0003003","HP:0002664","HP:0005275"], "evidence":"PCS", "pubmed":1651174,  "effect":"causal"},
            "MLH1":  {"hpo_terms": ["HP:0003003","HP:0002664","HP:0005275"], "evidence":"PCS", "pubmed":7951315,  "effect":"causal"},
            "MSH2":  {"hpo_terms": ["HP:0003003","HP:0002664","HP:0200008"], "evidence":"PCS", "pubmed":7951315,  "effect":"causal"},
            "MSH6":  {"hpo_terms": ["HP:0003003","HP:0002664","HP:0200008"], "evidence":"PCS", "pubmed":9590181,  "effect":"causal"},
            "KRAS":  {"hpo_terms": ["HP:0003003","HP:0002664","HP:0002027"], "evidence":"PCS", "pubmed":17671637, "effect":"causal"},
            "BRAF":  {"hpo_terms": ["HP:0003003","HP:0002664","HP:0002027"], "evidence":"PCS", "pubmed":15897572, "effect":"causal"},
        }
    },
    "Parkinson Disease": {
        "hpo_id": "HP:0001300", "omim": "168600", "category": "Neurological",
        "description": "Progressive neurodegenerative disorder of dopaminergic neurons",
        "genes": {
            "SNCA":  {"hpo_terms": ["HP:0001300","HP:0001288","HP:0002072"], "evidence":"PCS", "pubmed":9197268,  "effect":"causal"},
            "LRRK2": {"hpo_terms": ["HP:0001300","HP:0001288","HP:0000726"], "evidence":"PCS", "pubmed":15541309, "effect":"causal"},
            "PINK1": {"hpo_terms": ["HP:0001300","HP:0001288","HP:0002074"], "evidence":"PCS", "pubmed":15029634, "effect":"causal"},
            "PRKN":  {"hpo_terms": ["HP:0001300","HP:0001288","HP:0002072"], "evidence":"PCS", "pubmed":9560156,  "effect":"causal"},
            "GBA":   {"hpo_terms": ["HP:0001300","HP:0001288","HP:0002527"], "evidence":"PCS", "pubmed":19846850, "effect":"risk"},
            "DJ1":   {"hpo_terms": ["HP:0001300","HP:0001288","HP:0002527"], "evidence":"PCS", "pubmed":12564870, "effect":"causal"},
        }
    },
    "Alzheimer Disease": {
        "hpo_id": "HP:0002511", "omim": "104300", "category": "Neurological",
        "description": "Progressive neurodegenerative dementia with amyloid plaques and tau tangles",
        "genes": {
            "APOE":  {"hpo_terms": ["HP:0002511","HP:0001260","HP:0000726"], "evidence":"PCS", "pubmed":8346443,  "effect":"risk"},
            "APP":   {"hpo_terms": ["HP:0002511","HP:0001260","HP:0001268"], "evidence":"PCS", "pubmed":1671712,  "effect":"causal"},
            "PSEN1": {"hpo_terms": ["HP:0002511","HP:0001260","HP:0001268"], "evidence":"PCS", "pubmed":7596406,  "effect":"causal"},
            "PSEN2": {"hpo_terms": ["HP:0002511","HP:0001260","HP:0001268"], "evidence":"PCS", "pubmed":7584071,  "effect":"causal"},
            "CLU":   {"hpo_terms": ["HP:0002511","HP:0001260","HP:0000726"], "evidence":"PCS", "pubmed":19734902, "effect":"risk"},
            "BIN1":  {"hpo_terms": ["HP:0002511","HP:0001260","HP:0000726"], "evidence":"PCS", "pubmed":19734902, "effect":"risk"},
        }
    },
    "Cystic Fibrosis": {
        "hpo_id": "HP:0002099", "omim": "219700", "category": "Respiratory",
        "description": "Autosomal recessive multisystem disease affecting ion transport",
        "genes": {
            "CFTR":  {"hpo_terms": ["HP:0002099","HP:0002110","HP:0001696"], "evidence":"PCS", "pubmed":2475911, "effect":"causal"},
        }
    },
    "Sickle Cell Disease": {
        "hpo_id": "HP:0007760", "omim": "603903", "category": "Hematological",
        "description": "Autosomal recessive hemoglobinopathy causing vaso-occlusive crises",
        "genes": {
            "HBB":    {"hpo_terms": ["HP:0007760","HP:0001903","HP:0002621"], "evidence":"PCS", "pubmed":13054797, "effect":"causal"},
            "BCL11A": {"hpo_terms": ["HP:0007760","HP:0001903","HP:0001744"], "evidence":"PCS", "pubmed":19196962, "effect":"modifier"},
        }
    },
    "Familial Hypercholesterolemia": {
        "hpo_id": "HP:0003124", "omim": "143890", "category": "Metabolic",
        "description": "Autosomal dominant disorder of LDL cholesterol metabolism",
        "genes": {
            "LDLR":   {"hpo_terms": ["HP:0003124","HP:0001695","HP:0001659"], "evidence":"PCS", "pubmed":9390937,  "effect":"causal"},
            "APOB":   {"hpo_terms": ["HP:0003124","HP:0001695","HP:0002621"], "evidence":"PCS", "pubmed":11073740, "effect":"causal"},
            "PCSK9":  {"hpo_terms": ["HP:0003124","HP:0001695","HP:0001513"], "evidence":"PCS", "pubmed":15654334, "effect":"causal"},
            "LDLRAP1":{"hpo_terms": ["HP:0003124","HP:0001695","HP:0001513"], "evidence":"PCS", "pubmed":12471198, "effect":"causal"},
        }
    },
    "Hypertension": {
        "hpo_id": "HP:0000822", "omim": "145500", "category": "Cardiovascular",
        "description": "Persistent elevation of arterial blood pressure",
        "genes": {
            "ACE":   {"hpo_terms": ["HP:0000822","HP:0001659","HP:0004378"], "evidence":"PCS", "pubmed":8723067,  "effect":"risk"},
            "AGT":   {"hpo_terms": ["HP:0000822","HP:0001695","HP:0004378"], "evidence":"PCS", "pubmed":8423673,  "effect":"risk"},
            "AGTR1": {"hpo_terms": ["HP:0000822","HP:0001695","HP:0001744"], "evidence":"PCS", "pubmed":9570186,  "effect":"risk"},
            "NOS3":  {"hpo_terms": ["HP:0000822","HP:0001695","HP:0000822"], "evidence":"PCS", "pubmed":9633806,  "effect":"risk"},
            "ADD1":  {"hpo_terms": ["HP:0000822","HP:0004378","HP:0001695"], "evidence":"PCS", "pubmed":9633807,  "effect":"risk"},
        }
    },
    "Type 1 Diabetes": {
        "hpo_id": "HP:0100651", "omim": "222100", "category": "Metabolic",
        "description": "Autoimmune destruction of pancreatic beta cells causing absolute insulin deficiency",
        "genes": {
            "HLA-DRB1":{"hpo_terms":["HP:0100651","HP:0000819","HP:0001250"], "evidence":"PCS", "pubmed":3879741,  "effect":"causal"},
            "INS":    {"hpo_terms": ["HP:0100651","HP:0000819","HP:0001508"], "evidence":"PCS", "pubmed":9753729,  "effect":"causal"},
            "PTPN22": {"hpo_terms": ["HP:0100651","HP:0002960","HP:0000819"], "evidence":"PCS", "pubmed":14647274, "effect":"risk"},
            "CTLA4":  {"hpo_terms": ["HP:0100651","HP:0002960","HP:0000819"], "evidence":"PCS", "pubmed":10619861, "effect":"risk"},
            "IL2RA":  {"hpo_terms": ["HP:0100651","HP:0002960","HP:0001903"], "evidence":"PCS", "pubmed":17554260, "effect":"risk"},
        }
    },
    "Rheumatoid Arthritis": {
        "hpo_id": "HP:0001370", "omim": "180300", "category": "Autoimmune",
        "description": "Systemic autoimmune disease causing symmetric inflammatory polyarthritis",
        "genes": {
            "HLA-DRB1":{"hpo_terms":["HP:0001370","HP:0001369","HP:0002960"], "evidence":"PCS", "pubmed":8893090,  "effect":"causal"},
            "PTPN22": {"hpo_terms": ["HP:0001370","HP:0002960","HP:0001369"], "evidence":"PCS", "pubmed":15208781, "effect":"risk"},
            "STAT4":  {"hpo_terms": ["HP:0001370","HP:0002960","HP:0001369"], "evidence":"PCS", "pubmed":17804836, "effect":"risk"},
            "TNFAIP3":{"hpo_terms": ["HP:0001370","HP:0002960","HP:0001369"], "evidence":"PCS", "pubmed":17804836, "effect":"risk"},
            "PADI4":  {"hpo_terms": ["HP:0001370","HP:0001369","HP:0002960"], "evidence":"PCS", "pubmed":12754510, "effect":"risk"},
        }
    },
    "Autism Spectrum Disorder": {
        "hpo_id": "HP:0000729", "omim": "209850", "category": "Neurological",
        "description": "Neurodevelopmental condition affecting social communication and behavior",
        "genes": {
            "SHANK3": {"hpo_terms": ["HP:0000729","HP:0001256","HP:0000718"], "evidence":"PCS", "pubmed":17558391, "effect":"causal"},
            "NLGN3":  {"hpo_terms": ["HP:0000729","HP:0001256","HP:0000752"], "evidence":"PCS", "pubmed":12669065, "effect":"causal"},
            "NRXN1":  {"hpo_terms": ["HP:0000729","HP:0001256","HP:0000718"], "evidence":"PCS", "pubmed":18179900, "effect":"causal"},
            "CHD8":   {"hpo_terms": ["HP:0000729","HP:0001256","HP:0001250"], "evidence":"PCS", "pubmed":24267886, "effect":"causal"},
            "PTEN":   {"hpo_terms": ["HP:0000729","HP:0001256","HP:0001250"], "evidence":"PCS", "pubmed":15805158, "effect":"causal"},
            "FMR1":   {"hpo_terms": ["HP:0000729","HP:0001256","HP:0000752"], "evidence":"PCS", "pubmed":2180487,  "effect":"causal"},
        }
    },
    "Schizophrenia": {
        "hpo_id": "HP:0100753", "omim": "181500", "category": "Neurological",
        "description": "Severe psychiatric disorder with psychosis, cognitive deficits, and negative symptoms",
        "genes": {
            "DISC1":  {"hpo_terms": ["HP:0100753","HP:0000708","HP:0000716"], "evidence":"PCS", "pubmed":11062471, "effect":"risk"},
            "DTNBP1": {"hpo_terms": ["HP:0100753","HP:0000708","HP:0001249"], "evidence":"PCS", "pubmed":12098051, "effect":"risk"},
            "NRG1":   {"hpo_terms": ["HP:0100753","HP:0000708","HP:0000716"], "evidence":"PCS", "pubmed":12089511, "effect":"risk"},
            "COMT":   {"hpo_terms": ["HP:0100753","HP:0000708","HP:0001327"], "evidence":"PCS", "pubmed":11385583, "effect":"risk"},
            "22q11":  {"hpo_terms": ["HP:0100753","HP:0000708","HP:0001250"], "evidence":"PCS", "pubmed":8290963,  "effect":"risk"},
        }
    },
    "Hemophilia A": {
        "hpo_id": "HP:0001836", "omim": "306700", "category": "Hematological",
        "description": "X-linked coagulation disorder due to Factor VIII deficiency",
        "genes": {
            "F8":    {"hpo_terms": ["HP:0001836","HP:0003645","HP:0000978"], "evidence":"PCS", "pubmed":2834511, "effect":"causal"},
        }
    },
    "Hemophilia B": {
        "hpo_id": "HP:0001836", "omim": "306900", "category": "Hematological",
        "description": "X-linked coagulation disorder due to Factor IX deficiency",
        "genes": {
            "F9":    {"hpo_terms": ["HP:0001836","HP:0003645","HP:0000978"], "evidence":"PCS", "pubmed":2433168, "effect":"causal"},
        }
    },
    "Marfan Syndrome": {
        "hpo_id": "HP:0001166", "omim": "154700", "category": "Connective Tissue",
        "description": "Autosomal dominant connective tissue disorder with cardiovascular, ocular, skeletal features",
        "genes": {
            "FBN1":   {"hpo_terms": ["HP:0001166","HP:0000518","HP:0001083"], "evidence":"PCS", "pubmed":1303228,  "effect":"causal"},
            "TGFBR1": {"hpo_terms": ["HP:0001166","HP:0001083","HP:0002616"], "evidence":"PCS", "pubmed":15731757, "effect":"causal"},
            "TGFBR2": {"hpo_terms": ["HP:0001166","HP:0001083","HP:0002616"], "evidence":"PCS", "pubmed":15731757, "effect":"causal"},
        }
    },
    "Huntington Disease": {
        "hpo_id": "HP:0002072", "omim": "143100", "category": "Neurological",
        "description": "Autosomal dominant neurodegenerative disorder with CAG repeat expansion",
        "genes": {
            "HTT":   {"hpo_terms": ["HP:0002072","HP:0000726","HP:0001300"], "evidence":"PCS", "pubmed":8458085, "effect":"causal"},
        }
    },
    "Duchenne Muscular Dystrophy": {
        "hpo_id": "HP:0003560", "omim": "310200", "category": "Muscular",
        "description": "X-linked progressive muscle degeneration due to dystrophin absence",
        "genes": {
            "DMD":   {"hpo_terms": ["HP:0003560","HP:0003701","HP:0003236"], "evidence":"PCS", "pubmed":2885835, "effect":"causal"},
            "SGCA":  {"hpo_terms": ["HP:0003560","HP:0003701","HP:0003557"], "evidence":"PCS", "pubmed":9285790, "effect":"causal"},
        }
    },
    "Phenylketonuria": {
        "hpo_id": "HP:0001360", "omim": "261600", "category": "Metabolic",
        "description": "Autosomal recessive disorder of phenylalanine metabolism",
        "genes": {
            "PAH":   {"hpo_terms": ["HP:0001360","HP:0001256","HP:0001903"], "evidence":"PCS", "pubmed":2897727, "effect":"causal"},
            "GCH1":  {"hpo_terms": ["HP:0001360","HP:0001256","HP:0003128"], "evidence":"PCS", "pubmed":7874130, "effect":"causal"},
        }
    },
    "Wilson Disease": {
        "hpo_id": "HP:0001410", "omim": "277900", "category": "Metabolic",
        "description": "Autosomal recessive copper metabolism disorder",
        "genes": {
            "ATP7B": {"hpo_terms": ["HP:0001410","HP:0001300","HP:0010278"], "evidence":"PCS", "pubmed":8298641, "effect":"causal"},
        }
    },
    "Gaucher Disease": {
        "hpo_id": "HP:0001394", "omim": "230800", "category": "Metabolic",
        "description": "Lysosomal storage disorder due to glucocerebrosidase deficiency",
        "genes": {
            "GBA":   {"hpo_terms": ["HP:0001394","HP:0001903","HP:0001744"], "evidence":"PCS", "pubmed":5646079, "effect":"causal"},
        }
    },
    "Fabry Disease": {
        "hpo_id": "HP:0003429", "omim": "301500", "category": "Metabolic",
        "description": "X-linked lysosomal storage disorder due to alpha-galactosidase A deficiency",
        "genes": {
            "GLA":   {"hpo_terms": ["HP:0003429","HP:0000822","HP:0001371"], "evidence":"PCS", "pubmed":4859836, "effect":"causal"},
        }
    },
    "Neurofibromatosis 1": {
        "hpo_id": "HP:0006746", "omim": "162200", "category": "Neurological",
        "description": "Autosomal dominant disorder causing neurofibromas and cafe-au-lait spots",
        "genes": {
            "NF1":   {"hpo_terms": ["HP:0006746","HP:0001250","HP:0009732"], "evidence":"PCS", "pubmed":2675589, "effect":"causal"},
        }
    },
    "Tuberous Sclerosis": {
        "hpo_id": "HP:0001250", "omim": "191100", "category": "Neurological",
        "description": "Autosomal dominant disorder causing hamartomas in multiple organ systems",
        "genes": {
            "TSC1":  {"hpo_terms": ["HP:0001250","HP:0001263","HP:0009733"], "evidence":"PCS", "pubmed":9291101, "effect":"causal"},
            "TSC2":  {"hpo_terms": ["HP:0001250","HP:0001263","HP:0009733"], "evidence":"PCS", "pubmed":8162069, "effect":"causal"},
        }
    },
    "Fragile X Syndrome": {
        "hpo_id": "HP:0000752", "omim": "300624", "category": "Neurological",
        "description": "X-linked intellectual disability due to FMR1 gene silencing",
        "genes": {
            "FMR1":  {"hpo_terms": ["HP:0000752","HP:0001249","HP:0000729"], "evidence":"PCS", "pubmed":2180487, "effect":"causal"},
        }
    },
    "Down Syndrome": {
        "hpo_id": "HP:0001249", "omim": "190685", "category": "Chromosomal",
        "description": "Trisomy 21 causing intellectual disability and characteristic phenotype",
        "genes": {
            "DYRK1A":{"hpo_terms": ["HP:0001249","HP:0001250","HP:0000486"], "evidence":"PCS", "pubmed":10079247, "effect":"causal"},
            "APP":   {"hpo_terms": ["HP:0001249","HP:0001250","HP:0002511"], "evidence":"PCS", "pubmed":8023167,  "effect":"risk"},
            "SOD1":  {"hpo_terms": ["HP:0001249","HP:0001250","HP:0000486"], "evidence":"PCS", "pubmed":1584799,  "effect":"risk"},
        }
    },
    "Turner Syndrome": {
        "hpo_id": "HP:0000786", "omim": "163950", "category": "Chromosomal",
        "description": "45,X monosomy causing short stature, gonadal dysgenesis in females",
        "genes": {
            "SHOX":  {"hpo_terms": ["HP:0000786","HP:0004322","HP:0000823"], "evidence":"PCS", "pubmed":9590291, "effect":"causal"},
        }
    },
    "Klinefelter Syndrome": {
        "hpo_id": "HP:0000744", "omim": "312865", "category": "Chromosomal",
        "description": "47,XXY karyotype causing hypogonadism and infertility in males",
        "genes": {
            "AR":    {"hpo_terms": ["HP:0000744","HP:0000823","HP:0008669"], "evidence":"PCS", "pubmed":9020839, "effect":"modifier"},
        }
    },
    "Lupus (SLE)": {
        "hpo_id": "HP:0002725", "omim": "152700", "category": "Autoimmune",
        "description": "Systemic autoimmune disease with multi-organ involvement and anti-nuclear antibodies",
        "genes": {
            "TREX1": {"hpo_terms": ["HP:0002725","HP:0002960","HP:0000822"], "evidence":"PCS", "pubmed":17660530, "effect":"causal"},
            "IRF5":  {"hpo_terms": ["HP:0002725","HP:0002960","HP:0001903"], "evidence":"PCS", "pubmed":15716159, "effect":"risk"},
            "BLK":   {"hpo_terms": ["HP:0002725","HP:0002960","HP:0001903"], "evidence":"PCS", "pubmed":18204447, "effect":"risk"},
            "PTPN22":{"hpo_terms": ["HP:0002725","HP:0002960","HP:0001903"], "evidence":"PCS", "pubmed":14647274, "effect":"risk"},
            "STAT4": {"hpo_terms": ["HP:0002725","HP:0002960","HP:0001903"], "evidence":"PCS", "pubmed":18204448, "effect":"risk"},
        }
    },
    "Celiac Disease": {
        "hpo_id": "HP:0002608", "omim": "212750", "category": "Autoimmune",
        "description": "Gluten-triggered autoimmune enteropathy in genetically susceptible individuals",
        "genes": {
            "HLA-DQ2":{"hpo_terms":["HP:0002608","HP:0002013","HP:0001903"], "evidence":"PCS", "pubmed":2002548,  "effect":"causal"},
            "HLA-DQ8":{"hpo_terms":["HP:0002608","HP:0002013","HP:0001903"], "evidence":"PCS", "pubmed":2002548,  "effect":"causal"},
            "MYO9B":  {"hpo_terms": ["HP:0002608","HP:0002013","HP:0002243"], "evidence":"PCS", "pubmed":15895351, "effect":"risk"},
        }
    },
    "Amyotrophic Lateral Sclerosis": {
        "hpo_id": "HP:0007354", "omim": "105400", "category": "Neurological",
        "description": "Progressive motor neuron disease causing muscle weakness and paralysis",
        "genes": {
            "SOD1":   {"hpo_terms": ["HP:0007354","HP:0001257","HP:0003326"], "evidence":"PCS", "pubmed":8271574,  "effect":"causal"},
            "TDP-43": {"hpo_terms": ["HP:0007354","HP:0001257","HP:0002500"], "evidence":"PCS", "pubmed":17215369, "effect":"causal"},
            "FUS":    {"hpo_terms": ["HP:0007354","HP:0001257","HP:0003326"], "evidence":"PCS", "pubmed":19251628, "effect":"causal"},
            "C9orf72":{"hpo_terms": ["HP:0007354","HP:0002354","HP:0003326"], "evidence":"PCS", "pubmed":21944778, "effect":"causal"},
        }
    },
}

# Disease category colour map (used in network visualisation)
CATEGORY_COLORS = {
    "Metabolic":       "#3B82F6",  # blue
    "Cardiovascular":  "#EF4444",  # red
    "Oncology":        "#8B5CF6",  # purple
    "Neurological":    "#F97316",  # orange
    "Autoimmune":      "#10B981",  # green
    "Hematological":   "#EC4899",  # pink
    "Respiratory":     "#06B6D4",  # cyan
    "Connective Tissue":"#84CC16", # lime
    "Muscular":        "#F59E0B",  # amber
    "Chromosomal":     "#6366F1",  # indigo
}

EFFECT_COLORS = {
    "causal":     "#DC2626",
    "risk":       "#F97316",
    "protective": "#22C55E",
    "modifier":   "#6366F1",
}
