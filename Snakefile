import os
import csv
import json
import pandas as pd 
import ast
import glob

configfile: "config.yaml"

# -------------------------
# Helpers
# -------------------------

TASK_ROOT       = config["task_root"]
SCRIPT_ROOT     = config["script_root"]
SEQS_CSV        = config["input"]["sequences"].format(task_root=TASK_ROOT)
SMILES_CSV      = config["input"]["smiles"].format(task_root=TASK_ROOT)
INIT_DIRNAME    = config["initial"]["out_dir"]
FINE_DIRNAME    = config['fine']['out_dir']
METHODS         = list(config["initial"]["methods"])
P_SUMMARY       = config["initial"]["per_protein_summary"]
P_SELECTED      = config["initial"]["per_protein_selected"]
GLOBAL_SUMMARY  = config.get("global_summary", None)

TOOLS = config["tools"]

def protein_dir(p):          return os.path.join(TASK_ROOT, p)
def init_dir(p):             return os.path.join(protein_dir(p), INIT_DIRNAME)
def fine_dir(p):             return os.path.join(protein_dir(p), FINE_DIRNAME)
def method_dir(p, m):        return os.path.join(init_dir(p), m)
def input_csv(p):            return os.path.join(init_dir(p), "inputs", "finish.token")
def method_fail_csv(p, m):   return os.path.join(method_dir(p, m), "failed_smiles.csv")
def p_summary_csv(p):        return os.path.join(init_dir(p), P_SUMMARY)
def p_selected_csv(p):       return os.path.join(init_dir(p), P_SELECTED)

# If proteins aren’t listed in the config, we’ll let Snakemake discover them
# by creating input.csv for each protein and reading protein names from sequences.csv.
PROTEINS = config.get("proteins", None)

def get_proteins():
    """Generate proteins list if not provided in config"""
    if PROTEINS is None or len(PROTEINS) == 0:
        # Your script logic here to generate proteins
        # For example:
        proteins = pd.read_csv(SEQS_CSV)['name']
        # Add your logic to populate the proteins list
        # e.g., read from a file, scan a directory, etc.
        return proteins
    return PROTEINS

def init_targets(wildcards):
    targets = []
    for protein in get_proteins():
        ckpt = checkpoints.make_input_csv.get(protein=protein)
        init_in_dir = ckpt.output[0]  # .../init_screening/inputs/
        init_out_dir = [os.path.join(TASK_ROOT, protein, INIT_DIRNAME, "GraphDTA"),
                        os.path.join(TASK_ROOT, protein, INIT_DIRNAME, "ColdDTA"),
                        os.path.join(TASK_ROOT, protein, INIT_DIRNAME, "HMSA"),
                        # os.path.join(TASK_ROOT, protein, INIT_DIRNAME, "EviDTI"),  # DISABLED
                        os.path.join(TASK_ROOT, protein, INIT_DIRNAME, "DrugLAMP"),
                        os.path.join(TASK_ROOT, protein, INIT_DIRNAME, "ConPLex")]

        # IMPORTANT: capture the FULL stem used by your rule, e.g. "input_part_0"
        parts = glob_wildcards(os.path.join(init_in_dir, "input_{i}.csv")).i

        # Map those stems to the .done markers under the OUTPUT dir
        for dir in init_out_dir:
            targets.extend([os.path.join(dir, f"prediction_{i}.csv") for i in parts])
    return targets

def get_box_str(protein_name, seq_file, name_col='name', box_col='docking box'):
    df = pd.read_csv(seq_file)
    cell = df.loc[df[name_col] == protein_name, box_col].iloc[0]
    # If the CSV stores the box as a string like "[25,33,57,22,22,22]" parse it:
    if isinstance(cell, str):
        vals = ast.literal_eval(cell)
    else:
        vals = list(cell)
    center = vals[:3]
    size   = vals[3:]
    return " ".join(map(str, center)), " ".join(map(str, size))


def vina_targets(_wc):
    targets = []
    for protein in get_proteins():
        ckpt = checkpoints.split_csv.get(protein=protein)
        vin_in_dir = ckpt.output[0]  # .../Vina/input/
        vin_out_dir = os.path.join(TASK_ROOT, protein, FINE_DIRNAME, "Vina", "output")

        # IMPORTANT: capture the FULL stem used by your rule, e.g. "input_part_0"
        parts = glob_wildcards(os.path.join(vin_in_dir, "{part}.csv")).part

        # Map those stems to the .done markers under the OUTPUT dir
        targets.extend([os.path.join(vin_out_dir, f"{p}.done") for p in parts])
    return targets

def get_boltz2_inputs(wildcards):
    # Force evaluation of the checkpoint
    checkpoint_output = checkpoints.write_boltz2_input.get(protein=wildcards.protein).output[0]

    ids = glob_wildcards(os.path.join(checkpoint_output, "{i}.yaml")).i

    return ids

def get_boltz2_batch_jobs(protein, batch_id, n_batches=40):
    """Get all job indices for a specific batch"""
    smiles_file = checkpoints.compile_per_protein_summary.get(protein=protein).output.selected
    n_smiles = len(pd.read_csv(smiles_file))

    # Calculate jobs per batch (distribute evenly)
    jobs = list(range(n_smiles))
    batch_size = (n_smiles + n_batches - 1) // n_batches  # ceiling division

    # Get jobs for this batch
    start_idx = batch_id * batch_size
    end_idx = min((batch_id + 1) * batch_size, n_smiles)

    return list(range(start_idx, end_idx))

def boltz2_targets(wildcards):
	"""Create targets for 40 batches instead of individual jobs"""
	targets = []
	n_batches = 40
	for p in get_proteins():
		smiles_file = checkpoints.compile_per_protein_summary.get(protein=p).output.selected
		n_smiles = len(pd.read_csv(smiles_file))
		# Create one target per batch instead of per SMILES
		for batch_id in range(n_batches):
			batch_size = (n_smiles + n_batches - 1) // n_batches
			if batch_id * batch_size < n_smiles:  # Only create batches that have jobs
				token_file = os.path.join(TASK_ROOT, p, FINE_DIRNAME, f"Boltz2/output/token/batch_{batch_id}.done")
				targets.append(token_file)
	return targets

def get_boltz2_tar_files(wildcards):
	"""Get all tar files for a specific protein"""
	tar_files = []
	n_batches = 40
	smiles_file = checkpoints.compile_per_protein_summary.get(protein=wildcards.protein).output.selected
	n_smiles = len(pd.read_csv(smiles_file))

	for batch_id in range(n_batches):
		batch_size = (n_smiles + n_batches - 1) // n_batches
		if batch_id * batch_size < n_smiles:  # Only include batches that have jobs
			tar_file = os.path.join(TASK_ROOT, wildcards.protein, FINE_DIRNAME, f"Boltz2/output/batch_{batch_id}.tar.gz")
			tar_files.append(tar_file)
	return tar_files

def af3_targets(wildcards):
	targets = []
	for p in get_proteins():
		smiles_file = checkpoints.compile_per_protein_summary.get(protein=p).output.selected
		n_smiles = len(pd.read_csv(smiles_file))
		for i in range(n_smiles):
			token_file = os.path.join(TASK_ROOT, p, FINE_DIRNAME, f"AF3/output/token/{i}.done")
			targets.append(token_file)
	return targets

def diffdock_targets(wildcards):
	targets = []
	for p in get_proteins():
		ckpt = checkpoints.split_csv.get(protein=p)
		vin_in_dir = ckpt.output[0]
		parts = glob_wildcards(os.path.join(vin_in_dir, "{part}.csv")).part
		targets.extend([os.path.join(TASK_ROOT, p, FINE_DIRNAME, "PBSA/DiffDock/output", f"{i}.done") for i in parts])
	return targets

def md_targets(wildcards):
	targets = []
	for p in get_proteins():
		smiles_file = checkpoints.compile_per_protein_summary.get(protein=p).output.selected
		n_smiles = len(pd.read_csv(smiles_file))
		for i in range(n_smiles):
			token_file = os.path.join(TASK_ROOT, p, FINE_DIRNAME, f"PBSA/PBSA/MD/{p}_{i}/token.done")
			targets.append(token_file)
	return targets

def pbsa_targets(wildcards):
	targets = []
	for p in get_proteins():
		smiles_file = checkpoints.compile_per_protein_summary.get(protein=p).output.selected
		n_smiles = len(pd.read_csv(smiles_file))
		for i in range(n_smiles):
			token_file = os.path.join(TASK_ROOT, p, FINE_DIRNAME, f"PBSA/PBSA/PBSA/{p}_{i}/token.done")
			targets.append(token_file)
	return targets

rule all:
    input:
        init_targets, 
        lambda wildcards: expand(
            TASK_ROOT+"/{protein}/"+INIT_DIRNAME+"/"+P_SELECTED, 
            protein=get_proteins(),
        ),
        lambda wildcards: expand(
            TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/AF3/prefold/{protein}.json',
            protein=get_proteins(),
        ),
        lambda wildcards: expand(
            "/home/ubuntu/{protein}/boltz2_tmp/boltz_results_{protein}/predictions/{protein}/confidence_{protein}_model_0.json",
            protein=get_proteins(),
        ),

        lambda wildcards: expand(
		TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/AF3/prefold/prefold.done',
		protein=get_proteins(),
        ),
        vina_targets, 
        lambda wildcards: expand(
               "/home/ubuntu/{protein}/boltz2_tmp/boltz_input.done", protein=get_proteins()
        ),
        lambda wildcards: expand(
             TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/AF3/input/af3_input.done', protein=get_proteins()
        ),
        # boltz2_targets,
        af3_targets,
        # diffdock_targets,
        # md_targets,
        # pbsa_targets,
        # lambda wildcards: expand(
        #         TASK_ROOT + "/{protein}/" + FINE_DIRNAME + "/Vina/results.csv", 
        #         protein=get_proteins()
        # ),
        lambda wildcards: expand(
              TASK_ROOT + "/{protein}/" + FINE_DIRNAME + "/AF3/summary.csv",
              protein=get_proteins()
        ),        
        # lambda wildcards: expand(
        #         TASK_ROOT + "/{protein}/" + FINE_DIRNAME + "/Boltz2/summary.csv",
        #         protein=get_proteins()
        # ),
        # lambda wildcards: expand(
        #         TASK_ROOT + "/{protein}/" + FINE_DIRNAME + "/PBSA/summary.csv",
        #         protein=get_proteins()
        # )


checkpoint make_input_csv:
    output:
        directory(TASK_ROOT+"/{protein}/"+INIT_DIRNAME+"/inputs/")
    input:
        seqs=SEQS_CSV,
        smiles=SMILES_CSV
    params:
        outdir=lambda wildcards: os.path.join(init_dir(wildcards.protein), "inputs"),
        make_input=TOOLS['make_input'],
        token=lambda wildcards: os.path.join(init_dir(wildcards.protein), "inputs/finish.token")
    shell:
        r"""
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate HMSA
        mkdir -p "{params.outdir}"
        python "{params.make_input}" \
          --sequences "{input.seqs}" \
          --smiles "{input.smiles}" \
          --protein "{wildcards.protein}" \
          --outdir "{params.outdir}"
        touch {params.token}
        """


# -------------------------
# Method runners (GraphDTA / HMSA / ColdDTA)
# Each method writes its own folder.
# Adjust CLIs as needed.
# -------------------------

rule run_graphdta:
    input: token=lambda wildcards: os.path.join(init_dir(wildcards.protein), "inputs/finish.token") 
    output:
        csv=TASK_ROOT+"/{protein}/" + INIT_DIRNAME + "/GraphDTA/prediction_{i}.csv",
        failed=TASK_ROOT+"/{protein}/" + INIT_DIRNAME + "/GraphDTA/failed_smiles_{i}.csv"
    params: outdir=lambda wc: method_dir(wc.protein, "GraphDTA"), exe=TOOLS["graphdta"], csv=lambda wc: os.path.join(init_dir(wc.protein), "inputs", f"input_{wc.i}.csv")
    shell:
        r"""
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate graphdta 
        mkdir -p "{params.outdir}"
        python "{params.exe}" \
            "{params.csv}" --model "/home/ubuntu/screening_workflow/algos/GraphDTA/model_para/model_GINConvNet_kiba.pt" \
            --model-type "gin" --batch-size "64" --output "{output.csv}" \
            --failed-smiles "{output.failed}" --smiles-column "SMILES"
        """

rule run_hmsa:
    input: token=lambda wildcards: os.path.join(init_dir(wildcards.protein), "inputs/finish.token")
    output:
        csv=TASK_ROOT+"/{protein}/" + INIT_DIRNAME + "/HMSA/prediction_{i}.csv"
    params: outdir=lambda wc: method_dir(wc.protein, "HMSA"), exe=TOOLS["hmsa"], csv=lambda wc: os.path.join(init_dir(wc.protein), "inputs", f"input_{wc.i}.csv")
    shell:
        r"""
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate HMSA
        mkdir -p "{params.outdir}"
        python "{params.exe}" --test_path "{params.csv}" \
            --preds_path "{output.csv}" --checkpoint_paths "/home/ubuntu/screening_workflow/algos/HMSA-DTI/model_para/model.pt" \
            --smiles_columns "SMILES"
        """

rule run_colddta:
    input:
        token=lambda wildcards: os.path.join(init_dir(wildcards.protein), "inputs/finish.token")
    output:
        csv=TASK_ROOT+"/{protein}/" + INIT_DIRNAME + "/ColdDTA/prediction_{i}.csv",
        failed=TASK_ROOT+"/{protein}/" + INIT_DIRNAME + "/ColdDTA/failed_smiles_{i}.csv"
    params:
        outdir=lambda wc: method_dir(wc.protein, "ColdDTA"), 
        exe=TOOLS["colddta"],
        csv=lambda wc: os.path.join(init_dir(wc.protein), "inputs", f"input_{wc.i}.csv")
    shell:
        r"""
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate cold
        mkdir -p "{params.outdir}"
        python "{params.exe}" --input "{params.csv}" --output "{output.csv}" \
            --failed-smiles "{output.failed}" --batch-size 64 --target-length 1000 \
            --checkpoint /home/ubuntu/screening_workflow/algos/coldDTA/model/epoch1297test_loss0.1798.pt
        """

# -------------------------
# New Method Runners (EviDTI / DrugLAMP / ConPLex)
# -------------------------

# -------------------------
# EviDTI - DISABLED
# Uncomment the rules below to re-enable EviDTI in the workflow
# -------------------------

# checkpoint extract_evidti_features:
#     input:
#         token=lambda wildcards: os.path.join(init_dir(wildcards.protein), "inputs/finish.token")
#     output:
#         directory(TASK_ROOT+"/{protein}/"+INIT_DIRNAME+"/EviDTI/features/")
#     params:
#         input_dir=lambda wc: os.path.join(init_dir(wc.protein), "inputs"),
#         feature_dir=lambda wc: os.path.join(init_dir(wc.protein), "EviDTI/features"),
#         protein_name=lambda wc: wc.protein,
#         extract_script=TOOLS["evidti_extract_features"]
#     shell:
#         r"""
#         source $(conda info --base)/etc/profile.d/conda.sh
#         conda activate evidti
#         export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
#         mkdir -p "{params.feature_dir}"
#         python "{params.extract_script}" \
#             --input-dir "{params.input_dir}" \
#             --output-dir "{params.feature_dir}" \
#             --protein-name "{params.protein_name}"
#         """

# rule run_evidti:
#     input:
#         token=lambda wildcards: os.path.join(init_dir(wildcards.protein), "inputs/finish.token"),
#         features=lambda wildcards: os.path.join(init_dir(wildcards.protein), "EviDTI/features/")
#     output:
#         csv=TASK_ROOT+"/{protein}/"+INIT_DIRNAME+"/EviDTI/prediction_{i}.csv",
#         failed=TASK_ROOT+"/{protein}/"+INIT_DIRNAME+"/EviDTI/failed_smiles_{i}.csv"
#     params:
#         outdir=lambda wc: method_dir(wc.protein, "EviDTI"),
#         exe=TOOLS["evidti_wrapper"],
#         csv=lambda wc: os.path.join(init_dir(wc.protein), "inputs", f"input_{wc.i}.csv"),
#         feature_dir=lambda wc: os.path.join(init_dir(wc.protein), "EviDTI/features"),
#         model=TOOLS["evidti_model"],
#         protein_name=lambda wc: wc.protein
#     shell:
#         r"""
#         source $(conda info --base)/etc/profile.d/conda.sh
#         conda activate evidti
#         export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
#         mkdir -p "{params.outdir}"
#         python "{params.exe}" \
#             --input "{params.csv}" \
#             --output "{output.csv}" \
#             --failed-smiles "{output.failed}" \
#             --feature-dir "{params.feature_dir}" \
#             --model-path "{params.model}" \
#             --protein-name "{params.protein_name}" \
#             --batch-size 32
#         """

rule run_druglamp:
    input:
        token=lambda wildcards: os.path.join(init_dir(wildcards.protein), "inputs/finish.token")
    output:
        csv=TASK_ROOT+"/{protein}/"+INIT_DIRNAME+"/DrugLAMP/prediction_{i}.csv",
        failed=TASK_ROOT+"/{protein}/"+INIT_DIRNAME+"/DrugLAMP/failed_smiles_{i}.csv"
    params:
        outdir=lambda wc: method_dir(wc.protein, "DrugLAMP"),
        exe=TOOLS["druglamp_wrapper"],
        csv=lambda wc: os.path.join(init_dir(wc.protein), "inputs", f"input_{wc.i}.csv"),
        checkpoint=TOOLS["druglamp_checkpoint"]
    shell:
        r"""
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate drug_lamp
        export MKL_THREADING_LAYER=GNU
        mkdir -p "{params.outdir}"
        python "{params.exe}" \
            --input "{params.csv}" \
            --output "{output.csv}" \
            --failed-smiles "{output.failed}" \
            --checkpoint "{params.checkpoint}" \
            --model DrugLAMP \
            --n-layer 30 \
            --device cuda \
            --batch-size 16
        """

rule run_conplex:
    input:
        token=lambda wildcards: os.path.join(init_dir(wildcards.protein), "inputs/finish.token")
    output:
        csv=TASK_ROOT+"/{protein}/"+INIT_DIRNAME+"/ConPLex/prediction_{i}.csv",
        failed=TASK_ROOT+"/{protein}/"+INIT_DIRNAME+"/ConPLex/failed_smiles_{i}.csv"
    params:
        outdir=lambda wc: method_dir(wc.protein, "ConPLex"),
        exe=TOOLS["conplex_wrapper"],
        csv=lambda wc: os.path.join(init_dir(wc.protein), "inputs", f"input_{wc.i}.csv"),
        model=TOOLS["conplex_model"],
        protein_name=lambda wc: wc.protein
    shell:
        r"""
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate conplex-dti
        mkdir -p "{params.outdir}"
        python "{params.exe}" \
            --input "{params.csv}" \
            --output "{output.csv}" \
            --failed-smiles "{output.failed}" \
            --model-path "{params.model}" \
            --protein-name "{params.protein_name}" \
            --device 0 \
            --batch-size 128
        """

# -------------------------
# Per-protein compile summary from method outputs
# -------------------------
checkpoint compile_per_protein_summary:
    input:
        init_targets
    output:
        summary=TASK_ROOT+"/{protein}/" + INIT_DIRNAME + '/' + P_SUMMARY, 
        selected=TASK_ROOT+"/{protein}/" + INIT_DIRNAME + '/' + P_SELECTED
    params:
        outdir=lambda wc: init_dir(wc.protein),
        graphdta_dir=lambda wc: method_dir(wc.protein, "GraphDTA"),
        hmsa_dir=lambda wc: method_dir(wc.protein, "HMSA"),
        colddta_dir=lambda wc: method_dir(wc.protein, "ColdDTA"),
        # evidti_dir=lambda wc: method_dir(wc.protein, "EviDTI"),  # DISABLED
        druglamp_dir=lambda wc: method_dir(wc.protein, "DrugLAMP"),
        conplex_dir=lambda wc: method_dir(wc.protein, "ConPLex"),
        exe=TOOLS["compile"],
        target_n=config["initial"]["n_keep"]
    shell:
        r"""
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate general
        mkdir -p "{params.outdir}"
        python "{params.exe}" \
            --graphdta-dir "{params.graphdta_dir}" \
            --hmsa-dir "{params.hmsa_dir}" \
            --colddta-dir "{params.colddta_dir}" \
            --druglamp-dir "{params.druglamp_dir}" \
            --conplex-dir "{params.conplex_dir}" \
            --target-n {params.target_n} \
            --summary "{output.summary}" \
            --selected "{output.selected}"
        """


# -------------------------
# Pre-fold the proteins with AF3 and Boltz2 to get the re-usable files
# -------------------------

rule write_prefold_af3:
    input:  
        csv=SEQS_CSV
    output: 
        json=TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/AF3/prefold/{protein}.json' 
    params:
        outdir=lambda wc: os.path.join(fine_dir(wc.protein), "AF3/prefold"),
        protein_name = lambda wc: wc.protein,
        exe=TOOLS["write_af3_protein"],
    shell:
        r"""
        mkdir -p "{params.outdir}"
        python "{params.exe}" \
            --output-dir "{params.outdir}" \
            --protein-name "{params.protein_name}" \
            --input-csv "{input.csv}" 
        """

rule write_prefold_boltz2:
    input:  
        csv=SEQS_CSV
    output: 
        yaml=TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/Boltz2/prefold/{protein}.yaml' 
    params:
        outdir_top = lambda wc: fine_dir(wc.protein),
        outdir_1=lambda wc: os.path.join(fine_dir(wc.protein), "Boltz2/"),
        outdir=lambda wc: os.path.join(fine_dir(wc.protein), "Boltz2/prefold"),
        protein_name = lambda wc: wc.protein,
        exe=TOOLS["write_boltz2_yaml"],
    shell:
        r"""
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate boltz
        mkdir -p "{params.outdir}"
        python "{params.exe}" \
            --output "{params.outdir}" \
            --protein-name "{params.protein_name}" \
            --protein-file "{input.csv}" \
	    --protein-only  
        """

rule prefold_af3:
    input:  
        json=TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/AF3/prefold/{protein}.json' 
    output: 
        TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/AF3/prefold/prefold.done'      
    params:
        outdir=lambda wc: os.path.join(fine_dir(wc.protein), "AF3/prefold"),
        af3_weight_dir="/shared/programs/af3_weights",
        af3_db_dir="/shared/programs/af3_data",
        protein_name = lambda wc: wc.protein,
        protein_name_lower = lambda wc: wc.protein.lower(),
        hmmer_dir="/home/ubuntu/Applications/hmmer/bin", 
        exe=TOOLS["af3"],
    shell:
        r"""
        set +eu
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate af3
        python -u "{params.exe}" \
            --json_path {input.json} \
            --model_dir {params.af3_weight_dir} \
            --db_dir {params.af3_db_dir} \
            --jackhmmer_binary_path {params.hmmer_dir}/jackhmmer \
            --hmmalign_binary_path {params.hmmer_dir}/hmmalign \
            --hmmbuild_binary_path {params.hmmer_dir}/hmmbuild \
            --hmmsearch_binary_path {params.hmmer_dir}/hmmsearch \
            --nhmmer_binary_path {params.hmmer_dir}/nhmmer \
            --output_dir {params.outdir}
        if [[ -f "{params.outdir}/{params.protein_name_lower}/{params.protein_name_lower}_summary_confidences.json" ]]; then
            touch "{output}"
        else
            echo "Expected AF3 output not found: {params.outdir}/{params.protein_name_lower}/{params.protein_name_lower}_summary_confidences.json" >&2
            exit 1
        fi
        """

rule prefold_boltz2:
    input:
        yaml=TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/Boltz2/prefold/{protein}.yaml'
    output:
        json="/home/ubuntu/{protein}/boltz2_tmp/boltz_results_{protein}/predictions/{protein}/confidence_{protein}_model_0.json"
    params:
        outdir="/home/ubuntu/{protein}/boltz2_tmp/",
        protein_name = lambda wc: wc.protein,
        exe=TOOLS["boltz2"],
    shell:
        r"""
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate boltz
        "{params.exe}" predict "{input.yaml}" \
            --out_dir={params.outdir} \
            --use_msa_server \
            --override
        """

# =========================================================================================
# FINE Screening
# =========================================================================================
# Write inputs for all methods
# -----------------------------------

# split the initial screening csv into multiple for VINA & Diffdock
checkpoint split_csv:
    input:
        selected=TASK_ROOT+"/{protein}/" + INIT_DIRNAME + '/' + P_SELECTED
    output:
        directory(TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/Vina/input/')
    params:
        outdir=lambda wc: os.path.join(fine_dir(wc.protein), "Vina/input/"),
        pdb_file=TASK_ROOT+"/Input/protein_file/{protein}/{protein}.pdb",
        protein_name = lambda wc: wc.protein,
        exe=TOOLS["split_csv"],
    shell:
        r"""
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate general
        mkdir -p {params.outdir}
        python "{params.exe}" \
            --protein-name={params.protein_name} \
            --pdb-file={params.pdb_file} \
            --smiles-file={input.selected} \
            --chunk-size=100 \
            --output-dir={params.outdir}
        """    

# write AF3 & Boltz2 input 
checkpoint write_af3_input:
    input:
        selected=TASK_ROOT+"/{protein}/" + INIT_DIRNAME + '/' + P_SELECTED,
        token=TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/AF3/prefold/prefold.done',
    output:
        token = TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/AF3/input/af3_input.done', 
        dir = directory(TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/AF3/input/')
    params:
        outdir=TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/AF3/input/' ,
        exe=TOOLS['write_af3_lig'],
        prefold=lambda wc: os.path.join(TASK_ROOT, wc.protein, FINE_DIRNAME, 'AF3/prefold/', wc.protein.lower(), f'{wc.protein.lower()}_data.json')
    shell:
        r"""
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate general
        echo {params.prefold}
        mkdir -p {params.outdir}
        python {params.exe} --output-dir {params.outdir} --input-json {params.prefold} --smiles-file {input.selected}
        touch {output.token}
        """

checkpoint write_boltz2_input:
    input:
        csv=SEQS_CSV,
        selected=TASK_ROOT+"/{protein}/" + INIT_DIRNAME + '/' + P_SELECTED,
        json="/home/ubuntu/{protein}/boltz2_tmp/boltz_results_{protein}/predictions/{protein}/confidence_{protein}_model_0.json"
    output:
        token = "/home/ubuntu/{protein}/boltz2_tmp/boltz_input.done",
        dir = directory("/home/ubuntu/{protein}/boltz2_tmp/input/")
    params:
        msa="/home/ubuntu/{protein}/boltz2_tmp/boltz_results_{protein}/msa/{protein}_0.csv",
        outdir="/home/ubuntu/{protein}/boltz2_tmp/input/",
        exe=TOOLS['write_boltz2_yaml'],
        protein_name = lambda wc: wc.protein
    shell:
        r"""
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate boltz
        mkdir -p {params.outdir}
        python {params.exe} --output {params.outdir} --msa {params.msa} \
             --smiles-path {input.selected} --protein-name {params.protein_name} --protein-file {input.csv}
        touch {output.token}
        """

# -----------------------------------
# Runing Screening
# -----------------------------------
rule vina:
    input:
        csv = TASK_ROOT + "/{protein}/" + FINE_DIRNAME + "/Vina/input/{part}.csv",
        seqs = SEQS_CSV,
    output:
        # use a unique marker file to avoid directory race conditions
        touch(TASK_ROOT + "/{protein}/" + FINE_DIRNAME + "/Vina/output/{part}.done"),
    params:
        exe = TOOLS["vina"],
        pdb_file = lambda wc: (
            TASK_ROOT + f"/Input/protein_file/{wc.protein}/{wc.protein}.pdb"
        ),
        box_center = lambda wc: get_box_str(wc.protein, SEQS_CSV)[0],
        box_size   = lambda wc: get_box_str(wc.protein, SEQS_CSV)[1],
        outdir = TASK_ROOT + "/{protein}/" + FINE_DIRNAME + "/Vina/output/{part}/",
    shell:
        r"""
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate vina        
        mkdir -p {params.outdir}
        python {params.exe} \
          --smiles {input.csv} \
          --pdb {params.pdb_file} \
          --box-center {params.box_center} \
          --box-size {params.box_size} \
          --output {params.outdir} \
	  --smiles-col "ligand_description" 
        touch {output}
        """
rule collect_vina:
    input:
        tokens=vina_targets
    output:
        csv=TASK_ROOT + "/{protein}/" + FINE_DIRNAME + "/Vina/results.csv"
    params:
        exe=TOOLS['get_vina_scores'],
        results_dir=TASK_ROOT + "/{protein}/" + FINE_DIRNAME + "/Vina/output/",
        indir=TASK_ROOT + "/{protein}/" + FINE_DIRNAME + "/Vina/input/",
        outdir=TASK_ROOT + "/{protein}/" + FINE_DIRNAME + "/Vina/"
    shell:
        r"""
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate general
        python {params.exe} \
           --vina-results-folder {params.results_dir} \
           --output-dir {params.outdir} \
           --input-dir {params.indir}
        """

# ensure i is numeric (accepts zero-padded like 0001)
wildcard_constraints:
    i=r"\d+",
    batch_id=r"\d+"

def job_batch(wc, size=100):
    # i starts at 1 → batches: 1–100 → 0, 101–200 → 1, etc.
    return (int(wc.i) - 1) // size

rule af3:
    input:
        token = TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/AF3/input/af3_input.done'
    output:
        touch(TASK_ROOT + "/{protein}/" + FINE_DIRNAME + "/AF3/output/token/{i}.done")      
    params:
        outdir=lambda wc: os.path.join(fine_dir(wc.protein), "AF3/output"),
        token_dir=lambda wc: os.path.join(fine_dir(wc.protein), "AF3/output/token"),
        af3_weight_dir="/shared/programs/af3_weights",
        af3_db_dir="/shared/programs/af3_data",
        json=lambda wc: os.path.join(TASK_ROOT, wc.protein,  FINE_DIRNAME, f"AF3/input/{wc.protein.lower()}_{wc.i}.json"),
        exe=TOOLS["af3"]
    group: "af3_jobs"
    shell:
        r"""
        set -euo pipefail
        source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
        conda activate af3
	
        mkdir -p {params.token_dir}	

        export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"
        export XLA_PYTHON_CLIENT_PREALLOCATE=true
        export XLA_CLIENT_MEM_FRACTION=0.95

        python {params.exe} \
            --json_path={params.json} \
            --model_dir={params.af3_weight_dir} \
            --db_dir={params.af3_db_dir} \
            --output_dir={params.outdir} \
            --norun_data_pipeline
        touch {output}
	"""

rule collect_af3:
    input:
        token=af3_targets
    output:
        csv=TASK_ROOT + "/{protein}/" + FINE_DIRNAME + "/AF3/summary.csv"
    params:
        exe=TOOLS['get_af3_scores'],
        outdir=TASK_ROOT + "/{protein}/" + FINE_DIRNAME + "/AF3/",
        results_dir=TASK_ROOT + "/{protein}/" + FINE_DIRNAME + "/AF3/output"
    shell:
        r"""
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate general        
        python {params.exe} \
            --af3-results-folder {params.results_dir} \
            --output-dir {params.outdir}
        """

rule boltz2:
    input:
        token = "/home/ubuntu/{protein}/boltz2_tmp/boltz_input.done"
    output:
        touch(TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/Boltz2/output/token/batch_{batch_id}.done'),
        tar=TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/Boltz2/output/batch_{batch_id}.tar.gz'
    params:
        final_outdir=lambda wc: os.path.join(fine_dir(wc.protein), "Boltz2/output/"),
        outdir_token=lambda wc: os.path.join(fine_dir(wc.protein), "Boltz2/output/token/"),
        protein_name = lambda wc: wc.protein,
        input_dir = "/home/ubuntu/{protein}/boltz2_tmp/input",
        selected = lambda wc: os.path.join(TASK_ROOT, wc.protein, INIT_DIRNAME, "selected.csv"),
        exe=TOOLS["boltz2"],
        batch_id=lambda wc: int(wc.batch_id),
        n_batches=40
    group: "boltz2_jobs"
    shell:
        r"""
        source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
        conda activate boltz

        mkdir -p {params.outdir_token}

        # Get total number of jobs
        SMILES_FILE="{params.selected}"
        N_SMILES=$(tail -n +2 "$SMILES_FILE" | wc -l)

        # Calculate batch size and range
        BATCH_SIZE=$(( ($N_SMILES + {params.n_batches} - 1) / {params.n_batches} ))
        START_IDX=$(( {params.batch_id} * $BATCH_SIZE ))
        END_IDX=$(( $START_IDX + $BATCH_SIZE ))
        if [ $END_IDX -gt $N_SMILES ]; then
            END_IDX=$N_SMILES
        fi

        echo "Processing batch {params.batch_id}: jobs $START_IDX to $(($END_IDX - 1))"

        # Create a batch-level temporary directory
        BATCH_LOCAL_OUT=/tmp/boltz2_batch_{params.protein_name}_{params.batch_id}
        mkdir -p $BATCH_LOCAL_OUT

        # Process each job in the batch
        for i in $(seq $START_IDX $(($END_IDX - 1))); do
            YAML_FILE="{params.input_dir}/${{i}}.yaml"

            if [ ! -f "$YAML_FILE" ]; then
                echo "Warning: YAML file not found: $YAML_FILE"
                continue
            fi

            # Use local NVMe for temporary storage - subdirectory per job
            LOCAL_OUT=$BATCH_LOCAL_OUT/job_${{i}}
            mkdir -p $LOCAL_OUT
            
            echo "=== Job $i start: $(date +%s) ==="
            echo "Running job $i from batch {params.batch_id}"

            # Run prediction on local storage
            "{params.exe}" predict "$YAML_FILE" \
                --out_dir=$LOCAL_OUT \
                --override
            echo "=== Job $i end: $(date +%s) ==="
        done

        # After all jobs complete, compress relevant files
        echo "Compressing results for batch {params.batch_id}"
        cd $BATCH_LOCAL_OUT
        find . \( -name "affinity_*.json" -o -name "confidence*.json" -o -name "*model*.cif" \) -print0 | \
            tar -czf batch_{params.batch_id}.tar.gz --null -T -

        # Copy tar file to final output directory
        cp batch_{params.batch_id}.tar.gz {params.final_outdir}/

        # Clean up local batch directory
        rm -rf $BATCH_LOCAL_OUT

        touch {output[0]}
        """

rule collect_boltz2:
    input:
        token=boltz2_targets,
        tars=get_boltz2_tar_files
    output:
        csv=TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/Boltz2/summary.csv'
    params:
        exe=TOOLS['get_boltz2_scores'],
        results_dir=TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/Boltz2/output/',
        outdir=TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/Boltz2/',
        selected=TASK_ROOT+"/{protein}/" + INIT_DIRNAME + '/' + P_SELECTED
    shell:
        r"""
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate general

        # Decompress all tar files
        echo "Decompressing batch tar files..."
        for tar_file in {input.tars}; do
            if [ -f "$tar_file" ]; then
                echo "Extracting $tar_file"
                tar -xzf "$tar_file" -C {params.results_dir}
            else
                echo "Warning: tar file not found: $tar_file"
            fi
        done

        # Run the analysis on decompressed files
        echo "Running analysis..."
        python {params.exe} \
            --boltz-results-folder {params.results_dir} \
            --output-dir {params.outdir} \
            --smiles {params.selected}
        """


rule diffdock:
    input:  
        csv = TASK_ROOT + "/{protein}/" + FINE_DIRNAME + "/Vina/input/{part}.csv",
    output: 
        token=TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/PBSA/DiffDock/output/{part}.done'
    params:
        outdir=lambda wc: os.path.join(fine_dir(wc.protein), "PBSA/DiffDock/output/"),
        protein_name = lambda wc: wc.protein,
        exe_dir=TOOLS["diffdock"],
        config='/home/ubuntu/Applications/DiffDock/default_inference_args.yaml'
    shell:
        r"""
        mkdir -p {params.outdir}
        set +eu
        source conda
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate diffdock
        cd {params.exe_dir}
        python -m inference --config {params.config} \
            --protein_ligand_csv {input.csv} \
            --out_dir {params.outdir}
        touch {output.token}
        """

rule md:
    input:
        token=diffdock_targets,
        gro=TASK_ROOT+'/Input/protein_file/{protein}/{protein}.gro',
        top=TASK_ROOT+'/Input/protein_file/{protein}/system_EM.top'
    output:
        token=TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/PBSA/PBSA/MD/{protein}_{i}/token.done'
    params:
        exe=TOOLS['md'],
        sdf=TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/PBSA/DiffDock/output/{protein}_{i}/rank1.sdf',
        outdir=TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/PBSA/PBSA/MD/{protein}_{i}/',
        script_dir=SCRIPT_ROOT+"/pbsa/"
    shell:
        r"""
        mkdir -p {params.outdir}
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate gmxMMPBSA
        bash {params.exe} {params.sdf} {params.outdir} {input.gro} {input.top} {params.script_dir}
        touch {output.token}
        """

rule pbsa:
    input:
        token=TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/PBSA/PBSA/MD/{protein}_{i}/token.done'
    output:
        token=TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/PBSA/PBSA/PBSA/{protein}_{i}/token.done'
    params:
        xtc=TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/PBSA/PBSA/MD/{protein}_{i}/T298.xtc',
        tpr=TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/PBSA/PBSA/MD/{protein}_{i}/T298.tpr',
        top=TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/PBSA/PBSA/MD/{protein}_{i}/system.top',
        ndx=TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/PBSA/PBSA/MD/{protein}_{i}/index.ndx',
        exe=TOOLS["pbsa"],
        gmx=TOOLS["gmx"],
        outdir=TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/PBSA/PBSA/PBSA/{protein}_{i}/',
        dat=TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/PBSA/PBSA/PBSA/{protein}_{i}/FINAL_RESULTS_MMPBSA.dat',
        script_dir=SCRIPT_ROOT+"/pbsa/"
    shell:
        r"""
        set +eu
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate gmxMMPBSA        
	bash -c "source {params.gmx}"
        export PATH="/home/ubuntu/miniconda3/envs/gmxMMPBSA/bin:$PATH"
        mkdir -p {params.outdir}
        cd {params.outdir}
        env "PATH=$PATH" {params.exe} -O \
            -i {params.script_dir}/mmpbsa.in -cs {params.tpr} -ct {params.xtc} \
            -ci {params.ndx} -cg 1 13 -cp {params.top} -o {params.dat} -eo {params.outdir}/FINAL_RESULTS_MMPBSA.csv
        touch {output.token}
        """

rule collect_pbsa:
    input:
        token=pbsa_targets,
        selected=TASK_ROOT+"/{protein}/" + INIT_DIRNAME + '/' + P_SELECTED
    output:
        csv=TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/PBSA/summary.csv'
    params:
        exe1=TOOLS['get_pbsa_scores'],
        exe2=TOOLS['map_pbsa_smiles'],
        tmp_csv=TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/PBSA/tmp.csv',
        pbsa_dir=TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/PBSA/PBSA/PBSA/',
        outdir=TASK_ROOT+"/{protein}/" + FINE_DIRNAME + '/PBSA/'
    shell:
        r"""
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate general
        bash {params.exe1} {params.pbsa_dir} {params.outdir}
        python {params.exe2} --collected {params.tmp_csv} \
          --smiles_csv {input.selected} \
          --outdir {params.outdir}
        rm {params.tmp_csv}
        """
