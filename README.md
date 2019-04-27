# Raw data and MCMC output underlying plots and results in the McQueen paper

## Reference datasets *D*<sub>B</sub>
The combination of sequence data from the Los Alamos HIV sequence database, the Stanford drug resistance database and the HIV positive selection database, run through our quality control pipeline and used for our reference dataset at the time of performing inference is here:
* `reference_datasets/HLA_inference/combined_{protease,reverse_transcriptase}_references.fasta`

The viral sequence data taken from the Stanford drug resistance database, run through our quality control pipeline and used for drug associated selection inference in the paper is here:
* `reference_datasets/Drug_inference/{protease,reverse_transcriptase}_references.fasta`

## Query datasets *D*
The viral sequence data taken from the Stanford drug resistance database, run through our quality control pipeline and used for drug associated selection inference in the paper is here:
* `query_datasets/Drug_inference/{protease,reverse_transcriptase}_references.fasta`

Associated host drug regime information for these patients is here:
* `csv_files/drug_regime_data/{protease,revesrse_transcriptase}_

## Main Text
### Figures
* **Figure 1**: The underlying process and our approximation.
  - N/A
* **Figure 2**: Simulation results summary; Inference results for Simulation studies 1 and 2.
  - Rdata file containing median, mean, and quantiles displayed in Figure 2_a_ and 2_b_: `Rdata_files/simulation_study_1/dsim_closest_100_exp_omega_rec_0.Rdata`
  - Rdata file containing median, mean and quantiles displayed in Figure 2_c_: `Rdata_files/simulation_study_1/dsim_closest_100_exp_omega_rec_0_01.Rdata`
  - Rdata file containing median, mean and quantiles displayed in Figure 2_d_: `Rdata_files/simulation_study_2/dgen_sim_HLA_coeff.Rdata`
  - Rdata files containing the 100 ROC curves for each independent run using the methods detailed in the paper, and the average ROC curve for each method: `Rdata_files/simulation_study_2/ROC_curves/data_for_ROC_curve_1460_leaves.Rdata`.
* **Figure 3**: Drug associated selection analysis; Protease results summary.
  - Raw output from MCMC run: `mcmc_output_files/Drug_inference/protease_drugs_combined_{window,log}.txt.zip`
  - Rdata file, in nicer array format with removal of burn-in: `Rdata_files/Drug_inference/protease_drugs.Rdata`
* **Figure 4**: Selected results summaries.
  - Raw output from MCMC run: `mcmc_output_files/HLA_inference/reverse_transcriptase_combined_subtype_{B,C}_ref_combined_{window,log}.txt.zip`.
  - Rdata file, in nicer array format with removal of burn-in: `Rdata_files/HLA_inference/reverse_transcriptase_combined_subtype_{B,C}_ref.Rdata`.
* **Figure 5**:  HLA-B associated selection in reverse transcriptase.
  - Raw output from MCMC run: `mcmc_output_files/HLA_inference/reverse_transcriptase_combined_subtype_{B,C}_ref_combined_{window,log}.txt.zip`.
  - Rdata file, in nicer array format, with HLA labels and removal of burn-in: `Rdata_files/HLA_inference/reverse_transcriptase_combined_subtype_{B,C}_ref.Rdata`.

### Tables
* Table 1: Data used in this study: public databases.
  - N/A
* Table 2: Data used in this study: viral sequence data with associated host HLA information.
  - N/A
* Table 3: Drug associated selection analysis: sensitivity of inference.
  - Raw output from MCMC run: `mcmc_output_files/Drug_inference/{protease,reverse_transcriptase}_drugs_combined_{window,log}.txt.zip`
  - Rdata file, in nicer array format with removal of burn-in: `Rdata_files/Drug_inference/{protease,reverse_transcriptase}_drugs.Rdata`
* Table 4: Overlap between sites under selection and known HLA epitopes.
  - Raw output from MCMC run: `mcmc_output_files/{protease, reverse_transcriptase}_combined_subtype_{B,C}_ref_combined_{window, log}.txt.zip`.
  - Rdata file, in nicer array format, with HLA labels and removal of burn-in: `Rdata_files/HLA_inference/{protease, reverse_transcriptase}_combined_subtype_{B,C}_ref.Rdata`.
* Table 5: Codons in reverse transcriptase showing evidence for HLA associated selection. 
  - Raw output from MCMC run: `mcmc_output_files/reverse_transcriptase_combined_subtype_{B,C}_ref_combined_{window, log}.txt.zip`.
  - Rdata file, in nicer array format, with HLA labels and removal of burn-in: `Rdata_files/HLA_inference/reverse_transcriptase_combined_subtype_{B,C}_ref.Rdata`.
* Table 6: Codons in protease showing evidence for HLA associated selection. 
  - Raw output from MCMC run: `mcmc_output_files/protease_combined_subtype_{B,C}_ref_combined_{window, log}.txt.zip`.
  - Rdata file, in nicer array format, with HLA labels and removal of burn-in: `Rdata_files/HLA_inference/protease_combined_subtype_{B,C}_ref.Rdata`.

## Supplementary Materials
### Figures
* **Figure 1**: Simulation study 1: summary of the selection profiles used.
  - N/A
* **Figure 2**: Simulation study 1: parameter estimates, no recombination (_r_= 0)
  - Rdata file containing median, mean, and quantiles displayed: `Rdata_files/simulation_study_1/dsim_closest_100_exp_omega_rec_0.Rdata`
* **Figure 3**: Simulation Study 1: comparison of estimates obtained using the closest 10 sequences by Hamming distance to the closest 100 sequences for _r_= 0 and _r_= 0.01.
  - Rdata file containing median, mean, and quantiles displayed: `Rdata_files/simulation_study_1/dsim_closest_{10,100}_exp_omega_rec_{0,0_01}.Rdata`
* **Figure 4**: Simulation Study 1: parameter estimates, _r_= 0.01.
  - Rdata file containing median, mean, and quantiles displayed: `Rdata_files/simulation_study_1/dsim_closest_100_exp_omega_rec_0_01.Rdata`
* **Figure 5**:  Simulation Study 1: comparison of estimates obtained using the closest 100 sequences by Hamming distance to the sequences actually copied from in the simulation, _r_= 0.05. 
  - Rdata file containing median, mean, and quantiles displayed: 
  - `Rdata_files/simulation_study_1/dsim_closest_100_exp_omega_rec_0_05.Rdata`
* **Figure 6**: Simulation Study 2: parameter estimates.  
  - Rdata file containing median, mean, and quantiles displayed:
  - `Rdata_files/simulation_study_2/dgen_sim_HLA_coeff.Rdata`
* **Figure 7**: Simulation Study 2: sample size of 3000, results summary.
  - Rdata files containing _p_ values, medians, or estimated selection coefficient information used to generate ROC curves: `Rdata_files/simulation_study_2/dgen_3000_leaves/*`. The separate files contain results for the different methods and different subsets of the reference dataset:
    - `bhatt_and_fisher_list_3000_leaves.Rdata`: Fisher's Exact test _p_-values as in [Moore _et al._](https://science.sciencemag.org/content/sci/296/5572/1439.full.pdf) and Phylogeny aware _p_-values using the method of [Bhattacharya _et_al._](https://science.sciencemag.org/content/sci/315/5818/1583.full.pdf).
    - `helen_rate_list_3000_leaves.Rdata`: Escape rate estimates using the methods of [Fryer _et al._](https://journals.plos.org/plospathogens/article/file?id=10.1371/journal.ppat.1001196&type=printable).
    - `carlson_llk_list_and_p_3000_leaves.Rdata`: Using the methods of [Carlson _et al._](https://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1000225&type=printable).
    - `carlson_llk_list_and_p_OR_3000_leaves.Rdata`: Using the methods of [Carlson _et al._](https://jvi.asm.org/content/jvi/86/9/5230.full.pdf).
    - `dgen_sim_HLA_coeff_3000_leaves_sims_average.Rdata`: Using all available reference sequence data from public databases.
* **Figure 8**: Simulation Study 2: impact of reduced reference sequence set, results summary.
  - Rdata files containing _p_ values, medians, or estimated selection coefficient information used to generate ROC curves: `Rdata_files/simulation_study_2/dgen_1460_leaves/*`. The separate files contain results for the different methods and different subsets of the reference dataset:
    - `bhatt_and_fisher_list.Rdata`: Fisher's Exact test _p_-values as in [Moore _et al._](https://science.sciencemag.org/content/sci/296/5572/1439.full.pdf) and Phylogeny aware _p_-values using the method of [Bhattacharya _et_al._](https://science.sciencemag.org/content/sci/315/5818/1583.full.pdf).
    - `helen_rate_list.Rdata`: Escape rate estimates using the methods of [Fryer _et al._](https://journals.plos.org/plospathogens/article/file?id=10.1371/journal.ppat.1001196&type=printable).
    - `carlson_llk_list_and_p.Rdata`: Using the methods of [Carlson _et al._](https://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1000225&type=printable).
    - `carlson_llk_list_and_p_OR.Rdata`: Using the methods of [Carlson _et al._](https://jvi.asm.org/content/jvi/86/9/5230.full.pdf).
    - `dgen_sim_HLA_coeff_sims_average.Rdata`: Using all available reference sequence data from public databases.
    - `dgen_sim_HLA_coeff_combined_sims_average.Rdata`: As above, plus query sequences using a leave one out approach.
    - `dgen_sim_HLA_coeff_divide_10_sims_average.Rdata`: Using 1/10 of HIV sequence data from public databases.
    - `dgen_sim_HLA_coeff_combined_divide_10_sims_average.Rdata`: As above, plus query sequences using a leave one out approach.
    - `dgen_sim_HLA_coeff_divide_100_sims_average.Rdata`: Using 1/100 of HIV sequence data from public databases.
    - `dgen_sim_HLA_coeff_combined_divide_100_sims_average.Rdata`: As above, plus query sequences using a leave one out approach.
    - `dgen_sim_HLA_coeff_same_refs_sims_average.Rdata`: Using only query sequences in the reference set, using a leave one out approach.
  
* **Figure 9**: Simulation Study 3: example results for HLA-A alleles.
  - Rdata files containing quantile information (2.5%, 25%, 50%, 75%, 97.5%) for simulation study 3:
  `Rdata_files/simulation_study_3/dgen_sim_HLA_coeff_combined_botswana_{*}_{average,quantiles}.Rdata`
    - The separate folders contain results for different reference datasets whether we used a leave one out approach for the query datasets:
      - `sims`: The reference dataset that was used for inference in the paper (combining across all available public databases)
      - `add_queries_to_ref_set_sims`: As above, but also including the query sequences using a leave one out approach.
      - `same_refs_sims`: Only using the query sequences, using a leave one out approach.
      - `los_alamos_all_sims`: Include only sequence data taken from the Los Alamos sequence database in the reference dataset.
      - `los_alamos_all_add_queries_to_ref_set_sims`: As above, but also including the query sequences using a leave one out approach.
      - `los_alamos_remove_botswana_sims`: Include only sequence data taken from the Los Alamos sequence database in the reference dataset, and restrict to samples not from Botswana.
      - `los_alamos_remove_botswana_add_queries_to_ref_set_sims`: As above, but also including the query sequences using a leave one out approach.
      - `gold_standard_sims`: Use the sequences actually simulated from as the reference dataset.
      
* **Figure 10**: Simulation Study 3: example results for HLA-B alleles. 
  - As in Supplementary Figure 9.
* **Figure 11**: Simulation Study 3: example results for HLA-C alleles.
  - As in Supplementary Figure 9.
* **Figure 12**: Drug associated selection analysis: reverse transcriptase results summary.
  - Raw output from MCMC run: `mcmc_output_files/Drug_inference/reverse_transcriptase_drugs_combined_{window,log}.txt.zip`
  - Rdata file, in nicer array format with removal of burn-in: `Rdata_files/Drug_inference/reverse_transcriptase_drugs.Rdata`
* Table 4: Overlap between sites under selection and known HLA epitopes.
* **Figure 13**: Hierarchical clustering of drug selection profiles.
  - `Rdata_files/Drug_inference/{protease,reverse_transcriptase}_drugs.Rdata`
* **Figure 14**: Summary of HLA-A associated selection in protease. 
  - Raw output from MCMC run: `mcmc_output_files/protease_combined_subtype_{B,C}_ref_combined_{window, log}.txt.zip`.
  - Rdata file, in nicer array format, with HLA labels and removal of burn-in: `Rdata_files/HLA_inference/reverse_transcriptase_combined_subtype_{B,C}_ref.Rdata`.
* **Figure 15**: Summary of HLA-B associated selection in protease. 
  - Raw output from MCMC run: `mcmc_output_files/protease_combined_subtype_{B,C}_ref_combined_{window, log}.txt.zip`.
  - Rdata file, in nicer array format, with HLA labels and removal of burn-in: `Rdata_files/HLA_inference/reverse_transcriptase_combined_subtype_{B,C}_ref.Rdata`.
* **Figure 16**: Summary of HLA-C associated selection in protease. 
  - Raw output from MCMC run: `mcmc_output_files/protease_combined_subtype_{B,C}_ref_combined_{window, log}.txt.zip`.
  - Rdata file, in nicer array format, with HLA labels and removal of burn-in: `Rdata_files/HLA_inference/reverse_transcriptase_combined_subtype_{B,C}_ref.Rdata`.
* **Figure 17**: Summary of HLA-A associated selection in reverse transcriptase.
  - Raw output from MCMC run: `mcmc_output_files/reverse_transcriptase_combined_subtype_{B,C}_ref_combined_{window, log}.txt.zip`.
  - Rdata file, in nicer array format, with HLA labels and removal of burn-in: `Rdata_files/HLA_inference/reverse_transcriptase_combined_subtype_{B,C}_ref.Rdata`.
* **Figure 18**: Summary of HLA-C associated selection in reverse transcriptase.
  - Raw output from MCMC run: `mcmc_output_files/reverse_transcriptase_combined_subtype_{B,C}_ref_combined_{window, log}.txt.zip`.
  - Rdata file, in nicer array format, with HLA labels and removal of burn-in: `Rdata_files/HLA_inference/reverse_transcriptase_combined_subtype_{B,C}_ref.Rdata`.
* **Figure 19**: Estimated recombination probabilities between sites across protease and reverse transcriptase.
  - `mcmc_output_files/{protease,reverse_transcriptase}_combined_subtype_{B,C}_ref_combined_log.txt.zip`.
* **Figure 20**: The impact of HLA alleles on differentiation between HIV-1 subtypes.
  - `Rdata_files/HLA_inference/{protease,reverse_transcriptase}_selection_list_all_HLA.Rdata`

### Tables
* **Table 1**: Priors on parameters of the model.
  - N/A
* **Table 2**: The collection of ‘false positive’ associations are interrogated.
  - N/A
* **Table 3**: The collection of apparent false negative associations are interrogated. 
  - N/A
* **Table 4**: The collection of top-tier and second-tier sites in protease.
  - Raw output from MCMC run: `mcmc_output_files/protease_combined_subtype_{B,C}_ref_combined_{window, log}.txt.zip`.
  - Rdata file, in nicer array format, with HLA labels and removal of burn-in: `Rdata_files/HLA_inference/protease_combined_subtype_{B,C}_ref.Rdata`.
  - Summary information containing the median, lower 2.5% and 10% quantile information: `Rdata_files/HLA_inference/protease_selection_list_all_HLA.Rdata`.
* **Table 5**: The collection of top-tier and second-tier sites in reverse transriptase.
  - Raw output from MCMC run: `mcmc_output_files/reverse_transcriptase_combined_subtype_{B,C}_ref_combined_{window, log}.txt.zip`.
  - Rdata file, in nicer array format, with HLA labels and removal of burn-in: `Rdata_files/HLA_inference/reverse_transcriptase_combined_subtype_{B,C}_ref.Rdata`.
  - Summary information containing the median, lower 2.5% and 10% quantile information: `Rdata_files/HLA_inference/reverse_transcriptase_selection_list_all_HLA.Rdata`.
* **Table 6**: Overlap between studies, protease.
  - Summary information containing the median, lower 2.5% and 10% quantile information: `Rdata_files/HLA_inference/protease_selection_list_all_HLA.Rdata`.
  - Also, summarise where I got the Carlson information from.
* **Table 7**: Overlap between studies, reverse transcriptase.
  - Summary information containing the median, lower 2.5% and 10% quantile information: `Rdata_files/HLA_inference/protease_selection_list_all_HLA.Rdata`.
* **Table 8**: Table of _p_-values for differences between _H_<sub>B</sub> and _H_<sub>C</sub>) associated selection at sites that distinguish subtype B and subtype C in protease and reverse transcriptase.
  - TO DO
* **Table 9**: Summary of reference datasets for simulation study 3.
  - As in Supplementary Figure 9.

## Source datafile - to be available for download with manuscript.
These will be small in filesize so that they may be downloaded from the journal website and provide sufficient data to recreate the various plots.

## Main text
* **Figure 2** `source_data_files/Figure_2.Rdata` is a list containing all the information to create Figure 2. The list has a collection of sublists:
  * `simulation_study_1_closest_100_R_0`: information for simulation study 1 with the recombination rate between neighbouring sites set to 0 and using the closest 100 sequences by Hamming distance as references.
  * `simulation_study_1_closest_100_R_0_01`: information for simulation study 1 with the recombination rate between neighbouring sites set to 0.01 and using the closest 100 sequences by Hamming distance as references.
  * `simulation_study_2`: information for simulation study 2 estimations after simulating under a birth-death process.
  * `simulation_study_2_ROC`: ROC curve information.

  Each of the first three lists consists of the following entries:
    * `truth`: a list consisting of the true underlying host associated selection parameters (`HLA_esc_truth`) and reversion parameters `rev_truth`.
    * `omega_truth`: the true underlying dN/dS ratio for the codon model, for one of the 100 independent runs.
    * `omega_mean`: the mean dN/dS estimates for one of the 100 independent runs.
    * `omega_median`: the median dN/dS estimates for one of the 100 independent runs.
    * `omega_credible_95`: the 95% credible interval for dN/dS estimates for one of the 100 independent runs.
    * `omega_credible_50`: the 50% credible interval for dN/dS estimates for one of the 100 independent runs.
    * `R_truth`: the true underlying recombination probability between adjacent sites for the simulation study.
    * `mean_R`: the mean estimate of recombination between adjacent sites.
    * `median_R`: the median estimate of recombination between adjacent sites.
    * `credible_95_R`: the average 95% credible interval across the 100 independent runs.
    * `credible_50_R`: the average 50% credible interval across the 100 independent runs.
    * `mean_escape`: the mean host associated selection estimates.
    * `median_escape`: the median host associated selection estimates.
    * `coverage_95_escape`: the coverage of the 95% credible interval across the 100 independent runs.
    * `coverage_50_escape`: the coverage of the 50% credible interval across the 100 independent runs.
    * `credible_95_escape`: the average 95% credible interval across the 100 independent runs.
    * `credible_50_escape`: the average 50% credible interval across the 100 independent runs.
    * `mean_reversion`: the mean reversion estimates at each site.
    * `median_reversion`: the median reversion estimates at each site.
    * `coverage_95_reversion`: the coverage of the 95% credible interval across the 100 independent runs.
    * `coverage_50_reversion`: the coverage of the 50% credible interval across the 100 independent runs.
    * `credible_95_reversion`: the average 95% credible interval across the 100 independent runs.
    * `credible_50_reversion`: the average 50% credible interval across the 100 independent runs.

  The final list `simulation_study_2_ROC` consists of the following entries:
    * `bhattacharya_average`: average ROC curve using the methods of [Bhattacharya _et_al._](https://science.sciencemag.org/content/sci/315/5818/1583.full.pdf).
    * `phyloD_average`: average ROC curve using the methods of [Carlson _et al._](https://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1000225&type=printable).
    * `phyloD_false_positive`: rolling false positive rate as threshold is decreased.
    * `phyloD_true_positive`: rolling true positive rate as threshold is decreased.
    * `phyloD_OR_average`: average ROC curve using the methods of [Carlson _et al._](https://jvi.asm.org/content/jvi/86/9/5230.full.pdf).
    * `phyloD_OR_false_positive`: rolling false positive rate as threshold is decreased.
    * `phyloD_OR_true_positive`: rolling true positive rate as threshold is decreased.
    * `fishers_exact_average`: average ROC curve using the methods of [Moore _et al._](https://science.sciencemag.org/content/sci/296/5572/1439.full.pdf).
    * `fryer_average`: average ROC curve using the methods of [Fryer _et al._](https://journals.plos.org/plospathogens/article/file?id=10.1371/journal.ppat.1001196&type=printable).
    * `fryer_false_positive`: rolling false positive rate as threshold is decreased.
    * `fryer_true_positive`: rolling true positive rate as threshold is decreased.
    * `our_method_average`: average ROC curve using our method, McQueen.
    * `our_false_positive`: rolling false positive rate as threshold is decreased.
    * `our_true_positive`: rolling true positive rate as threshold is decreased.
  
* **Figure 3** `source_data_files/Figure_3.Rdata` contains the 2.5%, 10% and 50% quantiles with rows labelled by drug, and the median omega.
* **Figure 4** `source_data_files/Figure_4.Rdata` contains the 2.5%, 10% and 50% quantiles for HLA associated selection, and median omega for inference in reverse transcriptase. Information is contained in named lists with 2 elements, one from subtype B consensus and one from subtype C consensus. We also include the number of individuals harbouring each of the HLA types.
* **Figure 5** `source_data_files/Figure_5.Rdata` contains the 2.5%, 10% and 50% quantiles for HLA associated selection, and median omega for inference in reverse transcriptase. Information is contained in named lists with 2 elements, one from subtype B consensus and one from subtype C consensus. We also include the number of individuals harbouring each of the HLA types.

## Supplementary Materials
* **Figure 2** `source_data_files/Supplementary_Figure_2.Rdata` contains the following information for simulation study 1 with recombination probability set to 0 between neighbouring sites and using the closest 100 reference samples by Hamming distance:
  * `truth`: a list consisting of the true underlying host associated selection parameters (`HLA_esc_truth`) and reversion parameters `rev_truth`.
  * `omega_truth`: the true underlying dN/dS ratio for the codon model, for one of the 100 independent runs.
  * `omega_mean`: the mean dN/dS estimates for one of the 100 independent runs.
  * `omega_median`: the median dN/dS estimates for one of the 100 independent runs.
  * `omega_credible_95`: the 95% credible interval for dN/dS estimates for one of the 100 independent runs.
  * `omega_credible_50`: the 50% credible interval for dN/dS estimates for one of the 100 independent runs.
  * `R_truth`: the true underlying recombination probability between adjacent sites for the simulation study.
  * `mean_R`: the mean estimate of recombination between adjacent sites.
  * `median_R`: the median estimate of recombination between adjacent sites.
  * `credible_95_R`: the average 95% credible interval across the 100 independent runs.
  * `credible_50_R`: the average 50% credible interval across the 100 independent runs.
  * `mean_escape`: the mean host associated selection estimates.
  * `median_escape`: the median host associated selection estimates.
  * `coverage_95_escape`: the coverage of the 95% credible interval across the 100 independent runs.
  * `coverage_50_escape`: the coverage of the 50% credible interval across the 100 independent runs.
  * `credible_95_escape`: the average 95% credible interval across the 100 independent runs.
  * `credible_50_escape`: the average 50% credible interval across the 100 independent runs.
  * `mean_reversion`: the mean reversion estimates at each site.
  * `median_reversion`: the median reversion estimates at each site.
  * `coverage_95_reversion`: the coverage of the 95% credible interval across the 100 independent runs.
  * `coverage_50_reversion`: the coverage of the 50% credible interval across the 100 independent runs.
  * `credible_95_reversion`: the average 95% credible interval across the 100 independent runs.
  * `credible_50_reversion`: the average 50% credible interval across the 100 independent runs.
* **Figure 3** `source_data_files/Supplementary_Figure_3.Rdata` contains lists of information for output from simulation study 1:
  * `closest_10_R_0` 
  * `closest_100_R_0`
  * `closest_10_R_0_01`
  * `closest_100_R_0_01`
  
  Each contains the information as in `source_data_files/Supplementary_Figure_2.Rdata`.
* **Figure 4** `source_data_files/Supplementary_Figure_4.Rdata` contains the information as in `source_data_files/Supplementary_Figure_2.Rdata` for simulation study 1 with recombination probability set to 0.01 between neighbouring sites and using the closest 100 reference samples by Hamming distance.
* **Figure 5** `source_data_files/Supplementary_Figure_5.Rdata` contains two lists of information for output from simulation study 1 with recombination probability set to 0.05 between neighbouring sites:
  * `closest_100_R_0_05`
  * `closest_100_R_0_05_actually_copied_from`
  
  Each contains the information as in `source_data_files/Supplementary_Figure_2.Rdata`.
* **Figure 6** 
* **Figure 7** `source_data_files/Supplementary_Figure_7.Rdata` contains the collection of ROC curves for the compared methods when simulating trees with 3000 leaves.
* **Figure 8** `source_data_files/Supplementary_Figure_8.Rdata` contains the collection of ROC curves for the compared methods when simulating trees with 1460 leaves.
* **Figure 9** `source_data_files/Supplementary_Figure_9.Rdata` contains a collection of lists:
  * `our_reference_set`: All available sequence data from the Los Alamos HIV sequence database, the Stanford drug resistance database and the HIV positive selection database.
  * `our_reference_set_LOO`: All available sequence data from the Los Alamos HIV sequence database, the Stanford drug resistance database and the HIV positive selection database, plus all query sequences using a leave one out approach.
  * `los_alamos_reference_set`: Restrict to only Los Alamos sequence data in the reference dataset.
  * `los_alamos_reference_set_LOO`: Restrict to only Los Alamos sequence data in the reference dataset, plus all query sequences using a leave one out approach.
  * `los_alamos_reference_set_remove_botswana`: Restrict to only Los Alamos sequence data in the reference dataset and remove all sequences sampled from Botswana.
  * `los_alamos_reference_set_remove_botswana_LOO`: Restrict to only Los Alamos sequence data in the reference dataset and remove all sequences sampled from Botswana and include all query sequences using a leave one out approach.
  * `only_queries`: Use query sequences as the reference dataset using a leave one out approach.
  * `gold_standard_reference_set`: Use the sequences actually copied from as the reference dataset.
Each list contains the following data:
  * `median_escape`: The median escape estimates across the runs site by site for each HLA.
  * `coverage_50_escape`: The coverage of the 50% credible invervals site by site for each HLA.
  * `coverage_95_escape`: The coverage of the 95% credible intervals site by site for each HLA.
  * `RMSE`: The root mean squared error site by site for each HLA.
In addition, there is also `HLA_truth` which contains the true underlying selection coefficients for the parametric bootstrap.
* **Figure 10** `source_data_files/Supplementary_Figure_10.Rdata` contains the same information as in `source_data_files/Supplementary_Figure_9.Rdata`.
* **Figure 11** `source_data_files/Supplementary_Figure_11.Rdata` contains the same information as in `source_data_files/Supplementary_Figure_9.Rdata`.
* **Figure 12** `source_data_files/Supplementary_Figure_12.Rdata` contains the 2.5%, 10% and 50% quantiles with rows labelled by drug.
* **Figure 13** `source_data_files/Supplementary_Figure_13.Rdata` contains `drug_selection_median`; the median drug associated selection coefficients, and `trees`; native `ape` format (requires the R library `ape`) for the trees determined using `hclust`.
* **Figure 14** `source_data_files/Supplementary_Figure_14.Rdata` contains the 2.5%, 10% and 50% quantiles for HLA associated selection, and median omega for inference in protease. Information is contained in named lists with 2 elements, one from subtype B consensus and one from subtype C consensus. We also include the number of individuals harbouring each of the HLA types.
* **Figure 15** `source_data_files/Supplementary_Figure_15.Rdata` contains the 2.5%, 10% and 50% quantiles for HLA associated selection, and median omega for inference in protease. Information is contained in named lists with 2 elements, one from subtype B consensus and one from subtype C consensus. We also include the number of individuals harbouring each of the HLA types.
* **Figure 16** `source_data_files/Supplementary_Figure_16.Rdata` contains the 2.5%, 10% and 50% quantiles for HLA associated selection, and median omega for inference in protease. Information is contained in named lists with 2 elements, one from subtype B consensus and one from subtype C consensus. We also include the number of individuals harbouring each of the HLA types.
* **Figure 17** `source_data_files/Supplementary_Figure_17.Rdata` contains the 2.5%, 10% and 50% quantiles for HLA associated selection, and median omega for inference in reverse transcriptase. Information is contained in named lists with 2 elements, one from subtype B consensus and one from subtype C consensus. We also include the number of individuals harbouring each of the HLA types.
* **Figure 18** `source_data_files/Supplementary_Figure_18.Rdata` contains the 2.5%, 10% and 50% quantiles for HLA associated selection, and median omega for inference in reverse transcriptase. Information is contained in named lists with 2 elements, one from subtype B consensus and one from subtype C consensus. We also include the number of individuals harbouring each of the HLA types.
* **Figure 19** `source_data_files/Supplementary_Figure_19.Rdata` contains the mean, median, 2.5%, 25%, 75%, and 97.5% quantiles for probability of recombination between adjacent sites in protease and reverse transcriptase, away from subtype B and subtype C consensus viral sequence.
* **Figure 20** `source_data_files/Supplementary_Figure_20.Rdata` contains estimates of HLA frequency weighted selection away from subtype B and subtype C viral consensus:
  * `log_protease_from_B_by_subtype_B`: log base 10 of selection away from subtype B protease consensus by HLA frequency distribution of hosts harbouring subtype B viruses.
  * `log_protease_from_B_by_subtype_C`: log base 10 of selection away from subtype B protease consensus by HLA frequency distribution of hosts harbouring subtype C viruses.
  * `log_protease_from_C_by_subtype_B`: log base 10 of selection away from subtype C protease consensus by HLA frequency distribution of hosts harbouring subtype B viruses.
  * `log_protease_from_C_by_subtype_C`: log base 10 of selection away from subtype C protease consensus by HLA frequency distribution of hosts harbouring subtype C viruses.
  * `log_RT_from_B_by_subtype_B`: log base 10 of selection away from subtype B reverse transcriptase consensus by HLA frequency distribution of hosts harbouring subtype B viruses.
  * `log_RT_from_B_by_subtype_C`: log base 10 of selection away from subtype B reverse transcriptase consensus by HLA frequency distribution of hosts harbouring subtype C viruses.
  * `log_RT_from_C_by_subtype_B`: log base 10 of selection away from subtype C reverse transcriptase consensus by HLA frequency distribution of hosts harbouring subtype B viruses.
  * `log_RT_from_C_by_subtype_C`: log base 10 of selection away from subtype C reverse transcriptase consensus by HLA frequency distribution of hosts harbouring subtype C viruses.
We also included the collection of sites that differ between subtype B and subtype C in reverse transcriptase: `B_C_diffs_protease`, `B_C_diffs_RT`.
