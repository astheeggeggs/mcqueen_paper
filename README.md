# Raw data and MCMC output underlying plots and results in the McQueen paper

## Main Text
### Figures
* **Figure 1**: The underlying process and our approximation.
  - N/A
* **Figure 2**: Simulation results summary; Inference results for Simulation studies 1 and 2.
  - TO DO
* **Figure 3**: Drug associated selection analysis; Protease results summary.
  - Raw output from MCMC run: `mcmc_output_files/Drug_inference/protease_drugs_combined_{window,log}.txt.zip`
  - Rdata file, in nicer array format with removal of burn-in: `Rdata_files/Drug_inference/protease_drugs.Rdata`
* **Figure 4**: Selected results summaries.
  - Raw output from MCMC run: `mcmc_output_files/reverse_transcriptase_combined_subtype_{B,C}_ref_combined_{window,log}.txt.zip`.
  - Rdata file, in nicer array format with removal of burn-in: `Rdata_files/HLA_inference/reverse_transcriptase_combined_subtype_{B,C}_ref.Rdata`.
* **Figure 5**:  HLA-B associated selection in reverse transcriptase.
  - Raw output from MCMC run: `mcmc_output_files/reverse_transcriptase_combined_subtype_{B,C}_ref_combined_{window,log}.txt.zip`.
  - .Rdata file, in nicer array format, with HLA labels and removal of burn-in: `Rdata_files/HLA_inference/reverse_transcriptase_combined_subtype_{B,C}_ref.Rdata`.

### Tables
* Table 1: Data used in this study:  public databases.
  - N/A
* Table 2: Data used in this study: viral sequence data with associated host HLA information.
  - N/A
* Table 3: Drug associated selection analysis: sensitivity of inference.
  - Raw output from MCMC run: `mcmc_output_files/Drug_inference/{protease,reverse_transcriptase}_drugs_combined_{window,log}.txt.zip`
  - Rdata file, in nicer array format with removal of burn-in: `Rdata_files/Drug_inference/{protease,reverse_transcriptase}_drugs.Rdata`
* Table 4: Overlap between sites under selection and known HLA epitopes.
  - Raw output from MCMC run: `mcmc_output_files/{protease, reverse_transcriptase}_combined_subtype_{B,C}_ref_combined_{window, log}.txt.zip`.
  - .Rdata file, in nicer array format, with HLA labels and removal of burn-in: `Rdata_files/HLA_inference/{protease, reverse_transcriptase}_combined_subtype_{B,C}_ref.Rdata`.
* Table 5: Codons in reverse transcriptase showing evidence for HLA associated selection. 
  - Raw output from MCMC run: `mcmc_output_files/reverse_transcriptase_combined_subtype_{B,C}_ref_combined_{window, log}.txt.zip`.
  - .Rdata file, in nicer array format, with HLA labels and removal of burn-in: `Rdata_files/HLA_inference/reverse_transcriptase_combined_subtype_{B,C}_ref.Rdata`.
* Table 6: Codons in protease showing evidence for HLA associated selection. 
  - Raw output from MCMC run: `mcmc_output_files/protease_combined_subtype_{B,C}_ref_combined_{window, log}.txt.zip`.
  - .Rdata file, in nicer array format, with HLA labels and removal of burn-in: `Rdata_files/HLA_inference/protease_combined_subtype_{B,C}_ref.Rdata`.

## Supplementary Materials
### Figures
* **Figure 1**: Simulation study 1: summary of the selection profiles used.
  - N/A
* **Figure 2**: Simulation study 1: parameter estimates, no recombination (_r_= 0)
  - TO DO
* **Figure 3**: Simulation Study 1: comparison of estimates obtained using the closest 10 sequences by Hamming distance to the closest 100 sequences for _r_= 0 and _r_= 0.01.
  - TO DO
* **Figure 4**: Simulation Study 1: parameter estimates, _r_= 0.01.
  - TO DO
* **Figure 5**:  Simulation Study 1: comparison of estimates obtained using the closest 100 sequences by Hamming distance to the sequences actually copied from in the simulation, _r_= 0.05. 
  - TO DO
* **Figure 6**: Simulation Study 2: parameter estimates.  
  - TO DO
* **Figure 7**: Simulation Study 2: sample size of 3000, results summary.
  - .Rdata files containing _p_ values, medians, or estimated selection coefficient information used to generate ROC curves: `Rdata_files/simulation_study_2/dgen_3000_leaves/*`. The separate files contain results for the different methods and different subsets of the reference dataset:
    - `bhatt_and_fisher_list_3000_leaves.Rdata`: Fisher's Exact test _p_-values as in [Moore _et al._](https://science.sciencemag.org/content/sci/296/5572/1439.full.pdf) and Phylogeny aware _p_-values using the method of [Bhattacharya _et_al._](https://science.sciencemag.org/content/sci/315/5818/1583.full.pdf).
    - `helen_rate_list_3000_leaves.Rdata`: Escape rate estimates using the methods of [Fryer _et al._](https://journals.plos.org/plospathogens/article/file?id=10.1371/journal.ppat.1001196&type=printable).
    - `carlson_llk_list_and_p_3000_leaves.Rdata`: Using the methods of [Carlson _et al._](https://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1000225&type=printable).
    - `carlson_llk_list_and_p_OR_3000_leaves.Rdata`: Using the methods of [Carlson _et al._](https://jvi.asm.org/content/jvi/86/9/5230.full.pdf).
    - `dgen_sim_HLA_coeff_3000_leaves_sims_average.Rdata`: Using all available reference sequence data from public databases.
* **Figure 8**: Simulation Study 2: impact of reduced reference sequence set, results summary.
  - .Rdata files containing _p_ values, medians, or estimated selection coefficient information used to generate ROC curves: `Rdata_files/simulation_study_2/dgen_1460_leaves/*`. The separate files contain results for the different methods and different subsets of the reference dataset:
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
  - .Rdata files containing quantile information (2.5%, 25%, 50%, 75%, 97.5%) for simulation study 3:
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
  - TO DO
* **Figure 14**: Summary of HLA-A associated selection in protease. 
  - Raw output from MCMC run: `mcmc_output_files/protease_combined_subtype_{B,C}_ref_combined_{window, log}.txt.zip`.
  - .Rdata file, in nicer array format, with HLA labels and removal of burn-in: `Rdata_files/HLA_inference/reverse_transcriptase_combined_subtype_{B,C}_ref.Rdata`.
* **Figure 15**: Summary of HLA-B associated selection in protease. 
  - Raw output from MCMC run: `mcmc_output_files/protease_combined_subtype_{B,C}_ref_combined_{window, log}.txt.zip`.
  - .Rdata file, in nicer array format, with HLA labels and removal of burn-in: `Rdata_files/HLA_inference/reverse_transcriptase_combined_subtype_{B,C}_ref.Rdata`.
* **Figure 16**: Summary of HLA-C associated selection in protease. 
  - Raw output from MCMC run: `mcmc_output_files/protease_combined_subtype_{B,C}_ref_combined_{window, log}.txt.zip`.
  - .Rdata file, in nicer array format, with HLA labels and removal of burn-in: `Rdata_files/HLA_inference/reverse_transcriptase_combined_subtype_{B,C}_ref.Rdata`.
* **Figure 17**: Summary of HLA-A associated selection in reverse transcriptase.
  - Raw output from MCMC run: `mcmc_output_files/reverse_transcriptase_combined_subtype_{B,C}_ref_combined_{window, log}.txt.zip`.
  - .Rdata file, in nicer array format, with HLA labels and removal of burn-in: `Rdata_files/HLA_inference/reverse_transcriptase_combined_subtype_{B,C}_ref.Rdata`.
* **Figure 18**: Summary of HLA-C associated selection in reverse transcriptase.
  - Raw output from MCMC run: `mcmc_output_files/reverse_transcriptase_combined_subtype_{B,C}_ref_combined_{window, log}.txt.zip`.
  - .Rdata file, in nicer array format, with HLA labels and removal of burn-in: `Rdata_files/HLA_inference/reverse_transcriptase_combined_subtype_{B,C}_ref.Rdata`.
* **Figure 19**: Estimated recombination probabilities between sites across protease and reverse transcriptase.
  - `mcmc_output_files/{protease,reverse_transcriptase}_combined_subtype_{B,C}_ref_combined_log.txt.zip`.
* **Figure 20**: The impact of HLA alleles on differentiation between HIV-1 subtypes.
  - TO DO

### Tables
* **Table 1**: Priors on parameters of the model.
  - N/A
* **Table 2**: The collection of ‘false positive’ associations are interrogated.
  - N/A
* **Table 3**: The collection of apparent false negative associations are interrogated. 
  - N/A
* **Table 4**: The collection of top-tier and second-tier sites in protease.
  - Raw output from MCMC run: `mcmc_output_files/protease_combined_subtype_{B,C}_ref_combined_{window, log}.txt.zip`.
  - .Rdata file, in nicer array format, with HLA labels and removal of burn-in: `Rdata_files/HLA_inference/protease_combined_subtype_{B,C}_ref.Rdata`.
  - Summary information containing the median, lower 2.5% and 10% quantile information: TO DO.
* **Table 5**: The collection of top-tier and second-tier sites in reverse transriptase.
  - Raw output from MCMC run: `mcmc_output_files/reverse_transcriptase_combined_subtype_{B,C}_ref_combined_{window, log}.txt.zip`.
  - .Rdata file, in nicer array format, with HLA labels and removal of burn-in: `Rdata_files/HLA_inference/reverse_transcriptase_combined_subtype_{B,C}_ref.Rdata`.
  - Summary information containing the median, lower 2.5% and 10% quantile information: TO DO.
* **Table 6**: Overlap between studies, protease.
  - Summary information containing the median, lower 2.5% and 10% quantile information: TO DO.
  - Also, summarise where I got the Carlson information from.
* **Table 7**: Overlap between studies, reverse transcriptase.
  - Summary information containing the median, lower 2.5% and 10% quantile information: TO DO.
  - Also, summarise where I got the Carlson information from.
* **Table 8**: Table of _p_-values for differences between _H_<sub>B</sub> and _H_<sub>C</sub>) associated  selection at sites that distinguish subtype B and subtype C in protease and reverse transcriptase.
  - TO DO
* **Table 9**: Summary of reference datasets for simulation study 3.
  - As in Supplementary Figure 9.


