import requests
import sys
import json
import jsonlines
from langchain_openai import ChatOpenAI
from langchain_core.prompts import ChatPromptTemplate

def search_projects_by_keyword(keyword):
    # Encode the keyword with %20 for spaces
    encoded_keyword = "%20".join(keyword.split())
    
    base_url = "https://www.ebi.ac.uk/pride/ws/archive/v2/"
    search_url = f"{base_url}search/projects?keyword={encoded_keyword}&filter=project_submission_type%3D%3DPARTIAL%2Cproject_submission_type%3D%3DCOMPLETE&pageSize=100&sortDirection=DESC&sortFields=submission_date"
    
    response = requests.get(search_url)
    if response.status_code == 200:
        response_data = response.json()
        total_projects = response_data['page']['totalElements']
        print(f"Total number of projects related to '{keyword}': {total_projects}")
        
        projects_data = []
        
        for project in response_data['_embedded']['compactprojects']:
            project_details = {
                "project_title": project['title'],
                "project_accession": project['accession'],
                "project_description": project.get('projectDescription', 'N/A'),
                "sample_processing_protocol": project.get('sampleProcessingProtocol', 'N/A'),
                "data_processing_protocol": project.get('dataProcessingProtocol', 'N/A'),
                "keywords": project.get('keywords', []),
                "softwares": project.get('softwares', []),
                "instruments": project.get('instruments', []),
                "organisms": project.get('organisms', []),
                "organism_parts": project.get('organismParts', [])
            }
            # Fetch and append file URLs for each project
            project_details['file_urls'] = fetch_ftp_urls(project['accession'])
            
            projects_data.append(project_details)
        
        write_to_jsonl(projects_data)  # Write project data to JSONL file
        
    else:
        print("Failed to retrieve data")

def fetch_ftp_urls(project_accession):
    files_response = requests.get(f"https://www.ebi.ac.uk/pride/ws/archive/v2/files/byProject?accession={project_accession}")
    if files_response.status_code == 200:
        files_data = files_response.json()
        ftp_urls = []
        for file_info in files_data:
            file_category = file_info.get('fileCategory', {}).get('value', '').lower()
            if file_category in ['raw', 'mzml']:
                for location in file_info.get('publicFileLocations', []):
                    if location.get('name') == 'FTP Protocol':
                        ftp_urls.append(location.get('value'))
        ftp_urls.sort()
        return ftp_urls
    else:
        print("Failed to retrieve files data")
        return []

def write_to_jsonl(data):
    with open('output.jsonl', 'w') as f:
        for item in data:
            json.dump(item, f)
            f.write('\n')

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script_name.py keyword")
        sys.exit(1)
    
    keyword = sys.argv[1]
    search_projects_by_keyword(keyword)

llm = ChatOpenAI(model='gpt-4-turbo',api_key="add-key")

prompt = ChatPromptTemplate.from_messages([
    ("system", "You are a mass spectrometry and proteomics expert. You are tasked to determine mass spec parameters based off of experimental descriptions" +
     '''Summarize information from the following descriptions in list form for use in mass spectrometry database searches. 

Project Description: Circular RNAs (circRNAs) are covalently closed non-coding RNAs lacking the 5’ cap and the poly-A tail. Nevertheless, it has been demonstrated that certain circRNAs can undergo active translation. Therefore, aberrantly expressed circRNAs in human cancers could be an unexplored source of tumor-specific antigens, potentially mediating anti-tumor T cell responses. This study presents an immunopeptidomics workflow with a specific focus on generating a circRNA-specific protein fasta reference. The main goal of this workflow is to streamline the process of identifying and validating human leukocyte antigen (HLA) bound peptides potentially originating from circRNAs. We increased the analytical stringency of our workflow by retaining peptides identified independently by two mass spectrometry search engines and/or by applying a group-specific FDR for canonical-derived and circRNA-derived peptides. A subset of circRNA-derived peptides specifically encoded by the region spanning the back-splice junction (BSJ) were validated with targeted MS, and with direct Sanger sequencing of the respective source transcripts. Our workflow identified 54 unique BSJ-spanning circRNA-derived peptides in the immunopeptidome of melanoma and lung cancer samples. Our novel approach enlarges the catalog of source proteins that can be explored for immunotherapy.
Sample Processing Protocol: Purification of HLA-I peptides and LC-MS/MS Analysis: HLA-I immunoprecipitation was performed according to our previous established protocols6, 64. Shortly, protein-A Sepharose 4B (Pro-A) beads (Invitrogen) were used to purify W6/32 monoclonal antibodies from the supernatant of HB95 hybridoma cells (ATCC HB-95). After antibody crosslinking, Pro-A beads were used for immunoaffinity purification of HLA-I complexes from tissue or cell line lysates. HLA-I peptides were then purified using a C18 solid phase extraction (SPE) and dried using vacuum centrifugation (Concentrator plus, Eppendorf). Samples were stored at -80°C if not immediately submitted to mass spectrometry analysis. Finally, immunopeptides were re-suspended in 2% ACN and 0.1% FA (formic acid). iRT peptides (Biognosis, Schlieren, Switzerland) were spiked into in the samples (Biognosis) as indicated in Supplementary Table 7 and analyzed by LC-MS/MS. :Liquid chromatography and mass spectrometry (LC-MS): The LC-MS system consisted of an Easy-nLC 1200 coupled to Q Exactive HF-X mass spectrometer (ThermoFisher scientific, Bremen, Germany) or to Eclipse tribrid mass spectrometer (ThermoFisher Scientific, San Jose, USA). The peptides were eluted on a 450mm analytical column (8 m tip, 75 m ID) packed with ReproSil-Pur C18 (1.9 m particle size, 120 A pore size, Dr Maisch, GmbH) and separated at a flow rate of 250 nL/min as described64. For DDA measurements, the top 20 most abundant precursor ions selection was performed on the Q Exactive as described64. For DIA, the Eclipse tribrid mass spectrometer was used to samples ions. The cycle of acquisitions consists of a full MS scan from 300 to 1650 m/z (R = 120,000, ion accumulation time of 60 ms and normalized AGC of 250%) and 22 DIA MS/MS scans in the orbitrap. For each DIA MS/MS scan, a resolution of 30,000, a normalized AGC of 250%, and a stepped normalized collision energy (27, 30, and 32) were used. The maximum ion accumulation was set to auto, the fixed first mass was set to 200 m/z, and the overlap between consecutive MS/MS scans was 1 m/z as described79.
Data Processing Protocol: DDA MS search with group-specific FDR: MS-derived raw files resultant from three biological replicates of T1185B (Supplementary Table 7) were searched using MaxQuant51 (version 2.1.0.0) with a PSM FDR of 0.1 and Comet52 against the generated reference fasta file containing both UniProt and the trimmed circRNA-derived putative ORFs around the BSJ and initiated by the canonical start codon ATG. Outputs were then intersected by NewAnce (v1.6) setting a group-specific FDR of 0.03 for protein- and circRNA-derived peptides. The search was done setting a nonspecific protein digestion cleavage, no fixed modifications, methionine oxidation and protein N-term acetylation as variable modifications, and restricting the peptide length to 8-15 amino acids. Same approach was used for the cell line and tumor tissues of patient Mel-1 patient (two biological replicates and three different lymph node regions, respectively). Hybrid DIA MS search: DIA files corresponding to the immunopeptidomics samples of MG132 and IFN treatment of T1185B cells were searched using a hybrid DIA approach using two computational tools, FragPipe (v.20.0) with group-specific FDR calculation57, 60, 82 and Spectronaut (v.18.4)83, against the generated reference fasta file containing both UniProt and the trimmed circRNA-derived ORFs around the BSJ to which a list of common MS contaminant proteins were added. To increase the coverage of the spectral library, we assembled available DDA raw files from T1185B cells treated or not with IFN (from the PRIDE accession PXD0136496), the newly generated DDA data, together with the DIA files of T1185B cells treated with MG132 and the new IFN treatments (and their respective controls), as indicated in Supplementary Table 7. In FragPipe we applied a group-specific FDR threshold of 0.03 (MSFragger Group variable: Protein evidence from FASTA file) while in Spectronaut we applied a global peptide FDR threshold of 0.01. In both engines, the search was done by applying a FDR threshold of 1 for proteins, nonspecific protein digestion cleavage, no fixed modifications, methionine oxidation and protein N-term acetylation as variable modifications, and restricting the peptide length to 8-15 amino acids. Hybrid spectral libraries were then used to match and quantify peptides from the immunopeptidomics DIA data using a peptide precursor group-specific FDR of 0.03 for FragPipe or global FDR of 0.01 for Spectronaut. Default decoy generation methods were used for each MS search tool, reversed and mutated sequences for FragPipe and Spectronaut, respectively. Data analysis was performed using Fragpipe quantification values after overlapping identified sequences from both FragPipe and Spectronaut MS analysis tools. DDA-DIA MS search with group-specific FDR calculation: HLA-I and HLA-II raw files of the lung cancer cohort of 8 patients and 52 tumoral and healthy matched tissues62 (Supplementary Table 7) were downloaded from PRIDE PXD034772 and analyzed by FragPipe (v.19.2-build39 for HLA-I and v.20.0 for HLA-II immunopeptidomes) which supported group-specific FDR calculation. Spectral library generation was performed using the DDA immunopeptidomics data. The search was done setting a nonspecific protein digestion cleavage, no fixed modifications, methionine oxidation and protein N-term acetylation as variable modifications, a group-specific FDR threshold of 0.03 for peptides (MSFragger Group variable: Protein evidence from FASTA file), a FDR threshold of 1 for proteins, and restricting the peptide length to 8-15 or to 8-25 amino acids for HLA-I or HLA-II MS searches, respectively. Respective spectral libraries were then used to match and quantify peptides from the immunopeptidomics DIA data using a peptide precursor group-specific FDR of 0.03. DIA immunopeptidomics raw files were used for matching and quantification of peptides. Peptides from canonical and non-canonical circRNA groups were used to calculate the FDR separately because the score distributions are different. Pooling them together would result in underestimated FDR for the circRNA group. HLA-I and HLA-II library generation and peptide identification were performed in independent analysis.

Sample Processing in-brief:
- W6/32 antbody was used for immunoaffinity purification of HLA-I complexes
- Steps include C18 cleanup, drying, and resuspension
Mass spectrometer: Q Exactive HF-X mass spectrometer or Eclipse tribrid mass spectrometer
DDA processing: Q Exactive instrument, raw files resultant from three biological replicates of T1185B and tumor tissues of patient Mel-1 patient (two biological replicates and three different lymph node regions, respectively)
- Software: MaxQuant v2.1 and Comet
- Reference: both UniProt and the trimmed circRNA-derived putative ORFs around the BSJ and initiated by the canonical start codon ATG
- PSM_FDR: 0.1
- group-specific FDR: 0.03 for protein- and circRNA-derived peptides(NewAnce)
- protein digetsion: nonspecific
- Fixed modifications: none
- Variable modifications: methionine-oxidation and n-term acetylation
- Maximum peptide length = 15 
- Minimun peptide length = 8 
Hybrid-DIA processing: Eclipse tribrid, MG132 and IFN treatment of T1185B cells
- Software: FragPipe v20 (group-specific FDR) and Spectronaut
- Reference: both UniProt and the trimmed circRNA-derived putative ORFs around the BSJ and initiated by the canonical start codon ATG (with contaminants added)
- Spectral Library: available DDA raw files from T1185B cells treated or not with IFN (from the PRIDE accession PXD0136496), newly generated DDA data, together with the DIA files of T1185B cells treated with MG132 and the new IFN treatments (and their respective controls)
- FragPipe group-specific FDR threshold of 0.03 (MSFragger Group variable: Protein evidence from FASTA file) 
- Spectronaut global peptide FDR threshold of 0.01. 
- Protein-FDR threshold = 1
- protein digetsion: nonspecific
- Fixed modifications: none
- Variable modifications: methionine-oxidation and n-term acetylation
- Maximum peptide length = 15 
- Minimun peptide length = 8 
DDA-DIA MS search with group-specific FDR calculation: HLA-I and HLA-II raw files of the lung cancer cohort of 8 patients and 52 tumoral and healthy matched tissues  were downloaded from PRIDE PXD034772
- Software: FragPipe (v.19.2-build39 for HLA-I and v.20.0 for HLA-II immunopeptidomes)
- Spectral Library: DDA data
- FragPipe group-specific FDR threshold of 0.03 (MSFragger Group variable: Protein evidence from FASTA file) 
- protein digetsion: nonspecific
- Protein-FDR threshold = 1
- protein digetsion: nonspecific
- Fixed modifications: none
- Variable modifications: methionine-oxidation and n-term acetylation
- HLA-I:
- Maximum peptide length = 15 
- Minimun peptide length = 8 
- HLA-II:
- Maximum peptide length = 25 
- Minimun peptide length = 8 

Summarize information from the following descriptions in list form for use in mass spectrometry database searches. 

Project Description: Background: Patient derived organoids (PDOs) can be established from colorectal cancers as in vitro models to interrogate cancer biology and its clinical relevance. We applied mass spectrometry (MS) immunopeptidomics to investigate neoantigen presentation and whether this can be augmented through interferon gamma (IFN) or MEK-inhibitors. Methods: Four PDOs from chemotherapy refractory and one from a treatment naïve CRC were expanded to replicates with 100 million cells each, and HLA class I and class II peptide ligands were analysed by MS. Results: We identified an average of 9,936 unique peptides per PDO which compares favourably against published immunopeptidomics studies, suggesting high sensitivity. Loss of heterozygosity of the HLA locus was associated with low peptide diversity in one PDO. Peptides from genes without detectable expression by RNA-sequencing were rarely identified by MS. Only 3 out of 612 non-silent mutations encoded for neoantigens that were detected by MS. Treatment of four PDOs with IFN upregulated HLA class I expression and qualitatively changed the immunopeptidome, with increased presentation of IFN-inducible genes. HLA class II presented peptides increased on average 16-fold with IFN treatment. MEK-inhibitor treatment showed no consistent effect on class I or II HLA expression or the peptidome. Importantly, no additional HLA class I or II presented neoantigens became detectable with any treatment.  Conclusions: Only 3 out of 612 non-synonymous mutations encoded for neoantigens that were detectable by MS. Although MS has sensitivity limits and biases, and likely underestimated the true neoantigen burden, this established a lower bound of the percentage of non-silent mutations that encode for presented neoantigens, which may be as low as 0.5%. This could be a reason for the poor responses of non-hypermutated CRCs to immune checkpoint inhibitors. MEK-inhibitors recently failed to improve checkpoint-inhibitor efficacy in CRC and the observed lack of HLA upregulation or improved peptide presentation may explain this.
Sample Processing Protocol: Purification of HLA-I and HLA-II peptides Anti-pan-HLA-I and anti-HLA-II monoclonal antibodies were purified from the supernatant of HB95 (ATCC® HB-95™) and HB145 cells (ATCC® HB-145™), respectively, grown in CELLLine CL-1000 flasks (Sigma-Aldrich, St. Louis, Missouri) using Protein A-Sepharose 4B beads (Pro-A beads; Invitrogen, Carlsbad, California). Antibodies were cross-linked to Pro-A beads at a concentration of 5 mg of antibodies per 1 mL volume of beads with Dimethyl pimelimidate dihydrochloride (Sigma-Aldrich) in 0.2 M Sodium Borate buffer pH 9 (Sigma-Aldrich) at a final concentration of 20 mM for 30 minutes. The reaction was quenched by incubation with 0.2 M ethanolamine pH 8 (Sigma-Aldrich) for 2 hours. Cross-linked antibodies were kept at 4°C until use.  PDO cells were lysed in biological replicates of 100 million cells at 4°C for 1 hour in PBS containing 0.25% sodium deoxycholate (Sigma-Aldrich), 0.2 mM iodoacetamide (Sigma-Aldrich), 1 mM EDTA, 1:200 Protease Inhibitors Cocktail (Sigma-Aldrich), 1 mM Phenylmethylsulfonylfluoride (Roche, Basel, Switzerland), 1% octyl-beta-D glucopyranoside (Sigma-Alrich). Cell lysates were cleared by centrifugation with a table-top centrifuge (Eppendorf Centrifuge, Hamburg, Germany) at 4°C at 14,200 rpm for 50 min. For the sequential purification of HLA-I and HLA-II, two stacked 96-well single-use micro-plates (3µm glass fiber, 10µm polypropylene membranes; ref number 360063, Seahorse Bioscience, North Billerica, Massachusetts) were used. The first contained cross-linked beads with anti HLA-I antibodies and the second contained the anti HLA-II corss-linked beads. The lysates were loaded on plates by gravity. Then the Waters Positive Pressure-96 Processor (Waters, Milford, Massachusetts) was employed and the plates were washed separately with 4 times 2 mL of 150 mM sodium chloride (NaCl) (Carlo-Erba, Val de Reuil, France) in 20 mM Tris-HCl pH 8, 4 times 2 mL of 400 mM NaCl in 20 mM Tris-HCl pH 8 and again with 4 times 2 mL of 150 mM NaCl in 20 mM Tris-HCl pH 8. Finally, the plates were washed twice with 2 mL of 20 mM Tris-HCl pH 8. Each plate was stacked on top of a Sep-Pak tC18 100 mg Sorbent 96-well plate (ref number: 186002321, Waters) equilibrated with 1 mL of 80% acetonitrile (ACN) in 0.1 % TFA and with 2 mL of 0.1% TFA. The HLA and peptides were eluted with 500 µL 1% TFA into the Sep-Pak plate and then the Sep-Pak plates were washed with 2 mL of 0.1 % TFA. HLA-I and HLA-II peptides were eluted with 500 µL of 32% ACN in 0.1% TFA into a collection plate. Recovered peptides were dried using vacuum centrifugation (Concentrator plus Eppendorf) and stored at -20°C. LC-MS/MS analyses Prior to MS analysis, HLA-I and HLA-II peptide samples were re-suspended in 9 µL of 0.1 % FA  and 2/3 of the sample volume were placed in the Ultra HPLC autosampler. Each sample was measured in technical duplicates. HLA peptides were separated by a nanoflow HPLC (Proxeon Biosystems, Thermo Fisher Scientific, Odense) coupled on-line to a Q Exactive HFX mass spectrometers (Thermo Fisher Scientific, Bremen) with a nanoelectrospray ion source (Proxeon Biosystems). Self-packed 50 cm long (75 μm inner diameter) column with ReproSil-Pur C18-AQ 1.9 μm resin (Dr. Maisch GmbH, Ammerbuch-Entringen, Germany) in buffer A (0.5% acetic acid). HLA-I and HLA-II peptides were eluted with a linear gradient of 2–30% buffer B (80% ACN and 0.5% acetic acid) at a flow rate of 250 nl/min over 125 min  and 90 min, respectively. MS spectra were acquired from m/z = 300-1650 in the Orbitrap with a resolution of 60’000 (m/z = 200) and ion accumulation time of 80 ms. The auto gain control (AGC) was set to 3e6 ions. MS/MS spectra were acquired on 10 most abundant precursor ions (if present) with a resolution of 15,000 (m/z = 200), ion accumulation time of 120 ms and an isolation window of 1.2 m/z. The AGC was set to 2e5 ions, dynamic exclusion to 20 s and a normalized collision energy (NCE) of 27 was used for fragmentation. The peptide match option was disabled. No fragmentation was performed for HLA-I peptides in case of assigned precursor ion charge states of four and above, and for HLA-II peptides, in case of assigned precursor ion charge states of one, and also from six and above.
Data Processing Protocol: Peptide identification and Statistical analysis We searched the immunopeptidomics peak lists data with a PSM false discovery date of 1% with the MaxQuant platform (PMID 19029910) version 1.5.5.1 against a fasta file containing the human proteome (Homo_sapiens_UP000005640_9606, the reviewed part of UniProt, with no isoforms, including 21,026 entries downloaded in March 2017) and a  list of 247 frequently observed contaminants. In addition, MaxQuant version 1.5.9.4 was used for a search against a customized reference database were somatic mutations were included, with a PSM false discovery rate of 5%. Peptides with a length between 8 and 25 AA were allowed. The second peptide identification option in Andromeda was enabled. The enzyme specificity was set as unspecific and no protein FDR was set. The initial allowed mass deviation of the precursor ion was set to 6 ppm and the maximum fragment mass deviation was set to 20 ppm. Methionine oxidation and N-terminal acetylation were set as variable modifications. The ‘match between runs’ analysis was set separately for all replicates per PDO, and separately for HLA class I and HLA class II samples. For analysis of unique identified peptide sequences, we utilized a simple binary criterion of present or absent. A peptide was only defined as present if it was detected in both technical replicates. A peptide was then defined as present if it was detected in any of the biological replicates, to give the highest sensitivity of detection. Venn diagrams were generated by Venny 2.1.0. The raw MS intensity values were Log2 transformed. For differential expression analyses the Perseus computational platform (PMID 27348712) was used for ‘width normalization’, and missing values were imputed.

Sample Processing in-brief:
- Anti-pan-HLA-I and anti-HLA-II monoclonal antibodies were used for immunoaffinity purification of HLA I and II complexes
- Steps include C18 cleanup, drying, and resuspension
Mass spectrometer: Q Exactive HFX mass spectrometers
Data processing:
- Software: MaxQuant
- PSM_FDR: 0.01 for normal reference and 0.05 for custom database
- Reference: human proteome (Homo_sapiens_UP000005640_9606, the reviewed part of UniProt, with no isoforms, including 21,026 entries downloaded in March 2017) and a  list of 247 frequently observed contaminants and customized reference database with somatic mutations
- Mazimum peptide length = 25 
- Minimum peptide length = 8
- protein digetsion: nonspecific
- precursor mass tolerance: 6 ppm
- fragment mass tolerane: 20pm
- Fixed modifications: none
- Variable modifications: methionine-oxidation and n-term acetylation

Summarize information from the following descriptions in list form for use in mass spectrometry database searches. 

Project Description: Immune checkpoint inhibitor and adoptive lymphocyte transfer-based therapies have shown great therapeutic potential for cancers with high tumor mutation burden (TMB). Here, we employed mass spectrometry (MS)-based proteogenomic large-scale profiling to identify potential immunogenic human leukocyte antigen (HLA) Class ǀ-presented peptides in both melanoma, a “hot tumor” with high TMB, and EGFR mutant lung adenocarcinoma, a “cold tumor” with low TMB. We used cell line and patient-specific databases constructed using variants identified from whole-exome sequencing, as well as a de novo search algorithm from the PEAKS search algorithm to interrogate the mass spectrometry data of the Class I immunopeptidome. We identified 12 mutant neoantigens. Several classes of tumor-associated antigen-derived peptides were also identified. We constructed a cancer germline (CG) antigen database with 285 antigens and identified 42 Class I-presented CG antigens. We identified more than 1000 post-translationally modified (PTM) peptides representing 58 different PTMs. Our results suggest that PTMs play a critical role impacting HLA-binding affinity dramatically.  Finally, leveraging de novo search and an annotated lncRNA database, we developed a novel non-canonical peptide discovery pipeline to identify 44 lncRNA-derived peptides that are presented by Class I. We validated MS/MS spectra of select peptides using synthetic peptides and performed HLA Class I binding assays to demonstrate binding of select neo-peptides and lncRNA-derived peptides to Class I proteins. In summary, we provide direct evidence of HLA Class I presentation of a large number of mutant and tumor-associated peptides for potential use as vaccine or adoptive cell therapy in melanoma and EGFR mutant lung cancer.
Sample Processing Protocol: For melanoma and lung adenocarcinoma cell line HLA peptide enrichment, 2.0 × 108 cells were harvested in 4ml ice-cold lysis buffer (20mM Tris-HCl pH=8.5, 100mM NaCl, 1 mM EDTA. 1% triton X-100 supplemented with Halt 1:100 protease Inhibitor cocktail Cat. No 78430, Thermo Scientific). After 30min on ice, lysates were subjected to needle sonication for 30s. Approximately 30mg snap-frozen lung tumor tissue was homogenized in 4ml ice-cold lysis buffer for 30s at 4°C using the Qiagen TissueLyser II. Cell/tissue lysates were centrifuged for 2 hours at 4°C at 20,000g, and the supernatant used in subsequent experiments. HLA-peptides complexes were isolated by interacting with 0.5mg HLA Class ǀ pan antibody clone W6/32 (BioXcell, West Lebanon, NH) pre-coupled to 200μl slurry of protein A/G PLUS agarose resin (Santa Cruz Biotechnology) overnight at 4°C with constant rotation. Agarose beads were then washed three times with ice-cold lysis buffer (without triton and protease inhibitors), followed by two washes in ice-cold 20mM Tris-HCl (pH=8.5), then one wash in ice-cold HPLC grade water. Complexes were eluted four times with 0.15% trifluoroacetic acid (TFA) in water at room temperature and combined. To purify immunopeptides, HLA-peptide complexes were loaded on preconditioned 50mg C18 desalting columns (Sigma Millipore), then followed three 0.1% TFA in water washes. HLA peptides were then eluted with 40% acetonitrile (ACN) in 0.1% TFA. Purified peptides were lyophilized at -80°C for 2hours, then, reconstituted in 0.1%TFA, 2%ACN loading buffer for MS analysis.For the whole-cell proteome of H1975, Roughly 0.2mg of the cell lysate of H1975 were reduced and alkylated, and further digested by MS grade trypsin/lysC at 37°C for 16 hours. The tryptic peptides were fractionated by an off-line high-pH reverse phase liquid chromatography (LC) into 12 fractions, each of which was subjected to a 120 min gradient separation on a nanoLC and analyzed by an Orbitrap Q-Exactive HF using discovery mode.
Data Processing Protocol: Patient and cell line specific protein sequence databases were first generated. The aligned whole exome sequencing BAM files from patient blood (germline) or tumor tissue/cell lines were used to retrieve the variants call format (VCF) files using HaplotypeCaller (Ren et al., 2019), and the intermediate VCF files were further annotated by SnpEff (Cingolani et al., 2012)which filtered out only nonsynonymous on exome regions including SNVs and Indels. Similarly, BBduk was used to remove adapter sequences and low-quality reads from paired-end FASTQ files, which were then used as input for STAR-Fusion (Haas et al., 2019). The final VCF files were in silico translated to sample specific protein sequence libraries (FASTA files) using QUILTS (Ruggles et al., 2016), which were merged with refseq hg38 converted human proteome database.  The database search of MS raw files was carried out by PEAKS studio (Tran et al., 2019) (Bioinformatics Solutions) using the patient/cell line specific databases described above. In the PEAKS searching engine, no enzyme digestion was selected because HLA peptides are natural peptides without artificial digestion. Importantly, the unique PEAKS built-in functions, pan-PTMs including 650 different variable modifications and de novo searching were used. The precursor mass tolerance was set to 15 ppm and fragment ion tolerance was set to 0.5 Da. For the mutant neoepitope screening, the false discovery rate (FDR), estimated by decoy-fusion database, was chosen at 0.05. For the wt peptides screening (e.g., CG antigen, LncRNA and PTMs), FDR was set to 0.01. For the de novo search, we applied very stringent criteria, a) the average local confidence (ALC) score (substituted for FDR estimation) of each peptides must be >50%; b) the lowest %rank, predicted by NetMHCpan4.0, of each peptides against their corresponding HLA  molecules in the same sample must be <2.0. The detected raw peptide intensity was log2 transformed for further statistical analysis.

Sample Processing in-brief:
Immunopeptidome: melanoma and lung adenocarcinoma cell line
- W6/32 antbody was used for immunoaffinity purification of HLA-I complexes
- Steps include C18 cleanup, drying, and resuspension
Whole protoeme: cell lysate of H1975
- reduction/alkylation
- digestion with trypsin/LysC
- 12 fractions
Mass spectrometer: Orbitrap Q-Exactive HF using discovery mode
Data processing:
- Software: PEAKS
- Reference: human proteome merged with patient/cell specific databases
- protein digetsion: none
- pan-PTMs including 650 different variable modifications 
- de novo searching 
- precursor mass tolerance: 15ppm
- fragment mass tolerane: 0.5Da
- PSM_FDR: 0.01 for normal reference and 0.05 for custom database

Summarize information from the following descriptions in list form for use in mass spectrometry database searches. 

Project Description: Effective immunosurveillance of cancer requires the presentation of peptide antigens on major histocompatibility complex Class I (MHC-I). Recent developments in proteomics have improved the identification of peptides that are naturally presented by MHC-I, collectively known as the “immunopeptidome”. Current approaches to profile tumor immunopeptidomes have been limited to in vitro investigation, which fails to capture the in vivo repertoire of MHC-I peptides, or bulk tumor lysates, which are obscured by the lack of cancer cell-specific MHC-I isolation. To overcome these limitations, we report here the engineering of a Cre recombinase-inducible affinity tag into the endogenous mouse MHC-I gene and targeting of this allele to the KrasLSL-G12D/+; p53fl/fl (KP) mouse model (KP/KbStrep). This novel approach has allowed us to precisely isolate MHC-I peptides from autochthonous pancreatic ductal adenocarcinoma (PDAC) and lung adenocarcinoma (LUAD) in vivo. With this powerful analytical tool, we were able to profile the evolution of the LUAD immunopeptidome from the alveolar type 2 cell-of-origin through late-stage disease. Differential peptide presentation in LUAD is not driven by increased mRNA expression or translation rate and is likely driven by post-translational mechanisms. Vaccination of mice with peptides presented by LUAD in vivo provoked CD8 T-cell responses in naïve and tumor bearing mice. Many peptides unique to LUAD, including immunogenic peptides, exhibited very low expression of the cognate mRNA provoking reconsideration of antigen prediction pipelines that triage peptides according to transcript abundance. Beyond cancer, the KbStrep allele is compatible with a broad range of Cre-driver lines to explore antigen presentation in vivo in the pursuit of understanding basic immunology, infectious disease, and autoimmunity.
Sample Processing Protocol: pMHCs were isolated by affinity purification and purified prior to LC-MS/MS analysis. Exact methodological details can be found in the associated manuscript.
Data Processing Protocol: All mass spectra were analyzed with Proteome Discoverer (PD, version 2.5) and searched using Mascot (version 2.4) against the mouse SwissProt database (2021_03; 2021_02 for label free quantification) supplemented with a list of murine ORFs previously identified by ribosome profiling (www.sorfs.org). Peptides were searched with no enzyme and variable methionine oxidation. Peptide spectrum matches (PSMs) were further filtered according to the following criteria: ions score ≥ 15, search engine rank = 1 and results from technical replicates of each sample analysis were combined. Median retention time (RT) was calculated using the RT values of filtered PSMs from all replicates of a given sample. Label free quantitation was done using the Minora Feature Detector (precursor abundance values measured based on area under the curve) in PD with match between runs enabled and filtered for peptides with ions score ≥ 15 and search engine rank = 1. Abundances were averaged across technical and biological replicates. The lot number used for TMT 6plex (Thermo fisher) was UG285371.

Sample Processing in-brief:
- affinity purification of HLA-I complexes
- see manuscript for more info
Data processing:
- Software: Proteome Discoverer v2.5 and Mascot v 2.4
- Reference: mouse SwissProt database merged with ORFs from ribosome profiling
- protein digetsion: none
- Fixed modifications: none
- Variable modifications: methionine oxidation

Summarize information from the following descriptions in list form for use in mass spectrometry database searches. 

Project Description: T cell recognition of human leukocyte antigen (HLA)-presented tumor-associated peptides is central for cancer immune surveillance. Mass spectrometry (MS)-based immunopeptidomics represents the only unbiased method for directly identifying and characterizing naturally presented tumor-associated peptides, which represent a key prerequisite for the development of T cell-based immunotherapies. This study reports on the de novo implementation of ion mobility separation-based timsTOF MS for next-generation immunopeptidomics, enabling high-speed and sensitive detection of HLA peptides. A direct comparison of timsTOF-based with state-of-the-art immunopeptidomics using orbitrap technology showed significantly increased HLA peptide identifications from benign and malignant primary samples of solid tissue and hematological origin. First application of timsTOF immunopeptidomics for tumor antigen discovery enabled (i) the expansion of benign reference immunopeptidome databases with more than 150,000 HLA peptides from 94 primary benign tissue samples, (ii) the refinement of previously described tumor antigens, and (iii) the identification of a vast array of novel tumor antigens comprising low abundant neoepitopes that might serve as targets for future cancer immunotherapy development.
Sample Processing Protocol: HLA class I and II molecules were isolated using standard immunoaffinity purification as described (DOI: 10.1007/978-1-62703-218-6_12) using the pan-HLA class IÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Å“specific mAb W6/32, the pan-HLA class IIÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Å“specific mAb TuÃƒÅ’Ã‹â€ 39, and the HLA-DRÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Å“specific mAb L243 (produced in-house). HLA ligand extracts were analyzed in 3 technical replicates using the timsTOF Pro device (Bruker). In brief, samples were separated using high-performance liquid chromatography (nanoElute, Bruker) using an acclaim TM PepMap (thermo Fisher Scientific) and a 75Ãƒâ€šÃ‚Âµmx25cm Aurora emitter column (IonOpticks) and a gradient ranging from 0-60% over 55 min folowed by 5 min 95% acetonitril with 0.01% formic acid (total 60 min). Eluting peptides were analyzed on-line timsTOF Pro in PASEF mode and collision induced dissociation. Ten samples were additionally analyzed using prevously described methods(DOI: 10.1073/pnas.1416389112). In brief, peptide samples were separated by nanoflow high-performance liquid chromatography (Ultimate 3000 RSLC Nano UHPLC; Thermo Fisher Scientific) using a 75 ÃƒÅ½Ã‚Â¼m ÃƒÆ’Ã¢â‚¬â€� 2 cm PepMap Nanotrap and a 50 Ãƒâ€šÃ‚Âµm x 25 cm PepMap RSLC C18 column (Thermo Fisher Scientific) and a gradient ranging from 2.4% to 32.0% acetonitrile over the course of 90 minutes. Eluting peptides were analyzed in an on-line coupled Orbitrap Fusion Lumos mass spectrometer (Thermo Fisher Scientific) employing a top-speed collisional-induced dissociation or higher-energy collisional dissociation fragmentation.
Data Processing Protocol: Data processing was performed using Bioinformatics Solutions Inc. PEAKS Studio 10.6. Three DDA runs, i.e. three technical replicates, were co-processed. For precursor mass an error tolerance of 20.0 ppm (monoisotopic mass) and 0.02 for fragment ions was allowed. For enzymatic digest, "none" was selected in unspecific digest mode. As a variable PTM, methionine oxidation (+15.9949) was selected with maximum 3 PTMs allowed per peptide. As database the human reference proteome datase was selected/ downloaded prior from Uniprot (14.10.2020). A False-discovery-rate was estimated with a decoy-fusion. Filters post-processing included: 1% FDR for Peptides with a PTM Ascore >20; -10lgP 0 abd 0 unique peptides with significant peptides for proteins; De Novo score >80.

Sample Processing in-brief:
- Anti-pan-HLA-I (W6/32) and anti-HLA-II monoclonal antibodies were used for immunoaffinity purification of HLA I and II complexes as well as HLA-DRA
Mass spectrometer: timsTOF Pro in PASEF mode and collision induced dissociation and Orbitrap Fusion Lumos mass spectrometer
Data processing:
- Software: PEAKS v10.6
DDA processing: 
- precursor mass tolerance: 20ppm
- fragment mass tolerane: 0.02
- protein digetsion: none
- Fixed modifications: none
- Variable modificaions: methionine oxidation
- Reference: human proteome
- PSM_FDR: 1% 

Summarize information from the following descriptions in list form for use in mass spectrometry database searches. 

Project Description: Cysteine residues undergo various oxidative modifications, acting as key sensors of reactive oxygen species (ROS) and reactive nitrogen species (RNS). Given that ROS and RNS have known roles in many pathophysiological processes, numerous proteome-wide strategies to profile cysteine oxidation state have emerged in the last decade. Recent advancements to traditional redox profiling methods include incorporation of costly isotopic labeling reagents to allow for more quantitative assessment of oxidation states. These methods are typically carried out by using sequential thiol capping and reduction steps in order to label redox-sensitive cysteines, often found in di-cysteine motifs (‘CXXC’ or ‘CXXXC’). Tailored, pricy algorithms are commonly used to analyze redox-profiling datasets, the majority of which cannot accurately quantify site-of-labeling in redox-motifs; moreover, accurate quantification is confounded by excess labeling reagents during sample preparation. Here, we present a low-cost redox-profiling workflow using newly synthesized isotopic reagents compatible with SP3-bead technology, termed SP3-ROx, that allows for high throughput, rapid identification of redox-sensitive cysteines. We optimize cysteine labeling quantification using the FragPipe suite, an open source GUI for MSfragger-based search algorithm. Application of SP3-ROx to naive and activated T cells identifies redox-senstive cysteines, showcasing the utility of this workflow to study biological processes.
Sample Processing Protocol: The samples were analyzed by liquid chromatography tandem mass spectrometry using a Thermo Scientific™ Orbitrap Eclipse™ Tribrid™ mass spectrometer coupled with a High Field Asymmetric Waveform Ion Mobility Spectrometry (FAIMS) Interface. Peptides were resuspended in 5% formic acid and fractionated online using a 18cm long, 100 μM inner diameter (ID) fused silica capillary packed in-house with bulk C18 reversed phase resin (particle size, 1.9 μm; pore size, 100 Å; Dr. Maisch GmbH). The 70-minute water acetonitrile gradient was delivered using a Thermo Scientific™ EASY-nLC™ 1200 system at different flow rates (Buffer A: water with 3% DMSO and 0.1% formic acid and Buffer B: 80% acetonitrile with 3% DMSO and 0.1% formic acid). The detailed gradient includes 0 – 5 min from 3 % to 10 % at 300 nL/min, 5 – 64 min from 10 % to 50 % at 220 nL/min, and 64 – 70 min from 50 % to 95 % at 250 nL/min buffer B in buffer A. Data was collected with charge exclusion (1, 8,>8). Data was acquired using a Data-Dependent Acquisition (DDA) method comprising a full MS1 scan (Resolution = 120,000) followed by sequential MS2 scans (Resolution = 15,000) to utilize the remainder of the 1 second cycle time. Time between master scans was set 1 s. HCD collision energy of MS2 fragmentation was 30 %.
Data Processing Protocol: Raw data collected by LC-MS/MS and converted to mzML format with peakPicking for MS levels 1 and 2 using MSConvert (ProteoWizard release 3.0.20287) were searched using FragPipe GUI v16.0 with MSFragger (version 3.3) [3-7], Philosopher (version 4.0.0) and IonQuant (version 1.7.5) enabled. Precursor and fragment mass tolerance was set as 20 ppm. Missed cleavages were allowed up to 2. Peptide length was set 7 - 50 and peptide mass range was set 500 - 5000. Cysteine residues were searched with variable modifications at cysteine residues for carboxyamidomethylation (+57.02146), IPIAA-L (+463.2366), and IPIAA-H (+467.2529) labeling allowing for 3 max occurrences and all mods used in first search checked. Permissive IonQuant parameters allowed minimum scan/isotope numbers set to 1. PTM-prophet information was obtained from psm.tsv using ‘heavy’ and ‘light’ localizations scores.

Sample Processing in-brief:
Mass Spectrometer: Thermo Scientific™ Orbitrap Eclipse™ Tribrid™ mass spectrometer coupled with a High Field Asymmetric Waveform Ion Mobility Spectrometry (FAIMS) Interface
Data Processing:
- Software: FragPipe v16, Philosopher (version 4.0.0) and IonQuant (version 1.7.5) enabled
- Covert to mzML
- precursor mass tolerance: 20ppm
- fragment mass tolerane: 20ppm
- Maximum peptide length = 50 
- Minimun peptide length = 7  
- Fixed modifications: none
- Variable modifications: methionine oxidation, n-term acetylation, Cysteine carbidomethylation, IPIAA-Light, IPIAA-Heavy
- msfragger.table.var-mods=15.994900,M,true,3; 42.010600,[^,false,3; 57.021460,C,true,3; 463.2366,C,true,3; 467.2529,C,true,3
- IonQuant parameters allowed minimum scan/isotope numbers set to 1

Summarize information from the following descriptions in list form for use in mass spectrometry database searches. 

Project Description: Chemoproteomics investigates small molecule-protein interactions and has made significant progress in recent years. Despite its vast potential, the proteome-wide profiling of reactive cysteine ligandability remains a formidable task to adapt for high throughput applications. This is primarily due to a lack of platforms capable of achieving the desired depth using low sample input in 96- or 384-well plates. Here we have revamped the cysteine profiling platform to address the challenge with an eye toward performing high-throughput library screening in plates. By incorporating several changes including i) an 18-plex TMT sample multiplexing strategy, ii) a magnetic beads-based one-pot workflow, iii) a 10X higher capacity streptavidin resin, and iv) optimized mass spectrometry analyses, a plate-based platform was developed that enables routine interrogation of either ~18,000 or ~24,000 reactive cysteines based on starting amounts of 10 or 20 µg, respectively. We applied the platform to screen a library of 192 electrophiles in the native HEK293T proteome, mapping the ligandablity of 38,450 reactive cysteines from 8,274 human proteins. The significantly improved depth revealed many previously unknown reactive cysteines and cysteine-ligand interactions and led to the identification of an azepane-containing acrylamide which has preferential binding to cysteines in EGF-like domains. We further applied the platform to characterize new cellular targets of well-studied compounds and covalent drugs in three different human cell lines. We found that ARS-1620, a KRASG12C inhibitor, also binds to cysteine 140 of an off-target adenosine kinase ADK, inhibiting its kinase activity. The platform represents a major step forward to high throughput evaluation of reactive cysteines on a proteome-wide scale.
Sample Processing Protocol: Cell lysate was diluted to 2 µg/µL with lysis buffer. 5 or 10 µL lysate, containing 10 or 20 µg protein respectively, was loaded into each well in 96-well plate. The following steps are used for 20 µg protein input. Amount of each reagent should be adjusted accordingly depending on the input. 5 µL of compound solution in lysis buffer was added to a final concentration of 50 µM and incubated for 1 hour. 5 µL of DBIA solution in lysis buffer was added to a concentration of 500 µM and incubated in dark for 1 hour. Addition of 3 µL SP3 beads (1:1 mixture of hydrophobic and hydrophilic type, 50 mg/mL, Cat. #45152105050250 and Cat. #65152105050250,) was followed by addition of 30 µL ~98% ethanol supplemented with 20 mM DTT. Lysate-bead mixture was incubated for 10 minutes with mild shaking before placing the plate on a magnetic stand to aspirate the supernatant. Beads were washed once with 80% ethanol and resuspended in 30 µL lysis buffer supplemented with 20 mM IAA and incubated in dark for 30 minutes with vigorous shaking. 60 µL ~98% ethanol supplemented with 20 mM DTT was added to the mixture. Mild shaking was performed before two washes using 80 % ethanol. The aqueous phase was then removed followed by adding 30 µL 200 mM EPPS buffer (pH 8.5) containing 0.2 µg Lys-C. After 3-hour incubation at room temperature, 5 µL EPPS buffer containing 0.2 µg trypsin was added and incubated with beads at 37 °C overnight. To the mixture of digested peptides and beads, 11 µL acetonitrile and 4 µL TMT (10 µg/µL) reagent were sequentially added, followed by gentle mixing for 60 minutes at room temperature. Reaction was quenched by adding 7 µL 5% hydroxyl amine. Aqueous phase containing TMT-labeled peptides were combined and dried using a SpeedVac. The resulting sample was desalted using a 100-mg Sep-Pak column. Desalted TMT-labeled peptides were resuspended in 310 µL100 mM HEPES buffer (pH 7.4). Addition of 50 µL Pierce™ High Capacity Streptavidin Agarose (Cat. #20359) was followed by incubation of peptide-beads mixture at room temperature for 3 hours. Resulting mixture was then loaded on a Ultrafree-MC centrifugal filter (hydrophilic PTFE, 0.22 µm pore size) and centrifugated at 1000 g for 30 seconds. Beads were washed sequentially with 300 µL 100 mM HEPES (pH 7.4) with 0.05% NP-40 twice, 350 µL 100 mM HEPES (pH 7.4) three times and 400 µL water once. Peptide were eluted sequentially by 1) incubation with elution buffer (80% acetonitrile, 0.1% formic acid) for 20 minutes at room temperature; 2) incubation with elution buffer for 10 minutes at room temperature; 2) incubation with elution buffer for 10 minutes at 72 °C. The combined eluent was dried in a SpeedVac and desalted via StageTip prior to LC-MS/MS analysis.
Data Processing Protocol: Raw files of cysteine profiling were searched using the open-source Comet search engine (ver. 2019.01.5)54 with the Uniprot human proteome database (downloaded 11/24/2021) with contaminants and reverse decoy sequences appended. Precursor error tolerance was 50 p.p.m. and fragment error tolerance was 0.02 Da. Static modifications include Cys carboxyamidomethylation (+57.0215) and TMTpro18 (+304.2071) on Lys side chains and peptide N-termini. Methionine oxidation (+15.9949) and DBIA-modification on cysteine residues (+239.1628) were allowed as variable modifications.  Peptide spectral matches were filtered to a peptide false discovery rate (FDR) of <1% using linear discriminant analysis employing a target-decoy strategy55,56. Resulting peptides were further filtered to obtain a 1% protein FDR at the entire dataset level (including all plexes in an experiment)57. Cysteine-modified peptides were filtered for site localization using the AScorePro algorithm with a cutoff of 13 (P < 0.05) as previously described58,59. Overlapping peptide sequences generated from different charge states, retention times and tryptic termini were grouped together into a single entry. A single quantitative value was reported, and only unique peptides were reported. Reporter ion intensities were adjusted to correct for impurities during synthesis of different TMT reagents according to the manufacturer’s specifications. For quantification of each MS2 spectrum, a total sum signal-to-noise of all reporter ions of 180 (TMTPro 18-plex) was required. Peptide quantitative values were normalized so that the sum of the signal for all proteins in each channel was equal to account for sample loading differences (column normalization).

Sample Processing in-brief:
- Sample: HEK293T cell lysate, 10 or 20 µg protein per well in a 96-well plate
- Compound incubation: 50 µM for 1 hour
- DBIA incubation: 500 µM for 1 hour in dark
- SP3 beads addition: 3 µL of a 1:1 mixture of hydrophobic and hydrophilic beads
- Reduction and alkylation: 30 µL ~98% ethanol with 20 mM DTT followed by 30 µL lysis buffer with 20 mM IAA
- Digestion: 3 hours with Lys-C, overnight with trypsin at 37 °C
- TMT labeling: 11 µL ACN, 4 µL TMT reagent, 60 minutes at room temperature
- Quenching: 7 µL 5% hydroxylamine
- Streptavidin purification: 50 µL High Capacity Streptavidin Agarose, incubated 3 hours
- Elution: Sequential elutions at room temperature and 72°C with 80% ACN, 0.1% FA
- Desalting: StageTip prior to LC-MS/MS analysis
Mass Spectrometer:
- Instrument: Not specified in the summary, assume high-resolution LC-MS/MS
Data Processing in-brief:
- Software: Comet (ver. 2019.01.5)
- Database: Uniprot human proteome (downloaded 11/24/2021) with contaminants and reverse decoy sequences
- Precursor mass tolerance: 50 ppm
- Fragment mass tolerance: 0.02 Da
- Fixed modifications: Cys carboxyamidomethylation (+57.0215), TMTpro18 (+304.2071) on Lys and peptide N-termini
- Variable modifications: Methionine oxidation (+15.9949), DBIA-modification on cysteine residues (+239.1628)
- FDR: Peptide <1%, Protein 1% using linear discriminant analysis and target-decoy strategy
- Quantification: MS2 spectrum total sum S/N of all reporter ions of 180 required, normalized for sample loading differences (column normalization)
- Site localization: AScorePro algorithm with a cutoff of 13 (P < 0.05) for cysteine-modified peptides
- Data Adjustment: Reporter ion intensities adjusted for TMT reagent impurities

Summarize information from the following descriptions in list form for use in mass spectrometry database searches. 

Project Description: One key barrier to improving efficacy of personalized cancer immunotherapies that are dependent on the tumor antigenic landscape remains patient stratification. While patients with CD3+CD8+ T cell inflamed tumors typically show better response to immune checkpoint inhibitors, it is still unknown if the repertoire of HLA bound peptides presented in highly inflamed and non-inflamed tumors is substantially different. We surveyed 61 tumor regions and adjacent non-malignant lung tissues from eight lung cancer patients and performed deep antigen discovery combining immunopeptidomics, genomics, bulk and spatial transcriptomics and explored the heterogeneous expression and presentation of tumor (neo)antigens. Here we associated diverse immune cell populations with the immunopeptidome in CD3+CD8+ T cell excluded, infiltrated, highly- and lowly-inflamed tumors. We found evidence for lower immune-editing and higher presentation efficiency of tumor-associated antigens in lowly-inflamed and CD3+CD8+ T cell excluded tumors. This could have implications for the choice of combination therapies tailored to the patient’s mutanome and microenvironment.
Sample Processing Protocol: Immunoaffinity purification of HLA peptides We performed HLA immunoaffinity purification according to our previously established protocols. W6/32 and HB145 monoclonal antibodies were purified from the supernatants of HB95 (ATCC® HB-95™) and HB145 cells (ATCC® HB-145™) using protein-A sepharose 4B (Pro-A) beads (Invitrogen), and antibodies were then cross-linked to Pro-A beads. Snap-frozen tissue samples were lysed in lysis buffer and homogenized on ice in 3–5 short intervals of 5 s each using an Ultra Turrax homogenizer (IKA), and were left for 1h at 4 °C. The lysis buffer contained PBS with 0.25% sodium deoxycholate (Sigma Aldrich), 0.2 mM iodoacetamide (Sigma Aldrich), 1 mM EDTA, a 1:200 protease inhibitors cocktail (Sigma Aldrich), 1 mM phenylmethylsulfonylfluoride (Roche), and 1% octyl-beta-D glucopyranoside (Sigma Alrich). The lysates were then cleared by centrifugation at 75,600 x g in a high-speed centrifuge (Beckman Coulter, Avanti JXN-26 Series, JA-25.50 rotor) at 4 °C for 50 min. For HLA immunopurification, we employed the Waters Positive Pressure-96 Processor (Waters) and 96-well single-use filter micro-plates with 3 µm glass fibers and 25 µm polyethylene membranes (Agilent, 204495-100). The lysates were passed sequentially through the first plate with wells containing Pro-A beads for endogenous antibody depletion, then a second plate with wells containing anti-pan HLA-I antibody-cross-linked beads, and lastly through the third plate containing anti-pan HLA-II antibody-cross-linked beads. The HLA-I and HLA-II crosslinked beads were then washed separately in their plates with varying concentrations of salts using the processor. Finally, the beads were washed twice with 2 mL of 20 mM Tris-HCl pH 8. Sep-Pak tC18 100 mg Sorbent 96-well plates (Waters, ref no: 186002321) were used for the purification and concentration of HLA-I and HLA-II peptides. The C18 sorbents were conditioned, and the HLA complexes and bound peptides were directly eluted from the filter plates with 1% trifluoroacetic acid (TFA; Sigma Aldrich). After washing the C18 sorbents with 2 mL of 0.1% TFA, HLA-I peptides were eluted with 28% acetonitrile (ACN; Sigma Aldrich) in 0.1% TFA, and HLA-II peptides were eluted with 32% ACN in 0.1% TFA. Recovered HLA-I and -II peptides were dried using vacuum centrifugation (Concentrator plus, Eppendorf) and stored at −20 °C. Liquid chromatography–mass spectrometry (LC-MS/MS) analyses The LC-MS/MS system consisted of an Easy-nLC 1200 (Thermo Fisher Scientific) connected to a Q Exactive HF-X mass spectrometer (Thermo Fisher Scientific). Peptides were separated on a 450 mm analytical column (8 µm tip, 75 µm inner diameter, PicoTipTMEmitter, New Objective) packed with ReproSil-Pur C18 (1.9 µm particles, 120 Å pore size, Dr. Maisch GmbH). The separation was performed at a flow rate of 250 nL/min by a gradient of 0.1% formic acid (FA) in 80% ACN (solvent B) in 0.1% FA in water (solvent A). HLAIp were analyzed by the following gradient: 0–5 min (5% B); 5–85 min (5–35% B); 85–100 min (35–60% B); 100–105 min (60–95% B); 105–110 min (95% B); 110–115 min (95–2% B) and 115–125 min (2% B). HLAIIp were analyzed by the following gradient: 0–5 min (2–5% B); 5–65 min (5–30% B); 65–70 min (30–60% B); 70–75 min (60–95% B); 75–80 min (95% B), 80–85 min (95–2% B) and 85–90 min (2% B). For data dependent acquisition, full MS spectra were acquired in the Orbitrap from m/z = 300–1650 with a resolution of 60,000 (m/z = 200) and an ion accumulation time of 80 ms. The auto gain control (AGC) was set to 3e6 ions. MS/MS spectra were acquired in a data-dependent manner on the 20 most abundant precursor ions (if present) with a resolution of 15,000 (m/z = 200), an ion accumulation time of 120 ms, and an isolation window of 1.2 m/z. The AGC was set to 2e5 ions, the dynamic exclusion was set to 20 s, and a normalized collision energy (NCE) of 27 was used for fragmentation. No fragmentation was performed for HLAIp with assigned precursor ion charge states of four and above or for HLAIIp with an assigned precursor ion charge state of one, or six and above. The peptide match option was disabled. Data was acquired in an independent manner using the following parameters. The cycle of acquisition consisted of a full MS scan from 300-1650 m/z (R = 60,000 and ion accumulation time of 60 ms) and 21 DIA MS/MS scans in the Orbitrap. For each DIA MS/MS scan, a resolution of 30,000, an AGC of 3e6 and a ramping normalized collision energy (NCE = 25.5, 27 and 30) were used. The maximum ion accumulation was set to auto and the overlap between consecutive MS/MS scans was 1 m/z. For some peptide samples, only DDA measurements were available.
Data Processing Protocol: MS-based searches First, for each region, we search the corresponding MS raw file against the a personalized proteome reference as described above using Comet with the following setup: precursor mass tolerance 20 ppm, MS/MS fragment tolerance of 0.02 Da, peptide length of 8–15 when searching only HLA-I peptides and 8–25 for both HLA-I and HLA-II peptides, no fixed modifications, while methionine oxidation and phosphorylation on serine, threonine and tyrosine was included as a variable modification. A group specific 3% FDR (smoothing value: 5) for protein-coding, non-canonical sources and transposable elements) was calculated by NewAnce v1.7.1 (smoothing value: 5) as previously described. We next generated a single comprehensive reference database containing all the protein sequence sources of the detected personalized variant and non-canonical peptides from all patient tissue regions, and concatenate these to a generic GENCODE database. Then, Comet and NewAnce searches were run again against this database using the entire immunopeptidomics dataset, yet separately for HLA-I and HLA-II raw files with same parameters as above. The outputs of this search were used to create spectral libraries for targeted DIA analysis with Spectronaut. The spectral libraries were generated using the results from the Comet and NewAnce search by parsing the PSMs into the BGS generic format recommended by Spectronaut (version 14.6.2, Biognosys). The default settings were used except the following parameters. For the mass tolerance, the calibration search was set to ‘dynamic’ and a correction factor of 1 was used for MS and MS/MS. For identification, a FDR threshold of 0.01 and unspecific digestion rule were used. For spectral library filtering, only peptides with longer than 3 a.a and with at least 3 best fragments (up to 6) were considered. A deep learning assisted iRT regression and a R-square of 0.8 were used for iRT reference strategy and correlation score, respectively. The selection of fragment ions based on intensity was used. For targeted DIA based identification of the peptides, the library was matched against the corresponding immunopeptidomics DIA raw files with q-value cut-off of 0.01 and 1 respectively for precursor and protein. MS and MS/MS data were extracted using maximum intensity strategy within the given mass tolerance. A dynamic mass tolerance and a correction factor of 1 were used for both MS and MS/MS. The quantification was done at MS/MS level, using at least 3 fragments. The global normalization per run based on the median was switched off to avoid cross run normalization. Results from Spectronaut were exported in peptide-centered file formats: peptide quantity and peptide score. For more extended analysis of HLA-I peptides derived from non-canonical sources, we used the high-confidence database of translated nuORFs across tissues (nuORFdb) (concatenated with the human reference proteome; 323,848 entries, PA_nuORFdb_v1.0.fasta) and a reduced version of the above mentioned personalized references per patient, where the ORFs non-canonical sources were restricted to in silico methionine-to-stop translated transcript entries, resulting in fasta files with overall similar size per patient (ranging from 521,779 entries for patient 02672 to 599,300 entries for patient 02287). We then used the hybrid DIA approach with Spectronaut (version 16.3: The parameters are provided in Supplementary software 1). Peptide identification was performed by PulsarTM on DIA and DDA files (separately per patient) using unspecific digestion and with a peptide length from 8 – 15a.a. Acetylation at Protein N-term and Oxidation of methionine were considered as variable modifications. –b and y ions were used as main ion types. For non-canonical sources, we used the gencode annotation with priorities in the following order: lncRNAs, processed transcripts, pseudogenes, retained introns, noncanonical ORFs and ‘others’.

Sample Processing in-brief:
- Lysates passed sequentially to HLA I and HLA II pan antibodies 
- W6/32 and HB145 abtibodies used for HLA I immunoaffinity purfication
- C18 cleanup, peptide elution, and vacuum drying
Mass Spectrometer:
- Instrument:  Easy-nLC 1200 connected to a Q Exactive HF-X mass spectrometer 
Data Processing in-brief:
- Software:Comet 
- Reference: personalized proteome references
- Precursor mass tolerance: 20 ppm
- Fragment mass tolerance: 0.02 Da
- Maximum peptide length = 25 
- Minimun peptide length = 8  
- Fixed modifications: none
- Variable modifications: methionine oxidation, phosphorylation on serine, threonine and tyrosine
- group specific 3% FDR (smoothing value: 5) for protein-coding, non-canonical sources and transposable elements) was calculated by NewAnce v1.7.1
DIA processing:
- Spectral library: entire immunopeptidomics dataset, yet separately for HLA-I and HLA-II raw files with same parameters as above with Comet and NewAnce
- Reference: single comprehensive database containing all the protein sequence sources of the detected personalized variant and non-canonical peptides from all patient tissue regions
- DIA analysis with Spectronaut: PSMs into the BGS generic format recommended by Spectronaut,  the mass tolerance, the calibration search was set to ‘dynamic’ and a correction factor of 1 was used for MS and MS/MS
- FDR threshold of 0.01
- Protein digestion: unspecific 
Hybrid DIA processing:
-Software: Spectronaut
-Protein digestion: non-specific
- Maximum peptide length = 15 
- Minimun peptide length = 8 
- Fixed modifications: none
- Variable modifications: n-term acetylation, methionine oxidation
'''),
    ("user", "{input}")
])

chain = prompt | llm 



# Define the path to your output.jsonl file
jsonl_file = "output.jsonl"


# Open the output.jsonl file
with jsonlines.open(jsonl_file) as reader:
    # Iterate through each line in the file
    for obj in reader:
        # Extract project title, accession, project_description, sample_processing_protocol, and data_processing_protocol
        project_title = obj["project_title"]
        project_accession = obj["project_accession"]
        project_description = obj["project_description"]
        sample_processing_protocol = obj["sample_processing_protocol"]
        data_processing_protocol = obj["data_processing_protocol"]

        # Construct the query using the extracted text
        query = f'''Summarize information from the following descriptions in list form for use in mass spectrometry database searches. 

Project Description: {project_description}
Sample Processing Protocol: {sample_processing_protocol}
Data Processing Protocol: {data_processing_protocol}
'''

        # Invoke the chain for the current project description
        result = chain.invoke({"input": query})
        
        # Print the project title and accession number
        print("Project Title:", project_title)
        print("Accession:", project_accession)
        print()
        
        # Print or use the result as needed
        print(result.content)
        print()

