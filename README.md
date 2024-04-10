# pride_mine

Mine the PRIDE repository for projects matching keyword(s). \
Output descriptions are input for GPT-4 to generate concise summaries of information relevant for data searching. \
The FTP locations of relevant projects can be used in downstream processing.

![Model](pride_mine.png)


## To run
- Obtain an Open AI Key for line 77 of get_summaries.py
- python get_summaries.py "keywords"

## Output
Project tiles, accessions, and summaries are printed.
