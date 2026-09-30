# Restricted analysis data

The individual-level Busselton Health Study (BHS) analysis data required for the analyses in this repository are not publicly distributed.

BHS genotype, phenotype, pregnancy-history and linked health data are subject to data-access, governance and privacy restrictions.

The Fine-Gray analysis script expects an authorised local analysis dataset containing the variables required for the models.

### Follow-up time

The `follow_up_years` variable represents time from the 1994–95 BHS baseline survey to the relevant time-to-event endpoint. For participants who experienced incident CVD, follow-up ended at the first recorded CVD event. For participants without incident CVD, follow-up ended at the recorded date of death, where available, or at the administrative end of follow-up, defined as 27 years after baseline (end of 2022).

For the Fine-Gray analyses of specific first CVD outcomes, the CVD subtype of interest was treated as the event of interest, while alternative first CVD outcomes were treated as competing events. Participants without an incident CVD event were censored at their defined follow-up time.

The public repository does not contain:

- individual participant records
- participant identifiers
- genotype data
- individual polygenic risk score values
- linked hospitalisation records
- mortality records
- pregnancy-history records
- participant ID mapping files
Researchers interested in accessing BHS data should follow the relevant BHS data-access and governance procedures.
